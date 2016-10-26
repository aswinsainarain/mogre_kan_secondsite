import pandas as pd
import re

ptt = pd.DataFrame.from_csv('NC_000913.ptt', header = 2, sep = '\t', index_col = 5) # read ptt file
rnt = pd.DataFrame.from_csv('NC_000913.rnt', header = 2, sep = '\t', index_col = 5) # read rnt file

ann = pd.concat([ptt, rnt]) # concatenate them into a single table

# Build a dictionary with coordinates as keys and gene names as values
genes = dict()
for i in range(len(ann)):
    start, end = re.split(r'\.\.', ann.ix[i, 0])
    start = int(start)
    end = int(end)
    geneName = ann.ix[i, 4]
    for j in range(start, end+1):
        if j not in genes:
            genes[j] = geneName
        elif j in genes:
            genes[j] += '.' + geneName # If a gene is already assigned to the coordinate, add the new gene name on top of the existing gene name separated by a '.'

# len(re.split(r'\.', genes[16903])) # This is an example of a coordinate that maps muplitple genes hokC is embedded in mokC
# Remove coordinates that have overlapping genes
for coordinate in genes.keys():
    if len(re.split(r'\.', genes[coordinate])) > 1: # if there are mulitple genes mapping to a coordinate remove them
        del genes[coordinate]

# Get the 'new' gene lengths. GENES WITH OVERLAPPING REGIONS WILL BE SHORTER IN LENGTH (since the overlapping region has been removed).
geneLengths = dict()
for coordinate in genes.keys():
    if genes[coordinate] in geneLengths:
        geneLengths[genes[coordinate]] += 1
    else:
        geneLengths[genes[coordinate]] = 1

# Write these gene lengths to a file
with open('gene_lengths.txt', 'w') as o:
    for gene in sorted(geneLengths):
        write = gene + '\t' + str(geneLengths[gene]) + '\n'
        o.write(write)