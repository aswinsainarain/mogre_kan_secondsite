import collections
import subprocess
import re
import pandas as pd

# Annotations
ptt = pd.DataFrame.from_csv('NC_000913.ptt', header = 2, sep = '\t', index_col = 5) # read ptt file
#rnt = pd.DataFrame.from_csv('NC_000913.rnt', header = 2, sep = '\t', index_col = 5) # read rnt file

#ann = pd.concat([ptt, rnt]) # concatenate them into a single table
ann = ptt

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

# len(re.split(r'\.', genes[16903])) # This is an example of a coordinate that maps multiple genes hokC is embedded in mokC
# Remove coordinates that have overlapping genes
for coordinate in genes.keys():
    if len(re.split(r'\.', genes[coordinate])) > 1: # if there are mulitple genes mapping to a coordinate remove them
        del genes[coordinate]

# For each .readsPerNt
files = subprocess.Popen('ls *.readsPerNt', stdout = subprocess.PIPE, shell = True).communicate()
files = files[0].split()

for file in files:
    name = file.split(r'.')[0]
    output = name + '.readsPerGene'
    
    # Constructing the counts dictionary
    counts = collections.defaultdict(int)
    
    # Populating the dictionary counts with gene based read counts
    with open(file) as readsPerNt:
        for line in readsPerNt:
            a = line.split('\t')
            if int(a[0]) in genes:
                counts[genes[int(a[0])]] += int(a[1])  # Get the (value) gene name from the coordinate and add the number of reads mapped to that coordinate to a new dictionary called counts
    
    o = open (output, 'w')
    o.write('Gene\tReadCounts\n')
    for gene in sorted(counts):
        write = gene + '\t' + str(counts[gene]) + '\n'
        o.write(write)
    o.close()        
