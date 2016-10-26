import subprocess
import collections
import pandas as pd
import re


# Annotations
genes = collections.defaultdict()
with open('mg1655.annotation') as ann:
    seen = 0
    for line in ann:
        a = line.split()
        if seen == 1:
            for i in range(int(a[1]), int(a[2])+1):
                genes[i] = a[0] # coordinate is given value gene name
        seen = 1


snpFiles = subprocess.Popen('ls *.snp', stdout = subprocess.PIPE, shell = True).communicate()[0].split()
indelFiles = subprocess.Popen('ls *.indel', stdout = subprocess.PIPE, shell = True).communicate()[0].split()

# Snp Table Generation
snpTable = pd.DataFrame()
for snpFile in snpFiles: # For all the .snp files
    tempDict = collections.defaultdict() # reinitialize tempDict
    with open(snpFile) as snp:
        i = 0
        for line in snp:
            if i == 1:
                line = line.rstrip()
                a = re.split(r'\t', line)
                if int(a[1]) in genes:
                    geneName = genes[int(a[1])]
                else:
                    geneName = 'IGR'
                tempKey = geneName + ' ' + a[1] + ' ' + a[2] + ' -> ' + a[-1] 
                tempValue = re.split(r'%', a[6])[0]
                tempValue = float(tempValue)
                tempDict[tempKey] = tempValue
            i = 1
        tempDF = pd.DataFrame.from_dict(tempDict, orient = 'index')
    snpTable = pd.concat([snpTable, tempDF], axis = 1)

snpTable.columns = snpFiles
snpTable.to_csv('Final_snp_list.txt', sep = '\t', na_rep = 'NA', index_label = 'Annotation')

# snp table with cutoff 90% in at least one of the samples
retainRows = list()
for i in range(len(snpTable)):
    for varfreq in snpTable.ix[i]:
        if varfreq >= 90:
            retainRows.append(i)
            break

snpTable90 = snpTable.ix[retainRows]
snpTable90.to_csv('cutoff_90_snp_list.txt', sep = '\t', na_rep = 'NA', index_label = 'Annotation')



# snp table with cutoff 50% in at least one of the samples
retainRows = list()
for i in range(len(snpTable)):
    for varfreq in snpTable.ix[i]:
        if varfreq >= 50:
            retainRows.append(i)
            break

snpTable50 = snpTable.ix[retainRows]
snpTable50.to_csv('cutoff_50_snp_list.txt', sep = '\t', na_rep = 'NA', index_label = 'Annotation')

##################

# Indel Table Generation
indelTable = pd.DataFrame()
for indelFile in indelFiles: # For all the .snp files
    tempDict = collections.defaultdict() # reinitialize tempDict
    with open(indelFile) as indel:
        i = 0
        for line in indel:
            if i == 1:
                line = line.rstrip()
                a = re.split(r'\t', line)
                if int(a[1]) in genes:
                    geneName = genes[int(a[1])]
                else:
                    geneName = 'IGR'
                tempKey = geneName + ' ' + a[1] + ' ' + a[2] + ' -> ' + a[-1] 
                tempValue = re.split(r'%', a[6])[0]
                tempValue = float(tempValue)
                tempDict[tempKey] = tempValue
            i = 1
        tempDF = pd.DataFrame.from_dict(tempDict, orient = 'index')
    indelTable = pd.concat([indelTable, tempDF], axis = 1)

indelTable.columns = indelFiles
indelTable.to_csv('Final_indel_list.txt', sep = '\t', na_rep = 'NA', index_label = 'Annotation')

# snp table with cutoff 90% in at least one of the samples
retainRows = list()
for i in range(len(indelTable)):
    for varfreq in indelTable.ix[i]:
        if varfreq >= 90:
            retainRows.append(i)
            break

indelTable90 = indelTable.ix[retainRows]
indelTable90.to_csv('cutoff_90_indel_list.txt', sep = '\t', na_rep = 'NA', index_label = 'Annotation')



# indel table with cutoff 50% in at least one of the samples
retainRows = list()
for i in range(len(indelTable)):
    for varfreq in indelTable.ix[i]:
        if varfreq >= 50:
            retainRows.append(i)
            break

indelTable50 = indelTable.ix[retainRows]
indelTable50.to_csv('cutoff_50_indel_list.txt', sep = '\t', na_rep = 'NA', index_label = 'Annotation')
