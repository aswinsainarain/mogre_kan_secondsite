# -*- coding: utf-8 -*-

import pandas
import subprocess

# For each .readsPerGene file
files = subprocess.Popen('ls *.readsPerGene', stdout = subprocess.PIPE, shell = True).communicate()
files = files[0].split()

geneCounts = pandas.DataFrame()
colnames = list()
for file in files:
    colnames.append(file.split('.')[0])
    with open(file) as readsPerGene:
        temp = pandas.read_table(readsPerGene, index_col = 0)
        geneCounts = pandas.concat([geneCounts, temp], axis = 1, ignore_index = 1)
geneCounts.columns = colnames

geneCounts.to_csv('geneReadCounts.txt', sep = '\t', na_rep = '0', index_label = 'Genes')
