import collections
import subprocess
import re
import pandas as pd

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

# For each .readsPerNt
files = subprocess.Popen('ls *.sorted.sam', stdout = subprocess.PIPE, shell = True).communicate()
files = files[0].split()

'''
Mostly you'll see three different flags
16: read reverse strand
0: I don't know what I guess forward strand
4: Read unmapped
Descriptions of the flags can be found on this website:https://broadinstitute.github.io/picard/explain-flags.html
Or you could consult the sam format specification and then convert hexadecimal to decimal

You can use the following code on a samfile to figure out how many of these are there:
temp = collections.defaultdict(int)    
with open(file) as samfile:
    tempcounter = 0
    for line in samfile:
        tempcounter += 1
        a = re.split(r'\t', line)
        temp[int(a[1])] += 1

'''
for file in files:
    name = file.split(r'.')[0]
    output = name + '.perGeneFRCounts'  
    # Constructing the counts dictionary
    forward = collections.defaultdict(int)
    reverse = collections.defaultdict(int)
    # Populating the dictionary counts with gene based read counts
    with open(file) as samfile:
        for line in samfile:
            a = re.split(r'\t', line)   # a[3] is coordinate; a[1] is sam flag
            if int(a[4]) >= 20: # a[4] is mapping quality. If above 20 count whether it maps to forward or reverse strand
                if int(a[3]) in genes:
                    if int(a[1]) == 0: # If read maps to forward strand
                        forward[genes[int(a[3])]] += 1
                    if int(a[1]) == 16: # If read maps to reverse strand
                        reverse[genes[int(a[3])]] += 1
    forwardDF = pd.DataFrame.from_dict(forward, orient = 'index')
    reverseDF = pd.DataFrame.from_dict(reverse, orient = 'index')
    mergedDF = pd.concat([forwardDF, reverseDF], axis = 1)
    mergedDF = mergedDF.fillna(value = 0)
    mergedDF.columns = ['forward', 'reverse']
    mergedDF.to_csv(output, sep = '\t')