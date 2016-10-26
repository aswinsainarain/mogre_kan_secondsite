# For each sorted sam file
import subprocess
files = subprocess.Popen('ls *.sorted.sam', stdout = subprocess.PIPE, shell = True).communicate()
files = files[0].split()

import collections
for file in files:
    name = file.split(r'.')[0]
    output = name + '.readsPerNt'
    
    # Constructing the counts dictionary with genome coordinates
    counts = collections.defaultdict(int)
    
    # Populating the dictionary counts with position based read counts
    with open(file) as sam:
        for line in sam:
            a = line.split('\t')
            if int(a[4]) >= 20:   # mapping quality is set here
                counts[int(a[3])] += 1 # coordinate read count is incremented by 1
    o = open (output, 'w')
    for coordinate in sorted(counts):
        write = str(coordinate) + '\t' + str(counts[coordinate]) + '\n'
        o.write(write)
    o.close()
