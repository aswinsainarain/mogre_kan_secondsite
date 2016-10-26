import subprocess
files = subprocess.Popen('ls *.sorted.sam', stdout = subprocess.PIPE, shell = True).communicate()
files = files[0].split()

import collections
for file in files:
    mapq = collections.defaultdict(int)
    with open(file) as sam:
        for line in sam:
            a = line.split('\t')
            mapq[int(a[4])] += 1
    print file
    for quality in sorted(mapq):
        print(str(quality) + '\t' + str(mapq[quality]))
