import subprocess
import re

# For each .fastq
files = subprocess.Popen('ls ../trimmed_fastqs/*.fastq', stdout = subprocess.PIPE, shell = True).communicate()
files = files[0].split()

for file in files:
    print file + ':'
    a = re.compile(r'@')
    n = re.compile(r':N:')
    y = re.compile(r':Y:')
    Ncounts = 0
    Ycounts = 0
    with open(file) as fastq:
        for line in fastq:
            if re.match(a, line):
                if re.search(n, line):
                    Ncounts += 1
                if re.search(y, line):
                    Ycounts += 1
    print 'N counts = ' + str(Ncounts)
    print 'Y counts = ' + str(Ycounts)