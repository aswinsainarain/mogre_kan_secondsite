import subprocess
import re

sams = subprocess.Popen('ls ../*sorted.sam', stdout = subprocess.PIPE, shell = True).communicate()
sams = sams[0].split()

for sam in sams:
    output = re.split(r'[/.]', sam)[3] + '.unmapped.fasta'    
    with open(output, 'w') as o:
        with open(sam) as samFile:
            for line in samFile:
                a = re.split(r'\t', line)
                header = a[0]
                read = a[9]
                genome = a[2]
                if genome == '*':
                    o.write('>' + header + '\n')
                    o.write(read + '\n')
