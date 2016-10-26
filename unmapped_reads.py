# For each sorted sam file
import subprocess
files = subprocess.Popen('ls *sorted.sam', stdout = subprocess.PIPE, shell = True).communicate()
files = files[0].split()


for file in files:
    name = file.split(r'.')[0]

    # Populating the dictionary counts with position based read counts
    with open(file) as sam:
        print 'File:' + file
        count = 0
        count1 = 0
        for line in sam:
            a = line.split('\t')
            if a[2] == '*':    
                count += 1 # coordinate read count is incremented by 1 if coordinate field has *
            if count == 1 and count1 == 0:
                count1 = 1
                print 'Example Unmapped Read:\n' + line.rstrip()
        print 'Unmapped Reads:' + str(count) + '\n'
