# cutadapt version 1.9.dev1 with Python 2.7.4
# bwa Version: 0.7.5a-r405
# samtools Version: 1.3 (using htslib 1.3)
# VarScan version 2.3.7

#####################################
fastqGzFiles <-list.files(path = "fastq_gz_files/",  # Use this part of the script to extract all the fastq.gz files
                       pattern = ".fastq.gz$")

for (i in fastqGzFiles) {
 filename = strsplit(i, split = "[.]")[[1]][1]
 print(paste0("COMMAND: gzip -dc fastq_gz_files/", i, " > fastq_files/", filename, ".fastq"))
 system(paste0("gzip -dc fastq_gz_files/", i, " > fastq_files/", filename, ".fastq"))
}
######################################

fastqFiles <-list.files(path = "fastq_files/",   # all .fastq files
                          pattern = ".fastq$")

system("mkdir trimmed_fastqs")

# Cutadapt
# Cut reads shorter than 30 bases are discarded by the -m 30 flag
trimmedFastqs <- NULL
for (i in fastqFiles) {
  filename <- strsplit(i, split = "[.]")[[1]][1]
  trimmedFile <- paste0(filename, ".trimmed.fastq")
  trimmedFastqs <- c(trimmedFastqs, trimmedFile)
  print(paste0("COMMAND: cutadapt -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -b TCTCGTATGCCGTCTTCTGCTTG -m 30 -o trimmed_fastqs/", 
                trimmedFile, 
                " fastq_files/", 
                i))
   system(paste0("cutadapt -b GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -b TCTCGTATGCCGTCTTCTGCTTG -m 30 -o trimmed_fastqs/", 
                 trimmedFile, 
                 " fastq_files/", 
                 i))
}

trimmedFastqs <- list.files(path = "trimmed_fastqs/",
                            pattern = ".fastq$")


# For all fastq files align them, convert them to bam, generate pileup and call snps and indels using VarScan
# bwa index the genome
print("COMMAND: bwa index -a is NC_000913.fna")
system("bwa index -a is NC_000913.fna")
# bwa aln
# -q 30 flag trims ends based on Phred quality 30. There is some formula for choosing which bases to trim. Check the bwa manual.
for (i in trimmedFastqs) {
    filename <- strsplit(i, split = "[.]")[[1]][1]
    print(paste0("COMMAND: bwa aln -t 7 -q 30 NC_000913.fna trimmed_fastqs/", i, " > ", filename, ".sai"))
    system(paste0("bwa aln -t 7 -q 30 NC_000913.fna trimmed_fastqs/", i, " > ", filename, ".sai"))
}
# bwa samse
for (i in trimmedFastqs) {
    filename <- strsplit(i, split = "[.]")[[1]][1]
    print(paste0("COMMAND: bwa samse NC_000913.fna ", filename, ".sai trimmed_fastqs/", i, " > ", filename, ".sam"))
    system(paste0("bwa samse NC_000913.fna ", filename, ".sai trimmed_fastqs/", i, " > ", filename, ".sam"))
}
# samtools view -bS
for (i in trimmedFastqs) {
    filename <- strsplit(i, split = "[.]")[[1]][1]
    print(paste0("COMMAND: samtools view -bS ", filename, ".sam > ", filename, ".bam"))
    system(paste0("samtools view -bS ", filename, ".sam > ", filename, ".bam"))
}
# samtools sort
for (i in trimmedFastqs) {
    filename <- strsplit(i, split = "[.]")[[1]][1]
    print(paste0("COMMAND: samtools sort ", filename, ".bam -o ", filename, ".sorted.bam"))
    system(paste0("samtools sort ", filename, ".bam -o ", filename, ".sorted.bam"))
}
# samtools view
for (i in trimmedFastqs) {
    filename <- strsplit(i, split = "[.]")[[1]][1]
    print(paste0("COMMAND: samtools view ", filename, ".sorted.bam > ", filename, ".sorted.sam"))
    system(paste0("samtools view ", filename, ".sorted.bam > ", filename, ".sorted.sam"))
}
# samtools mpileup
for (i in trimmedFastqs) {
    filename <- strsplit(i, split = "[.]")[[1]][1]
    print(paste0("COMMAND: samtools mpileup -f NC_000913.fna -B ", filename, ".sorted.bam > ", filename, ".pileup"))
    system(paste0("samtools mpileup -f NC_000913.fna -B ", filename, ".sorted.bam > ", filename, ".pileup"))
}
# java -jar VarScan pileup2snp
for (i in trimmedFastqs) {
    filename <- strsplit(i, split = "[.]")[[1]][1]
    print(paste0("COMMAND: java -jar VarScan.v2.3.7.jar pileup2snp ", filename, ".pileup --min-avg-qual 30 > ", filename, ".snp"))
    system(paste0("java -jar VarScan.v2.3.7.jar pileup2snp ", filename, ".pileup --min-avg-qual 30 > ", filename, ".snp"))
}
# java -jar VarScan pileup2indel
for (i in trimmedFastqs) {
    filename <- strsplit(i, split = "[.]")[[1]][1]
    print(paste0("COMMAND: java -jar VarScan.v2.3.7.jar pileup2indel ", filename, ".pileup --min-avg-qual 30 > ", filename, ".indel"))
    system(paste0("java -jar VarScan.v2.3.7.jar pileup2indel ", filename, ".pileup --min-avg-qual 30 > ", filename, ".indel"))
}
