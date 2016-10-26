Some datasets required by the script are included in the datasets folder.
The file paths within these scripts will need to be changed since these were made to work for my folder arrangement.

The first script to use is analysis.R which uses the following scripts in order:
1. align_and_snp_call.R (performs alignment of reads and snp calling)
2. annotation_from_ptt_rnt.R (generates some annotation)
3. annotation_from_ptt_rnt.py (generates some annotation)
4. condense_snp_indel_info.py (condenses information of snps and indels across all samples into one table)
5. check_NY_reads.py (checks if there are any reads with mismatches in the barcode)
6. mapq_readcounts.py (gets counts of reads of each match quality)
7. unmapped_reads.py (gets number of unmapped reads)
8. read_count_per_nt.py (gets number of reads mapping to each nucleotide)
9. per_gene_coding_noncoding_read_counts.py (gets number of reads mapping to the + or – strand in each gene)
10. check_forward_reverse_corr.R (checks correlation between numbers of reads mapping to the + and – strand in genes)
11. read_count_per_gene.py (gets the number of reads mapping to each gene)
12. tabulate_per_gene_readcounts.py (tabulates the per gene read counts of all samples for analysis with edgeR)

The align_and_snp_call.R script calls the following programs (versions are mentioned in the script):
1. cutadapt
2. bwa
3. samtools
4. varscan

Once the per gene read counts are obtained, the edgeR_Rscript.R can be used. This script performs the basic edgeR analysis to give a list of differentially expressed genes, their fold changes and P values.

The remaining scripts perform more analyses, some require the edgeR_Rscript.R to be run first or it’s workspace opened first (this is mentioned in the script):
make_unmapped_fastas.py (makes fastas of unmapped reads - to check what these are)
analysis_of_operonic_and_non_operonic_genes.R (was used for generating Fig. S4)
correlations_of_raw_read_counts_Rscript.R (gets correlations between samples based on raw read counts Fig S3)
checking_edgeR_fold_changes_with_manually_calculated_fold_changes_Rscript.R (performs a crude manual caluculation of fold changes and compares to edgeR derived ones – sanity check)
create_topgo_annotation_file.R (creates annotation required for topGO to perform enrichment analyses)
topgo_analysis_Rscript.R (perform topGO enrichment analysis)
functions_Rscript.R (contains some functions I use in the scripts)
peter_etal_analysis_Rscript.R (analysis with Peter et. al. Dataset, used in figures S18 and figure 4)
tf_target_analysis_Rscript.R (analysis with TFs and their targets among the differentially expressed genes)

