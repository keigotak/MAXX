# MAXX
The MAXX script generates a tumor specific reference genome and an accompaning mutation index file from a list of muations, reference geneome, and gtf file.  Tumor specific reference genomes have been demonstatrated to dramatically enhance alignment of RNA-sequencing reads containing indel mutations, providing a comprehensive analysis of mutation allelic expression.  Since sequencing methods and mutation detection softwares are dramatically improving, the ability to detect indel mutations within RNA-sequencing data will continue to be an important aspect in performing mutation allelic expression analyses.
## Setup and Processing
### Requirements
Python: any version
### Running MAXX
MAXX is a python script, thus it is ran using the command line.  To run MAXX, first download the MAXX.py script.  Next open the terminal (MAC) or the command prompt (PC) and change your working directory to where the MAXX.py script is located. Once in the proper directory, enter "python MAXX.py -m Mutation_File_Path -f Reference_Genome_Path -g GTF_File_Path -s Sample_Name" into the terminal/command prompt and hit enter.  MAXX will take approximatly 2-3 minutes to run and output a tumor specific reference genome and an accompaning mutation index file.
## Required Input Paramters
### -m
The path to the mutation list file. This file cotains a list of mutations that will be used to generate the tumor specific reference genome.  For MAXX to run without an error, the mutation file must be tab delimted and follow the format demonstrated below.  It is also crucial that the Ref_allele and Alt_allele columns follow the vcf format.  Where insertions and deletions contain the previous reference nucleotide.  
```
Gene    Chr     Start   End     Ref_allele      Alt_allele
Gene1   1       222     222     C       T
Gene2   1       212     213     CA      C
Gene3   1       215     215     T       TG
Gene4   1       205     209     TGGAA   T
Gene4   1       201     201     A       AT
```
In addition to the six columns, there may be additional columns futher describing the mutation.  MAXX will ignore any additional columns, however they will be carried over to the mutation index file. 
### -f 
The path to the reference genome file, most likely will be Hg19 or GRCh38.  MAXX uses the input reference genome to obtain the genetic sequence of genes presented in the mutation list.  The obtained gene sequence is then altered to match the gene sequence exhibited by the tumor. For MAXX to run properly, the header line for each chromsome within the reference file must be labled as >chr + chromsome.  For example:
```
>chr1
```
### -g
The path to the GTF file that corresponds to the input reference geneome. MAXX uses the GTF file to identify the start and end positions for each mutated gene, these positions are then used to obtain the genetic sequence from the reference genome.  In order for MAXX to run properly, the GTF attribute "gene_name" must be positoned 5th in the attribue list. For example:  
```
chr1    HAVANA  gene    201     225     .       +       .       gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "Gene1"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "Gene1"; level 2; havana_gene "OTTHUMG00000000961.2";
``` 
### -s
This is simply the name of the sample that contains the mutations within the mutation list.  The -s paramter is used to label the output files, e.g. SAMPLE.fa and SAMPLE_Index.txt
## Output Files
### SAMPLE.fa 
A tumor specific reference geneome in the fasta format.  The new reference genome contains the wild type and mutant sequence for all mutated genes present in the mutation list.  This reference geneome is ideally used in conjuction with a RNA-sequencing aligner to increase alignment of mutant RNA-sequencing reads.  
### SAMPLE_Index.txt 
A file containing the new location of all input mutations within the tumor specific reference genome.  This index file is required by softwares such as Bam Read Count and Sam Tools to identify how many reads aligned to that particular nucleotide.
## Special Features of MAXX
1. Can make mutant gene sequences from genes that contain more than one mutation by using a nucleotide shifting algorithim.
2. Extends gene sequences placed within the tumor specific reference genome by 200 nucleotides on each end to provide proper RNA-sequencing aligment.
3. Produces a warning when the reference allele nucleotide within the mutation list doesn't correspond within the nucleotide within the reference genome. This ensures that the correct reference genomes and GTF files are being used. 
4. MAXX generated reference genomes are approximatly 70 times smaller than a typical Hg19 reference genome.  They also require much less comptaional power to run a RNA-sequencing aligner and maintain accuracy of RNA-sequencing alignment for the selected genes. 
5. Example input and output files are provided for testing and tweaking the MAXX.py script.  
## FAQs

