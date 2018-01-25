# MAXX
The MAXX.py script generates a tumor specific reference genome and an accompanying mutation index file from a list of mutations, reference genome, and gtf file.  Tumor specific reference genomes have been demonstrated to dramatically enhance alignment of RNA-sequencing reads containing indel mutations, providing an unbiased analysis of mutation allelic expression.  Since sequencing methods and mutation detection software are dramatically improving, the ability to detect indel mutations within RNA-sequencing data will continue to be an important aspect in performing mutation allelic expression analyses.
## Setup and Processing
### Requirements
Python: any version
### Running MAXX
MAXX is a python script; thus, it is ran using the command line.  To run MAXX, first download the MAXX.py script.  Next open the terminal (MAC) or the command prompt (PC) and change your working directory to where the MAXX.py script is located. Once in the proper directory, enter "python MAXX.py -m Mutation_File_Path -f Reference_Genome_Path -g GTF_File_Path -s Sample_Name" into the terminal/command prompt and hit enter.  MAXX will take approximately 2-3 minutes to run and output a tumor specific reference genome along with an accompanying mutation index file.
## Required Input Parameters
### -m Mutation_File_Path
The path to the mutation list file. This file contains a list of mutations that will be used to generate the tumor specific reference genome.  For MAXX to run without an error, the mutation file must be tab delimited and follow the format demonstrated below.  It is also crucial that the Ref_allele and Alt_allele columns follow the vcf format.  Where insertions and deletions contain the previous reference nucleotide.  
```
Gene    Chr     Start   End     Ref_allele      Alt_allele
Gene1   1       222     222     C       T
Gene2   1       212     213     CA      C
Gene3   1       215     215     T       TG
Gene4   1       205     209     TGGAA   T
Gene4   1       201     201     A       AT
```
In addition to the six columns, there may be additional columns further describing the mutation.  MAXX will ignore any additional columns; however, they will be carried over to the mutation index file. 
### -f Reference_Genome_Path
The path to the reference genome file, most likely will be Hg19 or GRCh38.  MAXX uses the input reference genome to obtain the genetic sequence of genes presented in the mutation list.  The obtained gene sequence is then altered to match the gene sequence exhibited by the tumor. For MAXX to run properly, the header line for each chromosome within the reference file must be labeled as >chr + chromosome.  For example:
```
>chr1
```
### -g GTF_File_Path
The path to the GTF file that corresponds to the input reference genome. MAXX uses the GTF file to identify the start and end positions for each mutated gene, these positions are then used to obtain the genetic sequence from the reference genome.  In order for MAXX to run properly, the GTF attribute "gene_name" must be positioned 5th in the attribute list. For example:  
```
chr1    HAVANA  gene    201     225     .       +       .       gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "Gene1"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "Gene1"; level 2; havana_gene "OTTHUMG00000000961.2";
``` 
### -s Sample_Name
This is simply the name of the sample in which the tumor specific reference genome is being created for.  The -s parameter is used to label the output files, e.g. SAMPLE.fa and SAMPLE_Index.txt
## Output Files
### SAMPLE.fa 
A tumor specific reference genome in the fasta format.  The new reference genome contains the wild type and mutant sequence for all mutated genes present in the mutation list.  This reference genome is ideally used in conjunction with a RNA-sequencing aligner to increase alignment accuracy of mutant RNA-sequencing reads.  
### SAMPLE_Index.txt 
A file containing the new positions of all input mutations, relative to the tumor specific reference genome.  This index file is required by softwares such as Bam Read Count and Sam Tools to identify how many reads aligned to that particular nucleotide, providing the information to calculate the RNA allele frequency.
## Special Features of MAXX
1. Can create mutant gene sequences for genes that contain more than one mutation by using a nucleotide shifting algorithm.
2. Extends gene sequences placed into the tumor specific reference genome by 200 nucleotides on each end to provide proper RNA-sequencing alignment.
3. Produces a warning when the reference allele nucleotide within the mutation list doesn't match the nucleotide within the reference genome. This ensures that the correct reference genomes and GTF files are being used. 
4. MAXX generated reference genomes are approximately 70 times smaller than a typical Hg19 reference genome.  They also require much less computational power to run a RNA-sequencing aligner.  Despite this decrease in computational storage and power, tumor specific reference genomes maintain a high RNA-sequencing alignment for mutated genes. 
5. Example input and output files are provided for testing and altering the MAXX.py script for indvidual needs.  
## FAQs
