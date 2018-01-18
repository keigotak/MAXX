# MAXX
The MAXX script generates a tumor specific reference genome and an accompaning mutation index file from a list of muations, reference geneome, and gtf file.  Tumor specific reference genomes have been demonstatrated to dramatically enhance alignment of RNA-sequencing reads containing indel mutations, providing a comprehensive analysis of tumor allelic expression.  Since sequencing methods and mutation detection softwares are dramatically improving, the ability to detect indel mutations within RNA-sequencing data will continue to be an important aspect in performing mutation allelic expression analyses.
## Setup and Processing
### Requirements
Python: any version
### Running MAXX
MAXX is a python script, thus it is ran using the command line.  To run MAXX, first download the MAXX.py script.  Next open the terminal (MAC) or the command prompt (PC) and change your working directory to where the MAXX.py script is located. Finally, enter "python MAXX.py -m Mutation_File_Path -f Reference_Genome_Path -g GTF_File_Path -s Sample_Name" into the terminal/command prompt and hit enter.  MAXX will take approximatly 2-3 minutes to run and output a tumor reference genome and an accompaning mutation index file.
## Required Input Paramters
### -m
The path to the mutation list that are going to be used to generate the tumor specific reference genome.  This file must be tab delimted and follow the format demonstrated below.
```
Gene    Chr     Start   End     Ref_allele      Alt_allele
Gene1   1       222     222     C       T
Gene2   1       212     213     CA      C
Gene3   1       215     215     T       TG
Gene4   1       205     209     TGGAA   T
Gene4   1       201     201     A       AT
```
However, there may be additional columns futher describing the mutation.
### -f 
The path to the reference genome where the nuceleotides will be altered to match the tumor's mutated sequence.  Note, MAXX will identify if the reference nucelotied presented in the mutation list does not match the reference nucelotide in the reference geneome.  As for the format of the header line for each chromsome, it has to contain only >chr + chromsome.  For example:
```
>chr1
```
### -g
The path to the GTF file that corresponds with the reference geneome. MAXX uses the GTF file to identify the start and end positions for each mutated gene.  In order for MAXX to run properly, the GTF attribute identifies must be formatted such as this.  
```
chr1    HAVANA  gene    201     225     .       +       .       gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "Gene1"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "Gene1"; level 2; havana_gene "OTTHUMG00000000961.2";
```
### -s
This is simply the name of the sample that contains the mutations within the mutation list.  The -s paramter is used to label the output files e.g. SAMPLE.fa and SAMPLE_Index.txt

## Output Files
### SAMPLE.fa 
The tumor specific reference geneome, containing the wild type and mutant sequence for all mutated genes presented in the mutation list.  The reference geneome is used with to align the RNA-sequencing data. 
### SAMPLE_Index.txt 
A file containing the name of the wild type and mutant alleles along with the positions of the reference nucelotide and mutant nucelodie.  The index file is required by softwares such as Bam Read Count and Sam Tools to identify how many reads aligned to that particular nucleotide.

## Special Features
1. Handles genes that have more than one mutation.
2. Extends gene sequences by additional 200 nucelotides on each end.
3. Identifies refernce nuceltoide mismatches between gtf reference geneome and
4. Condensed example input and output files.
5. relatively small size while mainting effectivness.

##FAQs



