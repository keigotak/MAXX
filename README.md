# MAXX
The MAXX script generates a tumor specific reference genome from a list of muations, reference geneome, and a gtf file.  Tumor specific reference genomes have demonstatrated to dramatically enhance RNA-sequencing reads to indel mutations.  With the increasing ability to detect indel mutations, tumor specific reference geneomes generated via MAXX provide the most comprehensive method of tumor allelic expression.
## Setup
### Requirements
Python: any version
### Getting Started
MAXX is a python script, thus it is ran using the command line.  To run MAXX, first download the MAXX.py script.  Next open the terminal (MAC) or the command prompt (PC) and change your working directory to where the MAXX.py script is located. Finally, enter "python MAXX.py -m Mutation_File -f Reference_Genome -g GTF_File -s Sample_Name" into the terminal/command prompt and hit enter.
## Input Paramters
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



