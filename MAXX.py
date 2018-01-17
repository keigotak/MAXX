import argparse
import sys
import re
import copy
import itertools
import operator

parser = argparse.ArgumentParser()
parser.add_argument("-m", dest="mutation_file", required=True, help="path to mutation file, must be formatted as the example.")
parser.add_argument("-f", dest="fasta_file", required=True, help="path to reference genome, in fasta format")
parser.add_argument("-g", dest="gtf_file", required=True, help="path to GTF file that corresponds to refence genome")
parser.add_argument("-s", dest="sample_name", required=True, help="name of sample being used")
args = parser.parse_args()

#mutation list to generate a tumor specific genome
mutation_file = open(args.mutation_file,'r')
header_line = mutation_file.readline()

#gtf file genome
gtf_file = open(args.gtf_file,'r')

#reference genome
fasta_file = open(args.fasta_file,'r')

#the name of sample
sample_name = args.sample_name


#I'm going to go ahead and organize all of the mutations first.

#A dictionary of all mutations
mutation_dict = {}

for line in mutation_file:
	line = line.replace("\"","")
	line_list = line.strip().split()
	line_list[2] = int(line_list[2])
	line_list[3] = int(line_list[3])
	if line_list[0] not in mutation_dict:
		mutation_dict[line_list[0]] = [line_list]
	else:
		mutation_dict[line_list[0]].append(line_list)
	
for mutated_gene in mutation_dict:
	mutation_dict[mutated_gene] = sorted(mutation_dict[mutated_gene], key=lambda x: x[2])

#Getting the gene position of each gene with a mutation
mutated_gene_position_dict = {}
for line in gtf_file:
	line = line.replace("\"","")
	line_list = line.strip().split("\t")
	print line_list
	if "#" not in line and line_list[2] == "gene":
		attribute_list = line_list[8].split(";")
		gene = attribute_list[4].replace(" gene_name ","")
		if gene in mutation_dict:
			temp_list = []
			temp_list.append(int(line_list[3]))
			temp_list.append(int(line_list[4]))
			mutated_gene_position_dict[gene] = temp_list


#Check to see if all mutated genes were accounted for and remove the mutations that weren't present in the gtf file
if len(mutated_gene_position_dict) == 0:
	print "Weird, no genes from the mutation file matched the genes from the gtf file.  Maybe the order of the attributes within your gtf file is in a different order than the one I was using?"
	sys.exit()
else:
	all_genes_accounted_for = True
	mutated_gene_removal_list = []
	for mutated_gene in mutation_dict:
		if mutated_gene not in mutated_gene_position_dict:
			all_genes_accounted_for = False
			print mutated_gene + " was not found within the gtf file. My best guess is that the gene names aren't consitent between the gtf file and the mutation file."
			mutated_gene_removal_list.append(mutated_gene)

	if all_genes_accounted_for == True:
		print "We found the position for all mutated genes!"
	else:
		for mutated_gene in mutated_gene_removal_list:
			mutation_dict.pop(mutated_gene, None)

#Lets load in that Fasta File/Reference Genome! By the way, it's important that the key has "chr" in it (e.g. chr1).  This beacuae we're going to assoicate the mutation with the chromosome. 

chr_dict = {}
sequence = ""
chr = ""
for line in fasta_file:
	line = line.replace("\"","")
	line = line.strip()
	if "#" not in line:
		if ">" in line:
			if chr == "":
				chr = line.replace(">","")
			else:
				chr_dict[chr] = sequence
				chr = line.replace(">","")
				sequence = ""	
		else:
			sequence = sequence + line.upper()
chr_dict[chr] = sequence

#Time to get the wild type sequence of the genes that are mutated.  Also I want to make sure we are using the right fasta file, by checking that the reference nucleotides match between the fasta file and our mutation data. 
tumor_ref_genome = open(sample_name + ".fa",'w')

for mutated_gene in mutation_dict:
	tumor_ref_genome.write(">" + mutated_gene + "\n")
	current_chromosome = chr_dict["chr" + mutation_dict[mutated_gene][0][1]]
        start_pos = mutated_gene_position_dict[mutated_gene][0] - 201
        end_pos = mutated_gene_position_dict[mutated_gene][1] + 200
        tumor_ref_genome.write(current_chromosome[start_pos:end_pos] + "\n")
	for mutation in mutation_dict[mutated_gene]:
        	if mutation[4] != current_chromosome[mutation[2]-1:mutation[3]]:
                	print "That's weird, the reference nucleoted of " + mutation[0] + " at postion " + str(mutation[2]) + " dosen't match the reference genome you're using."

#Now I'm going to generate the tumor's mutated gene sequence. Because the mutation sequences may be altereted due to mutation sifting, the original mutation dict is going to be copied. 
altered_mutation_dict = copy.deepcopy(mutation_dict)

for mutated_gene in mutation_dict:
	current_chromosome = chr_dict["chr" + mutation_dict[mutated_gene][0][1]]
	start_pos = mutated_gene_position_dict[mutated_gene][0] - 201
        end_pos = mutated_gene_position_dict[mutated_gene][1] + 200
	gene_sequence = current_chromosome[start_pos:end_pos]
	tumor_ref_genome.write(">" + mutated_gene + "_Mut\n")
	if len(mutation_dict[mutated_gene]) == 1:
		mutation = mutation_dict[mutated_gene][0]
		mut_start_pos = mutation[2] - start_pos - 0
                mut_end_pos = mutation[3] - start_pos - 0
                mutant_sequence = gene_sequence[:(mut_start_pos-1)] + mutation[5] + gene_sequence[mut_end_pos:]
                tumor_ref_genome.write(mutant_sequence + "\n")
	else:
		shift = 0
		mutation_count = 0
		mutant_sequence = gene_sequence
		for mutation in mutation_dict[mutated_gene]:
                        mut_start_pos = mutation[2] - start_pos - 0 + shift
                        mut_end_pos = mutation[3] - start_pos - 0 + shift
                        mutant_seqeunce = mutant_sequence[:(mut_start_pos-1)] + mutation[5] + mutant_sequence[mut_end_pos:]
                        altered_mutation_dict[mutated_gene][mutation_count][2] = mut_start_pos
			altered_mutation_dict[mutated_gene][mutation_count][3] = mut_end_pos
			
			if len(mutation[5]) > 1:
                        	shift += 1 * (len(mutation[5]) - 1)
                        if len(mutation[4]) > 1:
                        	shift += -1 * (len(mutation[4]) - 1)
                        mutation_count += 1

                tumor_ref_genome.write(mutant_sequence + "\n")



#Generate tumor mutation index file.  
tumor_index_file = open(sample_name + "_Index.txt",'w')

for mutated_gene in mutation_dict:
	mutation_count = 0
	wt_header = ">" + mutated_gene
	mut_header = ">" + mutated_gene + "_Mut"
	start_pos = mutated_gene_position_dict[mutated_gene][0]
	for mutation in mutation_dict[mutated_gene]:
		wt_start = mutation[2] - start_pos + 201
		wt_end = mutation[3] - start_pos + 201
		if len(mutation[4]) > 1:
			wt_start += 1
		tumor_index_file.write(wt_header + "\t" + str(wt_start) + "\t" + str(wt_end) + "\t" + mutation[4])
		for x in range(6, len(mutation)):
			tumor_index_file.write("\t" + mutation[x])
		tumor_index_file.write("\n")
		mutation = altered_mutation_dict[mutated_gene][mutation_count]
		mut_start = mutation[2]
                mut_end = mutation[3]
		if len(mutation[5]) > 1:
                        mut_start = mut_start + 1
                        mut_end = mut_end + len(mutation[5])-1
                if len(mutation[4]) > 1:
                	mut_end = mut_end-len(mutation[4])+1
		
		tumor_index_file.write(mut_header + "\t" + str(mut_start) + "\t" + str(mut_end) + "\t" + mutation[5])
		for x in range(6, len(mutation)):
                        tumor_index_file.write("\t" + mutation[x])
                tumor_index_file.write("\n")	
		mutation_count += 1
