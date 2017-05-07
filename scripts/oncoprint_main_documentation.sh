#===============================================================================
#
#         FILE: home/proj/MDW_genomics/steepale/oncoprint/scripts/oncoprint_main_documentation.txt
#
#        USAGE: for documentation purposes, scripts inside
#
#  DESCRIPTION:  This script serves as a step by step documentation script and development script for performing oncoprint array preparations for GenVisR in R
#                
# REQUIREMENTS:  ---
#        NOTES:  ---
#       AUTHOR:  Alec Steep, steepale@msu.edu
#  AFFILIATION:  Michigan State University (MSU), East Lansing, MI, United States
#				         USDA ARS Avian Disease and Oncology Lab (ADOL), East Lansing, MI, United States
#				         Technical University of Munich (TUM), Weihenstephan, Germany
#      VERSION:  1.0
#      CREATED:  2017.02.17
#     REVISION:  
#===============================================================================

# Permanent PROJECT DIRECTORY (TUM Cluster)
cd /home/proj/MDW_genomics/steepale/oncoprint

# Files with annotated variants
./data/somatic_snvs_final.txt
./data/somatic_indels_final.txt
# Annotated CNVs and LOH
./data/scnas_annotation.bed
# Chicken-human orthologs from ensembl
/home/proj/MDW_genomics/steepale/oncoprint/data/biomart_chicken-human_orthologs_2017_01_04.txt
# Annotated somatic gene fusions
/home/proj/MDW_genomics/steepale/gene_fusions/data/chimeric_junctions_filtered_annotated.txt

# Create reference files for array script
# Genes mutated via non-synonymous SNVs
grep -v "^#" /home/proj/MDW_genomics/steepale/pathway_analysis/results/somatic_snvs_final.txt | \
cut -f8 | sort | uniq > /home/users/a.steep/databases/samples/somatic_nonsyn_snv_mutated_genes.txt
# Genes mutated via non-synonymous INDELs
grep -v "^#" /home/proj/MDW_genomics/steepale/pathway_analysis/results/somatic_indels_final.txt | \
cut -f8 | sort | uniq > /home/users/a.steep/databases/samples/somatic_nonsyn_indel_mutated_genes.txt
# Genes mutated via gene fusions
(grep -v "^#" /home/proj/MDW_genomics/steepale/gene_fusions/data/chimeric_junctions_filtered_annotated.txt | cut -f2; \
grep -v "^#" /home/proj/MDW_genomics/steepale/gene_fusions/data/chimeric_junctions_filtered_annotated.txt | cut -f3) | \
sort | uniq > /home/users/a.steep/databases/samples/somatic_gene_fusion_genes_ensembl_gene_name.txt










# Create an array in python
python ./scripts/create_oncoprint_somatic_snv_array.py

# ./scripts/create_oncoprint_somatic_snv_array.py
###################################
import sys
import os
import numpy

# Reference files
tumor_birds_file = "/home/users/a.steep/databases/samples/tumor_sample_dnaseq_list_NNN-N_SN.txt"
germline_birds_file = "/home/users/a.steep/databases/samples/germline_sample_dnaseq_list_NNN-N_SN.txt"
snv_file = open("/home/proj/MDW_genomics/steepale/pathway_analysis/results/somatic_snvs_final.txt")
indel_file = open("/home/proj/MDW_genomics/steepale/pathway_analysis/results/somatic_indels_final.txt")
snv_genes_file = open("/home/users/a.steep/databases/samples/somatic_nonsyn_snv_mutated_genes.txt")
indel_genes_file = open("/home/users/a.steep/databases/samples/somatic_nonsyn_indel_mutated_genes.txt")
mut_genes_rename_file = open("./data/somatic_nonsyn_snv_indel_mutated_genes_renamed.txt", 'w')
snv_array_file = "./data/oncoprint_somatic_snv_array_geneid-x_tumors-y.txt"
indel_array_file = "./data/oncoprint_somatic_indel_array_geneid-x_tumors-y.txt"
snv_indel_file = open("/home/proj/MDW_genomics/steepale/pathway_analysis/results/somatic_snvs_and_indels_final.txt")

# Create a unique list of all tumor samples
tumors = set()
for tumor_bird in open(tumor_birds_file):
	tumor_bird = tumor_bird.rstrip()
	tumors.add(tumor_bird)
tumors = list(tumors)
tumors.sort()
#print(tumors)

# Create and fill null dictionary for gene_id to gene_symbol relationship
gene_id2gene_symbol = {}
# Create and fill null dictionary to say if gene in sample is mutated via SNV
snv_status = {}
for snv_gene in snv_genes_file:
	snv_gene = snv_gene.rstrip()
	snv_status[snv_gene] = []

# Create and fill null dictionary to say if gene in sample is mutated via INDEL
indel_status = {}
for indel_gene in indel_genes_file:
	indel_gene = indel_gene.rstrip()
	indel_status[indel_gene] = []

# Create a unique list of all gene_ids
gene_ids = set()
for snv_indel_line in snv_indel_file:
	snv_indel_line = snv_indel_line.rstrip()
	if snv_indel_line[0] != '#':
		snv_indel_col = snv_indel_line.split('\t')
		snv_indel_gene_id = snv_indel_col[7]
		snv_indel_gene_symbol = snv_indel_col[6]
		snv_indel_tsn = int(snv_indel_col[8])
		snv_indel_sample = snv_indel_col[9]
		gene_ids.add(snv_indel_gene_id)
		# Add to the dictionary of gene_id2gene_symbol
		if snv_indel_gene_symbol != 'NA':
			gene_id2gene_symbol[snv_indel_gene_id] = snv_indel_gene_symbol
		elif snv_indel_gene_symbol == 'NA':
			gene_id2gene_symbol[snv_indel_gene_id] = snv_indel_gene_id
gene_ids = list(gene_ids)
gene_ids.sort()
#print(gene_ids)

# Fill in dictionary to say if gene in sample is mutated via SNV
for snv_line in snv_file:
	snv_line = snv_line.rstrip()
	if snv_line[0] != '#':
		snv_col = snv_line.split('\t')
		snv_gene_id = snv_col[7]
		snv_gene_symbol = snv_col[6]
		snv_tsn = int(snv_col[8])
		snv_sample = snv_col[9]
		# Add to the dictionary of SNV mutated genes
		if snv_tsn == 1:
			snv_status[snv_gene_id].append(snv_sample)
		elif snv_tsn > 1:
			sample = []
			for n in range(snv_tsn):
				sample.append(n)
				sample[n] = snv_sample.split('|')[n]
				snv_status[snv_gene_id].append(sample[n])

# Fill in dictionary to say if gene in sample is mutated via INDEL
for indel_line in indel_file:
	indel_line = indel_line.rstrip()
	if indel_line[0] != '#':
		indel_col = indel_line.split('\t')
		indel_gene_id = indel_col[7]
		indel_gene_symbol = indel_col[6]
		indel_tsn = int(indel_col[8])
		indel_sample = indel_col[9]
		# Add to the dictionary of indel mutated genes
		if indel_tsn == 1:
			indel_status[indel_gene_id].append(indel_sample)
		elif indel_tsn > 1:
			sample = []
			for n in range(indel_tsn):
				sample.append(n)
				sample[n] = indel_sample.split('|')[n]
				indel_status[indel_gene_id].append(sample[n])

# Test if dictionary is accurate
#for gene_id, tumors in snv_status.items():
#	print('gene_id: ' + gene_id + '\t' + 'tumors: ' + str(tumors))

for gene_id in gene_ids:
	mut_genes_rename_file.write(gene_id2gene_symbol[gene_id] + '\n')

# Create and fill SNV mutation array
onco_snvs = numpy.zeros((len(gene_ids), len(tumors)))
for (gene_id_mp, gene_id) in enumerate(gene_ids):
	for (tumor_mp, tumor) in enumerate(tumors):
		for mut_gene_id, mut_tumors_list in snv_status.items():
			for mut_tumor in mut_tumors_list:
				if mut_gene_id == gene_id and mut_tumor == tumor:
					onco_snvs[gene_id_mp, tumor_mp] = 1

array_header = "738-1_S1"+'\t'+"741-1_S2"+'\t'+"756-3_S3"+'\t'+"766-1_S4"+'\t'+"777-3_S14"+'\t'+"787-2_S15"+'\t'+"788-1_S16"+'\t'+"794-1_S17"+'\t'+"798-1_S5"+'\t'+"833-1_S6"+'\t'+"834-2_2_S12"+'\t'+"834-2_S7"+'\t'+"835-1_S18"+'\t'+"841-3_S19"+'\t'+"842-2_2_S25"+'\t'+"842-2_S20"+'\t'+"855-1_S8"+'\t'+"863-1_S9"+'\t'+"884-2_S21"+'\t'+"901-2_2_S26"+'\t'+"901-2_S22"+'\t'+"906-1_S23"+'\t'+"911-1_2_S13"+'\t'+"911-1_S24"+'\t'+"918-3_S10"+'\t'+"927-2_S11"

numpy.savetxt(snv_array_file, onco_snvs, fmt='%i', delimiter='\t', newline='\n', header=array_header)

os.system('sed -i "s/^# //" ./data/oncoprint_somatic_snv_array_geneid-x_tumors-y.txt')

# Create and fill indel mutation array
onco_indels = numpy.zeros((len(gene_ids), len(tumors)))
for (gene_id_mp, gene_id) in enumerate(gene_ids):
	for (tumor_mp, tumor) in enumerate(tumors):
		for mut_gene_id, mut_tumors_list in indel_status.items():
			for mut_tumor in mut_tumors_list:
				if mut_gene_id == gene_id and mut_tumor == tumor:
					onco_indels[gene_id_mp, tumor_mp] = 1

array_header = "738-1_S1"+'\t'+"741-1_S2"+'\t'+"756-3_S3"+'\t'+"766-1_S4"+'\t'+"777-3_S14"+'\t'+"787-2_S15"+'\t'+"788-1_S16"+'\t'+"794-1_S17"+'\t'+"798-1_S5"+'\t'+"833-1_S6"+'\t'+"834-2_2_S12"+'\t'+"834-2_S7"+'\t'+"835-1_S18"+'\t'+"841-3_S19"+'\t'+"842-2_2_S25"+'\t'+"842-2_S20"+'\t'+"855-1_S8"+'\t'+"863-1_S9"+'\t'+"884-2_S21"+'\t'+"901-2_2_S26"+'\t'+"901-2_S22"+'\t'+"906-1_S23"+'\t'+"911-1_2_S13"+'\t'+"911-1_S24"+'\t'+"918-3_S10"+'\t'+"927-2_S11"

numpy.savetxt(indel_array_file, onco_indels, fmt='%i', delimiter='\t', newline='\n', header=array_header)

os.system('sed -i "s/^# //" ./data/oncoprint_somatic_indel_array_geneid-x_tumors-y.txt')

# Test for IKZF1
for (gene_id_mp, gene_id) in enumerate(gene_ids):
	for (tumor_mp, tumor) in enumerate(tumors):
		if gene_id == 'ENSGALG00000013086':
			print(onco_snvs[gene_id_mp, ])
			print(onco_indels[gene_id_mp, ])
###################################

# Create a reference file of the sample titles for gene fusions
grep -v "^#" /home/proj/MDW_genomics/steepale/gene_fusions/data/chimeric_junctions_filtered_annotated.txt | cut -f5 | sort | uniq | sed "s/;/\n/g" | sort | uniq >/home/users/a.steep/databases/samples/tumor_sample_rnaseq_list_gene_fusions_017NNN-N_redundent_1s.txt


# Create an array in python
python ./scripts/create_oncoprint_somatic_gene_fusion_array.py

# ./scripts/create_oncoprint_somatic_gene_fusion_array.py
###################################
import sys
import os
import numpy

# Reference files
tumor_birds_file = "/home/users/a.steep/databases/samples/tumor_sample_rnaseq_list_gene_fusions_017NNN-N_redundent_1s.txt"
#snv_file = open("/home/proj/MDW_genomics/steepale/pathway_analysis/results/somatic_snvs_final.txt")
#indel_file = open("/home/proj/MDW_genomics/steepale/pathway_analysis/results/somatic_indels_final.txt")
fusion_file = open("/home/proj/MDW_genomics/steepale/gene_fusions/data/chimeric_junctions_filtered_annotated.txt")
#snv_genes_file = open("/home/users/a.steep/databases/samples/somatic_nonsyn_snv_mutated_genes.txt")
#indel_genes_file = open("/home/users/a.steep/databases/samples/somatic_nonsyn_indel_mutated_genes.txt")
fusion_genes_file = open("/home/users/a.steep/databases/samples/somatic_gene_fusion_genes_ensembl_gene_name.txt")
#snv_array_file = "./data/oncoprint_somatic_snv_array_geneid-x_tumors-y.txt"
#indel_array_file = "./data/oncoprint_somatic_indel_array_geneid-x_tumors-y.txt"
fusion_array_file = "./data/oncoprint_somatic_gene_fusions_array_geneid-x_tumors-y.txt"

# Create a unique list of all tumor samples
tumors = set()
for tumor_bird in open(tumor_birds_file):
	tumor_bird = tumor_bird.rstrip()
	tumors.add(tumor_bird)
tumors = list(tumors)
tumors.sort()
#print(tumors)

# Create and fill null dictionary to say if gene in sample is mutated via gene fusion
# &
# Create a unique list of all fusion gene names
fusion_status = {}
fusion_genes = set()
for fusion_gene in fusion_genes_file:
	fusion_gene = fusion_gene.rstrip()
	fusion_genes.add(fusion_gene)
	fusion_status[fusion_gene] = []
fusion_genes = list(fusion_genes)
fusion_genes.sort()
#print(len(fusion_genes))

# Fill in dictionary to say if gene in sample is mutated via fusion
for fusion_line in fusion_file:
	fusion_line = fusion_line.rstrip()
	if fusion_line[0] != '#':
		fusion_col = fusion_line.split('\t')
		fusion_gene_nameA = fusion_col[1]
		fusion_gene_nameB = fusion_col[2]
		fusion_tsn = int(fusion_col[3])
		fusion_sample = fusion_col[4]
		# Add mutated samples to the dictionary of fusion mutated genes
		if fusion_tsn == 1:
			fusion_status[fusion_gene_nameA].append(fusion_sample)
			fusion_status[fusion_gene_nameB].append(fusion_sample)
		elif fusion_tsn > 1:
			sample = []
			for n in range(fusion_tsn):
				sample.append(n)
				sample[n] = fusion_sample.split(';')[n]
				fusion_status[fusion_gene_nameA].append(sample[n])
				fusion_status[fusion_gene_nameB].append(sample[n])

# Test if dictionary is accurate
#for gene, tumors in fusion_status.items():
#	print('gene: ' + gene + '\t' + 'tumors: ' + str(tumors))

# Create and fill fusion mutation array
onco_fusions = numpy.zeros((len(fusion_genes), len(tumors)))
for (gene_name_mp, gene_name) in enumerate(fusion_genes):
	for (tumor_mp, tumor) in enumerate(tumors):
		for mut_gene_name, mut_tumors_list in fusion_status.items():
			for mut_tumor in mut_tumors_list:
				if str(mut_gene_name) == str(gene_name) and str(mut_tumor) == str(tumor):
					onco_fusions[gene_name_mp, tumor_mp] = 1

array_header = "738-1_S"+'\t'+"741-1_S"+'\t'+"766-1_S"+'\t'+"777-3_S"+"794-1_S"+'\t'+"798-1_S"+'\t'+"798-1_2_S"+'\t'+"833-1_S"+'\t'+"835-1_S"+'\t'+"841-3_S"+'\t'+"842-2_1_S"+'\t'+"842-2_S"+'\t'+"855-1_S"+'\t'+"863-1_S"+'\t'+"884-2_S"+'\t'+"901-2_2_S"+'\t'+"906-1_S"+'\t'+"911-1_S"+'\t'+"911-1_2_S"+'\t'+"918-3_S"+'\t'+"927-2_S"

numpy.savetxt(fusion_array_file, onco_fusions, fmt='%i', delimiter='\t', newline='\n', header=array_header)

os.system('sed -i "s/^# //" ./data/oncoprint_somatic_gene_fusions_array_geneid-x_tumors-y.txt')

# Test for IKZF1
for (gene_id_mp, gene_id) in enumerate(fusion_genes):
	for (tumor_mp, tumor) in enumerate(tumors):
		if gene_id == 'TSG101':
			print(onco_fusions[gene_id_mp, ])


###################################
























# Generate a gene of interest file
# This file is generated based on frequency of mutations per gene

python ./scripts/genes_of_interest.py \
/home/proj/MDW_genomics/steepale/pathway_analysis/results/somatic_snvs_and_indels_final.txt \
./data/genes_of_interest_somatic_snvs_and_indels_final.int

# ./scripts/genes_of_interest.py
###################################
import sys
import os

# in- and out-file
infile = sys.argv[1]
outfile = open(sys.argv[2], 'w')

# Reference files
# Unique list of genes (either ensemble ID or gene symbol)
mut_genes_rename_file = "./data/somatic_nonsyn_snv_indel_mutated_genes_renamed.txt"

# Header for outfile
outfile.write('#GENE' + '\t' + 'TOTAL_MUTS_ACROSS_SAMPLES' + '\n')

# Read each line of the input file (collective filtered somatic SNVs and INDELs) and mutations per gene
for mut_gene in open(mut_genes_rename_file):
	mut_gene = mut_gene.rstrip()
	interest_genes_tsn = 0
	for in_line in open(infile):
		if in_line[0] != '#':
			in_line = in_line.rstrip()
			in_col = in_line.split('\t')
			in_chr = in_col[0]
			in_pos = in_col[1]
			in_symbol = in_col[6]
			in_geneid = in_col[7]
			in_tsn = int(in_col[8])
			if mut_gene == in_symbol or mut_gene == in_geneid:
				interest_genes_tsn = interest_genes_tsn + in_tsn
	outfile.write(mut_gene + '\t' + str(interest_genes_tsn) + '\n')
# Close the outfile
outfile.close()

###################################

# Sort the file and determine a threshold for which genes are interesting
(grep "^#" ./data/genes_of_interest_somatic_snvs_and_indels_final.int; \
grep -v "^#" ./data/genes_of_interest_somatic_snvs_and_indels_final.int | \
sort -k2,2nr) > ./data/genes_of_interest_somatic_snvs_and_indels_final.txt

# Send the files to macbook pro for oncoprint generation in R Studio
# From (/Users/Alec/Documents/Bioinformatics/MDV_Project/oncoprint)
# Transfer the SNV mutation status array
rsync -avp \
a.steep@barcelona.binfo.wzw.tum.de:/home/proj/MDW_genomics/steepale/oncoprint/data/oncoprint_somatic_snv_array_geneid-x_tumors-y.txt \
./data/
# Transfer the INDEL mutation status array
rsync -avp \
a.steep@barcelona.binfo.wzw.tum.de:/home/proj/MDW_genomics/steepale/oncoprint/data/oncoprint_somatic_indel_array_geneid-x_tumors-y.txt \
./data/
# Transfer mutated gene labels (x-axis for R matrices)
rsync -avp \
a.steep@barcelona.binfo.wzw.tum.de:/home/proj/MDW_genomics/steepale/oncoprint/data/somatic_nonsyn_snv_indel_mutated_genes_renamed.txt \
./data/
# Transfer the genes of interest file
rsync -avp \
a.steep@barcelona.binfo.wzw.tum.de:/home/proj/MDW_genomics/steepale/oncoprint/data/genes_of_interest_somatic_snvs_and_indels_final.txt \
./data/






