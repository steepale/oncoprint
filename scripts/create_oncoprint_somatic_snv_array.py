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

