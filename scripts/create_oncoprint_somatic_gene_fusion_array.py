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


