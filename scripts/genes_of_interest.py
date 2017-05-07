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

