# Script based on: http://bioconductor.org/packages/devel/bioc/vignettes/ComplexHeatmap/inst/doc/s8.oncoprint.html

setwd("/Users/Alec/Documents/Bioinformatics/MDV_Project/oncoprint")
#source("https://bioconductor.org/biocLite.R")
#biocLite("ComplexHeatmap")
library(ComplexHeatmap)

# Load the numpy generated SNV mutation array via a tab seperate file (contains header) into R
mat_snvs = read.table("./data/oncoprint_somatic_snv_array_geneid-x_tumors-y.txt", header = TRUE, stringsAsFactors=FALSE, sep = "\t")
somatic_mut_gene_labels = read.table("./data/somatic_nonsyn_snv_indel_mutated_genes_renamed.txt")
# Make the row names unique with make.names
rownames(mat_snvs) = make.names(somatic_mut_gene_labels$V1, unique=TRUE)
mat_snvs = as.matrix(mat_snvs)
mat_snvs[mat_snvs=="1"]<-"SNV"
mat_snvs[mat_snvs=="0"]<-""

# Load the numpy generated INDEL mutation array via a tab seperate file (contains header) into R
mat_indels = read.table("./data/oncoprint_somatic_indel_array_geneid-x_tumors-y.txt", header = TRUE, stringsAsFactors=FALSE, sep = "\t")
somatic_mut_gene_labels = read.table("./data/somatic_nonsyn_snv_indel_mutated_genes_renamed.txt")
# Make the row names unique with make.names
rownames(mat_indels) = make.names(somatic_mut_gene_labels$V1, unique=TRUE)
mat_indels = as.matrix(mat_indels)
mat_indels[mat_indels=="1"]<-"INDEL"
mat_indels[mat_indels=="0"]<-""

# Combine the SNV and INDEL Matrices
mat_final = matrix(paste(mat_snvs,mat_indels,sep=";"),nrow=nrow(mat_snvs),ncol=ncol(mat_snvs))
dimnames(mat_final) <- list(rownames(mat_indels), colnames(mat_indels))
mat_final <- replace(mat_final, mat_final == ";", "")
mat_final <- replace(mat_final, mat_final == "SNV;", "SNV")
mat_final <- replace(mat_final, mat_final == ";INDEL", "INDEL")

# Load the genes of interest file for all genes
genes_interest_counts = read.table("./data/genes_of_interest_somatic_snvs_and_indels_final.txt", header = FALSE, stringsAsFactors=FALSE, sep = "\t")

# Load gene of interest for specific pathways
# Calcium regulation pathway
#genes_interest_counts = read.table("./data/ca_regulation_somatic_snvs_indels.txt", header = FALSE, stringsAsFactors=FALSE, sep = "\t")
# Fanconi Anemia pathway
#genes_interest_counts = read.table("./data/fanconi_anemia_somatic_snvs_indels.txt", header = FALSE, stringsAsFactors=FALSE, sep = "\t")

colnames(genes_interest_counts) = c("GENE", "TOTAL_COHORT_MUTS")
# Genes with at least one mutation
genes_of_interest_int = genes_interest_counts[genes_interest_counts$TOTAL_COHORT_MUTS > 0, ]
# Alter the genes of interest length
#genes_of_interest = genes_of_interest[1:50]
# OR Selecti genes with gene name annotation
#genes_of_interest_int = genes_interest_counts[!grepl("ENSGALG", genes_interest_counts$GENE), ]
genes_of_interest = genes_of_interest_int[1:30, 'GENE']

# Final Matrix
mat_oncoprint = mat_final[genes_of_interest, ]

# Set the parameters for contructing the oncoprint figure
alter_fun = list(background = function(x, y, w, h, v) 
{grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
           gp = gpar(fill = "#CCCCCC", col = NA))}, 
SNV = function(x, y, w, h) 
{grid.rect(x, y, w*0.33, h*1.00, 
           gp = gpar(fill = "blue", col = NA))},
INDEL = function(x, y, w, h) 
{grid.rect(x, y, w*0.66, h*0.90, 
           gp = gpar(fill = "#008000", col = NA))})
#AMP = function(x, y, w, h) {
#        grid.rect(x, y, w-unit(0.5, "mm"), h*0.50, gp = gpar(fill = "orange", col = NA))},
#DEL = function(x, y, w, h) {
#        grid.rect(x, y, w-unit(0.5, "mm"), h*0.50, gp = gpar(fill = "red", col = NA))})

col = c("SNV" = "blue", "INDEL" = "#008000")

ht = oncoPrint(mat_oncoprint,
               get_type = function(x) strsplit(x, ";")[[1]],
               alter_fun = alter_fun, 
               col = col,
               show_pct = TRUE,
               pct_gp = gpar(fontsize = 6),
               show_column_names = TRUE,
               row_names_gp = gpar(fontsize = 6),
               column_names_gp = gpar(fontsize = 10),
               remove_empty_columns = FALSE,
               axis_gp = gpar(fontsize = 8),
               column_title = "Candidate Driver Gene Mutation Profiles Across MD Lymphomas",
               column_title_gp = gpar(fontsize = 18),
               heatmap_legend_param = list(title = "Alternations", at = c("SNV", "INDEL"), 
                                           labels = c("Non-Synonymous SNV", "Non-Synonymous INDEL"), 
                                           nrow = 1, title_position = "leftcenter"))

# Save the oncoprint for top 30 mutated genes
pdf("./figures/oncoprint_nonsyn_somatic_snvs_indels.pdf")
plot <- draw(ht, heatmap_legend_side = "bottom")     
dev.off()

# Save the oncoprint for Ca regulation
#pdf("./figures/ca_regulation_somatic_snvs_indels.pdf")
#plot <- draw(ht, heatmap_legend_side = "bottom")     
#dev.off()

# Save the oncoprint for Fanconi Anemia
#pdf("./figures/fanconi_anemia_somatic_snvs_indels.pdf")
#plot <- draw(ht, heatmap_legend_side = "bottom")     
#dev.off()

