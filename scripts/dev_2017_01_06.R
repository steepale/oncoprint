# Script based on: http://bioconductor.org/packages/devel/bioc/vignettes/ComplexHeatmap/inst/doc/s8.oncoprint.html

setwd("/Users/Alec/Documents/Bioinformatics/MDV_Project/oncoprint")
source("https://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
library(ComplexHeatmap)
#install.packages("seqminer")
#library(seqminer)

# Use seq miner to read VCF file and create a binary data matrix
# Download refFlat file from http://hgdownload.soe.ucsc.edu/goldenPath/galGal5/database/
# Date of file is 19-Feb-2017 (GenBank database)
# Filename: ./data/refFlat.txt.gz
# Load the refFlat file into R:
#geneFile = "./data/refFlat.txt.gz"
# Load the somatic SNV file into R
#vcfFile = "./data/somatic_snvs_tumors_vep.vcf.gz"
#cfh <- readVCFToListByGene(vcfFile, geneFile, geneName="IKZF1", "", "", "", c("GT"))



mat_snvs = read.table("/Users/Alec/Documents/Bioinformatics/MDV_Project/oncoprint/data/somatic_snvs_oncoprint.mat", header = TRUE, stringsAsFactors=FALSE, sep = "\t")
mat_snvs[is.na(mat_snvs)] = "NULL"
rownames(mat_snvs) = mat_snvs[, 1]
rownames(mat_snvs)
mat_snvs = mat_snvs[, -1]
colnames(mat_snvs)
mat_snvs = as.matrix(mat_snvs)

mat_indels = read.table("/Users/Alec/Documents/Bioinformatics/MDV_Project/oncoprint/data/somatic_indels_oncoprint.mat", header = TRUE, stringsAsFactors=FALSE, sep = "\t")
rownames(mat_indels) = mat_indels[, 1]
rownames(mat_indels)
mat_indels = mat_indels[, -1]
colnames(mat_indels)
mat_indels = as.matrix(mat_indels)

mat_amps = read.table("/Users/Alec/Documents/Bioinformatics/MDV_Project/oncoprint/data/somatic_amps_oncoprint.mat", header = TRUE, stringsAsFactors=FALSE, sep = "\t")
rownames(mat_amps) = mat_amps[, 1]
rownames(mat_amps)
mat_amps = mat_amps[, -1]
colnames(mat_amps)
mat_amps = as.matrix(mat_amps)

mat_dels = read.table("/Users/Alec/Documents/Bioinformatics/MDV_Project/oncoprint/data/somatic_dels_oncoprint.mat", header = TRUE, stringsAsFactors=FALSE, sep = "\t")
rownames(mat_dels) = mat_dels[, 1]
rownames(mat_dels)
mat_dels = mat_dels[, -1]
colnames(mat_dels)
mat_dels = as.matrix(mat_dels)

mat_loh = read.table("/Users/Alec/Documents/Bioinformatics/MDV_Project/oncoprint/data/somatic_loh_oncoprint.mat", header = TRUE, stringsAsFactors=FALSE, sep = "\t")
rownames(mat_loh) = mat_loh[, 1]
rownames(mat_loh)
mat_loh = mat_loh[, -1]
colnames(mat_loh)
mat_loh = as.matrix(mat_loh)

mat = matrix(paste(mat_snvs,mat_indels,mat_amps,mat_dels,mat_loh,sep=";"),nrow=nrow(mat_snvs),ncol=ncol(mat_snvs))
dimnames(mat) <- list(rownames(mat_indels), colnames(mat_indels))
mat <- replace(mat, mat == "NULL;NULL;NULL;NULL;NULL", "")
mat <- replace(mat, mat == "SNV;NULL;NULL;NULL;NULL", "SNV")
mat <- replace(mat, mat == "SNV;INDEL;NULL;NULL;NULL", "SNV;INDEL")
mat <- replace(mat, mat == "SNV;NULL;AMP;NULL;NULL", "SNV;AMP")
mat <- replace(mat, mat == "SNV;NULL;NULL;DEL;NULL", "SNV;DEL")
mat <- replace(mat, mat == "NULL;INDEL;NULL;NULL;NULL", "INDEL")
mat <- replace(mat, mat == "NULL;NULL;AMP;NULL;NULL", "AMP")
mat <- replace(mat, mat == "NULL;INDEL;AMP;NULL;NULL", "INDEL;AMP")
mat <- replace(mat, mat == "NULL;INDEL;NULL;DEL;NULL", "INDEL;DEL")
mat <- replace(mat, mat == "NULL;NULL;NULL;DEL;NULL", "DEL")
mat <- replace(mat, mat == "NULL;NULL;NULL;NULL;LOH", "LOH")
mat <- mat[-6,]
mat <- mat[-5,]
mat <- mat[-7,]

mat[is.na(mat)] = ""

alter_fun = list(background = function(x, y, w, h, v) 
        {grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
                   gp = gpar(fill = "#CCCCCC", col = NA))}, 
        SNV = function(x, y, w, h) 
                {grid.rect(x, y, w*0.33, h*1.00, 
                           gp = gpar(fill = "blue", col = NA))},
        INDEL = function(x, y, w, h) 
                {grid.rect(x, y, w*0.66, h*0.90, 
                   gp = gpar(fill = "#008000", col = NA))},
        AMP = function(x, y, w, h) {
                grid.rect(x, y, w-unit(0.5, "mm"), h*0.50, gp = gpar(fill = "orange", col = NA))},
        DEL = function(x, y, w, h) {
                grid.rect(x, y, w-unit(0.5, "mm"), h*0.50, gp = gpar(fill = "red", col = NA))})

col = c("SNV" = "blue", "INDEL" = "#008000", "AMP" = "orange", "DEL" = "red")

ht = oncoPrint(mat,
               get_type = function(x) strsplit(x, ";")[[1]],
               alter_fun = alter_fun, 
               col = col,
               show_pct = TRUE,
               pct_gp = gpar(fontsize = 9),
               show_column_names = TRUE,
               row_names_gp = gpar(fontsize = 9),
               column_names_gp = gpar(fontsize = 10),
               remove_empty_columns = FALSE,
               axis_gp = gpar(fontsize = 8),
               column_title = "Candidate Driver Gene Mutation Profiles Across MD Lymphomas",
               column_title_gp = gpar(fontsize = 18),
               heatmap_legend_param = list(title = "Alternations", at = c("SNV", "INDEL", "AMP", "DEL"), 
                                           labels = c("Non-Synonymous SNV", "Non-Synonymous INDEL", "Amplification", "Deletion"), 
                                           nrow = 1, title_position = "leftcenter"))
draw(ht, heatmap_legend_side = "bottom")

# Create a small heatmap for IKZF1

mat = matrix(paste(mat_snvs,mat_indels,mat_amps,sep=";"),nrow=nrow(mat_snvs),ncol=ncol(mat_snvs))
dimnames(mat) <- list(rownames(mat_indels), colnames(mat_indels))
mat <- replace(mat, mat == "NULL;NULL;NULL", "")
mat <- replace(mat, mat == "SNV;NULL;NULL", "SNV")
mat <- replace(mat, mat == "NULL;INDEL;NULL", "INDEL")
mat <- replace(mat, mat == "NULL;NULL;AMP", "")
mat <- replace(mat, mat == "NULL;INDEL;AMP", "INDEL;AMP")

mat[is.na(mat)] = ""

mat_IKZF1 <- mat[1:8, 1:26]
mat_IKZF1 <- mat_IKZF1[-1:-5, 1:26]
mat_IKZF1 <- mat_IKZF1[-8, 1:26]
mat_IKZF1[2,1] <- "LOW"
mat_IKZF1[2,14] <- "LOW"
mat_IKZF1[2,17] <- "AMP;LOW"
mat_IKZF1[2,18] <- "LOW"
mat_IKZF1[2,19] <- "LOW"
mat_IKZF1[2,22] <- "LOW"


rownames(mat_IKZF1)[1] = ""
rownames(mat_IKZF1)[3] = ""
mat_IKZF1[1,1:26] = ""
mat_IKZF1[3,1:26] = ""


alter_fun = list(background = function(x, y, w, h, v) 
{grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
           gp = gpar(fill = "#CCCCCC", col = NA))}, 
SNV = function(x, y, w, h) 
{grid.rect(x, y, w*0.33, h*1.00, 
           gp = gpar(fill = "blue", col = NA))},
INDEL = function(x, y, w, h) 
{grid.rect(x, y, w*0.66, h*0.90, 
           gp = gpar(fill = "#008000", col = NA))},
AMP = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.50, gp = gpar(fill = "orange", col = NA))},
LOW = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.50, gp = gpar(fill = "red", col = NA))})

col = c("SNV" = "blue", "INDEL" = "#008000", "AMP" = "orange", "LOW" = "red")

ht = oncoPrint(mat_IKZF1,
               get_type = function(x) strsplit(x, ";")[[1]],
               alter_fun = alter_fun, 
               col = col,
               show_pct = TRUE,
               pct_gp = gpar(fontsize = 14),
               show_column_names = TRUE,
               row_names_gp = gpar(fontsize = 14),
               show_row_barplot = FALSE,
               top_annotation = NULL,
               column_names_gp = gpar(fontsize = 12),
               remove_empty_columns = FALSE,
               axis_gp = gpar(fontsize = 8),
               column_title_gp = gpar(fontsize = 18),
               column_title = "IKZF1 Loss-of-Function: Binding Domain Mutations or Low Gene Expression", 
               heatmap_legend_param = list(title = "Alternations", at = c("SNV", "INDEL", "AMP", "LOW"), 
                                           labels = c("Non-Synonymous SNV", "Non-Synonymous INDEL", "Amplification", "Low Gene Expression"), 
                                           nrow = 1, title_position = "leftcenter"))
draw(ht, heatmap_legend_side = "bottom")










