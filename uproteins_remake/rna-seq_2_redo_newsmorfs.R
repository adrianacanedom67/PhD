library(DESeq2)
library(airway)
library(tidyverse)


counts_data <- read.csv("/home/adriana/uproteins/redo/differential_expression/isolates/two_strains_counts_filtered_updated_2_RENAMED.csv", check.names = FALSE)
# counts_data <- read_excel("/home/adriana/uproteins/isolates/two_strains/two_strains_counts_filtered_updated_2_RENAMED.xlsx")
head(counts_data)

#set the genes as rownames and remove the first column so that the rownames aren't numbers anymore
rownames(counts_data) <- counts_data[, 1]
counts_data <- counts_data[, -1]
head(counts_data)

#load column data - sample information
colData <- read.csv("/home/adriana/uproteins/isolates/two_strains/two_strains_sample_information.txt")
rownames(colData) <- colData$X

# remove the random X of the column names in counts_data
cleaned_counts_names <- gsub("^X", "", colnames(counts_data)) 
colnames(counts_data) <- cleaned_counts_names

#check whether the colnames from count_data match with colData rownames
all(colnames(counts_data) %in% rownames(colData))

#create DESeq DataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = colData, design = ~ strain)

#perform variance satabilizing transformation
vsd <- vst(dds, blind=TRUE) # normalized raw counts

pcadata <- plotPCA(vsd, intgroup="strain", returnData=TRUE)
plotPCA(vsd, intgroup = c("strain"), returnData = FALSE)

#run DESeq analysis
dds_results <- DESeq(dds)
resultsNames(dds_results) 

# this one is for h37rv vs 8462
res_h37rv_vs_8462 <- results(dds_results, name = "strain_h37rv_vs_8462")

# filter p value
res <- subset(res_h37rv_vs_8462, padj<0.05)
res$log2FoldChange > 1 # this boolean checks which rows have a log2foldchange >1 and returns T or F, but doesn't modify the dataframe itself
res <- res[res$log2FoldChange > 1, ] # filters out to include only rows with log2foldchange >2
#filter out Rvs
res <- res[!grepl("Rv", rownames(res)), ]
res <- res[!grepl("RV", rownames(res)), ]
#extract row names from res dataframe, and store the variance stabilized transformedcounts for the DE genes 
diff_1 <- vsd[rownames(res), ]




# for 94789 vs h37rv:
res_94789_vs_h37rv <- results(dds_results, name = "strain_94789_vs_h37rv", contrast=c("strain", "94789", "h37rv"))
# filter p value
res_2 <- subset(res_94789_vs_h37rv, padj<0.05)
res_2$log2FoldChange > 1
res_2 <- res_2[res_2$log2FoldChange > 1, ]
# filter out Rvs
res_2 <- res_2[!grepl("Rv", rownames(res_2)),]
res_2 <- res_2[!grepl("RV", rownames(res_2)),]

# Ensure there are no NA values in rownames(res_2) after filtering
if (any(is.na(rownames(res_2)))) {
  stop("There are NA values in the row names of res_2 after filtering")
}

diff_2 <- vsd[rownames(res_2)]


#concatenate diff_1 and diff_2 to compare h37rv with both strains
# diff_combined <- cbind(diff_1, diff_2)
concatenated_diff <- rbind(diff_1, diff_2)




# for 8462 vs 94789:
res_94789_vs_8462 <- results(dds_results, contrast=c("strain", "94789", "8462"))
res_3 <- subset(res_94789_vs_8462, padj<0.05)
res_3$log2FoldChange > 1
res_3 <- res_3[res_3$log2FoldChange > 1, ]
#filter out RVs
res_3 <- res_3[!grepl("Rv", rownames(res_3)),]
res_3 <- res_3[!grepl("RV", rownames(res_3)),]
diff_3 <- vsd[rownames(res_3)]


directory_path <- "/home/adriana/uproteins/isolates/two_strains/new_smorfs"

write.csv(as.data.frame(res_94789_vs_8462), file = paste0(directory_path, "res_94789_vs_8462.csv"), row.names = TRUE)
write.csv(as.data.frame(res_h37rv_vs_8462), file = paste0(directory_path, "res_h37rv_vs_8462.csv"), row.names = TRUE)
write.csv(as.data.frame(res_94789_vs_h37rv), file = paste0(directory_path, "res_94789_vs_h37rv.csv"), row.names = TRUE)

df <- as.data.frame(colData(dds)[,c("strain")])

#heatmap of the count matrix
library("pheatmap")
library("RColorBrewer")
library("grid")

# select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]

df <- as.data.frame(colData(dds)[,c("strain")])
# condition_colors <- brewer.pal(nlevels(df$condition), "Set1")]
colData_filtered <- colData[, !(names(colData) %in% c("type", "X"))]
colData_filtered$type <- NA
colData_filtered$X <- NA
colData_filtered <- colData_filtered[1:12]
colData_df <- as.data.frame(colData_filtered)
colData_df <- colData_df[!names(colData_df) %in% c("type", "X")]
colData_matrix <- t(as.matrix(colData_filtered))
# Check if dimensions match
if (ncol(assay(diff_1)) != ncol(colData_df)) {
  stop("Number of columns in expression matrix and annotation data don't match.")
}
my_colors <- colorRampPalette(c("blue", "white", "red"))(100)

condition_index <- which(names(colData) == "condition")
annotation_legend <- rep(FALSE, ncol(colData))
annotation_legend[condition_index] <- TRUE

colData_2 <- colData
colData_2 <- colData[, -which(names(colData) == "type")]
rownames(colData_2) <- gsub("_quant.sf", "", rownames(colData_2))
colData_2$X <- sub("_quant.sf", "", colData_2$X)
colData_2 



heatmap_rows <- rownames(colData_2)
heatmap_rows


heatmap_concatenated <- pheatmap(assay(concatenated_diff),
                                 main = "Reference vs Resistant strains",
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 show_rownames = TRUE,
                                 annotation_col = NULL,  # Remove column annotations
                                 fontsize = 7,
                                 color = my_colors,
                                 annotation_legend = FALSE,
                                 annotation_names_row = FALSE,
                                 annotation_names_col = FALSE,  # Remove column annotation names
                                 angle_col = 45,
                                 labels_col = heatmap_rows,
                                 scale = "none")


# library("ComplexHeatmap")
# complex_heat <- Heatmap(assay(diff), show_row_names=TRUE, show_column_names = FALSE)
# draw(complex_heat)

#heatmap of the sample-to-sample distances
# sampleDists <- dist(t(assay(dds)))
