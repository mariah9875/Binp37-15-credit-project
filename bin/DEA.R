
# load packages
library(DESeq2)
library(pheatmap)
library(plotly)
library(data.table)

# import data
countdata <- read.table("counts_default.csv")
gene_id <- read.csv("gene_id.txt")
gene_name <- read.csv("gene_name.txt")
gene_type <- read.csv("gene_type.txt")


## Data preparation
#---------------------------------------------------------------------
# how many lines 
nrow(gene_id)
nrow(gene_name)
nrow(gene_type)

# concatinate gene_id, gene_name, gene_type and make unique
gene_info <- data.frame(gene_id, gene_name, gene_type)
gene_info_unique <- unique(gene_info)
# *Delete columns
countdata <- data.frame(countdata$V1 ,countdata$V7, countdata$V8, countdata$V9, countdata$V10, countdata$V11, countdata$V12, countdata$V13, countdata$V14, countdata$V15, countdata$V16, countdata$V17)

# *Rename columns
countdata <- countdata %>% rename("gene_id"=countdata.V1, "DA026_S7"=countdata.V7, "DA027_S1"=countdata.V8, "DA028_S2"=countdata.V9, "DA035_S9"=countdata.V10, "DA036_S3"=countdata.V11, "DA050_S8"=countdata.V12, "DA056_S5"=countdata.V13, "DA057_S10"=countdata.V14, "DA058_S11"=countdata.V15, "DA059_S6"=countdata.V16, "DA061_S4"=countdata.V17)
gene_info <- gene_info_unique %>% rename("gene_id"=X.gene_id, "gene_name"=X.gene_name, "gene_type"=X.gene_type)
# delete first row
countdata = countdata[-1,]

# Remove characters and add gene_names to countdata
gene_info$gene_id <- gsub(";", "", gene_info$gene_id) 
gene_info$gene_name <- gsub(";", "", gene_info$gene_name)
gene_info$gene_type <- gsub(";", "", gene_info$gene_type)

#sort both files before merge together
gene_info <- gene_info[order(gene_info$gene_id),]
countdata <- countdata[order(countdata$gene_id),]
# check if identical gene_id
all(gene_info$gene_id == countdata$gene_id)

# merge countdata and gene_info
countdata_gene <- merge(gene_info,countdata, by="gene_id")

# gene_name as rownames
rownames(countdata_gene) <- make.unique(countdata_gene$gene_name)

# delete 3 columns
countdata_gene <- countdata_gene[-1]
countdata_gene <- countdata_gene[-1]
countdata_gene <- countdata_gene[-1]
countdata <- countdata_gene

# *Metadata
(condition <- factor(c("Chimp", "Chimp", "Human", "Human", "Human", "Human", "Human", "Human", "Chimp", "Chimp", "Chimp")))
(coldata <- data.frame(row.names=colnames(countdata),sample=colnames(countdata), condition))
#------------------------------------------------------------------------


#-----------------------------------------------------------------------------------
## DESeq2 analysis

# Make sure the counts are numeric
sapply(countdata, class)
countdata <- transform(countdata, DA026_S7=as.numeric(DA026_S7), DA027_S1=as.numeric(DA027_S1), DA028_S2=as.numeric(DA028_S2), DA035_S9=as.numeric(DA035_S9), DA036_S3=as.numeric(DA036_S3), DA050_S8=as.numeric(DA050_S8), DA056_S5=as.numeric(DA056_S5), DA057_S10=as.numeric(DA057_S10), DA058_S11=as.numeric(DA058_S11), DA059_S6=as.numeric(DA059_S6), DA061_S4=as.numeric(DA061_S4))
sapply(countdata, mode)

# DESeq object
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~condition)
# Filter out lowly expressed counts
dds <- dds[rowSums(counts(dds))>0.5,]

# *Set factor and run DESeq
dds$condition <- factor(dds$condition, levels = c('Chimp', 'Human'))
dds <- DESeq(dds)

# DESeq results object
res <- results(dds)
summary(res)
res

# Regularized log transformation for visualization
# Normalization
norm <- counts(dds, normalized=TRUE)
norm_dds <- norm
rld <- rlogTransformation(dds)

# Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="condition_Human_vs_Chimp", type="apeglm")

# Adjusted p-value, up and down-regulated genes
table(res$padj < 0.05)
table(res$padj < 0.05 & res$log2FoldChange > 1)
table(res$padj < 0.05 & res$log2FoldChange < -1)


# Print to table expressed, up and downregulated genes
up_regulated <- subset(res, res$padj < 0.05 & res$log2FoldChange > 1)
down_regulated <- subset(res, res$padj < 0.05 & res$log2FoldChange < -1)
all_genes <- subset(res)

write.table(up_regulated, file="up_regulated.txt", sep = "\t", quote = FALSE)
write.table(down_regulated, file= "down_regulated.txt", sep = "\t", quote = FALSE)
write.table(all_genes, file="all_regulated.txt", sep = "\t", quote = FALSE)

## Visualization
# Principal Component Analysis
plotPCA(rld) + theme_classic()

# MA plot
plotMA(res, alpha=0.05, ylim=c(-2,2))
# MA plot with shrunken log2 fold changes.
plotMA(resLFC, ylim=c(-2,2))

# Mean plot
norm_dds <- as.data.frame(norm_dds)
norm_dds_test <- norm_dds
# Mean of human and chimp
human_mean <- rowMeans(norm_dds_test[,subset(coldata, coldata$condition == "Human")$sample])
chimp_mean <- rowMeans(norm_dds_test[,subset(coldata, coldata$condition == "Chimp")$sample])

# Genes found in humans and chimp
human_genes_expressed <- norm_dds_test[,subset(coldata, coldata$condition == "Human")$sample]
human_genes_expressed <- human_genes_expressed[rowSums(human_genes_expressed)>0,]
chimp_genes_expressed <- norm_dds_test[,subset(coldata, coldata$condition == "Chimp")$sample]
chimp_genes_expressed <- chimp_genes_expressed[rowSums(chimp_genes_expressed)>0,]

# Mark up-regulate genes red and down-regulated genes blue
norm_dds_test$cex <- rep(1, nrow(norm_dds_test))
norm_dds_test[rownames(up_regulated),] <- 2
norm_dds_test[rownames(down_regulated),] <- 3
norm_dds_test$colour <- "black"
norm_dds_test$colour[norm_dds_test$cex==2] <- "red"
norm_dds_test$colour[norm_dds_test$cex==3] <- "blue"
# Plot with scatter plot
plot(log2(chimp_mean + 0.5),
     log2(human_mean + 0.5),
     col=norm_dds_test$colour,
     pch=16,
     xlab = "log2(mean Chimp)",
     ylab = "log2(mean Human)")
legend("bottomright",c("upregulated","not significant", "downregulated"),cex= .8, fill=c("red", "black", "blue"))


# Heatmaps
# select condition to be compared in heatmaps
coldata <- coldata[order(coldata$condition),]
# up and down-regulated genes
signdiff_up <- norm[rownames(subset(res, res$padj < 0.05 & res$log2FoldChange > 1)),]
signdiff_low <- norm[rownames(subset(res, res$padj < 0.05 & res$log2FoldChange < -1)),]

# Heatmap of top variance of genes
topVarGenes <- head(order(-rowVars(assay(rld))),20)
mat <- assay(rld)[topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("condition"),drop=FALSE])
pheatmap(mat, annotation_col = df)
# Heatmap of upregulation
pheatmap(log2(signdiff_up[,as.character(coldata$sample)]+0.5), cluster_cols = F, cluster_rows = F, annotation_col = df, fontsize = 8)
# Heatmap of downregulation
pheatmap(log2(signdiff_low[,as.character(coldata$sample)]+0.5), cluster_cols = F, cluster_rows = F, annotation_col = df, fontsize = 8)



# pie-chart of expressed small RNAs
# prepare data for pie chart
# Extract gene types from gene_info 
rownames(gene_info) <- make.unique(gene_info$gene_name)
type_all <- merge(gene_info, norm_dds, by=0)

type_upregulated <- merge(gene_info, signdiff_up, by=0)
type_downregulated <- merge(gene_info, signdiff_low, by=0)
write.table(type_upregulated, file="up_type.txt", sep = "\t", quote = FALSE)

type_human <- merge(gene_info, human_genes_expressed, by=0)
type_chimp <- merge(gene_info, chimp_genes_expressed, by=0)

freq_type <- (table(type_all$gene_type))
freq_type <- as.data.frame(freq_type)
freq_type <- freq_type %>% rename("geneType"=Var1)

# Ascendin order
freq_type <- freq_type[order(-freq_type$Freq),]
# Write to file all small RNAs type
write.table(freq_type, file="all_type.txt", sep = ",", quote = FALSE, row.names = F)

## *Pool together small groups of rna types
#--------------------------------------------------------------------------------------
## This pools together small groups of rna types
# make freq_type a character
freq_type$geneType <- as.character(freq_type$geneType)
# Pool together small goups of rna types into the group Others
freq_type$geneType[freq_type$geneType == "processed_pseudogene"] <- "pseudogenes"
freq_type$geneType[freq_type$geneType == "pseudogene"] <- "pseudogenes"
freq_type$geneType[freq_type$geneType == "IG_V_pseudogene"] <- "pseudogenes"
freq_type$geneType[freq_type$geneType == "polymorphic_pseudogene"] <- "pseudogenes"
freq_type$geneType[freq_type$geneType == "bidirectional_promoter_lncRNA"] <- "lncRNA"
freq_type$geneType[freq_type$geneType == "macro_lncRNA"] <- "lncRNA"
freq_type$geneType[freq_type$geneType == "Mt_rRNA"] <- "rRNA"
freq_type$geneType[freq_type$geneType == "rRNA_pseudogene"] <- "rRNA"
freq_type$geneType[freq_type$geneType == "Mt_tRNA"] <- "tRNA"
freq_type$geneType[freq_type$geneType == "unprocessed_pseudogene"] <- "Others"
freq_type$geneType[freq_type$geneType == "unitary_pseudogene"] <- "Others"
freq_type$geneType[freq_type$geneType == "transcribed_unprocessed_pseudogene"] <- "Others"
freq_type$geneType[freq_type$geneType == "transcribed_unitary_pseudogene"] <- "Others"
freq_type$geneType[freq_type$geneType == "transcribed_processed_pseudogene"] <- "Others"
freq_type$geneType[freq_type$geneType == "TR_V_gene"] <- "Others"
freq_type$geneType[freq_type$geneType == "TR_J_gene"] <- "Others"
freq_type$geneType[freq_type$geneType == "TEC"] <- "Others"
freq_type$geneType[freq_type$geneType == "3prime_overlapping_ncRNA"] <- "Others"
freq_type$geneType[freq_type$geneType == "processed_transcript"] <- "Others"
freq_type$geneType[freq_type$geneType == "ribozyme"] <- "Others"
freq_type$geneType[freq_type$geneType == "misc_RNA"] <- "Others"
freq_type$geneType[freq_type$geneType == "scaRNA"] <- "Others"
freq_type$geneType[freq_type$geneType == "snRNA"] <- "Others"
freq_type$geneType[freq_type$geneType == "snoRNA"] <- "Others"
freq_type$geneType[freq_type$geneType == "sense_overlapping"] <- "Others"
freq_type$geneType[freq_type$geneType == "sense_intronic"] <- "Others"
#----------------------------------------------------------------------------------------

# Delete duplicate rows
type <- tapply(freq_type$Freq, freq_type$geneType, FUN=sum)
z <- as.data.frame(type)

# make row names a column
z <- setDT(z, keep.rownames = TRUE)[]

# up-regulated types
freq_up <- (table(type_upregulated$gene_type))
freq_up <- as.data.frame(freq_up)
freq_up <- freq_up %>% rename("geneType"=Var1)
# down-regulated types
freq_down <- (table(type_downregulated$gene_type))
freq_down <- as.data.frame(freq_down)
freq_down <- freq_down %>% rename("geneType"=Var1)


# Plot pie chart with all types of expressed small RNAs
plot_ly(data = z, labels = ~rn, values= ~type,
          type = 'pie', sort= FALSE,
          marker = list(colors=colors, line = list(color="black", width = 0.5))) %>%
          layout(title = "Expressed small RNAs")
# up-regulated types pie chart 
plot_ly(data = freq_up, labels = ~geneType, values= ~Freq,
          type = 'pie', sort= FALSE,
          marker = list(colors=colors, line = list(color="black", width = 0.5))) %>%
          layout(title = "Up-regulated small RNAs")

# down-regulated types pie chart
plot_ly(data = freq_down, labels = ~geneType, values= ~Freq,
          type = 'pie', sort= FALSE,
          marker = list(colors=colors, line = list(color="black", width = 0.5))) %>%
          layout(title = "Down-regulated small RNAs")


