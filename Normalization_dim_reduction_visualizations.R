# Info --------------------------------------------------------------------
# Liam Hendrikse
# May 02, 2021
# Taylor lab


# Notes -------------------------------------------------------------------
# Code to generate visualizations and differential gene expression results 
# from Hendrikse et al, Nature, 2022


# Libraries ---------------------------------------------------------------
library(DESeq2)
library(plyr)
library(vsn)
library(umap)
library(ggplot2)
library(Rtsne)
library(gplots)
library(RColorBrewer)
library(DESeq2)
library(apeglm)
library(EnhancedVolcano)

reticulate::use_python("/opt/anaconda3/bin/python3.7", required = T)
library(reticulate)
py_config()


# Inputs ------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
count_path = args[1] #Path to output from STAR aligner
metadata_path = args[2] #Path to clinical metadata or sample annotations for samples in count_path


# Functions ---------------------------------------------------------------
plot_volc = function (res, sig = 0.05) {  
  res_sig=res[which(res$padj<sig),]
  res_ordered_sig <- res_sig[order(res_sig$log2FoldChange, decreasing = T),]
  
  #Output volcano plots
  max = max(res$log2FoldChange) + 0.25
  min = min(res$log2FoldChange) - 0.25
  ymax = max(-log10(res$padj), na.rm=TRUE) + 1
  title <- gsub("posterior SD: ", "", res@elementMetadata[[2]][3])

  #Enhanced volcano
  volc <- EnhancedVolcano(res,
                     lab = rownames(res),
                     x = 'log2FoldChange',
                     y = 'padj',
                     # selectLab = genes_highlight,
                     # colAlpha = c(ifelse(rownames(res) %in% genes_highlight, 0.8, 0.1)),
                     # drawConnectors = T,
                     # widthConnectors = 0.2,
                     # colConnectors = 'grey30',
                     FCcutoff = 1,
                     pCutoff = 0.05,
                     title = title,
                     subtitle = element_blank(),
                     caption  = element_blank(),
                     gridlines.major = F,
                     gridlines.minor = F,
                     #transcriptPointSize = 2,
                     legendPosition = "none",
                     xlim = c(min, max),
                     ylim = c(0, ymax))
  return(volc)

}


# Directories -------------------------------------------------------------
out_dir = "/outs/"
plot_dir = "/outs/plots/"
dir.create(out_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)
setwd(out_dir)


# Imports -----------------------------------------------------------------
counts = read.table(count_path, header=T, row.names=1, sep = "\t)
metadata = read.csv(metadata_path)
palette_plot = c("#FFF836", "#DBAF2B","#FFA920",  "#00F636","#D7EBCC","#007314", "white")
theme = theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),legend.title=element_blank(), legend.position = "none")


# Init DESeq object -------------------------------------------------------
dds = DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ 1 )

dds = dds[ rowSums(counts(dds)) > 1, ] #Optional filter
vsd = vst(dds, blind=F) #VST transformation and normalization


# Dimensionality reduction ------------------------------------------------

# PCA
plotPCA(vsd, intgroup=c("RNA_seq_subgroup_allMB_k4"))

# UMAP
metadata$RNA_seq_subgroup_allMB_k4 = factor(metadata$RNA_seq_subgroup_allMB_k4, 
                                            levels = c("Group3", "Group4"))

norm = vsd@assays@data@listData[[1]]
umap.custom = umap.defaults
umap.custom$n_neighbors = 12 

umap_obj = umap::umap(t(norm), config = umap.custom)

# Visualize clustering of subgroups on UMAP
gg_data = data.frame(UMAP1 = umap_obj$layout[,1] , 
                     UMAP2 = umap_obj$layout[,2] , 
                     Subgroup = metadata$RNA_seq_subgroup_allMB_k4)
ggplot(gg_data, aes(x= UMAP1, y = UMAP2, fill = group)) + 
  geom_point(colour = "black", shape = 21, size = 2) + 
  scale_fill_manual(values=palette_plot) +
  theme_classic() +
  theme

# Visualize gene expression of tumors on UMAP
gg_data = data.frame(UMAP1 = umap_obj$layout[,1] , 
                     UMAP2 = umap_obj$layout[,2] , 
                     Gene = umap_obj$data[,which(grepl("^OTX2_", colnames(umap_obj$data)))])
ggplot(gg_data, aes(x= UMAP1, y = UMAP2, color = Gene)) + 
  geom_point() + 
  scale_color_gradient2(midpoint=mean(gg_data$Gene), low="blue", mid="white",
                        high="red", space ="Lab" ) +
  #scale_color_gradient(low="lightgrey", high="blue") + 
  theme_classic()


# Gene expression correlation ---------------------------------------------

matrix_mod = as.matrix(t(umap_obj$data))
rownames(matrix_mod) = gsub("___.*", "", rownames(matrix_mod))
gene = as.numeric(matrix_mod["BARHL1",])
correlations = apply(matrix_mod,1,function(x){cor(gene,x)})

correlations = data.frame(Gene = names(correlations), 
                          Corr = as.numeric(correlations))
correlations = correlations[order(correlations$Corr, decreasing = T),]

# BARHL1 and DDX31 gene expression is expected to be correlated as they comprise a super enhancer in the developing cerebellum - positive control
plot(umap_obj$data[,which(grepl("^BARHL1_", colnames(umap_obj$data)))], 
     umap_obj$data[,which(grepl("^DDX31_", colnames(umap_obj$data)))])

# To confirm indivdual correlations
summary(lm(umap_obj$data[,which(grepl("^BARHL1_", colnames(umap_obj$data)))] ~ umap_obj$data[,which(grepl("^DDX31_", colnames(umap_obj$data)))]))


# Differential expression analysis ----------------------------------------

dds = DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ RNA_seq_subgroup_allMB_k4)
dds = dds[ rowSums(counts(dds)) > 10, ] 
dds = DESeq(dds, parallel=TRUE)
res = lfcShrink(dds, coef=2, type="apeglm")

# DEG visualization with volcano plot
volcano_plot = plot_volc(res, sig = 0.05)


sessionInfo()
rm(list=ls())




