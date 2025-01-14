#Download and Install packages
install.packages('BiocManager')
library(BiocManager)

install(c('tximport', 'DESeq2', 'biomaRt', 'pheatmap'))
install.packages('RColorBrewer')
install.packages('tidyverse')

#Library Packages
library(tximport)
library(DESeq2)
library(biomaRt)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)


#2.3: Import
sample_table=read_csv('https://raw.githubusercontent.com/sjcockell/mmb8052/main/practicals/practical_08/data/sample_table.csv')

#Generates a "files" vector containing addresses and names of each sample
files=pull(sample_table, Run) 
files=paste0('counts/', files, '/quant.sf')
names(files)=pull(sample_table, Run)

#Imports quant.sf files into a gene-annotated gene abundance estimate list (txi)
gene_map=read_csv('https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/extdata/gene_map.csv')
txi=tximport(files, type='salmon',tx2gene=gene_map, ignoreTxVersion=TRUE)


#2.4: Normalization and Statistics

#DESeq2() function run as its individual components

#Count list from TXimport converted into DESEq dataset labelling columns with sample IDs and specifying experimental groups (I.e. naïve, Allo2h and Allo24h for contrast)
dds=DESeqDataSetFromTximport(txi, colData = sample_table, design = ~ Group)
#Size factor estimation to correct for library size
dds=estimateSizeFactors(dds)
#Estimate genewise dispersions to model variability
dds=estimateDispersions(dds)
#Calculates Wald's negative binomial test
dds=nbinomWaldTest(dds)

#Generates Figure 1. Plots count dispersions, fitted mean and final counts after empirical Bayes shrinkage.
plotDispEsts(dds,
             xlab="Mean of Normalized Counts", 
             ylab="Dispersion")


#2.5: Quality Control

#Rlog transforms dataset, reducing the impact of very highly and lowly expressed genes on heatmap and PCA presentaion
rld=rlog(dds)

#Plots transformed data in a PCA plot with a colourblind suitable colour palette from RColourBrewer. Adjustments made using ggplot functions
plotPCA(rld, intgroup="Group")+
  geom_point(size=3, alpha=1, )+
  scale_color_brewer(palette="Dark2")+
  labs(colour="Sample Group")+
  ggtitle("Normalized Dataset PCA")+
  theme_classic()

#Works out Euclidean statistical distance between all samples and saves these figures in a matrix
sample_distance = dist(t(assay(rld)), method="euclidian")
sample_distance_matrix = as.matrix(sample_distance)

#Generates sequential color scale
heatmap_colours = colorRampPalette(brewer.pal(9, "YlOrRd"))(255)

#Colours column annotation with PCA colour palette
heatmap_annotation = data.frame(Group=colData(dds)[,c('Group')], row.names=rownames(colData(dds)))
col_colours = list(
  Group = c(
    "Allo24h" = "darkcyan",
    "Naive" = "darkslateblue",
    "Allo2h" = "darkorange3")) #Manual version of colourbrewer Dark2 scheme for colourblind

#Produces heatmap clustering by sample distance with Group annotation and no borders
pheatmap(sample_distance_matrix, width=1, height=1, 
         clustering_distance_rows=sample_distance,
         clustering_distance_cols=sample_distance,
         annotation_col=heatmap_annotation,
         annotation_row=heatmap_annotation,
         annotation_colors=col_colours,
         col=heatmap_colours,
         border_color=NA)


#2.6: Volcano Plots

#Done the same for all 3 contrasts excluding first line

#24h vs. Naive (24h)

#In the 3 character vector, Groups tells results to contrast based on the Group column, and the other two specify the numerator (24h) and denominator (naive)
results_table24h = results(dds, contrast= c('Group', 'Allo24h', 'Naive'))
#Send data to tibble for filtering
results_tibble24h = as_tibble(results_table24h, rownames='ensembl_gene_id')
#Filters any "NA" figures from data
filtered_results24h = filter(results_tibble24h, complete.cases(results_tibble24h))
#Adds a column constinaing absolute transformed log2FC (negative values converted to positive)
filtered_results24h = mutate(filtered_results24h, absFC=abs(log2FoldChange))
#Adds a column displaying TRUE if absFC < 1 and and FALSE if absFC > 1 and another column displaying TRUE if padj < 0.05 and  FALSE if padj < 0.05
filtered_results24h = mutate(filtered_results24h, SignificantPval=padj<0.05, SignificantLogFC=absFC>1)

#Uses these new columns to display type of significance of each gene's expression.
plot_results24h <- mutate(
  filtered_results24h,
  Significance = case_when(
    SignificantPval & SignificantLogFC ~ "Both",
    SignificantPval & !SignificantLogFC ~ "P-Value",
    !SignificantPval & SignificantLogFC ~ "logFC",
    TRUE ~ "Neither" #If neither = TRUE col displays "Neither"
  )
)

#Adds a column conaining the -log10 of padj for data transformation in volcano plot 
plot_results24h = mutate(plot_results24h, logPVal=-log10(padj))

#Volcano plot. Significance column is referenced in the aes(colour=) argument to colour points based on significance. X and Y intercepts introduced at Signficance cutoff points. X and Y axis limited to keep plot size consistent
ggplot(plot_results24h, aes(x=log2FoldChange, y=logPVal)) +
  geom_point(aes(colour=Significance), alpha=0.3, size=2)+
  scale_color_manual(values=c("darkcyan","darkslateblue","black", "darkorange3"))+ 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="black")+
  geom_vline(xintercept=c(1, -1), linetype="dashed", colour="black")+
  xlim(-15, 15)+
  ylim(-0.05, 200)+
  ggtitle("24h vs. Naïve")+
  theme_classic()


#2h vs. Naive (2h)

results_table2h = results(dds, contrast= c('Group', 'Allo2h', 'Naive'))
results_tibble2h = as_tibble(results_table2h, rownames='ensembl_gene_id')
filtered_results2h = filter(results_tibble2h, complete.cases(results_tibble2h))
filtered_results2h = mutate(filtered_results2h, absFC=abs(log2FoldChange))
filtered_results2h = mutate(filtered_results2h, SignificantPval=padj<0.05, SignificantLogFC=absFC>1)

plot_results2h <- mutate(
  filtered_results2h,
  Significance = case_when(
    SignificantPval & SignificantLogFC ~ "Both",
    SignificantPval & !SignificantLogFC ~ "P-Value",
    !SignificantPval & SignificantLogFC ~ "logFC",
    TRUE ~ "Neither"
  )
)

plot_results2h = mutate(plot_results2h, logPVal=-log10(padj))

ggplot(plot_results2h, aes(x=log2FoldChange, y=logPVal)) +
  geom_point(aes(colour=Significance), alpha=0.3, size=2)+
  scale_color_manual(values=c("darkcyan","darkslateblue","black", "darkorange3"))+ 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="black")+
  geom_vline(xintercept=c(1, -1), linetype="dashed", colour="black")+
  xlim(-14, 14)+
  ylim(-0.05, 200)+
  ggtitle("2h vs. Naïve")+
  theme_classic()


#24h vs. 2h (vs)

results_tablevs = results(dds, contrast= c('Group', 'Allo24h', 'Allo2h'))
results_tibblevs = as_tibble(results_tablevs, rownames='ensembl_gene_id')
filtered_resultsvs = filter(results_tibblevs, complete.cases(results_tibblevs))
filtered_resultsvs = mutate(filtered_resultsvs, absFC=abs(log2FoldChange))
filtered_resultsvs = mutate(filtered_resultsvs, SignificantPval=padj<0.05, SignificantLogFC=absFC>1)

plot_resultsvs <- mutate(
  filtered_resultsvs,
  Significance = case_when(
    SignificantPval & SignificantLogFC ~ "Both",
    SignificantPval & !SignificantLogFC ~ "P-Value",
    !SignificantPval & SignificantLogFC ~ "logFC",
    TRUE ~ "Neither"
  )
)

plot_resultsvs = mutate(plot_resultsvs, logPVal=-log10(padj))

ggplot(plot_resultsvs, aes(x=log2FoldChange, y=logPVal)) +
  geom_point(aes(colour=Significance), alpha=0.3, size=2)+
  scale_color_manual(values=c("darkcyan","darkslateblue","black", "darkorange3"))+ 
  geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="black")+
  geom_vline(xintercept=c(1, -1), linetype="dashed", colour="black")+
  xlim(-24, 24)+
  ylim(-0.05, 200)+
  ggtitle("24h vs. 2h")+
  theme_classic()


#2.7: Result Table Gene Annotation

#Accessing Ensembl mouse gene annotation release 108
ensembl108 = useEnsembl(biomart="ensembl", version=108)
ensembl108 = useDataset("mmusculus_gene_ensembl", mart=ensembl108)

#24h vs. Naive (24h)

#Selecting attributes to add to table
annotation = getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 
                                'start_position', 'end_position', 
                                'strand', 'gene_biotype', 'external_gene_name', 
                                'description'), 
                   filters = 'ensembl_gene_id', values = filtered_results24h$ensembl_gene_id, 
                   mart = ensembl108)

#Adding table of attributes to results table
annot_results24h = left_join(filtered_results24h, annotation)
annot_results24h = arrange(annot_results24h, padj)

#Selecting for significantly differentially expressed genes
degs24h = filter(annot_results24h, abs(log2FoldChange) > 1 & padj < 0.05)
write.csv(degs24h, "Data/24h_Naive.csv", row.names=FALSE)

#2h vs. Naive (2h)

annotation = getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 
                                'start_position', 'end_position', 
                                'strand', 'gene_biotype', 'external_gene_name', 
                                'description'), 
                   filters = 'ensembl_gene_id', values = filtered_results2h$ensembl_gene_id, 
                   mart = ensembl108)

annot_results2h = left_join(filtered_results2h, annotation)
annot_results2h = arrange(annot_results2h, padj)

degs2h = filter(annot_results2h, abs(log2FoldChange) > 1 & padj < 0.05)
write.csv(degs2h, "Data/2h_Naive.csv", row.names=FALSE)

#24h vs. 2h (vs)

annotation = getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 
                                'start_position', 'end_position', 
                                'strand', 'gene_biotype', 'external_gene_name', 
                                'description'), 
                   filters = 'ensembl_gene_id', values = filtered_resultsvs$ensembl_gene_id, 
                   mart = ensembl108)

annot_resultsvs = left_join(filtered_resultsvs, annotation)
annot_resultsvs = arrange(annot_resultsvs, padj)

degsvs = filter(annot_resultsvs, abs(log2FoldChange) > 1 & padj < 0.05)
write.csv(degsvs, "Data/24h_2h.csv", row.names=FALSE)


#2.8: Investigating Outliers

#Strongly upregulated(Gm14440/ENSMUSG00000078901) and downregulated at 24h(Gm2541/ENSMUSG00000089698) outliers selected from DESeq dataset to scrutinize consistency between samples (cols) and groups. drop=FALSE ensures resulting information is stored in a matrix
Outliers24h = counts(dds)[c("ENSMUSG00000089698", "ENSMUSG00000078901", "ENSMUSG00000107383", "ENSMUSG00000030854"), 
                          c("SRR7457551", "SRR7457552", "SRR7457561", "SRR7457562"),
                          drop=FALSE]

Outliers2h = counts(dds)[c("ENSMUSG00000089698", "ENSMUSG00000078901", "ENSMUSG00000107383", "ENSMUSG00000030854"), 
                         c("SRR7457553", "SRR7457554", "SRR7457555", "SRR7457556"),
                         drop=FALSE]

OutliersNaive = counts(dds)[c("ENSMUSG00000089698", "ENSMUSG00000078901", "ENSMUSG00000107383", "ENSMUSG00000030854"), 
                            c("SRR7457557", "SRR7457558", "SRR7457559", "SRR7457560"),
                            drop=FALSE]