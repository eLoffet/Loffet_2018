counts=read.csv("featcounts.txt", sep="", head=T, skip=1, row.names = "Geneid")
counts<-counts[,6:41]
count1<-as.matrix(counts)

#design experiment
library(DESeq2)
(condition <- factor(c(rep("d28HIO", 3), rep("d28HIOENS", 3), rep("spheroid", 3), rep("d14HIO", 3), rep("d14HIOENS", 2))))
(coldata <- data.frame(row.names=colnames(count1), condition))
dds=DESeqDataSetFromMatrix(countData = count1,colData = coldata,design =~condition)

#prefiltering
keep <- rowSums(counts(dds)) >= 10
keepdds <- dds[keep,]

#DESeq2
ddds <- DESeq(keepdds)
res_ddds=results(ddds)

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(ddds)


#pcaExplorer
library(pcaExplorer)
anno_df_biomart <- get_annotation(dds = ddds,
                                  biomart_dataset = "hsapiens_gene_ensembl",
                                  idtype = "ensembl_gene_id")
pcaExplorer(dds = ddds, rlt = rld,annotation=anno_df_biomart)

#pca
library(ggplot2)
pcaplot(rld, intgroup = "condition", ntop = 10000, returnData = FALSE, pcX = 1, pcY = 2, title = "PCA on 10000 genes", text_labels = FALSE, point_size = 6, ellipse = TRUE, ellipse.prob = 0.95) +
  theme(axis.text.y   = element_text(size=12), axis.text.x   = element_text(size=12), axis.title.y  = element_text(size=14), axis.title.x  = element_text(size=14),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_manual(values= c("#C9283E", "#0284A8", "#820333", "#02547D", "#FFBC67")) +
  geom_point(shape=1, size=6, colour="black")
  
#correlation matrix heatmap
library(RColorBrewer)
mycols <-(c("#C9283E", "#0284A8", "#820333", "#02547D", "#FFBC67"))
library(gplots)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          colRow=mycols[condition], colCol=mycols[condition],
          cexRow=1.2, cexCol=1.2,
          col=colorpanel(100, "dark blue", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(15, 15))
dev.off()

#get results with schrunken log fold change
reslfc1 <- lfcShrink(dds=ddds, contrast=c("condition","d28HIOENS","d28HIO"), alpha=0.05)
reslfc2 <- lfcShrink(dds=ddds, contrast=c("condition","d14HIOENS","d14HIO"), alpha=0.05)
reslfc3 <- lfcShrink(dds=ddds, contrast=c("condition","d28HIO","d14HIO"), alpha=0.05)
reslfc4 <- lfcShrink(dds=ddds, contrast=c("condition","d28HIOENS","d14HIOENS"), alpha=0.05)

#annotation HGCN symbol
library(org.Hs.eg.db)
reslfc1$symbol <- mapIds(org.Hs.eg.db,
  keys=row.names(reslfc1),
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals="first")
reslfc2$symbol <- mapIds(org.Hs.eg.db,
  keys=row.names(reslfc2),
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals="first")
  reslfc3$symbol <- mapIds(org.Hs.eg.db,
  keys=row.names(reslfc3),
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals="first")
reslfc4$symbol <- mapIds(org.Hs.eg.db,
  keys=row.names(reslfc4),
  column="SYMBOL",
  keytype="ENSEMBL",
  multiVals="first")
  
#significant genes
res_lfc_signi1 <- subset(reslfc1, padj < 0.05)
res_lfc_signi2 <- subset(reslfc2, padj < 0.05)
res_lfc_signi3 <- subset(reslfc3, padj < 0.05)
res_lfc_signi4 <- subset(reslfc4, padj < 0.05)

#get up and down gene lists
res_lfc_down1 <- subset(res_lfc_signi1,log2FoldChange< -1.5)
res_lfc_down_ord1 <- res_lfc_down1[ order(res_lfc_down1$log2FoldChange), ]
res_lfc_down2 <- subset(res_lfc_signi2,log2FoldChange< -1.5)
res_lfc_down_ord2 <- res_lfc_down2[ order(res_lfc_down2$log2FoldChange), ]
res_lfc_down3 <- subset(res_lfc_signi3,log2FoldChange< -1.5)
res_lfc_down_ord3 <- res_lfc_down3[ order(res_lfc_down3$log2FoldChange), ]
res_lfc_down4 <- subset(res_lfc_signi4,log2FoldChange< -1.5)
res_lfc_down_ord4 <- res_lfc_down4[ order(res_lfc_down4$log2FoldChange), ]

res_lfc_up1 <- subset(res_lfc_signi1,log2FoldChange> 1.5)
res_lfc_up_ord1 <- res_lfc_up1[ order(res_lfc_up1$log2FoldChange, decreasing = TRUE), ]
res_lfc_up2 <- subset(res_lfc_signi2,log2FoldChange> 1.5)
res_lfc_up_ord2 <- res_lfc_up2[ order(res_lfc_up2$log2FoldChange, decreasing = TRUE), ]
res_lfc_up3 <- subset(res_lfc_signi3,log2FoldChange> 1.5)
res_lfc_up_ord3 <- res_lfc_up3[ order(res_lfc_up3$log2FoldChange, decreasing = TRUE), ]
res_lfc_up4 <- subset(res_lfc_signi4,log2FoldChange> 1.5)
res_lfc_up_ord4 <- res_lfc_up4[ order(res_lfc_up4$log2FoldChange, decreasing = TRUE), ]

#MA plots
plotMA(reslfc1, ylim=c(-8,8))
plotMA(reslfc2, ylim=c(-4,4))
plotMA(reslfc3, ylim=c(-5,5))
plotMA(reslfc4, ylim=c(-6,6))

#Venn diagramms
library(VennDiagram)
library(gplots)
venn.plot1 <- venn.diagram(list(genes_down_lfc1_d28HIOENS_vs_d28HIO$X, genes_down_lfc2_d14HIOENS_vs_d14HIO$X),NULL,cat.fontfamily = "sans",lty="blank", fill=c("red", "cornflower blue "), alpha=c(0.5,0.5), cex = 1.5, cat.default.pos = "outer", margin = 0.25, category.names=c("log2FoldChange < -1.5 HIOENSj28/HIOj28", "log2FoldChange < -1.5 HIOENSj14/HIOj14"))
grid.draw(venn.plot1)

venn.plot2 <- venn.diagram(list(genes_up_lfc1_d28HIOENS_vs_d28HIO$X, genes_up_lfc2_d14HIOENS_vs_d14HIO$X), NULL,cat.fontfamily = "sans",lty="blank", fill=c("red", "cornflower blue "), alpha=c(0.5,0.5), cex = 1.5, cat.default.pos = "outer", margin = 0.25, category.names=c("log2FoldChange > 1.5 HIOENSj28/HIOj28", "log2FoldChange > 1.5 HIOj14/HIOj14"))
grid.draw(venn.plot2)

venn.plot3 <- venn.diagram(list(genes_down_lfc3_d28HIO_vs_d14HIO$X, genes_down_lfc4_d28HIOENS_vs_d14HIOENS$X), NULL,cat.fontfamily = "sans",lty="blank", fill=c("red", "cornflower blue "), cat.default.pos = "outer", margin = 0.25, alpha=c(0.5,0.5), cex = 1.5, category.names=c("log2FoldChange < -1.5 HIOj28/HIOj14", "log2FoldChange < -1.5 HIOENSj28/HIOENSj14"))
grid.draw(venn.plot3)

venn.plot4 <- venn.diagram(list(genes_up_lfc3_d28HIO_vs_d14HIO$X, genes_up_lfc4_d28HIOENS_vs_d14HIOENS$X), NULL,cat.fontfamily = "sans",lty="blank", fill=c("red", "cornflower blue "), cat.default.pos = "outer", margin = 0.25, alpha=c(0.5,0.5), cex = 1.5, category.names=c("log2FoldChange > 1.5 HIOj28/HIOj14", "log2FoldChange > 1.5 HIOENSj28/HIOENSj14"))
grid.draw(venn.plot4)

#gene lists from venn diagramms
gdwn_1 <- genes_down_lfc1_d28HIOENS_vs_d28HIO[,c(1,8)]
gdwn_2 <- genes_down_lfc2_d14HIOENS_vs_d14HIO[,c(1,8)]
gup_1 <- genes_up_lfc1_d28HIOENS_vs_d28HIO[,c(1,8)]
gup_2 <- genes_up_lfc2_d14HIOENS_vs_d14HIO[,c(1,8)]
gdwn_3 <- genes_down_lfc3_d28HIO_vs_d14HIO[,c(1,8)]
gdwn_4 <- genes_down_lfc4_d28HIOENS_vs_d14HIOENS[,c(1,8)]
gup_3 <- genes_up_lfc3_d28HIO_vs_d14HIO[,c(1,8)]
gup_4 <- genes_up_lfc4_d28HIOENS_vs_d14HIOENS[,c(1,8)]
gdwn_A <- merge(gdwn_1,gdwn_2)
gup_A <- merge(gup_1,gup_2)
gdwn_B <- merge(gdwn_3,gdwn_4)
gup_B <- merge(gup_3,gup_4)

a <- venn(list(genes_down_lfc1_d28HIOENS_vs_d28HIO$X, genes_down_lfc2_d14HIOENS_vs_d14HIO$X), show.plot=FALSE)
str(a)
inters <- attr(a,"intersections")
gdwn_1_uniq <- as.data.frame(inters$A)
colnames(gdwn_1_uniq) <- "X"
gdwn_1_unique <- merge(gdwn_1_uniq,gdwn_1)
gdwn_2_uniq <- as.data.frame(inters$B)
colnames(gdwn_2_uniq) <- "X"
gdwn_2_unique <- merge(gdwn_2_uniq,gdwn_2)
b <- venn(list(genes_up_lfc1_d28HIOENS_vs_d28HIO$X, genes_up_lfc2_d14HIOENS_vs_d14HIO$X), show.plot=FALSE)
intersb <- attr(b,"intersections")
gup_1_uniq <- as.data.frame(intersb$A)
colnames(gup_1_uniq) <- "X"
gup_1_unique <- merge(gup_1_uniq,gup_1)
gup_2_uniq <- as.data.frame(intersb$B)
colnames(gup_2_uniq) <- "X"
gup_2_unique <- merge(gup_2_uniq,gup_2)
c <-venn(list(genes_down_lfc3_d28HIO_vs_d14HIO$X, genes_down_lfc4_d28HIOENS_vs_d14HIOENS$X), show.plot=FALSE)
intersc <- attr(c,"intersections")
gdwn_3_uniq <- as.data.frame(intersc$A)
colnames(gdwn_3_uniq) <- "X"
gdwn_3_unique <- merge(gdwn_3_uniq,gdwn_3)
gdwn_4_uniq <- as.data.frame(intersc$B)
colnames(gdwn_4_uniq) <- "X"
gdwn_4_unique <- merge(gdwn_4_uniq,gdwn_4)
d <-venn(list(genes_up_lfc3_d28HIO_vs_d14HIO$X, genes_up_lfc4_d28HIOENS_vs_d14HIOENS$X), show.plot=FALSE)
intersd <- attr(d,"intersections")
gup_3_uniq <- as.data.frame(intersd$A)
colnames(gup_3_uniq) <- "X"
gup_3_unique <- merge(gup_3_uniq,gup_3)
gup_4_uniq <- as.data.frame(intersd$B)
colnames(gup_4_uniq) <- "X"
gup_4_unique <- merge(gup_4_uniq,gup_4)

#volcano plots
library(calibrate)
with(reslfc1, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot organoïdes à 28 jours avec SNE vs sans SNE", xlim=c(-11,11),ylim=c(-0.5,60)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(reslfc1, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(reslfc1, abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(reslfc1, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
with(subset(reslfc1, padj<.05 & abs(log2FoldChange)>3.5), textxy(log2FoldChange, -log10(pvalue), labs=symbol, cex=.7))

with(reslfc2, plot(log2FoldChange, -log10(pvalue), pch=20, main="volcano plot organoïdes à 14 jours avec SNE vs sans SNE", xlim=c(-5.5,5),ylim=c(-0.5,14)))
with(subset(reslfc2, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(reslfc2, abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(reslfc2, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
with(subset(reslfc2, padj<.05 & abs(log2FoldChange)>1.5), textxy(log2FoldChange, -log10(pvalue), labs=symbol, cex=.7))

with(reslfc3, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot d28HIO vs d14HIO", xlim=c(-15,15),ylim=c(-0.5,300)))
with(subset(reslfc3, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(reslfc3, abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(reslfc3, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
with(subset(reslfc3, padj<.05 & abs(log2FoldChange)>3), textxy(log2FoldChange, -log10(pvalue), labs=symbol, cex=.6))

with(reslfc4, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot d28HIOENS vs d14HIOENS", xlim=c(-12,12),ylim=c(-0.5,220)))
with(subset(reslfc4, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(reslfc4, abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(reslfc4, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
with(subset(reslfc4, padj<.05 & abs(log2FoldChange)>3), textxy(log2FoldChange, -log10(pvalue), labs=symbol, cex=.6))

#box plots
#log10 transformation gotten from pcaExplorer
log10 = read.csv("log10_pseudocount_added_pcaexplorer.txt", sep=",", head=T, row.names = 1, colClasses=c("character",rep("numeric",14)))

log10_T <- log10[c("ENSG00000164458"),]
test_log10_T <- as.data.frame(t(log10_T))
test_log10_T$sample <- c(rep("d28 HIO",3),rep("d28 HIOENS",3),rep("spheroid",3),rep("d14 HIO",3), rep("d14 HIOENS", 2))
test_log10_T$sample<-factor(test_log10_T$sample, levels=c("spheroid", "d14 HIO", "d14 HIOENS","d28 HIO", "d28 HIOENS"))
p <- ggplot(test_log10_T , aes(x=sample, y=ENSG00000164458))
p +  theme(axis.text.x = element_text(face="bold", color=mycols2, size=12, angle=0))+
  theme(plot.title = element_text(size = 16, face = "bold"))+
  scale_y_continuous(limits=c(0,3))+
  scale_x_discrete(name = "condition",labels=c("sphéroïdes","HIO j14","HIO+ENS j14","HIO j28","HIO+ENS j28"))+
  geom_boxplot(fill = mycols2) + geom_jitter(width = 0.2)+
  ylab("normalized counts -log10 scale") + ggtitle("T (ENSG00000164458)")
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0))
