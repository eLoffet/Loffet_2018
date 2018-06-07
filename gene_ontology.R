res_lfc_signi1 <- subset(reslfc1, padj < 0.05)
res_lfc1_signi_df <- as.data.frame(res_lfc1_signi)
reslfc1_pourGO <- na.omit(res_lfc_signi1_df)
write.csv(reslfc1_pourGO, "reslfc1_pourGO.csv")

###get GO tab from https://toppgene.cchmc.org/ with reslfc1_pourGO gene list

res_lfc1_df <- as.data.frame(reslfc1)
de_df <- res_lfc1_df[res_lfc1_df$padj < .05 & !is.na(res_lfc1_df$padj) & !is.na(res_lfc1_df$symbol),]
de_symbols <- de_df$symbol
# extract background genes
bg_ids <- rownames(ddds)[rowSums(counts(ddds)) > 0]
bg_symbols <- mapIds(org.Hs.eg.db,
                     keys=bg_ids,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

toppgene_tout <- read.csv("GOtoppgene_reslfc1_de_df_test.txt", sep="\t", head=T)

genes <- reslfc1_pourGO[,c("symbol","log2FoldChange")]
colnames(genes) <- c("ID","logFC") #est-ce qu'il faut changer log2FC en log FC ?
colnames(terms) <- c("category","ID","term","adj_pval","genes")
terms$category <- gsub("GO: Cellular Component", "CC", terms$category)
terms$category <- gsub("GO: Biological Process", "BP", terms$category)
terms$category <- gsub("GO: Molecular Function", "MF", terms$category)

library(GOplot)
circle_dat_lfc1up1 <- circle_dat(terms, genes)
circ <- circle_dat_lfc1up1
reduced_circ <- reduce_overlap(circ, overlap = 0.75) # Reduce redundant terms with a gene overlap >= 0.75...
GOBubble(reduced_circ, labels = 2.8, ID=FALSE, title = 'Bubble plot',display = 'multiple',bg.col = T )

IDs <- c('GO:0002053', 'GO:0060537', 'GO:0048333', 'GO:0072359','GO:0030198')
GOCircle(circ, nsub = IDs)
