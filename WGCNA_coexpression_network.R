counts_wgcna=read.csv("featcounts5.txt", sep="", head=T, skip=1, row.names = "Geneid")
counts_wgcna<-counts_wgcna[,6:41]
counts_wgcna<-as.matrix(counts_wgcna)

#normalise data
library(limma)
RNAseq = counts_wgcna[apply(counts_wgcna,1,function(x) sum(x==0))<ncol(counts_wgcna)*0.8,]
RNAseq_voom = voom(RNAseq)$E

library(flashClust)
library(WGCNA)
library(ape)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:5000],])
t_WGCNA_matrix <- t(WGCNA_matrix)
s = abs(bicor(WGCNA_matrix))

#pick soft power threshold
powers_1 = c(c(1:10), seq(from = 12, to=30, by=2))
sft_1 = pickSoftThreshold(WGCNA_matrix, powerVector = powers_1, verbose = 5, networkType ="signed", corFnc= "bicor")
plot(sft_1$fitIndices[,1], -sign(sft_1$fitIndices[,3])*sft_1$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'),
     ylim=c(-1,1));
text(sft_1$fitIndices[,1], -sign(sft_1$fitIndices[,3])*sft_1$fitIndices[,2],
     labels=powers_1,cex=1,col='red'); abline(h=0.90,col='red')

softpower = 30

#adjacency matrix
adj= adjacency(WGCNA_matrix,type = "signed", power = softpower);

#turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(WGCNA_matrix,networkType = "signed", TOMType = "signed", power = softpower)
#TOM2=TOMsimilarity(adj, TOMType = "signed")

plot(TOM)
 
colnames(TOM) <- colnames(WGCNA_matrix)
rownames(TOM) <- colnames(TOM)
dissTOM=1-TOM

geneTree_1 = flashClust(as.dist(dissTOM),method="average")
plot(geneTree_1, main = "Gene clustering on TOM-based dissimilarity", xlab="", sub="",cex=0.3)

minModuleSize = 20
dynamicMods = cutreeDynamic(dendro = geneTree_1,  method="tree", minClusterSize = minModuleSize)

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods) 
#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
plotDendroAndColors(geneTree_1, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

#set the diagonal of the dissimilarity to NA 
diag(dissTOM) = NA

#Visualize the Tom plot
TOMplot(dissTOM^4, geneTree_1, as.character(dynamicColors))


#module extraction
module_colors= setdiff(unique(dynamicColors), "grey")
for (color in module_colors){
  module=colnames(TOM)[which(dynamicColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

module.order <- unlist(tapply(1:ncol(WGCNA_matrix),as.factor(dynamicColors),I))
m<-t(t(WGCNA_matrix[,module.order])/apply(WGCNA_matrix[,module.order],2,max))

heatmap(t(m),zlim=c(0,1),col=blueWhiteRed(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=dynamicColors[module.order])

# Calculate eigengenes
MEList = moduleEigengenes(WGCNA_matrix, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
     
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(WGCNA_matrix, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

plotDendroAndColors(geneTree_1, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

txt_files <- list.files(pattern = "*.txt") 
for(i in txt_files) { x <- read.table(i, stringsAsFactors=F)

#get HGCN symbols
library(org.Hs.eg.db)
 x$symbol <- mapIds(org.Hs.eg.db,
                                keys=x$V1,
                                column="SYMBOL",
                                keytype="ENSEMBL",
                                multiVals="first")
  assign(i,x)
}

#eigengenes
for (color in moduleColors){
  module=colnames(TOM)[which(mergedColors==color)]
  write.table(module, paste("module_eigen_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

txt_files_eigen <- list.files(pattern = "*.txt") 
for(i in txt_files_eigen) { x <- read.table(i, stringsAsFactors=F)

x$symbol <- mapIds(org.Hs.eg.db,
                   keys=x$V1,
                   column="SYMBOL",
                   keytype="ENSEMBL",
                   multiVals="first")
  assign(i,x)
}


module.order.eigen <- unlist(tapply(1:ncol(WGCNA_matrix),as.factor(moduleColors),I))

#heatmap eigengens
meigen<-t(t(WGCNA_matrix[,module.order.eigen])/apply(WGCNA_matrix[,module.order.eigen],2,max))
heatmap(t(meigen),zlim=c(0,1),col=blueWhiteRed(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=moduleColors[module.order.eigen])
