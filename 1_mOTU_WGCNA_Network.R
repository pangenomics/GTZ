

############################## Load Necessary Packages
library(WGCNA)
library(igraph)
library(DESeq2)

############################ Enable multithreading for WGCNA
enableWGCNAThreads()

######################################################################################
############ Load mOTU data and normalize by total abundance in a given sample #######
######################################################################################
motu <- read.table(file="mapping_3.bases.coverage.allMGs.filtered.0.5.COG0012.tsv", sep="\t", header=TRUE, row.names=1)
sums <- apply(motu, 2, sum)
motu2 <- scale(motu, scale=sums, center=F)

write.table(motu2, file="motu_abundances_sample_normalized.txt", sep="\t", quote=F)
#rowsums <- rowSums(motu)
#motu3 <- t(scale(t(motu), scale=rowsums, center=F))

## For Pearson's Correlation and clustering of motus
corr <- cor(t(motu2), method="pearson")
dist <- as.dist(1-corr)
clust <- stats::hclust(dist, method="complete")
#tree <- cutreeDynamic(clust, method="tree", cutHeight=1.5, minClusterSize=2)
tree <- cutree(clust, h=0.75)
color_tree <- labels2colors(tree)

# Cluster columns, in this case samples
dist <- as.dist(1 - cor(motu))
clust2 <- stats::hclust(dist, method="complete")
tree2 <- cutree(clust2, h=0.8)
color_treeR <- labels2colors(tree2)

#Plot sample dendrogram
jpeg(file="samples_cluster.jpg", quality=75, res=150, height=8, width=16, units="in")
plotDendroAndColors(clust2, color_treeR, hang = 0.05, addGuide = TRUE, guideHang = 0.05, cex.dendroLabels=0.6)
dev.off()

#Plot motu dendrogram
jpeg(file="motu_dendrogram.jpg", quality=75, res=150, height=20, width=50, units="in")
plotDendroAndColors(clust, color_tree, hang = 0.05, addGuide = TRUE, guideHang = 0.05, cex.dendroLabels=0.6)
dev.off()

# Plot heatmap with clustering dendrogram on top
library(RColorBrewer)
library(gplots)
colors <- colorRampPalette(brewer.pal(9,"Blues"))(200)
jpeg(file="motu_heatmap_dendrogram.jpg", quality=100, res=600, height=16, width=12, units="in")
heatmap.2(as.matrix(motu2), trace="none", Rowv = as.dendrogram(clust), dendrogram="both", Colv=as.dendrogram(clust2), scale="row", na.color="grey", col=colors, key=F, keysize=0.5, cexCol = 1, cexRow = 0.2, RowSideColors=color_tree, ColSideColors=color_treeR, labRow="", margins = c(9, 3), colRow=2)
dev.off()

#####################
######## End 
#####################


#################################################################
###################### Network ##################################
#################################################################

library(WGCNA)
###########################################################
####### Section 1: Pick thresholds and generate modules
###########################################################
# load data and normalize
motu <- read.table(file="mapping_3.bases.coverage.allMGs.filtered.0.5.COG0012.tsv", sep="\t", header=TRUE, row.names=1)
sums <- apply(motu, 2, sum)
motu2 <- scale(motu, scale=sums, center=F)

#Pick Soft Threshold
powers =c(c(1:10),seq(from = 12, to=20,by=2))
sft=pickSoftThreshold(t(motu2), powerVector=powers, verbose=5)

#Document Scale-free topology fits and power-connectivity relationship, output results as a figure
pdf(file="motu_SFT_Fit_DESeq_Final.pdf", height=5, width=9)
par(mfrow =c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signedR^2",type="n",main =paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main =paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()

##Power of 6 chosen based on criteria of finding lowest value for which R^2 > 0.8 in scale-free network test
pow <- 6
net = blockwiseModules(t(motu2), power= pow,TOMType ="signed", minModuleSize = 2, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, verbose = 3) 

#Get unique colors to associate with each module identified.
mergedColors = labels2colors(net$colors)  

#Combine module information (ortholog cluster ID, module number, module color) and output table
mod <- cbind(row.names(motu2), net$colors, mergedColors); colnames(mod) <- c("Cluster", "Module", "Module_Color")
write.table(mod, file="motu.modules", quote=F, sep="\t", row.names=F)

#Plot dendrogram
jpeg(file="combined.jpg", quality=75, res=150, height=15, width=30, units="in")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],dendroLabels = FALSE, hang = 0.05,addGuide = TRUE, guideHang = 0.05)
dev.off()

#calculate and output Topological Overlap Matrix (TOM, for later analyses)
TOM = TOMsimilarityFromExpr(t(motu2), power= pow)
colnames(TOM) = row.names(motu2)
write.table(TOM, file="motu.TOM", sep="\t")

###########################################################
####### Section 2: Calculate Eigengenes of Modules ########
###########################################################
modules <- names(sort(table(mergedColors), decreasing=T)[1:30])
Module_Colors <- read.table(file="motu.modules", sep="\t", header=TRUE, row.names=1)
mod <- data.frame(sort(table(Module_Colors[,2]), decreasing=T))
mergedColors <- Module_Colors[,2]

motu_ordered <- motu2[,order(colnames(motu2))]

depths = list()
splitter <- strsplit(colnames(motu2), "_")
for(i in 1:83) {
	tail <- tail(splitter[[i]], n=1)
	depths[i] <- tail
}
depths <- as.character(depths)
order <- order(depths, decreasing=T)

eig <- moduleEigengenes(t(motu2[,order]), as.character(Module_Colors[,2]), excludeGrey=TRUE)
data <- eig$eigengenes
colnames(data) <- gsub("ME", "", colnames(data))
data1 <- eig$averageExpr

jpeg(file="motu_eigengenes.jpg", quality=100, res=300, height=8, width=7, units="in")
par(mfrow=c(6,2), mar=c(0.5, 0.5, 1, 0.5), oma=c(2,0,1,0))
for(i in 1:12) {
#for(i in 1:(dim(data)[2]-1)) {
	color <- modules[i]
	axis <- c(1:83)
	#row.names(ag) <- axis
	title <- paste(c("Module ", i, "; ", color, "; ", "n=", mod[modules[i],]), collapse="")
	title <- paste(c("Module ", i, "; ", "n=", mod[modules[i],]), collapse="")
	print(title)
	plot(axis, data[,modules[i]], col=color, type="l", ylab=NA, lwd=2, main=title, cex.main=1.2, axes=F)
	#rect(xleft = c(0, 5.5, 11.5, 16.5, 23.5, 31.5), xright=c(2.5, 8.5, 13.5, 20.5, 27.5, 34.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(1, 1, 1, 1, 1, 1), col='grey90', border=NA)
	par(new=T)
	plot(axis, data[,modules[i]], col='steelblue', type="l", ylab=NA, lwd=2, main=title, cex.main=1.2, axes=F)
	box(lty=1, col="grey60")
}
dev.off()


jpeg(file="motu_eigengenes.jpg", quality=100, res=300, height=10, width=7, units="in")
par(mfrow=c(4,6), mar=c(0.5, 0.5, 1, 0.5), oma=c(2,3,1,0))
for(i in 1:24) {
#for(i in 1:(dim(data)[2]-1)) {
	color <- modules[i]
	axis <- c(1:83)
	#row.names(ag) <- axis
	title <- paste(c("Module ", i, "; ", color, "; ", "n=", mod[modules[i],]), collapse="")
	title <- paste(c("Module ", i, "; ", "n=", mod[modules[i],]), collapse="")
	print(title)
	plot(data[,modules[i]], axis, col=color, type="l", ylab=NA, lwd=2, main=title, cex.main=1.2, axes=F)
	#if(i %in% c(1, 7, 13, 19)) { axis(2, at = seq(from=1, to=83, by=1), labels = colnames(motu2)[order], cex.axis=0.5, las=2) }
	if(i %in% c(1, 7, 13, 19)) { colorbar.plot(-0.5, 40, c(1:83), horizontal=F) }
	par(new=T)
	plot(data[,modules[i]], axis, col='steelblue', type="l", ylab=NA, lwd=2, main=title, cex.main=1.2, axes=F)
	box(lty=1, col="grey60")

}
dev.off()


###########################################################
####### Section 3: Visualize Networks
###########################################################

########################################################################################
############################# Construct Networks using igraph package ##################
########################################################################################
############################## Load Necessary Packages
library(WGCNA)
library(igraph)
library(DESeq2)

############################ Enable multithreading for WGCNA
enableWGCNAThreads()

# load data from previous sections
TOM <- read.table("motu.TOM", header=T, row.names=1, sep="\t")
motu2 <- read.table("motu_abundances_sample_normalized.txt", sep="\t", header=T, row.names=1)
Module_Colors <- read.table(file="motu.modules", sep="\t", header=TRUE, row.names=1)

# initialize the datasets
modules <- names(sort(table(Module_Colors$Module_Color), decreasing=T))[1:20] 
row.names(TOM) <- colnames(TOM)
nodes <- row.names(motu2)[which(Module_Colors[,2] %in% modules)]
modTOM <- TOM[nodes, nodes]

# choose cutoffs for what edges to show in the network
#tom_w_signs <- TOM * sign(cor(t(motu2)))
#tomvals <- apply(tom_w_signs, 2, as.numeric)
#pos_cutoff <- mean(tomvals) + 2*sd(tomvals)
#neg_cutoff <- mean(tomvals) - 2*sd(tomvals)
#to_keep <- (tom_w_signs > pos_cutoff) + (tom_w_signs < neg_cutoff)
#tomfinal <- as.matrix(tom_w_signs * to_keep - diag(dim(to_keep)[1]))
#g  <- graph.adjacency(tomfinal*sign(tomfinal), weighted=T)
#df <- get.data.frame(g)

### Choose appropriate threshold for network representation
quant = quantile(apply(TOM, 2, as.numeric), 0.98)

cyt = exportNetworkToCytoscape(modTOM, weighted= TRUE, threshold = quant, nodeNames = nodes, nodeAttr = as.character(Module_Colors[nodes,2]))
#cyt = exportNetworkToCytoscape(modTOM,weighted= TRUE, threshold = quant, nodeNames = modProbes, nodeAttr = mergedColors[inModule])

newedges=(cyt$edgeData[,3])
edges <- cbind(as.character(cyt$edgeData[,1]), as.character(cyt$edgeData[,2]), cyt$edgeData[,3])
write.table(edges, file="motu.edges",  row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

######### Make Fruchterman-Reingold layout for graph
igr=read.graph("motu.edges", format="ncol", weights="yes")
set.seed(123)
l <- layout_with_kk(igr)

#########################################################################################
############# Color for cluster-based network ###########################################
#########################################################################################
nodes <- cbind(as.character(cyt$nodeData[,1]), as.character(cyt$nodeData[,3]))
write.table(nodes, file="motu.nodes", row.names=FALSE, col.names=c("node", "attr"), sep="\t", quote=FALSE)

x=read.table(file="motu.nodes", sep="\t", header=TRUE)
V(igr)$attr=as.character(x$attr[match(V(igr)$name, x$node)])
#V(igr)$color="grey75"
V(igr)$color=V(igr)$attr

mod_names <- c("Mod 1", "Mod 2", "Mod 3", "Mod 4", "Mod 5", "Mod 6", "Mod 7", "Mod 8", "Mod 9", "Mod 10", "Mod 11", "Mod 12", "Mod 13", "Mod 14", "Mod 15", "Mod 16", "Mod 17", "Mod 18", "Mod 19", "Mod 20")
#tiff(file="clusters.tiff", quality=75, res=100, height=30, width=30, units="in")
jpeg(file="motu_module_Network.0.99.jpg", quality=100, res=600, height=8, width=8, units="in")
#postscript(file="canon_module_Network.3.eps", height=2, width=2)
plot.igraph(igr,vertex.label=NA,layout=l, vertex.size=1.3, edge.color="grey85", vertex.frame.color=NA, edge.width=0.3)
legend("topleft", legend=mod_names, col=modules, pch=15, cex=0.7, bty="n")
dev.off()


#########################################################################################
############# Color by depth ###########################################
#########################################################################################

depth_file <- read.table(file="sample_depths.txt", sep="\t", header=T, row.names=1, stringsAsFactors=FALSE)
depth_file$depth <- as.character(depth_file$depth)
depths = as.character(unique(depth_file$depth))

motu <- read.table(file="mapping_3.bases.coverage.allMGs.filtered.0.5.COG0012.tsv", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)

depth_mat <- data.frame(t(depth_file))
data <- data.frame(motu[,row.names(depth_file)])

datalist <- list()
for(i in 1:length(depths)) {
	samples <- row.names(depth_file)[which(depth_file$depth == depths[i])]
	datalist[[i]] <- rowSums(motu[,samples])
}
data = do.call(cbind, datalist)
colnames(data) <- depths

most_abund <- list()
for(i in 1:dim(data)[1]) {
	d <- names(which.max(data[i,]))
	most_abund[i] <- d
}
most_abund <- as.character(most_abund)

node_info <- data.frame(cbind(as.character(row.names(data)), as.character(most_abund)))
#colnames(node_info) <- c("node", "depth")
node_info <- node_info[order(most_abund),]

############### Set color palette using RColorBrewer
library(RColorBrewer)
mypal <- c("green4", brewer.pal(9, "Set1")); mypal[7] = "gold"; mypal[4] <- "#54D126"

attr <- labels2colors(node_info[,2], colorSeq=mypal)  #  c("lightblue", "dodgerblue", "gold", "green4", "darkblue", "purple", "black"))
write.table(cbind(node_info, attr), file="motu.depth.nodes", row.names=FALSE, col.names=c("node", "depth", "attr"), sep="\t", quote=FALSE)

# initialize the datasets
modules <- names(sort(table(Module_Colors$Module_Color), decreasing=T))[1:20] 
row.names(TOM) <- colnames(TOM)
nodes <- row.names(motu2)[which(Module_Colors[,2] %in% modules)]
modTOM <- TOM[nodes, nodes]

### Choose appropriate threshold for network representation
quant = quantile(apply(TOM, 2, as.numeric), 0.9799)

cyt = exportNetworkToCytoscape(modTOM, weighted= TRUE, threshold = quant, nodeNames = nodes, nodeAttr = as.character(Module_Colors[nodes,2]))
#cyt = exportNetworkToCytoscape(modTOM,weighted= TRUE, threshold = quant, nodeNames = modProbes, nodeAttr = mergedColors[inModule])

newedges=(cyt$edgeData[,3])
edges <- cbind(as.character(cyt$edgeData[,1]), as.character(cyt$edgeData[,2]), cyt$edgeData[,3])
write.table(edges, file="motu.edges",  row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

######### Make Fruchterman-Reingold layout for graph
igr=read.graph("motu.edges", format="ncol", weights="yes")
set.seed(123)
l <- layout_with_kk(igr)

#########################################################################################
############# Color for cluster-based network ###########################################
#########################################################################################

x=read.table(file="motu.depth.nodes", sep="\t", header=TRUE, comment.char="")
V(igr)$attr=as.character(x$attr[match(V(igr)$name, x$node)])
#V(igr)$color="grey75"
V(igr)$color=V(igr)$attr

mod_names <- c("25m", "75m", "125m", "200m", "500m", "770m", "1000m")
jpeg(file="motu_module_Network.0.98_bydepth.jpg", quality=100, res=600, height=8, width=8, units="in")
plot.igraph(igr,vertex.label=NA,layout=l, vertex.size=1.5, edge.color="grey85", vertex.frame.color=NA, edge.width=0.3)
legend("topright", legend=mod_names, col=mypal, pch=15, cex=1.2, bty="n")
dev.off()






