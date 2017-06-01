


setwd("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\mapping_mOTU_5")

dataTable_map_motu <- read.table("mapping_5.bases.coverage.allMGs.filtered.0.5.COG0012.tsv2", sep="\t", header=TRUE, row.names=1)
newNames_1 <- colnames(dataTable_map_motu)
newNames2 <- as.vector(newNames_1[-grep("0045m", newNames_1)])
dataTable_map_motu <- dataTable_map_motu[,newNames2]


library(vegan)
bray_matrix_map_mOTU <- vegdist(t(dataTable_map_motu))
bray_matrix_map_mOTU_2 <- bray_matrix_map_mOTU
bray_matrix_map_mOTU_2[!is.finite(bray_matrix_map_mOTU_2)] <- 1.0
#save(bray_matrix_map_mOTU_2, file = "aloha_cat_bray_matrix_map.R")

library("ape")
#change names
samples2cruiseDepth <- read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\samples2cruiseDepth.txt", sep="\t", row.names=1)
bray_matrix_map_mOTU_mat <- as.matrix(bray_matrix_map_mOTU_2)
names = rownames(bray_matrix_map_mOTU_mat)

newNames =  samples2cruiseDepth[names,]
colnames(bray_matrix_map_mOTU_mat) = newNames
rownames(bray_matrix_map_mOTU_mat) = newNames


dataTable_motu_2 <- read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\samples2cruiseDepth_data.txt", sep="\t", row.names=1)
dataTable_motu_2 =  dataTable_motu_2[names,]

depthLabels <- dataTable_motu_2[2]
layerLabels <- dataTable_motu_2[3]

depth2Colors <- read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\depth2Colors.txt", sep="\t", row.names=1)
colorsDepth =  depth2Colors[as.vector(depthLabels[,1]),]
colorsDepth <-  as.vector(colorsDepth)

sampleOrderByDepth <- read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\sampleOrderByDepth.txt")
sampleOrderByDepth <- as.vector(t(sampleOrderByDepth))
sampleOrderByDepth <- as.vector(sampleOrderByDepth[-grep("0045m", sampleOrderByDepth)])

layer2number <- read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\layer2number.txt", sep="\t", row.names=1)
layerLabelsIndices <- as.vector(t(layerLabels))
layer2number_bySample <- as.numeric(as.vector(layer2number[layerLabelsIndices,]))

clustering_bray_matrix_map <- hclust(as.dist(bray_matrix_map_mOTU_mat))
plot(clustering_bray_matrix_map)

library(dendextend)
clustering_bray_dend <- as.dendrogram(clustering_bray_matrix_map)
clustering_bray_dend <- rotate(clustering_bray_dend, rev(sampleOrderByDepth))
#clustering_bray_dend <- rotate(clustering_bray_dend, sampleOrderByDepth)
clustering_bray_dend <- hang.dendrogram(clustering_bray_dend,hang_height=0.15)
clustering_bray_dend <- assign_values_to_branches_edgePar(clustering_bray_dend, value = 3, edgePar = "lwd")
clustering_bray_dend <- assign_values_to_branches_edgePar(clustering_bray_dend, value = "#202020", edgePar = "col")
	 
treeOrder <- order.dendrogram(clustering_bray_dend)
#layer2number_byTree <- layer2number_bySample[treeOrder]
#layer2Colors <- read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\layer2color_3.txt", sep="\t", row.names=1)
#colsBranches = as.vector(layer2Colors[c("Mixed","aboveDCM","SubDCM","200m","Meso"),])
#clustering_bray_dend <- color_branches(clustering_bray_dend, clusters=layer2number_byTree,col=colsBranches)

colsLeaves <- as.vector(depth2Colors[as.vector(depthLabels[,1])[treeOrder],])
clustering_bray_dend <- assign_values_to_leaves_edgePar(dend=clustering_bray_dend, value = colsLeaves, edgePar = "col")
labels_colors(clustering_bray_dend) <- "#00000000"

par(mar = c(12,3,3,3))
plot(clustering_bray_dend, main = "Clustered ALOHA metagenomes", horiz = FALSE,)

###################################
###################################
library("plotrix")

meta_data <- read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\HOT_2010_11_metadata_4_2.tab", sep="\t", header=TRUE, row.names=2)
HOT_2010_11_sampleNames <- colnames(bray_matrix_map_mOTU_mat)


meta_data_renamed <- meta_data[HOT_2010_11_sampleNames,]
pressure <- meta_data_renamed$press
pressure_cols <- color.scale(pressure)

depth <- meta_data_renamed$depth
depth_cols <- color.scale(depth)
depth <- log(depth)
depth_cols <- color.scale(depth, extremes=c("#ffff33","#377eb8"))

depth <- meta_data_renamed$depth_m
depth2Colors <- read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\depth2Colors.txt", sep="\t", row.names=1)
depth_cols <- depth2Colors[as.vector(depth),]
depth_cols <- as.vector(depth_cols)

surfaceLight <- meta_data_renamed$SW_5
surfaceLight <- (surfaceLight)
surfaceLight_cols <- color.scale(surfaceLight, extremes=c("#ffff33","#377eb8"))

month <- meta_data_renamed$month
month <- (month)
month_cols <- color.scale(month, extremes=c("#ffff33","#377eb8"))

year <- meta_data_renamed$year
year <- (year)
year_cols <- color.scale(year, extremes=c("#ffff33","#377eb8"))

julianDay <- meta_data_renamed$julian_day
julianDay <- (julianDay)
julianDay_cols <- color.scale(julianDay , extremes=c("#ffff33","#377eb8"))

DCM_position  <- meta_data_renamed$DCM_numeric
DCM_position <- (DCM_position)
DCM_position_cols <- color.scale(DCM_position , extremes=c("#ffff33","#377eb8"))

Mixed_layer_position  <- meta_data_renamed$Mixed_layer_num
Mixed_layer_position <- (Mixed_layer_position)
Mixed_layer_position_cols <- color.scale(Mixed_layer_position , extremes=c("#377eb8","#ffff33"))



oxygen <- meta_data_renamed$boxy
oxygen_cols <- color.scale(oxygen)
oxygen <- log(oxygen)
oxygen_cols <- color.scale(oxygen, extremes=c("#377eb8","#ffff33"))

salt <- meta_data_renamed$bsal
salt_cols <- color.scale(salt)
salt <- log(salt)
salt_cols <- color.scale(salt, extremes=c("#377eb8","#ffff33"))

phos <- meta_data_renamed$llp.phos
phos_cols <- color.scale(phos)
phos <- log(phos)
phos_cols <- color.scale(phos, extremes=c("#377eb8","#ffff33"))

nitrate <- meta_data_renamed$lln.nit
nitrate_cols <- color.scale(nitrate)
nitrate <- (log(nitrate))
nitrate_cols <- color.scale(nitrate, extremes=c("#377eb8","#ffff33"))

temp <- meta_data_renamed$temp
temp_cols <- color.scale(temp)
temp <- log(temp)
temp_cols <- color.scale(temp, extremes=c("#377eb8","#ffff33"))

GC_gene_cat <- meta_data_renamed$GC_gene_cat_95
GC_gene_cat_cols <- color.scale(GC_gene_cat)
GC_gene_cat <- (GC_gene_cat)
GC_gene_cat_cols <- color.scale(GC_gene_cat, extremes=c("#377eb8","#ffff33"))

N_ARSC_gene_cat <- meta_data_renamed$N.ARSC_gene_cat_95
N_ARSC_gene_cat_cols <- color.scale(N_ARSC_gene_cat)
N_ARSC_gene_cat <- (N_ARSC_gene_cat)
N_ARSC_gene_cat_cols <- color.scale(N_ARSC_gene_cat, extremes=c("#377eb8","#ffff33"))

layer <- meta_data_renamed$Layer_fineGrained
layer2Colors <- read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\layer2color_3.txt", sep="\t", row.names=1)
layer_cols =  layer2Colors[as.vector(layer),]
layer_cols <-  as.vector(layer_cols)


dataTable_genomeSizes <- read.table("C:\\work\\delong\\aloha_gene_catalogue\\genomeSizes\\all.microbe_census.out", sep="\t", row.names=1)
names <- rownames(dataTable_genomeSizes)
newNames = samples2cruiseDepth[names,1]
rownames(dataTable_genomeSizes) <- newNames
dataTable_genomeSizes_renamed <- dataTable_genomeSizes [HOT_2010_11_sampleNames,]
dataTable_genomeSizes_renamed_cols <- color.scale(dataTable_genomeSizes_renamed, extremes=c("#377eb8","#ffff33"))


#the_bars <- cbind(N_ARSC_gene_cat_cols, GC_gene_cat_cols, nitrate_cols, phos_cols,  oxygen_cols, salt_cols, temp_cols, depth_cols, layer_cols)
#colnames(the_bars) <- c("N-ARSC","GC","Oxygen","Salinity","Nitrate","Phosphate","Temperature","Depth","Community Cluster")
#the_bars_rev <- the_bars[,rev(colnames(the_bars))]

# surfaceLight_cols, month_cols, year_cols, julianDay_cols, DCM_position_cols
# "Surface Solar Irridiance", "Month", "Year", "Julian Day", "above/below DCM"

#the_bars <- cbind(layer_cols, DCM_position_cols, surfaceLight_cols, julianDay_cols, dataTable_genomeSizes_renamed_cols, N_ARSC_gene_cat_cols, GC_gene_cat_cols, phos_cols, nitrate_cols, oxygen_cols, salt_cols, temp_cols, depth_cols)
#colnames(the_bars) <- c("Clustering", "above/below DCM", "Surface Solar Irridiance", "Julian Day", "Avg. Genome Size","N-ARSC","GC","Phosphate","Nitrate","Oxygen","Salinity","Temperature","Depth")
#the_bars_rev <- the_bars[,rev(colnames(the_bars))]

#the_bars <- cbind(layer_cols,  N_ARSC_gene_cat_cols, GC_gene_cat_cols, nitrate_cols, phos_cols, oxygen_cols, salt_cols, temp_cols, Mixed_layer_position_cols, DCM_position_cols, depth_cols)
#colnames(the_bars) <- c("Clustering",  "N-ARSC", "GC","Nitrate","Phosphate","Oxygen","Salinity","Temperature","above/below Mixed Layer","above/below DCM","Depth")
#the_bars_rev <- the_bars[,rev(colnames(the_bars))]
the_bars <- cbind(layer_cols,  N_ARSC_gene_cat_cols, GC_gene_cat_cols, nitrate_cols, phos_cols, oxygen_cols, salt_cols, temp_cols, DCM_position_cols, depth_cols)
colnames(the_bars) <- c("Clustering",  "N-ARSC", "GC","Nitrate","Phosphate","Oxygen","Salinity","Temperature","above/below DCM","Depth")
the_bars_rev <- the_bars[,rev(colnames(the_bars))]

par(mar = c(3,6,3,28))
plot(clustering_bray_dend, main = "Clustered ALOHA metagenomes", horiz = TRUE,)
colored_bars(colors = the_bars_rev, dend = clustering_bray_dend, y_scale=1.0, y_shift=1.1, horiz = TRUE)


clustering_bray_dend <- assign_values_to_branches_edgePar(clustering_bray_dend, value = 2, edgePar = "lwd")
pdf(file = "figure1_1_top.pdf", 7, 5, bg = "transparent")
par(mar = c(3,6,3,18))
plot(clustering_bray_dend, main = "Clustered ALOHA metagenomes", horiz = TRUE,)
colored_bars(colors = the_bars_rev, dend = clustering_bray_dend, y_scale=0.8, y_shift=1.1, horiz = TRUE)
dev.off()
	 

	 
	 


#diversity boxplots

setwd("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\mapping_mOTU_5")
dataTable_map_motu <- read.table("mapping_5.bases.coverage.allMGs.filtered.0.5.COG0012.tsv2", sep="\t", header=TRUE, row.names=1)
depth2Colors  <- as.matrix(read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\depth2Colors.txt", sep="\t", row.names=1))

library("vegan")
library("ape")
library("reshape2")
bray_matrix_map_mOTU <- vegdist(t(dataTable_map_motu))
bray_matrix_map_mOTU_2 <- bray_matrix_map_mOTU
bray_matrix_map_mOTU_2[!is.finite(bray_matrix_map_mOTU_2)] <- 1.0

dataTable_map_motu_norm <- read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\samples2cruiseDepth.txt", sep="\t", row.names=1)
bray_matrix_map_mOTU_mat <- as.matrix(bray_matrix_map_mOTU_2)
names = rownames(bray_matrix_map_mOTU_mat)
newNames =  dataTable_map_motu_norm[names,]
colnames(bray_matrix_map_mOTU_mat) = newNames
rownames(bray_matrix_map_mOTU_mat) = newNames

dataTable_metadata <- read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\samples2cruiseDepth_data.txt", sep="\t", row.names=1)
dataTable_metadata =  dataTable_metadata[names,]

depthLabels <- dataTable_metadata[2]
layerLabels <- dataTable_metadata[3]


bray_matrix_map_mOTU_df <- melt(bray_matrix_map_mOTU_mat)
colnames(bray_matrix_map_mOTU_df) <- c("Sample1", "Sample2", "Distance")
bray_matrix_map_mOTU_df$depth1 <- dataTable_metadata[bray_matrix_map_mOTU_df$Sample1,2]
bray_matrix_map_mOTU_df$depth2 <- dataTable_metadata[bray_matrix_map_mOTU_df$Sample2,2]
	 
bray_matrix_map_mOTU_df_filter1 <- bray_matrix_map_mOTU_df[bray_matrix_map_mOTU_df$depth1==bray_matrix_map_mOTU_df$depth2,]
bray_matrix_map_mOTU_df_filter1_1 <- bray_matrix_map_mOTU_df_filter1[bray_matrix_map_mOTU_df_filter1$Sample1!=bray_matrix_map_mOTU_df_filter1$Sample2,]
bray_matrix_map_mOTU_df_filter1_1 <- bray_matrix_map_mOTU_df_filter1_1[which(bray_matrix_map_mOTU_df_filter1_1$Sample1 != "HOT224-0125m" & bray_matrix_map_mOTU_df_filter1_1$Sample1 != "HOT227-0125m" & bray_matrix_map_mOTU_df_filter1_1$Sample1 != "HOT233-0125m"  & bray_matrix_map_mOTU_df_filter1_1$Sample2 != "HOT224-0125m" & bray_matrix_map_mOTU_df_filter1_1$Sample2 != "HOT227-0125m" & bray_matrix_map_mOTU_df_filter1_1$Sample2 != "HOT233-0125m"),]

rows <- c(which(bray_matrix_map_mOTU_df$depth1=="25m" & bray_matrix_map_mOTU_df$depth2=="75m"), which(bray_matrix_map_mOTU_df$depth1=="75m" & bray_matrix_map_mOTU_df$depth2=="125m") , which(bray_matrix_map_mOTU_df$depth1=="125m" & bray_matrix_map_mOTU_df$depth2=="200m"),which(bray_matrix_map_mOTU_df$depth1=="200m" & bray_matrix_map_mOTU_df$depth2=="500m"), which(bray_matrix_map_mOTU_df$depth1=="500m" & bray_matrix_map_mOTU_df$depth2=="770m"), which(bray_matrix_map_mOTU_df$depth1=="770m" & bray_matrix_map_mOTU_df$depth2=="1000m"))

bray_matrix_map_mOTU_df_filter2 <-  bray_matrix_map_mOTU_df[rows,]

bray_matrix_map_mOTU_df_filter3 <- bray_matrix_map_mOTU_df_filter2

bray_matrix_map_mOTU_df_filter3 <- bray_matrix_map_mOTU_df_filter3[c(which(bray_matrix_map_mOTU_df_filter3$Sample2 == "HOT224-0125m"), which(bray_matrix_map_mOTU_df_filter3$Sample2 == "HOT227-0125m"), which(bray_matrix_map_mOTU_df_filter3$Sample2 == "HOT233-0125m")),]

#bray_matrix_map_mOTU_df_filter3 <- bray_matrix_map_mOTU_df_filter3[c(which(bray_matrix_map_mOTU_df_filter3$Sample1 == "HOT224-0125m"), which(bray_matrix_map_mOTU_df_filter3$Sample1 == "HOT227-0125m"), which(bray_matrix_map_mOTU_df_filter3$Sample1 == "HOT233-0125m"),  which(bray_matrix_map_mOTU_df_filter3$Sample2 == "HOT224-0125m"), which( bray_matrix_map_mOTU_df_filter3$Sample2 == "HOT227-0125m"), which(bray_matrix_map_mOTU_df_filter3$Sample2 == "HOT233-0125m")),]

bray_matrix_map_mOTU_df_filter2 <-  bray_matrix_map_mOTU_df_filter2[which(bray_matrix_map_mOTU_df_filter2$Sample1 != "HOT224-0125m" & bray_matrix_map_mOTU_df_filter2$Sample1 != "HOT227-0125m" & bray_matrix_map_mOTU_df_filter2$Sample1 != "HOT233-0125m"  & bray_matrix_map_mOTU_df_filter2$Sample2 != "HOT224-0125m" & bray_matrix_map_mOTU_df_filter2$Sample2 != "HOT227-0125m" & bray_matrix_map_mOTU_df_filter2$Sample2 != "HOT233-0125m"),]


library("ggplot2")

bray_matrix_map_mOTU_df_filter1_1$depth  <- (factor((bray_matrix_map_mOTU_df_filter1_1$depth1), levels=c("1000m", "770m", "500m", "200m", "125m", "75m", "45m", "25m")))

#depth_list = rev(c("25m", "45m", "75m", "125m", "200m", "500m", "770m", "1000m"))
#colsDepth = as.vector(depth2Colors[depth_list, 1])
greyCols = c("#202020","#202020","#202020","#202020","#202020","#202020","#202020","#202020")

depth_list_1 = rev(c("25m",  "75m", "125m", "200m", "500m", "770m"))
colsDepth_1 = as.vector(depth2Colors[depth_list_1, 1])
depth_list_2 = rev(c("75m", "125m", "200m", "500m", "770m", "1000m"))
colsDepth_2 = as.vector(depth2Colors[depth_list_2, 1])

bray_matrix_map_mOTU_df_filter2$depth  <- (factor((bray_matrix_map_mOTU_df_filter2$depth1), levels=c("770m", "500m", "200m", "125m", "75m", "25m")))

bray_matrix_map_mOTU_df_filter3$depth  <- (factor((bray_matrix_map_mOTU_df_filter3$depth1), levels=c("770m", "500m", "200m", "125m", "75m", "25m")))
color3 = c("#fb8072","#fb8072","#fb8072","#fb8072","#fb8072","#fb8072","#fb8072","#fb8072")

xlabels = c("770m vs 1000m","500m vs 770m","200m vs 500m","125m vs 200m","75m vs 125m","25m vs 75m")

p <- ggplot(bray_matrix_map_mOTU_df_filter2, aes(x=depth, y=Distance, fill=depth, color=depth )) + scale_color_manual(values=c(colsDepth_2)) + scale_fill_manual(values=c(colsDepth_1)) + geom_boxplot(outlier.shape = NA, lwd=3, fatten =0.75) + geom_jitter(position=position_jitter(0.3), cex=5, shape=21) + labs(x="Depth",y="Distance") 
p2 <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_flip() + scale_x_discrete("Depth", waiver(), xlabels)
p2

p3 <- p2 + geom_jitter(data=bray_matrix_map_mOTU_df_filter3,  position=position_jitter(0.3), cex=5, shape=24, fill = "#e41a1c", colour = "#ff7f00")  
p3

p3 <- p2 + geom_jitter(data=bray_matrix_map_mOTU_df_filter3,  position=position_jitter(0.3), cex=5, shape=24, fill = "#d0d0d0", colour = "#000000") 
p3

ggsave("mOTU_intraDepth_brayCurtis.pdf", bg = "transparent",  width = 10, height = 8)
ggsave("mOTU_intraDepth_brayCurtis_bg.pdf", bg = "white",  width = 10, height = 8)
ggsave("mOTU_intraDepth_brayCurtis.svg", bg = "transparent",  width = 10, height = 8)
ggsave("mOTU_intraDepth_brayCurtis_bg.svg", bg = "white",  width = 10, height = 8)




bray_matrix_map_mOTU_df_filter2$depth  <- (factor((bray_matrix_map_mOTU_df_filter2$depth1), levels=c("770m", "500m", "200m", "125m", "75m", "25m")))

depth_list_1 = rev(c("25m",  "75m", "125m", "200m", "500m", "770m"))
colsDepth_1 = as.vector(depth2Colors[depth_list_1, 1])

depth_list_2 = rev(c("75m", "125m", "200m", "500m", "770m", "1000m"))
colsDepth_2 = as.vector(depth2Colors[depth_list_2, 1])
greyCols = c("#202020","#202020","#202020","#202020","#202020","#202020","#202020","#202020")

xlabels = c("770m vs 1000m","500m vs 770m","200m vs 500m","125m vs 200m","75m vs 125m","25m vs 75m")

p <- ggplot(bray_matrix_map_mOTU_df_filter2, aes(x=depth, y=Distance, fill=depth, color=depth )) + scale_color_manual(values=c(colsDepth_2)) + scale_fill_manual(values=c(colsDepth_1)) + geom_boxplot(outlier.shape = NA, lwd=3, fatten =0.5) + geom_jitter(position=position_jitter(0.3), cex=5, shape=21) + labs(x="Depth",y="Distance (BC)") 
p2 <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_flip() + scale_x_discrete("Depth", waiver(), xlabels)
p2

ggsave("mOTU_betweenDepth_brayCurtis.pdf", bg = "transparent",  width = 10, height = 8)
ggsave("mOTU_betweenDepth_brayCurtis_bg.pdf", bg = "white",  width = 10, height = 8)




bray_matrix_map_mOTU_df_filter2$depth  <- (factor((bray_matrix_map_mOTU_df_filter2$depth1), levels=c("770m", "500m", "200m", "125m", "75m", "25m")))

depth_list_1 = rev(c("25m",  "75m", "125m", "200m", "500m", "770m"))
colsDepth_1 = as.vector(depth2Colors[depth_list_1, 1])

depth_list_2 = rev(c("75m", "125m", "200m", "500m", "770m", "1000m"))
colsDepth_2 = as.vector(depth2Colors[depth_list_2, 1])
greyCols = c("#202020","#202020","#202020","#202020","#202020","#202020","#202020","#202020")

xlabels = c("770m vs 1000m","500m vs 770m","200m vs 500m","125m vs 200m","75m vs 125m","25m vs 75m")

p <- ggplot(bray_matrix_map_mOTU_df_filter2, aes(x=depth, y=Distance, fill=depth, color=depth )) + scale_color_manual(values=c(colsDepth_2)) + scale_fill_manual(values=c(colsDepth_1)) + geom_boxplot(outlier.shape = NA, lwd=2, fatten =0.75) + geom_jitter(position=position_jitter(0.3), cex=5, shape=21) + labs(x="Depth",y="Distance (BC)") 
p2 <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_flip() + scale_x_discrete("Depth", waiver(), xlabels)
p2

ggsave("mOTU_betweenDepth_brayCurtis_v2.pdf", bg = "transparent",  width = 10, height = 8)
ggsave("mOTU_betweenDepth_brayCurtis_v2_bg.pdf", bg = "white",  width = 10, height = 8)


p3 <- p2 + geom_jitter(data=bray_matrix_map_mOTU_df_filter3,  position=position_jitter(0.3), cex=5, shape=21, fill = "#fb8072", colour = "#fb8072") 
p3

ggsave("mOTU_betweenDepth_brayCurtis_v2_plus.pdf", bg = "transparent",  width = 10, height = 8)
ggsave("mOTU_betweenDepth_brayCurtis_v2_plus_bg.pdf", bg = "white",  width = 10, height = 8)




setwd("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\mapping_mOTU_5")

library("ggplot2")

bray_matrix_map_mOTU_df_filter1_1$depth  <- (factor((bray_matrix_map_mOTU_df_filter1_1$depth1), levels=c("1000m", "770m", "500m", "200m", "125m", "75m", "45m", "25m")))

depth_list = rev(c("25m", "45m", "75m", "125m", "200m", "500m", "770m", "1000m"))
colsDepth = as.vector(depth2Colors[depth_list, 1])
greyCols = c("#202020","#202020","#202020","#202020","#202020","#202020","#202020","#202020")

bray_matrix_map_mOTU_df_filter3$depth  <- (factor((bray_matrix_map_mOTU_df_filter3$depth1), levels=c("1000m", "770m", "500m", "200m", "125m", "75m", "45m", "25m")))
color3 = c("#fb8072","#fb8072","#fb8072","#fb8072","#fb8072","#fb8072","#fb8072","#fb8072")


p <- ggplot(bray_matrix_map_mOTU_df_filter1_1, aes(x=depth, y=Distance, fill=depth, color=depth )) + scale_color_manual(values=c(greyCols)) + scale_fill_manual(values=c(colsDepth)) + geom_boxplot(outlier.shape = NA, lwd=2, fatten =0.75) + geom_jitter(position=position_jitter(0.3), cex=5, shape=21) + labs(x="Depth",y="Distance (BC)") 
p2 <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_flip()
p2


#no 45m
setwd("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\mapping_mOTU_5")

library("ggplot2")

bray_matrix_map_mOTU_df_filter1_2 <- bray_matrix_map_mOTU_df_filter1_1[bray_matrix_map_mOTU_df_filter1_1$depth1!="45m",]

bray_matrix_map_mOTU_df_filter1_2$depth  <- (factor((bray_matrix_map_mOTU_df_filter1_2$depth1), levels=c("1000m", "770m", "500m", "200m", "125m", "75m", "25m")))

depth_list = rev(c("25m",  "75m", "125m", "200m", "500m", "770m", "1000m"))
colsDepth = as.vector(depth2Colors[depth_list, 1])
greyCols = c("#202020","#202020","#202020","#202020","#202020","#202020","#202020","#202020")

bray_matrix_map_mOTU_df_filter3$depth  <- (factor((bray_matrix_map_mOTU_df_filter3$depth1), levels=c("1000m", "770m", "500m", "200m", "125m", "75m",  "25m")))
color3 = c("#fb8072","#fb8072","#fb8072","#fb8072","#fb8072","#fb8072","#fb8072","#fb8072")


p <- ggplot(bray_matrix_map_mOTU_df_filter1_2, aes(x=depth, y=Distance, fill=depth, color=depth )) + scale_color_manual(values=c(greyCols)) + scale_fill_manual(values=c(colsDepth)) + geom_boxplot(outlier.shape = NA, lwd=1, fatten =1.5) + geom_jitter(position=position_jitter(0.3), cex=5, shape=21) + labs(x="Depth",y="Distance (BC)") 
p2 <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_flip()
p2

ggsave("mOTU_intraDepth_brayCurtis.pdf", bg = "transparent",  width = 10, height = 8)
ggsave("mOTU_intraDepth_brayCurtis_bg.pdf", bg = "white",  width = 10, height = 8)
