
setwd("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\mapping_2_kegg\\nitrate_metabolism")


dataTable_map_kegg <- read.table("KO_nitrate_metabolism.combined.bases.coverage.mOTUnorm.tsv", sep="\t", header=TRUE, row.names=1)
newNames_1 <- colnames(dataTable_map_kegg)
newNames2 <- as.vector(newNames_1[-grep("0045m", newNames_1)])
dataTable_map_kegg <- dataTable_map_kegg[,newNames2]

dataTable_map_kegg_2 <- dataTable_map_kegg

library("ape")
#change names
samples2cruiseDepth <- read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\samples2cruiseDepth.txt", sep="\t", row.names=1)
names = colnames(dataTable_map_kegg)

newNames =  samples2cruiseDepth[names,]
colnames(dataTable_map_kegg_2) = newNames

dataTable_samples2cruise_2 <- read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\samples2cruiseDepth_data.txt", sep="\t", row.names=1)
dataTable_samples2cruise_2 =  dataTable_samples2cruise_2[names,]

depthLabels <- dataTable_samples2cruise_2[2]
layerLabels <- dataTable_samples2cruise_2[3]

rownames(depthLabels) <- dataTable_samples2cruise_2[,1]
rownames(layerLabels) <- dataTable_samples2cruise_2[,1]


depth2Colors  <- as.matrix(read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\depth2Colors.txt", sep="\t", row.names=1))
library("ggplot2")
#depth_list = rev(c("25m", "75m", "125m", "200m", "500m", "770m", "1000m"))
depth_list = (c("25m", "75m", "125m", "200m", "500m", "770m", "1000m"))
colsDepth = as.vector(depth2Colors[depth_list, 1])

library(reshape2)
library(ggplot2)
df_dataTable_map_kegg_2 <- melt(t(dataTable_map_kegg_2))
colnames(df_dataTable_map_kegg_2) <- c("sample","KO","abundance")

df_dataTable_map_kegg_2$depth <- as.vector(depthLabels[df_dataTable_map_kegg_2$sample,])
df_dataTable_map_kegg_2$depth  <- (factor((df_dataTable_map_kegg_2$depth), levels=rev(c("1000m", "770m", "500m", "200m", "125m", "75m", "45m", "25m"))))

df_dataTable_map_kegg_2_melt <- melt(df_dataTable_map_kegg_2, id = c( "sample","KO", "depth"))
df_dataTable_map_kegg_3 <- dcast(df_dataTable_map_kegg_2_melt, KO + depth ~ variable, mean)

p <- ggplot(df_dataTable_map_kegg_3, aes(KO, abundance))  +  geom_bar(aes(fill = depth), position = "dodge", stat="identity") + scale_fill_manual(values=c(colsDepth)) + labs(x="Depth",y="Abundance") 
p2 <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p2


p <- ggplot(df_dataTable_map_kegg_3, aes(depth, abundance))  +  geom_bar(aes(fill = depth), position = "dodge", stat="identity") + scale_fill_manual(values=c(colsDepth)) + labs(x="Depth",y="Genome Size") 
p <- p + facet_grid(. ~ KO)
p2 <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + coord_cartesian(ylim = c(0, 1.5)) + geom_hline(yintercept = 1.0, colour="#808080") + geom_hline(yintercept = 0.1, colour="#808080")
p2

ggsave("nitrogenMetabolism_allKOs_depth.svg", bg = "transparent",  width = 80, height = 4, limitsize = FALSE)
