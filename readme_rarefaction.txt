library(vegan)

setwd("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\mapping_mOTU_5\\rarefaction")

dataTable_map_motu <- read.table("mapping_5.insert.rawCounts.allMGs.uniq.filtered.0.5.COG0012.tsv2", sep="\t", header=TRUE, row.names=1)
samples2cruiseDepth <- read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\samples2cruiseDepth.txt", sep="\t", row.names=1)
names = colnames(dataTable_map_motu)
newNames =  as.vector(samples2cruiseDepth[names,])
colnames(dataTable_map_motu) = newNames


#min(colSums(dataTable_map_motu))
# 1955

dataTable_map_motu_rarefied <- rrarefy(t(dataTable_map_motu), min(colSums(dataTable_map_motu)))
#rarecurve(t(dataTable_map_motu))
rowSums(dataTable_map_motu_rarefied)

layer2Colors <- as.matrix(read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\layer2color_3.txt", sep="\t", row.names=1))
#colsBranches = as.vector(layer2Colors[c("Mixed","aboveDCM","SubDCM","200m","Meso"),])

meta_data <- read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\HOT_2010_11_metadata_4_2.tab", sep="\t", header=TRUE, row.names=2)




speciesAccCurve_all <- specaccum(dataTable_map_motu_rarefied, method="random")
plot(speciesAccCurve_all)


dataTable_map_motu_rarefied

samples0025m <- as.vector(newNames[grep("0025m", newNames)])
samples0045m <- as.vector(newNames[grep("0045m", newNames)])
samples0075m <- as.vector(newNames[grep("0075m", newNames)])

samples0125m <- as.vector(newNames[grep("0125m", newNames)])
samples0200m <- as.vector(newNames[grep("0200m", newNames)])
samples0500m <- as.vector(newNames[grep("0500m", newNames)])
samples0770m <- as.vector(newNames[grep("0770m", newNames)])
samples1000m <- as.vector(newNames[grep("1000m", newNames)])


samples_SRF <- c(samples0025m, samples0045m, samples0075m)
dataTable_map_motu_rarefied_SRF <- dataTable_map_motu_rarefied[samples_SRF,]
dataTable_map_motu_rarefied_SRF = dataTable_map_motu_rarefied_SRF[, colSums(dataTable_map_motu_rarefied_SRF)!=0] 
speciesAccCurve_SRF <- specaccum(dataTable_map_motu_rarefied_SRF, method="random")
#plot(speciesAccCurve_SRF)
plot(speciesAccCurve_SRF, add=TRUE, col=layer2Colors["Mixed",1])


dataTable_map_motu_rarefied_0200m <- dataTable_map_motu_rarefied[samples0200m,]
dataTable_map_motu_rarefied_0200m = dataTable_map_motu_rarefied_0200m[, colSums(dataTable_map_motu_rarefied_0200m)!=0] 
speciesAccCurve_0200m <- specaccum(dataTable_map_motu_rarefied_0200m, method="random")
#plot(speciesAccCurve_0200m)
plot(speciesAccCurve_0200m, add=TRUE, col=layer2Colors["200m",1])

samples0125m_subDCM <- c("HOT237-0125m", "HOT232-0125m", "HOT236-0125m", "HOT238-0125m", "HOT234-0125m","HOT225-0125m","HOT229-0125m","HOT231-0125m")
dataTable_map_motu_rarefied_0125m_subDCM <- dataTable_map_motu_rarefied[samples0125m_subDCM,]
dataTable_map_motu_rarefied_0125m_subDCM = dataTable_map_motu_rarefied_0125m_subDCM[, colSums(dataTable_map_motu_rarefied_0125m_subDCM)!=0] 
speciesAccCurve_0125m_subDCM <- specaccum(dataTable_map_motu_rarefied_0125m_subDCM, method="random")
#plot(speciesAccCurve_0125m_subDCM)
plot(speciesAccCurve_0125m_subDCM, add=TRUE, col=layer2Colors["SubDCM",])


samples_MESO <- c(samples0500m, samples0770m, samples1000m)
dataTable_map_motu_rarefied_MESO <- dataTable_map_motu_rarefied[samples_MESO,]
dataTable_map_motu_rarefied_MESO = dataTable_map_motu_rarefied_MESO[, colSums(dataTable_map_motu_rarefied_MESO)!=0] 
speciesAccCurve_MESO <- specaccum(dataTable_map_motu_rarefied_MESO, method="random")
#plot(speciesAccCurve_MESO)
plot(speciesAccCurve_MESO, add=TRUE, col=layer2Colors["Meso",1])

plot(speciesAccCurve_all)
plot(speciesAccCurve_SRF, add=TRUE, col=layer2Colors["Mixed",1])
plot(speciesAccCurve_0200m, add=TRUE, col=layer2Colors["200m",1])
plot(speciesAccCurve_0125m_subDCM, add=TRUE, col=layer2Colors["SubDCM",])
plot(speciesAccCurve_MESO, add=TRUE, col=layer2Colors["Meso",1])






#C:\work\delong\aloha_gene_catalogue\analysis\mapping_table_mOTUnorm\rarefaction

setwd("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\mapping_mOTU_5\\rarefaction")


#min(colSums(dataTable_map_motu))
# 1955

dataTable_map_motu_rarefied <- rrarefy(t(dataTable_map_motu), min(colSums(dataTable_map_motu)))
#rarecurve(t(dataTable_map_motu))
rowSums(dataTable_map_motu_rarefied)

layer2Colors <- as.matrix(read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\layer2color_3.txt", sep="\t", row.names=1))
#colsBranches = as.vector(layer2Colors[c("Mixed","aboveDCM","SubDCM","200m","Meso"),])

depth2Colors  <- as.matrix(read.table("C:\\work\\delong\\aloha_gene_catalogue\\analysis\\depth2Colors.txt", sep="\t", row.names=1))

speciesAccCurve_all <- specaccum(dataTable_map_motu_rarefied, method="random")
plot(speciesAccCurve_all)

dataTable_map_motu_rarefied

samples0025m <- as.vector(newNames[grep("0025m", newNames)])
samples0045m <- as.vector(newNames[grep("0045m", newNames)])
samples0075m <- as.vector(newNames[grep("0075m", newNames)])
samples0125m <- as.vector(newNames[grep("0125m", newNames)])
samples0200m <- as.vector(newNames[grep("0200m", newNames)])
samples0500m <- as.vector(newNames[grep("0500m", newNames)])
samples0770m <- as.vector(newNames[grep("0770m", newNames)])
samples1000m <- as.vector(newNames[grep("1000m", newNames)])

samples_SRF <- c(samples0025m, samples0045m, samples0075m)
dataTable_map_motu_rarefied_SRF <- dataTable_map_motu_rarefied[samples_SRF,]
dataTable_map_motu_rarefied_SRF = dataTable_map_motu_rarefied_SRF[, colSums(dataTable_map_motu_rarefied_SRF)!=0] 
speciesAccCurve_SRF <- specaccum(dataTable_map_motu_rarefied_SRF, method="random")
#plot(speciesAccCurve_SRF)
plot(speciesAccCurve_SRF, add=TRUE, col=layer2Colors["Mixed",1])

dataTable_map_motu_rarefied_0200m <- dataTable_map_motu_rarefied[samples0200m,]
dataTable_map_motu_rarefied_0200m = dataTable_map_motu_rarefied_0200m[, colSums(dataTable_map_motu_rarefied_0200m)!=0] 
speciesAccCurve_0200m <- specaccum(dataTable_map_motu_rarefied_0200m, method="random")
#plot(speciesAccCurve_0200m)
plot(speciesAccCurve_0200m, add=TRUE, col=layer2Colors["200m",1])

samples0125m_subDCM <- c("HOT237-0125m", "HOT232-0125m", "HOT236-0125m", "HOT238-0125m", "HOT234-0125m","HOT225-0125m","HOT229-0125m","HOT231-0125m")
dataTable_map_motu_rarefied_0125m_subDCM <- dataTable_map_motu_rarefied[samples0125m_subDCM,]
dataTable_map_motu_rarefied_0125m_subDCM = dataTable_map_motu_rarefied_0125m_subDCM[, colSums(dataTable_map_motu_rarefied_0125m_subDCM)!=0] 
speciesAccCurve_0125m_subDCM <- specaccum(dataTable_map_motu_rarefied_0125m_subDCM, method="random")
#plot(speciesAccCurve_0125m_subDCM)
plot(speciesAccCurve_0125m_subDCM, add=TRUE, col=layer2Colors["SubDCM",])

samples_MESO <- c(samples0500m, samples0770m, samples1000m)
dataTable_map_motu_rarefied_MESO <- dataTable_map_motu_rarefied[samples_MESO,]
dataTable_map_motu_rarefied_MESO = dataTable_map_motu_rarefied_MESO[, colSums(dataTable_map_motu_rarefied_MESO)!=0] 
speciesAccCurve_MESO <- specaccum(dataTable_map_motu_rarefied_MESO, method="random")
#plot(speciesAccCurve_MESO)
plot(speciesAccCurve_MESO, add=TRUE, col=layer2Colors["Meso",1])

plot(speciesAccCurve_all)
plot(speciesAccCurve_SRF, add=TRUE, col=layer2Colors["Mixed",1])
plot(speciesAccCurve_0200m, add=TRUE, col=layer2Colors["200m",1])
plot(speciesAccCurve_0125m_subDCM, add=TRUE, col=layer2Colors["SubDCM",])
plot(speciesAccCurve_MESO, add=TRUE, col=layer2Colors["Meso",1])



dataTable_map_motu_rarefied_0125m <- dataTable_map_motu_rarefied[samples0125m,]
dataTable_map_motu_rarefied_0125m = dataTable_map_motu_rarefied_0125m[, colSums(dataTable_map_motu_rarefied_0125m)!=0] 
speciesAccCurve_0125m <- specaccum(dataTable_map_motu_rarefied_0125m, method="random")
#plot(speciesAccCurve_0125m)
plot(speciesAccCurve_0125m, col=depth2Colors["125m",1])

samples0125m_subDCM <- c("HOT237-0125m", "HOT232-0125m", "HOT236-0125m", "HOT238-0125m", "HOT234-0125m","HOT225-0125m","HOT229-0125m","HOT231-0125m")
dataTable_map_motu_rarefied_0125m_subDCM <- dataTable_map_motu_rarefied[samples0125m_subDCM,]
dataTable_map_motu_rarefied_0125m_subDCM = dataTable_map_motu_rarefied_0125m_subDCM[, colSums(dataTable_map_motu_rarefied_0125m_subDCM)!=0] 
speciesAccCurve_0125m_subDCM <- specaccum(dataTable_map_motu_rarefied_0125m_subDCM, method="random")
#plot(speciesAccCurve_0125m_subDCM)
plot(speciesAccCurve_0125m_subDCM, col=layer2Colors["SubDCM",], xlim=c(0,12))

dataTable_map_motu_rarefied_0200m <- dataTable_map_motu_rarefied[samples0200m,]
dataTable_map_motu_rarefied_0200m = dataTable_map_motu_rarefied_0200m[, colSums(dataTable_map_motu_rarefied_0200m)!=0] 
speciesAccCurve_0200m <- specaccum(dataTable_map_motu_rarefied_0200m, method="random")
#plot(speciesAccCurve_0200m)
plot(speciesAccCurve_0200m, add=TRUE, col=layer2Colors["200m",1])

dataTable_map_motu_rarefied_0025m  <- dataTable_map_motu_rarefied[samples0025m,]
dataTable_map_motu_rarefied_0025m = dataTable_map_motu_rarefied_0025m[, colSums(dataTable_map_motu_rarefied_0025m)!=0] 
speciesAccCurve_0025m <- specaccum(dataTable_map_motu_rarefied_0025m, method="random")
plot(speciesAccCurve_0025m, add=TRUE, col=depth2Colors["25m",1])

dataTable_map_motu_rarefied_0075m <- dataTable_map_motu_rarefied[samples0075m,]
dataTable_map_motu_rarefied_0075m = dataTable_map_motu_rarefied_0075m[, colSums(dataTable_map_motu_rarefied_0075m)!=0] 
speciesAccCurve_0075m <- specaccum(dataTable_map_motu_rarefied_0075m, method="random")
#plot(speciesAccCurve_0075m)
plot(speciesAccCurve_0075m, add=TRUE, col=depth2Colors["75m",1])


dataTable_map_motu_rarefied_0500m <- dataTable_map_motu_rarefied[samples0500m,]
dataTable_map_motu_rarefied_0500m = dataTable_map_motu_rarefied_0500m[, colSums(dataTable_map_motu_rarefied_0500m)!=0] 
speciesAccCurve_0500m <- specaccum(dataTable_map_motu_rarefied_0500m, method="random")
#plot(speciesAccCurve_0500m)
plot(speciesAccCurve_0500m, add=TRUE, col=depth2Colors["500m",1])

dataTable_map_motu_rarefied_0770m <- dataTable_map_motu_rarefied[samples0770m,]
dataTable_map_motu_rarefied_0770m = dataTable_map_motu_rarefied_0770m[, colSums(dataTable_map_motu_rarefied_0770m)!=0] 
speciesAccCurve_0770m <- specaccum(dataTable_map_motu_rarefied_0770m, method="random")
#plot(speciesAccCurve_0770m)
plot(speciesAccCurve_0770m, add=TRUE, col=depth2Colors["770m",1])

dataTable_map_motu_rarefied_1000m <- dataTable_map_motu_rarefied[samples1000m,]
dataTable_map_motu_rarefied_1000m = dataTable_map_motu_rarefied_1000m[, colSums(dataTable_map_motu_rarefied_1000m)!=0] 
speciesAccCurve_1000m <- specaccum(dataTable_map_motu_rarefied_1000m, method="random")
#plot(speciesAccCurve_1000m)
plot(speciesAccCurve_1000m, add=TRUE, col=depth2Colors["1000m",1])

plot(speciesAccCurve_0200m, col=layer2Colors["200m",1], xlim=c(1,12), ylim=c(0,1300),ylab="# of species", xlab="# of samples",lwd = 4)
plot(speciesAccCurve_0125m_subDCM, add=TRUE, col=layer2Colors["SubDCM",], xlim=c(1,12),lwd = 4)
plot(speciesAccCurve_0025m, add=TRUE, col=depth2Colors["25m",1],lwd = 4)
plot(speciesAccCurve_0075m, add=TRUE, col=depth2Colors["75m",1],lwd = 4)
plot(speciesAccCurve_0500m, add=TRUE, col=depth2Colors["500m",1],lwd = 4)
plot(speciesAccCurve_0770m, add=TRUE, col=depth2Colors["770m",1],lwd = 4)
plot(speciesAccCurve_1000m, add=TRUE, col=depth2Colors["1000m",1],lwd = 4)


pdf(file = "speciesAccCurve_depth.pdf", 14, 10, bg = "transparent")
plot(speciesAccCurve_0200m, col=layer2Colors["200m",1], xlim=c(1,12), ylim=c(0,1300),ylab="# of species", xlab="# of samples",lwd = 4)
plot(speciesAccCurve_0125m_subDCM, add=TRUE, col=layer2Colors["SubDCM",], xlim=c(1,12),lwd = 4)
plot(speciesAccCurve_0025m, add=TRUE, col=depth2Colors["25m",1],lwd = 4)
plot(speciesAccCurve_0075m, add=TRUE, col=depth2Colors["75m",1],lwd = 4)
plot(speciesAccCurve_0500m, add=TRUE, col=depth2Colors["500m",1],lwd = 4)
plot(speciesAccCurve_0770m, add=TRUE, col=depth2Colors["770m",1],lwd = 4)
plot(speciesAccCurve_1000m, add=TRUE, col=depth2Colors["1000m",1],lwd = 4)
dev.off()







speciesNum_0025m <- specnumber(dataTable_map_motu_rarefied_0025m)
speciesNum_0075m <- specnumber(dataTable_map_motu_rarefied_0075m)
speciesNum_0125m <- specnumber(dataTable_map_motu_rarefied_0125m)
speciesNum_0200m <- specnumber(dataTable_map_motu_rarefied_0200m)
speciesNum_0500m <- specnumber(dataTable_map_motu_rarefied_0500m)
speciesNum_0770m <- specnumber(dataTable_map_motu_rarefied_0770m)
speciesNum_1000m <- specnumber(dataTable_map_motu_rarefied_1000m)

boxplot(speciesNum_0025m, speciesNum_0075m, speciesNum_0125m, speciesNum_0200m, speciesNum_0500m, speciesNum_0770m, speciesNum_1000m)


library("ggplot2")
speciesNum <- specnumber(dataTable_map_motu_rarefied)

depth <- as.matrix(meta_data$depth_m)
rownames(depth) <- rownames(meta_data)
colnames(depth) <- c("depth")
#colnames(depth) <- c("depth")
speciesNumDF <- as.data.frame(speciesNum)

sampleNames <- rownames(speciesNumDF)
speciesNumDF$depth <- depth[sampleNames,]

speciesNumDF2 <- speciesNumDF
speciesNumDF <- speciesNumDF2[speciesNumDF2$depth != "45m",]

samples0125m_aboveDMC <- c("HOT224-0125m", "HOT227-0125m", "HOT233-0125m")
speciesNumDF_aboveDCM <- speciesNumDF[samples0125m_aboveDMC,]
speciesNumDF_aboveDCM[samples0125m_aboveDMC,"depth"] <- "125m"
colnames(speciesNumDF_aboveDCM) <- c("speciesNum","depth")

samplesNames <- row.names(speciesNumDF)
samplesNames_filtered <-samplesNames [! samplesNames %in% c("HOT224-0125m", "HOT227-0125m", "HOT233-0125m")]


speciesNumDF3 <- speciesNumDF[samplesNames_filtered, ]

speciesNumDF <- speciesNumDF3

#speciesNumDF$depth  <- (factor((speciesNumDF$depth), levels=c("25m", "45m", "75m", "125m", "200m", "500m", "770m", "1000m")))
#speciesNumDF$depth  <- (factor((speciesNumDF$depth), levels=c("1000m", "770m", "500m", "200m", "125m", "75m", "45m", "25m")))
speciesNumDF$depth  <- (factor((speciesNumDF$depth), levels=c("1000m", "770m", "500m", "200m", "125m", "75m", "25m")))


#depth = rev(c("25m", "45m", "75m", "125m", "200m", "500m", "770m", "1000m"))
depth = rev(c("25m", "75m", "125m", "200m", "500m", "770m", "1000m"))
colsDepth = as.vector(depth2Colors[depth, 1])
greyCols = c("#202020","#202020","#202020","#202020","#202020","#202020","#202020","#202020")


p <- ggplot(speciesNumDF, aes(x=depth, y=speciesNum, color=depth)) + scale_color_manual(values=c(colsDepth)) + geom_jitter(position=position_jitter(0.3), cex=5) + labs(x="Depth",y="Number of Species") 
p2 <- p + geom_jitter(data=speciesNumDF_aboveDCM,  position=position_jitter(0.3), cex=4, shape=24, fill = "#d0d0d0", colour = "#000000")
#p2 <- p
p3 <- p2 + theme_bw()
p3

ggsave("speciesRichness.svg", bg = "transparent",  width = 10, height = 8)
ggsave("speciesRichness_bg.svg", bg = "white",  width = 10, height = 8)

p <- ggplot(speciesNumDF, aes(x=depth, y=speciesNum, fill=depth, color=depth )) + scale_color_manual(values=c(greyCols)) + scale_fill_manual(values=c(colsDepth)) + geom_boxplot(outlier.shape = NA, lwd=2, fatten=0.75) + geom_jitter(position=position_jitter(0.3), cex=5, shape=21) + labs(x="Depth",y="Number of Species")
p2 <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_flip()
p2
p3 <- p2 + geom_jitter(data=speciesNumDF_aboveDCM,  position=position_jitter(0.3), cex=5, shape=24, fill = "#d0d0d0", colour = "#000000")
p3

ggsave("speciesRichness_boxplot.svg", bg = "transparent",  width = 9, height = 8)
ggsave("speciesRichness_boxplot_bg.svg", bg = "white",  width = 9, height = 8)


#colored rarecurve
sampleNames <- colnames(dataTable_map_motu)
samples_depth <- depth[sampleNames,]
samples_colors <- depth2Colors[samples_depth,]

rarecurve(t(dataTable_map_motu),col=c(samples_colors), cex = 0.0)








p <- ggplot(speciesNumDF, aes(x=depth, y=speciesNum, fill=depth, color=depth )) + scale_color_manual(values=c(greyCols)) + scale_fill_manual(values=c(colsDepth)) + geom_boxplot(outlier.shape = NA, lwd=2, fatten=0.75) + geom_jitter(position=position_jitter(0.3), cex=5, shape=21) + labs(x="Depth",y="Number of Species")
p2 <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_flip()
p2




