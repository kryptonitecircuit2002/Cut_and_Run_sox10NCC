library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Drerio.UCSC.danRer10.refGene)
library(org.Dr.eg.db)
txdb <- TxDb.Drerio.UCSC.danRer10.refGene

#Annotate h2az peaks
peakAnno.h2 <- annotatePeak(all_peaks.h2,
                            tssRegion = c(-1000, 1000),
                            TxDb = txdb,
                            annoDb = "org.Dr.eg.db")
plotDistToTSS(peakAnno.h2) + ggtitle("Distribution of H2A.Z binding loci realtive to TSS")
H2.annotation <- as.data.frame(peakAnno.h2@anno)
write.csv(H2.annotation, file = "All_peaks_h2.csv")

#Annotate h3.3 peaks
peakAnno.h3 <- annotatePeak(all_peaks.h3,
                            tssRegion = c(-1000, 1000),
                            TxDb = txdb,
                            annoDb = "org.Dr.eg.db")
H3.annotation <- as.data.frame(peakAnno.h3@anno)
write.csv(H3.annotation, file = "All_peaks_h3.csv")
plotDistToTSS(peakAnno.h3) + ggtitle("Distribution of H3.3 binding loci realtive to TSS")

#Annotate h3k27me3 peaks
peakAnno.h3k27 <- annotatePeak(all_peaks.h3K27,
                            tssRegion = c(-1000, 1000),
                            TxDb = txdb,
                            annoDb = "org.Dr.eg.db")
plotDistToTSS(peakAnno.h3k27) + ggtitle("Distribution of H3K27 loci realtive to TSS")

#Annotate h3k4me3 peaks
peakAnno.h3k4 <- annotatePeak(all_peaks.h3K4,
                               tssRegion = c(-1000, 1000),
                               TxDb = txdb,
                               annoDb = "org.Dr.eg.db")
plotDistToTSS(peakAnno.h3k4) + ggtitle("Distribution of H3K4 loci realtive to TSS")

plotAnnoPie(peakAnno.h2, main = "Peak Distribution of H2A.Z binding regions")
plotAnnoPie(peakAnno.h3, main = "Peak Distribution of H3.3 binding regions")
plotAnnoPie(peakAnno.h3k27, main = "Peak Distribution of regions with H3K27me3") 
plotAnnoPie(peakAnno.h3k4, main = "Peak Distribution of regions with H3K4me3") 

