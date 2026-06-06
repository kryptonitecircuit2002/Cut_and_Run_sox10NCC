library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Drerio.UCSC.danRer10.refGene)
library(org.Dr.eg.db)

h2h3.peaks <- read.csv("....Common_peaks_h2azh33.csv")
##Run GO
entrezids <- h2h3.peaks$geneId %>%
  as.character() %>%
  unique()
ego_h2h3_peaks <- enrichGO(gene = entrezids,
                         keyType = "ENTREZID",
                         OrgDb = org.Dr.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05,
                         readable = TRUE)
cluster_summary <- data.frame(ego_h2h3_peaks)
write.csv(cluster_summary, "GO_down_dar_H2.csv")
processes <- c("cell fate commitment", "cell fate specification", "neural crest cell differentiation", "neural crest cell development", "retinoic acid receptor signaling pathway", "regulation of retinoic acid receptor signaling pathway",
               "muscle structure development", "cranial skeletal system development", "Wnt signaling pathway", "embryonic cranial skeleton morphogenesis", "cell-cell signaling by wnt")
dotplot(ego_h2h3_peaks, showCategory = processes)
