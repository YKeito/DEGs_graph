
#CY15_time
CY15_1h_FDR005 <- allRNASeq[allRNASeq$CY15_1h_q_value < 0.05, ]
write.table(CY15_1h_FDR005, "Nakano_RNAseq/network_analysis/base/genes_set/CY_FDR005_DEGs/CY15_1h_FDR005.txt",append=F, quote = F, sep = "\t", row.names = T)
CY15_3h_FDR005 <- allRNASeq[allRNASeq$CY15_3h_q_value < 0.05, ]
write.table(CY15_3h_FDR005, "Nakano_RNAseq/network_analysis/base/genes_set/CY_FDR005_DEGs/CY15_3h_FDR005.txt",append=F, quote = F, sep = "\t", row.names = T)

CY15_12h_FDR005 <- allRNASeq[allRNASeq$CY15_12h_q_value < 0.05, ]
write.table(CY15_12h_FDR005, "Nakano_RNAseq/network_analysis/base/genes_set/CY_FDR005_DEGs/CY15_12h_FDR005.txt",append=F, quote = F, sep = "\t", row.names = T)

CY15_24h_FDR005 <- allRNASeq[allRNASeq$CY15_24h_q_value < 0.05, ]
write.table(CY15_24h_FDR005, "Nakano_RNAseq/network_analysis/base/genes_set/CY_FDR005_DEGs/CY15_24h_FDR005.txt",append=F, quote = F, sep = "\t", row.names = T)

CY15_48h_FDR005 <- allRNASeq[allRNASeq$CY15_48h_q_value < 0.05, ]
write.table(CY15_48h_FDR005, "Nakano_RNAseq/network_analysis/base/genes_set/CY_FDR005_DEGs/CY15_48h_FDR005.txt",append=F, quote = F, sep = "\t", row.names = T)


#CY16_time
CY16_1h_FDR005 <- allRNASeq[allRNASeq$CY16_1h_q_value < 0.05, ]
write.table(CY16_1h_FDR005, "Nakano_RNAseq/network_analysis/base/genes_set/CY_FDR005_DEGs/CY16_1h_FDR005.txt",append=F, quote = F, sep = "\t", row.names = T)
CY16_3h_FDR005 <- allRNASeq[allRNASeq$CY16_3h_q_value < 0.05, ]
write.table(CY16_3h_FDR005, "Nakano_RNAseq/network_analysis/base/genes_set/CY_FDR005_DEGs/CY16_3h_FDR005.txt",append=F, quote = F, sep = "\t", row.names = T)

CY16_12h_FDR005 <- allRNASeq[allRNASeq$CY16_12h_q_value < 0.05, ]
write.table(CY16_12h_FDR005, "Nakano_RNAseq/network_analysis/base/genes_set/CY_FDR005_DEGs/CY16_12h_FDR005.txt",append=F, quote = F, sep = "\t", row.names = T)

CY16_24h_FDR005 <- allRNASeq[allRNASeq$CY16_24h_q_value < 0.05, ]
write.table(CY16_24h_FDR005, "Nakano_RNAseq/network_analysis/base/genes_set/CY_FDR005_DEGs/CY16_24h_FDR005.txt",append=F, quote = F, sep = "\t", row.names = T)

CY16_48h_FDR005 <- allRNASeq[allRNASeq$CY16_48h_q_value < 0.05, ]
write.table(CY16_48h_FDR005, "Nakano_RNAseq/network_analysis/base/genes_set/CY_FDR005_DEGs/CY16_48h_FDR005.txt",append=F, quote = F, sep = "\t", row.names = T)
#CY20_time
CY20_1h_FDR005 <- allRNASeq[allRNASeq$CY20_1h_q_value < 0.05, ]
write.table(CY20_1h_FDR005, "Nakano_RNAseq/network_analysis/base/genes_set/CY_FDR005_DEGs/CY20_1h_FDR005.txt",append=F, quote = F, sep = "\t", row.names = T)
CY20_3h_FDR005 <- allRNASeq[allRNASeq$CY20_3h_q_value < 0.05, ]
write.table(CY20_3h_FDR005, "Nakano_RNAseq/network_analysis/base/genes_set/CY_FDR005_DEGs/CY20_3h_FDR005.txt",append=F, quote = F, sep = "\t", row.names = T)

CY20_12h_FDR005 <- allRNASeq[allRNASeq$CY20_12h_q_value < 0.05, ]
write.table(CY20_12h_FDR005, "Nakano_RNAseq/network_analysis/base/genes_set/CY_FDR005_DEGs/CY20_12h_FDR005.txt",append=F, quote = F, sep = "\t", row.names = T)

CY20_24h_FDR005 <- allRNASeq[allRNASeq$CY20_24h_q_value < 0.05, ]
write.table(CY20_24h_FDR005, "Nakano_RNAseq/network_analysis/base/genes_set/CY_FDR005_DEGs/CY20_24h_FDR005.txt",append=F, quote = F, sep = "\t", row.names = T)

CY20_48h_FDR005 <- allRNASeq[allRNASeq$CY20_48h_q_value < 0.05, ]
write.table(CY20_48h_FDR005, "Nakano_RNAseq/network_analysis/base/genes_set/CY_FDR005_DEGs/CY20_48h_FDR005.txt",append=F, quote = F, sep = "\t", row.names = T)

####up, down regulation####
####FDR005####
#CY15#
#1h
CY15_1h_FDR005_up <- CY15_1h_FDR005[CY15_1h_FDR005$CY15_1h > 0, ]
CY15_1h_FDR005_down <- CY15_1h_FDR005[CY15_1h_FDR005$CY15_1h < 0, ]
write.table(CY15_1h_FDR005_up, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY15_1h_FDR005_up.txt",append=F, quote = F, sep = "\t", row.names = T)
write.table(CY15_1h_FDR005_down, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY15_1h_FDR005_down.txt",append=F, quote = F, sep = "\t", row.names = T)

#3h
CY15_3h_FDR005_up <- CY15_3h_FDR005[CY15_3h_FDR005$CY15_3h > 0, ]
CY15_3h_FDR005_down <- CY15_3h_FDR005[CY15_3h_FDR005$CY15_3h < 0, ]
write.table(CY15_3h_FDR005_up, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY15_3h_FDR005_up.txt",append=F, quote = F, sep = "\t", row.names = T)
write.table(CY15_3h_FDR005_down, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY15_3h_FDR005_down.txt",append=F, quote = F, sep = "\t", row.names = T)

#12h
CY15_12h_FDR005_up <- CY15_12h_FDR005[CY15_12h_FDR005$CY15_12h > 0, ]
CY15_12h_FDR005_down <- CY15_12h_FDR005[CY15_12h_FDR005$CY15_12h < 0, ]
write.table(CY15_12h_FDR005_up, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY15_12h_FDR005_up.txt",append=F, quote = F, sep = "\t", row.names = T)
write.table(CY15_12h_FDR005_down, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY15_12h_FDR005_down.txt",append=F, quote = F, sep = "\t", row.names = T)

#24h
CY15_24h_FDR005_up <- CY15_24h_FDR005[CY15_24h_FDR005$CY15_24h > 0, ]
CY15_24h_FDR005_down <- CY15_24h_FDR005[CY15_24h_FDR005$CY15_24h < 0, ]
write.table(CY15_24h_FDR005_up, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY15_24h_FDR005_up.txt",append=F, quote = F, sep = "\t", row.names = T)
write.table(CY15_24h_FDR005_down, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY15_24h_FDR005_down.txt",append=F, quote = F, sep = "\t", row.names = T)

#48h
CY15_48h_FDR005_up <- CY15_48h_FDR005[CY15_48h_FDR005$CY15_48h > 0, ]
CY15_48h_FDR005_down <- CY15_48h_FDR005[CY15_48h_FDR005$CY15_48h < 0, ]
write.table(CY15_48h_FDR005_up, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY15_48h_FDR005_up.txt",append=F, quote = F, sep = "\t", row.names = T)
write.table(CY15_48h_FDR005_down, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY15_48h_FDR005_down.txt",append=F, quote = F, sep = "\t", row.names = T)

####FDR005####
#CY16#
#1h
CY16_1h_FDR005_up <- CY16_1h_FDR005[CY16_1h_FDR005$CY16_1h > 0, ]
CY16_1h_FDR005_down <- CY16_1h_FDR005[CY16_1h_FDR005$CY16_1h < 0, ]
write.table(CY16_1h_FDR005_up, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY16_1h_FDR005_up.txt",append=F, quote = F, sep = "\t", row.names = T)
write.table(CY16_1h_FDR005_down, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY16_1h_FDR005_down.txt",append=F, quote = F, sep = "\t", row.names = T)

#3h
CY16_3h_FDR005_up <- CY16_3h_FDR005[CY16_3h_FDR005$CY16_3h > 0, ]
CY16_3h_FDR005_down <- CY16_3h_FDR005[CY16_3h_FDR005$CY16_3h < 0, ]
write.table(CY16_3h_FDR005_up, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY16_3h_FDR005_up.txt",append=F, quote = F, sep = "\t", row.names = T)
write.table(CY16_3h_FDR005_down, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY16_3h_FDR005_down.txt",append=F, quote = F, sep = "\t", row.names = T)

#12h
CY16_12h_FDR005_up <- CY16_12h_FDR005[CY16_12h_FDR005$CY16_12h > 0, ]
CY16_12h_FDR005_down <- CY16_12h_FDR005[CY16_12h_FDR005$CY16_12h < 0, ]
write.table(CY16_12h_FDR005_up, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY16_12h_FDR005_up.txt",append=F, quote = F, sep = "\t", row.names = T)
write.table(CY16_12h_FDR005_down, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY16_12h_FDR005_down.txt",append=F, quote = F, sep = "\t", row.names = T)

#24h
CY16_24h_FDR005_up <- CY16_24h_FDR005[CY16_24h_FDR005$CY16_24h > 0, ]
CY16_24h_FDR005_down <- CY16_24h_FDR005[CY16_24h_FDR005$CY16_24h < 0, ]
write.table(CY16_24h_FDR005_up, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY16_24h_FDR005_up.txt",append=F, quote = F, sep = "\t", row.names = T)
write.table(CY16_24h_FDR005_down, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY16_24h_FDR005_down.txt",append=F, quote = F, sep = "\t", row.names = T)

#48h
CY16_48h_FDR005_up <- CY16_48h_FDR005[CY16_48h_FDR005$CY16_48h > 0, ]
CY16_48h_FDR005_down <- CY16_48h_FDR005[CY16_48h_FDR005$CY16_48h < 0, ]
write.table(CY16_48h_FDR005_up, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY16_48h_FDR005_up.txt",append=F, quote = F, sep = "\t", row.names = T)
write.table(CY16_48h_FDR005_down, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY16_48h_FDR005_down.txt",append=F, quote = F, sep = "\t", row.names = T)

####FDR005####
#CY20#
#1h
CY20_1h_FDR005_up <- CY20_1h_FDR005[CY20_1h_FDR005$CY20_1h > 0, ]
CY20_1h_FDR005_down <- CY20_1h_FDR005[CY20_1h_FDR005$CY20_1h < 0, ]
write.table(CY20_1h_FDR005_up, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY20_1h_FDR005_up.txt",append=F, quote = F, sep = "\t", row.names = T)
write.table(CY20_1h_FDR005_down, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY20_1h_FDR005_down.txt",append=F, quote = F, sep = "\t", row.names = T)

#3h
CY20_3h_FDR005_up <- CY20_3h_FDR005[CY20_3h_FDR005$CY20_3h > 0, ]
CY20_3h_FDR005_down <- CY20_3h_FDR005[CY20_3h_FDR005$CY20_3h < 0, ]
write.table(CY20_3h_FDR005_up, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY20_3h_FDR005_up.txt",append=F, quote = F, sep = "\t", row.names = T)
write.table(CY20_3h_FDR005_down, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY20_3h_FDR005_down.txt",append=F, quote = F, sep = "\t", row.names = T)

#12h
CY20_12h_FDR005_up <- CY20_12h_FDR005[CY20_12h_FDR005$CY20_12h > 0, ]
CY20_12h_FDR005_down <- CY20_12h_FDR005[CY20_12h_FDR005$CY20_12h < 0, ]
write.table(CY20_12h_FDR005_up, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY20_12h_FDR005_up.txt",append=F, quote = F, sep = "\t", row.names = T)
write.table(CY20_12h_FDR005_down, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY20_12h_FDR005_down.txt",append=F, quote = F, sep = "\t", row.names = T)

#24h
CY20_24h_FDR005_up <- CY20_24h_FDR005[CY20_24h_FDR005$CY20_24h > 0, ]
CY20_24h_FDR005_down <- CY20_24h_FDR005[CY20_24h_FDR005$CY20_24h < 0, ]
write.table(CY20_24h_FDR005_up, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY20_24h_FDR005_up.txt",append=F, quote = F, sep = "\t", row.names = T)
write.table(CY20_24h_FDR005_down, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY20_24h_FDR005_down.txt",append=F, quote = F, sep = "\t", row.names = T)

#48h
CY20_48h_FDR005_up <- CY20_48h_FDR005[CY20_48h_FDR005$CY20_48h > 0, ]
CY20_48h_FDR005_down <- CY20_48h_FDR005[CY20_48h_FDR005$CY20_48h < 0, ]
write.table(CY20_48h_FDR005_up, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY20_48h_FDR005_up.txt",append=F, quote = F, sep = "\t", row.names = T)
write.table(CY20_48h_FDR005_down, "Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown/CY20_48h_FDR005_down.txt",append=F, quote = F, sep = "\t", row.names = T)

