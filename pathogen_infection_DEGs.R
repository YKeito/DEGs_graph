Botrytis_cinerea <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/pathogen infection/B.cinerea.txt", sep = "\t", header = T, stringsAsFactors = F)
PstDC3000 <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/pathogen infection/pseudomonas syringae.txt", sep = "\t", header = T, stringsAsFactors = F)
alltarget <- list(Botrytis_cinerea$B.cenerera, PstDC3000$pseudomonas.syringae)

filename <- list.files("~/Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown", pattern = ".txt", full.names = T)
object_name <- c()
CY <- c()
CY_Botrytis <- c()
CY_PstDC3000 <- c()
temp <- rep(c("CY15", "CY16", "CY20"), each =10)
temp2 <- rep(rep(c("12h", "1h","24h", "3h","48h"), each = 2), times = 3)
temp3 <- rep(c("down", "up"), times = 15)
n <- 1
for (n in 1:length(filename)){
    object_name <- c(object_name, paste0(temp[n], "_", temp2[n], "_FDR005_", temp3[n]))
    test <- assign(object_name[n], read.table(filename[n],  header=T, sep="\t", stringsAsFactors = F))
    CY <- c(CY, nrow(test))
    names(CY)[n] <- object_name[n]
    CY_Botrytis <- c(CY_Botrytis, length(intersect(rownames(test), Botrytis_cinerea$B.cenerera)))
    names(CY_Botrytis)[n] <- object_name[n]
    CY_PstDC3000 <- c(CY_PstDC3000, length(intersect(rownames(test), PstDC3000$pseudomonas.syringae)))
    names(CY_PstDC3000)[n] <-  object_name[n]
    print(n)
}

#CYall FDR005
#Botrytis_cinerea
Botrytis_DEGs <- data.frame(expression_change = rep(c("01other", "02up", "03down"), times = 15),
                            Numgenes = c(c(sum(CY[4], CY[3])-sum(CY_Botrytis[4], CY_Botrytis[3])), CY_Botrytis[4], CY_Botrytis[3], 
                                         c(sum(CY[8], CY[7])-sum(CY_Botrytis[8], CY_Botrytis[7])), CY_Botrytis[8], CY_Botrytis[7], 
                                         c(sum(CY[2], CY[1])-sum(CY_Botrytis[2], CY_Botrytis[1])), CY_Botrytis[2], CY_Botrytis[1], 
                                         c(sum(CY[6], CY[5])-sum(CY_Botrytis[6], CY_Botrytis[5])), CY_Botrytis[6], CY_Botrytis[5], 
                                         c(sum(CY[10], CY[9])-sum(CY_Botrytis[10], CY_Botrytis[9])), CY_Botrytis[10], CY_Botrytis[9], 
                                         c(sum(CY[14], CY[13])-sum(CY_Botrytis[14], CY_Botrytis[13])), CY_Botrytis[14], CY_Botrytis[13], 
                                         c(sum(CY[18], CY[17])-sum(CY_Botrytis[18], CY_Botrytis[17])), CY_Botrytis[18], CY_Botrytis[17], 
                                         c(sum(CY[12], CY[11])-sum(CY_Botrytis[12], CY_Botrytis[11])), CY_Botrytis[12], CY_Botrytis[11], 
                                         c(sum(CY[16], CY[15])-sum(CY_Botrytis[16], CY_Botrytis[15])), CY_Botrytis[16], CY_Botrytis[15], 
                                         c(sum(CY[20], CY[19])-sum(CY_Botrytis[20], CY_Botrytis[19])), CY_Botrytis[20], CY_Botrytis[19], 
                                         c(sum(CY[24], CY[23])-sum(CY_Botrytis[24], CY_Botrytis[23])), CY_Botrytis[24], CY_Botrytis[23], 
                                         c(sum(CY[28], CY[27])-sum(CY_Botrytis[28], CY_Botrytis[27])), CY_Botrytis[28], CY_Botrytis[27], 
                                         c(sum(CY[22], CY[21])-sum(CY_Botrytis[22], CY_Botrytis[21])), CY_Botrytis[22], CY_Botrytis[21], 
                                         c(sum(CY[26], CY[25])-sum(CY_Botrytis[26], CY_Botrytis[25])), CY_Botrytis[26], CY_Botrytis[25], 
                                         c(sum(CY[30], CY[29])-sum(CY_Botrytis[30], CY_Botrytis[29])), CY_Botrytis[30], CY_Botrytis[29])
                            ,
                            Category = rep(rep(c("011h", "023h", "0312h", "0424h", "0548h"), each = 3), times = 3),
                            CY = rep(c("CY15", "CY16", "CY20"), each = 15),
                            stringsAsFactors = F
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()

levels(Botrytis_DEGs$Category) <- c("1h", "3h", "12h", "24h", "48h")
levels(Botrytis_DEGs$expression_change) <- c("other", "up", "down")
g <- ggplot(
  Botrytis_DEGs,
  aes (
    x = Category,             # 遺伝子別でグルーピング
    y = Numgenes,
    fill = expression_change       #色塗りつぶしている
    #color = color       #淵枠の色
  )
)
g <- g + geom_bar(stat = "identity")
g <- g +  theme_bw()
g <- g + scale_fill_manual(values = c("grey80", "#F8766D", "#00BFC4"))
#g <- g + scale_colour_identity(guide = "legend")
g <- g + ylab("Number of DEGs")
g <- g + xlab("time-course")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
CYall <- g + facet_grid(CY ~ .)
CYall <- CYall + theme(legend.position = 'none')
plot(CYall)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/pathogen_infection/Botrytis_cinerea.png", plot = CYall)
#PstDC3000
PstDC3000_DEGs <- data.frame(expression_change = rep(c("01other", "02up", "03down"), times = 15),
                             Numgenes = c(c(sum(CY[4], CY[3])-sum(CY_PstDC3000[4], CY_PstDC3000[3])), CY_PstDC3000[4], CY_PstDC3000[3], 
                                          c(sum(CY[8], CY[7])-sum(CY_PstDC3000[8], CY_PstDC3000[7])), CY_PstDC3000[8], CY_PstDC3000[7], 
                                          c(sum(CY[2], CY[1])-sum(CY_PstDC3000[2], CY_PstDC3000[1])), CY_PstDC3000[2], CY_PstDC3000[1], 
                                          c(sum(CY[6], CY[5])-sum(CY_PstDC3000[6], CY_PstDC3000[5])), CY_PstDC3000[6], CY_PstDC3000[5], 
                                          c(sum(CY[10], CY[9])-sum(CY_PstDC3000[10], CY_PstDC3000[9])), CY_PstDC3000[10], CY_PstDC3000[9], 
                                          c(sum(CY[14], CY[13])-sum(CY_PstDC3000[14], CY_PstDC3000[13])), CY_PstDC3000[14], CY_PstDC3000[13], 
                                          c(sum(CY[18], CY[17])-sum(CY_PstDC3000[18], CY_PstDC3000[17])), CY_PstDC3000[18], CY_PstDC3000[17], 
                                          c(sum(CY[12], CY[11])-sum(CY_PstDC3000[12], CY_PstDC3000[11])), CY_PstDC3000[12], CY_PstDC3000[11], 
                                          c(sum(CY[16], CY[15])-sum(CY_PstDC3000[16], CY_PstDC3000[15])), CY_PstDC3000[16], CY_PstDC3000[15], 
                                          c(sum(CY[20], CY[19])-sum(CY_PstDC3000[20], CY_PstDC3000[19])), CY_PstDC3000[20], CY_PstDC3000[19], 
                                          c(sum(CY[24], CY[23])-sum(CY_PstDC3000[24], CY_PstDC3000[23])), CY_PstDC3000[24], CY_PstDC3000[23], 
                                          c(sum(CY[28], CY[27])-sum(CY_PstDC3000[28], CY_PstDC3000[27])), CY_PstDC3000[28], CY_PstDC3000[27], 
                                          c(sum(CY[22], CY[21])-sum(CY_PstDC3000[22], CY_PstDC3000[21])), CY_PstDC3000[22], CY_PstDC3000[21], 
                                          c(sum(CY[26], CY[25])-sum(CY_PstDC3000[26], CY_PstDC3000[25])), CY_PstDC3000[26], CY_PstDC3000[25], 
                                          c(sum(CY[30], CY[29])-sum(CY_PstDC3000[30], CY_PstDC3000[29])), CY_PstDC3000[30], CY_PstDC3000[29]
                                          )
                             ,
                             Category = rep(rep(c("011h", "023h", "0312h", "0424h", "0548h"), each = 3), times = 3),
                             CY = rep(c("CY15", "CY16", "CY20"), each = 15),
                             stringsAsFactors = F
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()

levels(PstDC3000_DEGs$Category) <- c("1h", "3h", "12h", "24h", "48h")
levels(PstDC3000_DEGs$expression_change) <- c("other", "up", "down")
g <- ggplot(
  PstDC3000_DEGs,
  aes (
    x = Category,             # 遺伝子別でグルーピング
    y = Numgenes,
    fill = expression_change       #色塗りつぶしている
    #color = color       #淵枠の色
  )
)
g <- g + geom_bar(stat = "identity")
g <- g +  theme_bw()
g <- g + scale_fill_manual(values = c("grey80", "#F8766D", "#00BFC4"))
#g <- g + scale_colour_identity(guide = "legend")
g <- g + ylab("Number of DEGs")
g <- g + xlab("time-course")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
CYall <- g + facet_grid(CY ~ .)
CYall <- CYall + theme(legend.position = 'none')
plot(CYall)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/pathogen_infection/PstDC3000.png", plot = CYall)
