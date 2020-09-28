MeJA_DEGs <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/MeJA_DEGs.txt", sep = "\t", header = T, stringsAsFactors = F)
filename <- list.files("~/Nakano_RNAseq/network_analysis/base/genes_set/CY_time_updown", pattern = ".txt", full.names = T)
object_name <- c()
CY <- c()
CY_MeJA <- c()
temp <- rep(c("CY15", "CY16", "CY20"), each =10)
temp2 <- rep(rep(c("12h", "1h","24h", "3h","48h"), each = 2), times = 3)
temp3 <- rep(c("down", "up"), times = 15)
n <- 1
for (n in 1:length(filename)){
  object_name <- c(object_name, paste0(temp[n], "_", temp2[n], "_FDR005_", temp3[n]))
  test <- assign(object_name[n], read.table(filename[n],  header=T, sep="\t", stringsAsFactors = F))
  CY <- c(CY, nrow(test[n]))
  names(CY)[n] <- object_name[n]
  CY_MeJA <- c(CY_MeJA, length(intersect(rownames(test[n]), MeJA_DEGs$AGI)))
  names(CY_MeJA)[n] <- object_name[n]
  print(n)
}


MeJA_DEGs <- data.frame(expression_change = rep(c("01other", "02up", "03down"), times = 15),
                            Numgenes = c(c(sum(CY[4], CY[3])-sum(CY_MeJA[4], CY_MeJA[3])), CY_MeJA[4], CY_MeJA[3], 
                                         c(sum(CY[8], CY[7])-sum(CY_MeJA[8], CY_MeJA[7])), CY_MeJA[8], CY_MeJA[7], 
                                         c(sum(CY[2], CY[1])-sum(CY_MeJA[2], CY_MeJA[1])), CY_MeJA[2], CY_MeJA[1], 
                                         c(sum(CY[6], CY[5])-sum(CY_MeJA[6], CY_MeJA[5])), CY_MeJA[6], CY_MeJA[5], 
                                         c(sum(CY[10], CY[9])-sum(CY_MeJA[10], CY_MeJA[9])), CY_MeJA[10], CY_MeJA[9], 
                                         c(sum(CY[14], CY[13])-sum(CY_MeJA[14], CY_MeJA[13])), CY_MeJA[14], CY_MeJA[13], 
                                         c(sum(CY[18], CY[17])-sum(CY_MeJA[18], CY_MeJA[17])), CY_MeJA[18], CY_MeJA[17], 
                                         c(sum(CY[12], CY[11])-sum(CY_MeJA[12], CY_MeJA[11])), CY_MeJA[12], CY_MeJA[11], 
                                         c(sum(CY[16], CY[15])-sum(CY_MeJA[16], CY_MeJA[15])), CY_MeJA[16], CY_MeJA[15], 
                                         c(sum(CY[20], CY[19])-sum(CY_MeJA[20], CY_MeJA[19])), CY_MeJA[20], CY_MeJA[19], 
                                         c(sum(CY[24], CY[23])-sum(CY_MeJA[24], CY_MeJA[23])), CY_MeJA[24], CY_MeJA[23], 
                                         c(sum(CY[28], CY[27])-sum(CY_MeJA[28], CY_MeJA[27])), CY_MeJA[28], CY_MeJA[27], 
                                         c(sum(CY[22], CY[21])-sum(CY_MeJA[22], CY_MeJA[21])), CY_MeJA[22], CY_MeJA[21], 
                                         c(sum(CY[26], CY[25])-sum(CY_MeJA[26], CY_MeJA[25])), CY_MeJA[26], CY_MeJA[25], 
                                         c(sum(CY[30], CY[29])-sum(CY_MeJA[30], CY_MeJA[29])), CY_MeJA[30], CY_MeJA[29])
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

  levels(MeJA_DEGs$Category) <- c("1h", "3h", "12h", "24h", "48h")
  levels(MeJA_DEGs$expression_change) <- c("other", "up", "down")
  g <- ggplot(
    MeJA_DEGs,
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