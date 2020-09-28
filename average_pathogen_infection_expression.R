Botrytis_cinerea <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/pathogen infection/B.cinerea.txt", sep = "\t", header = T, stringsAsFactors = F)
PstDC3000 <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/pathogen infection/pseudomonas syringae.txt", sep = "\t", header = T, stringsAsFactors = F)
alltarget <- list(Botrytis_cinerea$B.cenerera, PstDC3000$pseudomonas.syringae)

expression_velue <- list(CY15_1h_FDR005$CY15_1h, CY15_3h_FDR005$CY15_3h, CY15_12h_FDR005$CY15_12h, CY15_24h_FDR005$CY15_24h, CY15_48h_FDR005$CY15_48h,
                      CY16_1h_FDR005$CY16_1h, CY16_3h_FDR005$CY16_3h, CY16_12h_FDR005$CY16_12h, CY16_24h_FDR005$CY16_24h, CY16_48h_FDR005$CY16_48h,
                      CY20_1h_FDR005$CY20_1h, CY20_3h_FDR005$CY20_3h, CY20_12h_FDR005$CY20_12h, CY20_24h_FDR005$CY20_24h, CY20_48h_FDR005$CY20_48h
                      )
AGI <- list(rownames(CY15_1h_FDR005), rownames(CY15_3h_FDR005), rownames(CY15_12h_FDR005), rownames(CY15_24h_FDR005), rownames(CY15_48h_FDR005),
            rownames(CY16_1h_FDR005), rownames(CY16_3h_FDR005), rownames(CY16_12h_FDR005), rownames(CY16_24h_FDR005), rownames(CY16_48h_FDR005),
            rownames(CY20_1h_FDR005), rownames(CY20_3h_FDR005), rownames(CY20_12h_FDR005), rownames(CY20_24h_FDR005), rownames(CY20_48h_FDR005)
            )
CY_Botrytis <- c()
CY_PstDC3000 <- c()
n <- 1
for (n in 1:length(AGI)){
  t <- intersect(AGI[[n]], Botrytis_cinerea$B.cenerera)
  ex <- expression_velue[[n]]
  data <- ex[match(t, AGI[[n]])]
  CY_Botrytis <- c(CY_Botrytis, mean(data))
  
  t <- intersect(AGI[[n]], PstDC3000$pseudomonas.syringae)
  ex <- expression_velue[[n]]
  data <- ex[match(t, AGI[[n]])]
  CY_PstDC3000 <- c(CY_PstDC3000, mean(data))
  print(n)
}


#CYall FDR005
#Botrytis_cinerea
pathogen_mean <- data.frame(value = c(CY_Botrytis, CY_PstDC3000),
                            Category = rep(c("011h", "023h", "0312h", "0424h", "0548h"), times = 6),
                            CY = rep(rep(c("CY15", "CY16", "CY20"), each = 5), times = 2),
                            pathogen = rep(c("Botrytis_cinerea", "PstDC3000"), each = 15),
                            stringsAsFactors = F
)
library(ggplot2)
library(reshape2)
library(ggsci)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()

levels(pathogen_mean$Category) <- c("1h", "3h", "12h", "24h", "48h")

g <- ggplot(
  pathogen_mean,
  aes (
    x = Category,             # 遺伝子別でグルーピング
    y = value,
    fill =  pathogen      #色塗りつぶしている
    #color = color       #淵枠の色
  )
)
g <- g + geom_bar(stat = "identity", position = "dodge")
g <- g +  theme_bw()
g <- g + scale_fill_manual(values = c("black", "grey80"))
g <- g + ylab("Exprssion average")
g <- g + xlab("time-course")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
CYall <- g + facet_grid(CY ~ .)
plot(CYall)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/pathogen_infection/average_expression.png", plot = CYall)