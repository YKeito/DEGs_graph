#test[, 9] #BetweennessCentrality
#test[, 42] #degree


filename <- list.files("~/Nakano_RNAseq/network_analysis/base/mastercluster_info/", pattern = ".csv", full.names = T)
object_name <- c()
CY <- c()
temp <- rep(c("CY15", "CY16", "CY20"), each =5)
temp2 <- rep(c("12h", "1h","24h", "3h","48h"), times = 3)
temp3 <- rep(c("FDR005"), times = 15)
n <- 1
for (n in 1:length(filename)){
  object_name <- c(object_name, paste0(temp[n], "_", temp2[n], "_", temp3[n], "_masterclusterinfo"))
  test <- assign(object_name[n], read.table(filename[n],  header=T, sep=",", stringsAsFactors = F))
  g <- data.frame(AGI = test[, 45],
                  BetweennessCentrality = test[, 9],
                  degree = test[, 42],
                  stringsAsFactors = F
                  )
  gg <- ggplot(g, aes(x = BetweennessCentrality, y = degree))
  gg <- gg + geom_point()
  gg <- gg + theme_bw()
  gg <- gg + scale_color_identity(guide = "legend")
  gg <- gg + theme(legend.position = 'none')
  plot(gg)
  File <- paste0("~/Nakano_RNAseq/network_analysis/results_Fig/mastercluster_info/", object_name[n],"_degree_betweennesscentrality.png")
  ggsave(file = File, plot = gg)
}










sample <- rep("grey20", times = length(subcluster1_node$name))
names(sample) <- subcluster1_node$name
sample[match(subcluster1_node$name[subcluster1_node$TF == "TF"], names(sample))] <- "tomato"

g <- data.frame(AGI = subcluster1_node$name,
                degree = subcluster1_node$Degree,
                BetweennessCentrality = subcluster1_node$BetweennessCentrality,
                color = sample,
                stringsAsFactors = F
)

gg <- ggplot(g, aes(x = BetweennessCentrality, y = degree, color = sample))
gg <- gg + geom_point()
gg <- gg +  theme_bw()
gg <- gg +scale_color_identity(guide = "legend")
gg <- gg + theme(legend.position = 'none')
plot(gg)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/subcluster1/degree_betweennesscentrality.png", plot = gg)
