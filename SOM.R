# self oraganising map (single cell only)

library("kohonen")
library("tempR")
setwd("C:/Users/user/Desktop/single cell/pca/output/2199 genes")


raw = as.data.frame(read.table("C:/Users/user/Desktop/single cell/DATA/raw(2199).txt", sep = "\t",
                               header = T, stringsAsFactors = FALSE))

# common_data = as.data.frame(read.table("C:/Users/user/Desktop/single cell/pca/data/common_filtering(20190207).txt", sep = "\t",
#                                        header = T, stringsAsFactors = FALSE))
# common_list = common_data$Gene
# common_list = c("sample", "color", common_list)
# raw = raw[,common_list]
# row.names(raw) = raw[,1]


# zscore = as.data.frame(raw[2:nrow(raw), ], stringsAsFactors = FALSE)
# rownames(zscore) = zscore[,1]
# zscore = zscore[,-1]
# 
# zscore = as.data.frame(sapply(zscore, as.numeric))
# 
# write.table(zscore, "zscore.txt", sep = "\t")
name = raw[,1]
zscore = as.data.frame(raw[, 3:ncol(raw)], stringsAsFactors = FALSE)
rownames(zscore) = name
zscore = as.data.frame(sapply(zscore, as.numeric))
# zscore = sapply(zscore, function(x){log2(x)})
# is.na(zscore) = sapply(zscore, is.infinite)
zscore = as.data.frame(zscore)
zscore[195,] = apply(zscore[1:74,], 2, mean, na.rm = T)
zscore[196,] = apply(zscore[75:121,], 2, mean, na.rm = T)
raw[195:196,2] = "forestgreen"

zscore = as.matrix(zscore)
color = unlist(raw[,2])
#SOM
som_grid <- somgrid(14, 13, topo="hexagonal")
som_model <- som(zscore, 
                 grid=som_grid, 
                 rlen=500, 
                 alpha=c(0.01,0.01), 
                 keep.data = TRUE,
                 maxNA.fraction = 5)


coolBlueHotRed <- function(n, alpha = 1) {rainbow(n, end=4/6, alpha=alpha)[n:1]}
#par(mar=c(7,8,7,8))

png("quaility.png", width = 400, height = 300)
  plot(som_model, type = "quality", main = "Quality", keepMargins = TRUE, shape = "straight",codeRendering = "segments", 
       border = "black")
  #Training progress for SOM
dev.off()

png("count.png", width = 400, height = 300)
plot(som_model, type = "count", main = "Node Counts", shape = "straight", border = "WHITE") # Node count plot
dev.off()

png("distance.png", width = 400, height = 300)
plot(som_model, type = "dist.neighbours", main = "SOM neighbour distances", shape = "straight", border = "WHITE") # U-matrix visualisation
dev.off()
#plot(som_model, type = "codes") # Weight Vector View

#clustering of grids
codes = as.data.frame(som_model$codes)
som_cluster <- cutree(hclust(dist(codes)), 4)
cluster_assignment = som_cluster[som_model$unit.classif]

#drawing property heatmap
# setwd("C:/Users/user/Desktop/single cell/som/png")
# for(i in 1:3197){
#   
#   png(paste0("property", 5 , ".png"), width = 500, height = 300)
#   plot(som_model, type = "property", property = getCodes(som_model)[,5], main = colnames(getCodes(som_model))[5], palette.name = coolBlueHotRed)
#   dev.off()
#   
# }
#write.table(som_model$codes, "codes.txt", sep = "\t")

label = colnames(zscore)
# color[color == "darkblue"] = "blue3"
# color[color == "orange"] = "yellow1"
# edge = color
# edge[color == "yellow1"] = "orange"
# edge[color == "blue3"] = "gray48"
# edge[color == "grey"] = "black"
# edge[color == "red"] = "red4"
pretty_palette = c("#FFF0FD", "#FFFFBA","#BAE1FF", "#c6e5d9", "#fcdfc9")

plot(som_model, type="mapping", main = "Clusters", keepMargins = TRUE, shape = "straight", pchs = 21,
     cex = 0.5, col = "black", bg = color, bgcol = pretty_palette[som_cluster], codeRendering = "segments", 
     border = "grey", classif = som_model$unit.classif, labels = name)
png("cluster5(filtering)20190309.png", width = 500, height = 400)
plot(som_model, type="mapping", main = "Clusters", keepMargins = TRUE, shape = "straight", pchs = 21, cex.main = 2,
     cex = 1.8, col = "black", bg = color, bgcol = pretty_palette[som_cluster], codeRendering = "segments", 
     border = "grey", classif = som_model$unit.classif)

garbage = dev.off()
mean(som_model$distances)
#add.cluster.boundaries(som_model, som_cluster, lwd =5, col = "white")
#Kohonen heatmap creation 

unit = cbind(label, som_model$unit.classif, cluster_assignment)
colnames(unit) = c("Samples", "Grid No", "Cluster Group")

#write.table(data,"cluster_result.txt", sep = "\t")
write.table(unit, "unit.txt", sep = "\t", row.names = FALSE)
write.table(codes, "code.txt", sep = "\t", row.names = TRUE)
saveRDS(som_model, "som_model.rds")
som_model = readRDS("som_model.rds")
