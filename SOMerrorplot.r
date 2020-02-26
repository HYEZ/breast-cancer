# SOM error plot 

setwd("C:/Users/user/Desktop/single cell/pca/output/2199 genes")
raw = read.table("PCA-SOM closet distance.txt", sep = "\t", stringsAsFactors = T, header = T, na.strings = c("NA", "#N/A", "NA!", ""))

png("pca-som closest distance.png", width = 600, height = 400)

plot(raw$error, type = "l", xaxt = "n", xlab = "", main = "Avg of closest distance", ylab = "", cex.main = 2, cex.axis = 1.5)
min = which.min(raw$error) #min components
points(min, min(raw$error), pch = 1, col = "red", cex = 3)
#text((min+70), min(raw$error),labels = paste0(raw$row[min], "x",raw$col[min]), col = "red", cex = 1.5)
axis(1, at = seq(1, 1053, by = 100),  labels = raw$row[(seq(1,1053, by = 100))], las = 1, cex.axis = 1.3, tick = T, line = 0)
mtext("row", las = 1, side = 2, at = 2400, cex = 1.5)
axis(1, at = seq(1, 1053, by = 100),  labels = raw$col[(seq(1,1053, by = 100))], las = 1, cex.axis = 1.3, tick = F, line = 1)
mtext("col", las = 1, side = 2, at = 1600, cex = 1.5)
mtext("# of SOM grids", las = 1, side = 1, line = 4, cex = 1.8)
trash = dev.off()