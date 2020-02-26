##pca-som


library(FactoMineR)
library(missMDA)
library(factoextra)
library(corrplot)
library(ggplot2)
library(kohonen)
library(tempR)

#setwd("")
#data = read.table("", header=T, sep="\t", stringsAsFactors  = FALSE, na.strings=c("na", "NA", "NaN", ""))
## data: row-sample, col-gene


setwd("C:/Users/user/Desktop/ÇýÁö")
data = read.csv("CCLE_depMap_19Q1_TPM.csv", header=T, sep=",", 
                stringsAsFactors  = FALSE, na.strings=c("na", "NA", "NaN", ""))
order = names(data) # header¸¸ ÃßÃâ

t_data = read.table("TNBC_list.txt", header=T, sep="\t", stringsAsFactors  = FALSE, na.strings=c("na", "NA", "NaN", ""))
nt_data = read.table("nonTNBC_list.txt", header=T, sep="\t", stringsAsFactors  = FALSE, na.strings=c("na", "NA", "NaN", ""))



bc_data = matrix(nrow=0, ncol=ncol(data))
for(i in t_data$ID) {
  subtype = c("tnbc")
  bc_data= rbind(bc_data, cbind(subset(data, ID == i), subtype)) 
}


for(i in nt_data$depMapID) {
  subtype = c("non-tnbc")
  if(nrow(subset(data, ID == i)) == 0 ) {
    next 
  }
  bc_data= rbind(bc_data, cbind(subset(data, ID == i), subtype)) 
}


#bc_data = read.csv("breast_cell_line.csv", header=T, sep=",", 
                   stringsAsFactors  = FALSE, na.strings=c("na", "NA", "NaN", ""))
#bc_order = names(bc_data)


# data process 

raw = suppressWarnings(as.data.frame(sapply(bc_data[,], as.numeric)))
row.names(raw) = raw$ID

hh = 350
ww = 400
raw = as.numeric(bc_data)
#pre PCA for missing values
nb = estim_ncpPCA(bc_data[,3:ncol(bc_data)], ncp.max = 5)
res.comp = imputePCA(bc_data[,3:ncol(bc_data)], ncp=nb$ncp)

#PCA
pca_res = PCA(res.comp$completeObs, scale.unit=FALSE, graph = FALSE, ncp = 5)
#eigen value
pca.eig = pca_res$eig

# Results for Variables
res.var <- get_pca_var(pca_res)
res.var$coord          # Coordinates
plot(res.var$contrib)
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

# Results for individuals
res.ind <- get_pca_ind(pca_res)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation

# export PCA results
write.table(pca.eig, file = "pca.eig.txt", sep="\t")
write.table(res.var$coord, file = "pca_var_coord.txt", sep="\t")
write.table(res.var$contrib, file = "pca_var_contrib.txt", sep="\t")
write.table(res.var$cos2, file = "pca_var_cos2.txt", sep="\t")
write.table(res.ind$coord, file = "pca_ind_coord.txt", sep="\t")
write.table(res.ind$contrib, file = "pca_ind_contrib.txt", sep="\t")
write.table(res.ind$cos2, file = "pca_ind_cos2.txt", sep="\t")

png = paste0("screeplot.png")
png(filename=png, height=450, width=520, bg="white") # 
bar=barplot(pca_res$eig[1:10,2],main="Variances",names.arg=1:10, cex.names = 1.8, cex.lab = 1.8, cex.main = 2, cex.axis = 1.8,
            ylab = "Percentage of variances", xlab = "Principal Compnents", ylim = c(0,100))
label = round(pca_res$eig[1:10,2])
lines(x = bar, y = label, col = "blue", type = "o")
text(bar, pca_res$eig[1:10,2]+3, labels = paste0(label, "%"), cex = 1.5, font = 4, col = "black")
garbage = dev.off() 

##cumulative percentage of variance
png = paste0("cumulative percentage of variance.png")
png(filename=png, height=450, width=520, bg="white") # 
bar=barplot(pca_res$eig[1:10,3],main="Variances", names.arg=1:10,cex.names = 1.8, cex.lab = 1.8, cex.main = 2, cex.axis = 1.8,
            ylab = "Percentage of variances", xlab = "Principal Compnents", ylim = c(0,100))
label = round(pca_res$eig[1:10,3])
lines(x = bar, y = label, col = "blue", type = "o")
text(bar, pca_res$eig[1:10,3]-3, labels = paste0(label, "%"), cex = 1.5, font = 4, col = "black")
garbage = dev.off() 

ind = as.data.frame(pca_res$ind)
color = as.character(raw$color)
pca_result = cbind(ind, color)
#theme
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             strip.background=element_blank(),axis.text.x=element_text(colour="black", size = 18),
             axis.text.y=element_text(colour="black", size = 18),axis.ticks=element_line(colour="black"), axis.title =element_text(size = 18),
             plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
             plot.margin=unit(c(1,1,1,1),"line"))
percentage <- round(pca_res$eig[1:10,2], 2)
percentage <- paste( colnames(ind), "(", paste( as.character(percentage), "%", ")", sep="") )

# drawing individual plot 
png(filename="individual plot.png", height=400, width=450, bg="white")
p<-ggplot(pca_som,aes(x = coord.Dim.1, y = coord.Dim.2, bg = color))
p<-p+ theme + xlab(percentage[1]) + ylab(percentage[2]) + ggtitle("PCA result")
p
garbage = dev.off()




#########################################################################################
#SOM 

#Calculate SOM closest distance

i =  1 
for(i in 1:){
  j = 1
  while((i*j) <= ){
    som_grid = somgrid(i, j, topo = "hexagonal")
    som_model <- som(data, 
                     grid=som_grid, 
                     rlen=500, 
                     alpha=c(0.01,0.01), 
                     keep.data = TRUE,
                     maxNA.fraction = 5)
    print(paste0("done: ", i,"  ", j))
    if(i == 1 & j == 1){
      out = c(i,j,mean(som_model$distances))
      names(out) = c("row", "col", "error")
    }
    else{
      str = c(i,j,mean(som_model$distances))
      out = rbind(out,str)
    }
    j = j+1
  }
}

write.table(out, "closest distance.txt", sep = "\t", row.names = F)

## check minimum closest distance grids counts


#SOM
row = rownames(data)
col = colnames(data)
som.data = as.matrix(data)
rownamse(data) = row
colnames(data) = col


som_grid <- somgrid(, , topo="hexagonal")
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
plot(som_model, type = "code", main = "Node Counts", shape = "straight", border = "WHITE") # Node count plot
dev.off()

png("distance.png", width = 400, height = 300)
plot(som_model, type = "dist.neighbours", main = "SOM neighbour distances", shape = "straight", border = "WHITE") # U-matrix visualisation
dev.off()
#plot(som_model, type = "codes") # Weight Vector View

#clustering of grids
codes = as.data.frame(som_model$codes)
som_cluster <- cutree(hclust(dist(codes)), 3)
cluster_assignment = som_cluster[som_model$unit.classif]

label = rownames(data)

pretty_palette = c("#FFF0FD", "#FFFFBA","#BAE1FF", "#c6e5d9", "#fcdfc9")

plot(som_model, type="mapping", main = "Clusters", keepMargins = TRUE, shape = "straight", pchs = 21,
     cex = 0.5, col = "black", bg = color, bgcol = pretty_palette[som_cluster], codeRendering = "segments", 
     border = "grey", classif = som_model$unit.classif, labels = name)
png("cluster3(filtering)20190309.png", width = 400, height = 500)
plot(som_model, type="mapping", main = "Clusters", keepMargins = TRUE, shape = "straight", pchs = 21, cex.main = 3,
     cex = 1, col = "black", bg = color, bgcol = pretty_palette[som_cluster], codeRendering = "segments", 
     border = "grey", classif = som_model$unit.classif)

garbage = dev.off()
mean(som_model$distances) 
#add.cluster.boundaries(som_model, som_cluster, lwd =5, col = "white")
#Kohonen heatmap creation 

unit = cbind(label, som_model$unit.classif, cluster_assignment)
colnames(unit) = c("Samples", "Grid No", "Cluster Group")


write.table(unit, "unit2.txt", sep = "\t", row.names = FALSE)
write.table(codes, "code.txt", sep = "\t", row.names = TRUE)
saveRDS(som_model, "som_model.rds")

##PCA-SOM result
cluster = factor(cluster_assignment)
png(filename="PCA-SOM result.png", height=400, width=450, bg="white")
p<-ggplot(pca_result,aes(x = coord.Dim.1, y = coord.Dim.2, shape = cluster))+
  scale_shape_manual(values = c(21,22,25,23,24))
p<-p+geom_point(size = 4, bg = color, col = "black")+ theme + xlab(percentage[1]) + ylab(percentage[2]) + 
  ggtitle("PCA-SOM result")
p
garbage = dev.off()

