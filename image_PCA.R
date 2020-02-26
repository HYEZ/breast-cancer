#ref.
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/




#install.packages("package name") 

###################################################
library(FactoMineR)
library(missMDA)
library(factoextra)
library(corrplot)
library(ggplot2)

setwd("C:/Users/user/Desktop/혜지")
data = read.csv("CCLE_depMap_19Q1_TPM.csv", header=T, sep=",", 
                  stringsAsFactors  = FALSE, na.strings=c("na", "NA", "NaN", ""))
order = names(data) # header만 추출

bc_data = read.csv("breast_cell_line.csv", header=T, sep=",", 
                   stringsAsFactors  = FALSE, na.strings=c("na", "NA", "NaN", ""))
bc_order = names(bc_data)


for(i in data) {
  
}




raw = data[,3:ncol(data)]

raw = as.data.frame(sapply(raw, as.numeric))
raw[195,] = apply(raw[1:74,], 2, mean, na.rm = T) # 2: 세로(열), 1:가로
raw[196,] = apply(raw[75:121,], 2, mean, na.rm = T) 
raw$sample = c(data$sample, "A549_avg", "H1437_avg") # $뒤에는 피쳐 이름
raw$color = c(data$color, "forestgreen", "forestgreen")
raw = raw[,order]
row.names(raw) = raw$sample              
hh = 350
ww = 400

#PCA
#pre PCA for missing values
nb = estim_ncpPCA(raw[,3:ncol(data)], ncp.max = 5)
res.comp = imputePCA(raw[,3:ncol(data)], ncp=nb$ncp)

#PCA
pca_res = PCA(res.comp$completeObs, scale.unit=FALSE, graph = FALSE, ncp = 5)
#eigen value
pca.eig = pca_res$eig




#draw scree plot
par(mar=c(7,8,7,8))
png = paste0(name, "_eigenvalue(single).png")
png(filename=png, height=hh, width=ww, bg="white") # 
fviz_eig(pca_res, addlabels=TRUE, ylim=c(0,30))
garbage = dev.off() 


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

write.table(pca.eig, file = "pca.eig(single).txt", sep="\t")

write.table(res.var$coord, file = "pca_var_coord(0124).txt", sep="\t")
write.table(res.var$contrib, file = "pca_var_contrib(single).txt", sep="\t")
write.table(res.var$cos2, file = "pca_var_cos2(single).txt", sep="\t")

write.table(res.ind$coord, file = "pca_ind_coord(a549).txt", sep="\t")
write.table(res.ind$contrib, file = "pca_ind_contrib(single).txt", sep="\t")
write.table(res.ind$cos2, file = "pca_ind_cos2(0124).txt", sep="\t")
par(mar=c(0,0,0,0))
par(mar=c(7,8,7,8))
#same with scree plot
png = paste0("filtering", "_barplot(2019309).png")
png(filename=png, height=450, width=520, bg="white") # 
bar=barplot(pca_res$eig[1:10,2],main="Variances",names.arg=1:10, cex.names = 1.8, cex.lab = 1.8, cex.main = 2, cex.axis = 1.8,
        ylab = "Percentage of variances", xlab = "Principal Compnents", ylim = c(0,100))
label = round(pca_res$eig[1:10,2])
lines(x = bar, y = label, col = "blue", type = "o")
text(bar, pca_res$eig[1:10,2]+3, labels = paste0(label, "%"), cex = 1.5, font = 4, col = "black")
garbage = dev.off() 

##cumulative percentage of variance
png = paste0("raw", "_cumulative percentage of variance(2019309).png")
png(filename=png, height=hh, width=ww, bg="white") # 
bar=barplot(pca_res$eig[1:10,3],main="Variances", names.arg=1:10,
            ylab = "Percentage of variances", xlab = "Principal Compnents", ylim = c(0,100))
label = round(pca_res$eig[1:10,3])
lines(x = bar, y = label, col = "blue", type = "o")
text(bar, pca_res$eig[1:10,3]-3, labels = paste0(label, "%"), cex = 0.8, font = 4, col = "black")
garbage = dev.off() 
# draw coord of variables
fviz_pca_var(pca_res, col.var = "black")
#filtering results: cos2>=0.9
fviz_pca_var(pca_res, select.var=list(cos2=0.9))
#filtering results: the highest 50 cos2
fviz_pca_var(pca_res, select.var= list(cos2 = 50))

#draw contributions of variables: top
fviz_contrib(pca_res, choice = "var", axes = 1, top = 50)
fviz_contrib(pca_res, choice = "var", axes = 2, top = 50)


# individual and PCs
    pca.ind = pca_res$ind$coord
    write.table(pca.ind, file = "pca.ind(single).txt", sep="\t")

    pca.dimdesc = dimdesc(pca_res, axes=1:2)
    write.table(pca.dimdesc$Dim.1, file = "pca.dimdesc1.txt", sep="\t")
    write.table(pca.dimdesc$Dim.2, file = "pca.dimdesc2.txt", sep="\t")
    
res.var$cos2
head(res.var$cos2)

fviz_pca_ind(pca_res)
corrplot(res.var$cos2[1:20,], is.corr = FALSE)
color = as.character(raw$color)

png = paste0("pca(raw)", "20190309.png")
png(filename=png, height=400, width=400, bg="white") # 
p<-ggplot(ind,aes(x = coord.Dim.1, y = coord.Dim.2))
p<-p+geom_point(size = 5, bg = color, col = "black", pch = 21) + theme + xlab(percentage[1]) + ylab(percentage[2])+ 
  ggtitle("PCA result")+
  geom_point(data = pca_som[which(pca_som$color == "yellow1"),], size = 5, shape = 21, col = "orange", fill = "yellow1")
p
garbage = dev.off() 
ind = as.data.frame(pca_res$ind)

#p<-ggplot(ind,aes(x = coord.Dim.1, y=coord.Dim.2, color = color))
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             strip.background=element_blank(),axis.text.x=element_text(colour="black", size = 18),
             axis.text.y=element_text(colour="black", size = 18),axis.ticks=element_line(colour="black"), axis.title =element_text(size = 18),
             plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
             plot.margin=unit(c(1,1,1,1),"line"))
percentage <- round(pca_res$eig[1:10,2], 2)
percentage <- paste( colnames(ind), "(", paste( as.character(percentage), "%", ")", sep="") )
cluster = factor(cluster_assignment)

color[color == "yellow1"] = "Bulk(A549,H1437)"
color[color == "blue3"] = "Single(H1437)"
color[color == "gray48"] = "Bulk(LUAD)"
color[color == "red"] = "Single(A549)"
color[color == "forestgreen"] =  "Single.Avg(A549,H1437)"

pca_som = cbind(ind, color)
pca_som$color = color
png(filename="cluster5_pca_som(20190304).png", height=400, width=450, bg="white")
p<-ggplot(pca_som,aes(x = coord.Dim.1, y = coord.Dim.2, shape = cluster))+
  scale_shape_manual(values = c(21,22,25,23,24))
p<-p+geom_point(size = 4, bg = color, col = "black")+ theme + xlab(percentage[1]) + ylab(percentage[2]) + 
  ggtitle("PCA-SOM result")
  #scale_fill_manual(values=c("#f9d9f5", "#fcfc55","#abd9fc", "#a6ddc8", "#f9d6bb"))
p
garbage = dev.off()
#color: scale_fill_manual(values=c("#f9d9f5", "#fcfc55","#abd9fc", "#a6ddc8", "#f9d6bb"))