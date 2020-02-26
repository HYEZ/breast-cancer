#calculate som error

library("kohonen")
library("tempR")
setwd("C:/Users/user/Desktop/single cell/pca/output/2199 genes")

raw = as.data.frame(read.table("C:/Users/user/Desktop/single cell/DATA/raw(2199).txt", sep = "\t",
                               header = T, stringsAsFactors = FALSE))

name = raw[,1]
data = as.data.frame(raw[, 3:ncol(raw)], stringsAsFactors = FALSE)
rownames(data) = name
data = as.data.frame(sapply(data, as.numeric))
data = as.data.frame(data)
data[195,] = apply(data[1:74,], 2, mean, na.rm = T)
data[196,] = apply(data[75:121,], 2, mean, na.rm = T)
raw[195:196,2] = "forestgreen"

data = as.matrix(data)
color = unlist(raw[,2])

#SOM
i =  1 
for(i in 2:194){
  j = 1
while((i*j) <= 194){
  som_grid = somgrid(i, j, topo = "hexagonal")
  som_model <- som(data, 
                   grid=som_grid, 
                   rlen=500, 
                   alpha=c(0.01,0.01), 
                   keep.data = TRUE,
                   maxNA.fraction = 5)
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

write.table(out, "out2.txt", sep = "\t", row.names = F)
