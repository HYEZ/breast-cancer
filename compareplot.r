setwd("C:/Users/user/Desktop/single cell/pca/data/output")
data = read.table("C:/Users/user/Desktop/single cell/pca/data/rescale.txt",
                  sep = "\t", header = T)

A549_single = data[,1:75]
H1437_single = data[,76:122]
LUAD = data[,123:ncol(data)]

A549_bulk = LUAD[, grep("A549", colnames(LUAD))]
row.names(A549_single) = A549_single$gene
row.names(H1437_single) = A549_single$gene
row.names(LUAD) = A549_single$gene
A549_single = A549_single[,-1]

H1437_bulk = LUAD[, grep("H1437",  colnames(LUAD))]


A549_log = log2(A549_single)
is.na(A549_log) = sapply(A549_log, is.infinite)

H1437_log = log2(H1437_single)
is.na(H1437_log) = sapply(H1437_log, is.infinite)

A549_avg = apply(A549_log, 1, function(x){mean(x, na.rm = T)})
H1437_avg = apply(H1437_log, 1, function(x){mean(x, na.rm = T)})

A549_bulk_log = log2(A549_bulk)
is.na(A549_bulk_log) = sapply(A549_bulk_log, is.infinite)

H1437_bulk_log = log2(H1437_bulk)
is.na(H1437_bulk_log) = sapply(H1437_bulk_log, is.infinite)

A549_log_cor = cor(A549_avg, A549_bulk_log, use = "complete.obs", method = "pearson")
H1437_log_cor = cor(H1437_avg, H1437_bulk_log, use = "complete.obs", method = "pearson")



#draw compare plot
png("rescale_log2_A549(single vs bulk).png", width = 300, height = 300)
plot(x = A549_avg, y = A549_bulk_log,  main = "A549", xlab = "Avg of single cell",
     ylab = "CCLE RNA seq", xlim = c(-25,0), ylim = c(-25,0))
abline(lm(A549_bulk_log ~ A549_avg))
mtext(paste0("PCC: ", round(A549_log_cor,2)), side = 1, line = -2, adj = 0.98)
trash = dev.off()

png("rescale_log2_H1437(single vs bulk).png", width = 300, height = 300)
plot(x = H1437_avg, y = H1437_bulk_log, main = "H1437", xlab = "Avg of single cell",
     ylab = "CCLE RNA seq", xlim = c(-25,0), ylim = c(-25,0)) 
abline(lm(H1437_bulk_log ~ H1437_avg))
mtext(paste0("PCC: ", round(H1437_log_cor,2)), side = 1, line = -2, adj = 0.98)
trash = dev.off()

#Export dataset
A549_sum = cbind(A549_avg, A549_bulk_log)
H1437_sum = cbind(H1437_avg, H1437_bulk_log)

A549_sum = as.data.frame(A549_sum)
A549_sum = na.omit(A549_sum)
H1437_sum = as.data.frame(H1437_sum)
H1437_sum = na.omit(H1437_sum)

#write.csv(A549_sum, "A549(single_avg vs bulk).csv")
#write.csv(H1437_sum, "H1437(single_avg vs bulk).csv")
A549_cutoff = lm(A549_bulk_log-3 ~ A549_avg)
A549_sub = A549_sum[A549_sum$A549_bulk_log > predict(A549_cutoff, newdata = data.frame(A549_avg = A549_sum$A549_avg)),]
A549_filtering = A549_sum[A549_sum$A549_bulk_log <= predict(A549_cutoff, newdata = data.frame(A549_avg = A549_sum$A549_avg)),]
A549_sub_cor = cor(A549_sub$A549_avg, A549_sub$A549_bulk_log, use = "complete.obs", method = "pearson")

H1437_cutoff = lm(H1437_bulk_log-3 ~ H1437_avg)
H1437_sub = H1437_sum[H1437_sum$H1437_bulk_log > predict(H1437_cutoff, newdata = data.frame(H1437_avg = H1437_sum$H1437_avg)),]
H1437_filtering = H1437_sum[H1437_sum$H1437_bulk_log <= predict(H1437_cutoff, newdata = data.frame(H1437_avg = H1437_sum$H1437_avg)),]
H1437_sub_cor = cor(H1437_sub$H1437_avg, H1437_sub$H1437_bulk_log, use = "complete.obs", method = "pearson")

#draw compare plot + filtering
png("rescale_log2_A549(single vs bulk)_color.png", width = 300, height = 300)
plot(x = A549_avg, y = A549_bulk_log,  main = "A549", xlab = "Avg of single cell",
     ylab = "CCLE RNA seq", xlim = c(-25,0), ylim = c(-25,0), col = "grey48") + points(x = A549_sub$A549_avg, y = A549_sub$A549_bulk_log, col = "black")
abline(lm(A549_bulk_log -3 ~ A549_avg), col = "grey48")
abline(lm(A549_bulk_log ~ A549_avg))

#mtext(paste0("PCC: ", round(A549_log_cor,2)), side = 1, line = -2, adj = 0.98)
trash = dev.off()

png("rescale_log2_H1437(single vs bulk)_color.png", width = 300, height = 300)
plot(x = H1437_avg, y = H1437_bulk_log, main = "H1437", xlab = "Avg of single cell",
     ylab = "CCLE RNA seq", xlim = c(-25,0), ylim = c(-25,0), col = "grey48") + points(x = H1437_sub$H1437_avg, y = H1437_sub$H1437_bulk_log, col = "black")
abline(lm(H1437_bulk_log -3 ~ H1437_avg), col = "grey48")
abline(lm(H1437_bulk_log ~ H1437_avg))
#mtext(paste0("PCC: ", round(H1437_log_cor,2)), side = 1, line = -2, adj = 0.98)
trash = dev.off()

common = merge(A549_sub, H1437_sub, by = "row.names")
write.table(common, "common_filtering.txt", sep = "\t")
write.table(A549_filtering,"cut(A549).txt", sep = "\t")
write.table(A549_sub, "after filtering(A549).txt", sep = "\t")

write.table(H1437_filtering, "cut(H1437).txt", sep = "\t")
write.table(H1437_sub, "after filtering(H1437).txt", sep = "\t")

#after filtering plot

png("rescale_log2_A549(single vs bulk)_filtering.png", width = 300, height = 300)
plot(x = A549_sub$A549_avg, y = A549_sub$A549_bulk_log, main = "A549", xlab = "Avg of single cell",
     ylab = "CCLE RNA seq", xlim = c(-25,0), ylim = c(-25,0))
abline(lm(A549_sub$A549_bulk_log ~ A549_sub$A549_avg))
mtext(paste0("PCC: ", round(A549_sub_cor,2)), side = 1, line = -2, adj = 0.98)
trash = dev.off()


png("rescale_log2_H1437(single vs bulk)_filtering.png", width = 300, height = 300)
plot(x = H1437_sub$H1437_avg, y = H1437_sub$H1437_bulk_log, main = "H1437", xlab = "Avg of single cell",
     ylab = "CCLE RNA seq", xlim = c(-25,0), ylim = c(-25,0))
abline(lm(H1437_sub$H1437_bulk_log ~ H1437_sub$H1437_avg))
mtext(paste0("PCC: ", round(H1437_sub_cor,2)), side = 1, line = -2, adj = 0.98)
trash = dev.off()


common_cor_a549 = cor(common$A549_avg, common$A549_bulk_log, use = "complete.obs", method = "pearson")
common_cor_h1437 = cor(common$H1437_avg, common$H1437_bulk_log, use = "complete.obs", method = "pearson")


write.table(H1437_filtering, "H1437_filter.txt", sep = "\t")

# 
# A549_cor = cor(A549_avg, A549_bulk, use = "complete.obs", method = "pearson")
# H1437_cor = cor(H1437_avg, H1437_bulk, use = "complete.obs", method = "pearson")
# 
# png("RAW_A549(single vs bulk).png", width = 300, height = 300)
# plot(x = A549_avg, y = A549_bulk, main = "A549 single cell avg vs bulk cell" )
# abline(lm(A549_bulk ~ A549_avg))
# mtext(paste0("PCC: ", round(A549_cor,2)), side = 1, line = -2, adj = 0.98)
# trash = dev.off()
# 
# png("RAW_H1437(single vs bulk).png", width = 300, height = 300)
# plot(x = H1437_avg, y = H1437_bulk, main = "H1437 single cell avg vs bulk cell")
# abline(lm(H1437_bulk ~ H1437_avg))
# mtext(paste0("PCC: ", round(H1437_cor,2)), side = 1, line = -2, adj = 0.98)
# trash = dev.off()
