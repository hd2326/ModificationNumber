files <- list.files(path = "./event/", pattern = ".tsv", full.names = T)
table <- list()
for (f in files){
  t <- read.delim(f, header = T, stringsAsFactors = F)
  n <- sapply(strsplit(unlist(strsplit(f, split = ".tsv")), split = "./event//"), function(x) x[2])
  if(nrow(t) > 0 & sum(is.na(t)) == 0) table[[n]] <- t[, c("position", "event_index", "event_level_mean", "model_kmer", "cigar_string", "cigar_basecalled_kmer", "basecall_index")]}
#event table

system("samtools view ./fastq/primer_extension.sorted.bam | awk '{ print $1, $3, $5}' > ./fastq/primer_extension_sam.txt")
sam <- read.delim("./fastq/primer_extension_sam.txt", header = F, stringsAsFactors = F, sep = " ", quote = "")
colnames(sam) <- c("id", "contig", "mapq")
sam <- sam[sam$id %in% names(table), ]
#sam file

fastq <- read.delim("./fastq/2017_06_05_CAM_ONT_BrdU_PCRs_primer_ext_basecalled_guppy3.1.5.fastq", header = F, stringsAsFactors = F, sep = " ", quote = "")
fastq <- data.frame(id=unlist(lapply(strsplit(fastq$V1[seq(1, nrow(fastq), by = 4)], split = "@"), function(x) x[2])),
                    seq=fastq$V1[seq(2, nrow(fastq), by = 4)],
                    q=fastq$V1[seq(4, nrow(fastq), by = 4)], stringsAsFactors = F)
fastq <- fastq[fastq$id %in% names(table), ]
q <- structure(lapply(fastq$q, function(x) utf8ToInt(x)), names=fastq$id)#q-score
seq <- structure(lapply(fastq$seq, function(x) unlist(strsplit(x, split = ""))), names = fastq$id)#sequence
#fastq file

contig <- structure(sam$contig, names=sam$id)#contig
strand <- structure(unlist(lapply(table, function(x){
  s <- ""
  if (sum(x$event_index != sort(x$event_index, decreasing = T)) == 0) s <- "R"
  else if (sum(x$event_index != sort(x$event_index, decreasing = F)) == 0) s <- "F"
  else s <- "NA"
  s
})), names=names(table))#strand distribution
mapq <- structure(sam$mapq, names=sam$id)#mapq score
map2 <- structure(duplicated(sam$id, fromLast = T) + duplicated(sam$id, fromLast = F) == 0, names=sam$id)#no sec/sup mapping
par(mfrow = c(2, 3), mar = c(5, 5, 5, 2))
barplot(table(contig), ylab = "#Reads", cex.names = 2, cex.axis = 2, cex.lab = 2, names.arg = c("B", "E", "F", "I", "T"), main = "Contig", cex.main = 2)
barplot(table(strand), ylab = "#Reads", cex.names = 2, cex.axis = 2, cex.lab = 2, main = "Strand", cex.main = 2)
hist(mapq, xlab = "", ylab = "#Reads", cex.axis = 2, cex.lab = 2, main = "MAPQ", cex.main = 2)
barplot(table(map2), ylab = "#Reads", cex.names = 2, cex.axis = 2, cex.lab = 2, main = "No Sec/Sup Mapping", cex.main = 2)
plot(density(unlist(lapply(q, function(x) mean(x-33)))), xlab = "", cex.axis = 2, cex.lab = 2, main = "Average Q-Score", cex.main = 2)
plot(density(unlist(lapply(seq, function(x) length(x)))), xlab = "", cex.axis = 2, cex.lab = 2, main = "Read Length", cex.main = 2)
#basic qc

ind <- Reduce(intersect, list(names(contig[contig == "T" | contig == "FdU" | contig == "BrdU" | contig == "IdU" | contig == "EdU"]),
                              names(strand[strand != "NA"]),
                              names(mapq[mapq == 60]),
                              names(map2[map2])))
table <- table[match(ind, names(table))]
sam <- sam[match(ind, sam$id), ]
fastq <- fastq[match(ind, fastq$id), ]
q <- q[match(ind, names(q))]
seq <- seq[match(ind, names(seq))]
strand <- strand[match(ind, names(strand))]
contig <- contig[match(ind, names(contig))]
for (i in 1:length(table)){
  table[[i]]$base <- seq[[i]][table[[i]]$basecall_index+1]
  table[[i]]$q <- unlist(lapply(table[[i]]$basecall_index, function(i, q) min(q[(i+1):(i-4)]), q=q[[i]]))
  table[[i]]$strand <- rep(strand[i], nrow(table[[i]]))
  table[[i]]$contig <- rep(contig[i], nrow(table[[i]]))}
save(table, file = "primer_extension_table.rda")
#filter/reform/save event table

load("primer_extension_table.rda")
table <- do.call(rbind, table)
table <- table[table$model_kmer != "NNNNNN", ]
ref <- "GAATTGGGCCCGCTCAGCAGACACAGAGCCTGAGCATCGCCGCGGAC"
position <- 0:41
strand <- structure(lapply(position, function(p, table) table$strand[table$position == p], table=table), names=position)
contig <- structure(lapply(position, function(p, table) table$contig[table$position == p], table=table), names=position)
sig <- structure(lapply(position, function(p, contig, table){
  structure(lapply(unique(unlist(contig)), function(c, p, table) table$event_level_mean[table$position == p & table$contig == c & table$strand == "F"], p=p, table=table), names = unique(unlist(contig)))}, contig=contig, table=table), names = position)
par(mfrow = c(6, 7), mar = c(5, 5, 4, 1))
for (p in as.character(position)){
  s <- sig[[p]]
  median <- structure(lapply(s, function(x) median(x)), names=names(s))
  mad <- structure(lapply(s, function(x) mad(x)), names=names(s))
  s <- structure(lapply(names(s), function(x, s, median, mad) s[[x]][s[[x]] > median[[x]]-2*mad[[x]] & s[[x]] < median[[x]]+2*mad[[x]]], s=s, median=median, mad=mad), names=names(s))
  plot(NULL, NULL, xlab = "pA", ylab = "Density", xlim = range(s), ylim = c(0, 0.3), cex.axis = 2, cex.lab = 2,
       main = paste(p, substr(ref, as.integer(p)+1, as.integer(p)+6)), col.main = (p %in% 25:36)+1, cex.main = 3)
  for (c in unique(unlist(contig))) if (length(s[[c]]) > 500) lines(density(s[[c]]), col = match(c, unique(unlist(contig))), lwd = 2)
  legend("topleft", legend = unlist(lapply(unique(unlist(contig)), function(c, s) paste(c, length(s[[c]]), sep = ":"), s=s)), text.col = 1:8, bty = "n", border = NA)}
#position distribution

library(MixtureInf)
load("primer_extension_table.rda")
table <- do.call(rbind, table)
table <- table[table$model_kmer != "NNNNNN" & table$strand == "F", c("position", "event_level_mean", "contig")]
base <- c("T", "EdU", "FdU", "BrdU", "IdU")
sig <- structure(lapply(base, function(b, table) table$event_level_mean[table$position == 27 & table$contig == b], table=table), names=base)
sig <- lapply(sig, function(x, median, mad) x[x > median-2*mad & x < median+2*mad], median=median(unlist(sig)), mad=mad(unlist(sig)))
sig <- sig[which(unlist(lapply(sig, function(x) length(x) > 500)))]
fit <- structure(lapply(1:9, function(i, sig) emtest.norm(sig, m0 = i), sig=unlist(sig)), names=1:9)
save(sig, fit, file = "primer_extension_fit_27.rda")
#position 27, GCCTGA, original paper showcase
sig <- structure(lapply(base, function(b, table) table$event_level_mean[table$position == 34 & table$contig == b], table=table), names=base)
sig <- lapply(sig, function(x, median, mad) x[x > median-2*mad & x < median+2*mad], median=median(unlist(sig)), mad=mad(unlist(sig)))
sig <- sig[which(unlist(lapply(sig, function(x) length(x) > 500)))]
fit <- structure(lapply(1:9, function(i, sig) emtest.norm(sig, m0 = i), sig=unlist(sig)), names=1:9)
save(sig, fit, file = "primer_extension_fit_34.rda")
#position 34, CATCGC
#mixture model

par(mfrow = c(2, 3), mar = c(5, 8, 3, 4))
load("primer_extension_fit_27.rda")
p <- unlist(lapply(fit, function(x) -log10(x$"P-values")))
plot(names(fit), p, type = "b", xlab = "#Components", ylab = "-log10(p-value)",
     col = c(1, 1, 1, 1, 1, 1, 2, 1, 1), cex = 2, cex.axis = 2, cex.lab = 2, main = "")
#determine #components
barplot(fit$"7"[[1]][1, ], xlab = "Component", ylab = "Propotion", col = c(2, 2, 8, 8, 8, 8, 2),
        names.arg = 1:7, cex.axis = 2, cex.lab = 2, cex.names = 2,  main = "GCCTGA", cex.main = 4)
#filter components
p <- matrix(NA, length(sig), length(sig), dimnames = list(names(sig), names(sig)))
for (i in 1:length(sig)) for (j in 1:length(sig)) if (i > j) p[i, j] <- wilcox.test(sig[[i]], sig[[j]])$p.value
image(-log10(p), xaxt = "n", yaxt = "n", col = colorRampPalette(c("Grey", "Red"))(100), main = "")
axis(side = 1, at = seq(0, 1, length.out = nrow(p)), labels = rownames(p), tick = F, las = 2, cex.axis = 2)
axis(side = 2, at = seq(0, 1, length.out = ncol(p)), labels = colnames(p), tick = F, las = 2, cex.axis = 2)
for (i in 1:nrow(p)) for (j in 1:ncol(p)) if (i > j) text(seq(0, 1, length.out = nrow(p))[i], seq(0, 1, length.out = ncol(p))[j], labels = signif(-log10(p[i, j]), 2), cex = 1.5)
legend("topleft", legend = "-log10(p-value)", bty = "n", border = NA, cex = 2)
#signal comparison

load("primer_extension_fit_34.rda")
p <- unlist(lapply(fit, function(x) -log10(x$"P-values")))
plot(names(fit), p, type = "b", xlab = "#Components", ylab = "-log10(p-value)",
     col = c(1, 1, 1, 1, 1, 1, 2, 1, 1), cex = 2, cex.axis = 2, cex.lab = 2, main = "")
#determine #components
barplot(fit$"7"[[1]][1, ], xlab = "Component", ylab = "Propotion", col = c(2, 2, 8, 8, 8, 8, 2),
        names.arg = 1:7, cex.axis = 2, cex.lab = 2, cex.names = 2,  main = "CATCGC", cex.main = 4)
#filter components
p <- matrix(NA, length(sig), length(sig), dimnames = list(names(sig), names(sig)))
for (i in 1:length(sig)) for (j in 1:length(sig)) if (i > j) p[i, j] <- wilcox.test(sig[[i]], sig[[j]])$p.value
image(-log10(p), xaxt = "n", yaxt = "n", col = colorRampPalette(c("Grey", "Red"))(100), main = "")
axis(side = 1, at = seq(0, 1, length.out = nrow(p)), labels = rownames(p), tick = F, las = 2, cex.axis = 2)
axis(side = 2, at = seq(0, 1, length.out = ncol(p)), labels = colnames(p), tick = F, las = 2, cex.axis = 2)
for (i in 1:nrow(p)) for (j in 1:ncol(p)) if (i > j) text(seq(0, 1, length.out = nrow(p))[i], seq(0, 1, length.out = ncol(p))[j], labels = signif(-log10(p[i, j]), 2), cex = 1.5)
legend("topleft", legend = "-log10(p-value)", bty = "n", border = NA, cex = 2)
#signal comparison
#determine optimal number of components                        
 
library(mixtools)
layout(matrix(c(1, 1, 1, 5, 5, 5,
                2, 3, 4, 6, 7, 8), nrow = 2, ncol = 6, byrow = T), widths = rep(2, 6), heights = c(6, 2))
load("primer_extension_fit_27.rda")
par(mar = c(6, 6, 8, 6))
plot(density(unlist(sig)), ylim = c(0, 0.25), xlab = "pA", cex.axis = 2, cex.lab = 2, lwd = 3, main = "GCCTGA", cex.main = 5)
lines(density(unlist(rnormmix(50000, lambda = fit$"7"[[1]][1, ], mu = fit$"7"[[1]][2, ], sigma = sqrt(fit$"7"[[1]][3, ])))), lwd = 3, col = 2)
lines(density(unlist(rnormmix(50000, lambda = fit$"7"[[1]][1, 3:6], mu = fit$"7"[[1]][2, 3:6], sigma = sqrt(fit$"7"[[1]][3, 3:6])))), lwd = 3, col = 3)
legend("topright", legend = c("empirical", "fit", "fit refine"), col = 1:3, lty = 1, cex = 2, lwd = 3, bty = "n", border = NA)
par(new = T)
plot(NA, NA, axes = F, xlab = "", ylab = "", ylim = c(0, 0.4), xlim = range(unlist(sig)))
axis(side = 4, cex.axis = 2)
for (i in 1:length(sig)) lines(density(sig[[i]]), lty = 2, lwd = 3, col = i+3)
legend("topleft", legend = names(sig), col = 4:8, lty = 2, cex = 2, lwd = 3,  bty = "n", border = NA)
#signal distribution
i1 <- order(fit$"7"[[1]][2, 3:6])
i2 <- order(c(median(sig$T), median(sig$EdU), median(sig$FdU), median(c(sig$BrdU, sig$IdU))))
par(mar = c(5, 5, 1, 1))
plot(fit$"7"[[1]][2, 3:6][i1], c(median(sig$T), median(sig$EdU), median(sig$FdU), median(c(sig$BrdU, sig$IdU)))[i2],
     xlim = c(90, 97), ylim = c(90, 97), pch = 16, xlab = "fitted", ylab = "median", cex.axis = 2, cex.lab = 2)
text(fit$"7"[[1]][2, 3:6][i1], c(median(sig$T), median(sig$EdU), median(sig$FdU), median(c(sig$BrdU, sig$IdU)))[i2], labels = c("T", "EdU", "FdU", "Br/IdU")[i2], cex = 1.5)
abline(0, 1, NULL, NULL, col = 2, lty = 4, lwd = 2)
#mean/median
plot(sqrt(fit$"7"[[1]][3, 3:6])[i1], c(mad(sig$T), mad(sig$EdU), mad(sig$FdU), mad(c(sig$BrdU, sig$IdU)))[i2],
     xlim = c(0, 2.5), ylim = c(0, 2.5), pch = 16, xlab = "fitted", ylab = "mad", cex.axis = 1.5, cex.lab = 2)
text(sqrt(fit$"7"[[1]][3, 3:6])[i1], c(mad(sig$T), mad(sig$EdU), mad(sig$FdU), mad(c(sig$BrdU, sig$IdU)))[i2], labels = c("T", "IdU", "FdU", "Br/EdU")[i2], cex = 1.5)
abline(0, 1, NULL, NULL, col = 2, lty = 4, lwd = 2)
#sd/mad
plot((fit$"7"[[1]][1, 3:6]/sum(fit$"7"[[1]][1, 3:6]))[i1], (c(length(sig$T), length(sig$EdU), length(sig$FdU), length(c(sig$BrdU, sig$IdU)))/length(unlist(sig)))[i2],
     xlim = c(0.05, 0.5), ylim = c(0.05, 0.5), pch = 16, xlab = "fitted", ylab = "porpotion", cex.axis = 1.5, cex.lab = 2)
text(fit$"7"[[1]][1, 3:6][i1], (c(length(sig$T), length(sig$EdU), length(sig$FdU), length(c(sig$BrdU, sig$IdU)))/length(unlist(sig)))[i2], labels = c("T", "EdU", "FdU", "Br/IdU")[i2], cex = 1.5)
abline(0, 1, NULL, NULL, col = 2, lty = 4, lwd = 2)
#alpha

load("primer_extension_fit_34.rda")
par(mar = c(6, 6, 8, 6))
plot(density(unlist(sig)), ylim = c(0, 0.15), xlab = "pA", cex.axis = 2, cex.lab = 2, lwd = 3, main = "CATCGC", cex.main = 5)
lines(density(unlist(rnormmix(50000, lambda = fit$"7"[[1]][1, ], mu = fit$"7"[[1]][2, ], sigma = sqrt(fit$"7"[[1]][3, ])))), lwd = 3, col = 2)
lines(density(unlist(rnormmix(50000, lambda = fit$"7"[[1]][1, 3:6], mu = fit$"7"[[1]][2, 3:6], sigma = sqrt(fit$"7"[[1]][3, 3:6])))), lwd = 3, col = 3)
legend("topright", legend = c("empirical", "fit", "fit refine"), col = 1:3, lty = 1, cex = 2, lwd = 3, bty = "n", border = NA)
par(new = T)
plot(NA, NA, axes = F, xlab = "", ylab = "", ylim = c(0, 0.2), xlim = range(unlist(sig)))
axis(side = 4, cex.axis = 2)
for (i in 1:length(sig)) lines(density(sig[[i]]), lty = 2, lwd = 3, col = i+3)
legend("topleft", legend = names(sig), col = 4:8, lty = 2, cex = 2, lwd = 3,  bty = "n", border = NA)
#signal distribution
i1 <- order(fit$"7"[[1]][2, 3:6])
i2 <- order(c(median(sig$T), median(sig$EdU), median(sig$FdU), median(c(sig$BrdU, sig$IdU))))
par(mar = c(5, 5, 1, 1))
plot(fit$"7"[[1]][2, 3:6][i1], c(median(sig$T), median(sig$EdU), median(sig$FdU), median(c(sig$BrdU, sig$IdU)))[i2],
     xlim = c(112, 128), ylim = c(112, 128), pch = 16, xlab = "fitted", ylab = "median", cex.axis = 1.5, cex.lab = 2)
text(fit$"7"[[1]][2, 3:6][i1], c(median(sig$T), median(sig$EdU), median(sig$FdU), median(c(sig$BrdU, sig$IdU)))[i2], labels = c("T", "EdU", "FdU", "Br/IdU")[i2], cex = 1.5)
abline(0, 1, NULL, NULL, col = 2, lty = 4, lwd = 2)
#mean/median
plot(sqrt(fit$"7"[[1]][3, 3:6])[i1], c(mad(sig$T), mad(sig$EdU), mad(sig$FdU), mad(c(sig$BrdU, sig$IdU)))[i2],
     xlim = c(0.5, 3.5), ylim = c(0.5, 3.5), pch = 16, xlab = "fitted", ylab = "mad", cex.axis = 1.5, cex.lab = 2)
text(sqrt(fit$"7"[[1]][3, 3:6])[i1], c(mad(sig$T), mad(sig$EdU), mad(sig$FdU), mad(c(sig$BrdU, sig$IdU)))[i2], labels = c("T", "EdU", "FdU", "Br/IdU")[i2], cex = 1.5)
abline(0, 1, NULL, NULL, col = 2, lty = 4, lwd = 2)
#sd/mad
plot((fit$"7"[[1]][1, 3:6]/sum(fit$"7"[[1]][1, 3:6]))[i1], (c(length(sig$T), length(sig$EdU), length(sig$FdU), length(c(sig$BrdU, sig$IdU)))/length(unlist(sig)))[i2],
     xlim = c(0.1, 0.4), ylim = c(0.1, 0.4), pch = 16, xlab = "fitted", ylab = "propotion", cex.axis = 1.5, cex.lab = 2)
text(fit$"7"[[1]][1, 3:6][i1], (c(length(sig$T), length(sig$EdU), length(sig$FdU), length(c(sig$BrdU, sig$IdU)))/length(unlist(sig)))[i2], labels = c("T", "EdU", "FdU", "Br/IdU")[i2], cex = 1.5)
abline(0, 1, NULL, NULL, col = 2, lty = 4, lwd = 2)
#alpha
#visualize fitting results  
  
load("primer_extension_table.rda")
read <- unlist(lapply(table, function(x) sum(25:36 %in% x$position) == 12))
contig <- unlist(lapply(table[read], function(x) x$contig[1]))
sig_25_30 <- lapply(table[read], function(x) structure(unlist(lapply(25:30, function(p, x) mean(x$event_level_mean[x$position == p]), x=x)), names=25:30))
sig_25_30 <- do.call(rbind, sig_25_30)
sig_31_36 <- lapply(table[read], function(x) structure(unlist(lapply(31:36, function(p, x) mean(x$event_level_mean[x$position == p]), x=x)), names=31:36))
sig_31_36 <- do.call(rbind, sig_31_36)
sig_25_36 <- cbind(sig_25_30, sig_31_36)
#clustering reads
                                                                     
library(dendextend)
library(cluster)
layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = F), heights = c(4, 1))
pam <- pam(dist(sig_25_36), k = 2)$clustering
hclust <- hclust(dist(sig_25_36[pam == 1, ]), method = "ward.D")
cluster <- cutree(hclust, 5)
prop <- lapply(1:5, function(i, cluster, contig) table(contig[cluster == i]), contig=contig[pam == 1], cluster=cluster)
prop <- apply(do.call(cbind, prop), 2, prop.table)[, c(4, 1, 2, 3, 5)]
colnames(prop) <- paste("Cluster", 1:5, sep = "")
dend <- as.dendrogram(hclust)
labels(dend) <- ""
dend <- color_branches(dend, k = 5, col = rainbow(5))
par(mar = c(1, 1.5, 4, 1.5))
plot(dend, axes = F, main = "Pos 25-36", cex.main = 2)
legend("topleft", legend = c("BrdU", "EdU", "FdU", "IdU", "T"), fill = 1:5, cex = 1.5, bty = "n", border = NA)
legend("topright", legend = paste("Cluster", 1:5, sep = ""), fill = rainbow(5), cex = 1.5, bty = "n", border = NA)
par(mar = c(1, 1, 0, 1))
barplot(rep(1, sum(pam == 1)), space = 0, border = NA, col = as.factor(contig)[pam == 1][hclust$order], axes = F)
par(mar = c(6, 2, 4, 8))
image(prop, axes = F, main = "")
axis(side = 1, at = seq(0, 1, length.out = nrow(prop)), labels = rownames(prop), tick = F, las = 2, cex.axis = 1.5)
axis(side = 4, at = seq(0, 1, length.out = ncol(prop)), labels = colnames(prop), tick = F, las = 2, cex.axis = 1.5)
for (i in 1:nrow(prop)) for (j in 1:ncol(prop)) text(seq(0, 1, length.out = nrow(prop))[i],
                                                     seq(0, 1, length.out = ncol(prop))[j], labels = signif(prop[i, j], 2), cex = 1.5)
dev.off()
#clustering visualization
