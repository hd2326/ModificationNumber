files <- list.files(pattern = "_table.rda")
table <- structure(lapply(files, function(f){
  load(f)
  table <- do.call(rbind, table[1:min(5000, length(table))])
  table[table$model_kmer != "NNNNNN", ]
}), names=unlist(lapply(strsplit(files, split = "_"), function(x) x[1])))
ref <- read.delim("ec_16S.fa", stringsAsFactors=FALSE)
ref <- paste(ref[, 1], collapse="")
position <- 0:(nchar(ref)-5)
sig <- structure(lapply(position, function(p, table) structure(lapply(table, function(x, p) x$event_level_mean[x$position == p], p=p), names=names(table)), table=table), names=position)

par(mfrow = c(6, 7), mar = c(5, 5, 4, 1))
for (p in as.character(position)){
  s <- sig[[p]]
  plot(NULL, NULL, xlab = "pA", ylab = "Density", xlim = range(unlist(s)), ylim = c(0, 0.5), cex.axis = 2, cex.lab = 2,
       main = paste(p, gsub("T", "U", substr(ref, as.integer(p)+1, as.integer(p)+5))), cex.main = 3, col.main = p%in%c(511:515, 522:526)+1)
  if (sum(unlist(lapply(s, function(ss) length(ss) == 0))) == 0){
    median <- structure(lapply(s, function(x) median(x)), names=names(s))
    mad <- structure(lapply(s, function(x) mad(x)), names=names(s))
    s <- structure(lapply(names(s), function(x, s, median, mad) s[[x]][s[[x]] > median[[x]]-2*mad[[x]] & s[[x]] < median[[x]]+2*mad[[x]]], s=s, median=median, mad=mad), names=names(s))
    for (i in 1:length(s)) if (length(s[[i]]) > 500) lines(density(s[[i]]), col = i, lwd = 2)}
  legend("topleft", legend = unlist(lapply(1:length(s), function(i, s) paste(names(s)[i], length(s[[i]]), sep = ":"), s=s)), text.col = 1:8, bty = "n", border = NA)}
#position distribution

library(MixtureInf)
s <- sig$"515"
median <- structure(lapply(s, function(x) median(x)), names=names(s))
mad <- structure(lapply(s, function(x) mad(x)), names=names(s))
s <- structure(lapply(names(s), function(x, s, median, mad) s[[x]][s[[x]] > median[[x]]-2*mad[[x]] & s[[x]] < median[[x]]+2*mad[[x]]], s=s, median=median, mad=mad), names=names(s))
fit <- structure(lapply(1:4, function(i, s) emtest.norm(s, m0 = i), s=unlist(s)), names=1:4)
save(s, fit, file = "rrna_fit_515.rda")
#position 515, Tgcca
s <- sig$"523"
median <- structure(lapply(s, function(x) median(x)), names=names(s))
mad <- structure(lapply(s, function(x) mad(x)), names=names(s))
s <- structure(lapply(names(s), function(x, s, median, mad) s[[x]][s[[x]] > median[[x]]-2*mad[[x]] & s[[x]] < median[[x]]+2*mad[[x]]], s=s, median=median, mad=mad), names=names(s))
fit <- structure(lapply(1:4, function(i, s) emtest.norm(s, m0 = i), s=unlist(s)), names=1:4)
save(s, fit, file = "rrna_fit_523.rda")
#position 523, gccGc

par(mfrow = c(1, 2), mar = c(5, 5, 5, 5))
load("rrna_fit_515.rda")
fit <- fit[1:4]
p <- unlist(lapply(fit, function(x) -log10(x$"P-values")))
p[!is.finite(p)] <- -log10(.Machine$double.xmin)
plot(names(fit), p, type = "b", xlab = "#Components", ylab = "-log10(p-value)",
     col = c(1, 2, 1, 1), cex = 2, cex.axis = 2, cex.lab = 2, main = "UGCCA", cex.main = 3)

load("rrna_fit_523.rda")
fit <- fit[1:4]
p <- unlist(lapply(fit, function(x) -log10(x$"P-values")))
p[!is.finite(p)] <- -log10(.Machine$double.xmin)
plot(names(fit), p, type = "b", xlab = "#Components", ylab = "-log10(p-value)",
     col = c(1, 2, 1, 1), cex = 2, cex.axis = 2, cex.lab = 2, main = "GCCGC", cex.main = 3)
#determine optimal number of components

par(mfrow = c(1, 2), mar = c(5, 5, 5, 2))
load("rrna_fit_515.rda")
plot(density(unlist(s)), xlab = "pA", ylab = "Density", xlim = range(unlist(s)), ylim = c(0, 0.3), lwd = 3, cex.axis = 2, cex.lab = 2, main = "UGCCA", cex.main = 3)
for (i in 1:length(s)) lines(density(s[[i]]), col = i+1, lwd = 3)
for (i in 1:2) lines(seq(0, 200, length.out = 1000), dnorm(seq(0, 200, length.out = 1000), mean = fit$"2"[[1]]["mean", i], sd = sqrt(fit$"2"[[1]]["variance", i])), col = 4+i, lty = 2, lwd = 3)
legend("topleft", legend = c("All", unlist(lapply(1:length(s), function(i, s) paste(names(s)[i], ":", length(s[[i]]), " ", signif(length(s[[i]])/length(unlist(s))*100, 2), "%", sep = ""), s=s))), text.col = 1:4, bty = "n", border = NA)
legend("topright", legend = c("C1, 60%", "C2, 40%"), text.col = 5:6, bty = "n", border = NA)
legend("bottomright", legend = c("empirical", "fitted"), lty = c(1, 2), lwd = 3, bty = "n", border = NA)
load("rrna_fit_523.rda")
plot(density(unlist(s)), xlab = "pA", ylab = "Density", xlim = range(unlist(s)), ylim = c(0, 0.4), lwd = 3, cex.axis = 2, cex.lab = 2, main = "GCCGC", cex.main = 3)
for (i in 1:length(s)) lines(density(s[[i]]), col = i+1, lwd = 3)
for (i in 1:2) lines(seq(0, 200, length.out = 1000), dnorm(seq(0, 200, length.out = 1000), mean = fit$"2"[[1]]["mean", i], sd = sqrt(fit$"2"[[1]]["variance", i])), col = 4+i, lty = 2, lwd = 3)
legend("topleft", legend = c("All", unlist(lapply(1:length(s), function(i, s) paste(names(s)[i], ":", length(s[[i]]), " ", signif(length(s[[i]])/length(unlist(s))*100, 2), "%", sep = ""), s=s))), text.col = 1:4, bty = "n", border = NA)
legend("topright", legend = c("C1, 30%", "C2, 70%"), text.col = 5:6, bty = "n", border = NA)
legend("bottomright", legend = c("empirical", "fitted"), lty = c(1, 2), lwd = 3, bty = "n", border = NA)
#visualize fitting results                   
                   
files <- list.files(pattern = "_table.rda")
table <- structure(lapply(files, function(f){
  load(f)
  table <- table[1:min(5000, length(table))]
  lapply(table, function(t) t[t$model_kmer != "NNNNNN", ])
}), names=unlist(lapply(strsplit(files, split = "_"), function(x) x[1])))
table <- do.call(c, table)
read <- unlist(lapply(table, function(x) sum(c(511:515, 522:526) %in% x$position) == 10))
contig <- sapply(strsplit(names(read)[read], split = "[.]"), function(x) x[1])
sig_511_515 <- lapply(table[read], function(x) structure(unlist(lapply(511:515, function(p, x) mean(x$event_level_mean[x$position == p]), x=x)), names=511:515))
sig_511_515 <- do.call(rbind, sig_511_515)
sig_511_515[is.na(sig_511_515)] <- 0
sig_522_526 <- lapply(table[read], function(x) structure(unlist(lapply(522:526, function(p, x) mean(x$event_level_mean[x$position == p]), x=x)), names=522:526))
sig_522_526 <- do.call(rbind, sig_522_526)
sig_522_526[is.na(sig_522_526)] <- 0
sig_all <- cbind(sig_511_515, sig_522_526)
#clustering reads

library(dendextend)
library(cluster)
layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = F), heights = c(4, 1))
pca <- prcomp(sig_all, scale. = T, center = T)
hclust <- hclust(dist(pca$x), method = "ward.D2")
cluster <- cutree(hclust, 3)
prop <- lapply(1:3, function(i, cluster, contig) table(contig[cluster == i]), contig=contig, cluster=cluster)
prop <- apply(do.call(cbind, prop), 2, prop.table)[, c(3, 2, 1)]
colnames(prop) <- paste("Cluster", 1:3, sep = "")
dend <- as.dendrogram(hclust)
labels(dend) <- ""
dend <- color_branches(dend, k = 3, col = rainbow(3))
par(mar = c(1, 1.5, 4, 1.5))
plot(dend, axes = F, main = "Pos 511-515, 522-526", cex.main = 2)
legend("topleft", legend = c("m7G", "Native", "Psi516"), fill = 1:3, cex = 1.5, bty = "n", border = NA)
legend("topright", legend = paste("Cluster", 1:3, sep = ""), fill = rainbow(3), cex = 1.5, bty = "n", border = NA)
par(mar = c(1, 1, 0, 1))
barplot(rep(1, nrow(sig_all)), space = 0, border = NA, col = as.factor(contig)[hclust$order], axes = F)
par(mar = c(6, 2, 4, 8))
image(prop, axes = F, main = "")
axis(side = 1, at = seq(0, 1, length.out = nrow(prop)), labels = rownames(prop), tick = F, las = 2, cex.axis = 1.5)
axis(side = 4, at = seq(0, 1, length.out = ncol(prop)), labels = colnames(prop), tick = F, las = 2, cex.axis = 1.5)
for (i in 1:nrow(prop)) for (j in 1:ncol(prop)) text(seq(0, 1, length.out = nrow(prop))[i],
                                                     seq(0, 1, length.out = ncol(prop))[j], labels = signif(prop[i, j], 2), cex = 1.5)
dev.off()
#visualize clustering results

sig <- sig$"515"
median <- structure(lapply(sig, function(x) median(x)), names=names(sig))
mad <- structure(lapply(sig, function(x) mad(x)), names=names(sig))
sig <- structure(lapply(names(sig), function(x, sig, median, mad) sig[[x]][sig[[x]] > median[[x]]-2*mad[[x]] & sig[[x]] < median[[x]]+2*mad[[x]]], sig=sig, median=median, mad=mad), names=names(sig))
fit <- structure(lapply(c(100, 1000, 2000), function(s, sig){
  structure(lapply(c(0.01, 0.05, 0.1, 0.25, 0.5), function(p, s, sig){
    structure(lapply(1:10, function(i, p, s, sig){
      s1 <- sample(sig$"Psi516", s*p)
      s2 <- sample(sig$"native", s*(1-p))
      fit <- structure(lapply(1:3, function(i, s) emtest.norm(s, m0 = i), s=c(s1, s2)), names=1:3)
      list(psi516=s1, native=s2, fit=fit)
    }, p=p, s=s, sig=sig), names = 1:10)
  }, s=s, sig=sig), names=c(0.01, 0.05, 0.1, 0.25, 0.5))
}, sig=sig), names=c(100, 1000, 2000))
save(fit, file = "rrna_fit_down-sample.rda")
#down-sample analysis

alpha <- lapply(fit, function(f) t(do.call(rbind, lapply(f, function(ff) do.call(cbind, lapply(ff, function(fff) fff$"fit"$"2"[[1]][1, ]))[2, ]))))
mean1 <- lapply(fit, function(f) t(do.call(rbind, lapply(f, function(ff) do.call(cbind, lapply(ff, function(fff) fff$"fit"$"2"[[1]][2, ]))[1, ]))))
mean2 <- lapply(fit, function(f) t(do.call(rbind, lapply(f, function(ff) do.call(cbind, lapply(ff, function(fff) fff$"fit"$"2"[[1]][2, ]))[2, ]))))
par(mfrow = c(2, 3), mar = c(5, 5, 5, 1))
for (i in 1:3){
  boxplot(alpha[[i]], ylim = c(0, 1), ylab = "Fraction", cex.lab = 2, axes = F, main = paste(names(alpha)[i], "Total Reads"), cex.main = 1.5)
  axis(side = 1, at = 1:5, labels = c(0.01, 0.05, 0.1, 0.25, 0.5), cex.axis = 1.5, las = 2)
  axis(side = 2, at = seq(0, 1, length.out = 6), labels = seq(0, 1, length.out = 6), cex.axis = 1.5)
  for (j in c(0.01, 0.05, 0.1, 0.25, 0.5)) abline(j, 0, NULL, NULL, col = 2, lty = 4)}
for (i in 1:3){
  boxplot(mean1[[i]], ylim = c(90, 104), ylab = "pA", cex.lab = 2, axes = F, main = paste(names(alpha)[i], "Total Reads"), cex.main = 1.5)
  boxplot(mean2[[i]], ylim = c(90, 104), axes = F, add = T, border = 4, main = "")
  axis(side = 1, at = 1:5, labels = c(0.01, 0.05, 0.1, 0.25, 0.5), cex.axis = 1.5)
  axis(side = 2, at = seq(90, 104, length.out = 8), labels = seq(90, 104, length.out = 8), cex.axis = 1.5, las = 2)
  abline(median(sig$"native"), 0, NULL, NULL, col = 1, lty = 4)
  abline(median(sig$"Psi516"), 0, NULL, NULL, col = 4, lty = 4)
  legend("topright", legend = c("Component #1", "Component #2"), fill = c(1, 4), bty = "n", border = NA)}
##visualize down-sample results
