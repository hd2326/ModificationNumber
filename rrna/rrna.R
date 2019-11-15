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
fit <- structure(lapply(1:5, function(i, s) emtest.norm(s, m0 = i), s=unlist(s)), names=1:5)
save(s, fit, file = "rrna_fit_515.rda")
#position 515, Tgcca
s <- sig$"523"
median <- structure(lapply(s, function(x) median(x)), names=names(s))
mad <- structure(lapply(s, function(x) mad(x)), names=names(s))
s <- structure(lapply(names(s), function(x, s, median, mad) s[[x]][s[[x]] > median[[x]]-2*mad[[x]] & s[[x]] < median[[x]]+2*mad[[x]]], s=s, median=median, mad=mad), names=names(s))
fit <- structure(lapply(1:5, function(i, s) emtest.norm(s, m0 = i), s=unlist(s)), names=1:5)
save(s, fit, file = "rrna_fit_523.rda")
#position 523, gccGc

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
