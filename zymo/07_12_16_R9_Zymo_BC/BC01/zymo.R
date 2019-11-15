files <- list.files(path = "./event/BC01/", pattern = ".tsv", full.names = T)
table <- list()
for (f in files){
  t <- read.delim(f, header = T, stringsAsFactors = F)
  n <- sapply(strsplit(unlist(strsplit(f, split = ".tsv")), split = "./event/BC01//"), function(x) x[2])
  if(nrow(t) > 0 & sum(is.na(t)) == 0) table[[n]] <- t[, c("position", "event_index", "event_level_mean", "model_kmer", "cigar_string", "cigar_basecalled_kmer", "basecall_index")]}
#event table

system("samtools view ./fastq/BC01.sorted.bam | awk '{ print $1, $3, $5}' > ./fastq/BC01_sam.txt")
sam <- read.delim("./fastq/BC01_sam.txt", header = F, stringsAsFactors = F, sep = " ", quote = "")
colnames(sam) <- c("id", "contig", "mapq")
sam <- sam[sam$id %in% names(table), ]
#sam file

fastq <- read.delim("./fastq/BC01.fastq", header = F, stringsAsFactors = F, sep = " ", quote = "")
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
map2 <- structure(duplicated(sam$id, fromLast = T) + duplicated(sam$id, fromLast = F) == 0, names=sam$id)#no secondary mapping
par(mfrow = c(2, 3), mar = c(5, 5, 5, 2))
barplot(table(contig), ylab = "#Reads", cex.names = 2, cex.axis = 2, cex.lab = 2,  main = "Contig", cex.main = 2)
barplot(table(strand), ylab = "#Reads", cex.names = 2, cex.axis = 2, cex.lab = 2, main = "Strand", cex.main = 2)
hist(mapq, xlab = "", ylab = "#Reads", cex.axis = 2, cex.lab = 2, main = "MAPQ", cex.main = 2)
barplot(table(map2), ylab = "#Reads", cex.names = 2, cex.axis = 2, cex.lab = 2, main = "No Sec/Sup Mapping", cex.main = 2)
plot(density(unlist(lapply(q, function(x) mean(x-33)))), xlab = "", cex.axis = 2, cex.lab = 2, main = "Average Q-Score", cex.main = 2)
plot(density(unlist(lapply(seq, function(x) length(x)))), xlab = "", cex.axis = 2, cex.lab = 2, main = "Read Length", cex.main = 2)
dev.off()
#basic qc

ind <- Reduce(intersect, list(names(contig[contig == "ZYMO"]),
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
save(table, file = "BC01_table.rda")
#filter/reform/save event table

load("BC01_table.rda")
table <- do.call(rbind, table)
table <- table[table$model_kmer != "NNNNNN", ]

kmer <- unique(table$model_kmer)
position <- structure(lapply(kmer, function(k, table) table$position[table$model_kmer == k], table=table), names = kmer)
sig <- structure(lapply(kmer, function(k, table) table$event_level_mean[table$model_kmer == k], table=table), names=kmer)
strand <- structure(lapply(kmer, function(k, table) table$strand[table$model_kmer == k], table=table), names=kmer)
par(mfrow = c(10, 10), mar = c(4, 4, 4, 1))
for (k in kmer){
  p <- table(position[[k]])
  plot(NA, NA, xlab = "pA", ylab = "Density", xlim = range(unlist(sig[[k]])), ylim = c(0, 0.5), main = k)
  for (i in 1:length(p)){
    f <- as.character(position[[k]]) == names(p)[i] & strand[[k]] == "F"
    r <- as.character(position[[k]]) == names(p)[i] & strand[[k]] == "R"
    if(sum(f) > 1000) lines(density(sig[[k]][f]), lty = 1, col = i)
    if(sum(r) > 1000) lines(density(sig[[k]][r]), lty = 4, col = i)}
  legend("topleft", legend = c("F", "R"), lty = c(1, 4), bty = "n", border = NA)
  legend("topright", legend = names(p), fill = 1:8, bty = "n", border = NA)}
#kmer distribution

library(e1071)
ref <- as.character(read.delim("../zymo.fa", stringsAsFactors = F))
position <- unique(table$position)
strand <- structure(lapply(position, function(p, table) table$strand[table$position == p], table=table), names=position)
sig <- structure(lapply(position, function(p, table) list(f=table$event_level_mean[table$position == p & table$strand == "F"],
                                                          r=table$event_level_mean[table$position == p & table$strand == "R"]), table=table), names=position)
mean <- structure(lapply(sig, function(x) list(f=mean(x$f), r=mean(x$r))), names=names(sig))
median <- structure(lapply(sig, function(x) list(f=median(x$f), r=median(x$r))), names=names(sig))
sd <- structure(lapply(sig, function(x) list(f=sd(x$f), r=sd(x$r))), names=names(sig))
mad <- structure(lapply(sig, function(x) list(f=mad(x$f), r=mad(x$r))), names=names(sig))
z <- structure(lapply(names(sig), function(x, sig, median, mad) list(f=(sig[[x]]$f-median[[x]]$f)/mad[[x]]$f,
                                                                     r=(sig[[x]]$r-median[[x]]$r)/mad[[x]]$r), sig=sig, median=median, mad=mad), names=names(sig))
kurtosis <- structure(lapply(sig, function(x) list(f=kurtosis(x$f[x$f>median(x$f)-3*mad(x$f) & x$f<median(x$f)+2*mad(x$f)]),
                                                   r=kurtosis(x$r[x$r>median(x$r)-3*mad(x$r) & x$r<median(x$r)+2*mad(x$r)]))), names = names(sig))
skewness <- structure(lapply(sig, function(x) list(f=skewness(x$f[x$f>median(x$f)-3*mad(x$f) & x$f<median(x$f)+2*mad(x$f)]),
                                                   r=skewness(x$r[x$r>median(x$r)-3*mad(x$r) & x$r<median(x$r)+2*mad(x$r)]))), names = names(sig))
m <- structure(lapply(position, function(p, table) list(f=unlist(lapply(gregexpr("M", table$cigar_string[table$position == p & table$strand == "F"]), function(x) length(as.integer(x)))),
                                                        r=unlist(lapply(gregexpr("M", table$cigar_string[table$position == p & table$strand == "R"]), function(x) length(as.integer(x))))), table=table), names = position)
q <- structure(lapply(position, function(p, table) list(f=table$q[table$position == p & table$strand == "F"],
                                                        r=table$q[table$position == p & table$strand == "R"]), table=table), names = position)

library(Biostrings)
par(mfrow = c(10, 10), mar = c(4, 4, 4, 1))
for (p in as.character(sort(position, decreasing = F))){
  plot(NA, NA, xlab = "pA", ylab = "Density", xlim = range(unlist(sig[[p]])), ylim = c(0, 0.5), main = paste("Position", p))
  kmer <- substr(ref, as.integer(p)+1, as.integer(p)+6)
  if (length(sig[[p]]$f) > 1000){
    lines(density(sig[[p]]$f), col = 4, lty = 1)
    lines(seq(0, 200, length.out = 1000), dnorm(seq(0, 200, length.out = 1000), mean = mean[[p]]$f, sd = sd[[p]]$f), col = 4, lty = 2)
    lines(seq(0, 200, length.out = 1000), dnorm(seq(0, 200, length.out = 1000), mean = median[[p]]$f, sd = mad[[p]]$f), col = 4, lty = 3)}
  if (length(sig[[p]]$r) > 1000){
    lines(density(sig[[p]]$r), col = 6, lty = 1)
    lines(seq(0, 200, length.out = 1000), dnorm(seq(0, 200, length.out = 1000), mean = mean[[p]]$r, sd = sd[[p]]$r), col = 6, lty = 2)
    lines(seq(0, 200, length.out = 1000), dnorm(seq(0, 200, length.out = 1000), mean = median[[p]]$r, sd = mad[[p]]$r), col = 6, lty = 3)}
  legend("bottomright", legend = c("Forward", "Reverse"), fill = c(4, 6), bty = "n", border = NA)
  legend("bottomleft", legend = paste(c(kmer, as.character(reverseComplement(DNAString(kmer)))), c(length(sig[[p]]$f), length(sig[[p]]$r)), sep = ":"), text.col = c(4, 6), bty = "n", border = NA)
  legend("topleft", legend = c(paste("K:", signif(kurtosis[[p]]$f, 2), " S:", signif(skewness[[p]]$f, 2), sep = ""),
                               paste("K:", signif(kurtosis[[p]]$r, 2), " S:", signif(skewness[[p]]$r, 2), sep = "")), text.col = c(4, 6), bty = "n", border = NA)
  legend("topright", legend = c("empirical", "mean/sd", "median/mad"), lty = 1:3, bty = "n", border = NA)}
#position distribution

library(vioplot)
par(mfcol = c(2, 5), mar = c(5, 5, 2, 2))
smoothScatter(unlist(lapply(mean, function(x) x$f)), unlist(lapply(median, function(x) x$f)), nrpoints = 0,
              xlab = "mean", ylab = "median", cex.axis = 2, cex.lab = 2, colramp = colorRampPalette(c("white", 4)), main = "")
abline(0, 1, NULL, NULL, col = 2)
smoothScatter(unlist(lapply(mean, function(x) x$r)), unlist(lapply(median, function(x) x$r)), nrpoints = 0,
              xlab = "mean", ylab = "median", cex.axis = 2, cex.lab = 2, colramp = colorRampPalette(c("white", 6)), main = "")
abline(0, 1, NULL, NULL, col = 2)
legend("topleft", legend = c("F", "R"), fill = c(4, 6), cex = 2, bty = "n", border = NA)

smoothScatter(unlist(lapply(sd, function(x) x$f)), unlist(lapply(mad, function(x) x$f)), nrpoints = 0,
              xlab = "sd", ylab = "mad", cex.axis = 2, cex.lab = 2, colramp = colorRampPalette(c("white", 4)), main = "")
abline(0, 1, NULL, NULL, col = 2)
smoothScatter(unlist(lapply(sd, function(x) x$r)), unlist(lapply(mad, function(x) x$r)), nrpoints = 0,
              xlab = "sd", ylab = "mad", cex.axis = 2, cex.lab = 2, colramp = colorRampPalette(c("white", 6)), main = "")
abline(0, 1, NULL, NULL, col = 2)

smoothScatter(unlist(lapply(skewness, function(x) x$f)), unlist(lapply(kurtosis, function(x) x$f)), nrpoints = 0,
              xlab = "skewness", ylab = "kurtosis", cex.axis = 2, cex.lab = 2, colramp = colorRampPalette(c("white", 4)), main = "")
smoothScatter(unlist(lapply(skewness, function(x) x$r)), unlist(lapply(kurtosis, function(x) x$r)), nrpoints = 0,
              xlab = "skewness", ylab = "kurtosis", cex.axis = 2, cex.lab = 2, colramp = colorRampPalette(c("white", 6)), main = "")

vioplot(lapply(sort(unique(unlist(lapply(q, function(x) x$f)))), function(x, q, z) z[q == x], q=unlist(lapply(q, function(x) x$f)), z=unlist(lapply(z, function(x) x$f))),
        xlab = "q-score", ylab = "z-score", col = 4, names = sort(unique(unlist(lapply(q, function(x) x$f-33)))), cex.axis = 2, cex.lab = 2, cex.names = 2, main = "")
vioplot(lapply(sort(unique(unlist(lapply(q, function(x) x$r)))), function(x, q, z) z[q == x], q=unlist(lapply(q, function(x) x$r)), z=unlist(lapply(z, function(x) x$r))),
        xlab = "q-score", ylab = "z-score", col = 6, names = sort(unique(unlist(lapply(q, function(x) x$r-33)))), cex.axis = 2, cex.lab = 2, cex.names = 2, main = "")
vioplot(lapply(1:6, function(x, m, z) z[m == x], m=unlist(lapply(m, function(x) x$f)), z=unlist(lapply(z, function(x) x$f))),
        xlab = "#M", ylab = "z-score", col = 4, cex.axis = 2, cex.lab = 2, cex.names = 2, main = "")
vioplot(lapply(1:6, function(x, m, z) z[m == x], m=unlist(lapply(m, function(x) x$r)), z=unlist(lapply(z, function(x) x$r))),
        xlab = "#M", ylab = "z-score", col = 6, cex.axis = 2, cex.lab = 2, cex.names = 2, main = "")
dev.off()
#signal statistics

load("BC01_table.rda")
table <- do.call(rbind, table)
table <- table[table$model_kmer != "NNNNNN", c("position", "event_level_mean", "strand")]
position <- unique(table$position)
strand <- structure(lapply(position, function(p, table) table$strand[table$position == p], table=table), names=position)
sig <- structure(lapply(position, function(p, table) list(f=table$event_level_mean[table$position == p & table$strand == "F"],
                                                          r=table$event_level_mean[table$position == p & table$strand == "R"]), table=table), names=position)
median_bc01 <- structure(lapply(sig, function(x) list(f=median(x$f), r=median(x$r))), names=names(sig))
mad_bc01 <- structure(lapply(sig, function(x) list(f=mad(x$f), r=mad(x$r))), names=names(sig))
par(mfcol = c(1, 2), mar = c(5, 5, 5, 2))
smoothScatter(unlist(lapply(median_bc01, function(x) x$f)), unlist(lapply(mad_bc01, function(x) x$f)), nrpoints = 0,
              xlab = "median", ylab = "mad", cex.axis = 2, cex.lab = 2, colramp = colorRampPalette(c("white", 4)), main = "C, F", cex.main = 2)
smoothScatter(unlist(lapply(median_bc01, function(x) x$r)), unlist(lapply(mad_bc01, function(x) x$r)), nrpoints = 0,
              xlab = "median", ylab = "mad", cex.axis = 2, cex.lab = 2, colramp = colorRampPalette(c("white", 6)), main = "C, R", cex.main = 2)
#median-mad distribution
