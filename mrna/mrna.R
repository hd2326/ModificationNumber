files <- list.files(path = "./event/", pattern = ".tsv", full.names = T)
table <- list()
pb <- txtProgressBar(min = 0, max = length(files), style = 3)
for (f in files){
  setTxtProgressBar(pb, grep(f, files))
  t <- read.delim(f, header = T, stringsAsFactors = F)
  n <- sapply(strsplit(unlist(strsplit(f, split = ".tsv")), split = "./event//"), function(x) x[2])
  if(nrow(t) > 0 & sum(is.na(t)) == 0 & sum(t$contig != "1") == 0){
    table[[n]] <- t[, c("position", "event_index", "event_level_mean", "model_kmer")]}}
#event table

system("samtools view ./fastq/UCSC_Run1_20170907_DirectRNA.pass.dedup.sorted.bam | awk '{ print $1, $3, $5}' > ./fastq/mrna_sam.txt")
sam <- read.delim("./fastq/mrna_sam.txt", header = F, stringsAsFactors = F, sep = " ", quote = "")
colnames(sam) <- c("id", "contig", "mapq")
sam <- sam[sam$id %in% names(table), ]
#sam file

fastq <- read.delim("./fastq/UCSC_Run1_20170907_DirectRNA.pass.dedup.fastq", header = F, stringsAsFactors = F, sep = " ", quote = "")
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
  s})), names=names(table))#strand distribution
mapq <- structure(sam$mapq, names=sam$id)#mapq score
map2 <- structure(duplicated(sam$id, fromLast = T) + duplicated(sam$id, fromLast = F) == 0, names=sam$id)#no sec/sup mapping
par(mfrow = c(2, 3), mar = c(5, 5, 5, 2))
barplot(table(contig), ylab = "#Reads", cex.names = 2, cex.axis = 2, cex.lab = 2, names.arg = names(table(contig)), main = "Contig", cex.main = 2)
barplot(table(strand), ylab = "#Reads", cex.names = 2, cex.axis = 2, cex.lab = 2, main = "Strand", cex.main = 2)
hist(mapq, xlab = "", ylab = "#Reads", cex.axis = 2, cex.lab = 2, main = "MAPQ", cex.main = 2)
barplot(table(map2), ylab = "#Reads", cex.names = 2, cex.axis = 2, cex.lab = 2, main = "No Sec/Sup Mapping", cex.main = 2)
plot(density(unlist(lapply(q, function(x) mean(x-33)))), xlab = "", cex.axis = 2, cex.lab = 2, main = "Average Q-Score", cex.main = 2)
plot(density(unlist(lapply(seq, function(x) length(x)))), xlab = "", cex.axis = 2, cex.lab = 2, main = "Read Length", cex.main = 2)
#basic qc

ind <- Reduce(intersect, list(names(strand[strand != "NA"]),
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
  table[[i]]$strand <- rep(strand[i], nrow(table[[i]]))
  table[[i]]$contig <- rep(contig[i], nrow(table[[i]]))}
#filter/reform/save event table

con <- file("/home/ubuntu/refGenome/Hs/fa/Homo_sapiens.GRCh38.dna.fa", "r", blocking = F)
ref <- readLines(con)
close(con)
chr <- grep(">", ref)
ref <- ref[(chr[1]+1):(chr[2]-1)]
ref <- paste(ref, collapse = "")
#reference sequence

kmer <- function(k){
  return(apply(expand.grid(lapply(1:k, function(i, base) base, base=c("A", "C", "G", "T"))), 1, function(x) paste(rev(x), collapse = "")))}
#kmer generator

table <- lapply(table, function(x, ref){
  if (sum(table$event_index != sort(table$event_index, decreasing=F)) == 0){
    ref <- substr(ref, table$position[1]-1, table$position[nrow(table)]+7)
    table <- table[, c("position", "event_level_mean")]
    table$seq <- unlist(lapply(table$position-table$position[1]+1, function(i, ref) substr(ref, i, i+8), ref=ref))
    table}
}, ref=ref)
save(table, file = "./mrna_table_I.rda")
#forward strand events

table <- do.call(rbind, table[unlist(lapply(table, function(x) !is.null(x)))])
sig_k <- list(sig_7_2_0=structure(lapply(kmer(7), function(k, table) table[substr(table$seq, 1, 7) == k, c("position", "event_level_mean")], table=table), names=kmer(7)),
              sig_7_1_1=structure(lapply(kmer(7), function(k, table) table[substr(table$seq, 2, 8) == k, c("position", "event_level_mean")], table=table), names=kmer(7)),
              sig_7_0_2=structure(lapply(kmer(7), function(k, table) table[substr(table$seq, 3, 9) == k, c("position", "event_level_mean")], table=table), names=kmer(7)),
              sig_6_1_0=structure(lapply(kmer(6), function(k, table) table[substr(table$seq, 2, 7) == k, c("position", "event_level_mean")], table=table), names=kmer(6)),
              sig_6_0_1=structure(lapply(kmer(6), function(k, table) table[substr(table$seq, 3, 8) == k, c("position", "event_level_mean")], table=table), names=kmer(6)),
              sig_5_0_0=structure(lapply(kmer(5), function(k, table) table[substr(table$seq, 3, 7) == k, c("position", "event_level_mean")], table=table), names=kmer(5)),
              sig_4_1_0=structure(lapply(kmer(4), function(k, table) table[substr(table$seq, 4, 7) == k, c("position", "event_level_mean")], table=table), names=kmer(4)),
              sig_4_0_1=structure(lapply(kmer(4), function(k, table) table[substr(table$seq, 3, 6) == k, c("position", "event_level_mean")], table=table), names=kmer(4)),
              sig_3_2_0=structure(lapply(kmer(3), function(k, table) table[substr(table$seq, 5, 7) == k, c("position", "event_level_mean")], table=table), names=kmer(3)),
              sig_3_1_1=structure(lapply(kmer(3), function(k, table) table[substr(table$seq, 4, 6) == k, c("position", "event_level_mean")], table=table), names=kmer(3)),
              sig_3_0_2=structure(lapply(kmer(3), function(k, table) table[substr(table$seq, 3, 5) == k, c("position", "event_level_mean")], table=table), names=kmer(3)))
#kmer positional signal mad

k <- sapply(strsplit(names(sig_k), split="_"), function(x) as.integer(x[2]))
d_k <- lapply(sig_k, function(sig) unlist(lapply(sig, function(s) mad(s$event_level_mean))))
par(mfrow = c(1, 2), mar = c(5, 5, 2, 2))
plot(NA, NA, xlim = c(0, 6), ylim = c(0, 2), xlab = "log2(mad)", ylab = "p")
for (i in 1:length(d_k)) lines(density(log2(d_k[[i]]), na.rm = T), col = k[i]+1)
legend("topright", legend = paste("k=", 7:3, sep = ""), fill = 8:4, bty = "n", border = NA)
plot(NA, NA, xlim = c(0, 6), ylim = c(0, 2), xlab = "log2(mad)", ylab = "p")
for (i in 1:length(d_k)) lines(density(log2(d_k[[i]]), na.rm = T), col = rainbow(12)[i])
legend("topright", legend = names(d_k), fill = rainbow(12), bty = "n", border = NA)
                                                 
d_k <- lapply(sig_k, function(sig) unlist(lapply(sig, function(s) mad(s$event_level_mean))))
d <- matrix(1, length(d_k), length(d_k), dimnames = list(names(d_k), names(d_k)))
for (i in 1:length(d_k)) for (j in 1:length(d_k)) d[i, j] <- ks.test(d_k[[i]], d_k[[j]])$statistic
image(d, axes = F, main = "NA12878 mRNA")
axis(side = 1, at = seq(0, 1, length.out = nrow(d)), labels = rownames(d), tick = F, las = 2, cex.axis = 0.5)
axis(side = 2, at = seq(0, 1, length.out = ncol(d)), labels = colnames(d), tick = F, las = 2, cex.axis = 0.5)
for (i in 1:nrow(d)) for (j in 1:ncol(d)) text(seq(0, 1, length.out = nrow(d))[i],
                                               seq(0, 1, length.out = ncol(d))[j], labels = signif(d[i, j], 2), cex = 0.5)
                                                 
#determining optimal number of k

pb <- txtProgressBar(min = 0, max = 5*4^4, style = 3)
sig_p <- structure(lapply(1:5, function(i, table, x, pb){
  structure(lapply(x, function(xx, table, i, pb){
    setTxtProgressBar(pb, 4^4*(i-1)+grep(xx, x))
    list("A"=table$event_level_mean[substr(table$seq, 3, 7) == paste(substr(xx, 0, (i-1)), "A", substr(xx, i, 4), sep = "")],
         "C"=table$event_level_mean[substr(table$seq, 3, 7) == paste(substr(xx, 0, (i-1)), "C", substr(xx, i, 4), sep = "")],
         "G"=table$event_level_mean[substr(table$seq, 3, 7) == paste(substr(xx, 0, (i-1)), "G", substr(xx, i, 4), sep = "")],
         "T"=table$event_level_mean[substr(table$seq, 3, 7) == paste(substr(xx, 0, (i-1)), "T", substr(xx, i, 4), sep = "")])
  }, table=table, i=i, pb=pb), names=x)
}, table=table, x=kmer(4), pb=pb), names=1:5)
close(pb)
d_p <- lapply(sig_p, function(s){
  structure(lapply(s, function(ss){
    dist <- matrix(0, 4, 4, dimnames = list(c("A", "C", "G", "T"), c("A", "C", "G", "T")))
    for (x in rownames(dist)) for (y in colnames(dist)) dist[x, y] <- median(ss[[x]]) - median(ss[[y]])
    dist}), names=names(s))})
names(d_p) <- sub("sig", "d", names(sig_p))
#5mer nucleotide positional median difference

par(mfrow = c(2, 3), mar = c(4, 4, 2, 2))
plot(NA, NA, xlab = "kmer median difference", ylab = "p", xlim = c(0, 50), ylim = c(0, 0.6), main = "")
for (i in 1:length(d_p)) lines(density(unlist(lapply(d_p[[i]], function(x) abs(x[upper.tri(x, diag = F)])))), col = i)
legend("topright", legend = paste("Pos:", 1:5, sep = ""), fill = 1:5, bty = "n", border = NA)
dd <- lapply(d_p, function(x) do.call(rbind, lapply(x, function(xx) colMeans(xx))))
for (i in 1:length(dd)){
  plot(NA, NA, xlab = "mean base median difference", ylab = "p", xlim = c(-25, 25), ylim = c(0, c(1, 0.2, 0.2, 0.4, 0.4, 1)[i]), main = "")
  for (j in 1:4) lines(density(dd[[i]][, j]), col = i, lty = j)}
legend("topright", legend = c("A", "C", "G", "T"), lty = 1:4, bty = "n", border = NA)
#visualization
