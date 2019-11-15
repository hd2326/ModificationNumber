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
fit <- structure(lapply(5:9, function(i, sig) emtest.norm(sig, m0 = i), sig=unlist(sig)), names=3:7)
save(sig, fit, file = "primer_extension_fit_27.rda")
#position 27, GCCTGA, original paper showcase
sig <- structure(lapply(base, function(b, table) table$event_level_mean[table$position == 34 & table$contig == b], table=table), names=base)
sig <- lapply(sig, function(x, median, mad) x[x > median-2*mad & x < median+2*mad], median=median(unlist(sig)), mad=mad(unlist(sig)))
sig <- sig[which(unlist(lapply(sig, function(x) length(x) > 500)))]
fit <- structure(lapply(3:7, function(i, sig) emtest.norm(sig, m0 = i), sig=unlist(sig)), names=3:7)
save(sig, fit, file = "primer_extension_fit_34.rda")
#position 34, CATCGC
#mixture model

load("primer_extension_table.rda")
read <- unlist(lapply(table, function(x) sum(25:36 %in% x$position) == 12))
contig <- unlist(lapply(table[read], function(x) x$contig[1]))
sig_25_30 <- lapply(table[read], function(x) structure(unlist(lapply(25:30, function(p, x) mean(x$event_level_mean[x$position == p]), x=x)), names=25:30))
sig_25_30 <- do.call(rbind, sig_25_30)
sig_31_36 <- lapply(table[read], function(x) structure(unlist(lapply(31:36, function(p, x) mean(x$event_level_mean[x$position == p]), x=x)), names=31:36))
sig_31_36 <- do.call(rbind, sig_31_36)
sig_25_36 <- cbind(sig_25_30, sig_31_36)
#clustering reads
