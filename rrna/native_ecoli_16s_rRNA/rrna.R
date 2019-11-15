files <- list.files(path = "./event/native_ecoli_16s_rRNA/", pattern = ".tsv", full.names = T)
table <- list()
for (f in files){
  t <- read.delim(f, header = F, stringsAsFactors = F)
  n <- sapply(strsplit(unlist(strsplit(f, split = ".sm.forward.tsv")), split = "./event/native_ecoli_16s_rRNA//"), function(x) x[2])
  if(nrow(t) > 0 & sum(is.na(t)) == 0){
    t <- t[, c(2, 6, 13, 14, 16)]
    colnames(t) <- c("position", "event_index", "probability", "event_level_mean", "model_kmer")
    table[[n]] <- t[t$probability >= 0.9, ]}}
#event table

system("samtools view ./fastq/native_ecoli_16s_rRNA.sorted.bam | awk '{ print $1, $3, $5}' > ./fastq/native_ecoli_16s_rRNA_sam.txt")
sam <- read.delim("./fastq/native_ecoli_16s_rRNA_sam.txt", header = F, stringsAsFactors = F, sep = " ", quote = "")
colnames(sam) <- c("id", "contig", "mapq")
sam$id <- unlist(strsplit(sam$id, split = "_Basecall_1D_template"))
sam <- sam[sam$id %in% names(table), ]
#sam file

fastq <- read.delim("./fastq/native_ecoli_16s_rRNA.fastq", header = F, stringsAsFactors = F, sep = " ", quote = "")
fastq <- data.frame(id=unlist(lapply(strsplit(fastq$V1[seq(1, nrow(fastq), by = 4)], split = "@"), function(x) x[2])),
                    seq=fastq$V1[seq(2, nrow(fastq), by = 4)],
                    q=fastq$V1[seq(4, nrow(fastq), by = 4)], stringsAsFactors = F)
fastq$id <- unlist(strsplit(fastq$id, split = "_Basecall_1D_template"))
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
barplot(table(contig), ylab = "#Reads", cex.names = 2, cex.axis = 2, cex.lab = 2, main = "Contig", cex.main = 2)
barplot(table(strand), ylab = "#Reads", cex.names = 2, cex.axis = 2, cex.lab = 2, main = "Strand", cex.main = 2)
hist(mapq, xlab = "", ylab = "#Reads", cex.axis = 2, cex.lab = 2, main = "MAPQ", cex.main = 2)
barplot(table(map2), ylab = "#Reads", cex.names = 2, cex.axis = 2, cex.lab = 2, main = "No Sec/Sup Mapping", cex.main = 2)
plot(density(unlist(lapply(q, function(x) mean(x-33)))), xlab = "", cex.axis = 2, cex.lab = 2, main = "Average Q-Score", cex.main = 2)
plot(density(unlist(lapply(seq, function(x) length(x)))), xlab = "", cex.axis = 2, cex.lab = 2, main = "Read Length", cex.main = 2)
#basic qc

ind <- Reduce(intersect, list(names(contig[contig == "ecoli_MRE600"]),
                              names(strand[strand == "F"]),
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
save(table, file = "native_ecoli_16s_rRNA_table.rda")
#filter/reform/save event table
