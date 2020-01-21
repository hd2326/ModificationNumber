#Code for data pre-processing

/home/ubuntu/nanopolish/nanopolish index -d /home/ubuntu/mrna_na12878/ucsc_run1/fast5/ UCSC_Run1_20170907_DirectRNA.pass.dedup.fastq &> index.log
#indexing .fastq

minimap2 -ax map-ont /home/ubuntu/refGenome/Hs/fa/Homo_sapiens.GRCh38.dna.fa UCSC_Run1_20170907_DirectRNA.pass.dedup.fastq -o UCSC_Run1_20170907_DirectRNA.pass.dedup.sam -t 32 &> minimap2.log
#alignment .fastq

samtools view -b UCSC_Run1_20170907_DirectRNA.pass.dedup.sam > UCSC_Run1_20170907_DirectRNA.pass.dedup.bam
#.sam to .bam

samtools sort -o  UCSC_Run1_20170907_DirectRNA.pass.dedup.sorted.bam UCSC_Run1_20170907_DirectRNA.pass.dedup.bam
#sorting .bam

samtools index UCSC_Run1_20170907_DirectRNA.pass.dedup.sorted.bam
#indexing .bam

/home/ubuntu/nanopolish/nanopolish eventalign -o /home/ubuntu/mrna_na12878/ucsc_run1/event/ -r UCSC_Run1_20170907_DirectRNA.pass.dedup.fastq -b UCSC_Run1_20170907_DirectRNA.pass.dedup.sorted.bam -g /home/ubuntu/refGenome/Hs/fa/Homo_sapiens.GRCh38.dna.fa -t 32 --cigar_output --scale-events &> eventalign.log
#event table generation
