#Code for data pre-processing

for f in BC*.fastq; do /home/ubuntu/nanopolish/nanopolish index -d /brdu_controls_from_muller/yeast/fast5/ $f &> ${f%.fastq}_index.log ; done
#indexing .fastq

for f in BC*.fastq; do minimap2 -ax map-ont /home/ubuntu/refGenome/Sc/fa/Saccharomyces_cerevisiae.R64-1-1.dna.fa $f -o ${f%.fastq}.sam -t 32 &> ${f%.fastq}_minimap2.log; done
#alignment .fastq

for f in BC*.sam; do samtools view -b $f > ${f%.sam}.bam; done
#.sam to .bam

for f in BC*.bam; do samtools sort -o  ${f%.bam}.sorted.bam $f; done
#sorting .bam

for f in BC*.sorted.bam; do samtools index $f; done
#indexing .bam

for f in BC*.sorted.bam; do /home/ubuntu/nanopolish/nanopolish eventalign -o /brdu_controls_from_muller/yeast/event/${f%.sorted.bam}/ -r ${f%.sorted.bam}.fastq -b $f -g /home/ubuntu/refGenome/Sc/fa/Saccharomyces_cerevisiae.R64-1-1.dna.fa -t 32 --cigar_output --scale-events &> ${f%.sorted.bam}_eventalign.log;  done
#event table generation
