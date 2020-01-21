#Code for data pre-processing

/home/ubuntu/nanopolish/nanopolish extract -v -r -q -t template -o 07_12_16_R9_Zymo_BC.fastq /home/ubuntu/zymo/07_12_16_R9_Zymo_BC/ &> extract.log
#extract .fastq from .fast5

porechop -i /home/ubuntu/zymo/07_12_16_R9_Zymo_BC/fastq/07_12_16_R9_Zymo_BC.fastq -b /home/ubuntu/zymo/07_12_16_R9_Zymo_BC/fastq/ --untrimmed -t 12 &> porechop.log
#demultiplexing .fastq

for f in BC*.fastq; do /home/ubuntu/nanopolish/nanopolish index -d /home/ubuntu/zymo/07_12_16_R9_Zymo_BC/fast5/ $f &> ${f%.fastq}_index.log ; done
#indexing .fastq

for f in BC*.fastq; do minimap2 -ax map-ont /home/ubuntu/zymo/zymo.fa $f -o ${f%.fastq}.sam &> ${f%.fastq}_minimap2.log; done
#alignment .fastq

for f in BC*.sam; do samtools view -b $f > ${f%.sam}.bam; done
#.sam to .bam

for f in BC*.bam; do samtools sort -o  ${f%.bam}.sorted.bam $f; done
#sorting .bam

for f in BC*.sorted.bam; do samtools index $f; done
#indexing .bam

for f in BC*.sorted.bam; do /home/ubuntu/nanopolish/nanopolish eventalign -o /home/ubuntu/zymo/07_12_16_R9_Zymo_BC/event/${f%.sorted.bam}/ -r ${f%.sorted.bam}.fastq -b $f -g /home/ubuntu/zymo/zymo.fa -t 18 --cigar_output --scale-events &> ${f%.sorted.bam}_eventalign.log;  done
#event table generation
