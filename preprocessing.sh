#Code for data pre-processing

nanopolish=/home/ubuntu/nanopolish/nanopolish
porechop=/home/ubuntu/anaconda3/bin/porechop
minimap2=/home/ubuntu/minimap2/minimap2
samtools=/usr/bin/samtools
#path to nanopolish, porechop, minimap2 and samtools

$nanopolish extract -v -r -q -t template -o /path/to/original/fastq/file /path/to/fast5/file/folder/ &> extract.log
#extract .fastq from .fast5

$porechop -i /path/to/original/fastq/file -b /path/to/demultiplexed/fastq/file/folder/ --untrimmed -t 12 &> porechop.log
#demultiplexing .fastq

$nanopolish index -d /path/to/fast5/file/folder/ /path/to/un-indexed/fastq/file &> index.log
#indexing .fastq

$minimap2 -ax map-ont /path/to/reference/genome/file /path/to/indexed/fastq/file -o /path/to/sam/file -t 32 &> minimap2.log
#alignment .fastq

$samtools view -b /path/to/sam/file > /path/to/bam/file
#.sam to .bam

$samtools sort -o  /path/to/sorted/bam/file /path/to/bam/file
#sorting .bam

$samtools index /path/to/sorted/bam/file
#indexing .bam

$nanopolish eventalign -o /path/to/event/table/folder/ -r /path/to/indexed/fastq/file -b /path/to/sorted/bam/file -g /path/to/reference/genome/file -t 32 --cigar_output --scale-events &> eventalign.log
#event table generation
