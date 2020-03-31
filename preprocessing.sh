nanopolish=/home/ubuntu/nanopolish/nanopolish
minimap2=/home/ubuntu/minimap2/minimap2
samtools=/usr/bin/samtools
#path to nanopolish, minimap2 and samtools

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
