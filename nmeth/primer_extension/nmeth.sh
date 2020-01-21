#Code for data pre-processing

/home/ubuntu/nanopolish/nanopolish index -d /brdu_controls_from_muller/primer_extension/fast5/ 2017_06_05_CAM_ONT_BrdU_PCRs_primer_ext_basecalled_guppy3.1.5.fastq &> index.log
#indexing .fastq

minimap2 -ax map-ont /brdu_controls_from_muller/primer_extension/ref.fa /brdu_controls_from_muller/primer_extension/fastq/2017_06_05_CAM_ONT_BrdU_PCRs_primer_ext_basecalled_guppy3.1.5.fastq -o primer_extension.sam &> minimap2.log
#alignment .fastq

samtools view -b primer_extension.sam > primer_extension.bam
#.sam to .bam

samtools sort -o primer_extension.sorted.bam primer_extension.bam
#sorting .bam

samtools index primer_extension.sorted.bam
#indexing .bam

/home/ubuntu/nanopolish/nanopolish eventalign -o /brdu_controls_from_muller/primer_extension/event/ -r 2017_06_05_CAM_ONT_BrdU_PCRs_primer_ext_basecalled_guppy3.1.5.fastq -b primer_extension.sorted.bam -g /brdu_controls_from_muller/primer_extension/ref.fa -t 18 --cigar_output --scale-events &> eventalign.log
#event table generation
