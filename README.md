Scripts in this repository can be used to reproduce all results in manuscript "Gaussian Mixture Model-Based Unsupervised Nucleotide Modification Number Detection Using Nanopore Sequencing Readouts".

Scripts were organized in folders based on the datasets they were applied to. For each folder, the .sh script was used for data pre-processing, and the .R script was used for the actual analysis.

In total, including
1. The Zymo native synthesized oligo nanopore sequencing dataset (./zymo/) [1].
2. Thymidine analogs-containing primer extension dataset (./nmeth/primer_extension/) [2]
3. Native yeast genomic DNA nanopore sequencing datasets (./nmeth/yeast/) [2]
4. NA12878 cell line mRNA dataset (./mrna/) [3].
5. E.coli 16S rRNA nanopore sequencing dataset [4], containing the following 3 sub-datasets. Please note the event tables were provided by the authors of the original study therefore data pre-processing is not applicable here.
    5a. Native strain (./rrna/).
    5b. Pseudouridine-deficient (Psi516) strain (./rrna/).
    5c. m7G-deficient (m7G) strain (./rrna/).

References:
1. Rand, A. C., Jain, M., Eizenga, J. M., Musselman-Brown, A., Olsen, H. E., Akeson, M., & Paten, B. (2017). Mapping DNA methylation with high-throughput nanopore sequencing. Nature methods, 14(4), 411.
2. Mueller, C. A., Boemo, M. A., Spingardi, P., Kessler, B. M., Kriaucionis, S., Simpson, J. T., & Nieduszynski, C. A. (2019). Capturing the dynamics of genome replication on individual ultra-long nanopore sequence reads. Nature methods, 16(5), 429.
3. Workman, R. E., Tang, A. D., Tang, P. S., Jain, M., Tyson, J. R., Razaghi, R., ... & Sadowski, N. (2019). Nanopore native RNA sequencing of a human poly (A) transcriptome. Nature methods, 16(12), 1297-1305.
4. Smith, A. M., Jain, M., Mulroney, L., Garalde, D. R., & Akeson, M. (2019). Reading canonical and modified nucleobases in 16S ribosomal RNA using nanopore native RNA sequencing. PloS one, 14(5), e0216709.
