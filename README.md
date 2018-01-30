# U6L1
## Description
This script takes the aligned bam file to a repeat-masked reference and the annotation files together with the expression level of each gene to output the genes fused with U6 which could not be aligned to the reference genome inlcuding U6-L1 fusion reads, together with the number of supporting reads and read names

## Required Resources
Python:		https://www.python.org <br />
Samtools: http://samtools.sourceforge.net/ <br />
FLASH: https://ccb.jhu.edu/software/FLASH/ <br />
bwa: https://github.com/lh3/bwa <br />
picard tools: https://broadinstitute.github.io/picard/

## Method
The process of how the fusion reads are extracted:
![Alt text](https://github.com/mills-lab/U6L1/blob/master/pipeline.png)

