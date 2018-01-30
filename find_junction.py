#This simple script takes a fastq file or a sam file to output the reads that contains certain junction motif that you would like to search with exact match
#Output of this script is a pure text file with the reads (read content only for fastq files, and line of read for sam files) that contains the junction as an input from the commmand
#Use "wc -l OUTPUT.TXT", you could see the number of reads containing certain motif

import os
import sys
motif = sys.argv[1] #motif sequence, ex. "ATCATGAC"
input = sys.argv[2] #input file, fastq or sam
output = sys.argv[3] #output text file
os.system(f"grep {motif} {input} > {output}")
