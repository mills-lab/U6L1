
# This script provides the final result from the aligned bam file.
# The read.fastq.gz file was sorted by using this command "zcat read.fastq.gz |paste - - - -| sort -k1,1 -S 20G | tr '\t' '\n'| gzip > read.sort.fastq.gz "
import sys
import os
import re
import gzip
# Input files
sample = sys.argv[1] #sample name and the name of the working directory
sambam_dir = sys.argv[2] #directory of the original full sam and bam files
sam_name = sys.argv[3] #name for sam file, XXX.sam
bam_name = sys.argv[4] #name for sorted bam file, XXX.bam
fastq_file_1 = sys.argv[5] #sorted fastq file 1
fastq_file_2 = sys.argv[6] #sorted fastq file 2
reference = sys.argv[7] #indexed masked reference
non_masked_ref = = sys.argv[8] #indexed non-masked reference
working_dir = sys.argv[9] #working directory
gtf_input = sys.argv[10] #whole sorted annotation file
gtf_file = sys.argv[11] #gtf_file includes only gene annotation
flash_dir = sys.argv[12] #where flash installed 
fpkm_file = sys.argv[13] #expression level input from the result of cuffliks
picard_dir = sys.argv[14] #directory for picard tools
whole_sam_input = f"{sambam_dir}/{sample}/{sam_name}"
whole_bam_input = f"{sambam_dir}/{sample}/{bam_name}"
os.system(f"mkdir {working_dir}"
os.system(f"mkdir {working_dir}/flash_result")
os.system(f"mkdir {working_dir}/alignment")
flash_install_dir = f"{working_dir}/flash_result/"
alignment_dir = f"{working_dir}s/alignment/"
sam_search_fused=["RNU6-1|ucsc|107"]
gtf_search_fused_li=["RNU6","U6"]

# output files
os.system(f"mkdir {working_dir}/sub_fastq")
sub_fastq_1 = f"{working_dir}/sub_fastq/{sample}_U6_1.fastq"
sub_fastq_2 = f"{working_dir}/sub_fastq/{sample}_U6_2.fastq"
sub_fastq_LINE_1 = f"{working_dir}/sub_fastq/{sample}_U6_LINE_1.fastq"
sub_fastq_LINE_2 = f"{working_dir}/sub_fastq/{sample}_U6_LINE_2.fastq"
sub_fastq_Alu_1 = f"{working_dir}/sub_fastq/{sample}_U6_Alu_1.fastq"
sub_fastq_Alu_2 = f"{working_dir}/sub_fastq/{sample}_U6_Alu_2.fastq"
sub_fastq_SVA_1 = f"{working_dir}/sub_fastq/{sample}_U6_SVA_1.fastq"
sub_fastq_SVA_2 = f"{working_dir}/sub_fastq/{sample}_U6_SVA_2.fastq"
os.system(f"mkdir {working_dir}/results")
single_matched_out=f"{working_dir}/results/single_matched_{sample}.txt"
paired_matched_out=f"{working_dir}/results/paired_matched_{sample}.txt"
single_nonmatched_out=f"{working_dir}/results/single_nonmatched_{sample}.txt"
paired_nonmatched_out=f"{working_dir}/results/paired_nonmatched_{working_dir}.txt"
matched_out = f"{working_dir}/results/matched_{sample}.txt"
nonmatched_out = f"{working_dir}/results/nonmatched_{sample}.txt"
os.system(f"mkdir {working_dir}/read_names")
read_name_other=f"{working_dir}/read_names/read_name_other_{sample}"
read_name_LINE=f"{working_dir}/read_names/read_name_LINE_{sample}"
read_name_Alu=f"{working_dir}/read_names/read_name_Alu_{sample}"
read_name_SVA=f"{working_dir}/read_names/read_name_SVA_{sample}"
read_name_other_sort = f"{read_name_other}_sort"
read_name_LINE_sort = f"{read_name_LINE}_sort"
read_name_Alu_sort = f"{read_name_Alu}_sort"
read_name_SVA_sort = f"{read_name_SVA}_sort"
read_name_all = f"{working_dir}/read_names/read_name_all_{sample}"
read_name_all_sort = f"{read_name_all}_sort"
sno_outfile = f"{working_dir}/results/results_sep_sno_{sample}"
no_sno_outfile = f"{working_dir}/results/results_sep_non_sno_{sample}"

# Get the read names have one end mapped to U6, the other end mapped to L1
# Sort the read names in the read name list
sam_inf = open(whole_sam_input)
other_read_li = []
LINE_read_li = []
Alu_read_li = []
SVA_read_li = []
for sam_line in sam_inf:
    if re.match("^@",sam_line):
        continue
    else:
        element = sam_line.rstrip().split("\t")
        chrm1 = element[2]
        chrm2 = element[6]
        if chrm2 == "=":
            continue
        else:
            if chrm1 in sam_search_fused or chrm2 in sam_search_fused:
                read_name = element[0]
                if chrm1 == "HUMTNL22" or chrm2 == "HUMTNL22":
                    LINE_read_li.append(read_name)
                elif chrm1 == "AluY" or chrm2 == "AluY":
                    Alu_read_li.append(read_name)
                elif chrm1 == "SVAA" or chrm2 == "SVAA":
                    SVA_read_li.append(read_name)
                else:
                    other_read_li.append(read_name)
            else:
                continue
sam_inf.close()
other_read_li = sorted(set(other_read_li))
LINE_read_li = sorted(set(LINE_read_li))
Alu_read_li = sorted(set(Alu_read_li))
SVA_read_li = sorted(set(SVA_read_li))
print (len(other_read_li))
print (len(LINE_read_li))
print (len(Alu_read_li))
print (len(SVA_read_li))
sam_inf.close()

# Output the read names and sort them in the same way as the fastq file is sorted, push the sorted read names into a new list for read extraction from the fastq files
def read_name_output(readname_list,read_name_outfile,read_name_outfile_sort):
    read_name_outf = open(read_name_outfile,"w")
    print (readname_list,sep="\n",file=read_name_outf)
    read_name_outf.close()
    os.system(f"sort {read_name_outfile}>{read_name_outfile_sort}")
    read_name_list=[]
    read_name_inf = open(read_name_outfile_sort)
    for line in read_name_inf:
        read_name_list.append(line.rstrip())
    read_name_inf.close()
    return read_name_list
read_name_other_li=read_name_output(other_read_li,read_name_other,read_name_other_sort)
read_name_LINE_li=read_name_output(LINE_read_li,read_name_LINE,read_name_LINE_sort)
read_name_Alu_li=read_name_output(Alu_read_li,read_name_Alu,read_name_Alu_sort)
read_name_SVA_li=read_name_output(SVA_read_li,read_name_SVA,read_name_SVA_sort)

def add_identity(datatype,infile):
    inf = open(infile)
    for line in inf:
        read_name_all_li.append([line.rstrip(),datatype])
    inf.close()
read_name_all_li = []
add_identity("other",read_name_other_sort)
add_identity("LINE",read_name_LINE_sort)
add_identity("Alu",read_name_Alu_sort)
outf = open(read_name_all,"w")
for names in read_name_all_li:
    print (names,sep="\t",file=outf)
outf.close()
os.system(f"sort -k1,1 {read_name_all}>{read_name_all_sort}")
read_name_all_li=[]
inf = open(read_name_all_sort)
for line in inf:
    read_name_all_li.append([line.rstrip().split("\t")[0],line.rstrip().split("\t")[1]])
inf.close()

# Get all the reads in the read name list from the original fastq file
def get_reads_fastq(fastq_file,read_li,out_fastq_other,out_fastq_LINE,out_fastq_Alu):
    limit = len(read_li)
    print limit
    fastq_inf = gzip.open(fastq_file)
    other_out = open(out_fastq_other,"w")
    LINE_out = open(out_fastq_LINE,"w")
    Alu_out =open(out_fastq_Alu,"w")
    if read_li == []:
        fastq_inf.close()
        fastq_out.close()
    else:
        i = 0
        while True:
            line1 = fastq_inf.readline().rstrip()
            if not line1: 
                break
            line2 = fastq_inf.readline().rstrip()
            line3 = fastq_inf.readline().rstrip()
            line4 = fastq_inf.readline().rstrip()
            name = line1.split(" ")[0]
            if name == "@"+read_li[i][0]:
                if read_li[i][1] == "other":
                    print ([name,line2,line3,line4],sep="\n",file=other_out)
                elif read_li[i][1] == "LINE":
                    print ([name,line2,line3,line4],sep="\n",file=LINE_out)
                elif read_li[i][1] == "Alu":
                    print ([name,line2,line3,line4],sep="\n",file=Alu_out)
                i += 1
                if i >= limit:
                    break
                if read_li[i][0] == read_li[i-1][0]:
                    i+=1
                
            else:
                continue
    fastq_inf.close()
    other_out.close()
    LINE_out.close()
    Alu_out.close()
get_reads_fastq(fastq_file_1,read_name_all_li,sub_fastq_1,sub_fastq_LINE_1,sub_fastq_Alu_1)
get_reads_fastq(fastq_file_2,read_name_all_li,sub_fastq_2,sub_fastq_LINE_2,sub_fastq_Alu_2)

# Run Flash to get reads that could be re-constructed and paired reads
os.system(f"{flash_install_dir}/flash {sub_fastq_1} {sub_fastq_2} --max-mismatch-density=0.1 --max-overlap=70 -d {flash_dir} -o {sample}_other")
os.system(f"{flash_install_dir}/flash {sub_fastq_LINE_1} {sub_fastq_LINE_2} --max-mismatch-density=0.1 --max-overlap=70 -d {flash_dir} -o {sample}_LINE")
os.system(f"{flash_install_dir}/flash {sub_fastq_Alu_1} {sub_fastq_Alu_2} --max-mismatch-density=0.1 --max-overlap=70 -d {flash_dir} -o {sample}_Alu")
os.system(f"{flash_install_dir}/flash {sub_fastq_SVA_1} {sub_fastq_SVA_2} --max-mismatch-density=0.1 --max-overlap=70 -d {flash_dir} -o {sample}_SVA")

# Align single / paired reads to un-masked refrence
def single_alignemnt(data_type):
    os.system(f"bwa mem {non_masked_ref} {flash_dir}{sample}_{data_type}.extendedFrags.fastq > {alignment_dir}{sample}_{data_type}_single.sam")
    os.system(f"samtools view -bhS {alignment_dir}{sample}_{data_type}_single.sam > {alignment_dir}{sample}_{data_type}_single.bam")
    os.system(f"samtools sort {alignment_dir}{sample}_{data_type}_single.bam {alignment_dir}{sample}_{data_type}_single_sort")
    os.system(f"samtools rmdup {alignment_dir}{sample}_{data_type}_single_sort.bam {alignment_dir}{sample}_{data_type}_single_sort_rmdup_1.bam")
    os.system(f"java -jar {picard_dir}/picard.jar MarkDuplicates REMOVE_DUPLICATES=true `AS=true INPUT={alignment_dir}{sample}_{data_type}_single_sort_rmdup_1.bam OUTPUT={alignment_dir}{sample}_{data_type}_single_sort_rmdup.bam METRICS_FILE={alignment_dir}{sample}_{data_type}_single_metrics.txt")
    os.system(f"samtools index {alignment_dir}{sample}_{data_type}_single_sort_rmdup.bam")
    os.system(f"rm {alignment_dir}{sample}_{data_type}_single.bam")
single_alignemnt("other")
single_alignemnt("LINE")
single_alignemnt("Alu")
single_alignemnt("SVA")

# Separate the matched and unmatched files from each alignment sam files
# Get the list with read names,chromosome,position mapped and cigar of this read
def single_cigar_read(cigar):
    letters = re.split("\d",cigar)
    letters = filter(None,letters)
    numbers = re.split("\D",cigar)
    numbers = filter(None,numbers)
    if len(letters) == 1 and letters[0]=="M":
        return True
    else:
        return False
def single_separate(data_type):
    single_read_matched_li = []
    single_read_nonmatched_li=[]
    single_in = open(f"{alignment_dir}{sample}_{data_type}_single.sam")
    for single_line in single_in:
        if re.match("^@",single_line):
            continue
        else:
            cigar = single_line.rstrip().split("\t")[5]
            chrm = single_line.rstrip().split("\t")[2]
            pos = single_line.rstrip().split("\t")[3]
            read_name = single_line.rstrip().split("\t")[0]
            if single_cigar_read(cigar):
                single_read_matched_li.append([read_name,chrm,pos,cigar])
            else:
                single_read_nonmatched_li.append([read_name,chrm,pos,cigar])
    single_in.close()
    return single_read_matched_li,single_read_nonmatched_li
def match_length(cigar):
    letters = re.split("\d",cigar)
    letters = filter(None,letters)
    numbers = re.split("\D",cigar)
    numbers = filter(None,numbers)
    matched_length = 0
    pos = 0
    #print cigar
    for i in range(0,len(letters)):
        letter = letters[i]
        number = int(numbers[i])
        if letter == "M":
            matched_length += number
            pos += number
        else:
            pos += number
    return matched_length
single_read_matched_other_li,single_read_nonmatched_other_li=single_separate("other")
single_read_matched_LINE_li,single_read_nonmatched_LINE_li=single_separate("LINE")
single_read_matched_Alu_li,single_read_nonmatched_Alu_li=single_separate("Alu")
single_read_matched_SVA_li,single_read_nonmatched_SVA_li=single_separate("SVA")

# According to where the read maps, find the genes which contain this read, also, get rid of U6 reads in this step
def not_U6(chrm,start,end):
    U6 = os.popen(f"tabix {gtf_input} {chrm}:{start}-{end}")
    U6 = U6.read().rstrip().split("\n")
    for i in range(0,len(U6)):
        if gtf_search_fused_li[0] in U6[i] or gtf_search_fused_li[1] in U6[i]:
            return False
        else:
            continue
    if i == len(U6)-1:
        return True
def find_gene(read_list):
    gene_dict = {}
    gene_dict["no_gene"] = []
    for read in read_list:
        read_name = read[0]
        chrm = read[1]
        start=int(read[2])
        cigar=read[3]
        length = match_length(cigar)
        end = start+int(length)
        gene = os.popen(f"tabix {gtf_input} {chrm}:{start}-{end} | awk \'$3==\"gene\"\'")
        gene = gene.read().rstrip().split("\n")
        if cigar == "*":
            continue
        else:
            if len(gene) == 1:
                if gene == ['']:
                    gene_dict["no_gene"].append(read_name)
                else:
                    if not_U6(chrm,start,end):
                        gene_id = gene[0].split("\t")[8].split(";")[0].split("\"")[1]
                        if gene_id in gene_dict.keys():
                            gene_dict[gene_id].append(read_name)
                        else:
                            chrm_gene = gene[0].split("\t")[0]
                            start_gene = gene[0].split("\t")[3]
                            end_gene = gene[0].split("\t")[4]
                            gene_dict[gene_id] = [chrm_gene,start_gene,end_gene,read_name]
                    else:
                        continue
            else:
                for onegene in gene:
                    if not_U6(chrm,start,end):
                        gene_id = onegene.split("\t")[8].split(";")[0].split("\"")[1]
                        if gene_id in gene_dict.keys():
                            gene_dict[gene_id].append(read_name)
                        else:
                            chrm_gene = onegene.split("\t")[0]
                            start_gene = onegene.split("\t")[3]
                            end_gene = onegene.split("\t")[4]
                            gene_dict[gene_id] = [chrm_gene,start_gene,end_gene,read_name]
                    else:
                        continue
    for keys in gene_dict:
        sorted(set(gene_dict[keys]))
    print len(gene_dict)
    return gene_dict
single_other_matched_dict = find_gene(single_read_matched_other_li)
single_other_nonmatched_dict = find_gene(single_read_nonmatched_other_li)

# Count the number of reads belong to each gene, and normalize the number by FPKM called by tophat and cufflinks
def find_FPKM(keys):
    fpkm_line = os.popen(f"grep \"{keys}\" {fpkm_file}")
    fpkm_line= fpkm_line.read().rstrip().split("\t")
    fpkm = float(fpkm_line[9])
    ori_fpkm = fpkm
    if fpkm == 0:
        fpkm = 0.00001
    return ori_fpkm,fpkm
def count_number(genes_dict):
    result_li = []
    for keys in genes_dict:
        if keys != "no_gene":
            number = len(genes_dict[keys])-3
            fpkm_ori,fpkm=find_FPKM(keys)
            normalized_number = float(number)/float(fpkm)
            one_gene = [keys,str(number),normalized_number,fpkm_ori]
            for reads in genes_dict[keys]:
                one_gene.append(reads)
            result_li.append(one_gene)
    return result_li
single_count_other_matched=count_number(single_other_matched_dict)
single_count_other_nonmatched=count_number(single_other_nonmatched_dict)

# Calculate the FPKM for Mobile elements and count the number of reads
def FPKM_ME(ME_name,ME_length):
    read_count = os.popen(f"samtools view -F 4 {whole_bam_input} {ME_name} | wc -l")
    read_count = read_count.read()
    FPKM_ME = float(read_count)*1000000000/(float(ME_length)*float(total_read_count))
    return FPKM_ME
total_read_count = os.popen(f"samtools view -F 4 {whole_bam_input} | wc -l")
total_read_count = total_read_count.read()
FPKM_LINE=FPKM_ME("HUMTNL22",6024)
FPKM_Alu=FPKM_ME("AluY",282)
FPKM_SVA=FPKM_ME("SVAA",1387)
def count_number_ME(reads_li):
    read_name_li = []
    for i in reads_li:
        read_name_li.append(i[0])
    return len(set(sorted(read_name_li)))
single_number_LINE_matched = count_number_ME(single_read_matched_LINE_li)
single_number_Alu_matched = count_number_ME(single_read_matched_Alu_li)
single_number_SVA_matched = count_number_ME(single_read_matched_SVA_li)
single_number_LINE_nonmatched = count_number_ME(single_read_nonmatched_LINE_li)
single_number_Alu_nonmatched = count_number_ME(single_read_nonmatched_Alu_li)
single_number_SVA_matched = count_number_ME(single_read_matched_SVA_li)

# Output everything: gene name/ME name, number of reads,FPKM, counts normalized by FPKM, read names
#Need to add SVA parameters into the funtion if add SVA later
def output(infor_list,out_file,LINE,Alu,SVA,LINE_li,Alu_li,SVA_li):
    outf = open(out_file,"w")
    for items in infor_list:
        gene_name = items[0]
        number = items[1]
        nornalized_number = items[2]
        FPKM=items[3]
        chrm_gene = items[4]
        start_gene = items[5]
        end_gene = items[6]
        reads_li = []
        for t in range(7,len(items)):
            reads_li.append(items[t])
        print ([gene_name,str(chrm_gene),str(start_gene),str(end_gene),str(number)]+reads_li,sep="\t",file=outf)
    LINE_reads = []
    for t in single_read_matched_LINE_li:
        LINE_reads.append(t[0])
    print (["LINE","chrm_LINE","*","*",str(LINE),str(float(LINE)/float(FPKM_LINE)),str(FPKM_LINE)]+LINE_reads,sep="\t",file=outf)
    Alu_reads = []
    for t in single_read_matched_Alu_li:
        Alu_reads.append(t[0])
    print (["Alu","chrm_Alu","*","*",str(Alu),str(float(Alu)/float(FPKM_Alu)),str(FPKM_Alu)]+Alu_reads,sep="\t",file=outf)
    for t in single_read_matched_SVA_li:
        SVA_reads.append(t[0])
    print (["SVA","chrm_SVA","*","*",str(SVA),str(float(SVA)/float(FPKM_SVA)),str(FPKM_SVA)]+SVA_reads,sep="\t",file=outf)
    outf.close()
output(single_count_other_matched,single_matched_out,single_number_LINE_matched,single_number_Alu_matched,single_number_SVA_matched,single_read_matched_LINE_li,single_read_matched_Alu_li,single_read_matched_SVA_li)
output(single_count_other_nonmatched,single_nonmatched_out,single_number_LINE_nonmatched,single_number_Alu_nonmatched,single_number_SVA_nonmatched,single_read_nonmatched_LINE_li,single_read_nonmatched_Alu_li,single_read_nonmatched_SVA_li)

