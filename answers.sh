#! /usr/bin/env bash 

bed="$HOME/data/MOLB7621/data-sets/bed"
newbed="$HOME/data/MOLB7621/data-sets/new_bed"
problemset="$HOME/data/MOLB7621/problem-set-2"

# Question 1:  Use BEDtools intersect to identify 
# the size of the largest overlap between CTCF and H3K4me3 locations.

zcat $bed/encode.tfbs.chr22.bed.gz \
    | grep -w 'CTCF' | sort -k2,2n \
    > $newbed/CTCF_sorted.bed
zcat $bed/encode.h3k4me3.hela.chr22.bed.gz \
    | cut -f1,2,3 | sort -k2,2n \
    > $newbed/h3k4me3_sorted.bed

answer_1=$(bedtools intersect -wo -a $newbed/CTCF_sorted.bed \
    -b $newbed/h3k4me3_sorted.bed \
    | sort -k8,8rn| cut -f8| head -n1)
echo answer-1: $answer_1

# Question 2: Use BEDtools to calculate the GC content 
# of nucleotides 19,000,000 to 19,000,500 on chr22 of hg19 genome build. 
# Report the GC content as a fraction (e.g., 0.50).

fasta="$HOME/data/MOLB7621/data-sets/fasta"
misc="$HOME/data/MOLB7621/data-sets/misc"

# First, make a bed file with chr22 19000000 19000500

touch 1_line.bed
mv 1_line.bed $newbed
echo "chr22;19000000;19000500" >> $newbed/1_line.bed 
sed 's/;/\t/g' $newbed/1_line.bed > $newbed/oneline.bed
rm -f $newbed/1_line.bed

# Then use getfasta to extract the sequence 
# corresponding to the interval chr22 19000000 19000500
answer_2=$(bedtools getfasta -fi $fasta/chr22_hg19.fa \
    -bed $newbed/oneline.bed -fo 500nt.fa 
    total_number=$(cat $problemset/500nt.fa \
    | grep -v '>' | tr -d '\n'| wc -c) 
    GC_number=$(cat $problemset/500nt.fa \
    | grep -v '>' | sed 's/T//g' | sed 's/A//g' \
    | tr -d '\n'|  wc -c) 
    echo $GC_number/$total_number | bc -l)
echo answer-2: $answer_2

# Question 3
# Use BEDtools to identify the length of the CTCF ChIP-seq peak
# (i.e., interval) that has the largest mean signal in ctcf.hela.chr22.bg.gz

bedtools="$HOME/data/MOLB7621/data-sets/bedtools"

answer_3=$(bedtools merge -i $bedtools/ctcf.hela.chr22.bg.gz \
    -d 0 -c 4 -o mean | sort -k4,4rn | awk '{print $3 -$2}' \
    | head -n1)
echo answer-3: $answer_3

# Question 4
# Use BEDtools to identify the gene promoter 
# (defined as 1000 bp upstream of a TSS) with the highest 
# median signal in ctcf.hela.chr22.bg.gz. Report the gene name (e.g., 'ABC123')
# bedtools flank [OPTIONS] -i <BED/GFF/VCF> -g <GENOME> [-b or (-l and -r)]

genome="$HOME/data/MOLB7621/data-sets/genome"

answer_4=$(bedtools flank -i $bed/tss.hg19.chr22.bed.gz \
    -g $genome/hg19.genome -l 1000 -r 0 -s \
    > $newbed/genes.1kb.promoters.bed 
    bedtools sort -i $newbed/genes.1kb.promoters.bed \
    > $newbed/sorted_genes_1kb_promoters.bed
    bedtools map -c 4 -o median -a $newbed/sorted_genes_1kb_promoters.bed \
   -b $bedtools/ctcf.hela.chr22.bg.gz | sort -k7,7nr \
   | cut -f 4 | head -n1) 

echo answer-4: $answer_4

# Question 5
# Use BEDtools to identify the longest interval on chr22 that is not covered by genes.hg19.bed.gz. Report the interval like chr1:100-500.

cat $genome/hg19.genome | grep chr22 \
    | awk 'BEGIN {OFS="\t"} {print $1,0,$2}' > $newbed/whole_chr22.bed
answer_5=$(bedtools subtract -a $newbed/whole_chr22.bed \
    -b $bed/genes.hg19.bed.gz \
    | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$3-$2}' \
    | sort -k4,4rn | head -n1 \
    | awk 'BEGIN {OFS=""} {print $1,":",$2,"-",$3}')  
    
echo answer-5: $answer_5

# Question 6
# Use one or more BEDtools that we haven't covered in class. Be creative.

# Use bedtools 'closest' to identify the gene with the highest CTCF peak
# (the peaks are defined using bedtools merge)
# in its promoter (the promoter is defined by the 1000 bp upstream of 
# its transcription start sites) 
# File used: data-sets/bedtools/ctcf.hela.chr22.bg.gz

zcat $bed/genes.hg19.bed.gz | sort -k1,1 -k2,2n > $newbed/sorted_genes.hg19.bed
bedtools merge -i $bedtools/ctcf.hela.chr22.bg.gz -d 0 -c 4 -o median \
    | sort -k1,1 -k2,2n > $newbed/sorted_ctcf_chr22_peaks.bed

# Here I removed the 5th column of the sorted_genes.hg19.bed file to help
# me keep only the data that would help me select the right data next.

cat $newbed/sorted_genes.hg19.bed | cut -f1,2,3,4,6 \
    > $newbed/noexondata_sorted_genes.hg19.bed

answer_6=$(bedtools closest -D b -a $newbed/sorted_ctcf_chr22_peaks.bed \
    -b $newbed/noexondata_sorted_genes.hg19.bed \
    | awk 'BEGIN {OFS="\t"} ($10<0) && ($10 >= -1000) {print $0}' \
    | sort -k4,4rn | cut -f 8| head -n1)

echo answer-6: $answer_6

# Question 7
# How many genes on the chromosome 22 have at least one CTCF ChIP peak
# 1000 bp upstream of their transcription starting site and another CTCF
# ChIP peak 1000 bp downstream of their gene end?
# Or how many genes may potentially be contained in a chromatin loop
# characterized by anchors of the chromatin loop located at 1000 bp
# from the start and 1000 bp from the end of the gene?

bedtools closest -D b -a $newbed/sorted_ctcf_chr22_peaks.bed \
    -b $newbed/noexondata_sorted_genes.hg19.bed \
    | awk 'BEGIN {OFS="\t"} ($10 <= 1000) && ($10 >= -1000) {print $0}' \
    | sort -k8 | cut -f 8,10 | awk '($2>0) {print $0}' \
    > $newbed/downstream-CTCF

bedtools closest -D b -a $newbed/sorted_ctcf_chr22_peaks.bed \
    -b $newbed/noexondata_sorted_genes.hg19.bed \
    | awk 'BEGIN {OFS="\t"} ($10 <= 1000) && ($10 >= -1000) {print $0}' \
    | sort -k8 | cut -f 8,10 | awk '($2<0) {print $0}' \
    > $newbed/upstream-CTCF

awk '{print $1}' $newbed/upstream-CTCF | uniq -c \
    | awk '{print $2}' > $newbed/unique_upstream_CTCF 

awk '{print $1}' $newbed/downstream-CTCF | uniq -c \
    | awk '{print $2}'> $newbed/unique_downstream_CTCF

cat $newbed/unique_downstream_CTCF $newbed/unique_upstream_CTCF \
    > $newbed/both_up_down_CTCF

answer_7=$(cat $newbed/both_up_down_CTCF | sort | uniq -c \
    | awk '($1!=1) {print $0}' | wc -l) 

echo answer-7: $answer_7




