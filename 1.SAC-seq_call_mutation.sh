#!/bin/bash


bam_file=$1
if [[ $bam_file == *input* ]]; then 
    samtools mpileup --input-fmt-option 'filter=(flag == 99 || flag == 147)' -d 0 -Q 13 --reverse-del \
        -f ~/Genome/hg38_UCSC.fa ./mapping/${bam_file} | python ~/Genome/tools_custom/SAC-seq/scripts/SAC_seq_mpileup_T_input.py \
        > ./Mut_detection/${bam_file/sorted.bam/txt}
    samtools mpileup --input-fmt-option 'filter=(flag == 83 || flag == 163)' -d 0 -Q 13 --reverse-del \
        -f ~/Genome/hg38_UCSC.fa ./mapping/${bam_file} | python ~/Genome/tools_custom/SAC-seq/scripts/SAC_seq_mpileup_A_input.py \
        >> ./Mut_detection/${bam_file/sorted.bam/txt}
fi


if [[ $bam_file == *label* ]]; then
    samtools mpileup --input-fmt-option 'filter=(flag == 99 || flag == 147)' -d 0 -Q 13 --reverse-del \
        -f ~/Genome/hg38_UCSC.fa ./mapping/${bam_file} | python ~/Genome/tools_custom/SAC-seq/scripts/SAC_seq_mpileup_T_label.py \
        > ./Mut_detection/${bam_file/sorted.bam/txt}
    samtools mpileup --input-fmt-option 'filter=(flag == 83 || flag == 163)' -d 0 -Q 13 --reverse-del \
        -f ~/Genome/hg38_UCSC.fa ./mapping/${bam_file} | python ~/Genome/tools_custom/SAC-seq/scripts/SAC_seq_mpileup_A_label.py \
        >> ./Mut_detection/${bam_file/sorted.bam/txt}
fi


