#!/bin/bash

### note: format of .gz file: NHEK_Arsenic_72h_rp1.input.R1.fastq.gz NHEK_Arsenic_72h_rp1.input.R2.fastq.gz
### line 21: the file final location needs to be changed
### how to use: ./SAC-seq_mapping.sh NHEK_Arsenic_72h_rp1.input.R1.fastq.gz NHEK_Arsenic_72h_rp1.input.R2.fastq.gz
R1=$1
R2=$2

for i in ./log ./spike-in ./mapping ./mapping/STAR_log ./Mut_detection ./Mut_detection/Merged
do
    if [ ! -d $i ];then
        mkdir $i
    fi
done

cutadapt -j 15 -m 26 --max-n=0 -e 0.15 -q 20 --nextseq-trim=20 -O 6 --pair-filter=both \
    -a AGATCGGAAGAGCACACGTCTG -A AGATCGGAAGAGCGTCGTGT \
    -o ${R1/.fastq.gz/.trimmed_1.fastq.gz} -p ${R2/.fastq.gz/.trimmed_1.fastq.gz} \
    $R1 $R2 > ./log/${R1/.R1.fastq.gz/}.log

mv $R1 $R2 /home/yutaozhao/data/Yu-Ying

cutadapt -j 15 -m 26 -n 5 -O 12 \
    -g ACACGACGCTCTTCCGATCT -a AGATCGGAAGAGCGTCGTGT \
    -G ACACGACGCTCTTCCGATCT -A AGATCGGAAGAGCGTCGTGT \
    -o ${R1/.fastq.gz/.trimmed_2.fastq.gz} -p ${R2/.fastq.gz/.trimmed_2.fastq.gz} \
    ${R1/.fastq.gz/.trimmed_1.fastq.gz} ${R2/.fastq.gz/.trimmed_1.fastq.gz} >> ./log/${R1/.R1.fastq.gz/}.log
rm ${R1/.fastq.gz/.trimmed_1.fastq.gz} ${R2/.fastq.gz/.trimmed_1.fastq.gz}

clumpify.sh in=${R1/.fastq.gz/.trimmed_2.fastq.gz} in2=${R2/.fastq.gz/.trimmed_2.fastq.gz} \
    out=${R1/.fastq.gz/.dedupe.fastq.gz} out2=${R2/.fastq.gz/.dedupe.fastq.gz} dedupe 2>> ./log/${R1/.R1.fastq.gz/}.log
rm ${R1/.fastq.gz/.trimmed_2.fastq.gz} ${R2/.fastq.gz/.trimmed_2.fastq.gz}

cutadapt -j 15 -m 15 -U 11 -u -11 --rename='{id}_{r2.cut_prefix} {comment}' \
    -o ${R1/.fastq.gz/.barcoded.fastq.gz} -p ${R2/.fastq.gz/.barcoded.fastq.gz} \
    ${R1/.fastq.gz/.dedupe.fastq.gz} ${R2/.fastq.gz/.dedupe.fastq.gz} >> ./log/${R1/.R1.fastq.gz/}.log
rm ${R1/.fastq.gz/.dedupe.fastq.gz} ${R2/.fastq.gz/.dedupe.fastq.gz}

echo "SAC-seq spike_in mapping:" >> ./log/${R1/.R1.fastq.gz/}.log
bowtie2 -p 15 --no-unal --end-to-end -L 16 -N 1 --mp 5 --un-conc-gz ${R1/R1.fastq.gz/after_spike_in.fastq.gz} \
    -x ~/Genome/tools_custom/SAC-seq/spike_in -1 ${R1/.fastq.gz/.barcoded.fastq.gz} -2 ${R2/.fastq.gz/.barcoded.fastq.gz} \
    2>> ./log/${R1/.R1.fastq.gz/}.log | samtools sort -@ 15 --input-fmt-option "filter=[NM]<=10" \
    -O BAM -o ./spike-in/${R1/.R1.fastq.gz/}.spike_in.bam
rm ${R1/.fastq.gz/.barcoded.fastq.gz} ${R2/.fastq.gz/.barcoded.fastq.gz}

echo "SAC-seq 45S mapping:" >> ./log/${R1/.R1.fastq.gz/}.log
bowtie2 -p 15 --no-unal --end-to-end -L 16 -N 1 --mp 5 --un-conc-gz ${R1/R1.fastq.gz/after_rRNA.fastq.gz} \
    -x ~/Genome/hg38-45S/hg38_45S -1 ${R1/R1.fastq.gz/after_spike_in.fastq.1.gz} -2 ${R1/R1.fastq.gz/after_spike_in.fastq.2.gz} \
    >/dev/null 2>> ./log/${R1/.R1.fastq.gz/}.log
rm ${R1/R1.fastq.gz/after_spike_in.fastq.1.gz} ${R1/R1.fastq.gz/after_spike_in.fastq.2.gz}

ulimit -n 1000000
STAR --runThreadN 20 --genomeDir ~/Genome/hg38_UCSC \
    --readFilesIn ${R1/R1.fastq.gz/after_rRNA.fastq.1.gz} ${R1/R1.fastq.gz/after_rRNA.fastq.2.gz} \
    --readFilesCommand gunzip -c --alignEndsType Local --outFilterMatchNminOverLread 0.66 \
    --outFilterMatchNmin 15 --outFilterMismatchNmax 5 --outFilterMismatchNoverLmax 0.2 \
    --outFilterMultimapNmax 50 --outSAMmultNmax -1 --outReadsUnmapped None --limitBAMsortRAM 8000000000 \
    --outSAMtype BAM Unsorted --outFileNamePrefix ./mapping/${R1/.R1.fastq.gz/}_
cat ./mapping/${R1/.R1.fastq.gz/}_Log.final.out >> ./log/${R1/.R1.fastq.gz/}.log
mv ./mapping/${R1/.R1.fastq.gz/}_Log.final.out ./mapping/${R1/.R1.fastq.gz/}_Log.out \
    ./mapping/${R1/.R1.fastq.gz/}_Log.progress.out ./mapping/${R1/.R1.fastq.gz/}_SJ.out.tab ./mapping/STAR_log
rm ${R1/R1.fastq.gz/after_rRNA.fastq.1.gz} ${R1/R1.fastq.gz/after_rRNA.fastq.2.gz}

samtools sort -@ 10 --input-fmt-option 'filter=(flag == 99 || flag == 147 || flag == 83 || flag == 163 )' \
    -o ./mapping/${R1/.R1.fastq.gz/}.sorted.bam ./mapping/${R1/.R1.fastq.gz/}_Aligned.out.bam
rm ./mapping/${R1/.R1.fastq.gz/}_Aligned.out.bam





