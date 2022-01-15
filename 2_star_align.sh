#!/bin/bash


######################
# Begin work section #
######################


echo Starting Time is `date "+%Y-%m-%d %H:%M:%S"`
start=$(date +%s)

## main setting
cutadapt_out_dir=/cutadapt/out/dir
star_index_name=GRCh38
star_index=/star/index
star_out_dir=/star/out/dir
ncpus=24
sample=sampleName
sample_fa_file=$cutadapt_out_dir/$sample.R1.adapter3trim.collapse.fa


echo -e "\n$sample\n"

if [[ -d $star_out_dir/$sample ]];then rm -rf $star_out_dir/$sample;fi
mkdir -p $star_out_dir/$sample

echo -e "\nstar align -- genome mapping\n"
STAR --genomeDir $star_index --readFilesIn $sample_fa_file         \
    --outFileNamePrefix $star_out_dir/$sample/ \
    --runThreadN $ncpus --genomeLoad NoSharedMemory     \
    --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1    \
    --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.06              \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000         \
    --outSAMunmapped None --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
    --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --sjdbScore 1 --limitBAMsortRAM 10000000000 \
    --outWigType bedGraph --outWigStrand Stranded --outSAMstrandField None --outWigReferencesPrefix chr

samtools view -@ $ncpus -F 1548 -Shub $star_out_dir/$sample/Aligned.sortedByCoord.out.bam | samtools sort -T $star_out_dir/$sample -@ $ncpus -o $star_out_dir/$sample/$sample.$star_index_name.align.sorted.bam -

echo -e "\nsamtools index\n"
samtools index $star_out_dir/$sample/$sample.$star_index_name.align.sorted.bam

echo -e "\nsamtools flagstat\n"
samtools flagstat $star_out_dir/$sample/$sample.$star_index_name.align.sorted.bam > $star_out_dir/$sample/$sample.$star_index_name.align.sorted.flagstat.qc

wait

echo Ending Time is `date "+%Y-%m-%d %H:%M:%S"`
end=$(date +%s)
time=$(( ($end - $start) / 60 ))
echo Used Time is $time mins
