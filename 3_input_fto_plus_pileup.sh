#!/bin/bash


######################
# Begin work section #
######################


echo Starting Time is `date "+%Y-%m-%d %H:%M:%S"`
start=$(date +%s)

## main setting
ref_fa=/genome/sequence/genome.fa
star_index_name=GRCh38
ncpus=24


## input pileup setting
input_rep1_sample_name=input_rep1_name
input_rep2_sample_name=input_rep2_name
input_merge_sample_name=input_merge_name
input_bam_dir=/star/out/dir
input_pileup_dir=/pileup/out/dir

input_merge_sample_bam=$input_bam_dir/$input_merge_sample_name/$input_merge_sample_name.$star_index_name.align.sorted.bam
input_merge_sample_bam_plus=$input_bam_dir/$input_merge_sample_name/$input_merge_sample_name.$star_index_name.align.sorted.plus.bam
input_merge_sample_bam_minus=$input_bam_dir/$input_merge_sample_name/$input_merge_sample_name.$star_index_name.align.sorted.minus.bam


## input pileup start
if [[ ! -d $input_bam_dir/$input_merge_sample_name ]];then mkdir -p $input_bam_dir/$input_merge_sample_name;fi
samtools merge -f $input_merge_sample_bam $input_bam_dir/$input_rep1_sample_name/$input_rep1_sample_name.$star_index_name.align.sorted.bam $input_bam_dir/$input_rep2_sample_name/$input_rep2_sample_name.$star_index_name.align.sorted.bam
samtools index $input_merge_sample_bam

samtools view -hub -F 16 $input_merge_sample_bam | samtools sort -T $input_bam_dir -@ $ncpus -o $input_merge_sample_bam_plus -
samtools index $input_merge_sample_bam_plus
samtools view -hub -f 16 $input_merge_sample_bam | samtools sort -T $input_bam_dir -@ $ncpus -o $input_merge_sample_bam_minus -
samtools index $input_merge_sample_bam_minus

if [[ ! -d $input_pileup_dir/$input_merge_sample_name ]];then mkdir -p $input_pileup_dir/$input_merge_sample_name;fi
samtools mpileup -f $ref_fa $input_merge_sample_bam_plus -o $input_pileup_dir/$input_merge_sample_name/$input_merge_sample_name.plus.pileup
samtools mpileup -f $ref_fa $input_merge_sample_bam_minus -o $input_pileup_dir/$input_merge_sample_name/$input_merge_sample_name.minus.pileup

rm $input_merge_sample_bam_plus*
rm $input_merge_sample_bam_minus*

echo -e "\ninput pileup done\n"


## FTO+ pileup setting
ftop_rep1_sample_name=ftop_rep1_name
ftop_rep2_sample_name=ftop_rep2_name
ftop_merge_sample_name=ftop_merge_name
ftop_bam_dir=/star/out/dir
ftop_pileup_dir=/pileup/out/dir

ftop_merge_sample_bam=$ftop_bam_dir/$ftop_merge_sample_name/$ftop_merge_sample_name.$star_index_name.align.sorted.bam
ftop_merge_sample_bam_plus=$ftop_bam_dir/$ftop_merge_sample_name/$ftop_merge_sample_name.$star_index_name.align.sorted.plus.bam
ftop_merge_sample_bam_minus=$ftop_bam_dir/$ftop_merge_sample_name/$ftop_merge_sample_name.$star_index_name.align.sorted.minus.bam


## FTO+ pileup start
if [[ ! -d $ftop_bam_dir/$ftop_merge_sample_name ]];then mkdir -p $ftop_bam_dir/$ftop_merge_sample_name;fi
samtools merge -f $ftop_merge_sample_bam $ftop_bam_dir/$ftop_rep1_sample_name/$ftop_rep1_sample_name.$star_index_name.align.sorted.bam $ftop_bam_dir/$ftop_rep2_sample_name/$ftop_rep2_sample_name.$star_index_name.align.sorted.bam
samtools index $ftop_merge_sample_bam

samtools view -hub -F 16 $ftop_merge_sample_bam | samtools sort -T $ftop_bam_dir -@ $ncpus -o $ftop_merge_sample_bam_plus -
samtools index $ftop_merge_sample_bam_plus
samtools view -hub -f 16 $ftop_merge_sample_bam | samtools sort -T $ftop_bam_dir -@ $ncpus -o $ftop_merge_sample_bam_minus -
samtools index $ftop_merge_sample_bam_minus

if [[ ! -d $ftop_pileup_dir/$ftop_merge_sample_name ]];then mkdir -p $ftop_pileup_dir/$ftop_merge_sample_name;fi
samtools mpileup -f $ref_fa $ftop_merge_sample_bam_plus -o $ftop_pileup_dir/$ftop_merge_sample_name/$ftop_merge_sample_name.plus.pileup
samtools mpileup -f $ref_fa $ftop_merge_sample_bam_minus -o $ftop_pileup_dir/$ftop_merge_sample_name/$ftop_merge_sample_name.minus.pileup

rm $ftop_merge_sample_bam_plus*
rm $ftop_merge_sample_bam_minus*

echo -e "\nFTO+ pileup done\n"

echo Ending Time is `date "+%Y-%m-%d %H:%M:%S"`
end=$(date +%s)
time=$(( ($end - $start) / 60 ))
echo Used Time is $time mins
