#!/bin/bash


######################
# Begin work section #
######################


echo Starting Time is `date "+%Y-%m-%d %H:%M:%S"`
start=$(date +%s)

sample=sampleName

cutadapt -e 0.1 -n 1 -O 1 -q 10 -m 32 -a AGATCGGAAGAGCACACGTCT -A GATCGTCGGACTGTAGAACTC -o $sample.R1.adapter3trim.fastq.gz -p $sample.R2.adapter3trim.fastq.gz \
	$sample.R1.fastq.gz $sample.R2.fastq.gz > $sample.adapter3trim.log
zcat $sample.R1.adapter3trim.fastq.gz | fastx_collapser -Q33 -i - -o $sample.R1.adapter3trim.collapse.fa
zcat $sample.R2.adapter3trim.fastq.gz | fastx_collapser -Q33 -i - -o $sample.R2.adapter3trim.collapse.fa
cutadapt -e 0.1 -n 1 -O 1 -m 16 -u 5 -u -11 -o $sample.R1.clean.fa $sample.R1.adapter3trim.collapse.fa > $sample.R1.barcoder5n3trim.log
cutadapt -e 0.1 -n 1 -O 1 -m 16 -u 11 -u -5 -o $sample.R2.clean.fa $sample.R2.adapter3trim.collapse.fa > $sample.R2.barcoder5n3trim.log
rm $sample.R[12].adapter3trim.fastq.gz
rm $sample.R[12].adapter3trim.collapse.fa

echo -e "\nAll done\n"

echo Ending Time is `date "+%Y-%m-%d %H:%M:%S"`
end=$(date +%s)
time=$(( ($end - $start) / 60 ))
echo Used Time is $time mins
