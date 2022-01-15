#!/bin/bash


######################
# Begin work section #
######################


echo Starting Time is `date "+%Y-%m-%d %H:%M:%S"`
start=$(date +%s)

## main setting
ref_fa=/genome/sequence/genome.fa
chrom_sizes=/chrom/sizes/hg38.chrom.sizes
star_index_name=GRCh38
min_coverage=10
min_var_freq=0.05
validation=1
ncpus=24

## input merge pileup setting
input_pileup_dir=/pileup/out/dir
input_merge_sample_name=input_merge_name

input_plus_pileup=$input_pileup_dir/$input_merge_sample_name/$input_merge_sample_name.plus.pileup
input_minus_pileup=$input_pileup_dir/$input_merge_sample_name/$input_merge_sample_name..minus.pileup

## FTO+ merge pileup setting
ftop_pileup_dir=/pileup/out/dir
ftop_merge_sample_name=ftop_merge_name

ftop_plus_pileup=$ftop_pileup_dir/$ftop_merge_sample_name/$ftop_merge_sample_name.plus.pileup
ftop_minus_pileup=$ftop_pileup_dir/$ftop_merge_sample_name/$ftop_merge_sample_name.minus.pileup


## FTO- rep1 pileup setting
ftom_rep1_sample_name=ftom_rep1_name
ftom_rep1_bam_dir=/star/out/dir
ftom_rep1_varscan_dir=/varscan/out/dir

sample_rep1_bam=$ftom_rep1_bam_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.$star_index_name.align.sorted.bam
sample_rep1_bam_plus=$ftom_rep1_bam_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.$star_index_name.align.sorted.plus.bam
sample_rep1_bam_minus=$ftom_rep1_bam_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.$star_index_name.align.sorted.minus.bam

## FTO- rep1 pileup start
samtools view -hub -F 16 $sample_rep1_bam | samtools sort -T $ftom_rep1_bam_dir -@ $ncpus -o $sample_rep1_bam_plus -
samtools index $sample_rep1_bam_plus
samtools view -hub -f 16 $sample_rep1_bam | samtools sort -T $ftom_rep1_bam_dir -@ $ncpus -o $sample_rep1_bam_minus -
samtools index $sample_rep1_bam_minus

if [[ ! -d $ftom_rep1_varscan_dir/$ftom_rep1_sample_name ]];then mkdir -p $ftom_rep1_varscan_dir/$ftom_rep1_sample_name;fi
samtools mpileup -f $ref_fa $sample_rep1_bam_plus -o $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.plus.pileup
samtools mpileup -f $ref_fa $sample_rep1_bam_minus -o $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.minus.pileup

## FTO- rep1 vs input varscan
java -jar $(which varScan.jar) somatic $input_plus_pileup $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.plus.pileup $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.input.plus.varscan --min-var-freq $min_var_freq --min-coverage $min_coverage --validation $validation
java -jar $(which varScan.jar) somatic $input_minus_pileup $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.minus.pileup $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.input.minus.varscan --min-var-freq $min_var_freq --min-coverage $min_coverage --validation $validation

## FTO- rep1 vs FTO+ varscan
java -jar $(which varScan.jar) somatic $ftop_plus_pileup $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.plus.pileup $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.ftop.plus.varscan --min-var-freq $min_var_freq --min-coverage $min_coverage --validation $validation
java -jar $(which varScan.jar) somatic $ftop_minus_pileup $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.minus.pileup $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.ftop.minus.varscan --min-var-freq $min_var_freq --min-coverage $min_coverage --validation $validation

rm $sample_rep1_bam_plus*
rm $sample_rep1_bam_minus*
rm $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.plus.pileup
rm $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.minus.pileup

echo -e "\nFTO- rep1 pileup and varscan done\n"


## FTO- rep2 pileup setting
ftom_rep2_sample_name=ftom_rep2_name
ftom_rep2_bam_dir=/star/out/dir
ftom_rep2_varscan_dir=/varscan/out/dir

sample_rep2_bam=$ftom_rep2_bam_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.$star_index_name.align.sorted.bam
sample_rep2_bam_plus=$ftom_rep2_bam_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.$star_index_name.align.sorted.plus.bam
sample_rep2_bam_minus=$ftom_rep2_bam_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.$star_index_name.align.sorted.minus.bam

## FTO- rep2 pileup start
samtools view -hub -F 16 $sample_rep2_bam | samtools sort -T $ftom_rep2_bam_dir -@ $ncpus -o $sample_rep2_bam_plus -
samtools index $sample_rep2_bam_plus
samtools view -hub -f 16 $sample_rep2_bam | samtools sort -T $ftom_rep2_bam_dir -@ $ncpus -o $sample_rep2_bam_minus -
samtools index $sample_rep2_bam_minus

if [[ ! -d $ftom_rep2_varscan_dir/$ftom_rep2_sample_name ]];then mkdir -p $ftom_rep2_varscan_dir/$ftom_rep2_sample_name;fi
samtools mpileup -f $ref_fa $sample_rep2_bam_plus -o $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.plus.pileup
samtools mpileup -f $ref_fa $sample_rep2_bam_minus -o $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.minus.pileup

## FTO- rep2 vs input varscan
java -jar $(which varScan.jar) somatic $input_plus_pileup $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.plus.pileup $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.input.plus.varscan --min-var-freq $min_var_freq --min-coverage $min_coverage --validation $validation
java -jar $(which varScan.jar) somatic $input_minus_pileup $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.minus.pileup $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.input.minus.varscan --min-var-freq $min_var_freq --min-coverage $min_coverage --validation $validation

## FTO- rep2 vs FTO+ varscan
java -jar $(which varScan.jar) somatic $ftop_plus_pileup $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.plus.pileup $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.ftop.plus.varscan --min-var-freq $min_var_freq --min-coverage $min_coverage --validation $validation
java -jar $(which varScan.jar) somatic $ftop_minus_pileup $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.minus.pileup $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.ftop.minus.varscan --min-var-freq $min_var_freq --min-coverage $min_coverage --validation $validation

rm $sample_rep2_bam_plus*
rm $sample_rep2_bam_minus*
rm $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.plus.pileup
rm $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.minus.pileup

echo -e "\nFTO- rep2 pileup and varscan done\n"

rm input_plus_pileup
rm input_minus_pileup
rm ftop_plus_pileup
rm ftop_minus_pileup


## identification of sites
# output format: bed6
# column 4: a|b=c=d=e=f=g=h=i|j. Where a, m6A_site_ID; b, non-mutated counts (control); c, mutated counts (control); d, mutation rate (control); e, non-mutated counts (treatment); f, mutated counts (treatment); g, mutation rate (treatment); h, variant p value; i, somatic p value; j, DRACH motif
# column 5: delta mutation rate

# FTO- rep1 vs input
awk -F '\t' 'BEGIN {OFS="\t"} (ARGIND==1 && $3=="A") || (ARGIND==2 && $3=="T") {a=$6/($5+$6)*100;b=$10/($9+$10)*100;if($5+$6>5 && $9+$10>5 && a<5 && ((b-a>5 && $15<0.05) || ($10>=5 && b-a>10))) {if(ARGIND==1) {s="+";l=$1"_"$2"_"s} else {s="-";l=$1"_"$2"_"s};sub(/%/,"",$7);sub(/%/,"",$11);print $1,$2-1,$2,l"|"$5"="$6"="$7"="$9"="$10"="$11"="$14"="$15,b-a,s}}' \
	$ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.input.plus.varscan $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.input.minus.varscan | \
	slopBed -b 2 -i - -g $chrom_sizes | fastaFromBed -s -bedOut -fi $ref_fa -bed - | awk -F '\t' 'BEGIN {OFS="\t"} {a=toupper($7);gsub(/T/,"U",a);if(a~/[GAU][AG]AC[ACU]/) {m="DRACH"} else {m="nonDRACH"};if(m=="DRACH") print $1,$2+2,$3-2,$4"|"a,$5,$6}' | sort -k 1,1V -k 2,2n > $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.input.DRACH.bed
# FTO- rep2 vs input
awk -F '\t' 'BEGIN {OFS="\t"} (ARGIND==1 && $3=="A") || (ARGIND==2 && $3=="T") {a=$6/($5+$6)*100;b=$10/($9+$10)*100;if($5+$6>5 && $9+$10>5 && a<5 && ((b-a>5 && $15<0.05) || ($10>=5 && b-a>10))) {if(ARGIND==1) {s="+";l=$1"_"$2"_"s} else {s="-";l=$1"_"$2"_"s};sub(/%/,"",$7);sub(/%/,"",$11);print $1,$2-1,$2,l"|"$5"="$6"="$7"="$9"="$10"="$11"="$14"="$15,b-a,s}}' \
	$ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.input.plus.varscan $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.input.minus.varscan | \
	slopBed -b 2 -i - -g $chrom_sizes | fastaFromBed -s -bedOut -fi $ref_fa -bed - | awk -F '\t' 'BEGIN {OFS="\t"} {a=toupper($7);gsub(/T/,"U",a);if(a~/[GAU][AG]AC[ACU]/) {m="DRACH"} else {m="nonDRACH"};if(m=="DRACH") print $1,$2+2,$3-2,$4"|"a,$5,$6}' | sort -k 1,1V -k 2,2n > $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.input.DRACH.bed
# FTO- rep1 vs FTO+
awk -F '\t' 'BEGIN {OFS="\t"} (ARGIND==1 && $3=="A") || (ARGIND==2 && $3=="T") {a=$6/($5+$6)*100;b=$10/($9+$10)*100;if($5+$6>5 && $9+$10>5 && a<10 && ((b-a>5 && $15<0.1) || ($10>=5 && b-a>10))) {if(ARGIND==1) {s="+";l=$1"_"$2"_"s} else {s="-";l=$1"_"$2"_"s};sub(/%/,"",$7);sub(/%/,"",$11);print $1,$2-1,$2,l"|"$5"="$6"="$7"="$9"="$10"="$11"="$14"="$15,b-a,s}}' \
	$ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.ftop.plus.varscan $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.ftop.minus.varscan | \
	slopBed -b 2 -i - -g $chrom_sizes | fastaFromBed -s -bedOut -fi $ref_fa -bed - | awk -F '\t' 'BEGIN {OFS="\t"} {a=toupper($7);gsub(/T/,"U",a);if(a~/[GAU][AG]AC[ACU]/) {m="DRACH"} else {m="nonDRACH"};if(m=="DRACH") print $1,$2+2,$3-2,$4"|"a,$5,$6}' | sort -k 1,1V -k 2,2n > $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.ftop.DRACH.bed
# FTO- rep2 vs FTO+
awk -F '\t' 'BEGIN {OFS="\t"} (ARGIND==1 && $3=="A") || (ARGIND==2 && $3=="T") {a=$6/($5+$6)*100;b=$10/($9+$10)*100;if($5+$6>5 && $9+$10>5 && a<10 && ((b-a>5 && $15<0.1) || ($10>=5 && b-a>10))) {if(ARGIND==1) {s="+";l=$1"_"$2"_"s} else {s="-";l=$1"_"$2"_"s};sub(/%/,"",$7);sub(/%/,"",$11);print $1,$2-1,$2,l"|"$5"="$6"="$7"="$9"="$10"="$11"="$14"="$15,b-a,s}}' \
	$ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.ftop.plus.varscan $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.ftop.minus.varscan | \
	slopBed -b 2 -i - -g $chrom_sizes | fastaFromBed -s -bedOut -fi $ref_fa -bed - | awk -F '\t' 'BEGIN {OFS="\t"} {a=toupper($7);gsub(/T/,"U",a);if(a~/[GAU][AG]AC[ACU]/) {m="DRACH"} else {m="nonDRACH"};if(m=="DRACH") print $1,$2+2,$3-2,$4"|"a,$5,$6}' | sort -k 1,1V -k 2,2n > $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.ftop.DRACH.bed

## filtering of sites
# output format: bed6
# column 4: a=b=c=d=e. Where a, m6A_site_ID; b, mutation rate (rep1); c, mutation rate (rep2); d, p value min; e, DRACH motif
# column 5: mutation rate mean

# Alternative set 1 (without FTO+ controls):
intersectBed -wo -s -a $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.input.DRACH.bed -b $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.input.DRACH.bed | \
	awk -F '\t' 'BEGIN {OFS="\t"} {split($4,a,"|");split($10,b,"|");split(a[2],x,"=");split(b[2],y,"=");if(x[8]<y[8]) {p=x[8]} else {p=y[8]};m=($5+$11)/2;print $1,$2,$3,a[1]"="$5"="$11"="p"="a[3],m,$6}' > $ftom_rep1_varscan_dir/$ftom_rep1_sample_name.$ftom_rep2_sample_name.input.common.DRACH.bed
# Alternative set 2 (within any one of two FTO+ controls):
cat $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.ftop.DRACH.bed $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.ftop.DRACH.bed | \
	intersectBed -u -s -a $ftom_rep1_varscan_dir/$ftom_rep1_sample_name.$ftom_rep2_sample_name.input.common.DRACH.bed -b - > $ftom_rep1_varscan_dir/$ftom_rep1_sample_name.$ftom_rep2_sample_name.input.common.ftop.any.DRACH.bed
# Alternative set 3 (within both FTO+ controls):
intersectBed -u -s -a $ftom_rep1_varscan_dir/$ftom_rep1_sample_name/$ftom_rep1_sample_name.ftop.DRACH.bed -b $ftom_rep2_varscan_dir/$ftom_rep2_sample_name/$ftom_rep2_sample_name.ftop.DRACH.bed | \
	intersectBed -u -s -a $ftom_rep1_varscan_dir/$ftom_rep1_sample_name.$ftom_rep2_sample_name.input.common.DRACH.bed -b - > $ftom_rep1_varscan_dir/$ftom_rep1_sample_name.$ftom_rep2_sample_name.input.common.ftop.both.DRACH.bed

echo -e "\nAll done\n"

echo Ending Time is `date "+%Y-%m-%d %H:%M:%S"`
end=$(date +%s)
time=$(( ($end - $start) / 60 ))
echo Used Time is $time mins
