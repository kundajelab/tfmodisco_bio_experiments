#!/usr/bin/env bash

#copy over the new K562 dnase file
cp /oak/stanford/groups/akundaje/projects/atlas/dnase_processed/k562_dnase/hg38/cromwell-executions/atac/41bb41c3-c39b-4145-9305-58cc692ef2e4/call-reproducibility_overlap/execution/optimal_peak.narrowPeak.gz hg38_new_k562.narrowPeak.gz

#prep universal dnase negatives - for performance assessment
universal_negatives=/srv/scratch/annashch/deeplearning/gecco/encode_dnase/all_dnase_idr_peaks.sorted.bed
merged_universal_neg="merged_universal_neg.bed.gz"
cat $universal_negatives | mergeBed | gzip -c > $merged_universal_neg
bedtools intersect -sorted -a $merged_universal_neg -b $universal_negatives -wa -wb | gzip -c > merged_universal_neg_intersected.bed.gz
./take_best_peak.py | gzip -c > merged_universal_neg_representative_peaks.bed.gz

#Map the representative DNase negatives to hg38 using liftover
#download the liftover chain file if it doesn't exist already
[[ -e hg19ToHg38.over.chain.gz ]] || wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
#run the liftover command
/software/ucsc_tools/latest/liftOver merged_universal_neg_representative_peaks.bed.gz hg19ToHg38.over.chain.gz hg38_merged_universal_neg_representative_peaks.bed unmapped_merged_universal_neg_representative_peaks.bed
#compress and clean up
gzip -f hg38_merged_universal_neg_representative_peaks.bed
rm unmapped_merged_universal_neg_representative_peaks.bed
rm merged_universal_neg_representative_peaks.bed.gz

#take 2000bp around the summit of the k562 peaks
zcat hg38_new_k562.narrowPeak.gz | perl -lane 'print $F[0]."\t".(($F[1]+$F[9]))."\t".(($F[1]+$F[9]))' | bedtools slop -g hg38.chrom.sizes -b 1000 | perl -lane 'if ($F[2]-$F[1]==2000) {print $_}' | gzip -c > 2000bp_around_hg38_new_k562.narrowPeak.gz
#take 2000bp around the summit of the DNAse peaks
zcat hg38_merged_universal_neg_representative_peaks.bed.gz | perl -lane 'print $F[0]."\t".(($F[1]+$F[9]))."\t".(($F[1]+$F[9]))' | bedtools slop -g hg38.chrom.sizes -b 1000 | perl -lane 'if ($F[2]-$F[1]==2000) {print $_}' | gzip -c > 2000bp_around_hg38_merged_universal_neg_representative_peaks.bed.gz

#filter out any universal dnase that overlap anything in the k562 set 
intersectBed -v -a 2000bp_around_hg38_merged_universal_neg_representative_peaks.bed.gz -b 2000bp_around_hg38_new_k562.narrowPeak.gz -wa | gzip -c > nok562_2000bp_around_hg38_merged_universal_neg_representative_peaks.bed.gz
#clean up the intermediate file
rm 2000bp_around_hg38_merged_universal_neg_representative_peaks.bed.gz

