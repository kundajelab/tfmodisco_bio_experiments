#!/usr/bin/env bash

##filter the SPI1 bQTLs for the ones that are in the 1kb summits
#spi1_bqtls="/home/ktian/kundajelab/tfnet/results/nandi/bQTL/analysis/bQTL_all_SNPs/SPI1.txt"
#spi1_summits_and_flank_coords="/srv/scratch/ktian/kundajelab/tfmodisco_bio_experiments/results/nandi/SPI1/SPI1_GM12878_18_08_31/interpret.tsv"

bedtools intersect -a <(cat $1 | perl -lane 'if ($. > 1) {print $F[0]."\t".($F[1]-1)."\t".($F[1])."\t".$F[5]."\t".$F[6]."\t".$F[9]}') -b $2 -wa -wb


