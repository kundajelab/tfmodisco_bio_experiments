bedtools intersect -a $TFNET_ROOT/results/nandi/bQTL/analysis/bQTL_all_SNPs/JUND.tsv -b /srv/scratch/ktian/kundajelab/tfmodisco_bio_experiments/results/nandi/bQTL/bQTL_Peaks/unzip/LCL-JUND-human-EXP0-optimal_idr.narrowPeak > filter.bed
#bedtools intersect -a filter0.bed -b dnase.bed > filter.bed
cut -f 4- filter.bed | python $TFNET_ROOT/scripts/snp_to_bed.py 1000 > filter.tsv
head -n 3569 filter.tsv > filter_top.tsv
tail -n 3569 filter.tsv > filter_bot.tsv
cat filter_top.tsv filter_bot.tsv > filter_top_bot.tsv
cut -f 4- filter_top_bot.tsv > filter_top_bot.txt
echo header > header.txt
cat filter_top_bot.txt >> header.txt
mv header.txt filter_top_bot.txt 
#ln -s JUND_20k.tsv interpret.tsv
ln -s filter_top_bot.tsv interpret.tsv
bedtools getfasta -fi $TFNET_ROOT/genome/hg19.fa -bed interpret.tsv -fo interpret.fa

mkdir scores
python $TFNET_ROOT/scripts/gen_allele.py

mv interpret.fa interpret.fa0
mv interpret.fa* scores
ln -s scores/interpret.fa0 interpret.fa
mkdir logs

ln -s ../../../bQTL_ChipSeq_JUND_18_08_29/model_files .
#ln -s /home/ktian/kundajelab/tfnet/results/nandi/JUND/JUND_GM12878_refine_18_09_04/model_files .
#source $TFNET_ROOT/scripts/run_allele_deeplift.sh
source $TFNET_ROOT/scripts/run_allele_deeplift.sh 0

#python get_scores.py
