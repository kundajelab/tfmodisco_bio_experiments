cat $TFNET_ROOT/results/nandi/bQTL/analysis/bQTL_all_SNPs/SPI1_20k.txt | python $TFNET_ROOT/scripts/snp_to_bed.py 1000 > SPI1_20k.tsv
ln -s SPI1_20k.tsv interpret.tsv
bedtools getfasta -fi $TFNET_ROOT/genome/hg19.fa -bed interpret.tsv -fo interpret.fa

mkdir scores
python $TFNET_ROOT/scripts/gen_allele.py

mv interpret.fa interpret.fa0
mv interpret.fa* scores
ln -s scores/interpret.fa0 interpret.fa
mkdir logs

#ln -s ../../../bQTL_ChipSeq_SPI1_18_08_29/model_files .
ln -s /home/ktian/kundajelab/tfnet/results/nandi/SPI1/SPI1_GM12878_refine_18_09_04/model_files .
#source $TFNET_ROOT/scripts/run_allele_deeplift.sh
source $TFNET_ROOT/scripts/run_allele_deeplift.sh 0

#python get_scores.py
