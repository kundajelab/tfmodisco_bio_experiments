# create snp.txt:
# chrom	pos
# chr5	1295162
# chr5	1295228
# chr5	1295250

python $TFNET_ROOT/scripts/snp_to_bed.py 1000 < snp.txt > interpret.tsv
bedtools getfasta -fi $TFNET_ROOT/genome/hg19.fa -bed interpret.tsv -fo interpret.fa
mkdir scores
python $TFNET_ROOT/scripts/gen_allele.py

mv interpret.fa interpret.fa0
mv interpret.fa* scores
ln -s scores/interpret.fa0 interpret.fa
mkdir logs
source $TFNET_ROOT/scripts/run_allele_deeplift.sh
python /Users/kat/kundajelab/tfnet/scripts/gen_bis_scores.py
