cat $TFNET_ROOT/results/nandi/bQTL/analysis/bQTL_all_SNPs/SPI1_10k.txt | python $TFNET_ROOT/scripts/snp_to_bed.py 1000 > SPI1_10k.tsv
ln -s SPI1_10k.tsv interpret.tsv
bedtools getfasta -fi $TFNET_ROOT/genome/hg19.fa -bed interpret.tsv -fo interpret.fa
#link to model_files 
ln -s /home/ktian/kundajelab/tfnet/results/nandi/SPI1/SPI1_GM12878_refine_18_09_04/model_files .
python $TFNET_ROOT/scripts/run_pipeline.py --start 70 --end 80 --end-task 1
python get_scores.py
