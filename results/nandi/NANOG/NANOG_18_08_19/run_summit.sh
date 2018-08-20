#python /home/ktian/kundajelab/tfnet/scripts/pick_summit.py NANOG H1-hESC | bedtools sort > NANOG_summit.tsv
#bedtools getfasta -fi /home/ktian/kundajelab/tfnet/genome/hg19.fa -bed NANOG_summit.tsv -fo NANOG_summit.fa
#python /home/ktian/kundajelab/tfnet/scripts/run_deeplift.py model_files/record_1_ NANOG_summit.fa 1 &> logs/deeplift.log
python run_tfmodisco.py scores/hyp_scores_task_ NANOG_summit.fa NANOG_summit.tsv 1234 1 > logs/modisco.log 2>&1
