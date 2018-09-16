#$TFNET_ROOT/scripts/run_pipeline.py --tfs SPI1 --cells GM12878 --end-task 1 --start 10
$TFNET_ROOT/scripts/run_pipeline.py --tfs SPI1 --cells GM12878 --end-task 1 --start 80 --stride 1  --fdr 0.0001 --min-seqlets=11000
#0.0001 + min 11200  got 11762
#0.0005  got >20k
#0.0001  got  9555
#0.00007 got  6392
