#!/usr/bin/env bash

import os

ROOT_DIR   = os.getenv('TFNET_ROOT', "../../") 
genomeDir  = ROOT_DIR + "/genome/"

os.system("mkdir -p logs")

"""
Input:
-------

Output:
-------
splits/train.txt.gz
splits/valid.txt.gz
splits/test.txt.gz
"""
#os.system("prepare_data.py > logs/prepare.log 2>&1")

#os.system("momma_dragonn_train > logs/train.log 2>&1")

#os.system("mv model_files/record_1_*Json.json model_files/record_1_Json.json")
#os.system("mv model_files/record_1_*Weights.h5 model_files/record_1_Weights.h5")


#os.system("cat labels.txt | tail -n +2 | perl -lane 'if ($.%10==1) {print $F[0]}' | sed 's/:/\t/; s/-/\t/' > splits/subset.tsv")
#os.system("gunzip -c splits/test.txt.gz | sed 's/:/\t/; s/-/\t/' > splits/subset.tsv")

#os.system("bedtools getfasta -fi " + genomeDir + "hg19.fa -bed splits/subset.tsv -fo subset.fa")

os.system("python $TFNET_ROOT/scripts/run_deeplift.py model_files/record_1_ subset.fa 3 > logs/deeplift.log 2>&1")

os.system("cat subset.fa | grep -v '^>' > subset.txt")
os.system("python $TFNET_ROOT/scripts/run_tfmodisco.py ./rescale_conv_revealcancel_fc_multiref_10_task_ subset.txt 3 > logs/deeplift.log 2>&1")


