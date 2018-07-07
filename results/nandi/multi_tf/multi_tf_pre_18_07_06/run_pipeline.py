#!/usr/bin/env bash
# to run the modisco pipeline

import os

ROOT_DIR    = os.getenv('TFNET_ROOT', "../../") 
genomeDir   = ROOT_DIR + "/genome/"
resultDir   = ROOT_DIR + "/results/"
templateDir = resultDir + "/templates/"

num_task = 3

import sys
if len(sys.argv) > 2:
    print("Syntax: ", sys.argv[0] , " [step]")
    quit()

if len(sys.argv) == 2:
    step = int(sys.argv[1])
else:
    step = 1 # start with normal training

os.system("mkdir -p logs")

#0 prepare_data with union of positives (no background) and train a model
if step == 0:
    """
    os.system("cp -r " + templateDir + "/config .")
    os.system("ln -s " + templateDir + "/make_hdf5_yaml .")

    os.system("python prepare_data.py --no-bg > logs/pre_prepare.log 2>&1")
    os.system("cp -f config/hyperparameter_configs_list.yaml.from_scratch config/hyperparameter_configs_list.yaml")
    os.system("momma_dragonn_train > logs/pre_train.log 2>&1")
    """

    os.chdir("model_files")
    os.system("ln -s record_1_*Json.json  record_1_Json.json")
    os.system("ln -s record_1_*Weights.h5 record_1_Weights.h5")
    os.chdir("..")
    os.system("mv model_files model_files_pretrain")

    os.system("gunzip -c splits/test.txt.gz | sed 's/:/\t/; s/-/\t/' | sort -k1,1 -k2,2n > splits/subset_nobg.tsv")
    os.system("bedtools getfasta -fi " + genomeDir + "hg19.fa -bed splits/subset_nobg.tsv -fo subset_nobg.fa")
    os.system("mv splits splits_pretrain")

    os.system("mv runs_perf-metric-auROC.db model_files_pretrain/")
    os.system("cp -f config/hyperparameter_configs_list.yaml.pretrain config/hyperparameter_configs_list.yaml")
    print("step 0 pre_train done")

#1 prepare_data with background for the main training
if step >= 1:
    os.system("python prepare_data.py > logs/prepare.log 2>&1")
    print("step 1 prepare_data done")

#2 train to continue from pre-trained data
if step >= 2:
    os.system("momma_dragonn_train > logs/train.log 2>&1")
    os.chdir("model_files")
    os.system("ln -s record_1_*Json.json  record_1_Json.json")
    os.system("ln -s record_1_*Weights.h5 record_1_Weights.h5")
    os.chdir("..")
    print("step 2 training done")

#3 deeplift
if step >= 3:
    # use 1 mod 10 as subset on which to run deeplift
    # os.system("cat labels.txt | tail -n +2 | perl -lane 'if ($.%10==1) {print $F[0]}' | sed 's/:/\t/; s/-/\t/' > splits/subset.tsv")

    # use the test set as the subset for deeplift
    # os.system("gunzip -c splits/test.txt.gz | perl -lane 'if ($.%2==1) {print}' | sed 's/:/\t/; s/-/\t/' | sort -k1,1 -k2,2n > splits/subset.tsv") # select half of testset
    #os.system("gunzip -c splits/test.txt.gz | sed 's/:/\t/; s/-/\t/' | sort -k1,1 -k2,2n > splits/subset.tsv")
    #os.system("bedtools getfasta -fi " + genomeDir + "hg19.fa -bed splits/subset.tsv -fo subset.fa")

    os.system("python $TFNET_ROOT/scripts/run_deeplift.py model_files/record_1_ subset_nobg.fa 3 > logs/deeplift.log 2>&1")
    print("step 3 deeplift done")

#4 modisco
if step >= 4:
    os.system("python $TFNET_ROOT/scripts/run_tfmodisco.py scores/hyp_scores_task_ subset_nobg.txt 3 > logs/modisco.log 2>&1")

"""
'''
Input for prepare_data.py
-------
list of positive  peaks for each task
list of ambiguous peaks for each task
fasta of the whole genome
sizes of the chromosomes

Step 1: prepare_data.py
-----------------------------------------------------------------------
'''
# run prepare_data.py (executable), which calls label_regions script to create labeled bins

'''
Outputs of prepare_data.py, inputs for momma_dragonn_train
-------
inputs.fa           # the input sequences for momma_dragonn_train
splits/train.txt.gz # training set
splits/valid.txt.gz # validation set
splits/test.txt.gz  # testing set


Step 2: momma_dragonn_train
-----------------------------------------------------------------------
'''
# train momma dragonn model

# re-name momma dragonn model files

'''
Outputs of momma_dragonn_train (after rename)
------
model_files/record_1_Json.json  # model architecture
model_files/record_1_Weights.h5 # model weights

Additional input for deeplift
------
subset.fa      # fasta file for the sequences. either 10% of inputs.fa, or test set (chr1)

'''

'''
Step 3: run_deeplift.py
-----------------------------------------------------------------------
'''
# run deeplift

'''
Output from deeplift
------
importance scores
scores/hyp_scores_task_0.npy
...
scores/hyp_scores_task_N.npy
where N= number of tasks

Additional input for modisco:
------
subset.txt    # sequences corresponding to subset.fa, excluding the fasta header lines with '>'

Step 4: run_tfmodisco.py
-----------------------------------------------------------------------
'''
# run tf modisco

"""

