#!/usr/bin/env bash
# to run the modisco pipeline

import os

ROOT_DIR   = os.getenv('TFNET_ROOT', "../../") 
genomeDir  = ROOT_DIR + "/genome/"

os.system("mkdir -p logs")

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
os.system("python prepare_data.py > logs/prepare.log 2>&1")

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
os.system("momma_dragonn_train > logs/train.log 2>&1")

# re-name momma dragonn model files
os.system("mv model_files/record_1_*Json.json model_files/record_1_Json.json")
os.system("mv model_files/record_1_*Weights.h5 model_files/record_1_Weights.h5")

'''
Outputs of momma_dragonn_train (after rename)
------
model_files/record_1_Json.json  # model architecture
model_files/record_1_Weights.h5 # model weights

Additional input for deeplift
------
subset.fa      # fasta file for the sequences. either 10% of inputs.fa, or test set (chr1)

'''

# use 1 mod 10 as subset on which to run deeplift
#os.system("cat labels.txt | tail -n +2 | perl -lane 'if ($.%10==1) {print $F[0]}' | sed 's/:/\t/; s/-/\t/' > splits/subset.tsv")

# use the test set as the subset for deeplift
os.system("gunzip -c splits/test.txt.gz | sed 's/:/\t/; s/-/\t/' > splits/subset.tsv")
os.system("bedtools getfasta -fi " + genomeDir + "hg19.fa -bed splits/subset.tsv -fo subset.fa")

'''
Step 3: run_deeplift.py
-----------------------------------------------------------------------
'''
# run deeplift
os.system("python $TFNET_ROOT/scripts/run_deeplift.py model_files/record_1_ subset.fa 3 > logs/deeplift.log 2>&1")

#os.system("cat subset.fa | grep -v '^>' > subset.txt")  # modisco will take fasta files

'''
Output from deeplift
------
importance scores
./scores/hyp_scores_task_0.npy
...
./scores/hyp_scores_task_N.npy
where N= number of tasks

Additional input for modisco:
------
subset.txt    # sequences corresponding to subset.fa, excluding the fasta header lines with '>'

Step 4: run_tfmodisco.py
-----------------------------------------------------------------------
'''
# run tf modisco
os.system("python $TFNET_ROOT/scripts/run_tfmodisco.py ./scores/hyp_scores_task_ subset.fa 3 > logs/modisco.log 2>&1")


