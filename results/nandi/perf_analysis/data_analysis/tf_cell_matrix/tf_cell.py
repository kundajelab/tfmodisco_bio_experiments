#!/usr/bin/env python


import logging
'''
logging.basicConfig(
        format='%(asctime)s %(levelname)-5s %(message)s',
        level=logging.DEBUG,
        datefmt='%Y-%m-%d %H:%M:%S')
'''


import sys
logging.info(" ".join(sys.argv))

import os


"""
TFNET_ROOT is set to the root of TFNET dir
 +-- scripts       such as label_region, prepare_data.py, run_deeplift.py, run_modisco.py, run_pipeline.py etc
 +-- genome        hg19.fa hg19.chrom.sizes
 +-- ENCODE_data   intervals
 +-- results       results of the pipeline
      +-- ZNF143
      +-- CTCF
      +-- SIX5
"""
ROOT_DIR   = os.getenv('TFNET_ROOT', "../../") 
scriptDir  = ROOT_DIR + "/scripts/"
dataDir    = ROOT_DIR + "/ENCODE_data/"
genomeDir  = ROOT_DIR + "/genome/"
resultsDir = "./"
logDir     = resultsDir + "log/"
#tmpDir     = "./tmp/"

import sys

positives = []
ambiguous = []

# must run from resultsDir

# loop through the TFs

import glob
#tf_files = glob.glob(dataDir + "*-" + tf + "-human-*-optimal*")
all_files = glob.glob(dataDir + "*optimal*")

count = 0
task_list = []
for path_name in all_files:
    fn = os.path.basename(path_name)
    fn_list = fn.split('-')
    exp  = fn_list[-2]
    tf   = fn_list[-4]
    cell = '-'.join(fn_list[:-4])
    task_list.append([cell, tf, exp])
    #print(path_name)
    count = count + 1

#print(task_list)
print("total experiments=", count)

from sets import Set

cell_set = Set([])
tf_set   = Set([])
for cell, tf, exp in task_list:
    cell_set.add(cell)
    tf_set.add(tf)

tf_count = len(tf_set)
cell_count = len(cell_set)

print("total # of tfs   = ", tf_count)
print("total # of cells = ", cell_count)

import pandas as pd
import numpy as np

I = pd.Index(tf_set) #, name="rows")
C = pd.Index(cell_set) #, name="cols")
df = pd.DataFrame(data=np.zeros((tf_count, cell_count),dtype=int),
                  index=I, columns=C)

for cell, tf, exp in task_list:
    df.loc[tf, cell] += 1

#dfs = df.iloc[:5, :5]
dfs = df

dfs.loc['col_sum']= dfs.sum()

#print(dfs)
df1 = dfs.sort_values(by='col_sum', axis=1, ascending=False)
print("sorted columns")

df1.loc[:,'row_sum']= df1.sum(axis=1)
print(df1)

df2=df1.sort_values(by='row_sum', axis=0, ascending=False)

df3=df2.iloc[1:] # remove col_sum row
df3=df3.append(df2.loc['col_sum'])

print(df3)


df3.to_csv('tf_cell.csv')
'''
#trying to sort
print(df.loc['GM12878','ZNF143'])
df.reindex_axis(df.mean().sort_values(ascending=False ).index, axis=1)
df.reindex_axis(df.mean().sort_values(ascending=False).index, axis=0)
print(df)
'''


