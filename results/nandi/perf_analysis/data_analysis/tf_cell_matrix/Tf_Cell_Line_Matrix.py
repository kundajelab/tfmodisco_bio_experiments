
# coding: utf-8

# # Data for multi-cell-line Tf-binding models

# In[109]:


import os

dataDir = "/Users/kat/kundajelab/tfnet/ENCODE_data/"
namesFile = "/home/ktian/kundajelab/tfnet/metadata/filenames.txt"
os.chdir(dataDir)
os.system('pwd')

# generate a list of peak file names
os.system(u'ls > fileNames')


# In[112]:


tf_to_cell = {}
tf_set = set()
cell_set = set()

file_list = []
for fname in open(namesFile):
    fname = fname.rstrip()
    fname = fname.split('-')
    n = len(fname)
    peak_type = fname[-1]
    this_experiment = fname[-2]
    this_tf = fname[-4]
    # note: some use eGFP, if len(fname) != 5 and fname[1]=='eGFP'
    this_cell = '-'.join(fname[:-4])
    
    file_list.append([this_tf, this_cell, this_experiment])
    
    if peak_type == "optimal_idr.narrowPeak.gz":
        tf_set.add(this_tf)
        cell_set.add(this_cell)
        
        if this_tf in tf_to_cell.keys():
            tf_to_cell[this_tf].append(this_cell)

        else:
            tf_to_cell[this_tf]=[this_cell]


# In[114]:


tf_set = list(tf_set)
cell_set = sorted(list(cell_set))

tf_count = len(tf_set)
cell_count = len(cell_set)

print(tf_count)
print(cell_count)


# In[118]:


import pandas as pd
import numpy as np

df = pd.DataFrame( data=np.zeros((tf_count,cell_count)), dtype=int,
                   index = tf_set, columns = cell_set)

for tf in tf_to_cell:
    for cell in tf_to_cell[tf]:
        df.loc[tf,cell] = 1



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


df3.to_csv('tf_cell1.csv')

