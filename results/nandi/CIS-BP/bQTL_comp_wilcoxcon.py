
# coding: utf-8

# # bQTL-motif distance vs p-value: TF-MoDISco vs CISBP
# 
# - Calculate the distance between each bQTL SNP to its nearest MoDISco motif and to CIS-BP motif.
# - Plot distance vs -log10(p-value)
# - Cap the distance at 1000bp and plot the distance vs -log10(p-value)

# In[5]:


from __future__ import print_function, division

import sys
import os
#import matplotlib as mpl
#mpl.use('Agg')

import numpy as np
import sys
import argparse


# In[6]:


from collections import OrderedDict
tf_list = [ 'JUND', 'RELA', 'SPI1', 'STAT1']


# In[7]:


import math
import numpy as np
modisco_dir = "/home/ktian/kundajelab/tfnet/results/nandi/"

big_str = "_big"
big_str = ""
big_str = "_full"

modisco_tsv_fns = [
    "JUND/JUND_GM12878_refine_18_09_04/modisco.run2/JUND_modisco_snp_dist_pval" + big_str + ".tsv",
    "RELA/RELA_GM12878_refine_18_09_04/modisco.run3/RELA_modisco_snp_dist_pval" + big_str + ".tsv",
    "SPI1/SPI1_GM12878_refine_18_09_04/modisco.run3/SPI1_modisco_snp_dist_pval" + big_str + ".tsv",
    "STAT1/STAT1_GM12878_refine_18_09_04/modisco.run2/STAT1_modisco_snp_dist_pval" + big_str + ".tsv",
]

import scipy.stats
for i, tf in enumerate(tf_list):


    snp_list_m = []
    with open(modisco_dir + modisco_tsv_fns[i]) as fh:
        for row in fh:
            item = row.split('\t')
            item[0] = float(item[0])
            item[1] = float(item[1])
            if item[0] == 0:
                item[0] = 0.1
            snp_list_m.append(item)

    snp_list_c = []
    with open("under_1ksummit/" + tf + "_cisbp" + big_str + ".tsv") as fh:
        for row in fh:
            item = row.split('\t')
            item[0] = float(item[0])
            item[1] = float(item[1])
            if item[0] == 0:
                item[0] = 0.1
            snp_list_c.append(item)
   
    m50 = [item[1] for item in snp_list_m if item[0] <= 50]
    c50 = [item[1] for item in snp_list_c if item[0] <= 50]

    #print(m50[:5], c50[:5])

    statistic, pval = scipy.stats.ranksums(m50, c50)

    print("%-5s: len < 50: modisco=%4d, cisbp=%4d" % (tf, len(m50), len(c50)), " wilcoxcon= % .6f, pval= % .6f" %(statistic, pval))
    print("     : modisco median=%f, cisbp median = %f\n" % (np.median(m50), np.median(c50)))


