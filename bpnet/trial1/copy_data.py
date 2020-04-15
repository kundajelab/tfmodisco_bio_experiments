import h5py
import numpy as np

f = h5py.File("/oak/stanford/groups/akundaje/avsec/basepair/data/processed/comparison/output/nexus,peaks,OSNK,0,10,1,FALSE,same,0.5,64,25,0.004,9,FALSE,[1,50],TRUE,FALSE,1/deeplift.imp_score.h5", 'r')

outf = h5py.File("trial1.deeplift.imp_score.h5")
outf.create_dataset("hyp_imp/Nanog/profile/wn", data=f["/hyp_imp/Nanog/profile/wn"])
outf.create_dataset("inputs/seq", data=f["inputs/seq"])

f.close()
outf.close()

