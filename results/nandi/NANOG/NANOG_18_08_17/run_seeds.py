import os
for sd in range(0, 10000, 1000):
    seed = 1234 + sd
    print("seed=%d" % seed)
    os.system("python run_tfmodisco.py scores/hyp_scores_task_ subset_nobg.fa subset_nobg.tsv " + str(seed) + " 1 > logs/modisco.log 2>&1")
    dirname = "modisco.run_seed_" + str(seed)
    os.system("mkdir " + dirname)
    os.system("mkdir " + dirname + "/logs")
    os.system("mv logs/modisco.log " + dirname + "/logs")
    os.system("mv figures " + dirname)
    os.system("mv results.hdf5 " + dirname)
    os.system("cp tfmodisco*.ipynb " + dirname)
