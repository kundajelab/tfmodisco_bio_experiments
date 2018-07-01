#!/usr/bin/env bash

import os


os.system("mkdir -f log")

"""
Input:
-------

Output:
-------
splits/train.txt.gz
splits/valid.txt.gz
splits/test.txt.gz
"""
os.system("prepare_data.py >& log/prepare.log")



os.system("momma_dragonn_train >& log/train.log")

os.system("$TFNET_ROOT/script
