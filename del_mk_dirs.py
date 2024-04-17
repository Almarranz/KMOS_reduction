#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 10:13:36 2023

@author: alvaromartinez
"""

# Delete-creates directories 
import os
import shutil
import numpy as np

period = '107'
# del_path = '/Users/alvaromartinez/Desktop/Phd/KMOS/p%s/'%(period)
del_path = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/ABC_reduction/young_candidates/'

# 
# letts = ['A','B','C','ABC']
# letts = ['D']
halfs = [1,2]
ifus = np.arange(1,25)


# for lett in letts:
#     del_folds =['p%s_%s/'%(period,lett),'p%s_%s_sky/'%(period,lett)]
#     for fold in del_folds:
#         shutil.rmtree(del_path + del_folds[0])
#         shutil.rmtree(del_path + del_folds[1])
#         os.mkdir((del_path + del_folds[0]))
#         os.mkdir((del_path + del_folds[1]))

# kmos_dir = ['end_products', 'log','tmp_products', 'bookkeeping']
# for lett in letts:
#     del_folds =['p%s_%s/'%(period,lett),'p%s_%s_sky/'%(period,lett)]
#     for kmos in kmos_dir:
#         # shutil.rmtree(del_path + del_folds[0])
#         # shutil.rmtree(del_path + del_folds[1])
#         os.mkdir((del_path + del_folds[0] + kmos))
#         os.mkdir((del_path + del_folds[1]+ kmos))

for ifu in ifus:
    os.mkdir((del_path +'ifu_%s'%(ifu)))
    for half in halfs:
        # os.mkdir((del_path +'ifu_%s/half_%s'%(ifu, half)))
        
        os.mkdir((del_path +'ifu_%s'%(ifu) + '/half_%s'%(half)))
    
        

# %%
from_file = '/Users/amartinez/Desktop/PhD/My_papers/Libralato/'
to_file = '/Users/amartinez/Desktop/PhD/My_papers/Libralato/revised_files_to_editor/'

name = 'ms_id_139573_knn25_rad50.0_mean.png'

shutil.copyfile(from_file + name, to_file + name)








