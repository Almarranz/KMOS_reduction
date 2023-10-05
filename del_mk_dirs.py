#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 10:13:36 2023

@author: alvaromartinez
"""

# Delete-creates directories 
import os
import shutil

period = '107'
del_path = '/Users/alvaromartinez/Desktop/Phd/KMOS/p%s/'%(period)


# 
# letts = ['A','B','C','ABC']
letts = ['D']

for lett in letts:
    del_folds =['p%s_%s/'%(period,lett),'p%s_%s_sky/'%(period,lett)]
    for fold in del_folds:
        shutil.rmtree(del_path + del_folds[0])
        shutil.rmtree(del_path + del_folds[1])
        os.mkdir((del_path + del_folds[0]))
        os.mkdir((del_path + del_folds[1]))

kmos_dir = ['end_products', 'log','tmp_products', 'bookkeeping']
for lett in letts:
    del_folds =['p%s_%s/'%(period,lett),'p%s_%s_sky/'%(period,lett)]
    for kmos in kmos_dir:
        # shutil.rmtree(del_path + del_folds[0])
        # shutil.rmtree(del_path + del_folds[1])
        os.mkdir((del_path + del_folds[0] + kmos))
        os.mkdir((del_path + del_folds[1]+ kmos))
    
        

