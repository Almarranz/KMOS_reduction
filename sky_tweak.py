#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 12:48:45 2023

@author: alvaromartinez
"""
# Makes the necesary files and subtracts the sky using the sky_tweak algorithm

import os
import numpy as np
import sys
import glob
from subprocess import call
import shutil

period = 'p105'
wave_int = input('Which interval? (A,B,C,D or ABC):')
# wave_int = 'B'#Wave length interval used in the reduction
sci_cubes = '/Users/alvaromartinez/Desktop/Phd/KMOS/p105/p105_ABC/end_products/2023-09-29T12:27:36/KMOS.2021-06-03T04:19:43.880_combine_OBs/'
# The sky path has to be all the way to the .fits file you want to use.
sky = '/Users/alvaromartinez/Desktop/Phd/KMOS/p105/p105_ABC_sky/end_products/2023-09-29T16:31:01/KMOS.2021-06-03T04:47:47.641_tpl/GC_K_Sky_SINGLE_CUBES_KMOS.2021-06-03T04:47:47.641.fits'
pruebas = '/Users/alvaromartinez/Desktop/Phd/KMOS/pruebas/'
folder ='/Users/alvaromartinez/Desktop/Phd/KMOS/%s/%s_%s/'%(period, period, wave_int)
ca = os.path.dirname(sci_cubes)
ca_ =os.path.basename(os.path.split(ca)[0])
folder_to_save = folder +ca_ +'/'
os.mkdir(folder_to_save)



# for file in glob.glob(sci_half1 + '*.fits', recursive=True):
count = 0 
for file in glob.glob(sci_cubes + '*SINGLE_CUBES*', recursive=True):
    count +=1
    with open(folder_to_save +'data.sof', 'w') as f:
        f.write('')
        f.close
    
    with(open(folder_to_save +'data.sof','a')) as f:
        f.write(file + ' CUBE_OBJECT\n')
        # f.write(file + ' SINGLE_CUBES\n')
    print('------\n',file)
   
    with(open(folder_to_save +'data.sof','a')) as f:
        f.write(sky+ ' CUBE_SKY\n')
    call(["esorex","kmos_sky_tweak", "data.sof"], cwd =folder_to_save)
    
    # call(["ls"], cwd =pruebas)
    os.rename(folder_to_save+ 'SKY_TWEAK.fits', folder_to_save + 'ob_%s_%s_%s_tweak.fits'%(period, wave_int,count))

# %%
with open(folder_to_save +'data_tweak.sof', 'w') as f:
     f.write('')
     f.close
for file in glob.glob(folder_to_save + '*tweak.fits', recursive=True):
    count +=1
    with(open(folder_to_save +'data_tweak.sof','a')) as f: 
        f.write(file + ' SINGLE_CUBES\n')

call(["esorex","kmos_combine", "--method=header", "data_tweak.sof"], cwd = folder_to_save)





# %%


# del_path = '/Users/alvaromartinez/Desktop/Phd/KMOS/p107_ABC/'



# letts = ['A','B','C','ABC']

# for lett in letts:
#     del_folds =['P107_%s/'%(lett),'P107_%s_sky/'%(lett)]
#     for fold in del_folds:
#         shutil.rmtree(del_path + del_folds[0])
#         shutil.rmtree(del_path + del_folds[1])
#         os.mkdir((del_path + del_folds[0]))
#         os.mkdir((del_path + del_folds[1]))
















