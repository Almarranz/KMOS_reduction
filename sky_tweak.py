#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 12:48:45 2023

@author: alvaromartinez
"""
# Makes the necesary files and subtracts the sky using the sjy_tweak algorithm

import os
import numpy as np
import sys
import glob
from subprocess import call

period = 'p105'
wave_int = 'C'#Wave length interval used in the reduction
sci_cubes = '/Users/alvaromartinez/Desktop/Phd/KMOS/p105_C/end_products/2023-05-24T12:54:32/KMOS.2021-06-03T04:19:43.880_combine_OBs/'
# The sky path has to be all the way to the .fits file you want to use.
sky = '/Users/alvaromartinez/Desktop/Phd/KMOS/p105_C_sky/end_products/2023-05-24T14:22:33/KMOS.2021-06-03T04:47:47.641_combine_OBs/GC_K_Sky_SINGLE_CUBES_KMOS.2021-06-03T04:47:47.641.fits'
pruebas = '/Users/alvaromartinez/Desktop/Phd/KMOS/pruebas/'
folder ='/Users/alvaromartinez/Desktop/Phd/KMOS/p105_%s/'%(wave_int)
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























