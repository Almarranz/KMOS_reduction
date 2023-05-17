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

period = 'p107'
wave_int = 'B'#Wave length interval used in the reduction
folder = '/Users/alvaromartinez/Desktop/Phd/KMOS/p107_ABC/p107_%s/'%(wave_int)
sci_cubes = '/Users/alvaromartinez/Desktop/Phd/KMOS/p107_ABC/p107_B/end_products/2023-05-15T11:21:27/KMOS.2021-06-06T03:07:34.684_combine_OBs/'
sky = '/Users/alvaromartinez/Desktop/Phd/KMOS/p107_ABC/p107_B_sky/end_products/2023-05-15T13:10:32/KMOS.2021-06-06T03:37:33.041_tpl/GC_K_Sky_SINGLE_CUBES_KMOS.2021-06-06T03:37:33.041.fits'
pruebas = '/Users/alvaromartinez/Desktop/Phd/KMOS/pruebas/'
folder ='/Users/alvaromartinez/Desktop/Phd/KMOS/p107_ABC/p107_B/'
# for file in glob.glob(sci_half1 + '*.fits', recursive=True):
count = 0 
for file in glob.glob(sci_cubes + '*SINGLE_CUBES*', recursive=True):
    count +=1
    with open(folder +'data.sof', 'w') as f:
        f.write('')
        f.close
    with(open(folder +'data.sof','a')) as f:
        f.write(file + ' CUBE_OBJECT\n')
        # f.write(file + ' SINGLE_CUBES\n')
    print('------\n',file)

    with(open(folder +'data.sof','a')) as f:
        f.write(sky+ ' CUBE_SKY\n')
    call(["esorex","kmos_sky_tweak", "data.sof"], cwd =folder)
    # call(["ls"], cwd =pruebas)
    os.rename(folder + 'SKY_TWEAK.fits', folder + 'ob_%s_%s_%s_tweak.fits'%(period, wave_int,count))
# %%
with open(folder +'data_tweak.sof', 'w') as f:
     f.write('')
     f.close
for file in glob.glob(folder + '*tweak.fits', recursive=True):
    count +=1
    with(open(folder +'data_tweak.sof','a')) as f: 
        f.write(file + ' SINGLE_CUBES\n')

call(["esorex","kmos_combine", "--method=header", "data_tweak.sof"], cwd = folder)




























