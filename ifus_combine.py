#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 10:05:49 2023

@author: alvaromartinez
"""
# Combines ifu from different periods in order to check any improvemt in the SNR.
# The different ifus have been individually alignged with GNS using Aladin.
# The chop of the ifus have been performed using the log at the esorex kmos_combine recipe

import os
import numpy as np
import sys
import glob
from subprocess import call
import astropy.units as u
from astropy.utils.data import download_file
from astropy.io import fits  # We use fits to open the actual data file
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
# %%plotting pa    metres
from matplotlib import rc
from matplotlib import rcParams
rcParams.update({'xtick.major.pad': '7.0'})
rcParams.update({'xtick.major.size': '7.5'})
rcParams.update({'xtick.major.width': '1.5'})
rcParams.update({'xtick.minor.pad': '7.0'})
rcParams.update({'xtick.minor.size': '3.5'})
rcParams.update({'xtick.minor.width': '1.0'})
rcParams.update({'ytick.major.pad': '7.0'})
rcParams.update({'ytick.major.size': '7.5'})
rcParams.update({'ytick.major.width': '1.5'})
rcParams.update({'ytick.minor.pad': '7.0'})
rcParams.update({'ytick.minor.size': '3.5'})
rcParams.update({'ytick.minor.width': '1.0'})
rcParams.update({'font.size': 20})
rcParams.update({'figure.figsize':(30,30)})
rcParams.update({
    "text.usetex": False,
    "font.family": "sans",
    "font.sans-serif": ["Palatino"]})
plt.rcParams["mathtext.fontset"] = 'dejavuserif'
rc('font',**{'family':'serif','serif':['Palatino']})
plt.rcParams.update({'figure.max_open_warning': 0})
# %%


pruebas = '/Users/amartinez/Desktop/PhD/KMOS/practice/'
aling = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/ifu_alignment/'

log_7 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p107_ABC/'
esorex_cube_7 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p107_ABC/COMBINE_SKY_TWEAK_mapping.fits'


log_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_ABC/'
esorex_cube_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_ABC/COMBINE_SKY_TWEAK_mapping.fits'


ifu_sel = 7 #TODO

dic_x5 = {}
dic_y5 = {}
for ifu in range(1,24):
    x_shifts5 = []
    y_shifts5 = []
    with open(log_5 + 'esorex.log', 'r') as f:
        for line in f:
            if 'IFU %s]'%(ifu) in line:
                ind1 = line.find('in x:')
                ind2 = line.find('in y:')
                ind_x = ind1 + 6
                ind_y = ind2 + 6
                x_shifts5.append(float(line[ind_x:ind_x+7]))
                y_shifts5.append(float(line[ind_y:ind_y+7]))
                # print (line)
                # break
            else:
                continue
        x_shifts5 = np.array(x_shifts5)
        y_shifts5 = np.array(y_shifts5)
        if len(x_shifts5)>0:
            dic_x5['ifu%s'%(ifu)] = np.asarray(np.rint(217 - x_shifts5),int)
            dic_y5['ifu%s'%(ifu)] = np.asarray(np.rint(433 + y_shifts5),int)
            # dic_x['ifu%s'%(ifu)] = np.asarray(np.rint(217 - x_shifts),int)
            # dic_y['ifu%s'%(ifu)] = np.asarray(np.rint(433 + y_shifts ),int)

dic_x7 = {}
dic_y7 = {}
for ifu in range(1,24):
    x_shifts7 = []
    y_shifts7 = []
    with open(log_7 + 'esorex.log', 'r') as f:
        for line in f:
            if 'IFU %s]'%(ifu) in line:
                ind1 = line.find('in x:')
                ind2 = line.find('in y:')
                ind_x = ind1 + 6
                ind_y = ind2 + 6
                x_shifts7.append(float(line[ind_x:ind_x+7]))
                y_shifts7.append(float(line[ind_y:ind_y+7]))
                # print (line)
                # break
            else:
                continue
        x_shifts7 = np.array(x_shifts7)
        y_shifts7 = np.array(y_shifts7)
        if len(x_shifts5)>0:
            dic_x7['ifu%s'%(ifu)] = np.asarray(np.rint(217 - x_shifts7),int)
            dic_y7['ifu%s'%(ifu)] = np.asarray(np.rint(433 + y_shifts7),int)
            # dic_x['ifu%s'%(ifu)] = np.asarray(np.rint(217 - x_shifts),int)
            # dic_y['ifu%s'%(ifu)] = np.asarray(np.rint(433 + y_shifts ),int)

# %%
ima_5 = fits.open(esorex_cube_5)
ima5 = ima_5[1].data
noise5 = ima_5[2].data


temp5 = np.zeros((3072,433,650))
# temp5[:] = np.nan
# temp_noise5 =  np.zeros((3072,433,650))
temp5[:,min(dic_y5['ifu%s'%(ifu_sel)])-27:max(dic_y5['ifu%s'%(ifu_sel)]),min(dic_x5['ifu%s'%(ifu_sel)]):max(dic_x5['ifu%s'%(ifu_sel)])+27] = ima5[:,min(dic_y5['ifu%s'%(ifu_sel)])-27:max(dic_y5['ifu%s'%(ifu_sel)]),min(dic_x5['ifu%s'%(ifu_sel)]):max(dic_x5['ifu%s'%(ifu_sel)])+27]
# temp_noise5[:,min(dic_y5['ifu%s'%(ifu_sel)])-27:max(dic_y5['ifu%s'%(ifu_sel)]),min(dic_x5['ifu%s'%(ifu_sel)]):max(dic_x5['ifu%s'%(ifu_sel)])+27] = noise5[:,min(dic_y5['ifu%s'%(ifu_sel)])-27:max(dic_y5['ifu%s'%(ifu_sel)]),min(dic_x5['ifu%s'%(ifu_sel)]):max(dic_x5['ifu%s'%(ifu_sel)])+27]

header1 = ima_5[0].header
header2 = ima_5[1].header
header3 = ima_5[2].header

new_hdul = fits.HDUList()
new_hdul.append(fits.PrimaryHDU(header=header1))
new_hdul.append(fits.ImageHDU(temp5, header=header2,name='DATA'))
# new_hdul.append(fits.ImageHDU(temp_noise5,header=header3, name='ERROR'))
# new_hdul.writeto(pruebas + 'cube_ifu%s_%s.fits'%(ifu_sel,'p105'),overwrite=True)


ima_7 = fits.open(esorex_cube_7)
ima7 = ima_7[1].data
noise7 = ima_7[2].data


temp7 = np.zeros((3072,433,650))
# temp_noise7 =  np.zeros((3072,433,650))
xp, yp = np.loadtxt(aling + 'ifu%s_xy_plus.txt'%(ifu_sel),unpack = True)
xp, yp = int(xp), int(yp)
# yp = -2
# xp = 10
y_d = min(dic_y7['ifu%s'%(ifu_sel)])-27
y_up = max(dic_y7['ifu%s'%(ifu_sel)])
x_d = min(dic_x7['ifu%s'%(ifu_sel)])
x_up = max(dic_x7['ifu%s'%(ifu_sel)])+27

if x_up + xp > ima7.shape[2] or y_up + yp > ima7.shape[1] :
    print('PAAAADDD')
    pad = int(max(x_up + xp - ima7.shape[2],y_up + yp - ima7.shape[1]))
    ima7 = np.pad(ima7,(0,pad), mode = 'constant') 
    ima7 = ima7[:-pad,:,:]
    
temp7[:, y_d : y_up, x_d : x_up] =ima7[:, y_d+yp  : y_up+yp, x_d + xp : x_up +xp ]
# temp7[:,min(dic_y7['ifu%s'%(ifu_sel)])-27:max(dic_y7['ifu%s'%(ifu_sel)]),min(dic_x7['ifu%s'%(ifu_sel)]):max(dic_x7['ifu%s'%(ifu_sel)])+27] = ima7[:,min(dic_y7['ifu%s'%(ifu_sel)])-27+yp:max(dic_y7['ifu%s'%(ifu_sel)])+yp,min(dic_x7['ifu%s'%(ifu_sel)])+xp:max(dic_x7['ifu%s'%(ifu_sel)])+27+ xp]
# temp_noise7[:,min(dic_y7['ifu%s'%(ifu_sel)])-27:max(dic_y7['ifu%s'%(ifu_sel)]),min(dic_x7['ifu%s'%(ifu_sel)]):max(dic_x7['ifu%s'%(ifu_sel)])+27] = noise7[:,min(dic_y7['ifu%s'%(ifu_sel)])-27:max(dic_y7['ifu%s'%(ifu_sel)]),min(dic_x7['ifu%s'%(ifu_sel)]):max(dic_x7['ifu%s'%(ifu_sel)])+27]
# sys.exit(161)
header1 = ima_7[0].header
header2 = ima_7[1].header
header3 = ima_7[2].header


new_hdul = fits.HDUList()
new_hdul.append(fits.PrimaryHDU(header=header1))
new_hdul.append(fits.ImageHDU(temp7, header=header2,name='DATA'))
# new_hdul.append(fits.ImageHDU(temp_noise7,header=header3, name='ERROR'))
# new_hdul.writeto(pruebas + 'cube_ifu%s_%s.fits'%(ifu_sel,'p107'),overwrite=True)
# %%
fig, ax = plt.subplots(1,1,figsize=(8,8))
ax.imshow(temp7[1550,:,:], vmin = -0.8e-20, vmax = 0.1e-16, origin = 'lower', cmap = 'Greys')
# %%
temp = np.mean([temp5,temp7],axis = 0)

new_hdul = fits.HDUList()
new_hdul.append(fits.PrimaryHDU(header=header1))
new_hdul.append(fits.ImageHDU(temp, header=header2,name='DATA'))
new_hdul.writeto(aling + 'cube_ifu%s_mean.fits'%(ifu_sel),overwrite=True)
