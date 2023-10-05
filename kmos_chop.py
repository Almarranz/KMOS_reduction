#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 13:41:51 2023

@author: alvaromartinez
"""

# Chops the KMOS images with the pourpose of a better alignment with GNS
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
plt.rcParams.update({'figure.max_open_warning': 0})# 
# %%
period = 'p107'
pruebas = '/Users/alvaromartinez/Desktop/Phd/KMOS/pruebas/'

if period == 'p107':
    log = '/Users/alvaromartinez/Desktop/Phd/KMOS/p107/p107_ABC/log/kmos_combine_1/2023-10-02T19:13:38.496/'
    esorex_ima = '/Users/alvaromartinez/Desktop/Phd/KMOS/p107/p107_ABC/end_products/2023-10-02T18:20:11/KMOS.2021-06-06T03:07:34.684_combine_OBs/comoving_group_mosaic_K_Half1_COMBINED_IMAGE_mapping.fits'

elif period == 'p105':
    log = '/Users/alvaromartinez/Desktop/Phd/KMOS/p105/p105_ABC/log/kmos_combine_1/2023-09-29T15:10:04.845/'
    esorex_ima = '/Users/alvaromartinez/Desktop/Phd/KMOS/p105/p105_ABC/end_products/2023-09-29T12:27:36/KMOS.2021-06-03T04:19:43.880_combine_OBs/comoving_group_mosaic_K_Half1_COMBINED_IMAGE_mapping.fits'

# ima = fits.getdata(pruebas + 'one_ABC.fits')
# ima = np.squeeze(ima)


# ima = fits.getdata(esorex_ima)
ima_ = fits.open(esorex_ima)
ima = ima_[1].data
noise = ima_[2].data

fig, ax = plt.subplots(1,1,figsize=(8,8))
ax.imshow(ima, vmin = -0.8e-20, vmax = 0.1e-16, origin = 'lower', cmap = 'Greys')
ax.axvline(650-27*4, color = 'r', ls = 'dashed', alpha = 0.4)
ax.set_xlim(650-28*4,650)
ax.set_ylim(433-27*4,433)

# %%

dic_x = {}
dic_y = {}
for ifu in range(1,24):
    x_shifts = []
    y_shifts = []
    with open(log + 'esorex.log', 'r') as f:
        for line in f:
            if 'IFU %s]'%(ifu) in line:
                ind1 = line.find('in x:')
                ind2 = line.find('in y:')
                ind_x = ind1 + 6
                ind_y = ind2 + 6
                x_shifts.append(float(line[ind_x:ind_x+7]))
                y_shifts.append(float(line[ind_y:ind_y+7]))
                # print (line)
                # break
            else:
                continue
        x_shifts = np.array(x_shifts)
        y_shifts = np.array(y_shifts)
        if len(x_shifts)>0:
            dic_x['ifu%s'%(ifu)] = np.asarray(np.rint(217 - x_shifts),int)
            dic_y['ifu%s'%(ifu)] = np.asarray(np.rint(433 + y_shifts),int)
            # dic_x['ifu%s'%(ifu)] = np.asarray(np.rint(217 - x_shifts),int)
            # dic_y['ifu%s'%(ifu)] = np.asarray(np.rint(433 + y_shifts ),int)

# header1 = ima_[0].header
# header2 = ima_[1].header
# header3 = ima_[2].header

# ifu_sel = 7

# temp = np.zeros((433,650))
# temp_noise =  np.zeros((433,650))
# yp = -2
# xp = -1
# temp[min(dic_y['ifu%s'%(ifu_sel)])-27:max(dic_y['ifu%s'%(ifu_sel)]), min(dic_x['ifu%s'%(ifu_sel)]):max(dic_x['ifu%s'%(ifu_sel)])+27] = ima[min(dic_y['ifu%s'%(ifu_sel)])-27+yp:max(dic_y['ifu%s'%(ifu_sel)])+yp, min(dic_x['ifu%s'%(ifu_sel)])-xp: max(dic_x['ifu%s'%(ifu_sel)])+27-xp]
# temp_noise[min(dic_y['ifu%s'%(ifu_sel)])-27:max(dic_y['ifu%s'%(ifu_sel)]),min(dic_x['ifu%s'%(ifu_sel)]):max(dic_x['ifu%s'%(ifu_sel)])+27] = noise[min(dic_y['ifu%s'%(ifu_sel)])-27:max(dic_y['ifu%s'%(ifu_sel)]),min(dic_x['ifu%s'%(ifu_sel)]):max(dic_x['ifu%s'%(ifu_sel)])+27]


# new_hdul = fits.HDUList()
# new_hdul.append(fits.PrimaryHDU(header=header1))
# new_hdul.append(fits.ImageHDU(temp, header=header2,name='DATA'))
# new_hdul.append(fits.ImageHDU(temp_noise,header=header3, name='ERROR'))
# new_hdul.writeto(pruebas + 'ifu%s_%s_not_aligned.fits'%(ifu_sel,period),overwrite=True)
sys.exit('113')
# %%
ifu_sel = 7

temp = np.zeros((433,650))
temp_noise =  np.zeros((433,650))
temp[min(dic_y['ifu%s'%(ifu_sel)])-27:max(dic_y['ifu%s'%(ifu_sel)]),min(dic_x['ifu%s'%(ifu_sel)]):max(dic_x['ifu%s'%(ifu_sel)])+27] = ima[min(dic_y['ifu%s'%(ifu_sel)])-27:max(dic_y['ifu%s'%(ifu_sel)]),min(dic_x['ifu%s'%(ifu_sel)]):max(dic_x['ifu%s'%(ifu_sel)])+27]
temp_noise[min(dic_y['ifu%s'%(ifu_sel)])-27:max(dic_y['ifu%s'%(ifu_sel)]),min(dic_x['ifu%s'%(ifu_sel)]):max(dic_x['ifu%s'%(ifu_sel)])+27] = noise[min(dic_y['ifu%s'%(ifu_sel)])-27:max(dic_y['ifu%s'%(ifu_sel)]),min(dic_x['ifu%s'%(ifu_sel)]):max(dic_x['ifu%s'%(ifu_sel)])+27]

header1 = ima_[0].header
header2 = ima_[1].header
header3 = ima_[2].header

paras = ['NAXIS1',	'NAXIS2',	
         'CRPIX1',	'CRPIX2',	
        'CRVAL1',	'CRVAL2',	
        'CD1_1',	'CD1_2',	
        'CD2_1',	'CD2_2']
# if period = 'p105':
    # NAXIS1  = 650
    # NAXIS2  = 433
    # CRPIX1  = 313.73946895425433
    # CRPIX2  = 323.77945829814655
    # EQUINOX = 2000.0
    # CRVAL1  = 266.384440107102
    # CRVAL2  = -28.9444025
    # CTYPE1  = RA---TAN
    # CTYPE2  = DEC--TAN
    # RADECSYS= ICRS
    # CD1_1   = 7.3124987844013485E-6
    # CD1_2   = 2.6520556700209487E-5
    # CD2_1   = 2.6520556700209487E-5
    # CD2_2   = -7.3124987844013485E-6
    
# if period = 'p107':   
    # NAXIS1  = 650
    # NAXIS2  = 433
    # CRPIX1  = 319.6167163349949
    # CRPIX2  = 321.2200399618404
    # EQUINOX = 2000.0
    # CRVAL1  = 266.384440107102
    # CRVAL2  = -28.9444025
    # CTYPE1  = RA---TAN
    # CTYPE2  = DEC--TAN
    # RADECSYS= ICRS
    # CD1_1   = 7.132068353100868E-6
    # CD1_2   = 2.7043877610461897E-5
    # CD2_1   = 2.7043877610461897E-5
    # CD2_2   = -7.132068353100868E-6


new_para =[650,         433,
           319.6167163349949,  321.2200399618404,
           266.384440107102, -28.9444025,
           7.132068353100868E-6,   2.7043877610461897E-5
           ,2.7043877610461897E-5, -7.132068353100868E-6]
for i, para in enumerate(paras):
    print(header2[para])
    header2[para] = new_para[i]
    header3[para] = new_para[i]

new_hdul = fits.HDUList()
new_hdul.append(fits.PrimaryHDU(header=header1))
new_hdul.append(fits.ImageHDU(temp, header=header2,name='DATA'))
new_hdul.append(fits.ImageHDU(temp_noise,header=header3, name='ERROR'))
new_hdul.writeto(pruebas + 'ifu%s_%s.fits'%(ifu_sel,period),overwrite=True)


fig, ax = plt.subplots(1,1,figsize=(8,8))
ax.imshow(temp, vmin = -0.8e-20, vmax = 0.1e-16, origin = 'lower', cmap = 'Greys')










































