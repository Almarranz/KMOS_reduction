#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 13:41:51 2023

@author: alvaromartinez
"""

# CChop the selected ifus from each epoch. You have to open them in DS9 and 
# check the pixel offset and correct it by hand, by modifing 
# the parametres xp and yp
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
from astropy.wcs import WCS
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
period = 'p105'#TODO
reduction = 'ABC'

pruebas = '/Users/amartinez/Desktop/PhD/KMOS/practice/'
aling = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/ifu_alignment_%s/'%(reduction)
if period == 'p107':
    log = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p107_%s/'%(reduction)
    esorex_ima = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p107_%s/comoving_group_mosaic_K_Half1_COMBINED_IMAGE_mapping.fits'%(reduction)

elif period == 'p105':
    log = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/'%(reduction)
    esorex_ima = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/comoving_group_mosaic_K_Half1_COMBINED_IMAGE_mapping.fits'%(reduction)






ima_ = fits.open(esorex_ima)
ima = ima_[1].data
noise = ima_[2].data



header1 = ima_[0].header
header2 = ima_[1].header
header3 = ima_[2].header

paras = ['NAXIS1',	'NAXIS2',	
          'CRPIX1',	'CRPIX2',	
        'CRVAL1',	'CRVAL2',	
        'CD1_1',	'CD1_2',	
        'CD2_1',	'CD2_2']
if period =='p105':
    NAXIS1  = 650
    NAXIS2  = 433
    CRPIX1  = 313.73946895425433
    CRPIX2  = 323.77945829814655
    EQUINOX = 2000.0
    CRVAL1  = 266.384440107102
    CRVAL2  = -28.9444025

    CD1_1   = 7.3124987844013485E-6
    CD1_2   = 2.6520556700209487E-5
    CD2_1   = 2.6520556700209487E-5
    CD2_2   = -7.3124987844013485E-6
    
if period == 'p107':   
    NAXIS1  = 650
    NAXIS2  = 433
    CRPIX1  = 319.6167163349949
    CRPIX2  = 321.2200399618404
    EQUINOX = 2000.0
    CRVAL1  = 266.384440107102
    CRVAL2  = -28.9444025
    CD1_1   = 7.132068353100868E-6
    CD1_2   = 2.7043877610461897E-5
    CD2_1   = 2.7043877610461897E-5
    CD2_2   = -7.132068353100868E-6


new_para =[650,         433,
            CRPIX1,  CRPIX2,
            CRVAL1, CRVAL2,
            CD1_1,   CD1_2
            ,CD2_1, CD2_2]
for i, para in enumerate(paras):
    print(header2[para])
    header2[para] = new_para[i]
    header2[para] = new_para[i]

wcs = WCS(header2)
hdu_alig = fits.HDUList()
hdu_alig.append(fits.PrimaryHDU(header=header2))
hdu_alig.append(fits.ImageHDU(ima, header=header2,name='DATA'))
# new_hdul.append(fits.ImageHDU(temp_noise,header=header3, name='ERROR'))
hdu_alig.writeto(pruebas + 'ima_%s_aligned.fits'%(period),overwrite=True)

# sys.exit(127)

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


# sys.exit('125')
# %%

ifus =  [ 3,  5,  6,  7,  8,  9, 10, 11, 12, 13,15, 17,
         19, 21, 23]
halfs = [1,2]
# ifu_sel = 6#TODO
# half_ifu = 2#TODO
for ifu_sel in ifus:
    for half_ifu in halfs:
        yp =0#TODO Negative move ifus upward
        xp = 0#TODO Negative move ifus to the right
        # temp = np.zeros((433,650))
        # temp_noise =  np.zeros((433,650))
        if half_ifu == 0:
            y_d = int(min(dic_y['ifu%s'%(ifu_sel)])-27)
            y_up = int(max(dic_y['ifu%s'%(ifu_sel)]))
            x_d = int(min(dic_x['ifu%s'%(ifu_sel)]))
            x_up = int(max(dic_x['ifu%s'%(ifu_sel)])+27)
        elif half_ifu == 1:
            y_d = int(min(dic_y['ifu%s'%(ifu_sel)])-27)
            y_up = np.unique(dic_y['ifu%s'%(ifu_sel)])[1]
            x_d = int(min(dic_x['ifu%s'%(ifu_sel)]))
            x_up = int(max(dic_x['ifu%s'%(ifu_sel)])+27)
        elif half_ifu == 2:   
            y_d = np.unique(dic_y['ifu%s'%(ifu_sel)])[1]
            y_up = np.unique(dic_y['ifu%s'%(ifu_sel)])[3]
            x_d = int(min(dic_x['ifu%s'%(ifu_sel)]))
            x_up = int(max(dic_x['ifu%s'%(ifu_sel)])+27)
        temp = np.zeros((y_up-y_d, x_up-x_d))
        temp_noise =  np.zeros((y_up-y_d, x_up-x_d))
        if period == 'p105':
            temp = ima[y_d:y_up, x_d:x_up]
            # temp[217 : 325, 541 : 649] = ima[y_d:y_up, x_d:x_up]
            temp_noise = noise[y_d:y_up, x_d:x_up]
            ifu_header = wcs[y_d:y_up, x_d:x_up]
        elif period == 'p107':
            # print('temp dim = %s,%s'%(min(dic_y['ifu%s'%(ifu_sel)])-27 - max(dic_y['ifu%s'%(ifu_sel)]), min(dic_x['ifu%s'%(ifu_sel)])-(max(dic_x['ifu%s'%(ifu_sel)])+27)))
            # print('ima dim = %s, %s'%(min(dic_y['ifu%s'%(ifu_sel)])-27+yp -  (max(dic_y['ifu%s'%(ifu_sel)])+yp) ,min(dic_x['ifu%s'%(ifu_sel)])+xp-(max(dic_x['ifu%s'%(ifu_sel)])+27+xp)))
          
            pad = 0
            # sys.exit()
            if x_up + xp > ima.shape[1] or y_up + yp > ima.shape[0]:
                print('PAD X left or PAD Y up')
                pad = max(x_up + xp - ima.shape[1],y_up + yp - ima.shape[0])
                # ima = np.pad(ima,(0,pad), mode = 'constant') 
                ima = np.pad(ima,(0,pad), mode = 'minimum') 
                
            temp = ima[y_d+yp  : y_up+yp, x_d + xp : x_up +xp ]
            ifu_header = wcs[y_d+yp  : y_up+yp, x_d + xp : x_up +xp]    
                
            if y_d+yp < 0:
                print('PAD Y down')
                ima = np.pad(ima,(abs(y_d+yp),0),mode = 'minimum')
                temp = ima[0  : 54 , x_d + xp : x_up +xp ]
                ifu_header = wcs[0  : 54, x_d + xp : x_up +xp]
                
            
            
            # temp_noise[min(dic_y['ifu%s'%(ifu_sel)])-27:max(dic_y['ifu%s'%(ifu_sel)]),min(dic_x['ifu%s'%(ifu_sel)]):max(dic_x['ifu%s'%(ifu_sel)])+27] = noise[min(dic_y['ifu%s'%(ifu_sel)])-27+yp:max(dic_y['ifu%s'%(ifu_sel)])+yp,min(dic_x['ifu%s'%(ifu_sel)])+xp:max(dic_x['ifu%s'%(ifu_sel)])+27+xp]
            
            save_offsets = input('Do you want to save this offsets?(y/n)')
            # save_offsets = 'n'
            if save_offsets == 'y':
                if half_ifu == 0:
                    np.savetxt(aling + 'ifu%s_xy_plus.txt'%(ifu_sel),np.array([[xp,yp]]),fmt = '%.0f',header = 'xp, yp') 
                elif half_ifu != 0:
                    np.savetxt(aling + 'ifu%s_half%s_xy_plus.txt'%(ifu_sel,half_ifu),np.array([[xp,yp]]),fmt = '%.0f',header = 'xp, yp') 
        
        
        new_hdul = fits.HDUList()
        new_hdul.append(fits.PrimaryHDU(header=header2))
        new_hdul.append(fits.ImageHDU(temp, header=header2,name='DATA'))
        new_hdul[1].header.update(ifu_header.to_header())
        # new_hdul.append(fits.ImageHDU(temp_noise,header=header3, name='ERROR'))
        if half_ifu == 0:
            new_hdul.writeto(aling + 'ifu%s_%s.fits'%(ifu_sel,period),overwrite=True)
        elif half_ifu !=0:
            new_hdul.writeto(aling + 'ifu%s_half%s_%s.fits'%(ifu_sel,half_ifu,period),overwrite=True)
            
        
        
        fig, ax = plt.subplots(1,1,figsize=(8,8))
        ax.imshow(temp, vmin = -0.8e-20, vmax = 0.1e-16, origin = 'lower', cmap = 'Greys')
# %%

# This part roughly alings the image (not the cube) with GNS. 
# =============================================================================
# period = 'p105'#TODO
# pruebas = '/Users/amartinez/Desktop/PhD/KMOS/practice/'
# aling = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/ifu_alignment/'
# if period == 'p107':
#     log = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p107_ABC/'
#     esorex_ima = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p107_ABC/comoving_group_mosaic_K_Half1_COMBINED_IMAGE_mapping.fits'
# 
# elif period == 'p105':
#     log = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_ABC/'
#     esorex_ima = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_ABC/comoving_group_mosaic_K_Half1_COMBINED_IMAGE_mapping.fits'
# 
# imagen = fits.open(esorex_ima)
# cab0 = imagen[0].header
# cab1 = imagen[1].header
# cab2 = imagen[2].header
# paras = ['NAXIS1',	'NAXIS2',	
#           'CRPIX1',	'CRPIX2',	
#         'CRVAL1',	'CRVAL2',	
#         'CD1_1',	'CD1_2',	
#         'CD2_1',	'CD2_2']
# if period =='p105':
#     NAXIS1  = 650
#     NAXIS2  = 433
#     CRPIX1  = 313.73946895425433
#     CRPIX2  = 323.77945829814655
#     EQUINOX = 2000.0
#     CRVAL1  = 266.384440107102
#     CRVAL2  = -28.9444025
# 
#     CD1_1   = 7.3124987844013485E-6
#     CD1_2   = 2.6520556700209487E-5
#     CD2_1   = 2.6520556700209487E-5
#     CD2_2   = -7.3124987844013485E-6
#     
# if period == 'p107':   
#     NAXIS1  = 650
#     NAXIS2  = 433
#     CRPIX1  = 319.6167163349949
#     CRPIX2  = 321.2200399618404
#     EQUINOX = 2000.0
#     CRVAL1  = 266.384440107102
#     CRVAL2  = -28.9444025
#     CD1_1   = 7.132068353100868E-6
#     CD1_2   = 2.7043877610461897E-5
#     CD2_1   = 2.7043877610461897E-5
#     CD2_2   = -7.132068353100868E-6
# 
# 
# new_para =[650,         433,
#             CRPIX1,  CRPIX2,
#             CRVAL1, CRVAL2,
#             CD1_1,   CD1_2
#             ,CD2_1, CD2_2]
# for i, para in enumerate(paras):
#     print(cab1[para])
#     cab1[para] = new_para[i]
#     cab2[para] = new_para[i]
# 
# hdul = fits.HDUList()
# hdul.append(fits.PrimaryHDU(header=cab0))
# hdul.append(fits.ImageHDU(imagen[1].data, header=cab1,name='DATA'))
# hdul.writeto(pruebas + 'mapping_ima_%s_align.fits'%(period),overwrite=True)
# =============================================================================

# %%
# 

















