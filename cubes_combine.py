#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 12:08:37 2023

@author: alvaromartinez
"""

# We have fitted the sky in three different windows with different wavelehts.
# This scrips cuts the cubes according with  the wavelehts and combine then in a 
# single one

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
cubA = '/Users/alvaromartinez/Desktop/Phd/KMOS/p105_A/2023-05-23T13:53:08/'
# cubB = '/Users/alvaromartinez/Desktop/Phd/KMOS/p105_B/2023-05-23T17:51:34/'
cubB = '/Users/alvaromartinez/Desktop/Phd/KMOS/p105_B/2023-09-25T17:17:01/'
cubC = '/Users/alvaromartinez/Desktop/Phd/KMOS/p105_C/2023-05-24T12:54:32/'
cubABC= '/Users/alvaromartinez/Desktop/Phd/KMOS/p105_ABC/2023-05-20T22:05:16/'
name = 'COMBINE_SKY_TWEAK_mapping.fits'
pruebas = '/Users/alvaromartinez/Desktop/Phd/KMOS/pruebas/'
cubs =[cubA, cubB, cubC]

a_win = [1.975,1.987,1.993,2.010,2.041,2.060]
abc_win = [2.060, 2.269]
b_win = [2.269,2.291,2.308,2.335,2.360,2.379]
c_win = [2.416,2.440,2.445,2.475]

fitA = fits.open(cubA + name)
fitB = fits.open(cubB + name)
fitC = fits.open(cubC + name)
fitABC =  fits.open(cubABC + name)
lam_0 = fitA[1].header['CRVAL3']
step =  fitA[1].header['CDELT3']

lam_all =np.array([lam_0 + step*i for i in range(0,3072)])

# %%
# cut_A_up = np.where(abs(lam_all - a_win[-1]) == min(abs(lam_all - a_win[-1])))
# cut_ABC_down = cut_A_up
# cut_ABC_up = np.where(abs(lam_all - abc_win[-1]) == min(abs(lam_all - abc_win[-1])))
# cut_B_up = np.where(abs(lam_all - b_win[-1]) == min(abs(lam_all - b_win[-1])))
# cut_C_down = np.where(abs(lam_all - c_win[-1]) == min(abs(lam_all - b_win[-1])))


cut_A_up = np.where(abs(lam_all - a_win[-1]) == min(abs(lam_all - a_win[-1])))
cut_B_up = np.where(abs(lam_all - b_win[-1]) == min(abs(lam_all - b_win[-1])))


fit_all = np.zeros((fitC[1].data.shape[0],fitC[1].data.shape[1],fitC[1].data.shape[2]))
fit_all_noise = np.zeros((fitC[1].data.shape[0],fitC[1].data.shape[1],fitC[1].data.shape[2]))


fit_all[0:cut_A_up[0][0],:,:] = fitA[1].data[0:cut_A_up[0][0],:,:] 
fit_all[cut_A_up[0][0]:cut_B_up[0][0],:,:]= fitB[1].data[cut_A_up[0][0]:cut_B_up[0][0],:,:]
fit_all[cut_B_up[0][0]:,:,:] = fitC[1].data[cut_B_up[0][0]:,:,:]

fit_all_noise[0:cut_A_up[0][0],:,:] = fitA[2].data[0:cut_A_up[0][0],:,:] 
fit_all_noise[cut_A_up[0][0]:cut_B_up[0][0],:,:]= fitB[2].data[cut_A_up[0][0]:cut_B_up[0][0],:,:]
fit_all_noise[cut_B_up[0][0]:,:,:] = fitC[2].data[cut_B_up[0][0]:,:,:]

# A_wl = fitA[1].data[0:cut_A_up[0][0],:,:] 
# B_wl = fitB[1].data[cut_A_up[0][0]:cut_B_up[0][0],:,:]
# C_wl = fitC[1].data[cut_B_up[0][0]:,:,:]

# %%
ima, header1 = fits.getdata(cubA + name, header = True)
hh = fits.open(cubA + name)[0]
paras = ['NAXIS1',	'NAXIS2',	'CRPIX1',	'CRPIX2',	
        'CRVAL1',	'CRVAL2',	'CD1_1',	'CD1_2',	
        'CD2_1',	'CD2_2']
new_para =[650,         433,312.7286234667807,319.09363905364944,
           266.384440107102, -28.9444025,6.958225445320778E-6,2.6697138105672696E-5
           ,2.6697138105672696E-5, -6.958225445320778E-6]
for i, para in enumerate(paras):
    print(header1[para])
    header1[para] = new_para[i]
# I got this parameters for the header aligning with Aladin
# NAXIS1  = 650
# NAXIS2  = 433
# CRPIX1  = 312.7286234667807
# CRPIX2  = 319.09363905364944
# EQUINOX = 2000.0
# CRVAL1  = 266.384440107102
# CRVAL2  = -28.9444025
# CTYPE1  = RA---TAN
# CTYPE2  = DEC--TAN
# RADECSYS= ICRS
# CD1_1   = 6.958225445320778E-6
# CD1_2   = 2.6697138105672696E-5
# CD2_1   = 2.6697138105672696E-5
# CD2_2   = -6.958225445320778E-6
# %%

new_hdul = fits.HDUList()
new_hdul.append(fits.PrimaryHDU(header=header1))
new_hdul.append(fits.ImageHDU(fit_all, header=header1,name='DATA'))
new_hdul.append(fits.ImageHDU(fit_all_noise,header=header1, name='ERROR'))
new_hdul.writeto(pruebas + 'three_headers_2.fits',overwrite=True)
sys.exit()
# # %%
# # 
# # hdu2 = fits.ImageHDU([fitA[1].data])
# # hdu3 = fits.ImageHDU([fitA[2].data])
# # test_HDU = fits.HDUList([hdu2,hdu3])
# # hdu1 = fits.PrimaryHDU([hh])
# hdu2 = fits.ImageHDU()
# hdu3 = fits.ImageHDU()
# test_HDU = fits.HDUList([hh,hdu2,hdu3])
# test_HDU.writeto(pruebas + 'three_headers.fits', hh,overwrite=True)

# fits.update(pruebas + 'three_headers.fits',fitA[1].data,1)
# fits.update(pruebas + 'three_headers.fits',fitA[2].data,2)
# # %%
# fits.update(pruebas + 'three_headers.fits',header1,0)
# # %%
tri_hea = fits.open(pruebas + 'three_headers.fits')

# # fits.writeto(pruebas + 'three_headers.fits',test_HDU,header =header, overwrite=True)
# # %%
# %%
fits.writeto(pruebas + 'test_ABC_3.fits', fit_all,overwrite= True)
fits.append(pruebas + 'test_ABC_3.fits', header = header1)
fits.append(pruebas + 'test_ABC_3.fits', fit_all_noise)


tri_hea = fits.open(pruebas + 'test_ABC_3.fits')
sys.exit('90')
fits.writeto(pruebas + 'test_ABC.fits', fit_all, header = header, overwrite= True)
fit_one = np.zeros((1,fitC[1].data.shape[1],fitC[1].data.shape[2]))
fit_one = fit_all[1271:1272,:,:]
fits.writeto(pruebas + 'one_ABC.fits', fit_one, header = header, overwrite= True)


# %%

print(header1)





