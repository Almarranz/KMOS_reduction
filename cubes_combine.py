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
# %%
cubA = '/Users/alvaromartinez/Desktop/Phd/KMOS/p105_A/2023-05-23T13:53:08/'
cubB = '/Users/alvaromartinez/Desktop/Phd/KMOS/p105_B/2023-05-23T17:51:34/'
cubC = '/Users/alvaromartinez/Desktop/Phd/KMOS/p105_C/2023-05-24T12:54:32/'
name = 'COMBINE_SKY_TWEAK_mapping.fits'
pruebas = '/Users/alvaromartinez/Desktop/Phd/KMOS/pruebas/'
cubs =[cubA, cubB, cubC]

a_win = [1.975,1.987,1.993,2.010,2.041,2.060]
b_win = [2.269,2.291,2.308,2.335,2.360,2.379]
c_win = [2.416,2.440,2.445,2.475]
fitA = fits.open(cubA + name)
fitB = fits.open(cubB + name)
fitC = fits.open(cubC + name)
lam_0 = fitA[1].header['CRVAL3']
step =  fitA[1].header['CDELT3']

lam_all =np.array([lam_0 + step*i for i in range(0,3072)])
# %%
cut_A_up = np.where(abs(lam_all - a_win[-1]) == min(abs(lam_all - a_win[-1])))
cut_B_up = np.where(abs(lam_all - b_win[-1]) == min(abs(lam_all - b_win[-1])))
cut_C_down = np.where(abs(lam_all - c_win[-1]) == min(abs(lam_all - b_win[-1])))

fit_all = np.zeros((fitC[1].data.shape[0],fitC[1].data.shape[1],fitC[1].data.shape[2]))


fit_all[0:cut_A_up[0][0],:,:] = fitA[1].data[0:cut_A_up[0][0],:,:] 
fit_all[cut_A_up[0][0]:cut_B_up[0][0],:,:]= fitB[1].data[cut_A_up[0][0]:cut_B_up[0][0],:,:]
fit_all[cut_B_up[0][0]:,:,:] = fitC[1].data[cut_B_up[0][0]:,:,:]

# A_wl = fitA[1].data[0:cut_A_up[0][0],:,:] 
# B_wl = fitB[1].data[cut_A_up[0][0]:cut_B_up[0][0],:,:]
# C_wl = fitC[1].data[cut_B_up[0][0]:,:,:]

# %%
ima, header = fits.getdata(cubA + name, header = True)

paras = ['NAXIS1',	'NAXIS2',	'CRPIX1',	'CRPIX2',	
        'CRVAL1',	'CRVAL2',	'CD1_1',	'CD1_2',	
        'CD2_1',	'CD2_2']
new_para =[650,         433,312.7286234667807,319.09363905364944,
           266.384440107102, -28.9444025,6.958225445320778E-6,2.6697138105672696E-5
           ,2.6697138105672696E-5, -6.958225445320778E-6]
for i, para in enumerate(paras):
    print(header[para])
    header[para] = new_para[i]
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

fits.writeto(pruebas + 'test_ABC.fits', fit_all, header = header, overwrite= True)
fit_one = np.zeros((1,fitC[1].data.shape[1],fitC[1].data.shape[2]))
fit_one = fit_all[1271:1272,:,:]
fits.writeto(pruebas + 'one_ABC.fits', fit_one, header = header, overwrite= True)


# %%







