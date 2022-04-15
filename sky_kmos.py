#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 13:40:20 2022

@author: amartinez
"""

from astropy.io import fits
from astropy import units as u
import numpy as np
from matplotlib import pyplot as plt
from astropy.visualization import quantity_support
from matplotlib.colors import LogNorm
import glob
import os
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.data import get_pkg_data_fileobj
from astropy.utils.data import get_pkg_data_contents
import random

morralla ='/Users/amartinez/Desktop/morralla/'

pruebas='/Users/amartinez/Desktop/PhD/KMOS/practice/'
sky='/Users/amartinez/Desktop/PhD/KMOS/data/all_data_and_all_sky/end_products/2022-01-26T12:35:36/KMOS.2021-06-03T04:47:47.641_combine_OBs/'
sky_all='/Users/amartinez/Desktop/PhD/KMOS/data/all_data_and_all_sky/end_products/2022-04-13T08:43:35/KMOS.2021-06-03T04:47:47.641_tpl/'
sky_masked='/Users/amartinez/Desktop/PhD/KMOS/data/all_data_and_all_sky/end_products/sky_all_masked/'

# s_file = get_pkg_data_filename(sky_all+'GC_K_Sky_SCI_RECONSTRUCTED_KMOS.2021-06-03T04:47:47.641.fits')
#for some reason this only work if a punt the file in the same directory as the script
#but only when using spyder, in jupyter notebooks works just fine
period=['03','06']
half =[1,2]
for p in period:
    for h in half:
        s_file=get_pkg_data_filename('sky_%s_half%s.fits'%(p,h))

# %%
 
        s_header=fits.getheader('sky_%s_half%s.fits'%(p,h))
        
        hdu=fits.open(s_file)
        
        dic_ifu={}
        ifu_valid=[1,3,4,5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16,17, 19, 21, 23]#manually choosen from the fits file
        for i in ifu_valid:
        #     print(i)
            dic_ifu['ifu%s'%(i)]=hdu[i].data
        #  ifu15 also looks very bad in the target images...
        # dic_ifu['ifu15']=dic_ifu['ifu17']#ifu 15 looks very bad, in the sky, lets try to replace it with another one that looks better
        #%%
        sky_median_all=[]
        sky_std_all=[]
        sky_mean_all=[]
        for i in range(0,dic_ifu['ifu1'].shape[0]):
            sky_ifu=[]
            sky_ifu=[dic_ifu['ifu%s'%(j)][i,:,:] for j in ifu_valid]
            sky_ifu=np.array(sky_ifu)
            sky_median=np.nanmedian(sky_ifu)
            sky_std=np.nanstd(sky_ifu)
            sky_mean=np.nanmean(sky_ifu)
            
            sky_median_all.append(sky_median)
            sky_std_all.append(sky_std)
            sky_mean_all.append(sky_mean)
        #%%
        
        # Loops over every slice on the cube(2048x14x14)xlen(valid_ifu). 
        # On each slice(14X14) compares each pixel with the median for the same slice(same position) in the rest of the ifus, i.e., the median of (14X14) of the slice in position cube_deepth of every valid ifu. 
        # The values of the median of each slice are store in a variable called sky_median_all.
        # If the value of this pixel is biger > meidan + 3*std(of the values on all the same slice for each ifu), it is replaced by the median plus 1std*rdm_number.
        #
        for cube_depth in range(len(sky_median_all)):
            if np.isnan(sky_median_all[cube_depth])==False:
                for j in ifu_valid:
                    for x in range(14):
                        for y in range(14):
                            if (np.isnan(hdu[j].data[cube_depth,x,y])==True): 
                                hdu[j].data[cube_depth,x,y] = sky_std_all[cube_depth]*np.random.random() + sky_median_all[cube_depth]
                            elif (hdu[j].data[cube_depth,x,y] > sky_median_all[cube_depth] + 1*sky_std_all[cube_depth]):
                               hdu[j].data[cube_depth,x,y ]= sky_std_all[cube_depth]*np.random.random() + sky_median_all[cube_depth]
        #%%
        new_hdul = fits.HDUList()
        new_hdul.append(fits.PrimaryHDU(header=hdu[0].header))
        for i in range(1,len(hdu)):
            new_hdul.append(fits.ImageHDU(hdu[i].data,header=hdu[i].header))
        new_hdul.writeto(sky_masked+'sky_%s_half%s_masked.fits'%(p,h),overwrite=True)
        
        print('Saved fits file period %s half %s'%(p,h))
#%%


