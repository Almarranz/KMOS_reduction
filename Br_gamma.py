#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 12:26:51 2023

@author: amartinez
"""
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
from matplotlib.colors import LogNorm
from matplotlib.tri import Triangulation
import regions
from matplotlib.patches import Rectangle
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

# We are going to make a Br_g emision map
reduction = 'ABC'
pruebas = '/Users/amartinez/Desktop/PhD/KMOS/practice/'
aling = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/ifu_alignment_%s/'%('ABC')
log_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/'%('ABC')
esorex_cube_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/COMBINE_SKY_TWEAK_mapping.fits'%(reduction)
esorex_ima_5= '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/comoving_group_mosaic_K_Half1_COMBINED_IMAGE_mapping.fits'%(reduction)
ima = fits.open(esorex_cube_5)
mapa = WCS(ima[1].header).celestial
h0 = ima[0].header
h1 = ima[1].header

ima_1 = fits.open(esorex_ima_5)
h0_1 = ima_1[0].header
h1_1 = ima_1[1].header
# %%
paras = ['NAXIS1',	'NAXIS2',	
          'CRPIX1',	'CRPIX2',	
        'CRVAL1',	'CRVAL2',	
        'CD1_1',	'CD1_2',	
        'CD2_1',	'CD2_2']

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


new_para =[650,         433,
            CRPIX1,  CRPIX2,
            CRVAL1, CRVAL2,
            CD1_1,   CD1_2
            ,CD2_1, CD2_2]
for i, para in enumerate(paras):
    h1[para] = new_para[i]
    h1_1[para] = new_para[i]
    

# %%
temp5 = np.zeros((ima[1].data.shape[1],ima[1].data.shape[2]))
# %%
data = ima[1].data
brg = np.mean(data[1286:1288,:,:], axis =0) 
brg_mean = np.nanmean((data - brg),axis = 0)
# %%
hdul = fits.HDUList()
hdul.append(fits.PrimaryHDU(header = h0))
hdul.append(fits.ImageHDU(brg_mean*-1,h1, name ='Brg emission'))
hdul.writeto(pruebas  + 'brg_mean.fits', overwrite = True)
wcs = WCS(h1)
clus_o = np.loadtxt('/Users/amartinez/Desktop/PhD/ESO_proposals/KMOS/ID_AB.txt')
clus_young = np.loadtxt('/Users/amartinez/Desktop/PhD/ESO_proposals/KMOS/young_kmos_bg.txt', usecols=(0,1), unpack= True)

clus_xy = wcs.wcs_world2pix(clus_o[:,0],clus_o[:,1],1,1)
clus_young_xy = wcs.wcs_world2pix(clus_young[0],clus_young[1],1,1)
# %%

data = brg_mean*-1
fill_value = 0.0  # Replace with an appropriate fill value
data = np.nan_to_num(data, nan=fill_value, posinf=fill_value, neginf=fill_value)


region_name = pruebas + 'ds9.reg'
from regions import Regions
r = Regions.read(pruebas + 'ds9.reg',format='ds9')

# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
v_min =1
v_max = 0
data = brg_mean*-1
# im = ax.imshow(data*1e17, cmap='Greys', origin='lower',vmin=v_min,vmax = v_max)
im = ax.imshow(ima_1[1].data, cmap='Greys', origin='lower',vmin=0.001e-16,vmax = 0.1e-16)
for i in range(6):
    ax.plot(r[i].vertices.x,r[i].vertices.y, color = 'lime')
# %
xticks = ax.get_xticks()
yticks = ax.get_yticks()

x_t = np.arange(-100,900,100)
y_t = np.arange(-100,900,100)
coor_tiks = wcs.wcs_pix2world(x_t,y_t,1,1)
# ax.set_xticklabels(coor_tiks[0])
# ax.ax.set_xticklabels(coor_tiks[1])

new_x_tick_labels = np.round(coor_tiks[0],3)  # Replace with your desired labels
new_y_tick_labels =  np.round(coor_tiks[1],3)  # Replace with your desired labels

# Set the new tick labels for the x and y axes
ax.set_xticklabels(new_x_tick_labels)
ax.set_yticklabels(new_y_tick_labels)

# ax.tricontour(x,y,z, levels=niveles, colors='white')
# fig.colorbar(im,ax =ax,shrink=0.7)
# cbar = fig.colorbar(im, orientation='vertical')

# for p in range(15):
#     ax.add_patch(patch_list[p])
ax.scatter(clus_xy[0],clus_xy[1],s = 200,color ='lime', label = 'Co-moving group',edgecolor = 'k',lw =1)
ax.scatter(clus_young_xy[0],clus_young_xy[1],color = 'fuchsia',label = 'No CO lines')

lgnd = ax.legend()
for handle in lgnd.legend_handles:
    handle.set_sizes([200])

ax.add_patch(Rectangle((0, 435), 650, 400, alpha = 0.5, color = 'k',fill = False, lw =5,ls ='dashed',label = 'ID'))
# ax.legend()
ax.set_xlim(-50, 670)
ax.set_ylim(-50, 845)
ax.set_xlabel('Ra (°)')
ax.set_ylabel('Dec (°)')
plt.savefig(pruebas  + 'br_map.png', bbox_inches = 'tight')
# %%




# %
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
v_min =1
v_max = 0
data = brg_mean*-1
im = ax.imshow(data*1e17, cmap='Greys', origin='lower')
# im = ax.imshow(ima_1[1].data, cmap='Greys', origin='lower',vmin=v_min,vmax = v_max)

