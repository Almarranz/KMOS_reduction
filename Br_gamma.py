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
from regions import Regions
import IPython
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

# Enable automatic plotting mode
# IPython.get_ipython().run_line_magic('matplotlib', 'auto')
IPython.get_ipython().run_line_magic('matplotlib', 'inline')



# We are going to make a Br_g emision map
reduction = 'ABC'
pruebas = '/Users/amartinez/Desktop/PhD/KMOS/practice/'
aling = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/ifu_alignment_%s/'%('ABC')
log_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/'%('ABC')
esorex_cube_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/COMBINE_SKY_TWEAK_mapping.fits'%(reduction)
esorex_ima_5= '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/comoving_group_mosaic_K_Half1_COMBINED_IMAGE_mapping.fits'%(reduction)
cube_no_cont = pruebas + 'cube_no_cont_all.fits'

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

# %%
data = fits.open(cube_no_cont)[0].data

brg_mean = np.mean(data[1285:1291,:,:], axis =0) 
hdu_brg = fits.PrimaryHDU(data = brg_mean, header = h1)
hdu_brg.writeto(pruebas +  'brg_emission.fits', overwrite= True)

heI_mean = np.mean(data[711:716,:,:], axis = 0)
hdu_heI = fits.PrimaryHDU(data = heI_mean, header = h1)
hdu_heI.writeto(pruebas +  'heI_emission.fits', overwrite= True)

coI = np.arange(1969,1980)
coII = np.arange(2124,2134)
co = np.r_[coI,coII]
co_mean = np.mean(data[co,:,:], axis = 0)
hdu_co = fits.PrimaryHDU(data = co_mean, header = h1)
hdu_co.writeto(pruebas + 'co_absortion.fits', overwrite = True) 
# %%
# cube_mean = np.nanmedian(ima[1].data[:,:,:], axis = 0)
# hdu_cube = fits.PrimaryHDU(data = cube_mean, header = h1)
# hdu_cube.writeto(pruebas + 'cube_mean.fits', overwrite = True)
# %%
wcs = WCS(h1)
# clus_o = np.loadtxt('/Users/amartinez/Desktop/PhD/ESO_proposals/KMOS/ID_AB.txt')
clus_o = np.loadtxt('/Users/amartinez/Desktop/PhD/ESO_proposals/KMOS/p1_p113/ID_AB.txt')
# clus_young = np.loadtxt('/Users/amartinez/Desktop/PhD/ESO_proposals/KMOS/young_kmos_bg.txt', usecols=(0,1), unpack= True)

clus_xy = wcs.wcs_world2pix(clus_o[:,0],clus_o[:,1],1,1)
# clus_young_xy = wcs.wcs_world2pix(clus_young[0],clus_young[1],1,1)
# %%
ifu_sel_ls = np.arange(6,8)
half_ifus = [1,2]
reductions = ['ABC']
age = 'young_candidates'
# age = 'cluster_spectra'

clus_young_xy = np.loadtxt(pruebas + 'xy_young.txt')
# for  ifu_sel in ifu_sel_ls:
#     for half_ifu in half_ifus:
#         spec_fold = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/%s/ifu_%s/half_%s/'%(reductions[0],age,ifu_sel,half_ifu)

#         xy = np.loadtxt(spec_fold + 'xy_young_ifu%s_half%s.txt'%(ifu_sel, half_ifu),
#                            usecols=(0,1), unpack= True)
#         clus_young_xy.append(xy)
# sys.exit('163')    
# %%

r = Regions.read(pruebas + 'Brg_region.reg',format='ds9')

# %%
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
v_min =1
v_max = 0
# im = ax.imshow(data*1e17, cmap='Greys', origin='lower',vmin=v_min,vmax = v_max)
im = ax.imshow(ima_1[1].data, cmap='Greys', origin='lower',vmin=0.001e-16,vmax = 0.1e-16)
for i in range(17):
    ax.plot(r[i].vertices.x,r[i].vertices.y, color = 'white')
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
# %
# Set the new tick labels for the x and y axes
ax.set_xticklabels(new_x_tick_labels)
ax.set_yticklabels(new_y_tick_labels)
# sys.exit(193)
# %%
# ax.tricontour(x,y,z, levels=niveles, colors='white')
# fig.colorbar(im,ax =ax,shrink=0.7)
# cbar = fig.colorbar(im, orientation='vertical')

# for p in range(15):
#     ax.add_patch(patch_list[p])
ax.scatter(clus_xy[0],clus_xy[1],s = 200,color ='lime', label = 'Co-moving group',edgecolor = 'k',lw =1)
ax.scatter(clus_young_xy[:,0],clus_young_xy[:,1],color = 'fuchsia',label = 'No CO lines')

lgnd = ax.legend()
for handle in lgnd.legend_handles:
    handle.set_sizes([200])

ax.add_patch(Rectangle((0, 435), 650, 400, alpha = 0.5, color = 'k',fill = False, lw =5,ls ='dashed',label = 'ID'))
# ax.legend()
ax.set_xlim(-50, 670)
ax.set_ylim(-50, 845)
ax.set_xlabel('Ra (°)')
ax.set_ylabel('Dec (°)')
# plt.savefig(pruebas  + 'br_map.png', bbox_inches = 'tight')
# %%

fig, ax = plt.subplots(1, 1, figsize=(10, 10))
v_min =1
v_max = 0
im = ax.imshow(brg_mean, cmap='Greys', origin='lower',vmin=0.001e-17,vmax = 0.55e-17)
# im = ax.imshow(heI_mean, cmap='Greys', origin='lower',vmin=0e-17,vmax = 1e-17)
for i in range(17):
    ax.plot(r[i].vertices.x,r[i].vertices.y, color = 'lime')
# im = ax.imshow(ima_1[1].data, cmap='Greys', origin='lower',vmin=v_min,vmax = v_max)
fig.colorbar(im, ax=ax, orientation='vertical')



