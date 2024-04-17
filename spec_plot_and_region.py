#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 12:04:02 2023

@author: amartinez
"""

# Region file of young or old or all stars. Also independet ploting of the spectra

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
from matplotlib.patches import Circle
import IPython
import tkinter as tk
from tkinter import simpledialog
import shutil
from matplotlib.widgets import Slider
import random
import re


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
# Enable automatic plotting mode
IPython.get_ipython().run_line_magic('matplotlib', 'auto')
# IPython.get_ipython().run_line_magic('matplotlib', 'inline')

reduction = 'ABC'
pruebas = '/Users/amartinez/Desktop/PhD/KMOS/practice/'
aling = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/ifu_alignment_%s/'%('ABC')
log_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/'%('ABC')
esorex_cube_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/COMBINE_SKY_TWEAK_mapping.fits'%(reduction)
esorex_ima_5= '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/comoving_group_mosaic_K_Half1_COMBINED_IMAGE_mapping.fits'%(reduction)
log_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/'%('ABC')


pruebas = '/Users/amartinez/Desktop/PhD/KMOS/practice/'
reductions = ['ABC']
# reductions = ['tramos']
ifu_sel_ls = np.arange(1,24)
half_ifu_ls = [1,2]
# spec_fol = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/ABC_reduction/young_candidates/ifu_%s/half_%s'%(ifu_sel,half_ifu)

# cell_n = []
# for j in ifu_sel_ls:
#     for k in half_ifu_ls:
#         spec_fol = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/cluster_spectra/ifu_%s/half_%s'%(reductions[0],j,k)
#         cell_n.append(len(glob.glob(spec_fol +'/spec*.fits')))
# cells = sum(cell_n)
# 
age = 'young_candidates'
# age = 'cluster_spectra'
all_stars = []
x_y_full =[] 
for ifu_sel in ifu_sel_ls:
    colores = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])]
    for half_ifu in half_ifu_ls:
        
        spec_fol = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/%s/ifu_%s/half_%s/'%(reductions[0],age,ifu_sel,half_ifu)
        isdir = os.path.isdir(spec_fol)
        if not isdir or len(glob.glob(spec_fol +'/spec*')) == 0:
        
            print('No IFU %s'%(ifu_sel))
            continue
        else:
            all_stars.append(len(glob.glob(spec_fol +'/spec*')))
            for st in range(len(glob.glob(spec_fol +'/spec*'))):
                name = glob.glob(spec_fol +'/spec*.fits')[st]
                print(name)
                sp = fits.open(name)
                sp_hed = sp[0].header
                x_y_full.append((sp_hed['X_FULL']+1,sp_hed['Y_FULL']+1))
# np.savetxt(spec_fol + 'xy_young_ifu%s_half%s.txt'%(ifu_sel, half_ifu),np.array(x_y), fmt = '%.0f')
# np.savetxt(pruebas + 'xy_young.txt',np.array(x_y), fmt = '%.0f')
print(np.sum(all_stars))

# %%

ifu_sel_ls = np.arange(6,7)
half_ifu_ls = [2]



w1 = 2
w2 = 2



HeI = 2.058
COI =  2.29322
COII = 2.32246
COIII = 2.3525
Brg = 2.165
He = 2.12
HeII = 2.189
# H2 = 2.12
l_names = ['HeI', 'COI', 'COII','Br$\gamma$', 'He', 'HeII','C0III']
lines = [HeI, COI, COII,Brg, He, HeII,COIII]





# plt.draw()
# %%
dic_young = {}
for ifu_sel in ifu_sel_ls:
    colores = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])]
    for half_ifu in half_ifu_ls:
        aligned_cube = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/cubes/cube_ifu%s_half%s_mean.fits'%(reduction,ifu_sel,half_ifu)
        
        fig, ax = plt.subplots(1,2)
        ima = fits.open(aligned_cube)
        cube = ima[1].data
        cube_plot = np.nanmedian(cube, axis = 0)
        
        def update(val):
            vmin = vmin_slider.val
            vmax = vmax_slider.val
            im.set_clim(vmin, vmax)
            plt.draw()
        global im
        v_min =np.nanmin(cube_plot)
        v_max =np.nanmax(cube_plot)
        im = ax[1].imshow(cube_plot, cmap='Greys',label='overlays',origin='lower',vmin=v_min,vmax = v_max/5,alpha =1)
        # fig.colorbar(im, ax=ax[1], orientation='vertical')
        ax_vmin = plt.axes([0.6, 0.1, 0.15, 0.03], facecolor='lightgoldenrodyellow')
        vmin_slider = Slider(ax_vmin, 'vmin', v_min, v_max, valinit=v_min)

        # Add a slider for vmax
        ax_vmax = plt.axes([0.6, 0.05, 0.15, 0.03], facecolor='lightgoldenrodyellow')
        vmax_slider = Slider(ax_vmax, 'vmax', v_min, v_max, valinit=v_max)
        
        vmin_slider.on_changed(update)
        vmax_slider.on_changed(update)

       

        # plt.draw()
        difu_fol = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/%s/ifu_%s/half_%s/'%(reductions[0],'cluster_spectra',ifu_sel,half_ifu)
        spec_fol = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/%s/ifu_%s/half_%s/'%(reductions[0],age,ifu_sel,half_ifu)
        isdir = os.path.isdir(spec_fol)





        if not isdir or len(glob.glob(spec_fol +'/spec*')) == 0:
        
            print('No IFU %s'%(ifu_sel))
            continue
        else:
            all_stars.append(len(glob.glob(spec_fol +'/spec*')))
            difuse_flux = fits.open(difu_fol  + 'sky_ifu%s_half%s.fits'%(ifu_sel,half_ifu))[0].data
            for st in range(len(glob.glob(spec_fol +'/spec*'))):
                name = glob.glob(spec_fol +'/spec*.fits')[st]
                coord = re.findall(r'\d+', name)
                x,y = np.array(coord[4]).astype(int), np.array(coord[5]).astype(int)
                sc =ax[1].scatter(x,y, marker = 'x', s = 200, label = 'x,y=%s,%s'%(x,y),)
                print(sc.get_facecolors())
                ax[1].legend(fontsize = 8)
                print(name)
                print(x,y)
                ax[0].set_title(os.path.basename(name))
                sp = fits.open(name)
                cab = sp[0].header
                # sp_data = sp[0].data[0]
                sp_data = sp[0].data
                lam = np.array([cab['CRVAL1']+cab['CDELT1']*i for i in range(cab['NAXIS1'])])
                # d_spec = lam[0]
                # u_spec = lam[-1]
                # good = np.where((lam > d_spec) & (lam < u_spec))
                # lam = lam[good]
                ax[0].plot(lam,sp_data, label = 'x=%s\ny=%s'%(cab['X_FULL']+1,cab['Y_FULL']+1))

                for l in lines:
                    ax[0].axvline(l, color = 'grey', alpha = 0.5, ls = 'dashed')
                ax2 = ax[0].twiny()
                ax2.set_xticks(lines)
                tl = l_names
                ax2.set_xticklabels(tl)
                ax2.plot(lam,sp_data, alpha = 0)
                
                dic_young['s_%s_%s'%(x,y)] = sp_data
                # sig_w = 3
                # ind1 = np.where(abs(lam-He)== min(abs(lam-He)))
                # ind2 = np.where(abs(lam-Brg)== min(abs(lam-Brg)))
                # ax2.plot(lam,difuse_flux, color ='k', alpha = 0.3,zorder = 1, label ='Difuse emission')
                # ax2.axhline(np.nanmean(difuse_flux[ind1[0][0]:ind2[0][0]]), alpha = 0.3, color = 'k')
                # sig = np.nanstd(difuse_flux[ind1[0][0]:ind2[0][0]])
                # ax[0].axhspan(np.nanmean(difuse_flux[ind1[0][0]:ind2[0][0]])- sig*sig_w,
                #             np.nanmean(difuse_flux[ind1[0][0]:ind2[0][0]])+ sig*sig_w, color = 'k', alpha = 0.1)
                # ax[0].set_ylim(np.nanmean(difuse_flux[ind1[0][0]:ind2[0][0]])- sig*100,
                #             np.nanmean(difuse_flux[ind1[0][0]:ind2[0][0]])+ sig*100)
         
                # ax[0].legend()
plt.show()
# %%

leny = len(dic_young)
fig, ax = plt.subplots(leny, 1, sharex=True)
# Remove vertical space between axes
fig.subplots_adjust(hspace=0)

ls_y = list(dic_young.keys())
to_plot = np.where((lam>HeI)&(lam<COIII))
if leny >1:
    for i in range(leny):
        
        ax[i].plot(lam[to_plot],dic_young[ls_y[i]][to_plot])
        ax[i].set_xlim(HeI, COIII)
        for l in lines:
            ax[i].axvline(l, color = 'grey', alpha = 0.5, ls = 'dashed') 
    
    
    ax2 = ax[0].twiny()
    ax2.set_xlim(HeI, COIII)
    ax2.set_xticks(lines)
    tl = l_names
    ax2.set_xticklabels(tl)
    ax2.plot(lam[to_plot],dic_young[ls_y[0]][to_plot], alpha = 0)
else: 
    ax.plot(lam[to_plot],dic_young[ls_y[0]][to_plot])
    ax.set_xlim(HeI, COIII)
    for l in lines:
        ax.axvline(l, color = 'grey', alpha = 0.5, ls = 'dashed') 
    ax2 = ax.twiny()
    ax2.set_xlim(HeI, COIII)
    ax2.set_xticks(lines)
    tl = l_names
    ax2.set_xticklabels(tl)
    ax2.plot(lam[to_plot],dic_young[ls_y[0]][to_plot], alpha = 0)

# ax[0].set_xticks(lines)
# ax[0].plot(lam[to_plot],dic_young[ls_y[0]][to_plot], alpha = 0)
# tl = l_names
# ax[0].set_xticklabels(tl)

# %%
# Enable automatic plotting mode
IPython.get_ipython().run_line_magic('matplotlib', 'auto')
# IPython.get_ipython().run_line_magic('matplotlib', 'inline')
fig, ax = plt.subplots()
ax.plot(lam, sp_data[0])
ax.plot(lam, sp_data[1])
# plt.legend()
# %%
