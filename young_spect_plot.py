#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 13:51:29 2024

@author: amartinez
"""

# Plottim the spectra for the YSO with good SNR found in the whole KMOS fov

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
import numpy as np
from astropy.stats import sigma_clip
from numpy import mean

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
# IPython.get_ipython().run_line_magic('matplotlib', 'auto')
IPython.get_ipython().run_line_magic('matplotlib', 'inline')

# reduction = 'ABC'
reduction = 'tramos'

pruebas = '/Users/amartinez/Desktop/PhD/KMOS/practice/'
aling = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/ifu_alignment_%s/'%('ABC')
log_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/'%('ABC')
esorex_cube_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/COMBINE_SKY_TWEAK_mapping.fits'%(reduction)
esorex_ima_5= '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/comoving_group_mosaic_K_Half1_COMBINED_IMAGE_mapping.fits'%(reduction)
log_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/'%('ABC')


HeI = 2.058
COI =  2.29322
COII = 2.32246
COIII = 2.3525
Brg = 2.165
He = 2.12
HeII = 2.189
# H2 = 2.12
# l_names = ['HeI', '$^{12}$CO(2,0)', ' $^{12}$CO(3,1)\n','Br$\gamma$', 'He', 'HeII','$^{12}$CO(4,2)']
# lines = [HeI, COI, COII,Brg, He, HeII,COIII]

l_names = [ '$^{12}$CO(2,0)', ' $^{12}$CO(3,1)\n','Br$\gamma$', 'He','$^{12}$CO(4,2)']
lines = [COI, COII,Brg, He,COIII]

age = 'young_candidates'
# reductions = ['ABC']
reductions = ['tramos']
ls_spec = np.loadtxt(pruebas + 'ls_young.txt')
# ls_spec = np.delete(ls_spec,1, axis =0)
colorines = ['#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2','#7f7f7f','#8c564b', '#e377c2','#7f7f7f']

fig, ax = plt.subplots(len(ls_spec), 2, figsize = (30,18))
fig.subplots_adjust(hspace=0)

lim_d, lim_u = 2.15, COIII +0.01
# lim_d, lim_u = 2, COIII +0.01#!!!

norm_0, norm_1 = HeI, He
d_norm, u_norm =  2.2,2.28 #!!!
snr_d, snr_u = 2.2356, 2.2661#!!!
dic_yso = {}
dic_yso_snr = {}
for i in range(len(ls_spec)):
    ax[i,0].set_xlim(lim_d, lim_u)
    ifu_sel = ls_spec[i][0]
    half_ifu = ls_spec[i][1]
    x, y = ls_spec[i][2], ls_spec[i][3]
    spec_fol = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/%s/ifu_%.0f/half_%.0f/'%(reductions[0],age,ifu_sel,half_ifu)
    spec = fits.open(spec_fol + 'spec_ifu%.0f_half%.0f_%.0f_%.0f.fits'%(ifu_sel, half_ifu, x, y))
    cab = spec[0].header
    sp_data = spec[0].data
    lam = np.array([cab['CRVAL1']+cab['CDELT1']*i for i in range(cab['NAXIS1'])])
    dic_yso_snr['Y%s'%(i+1)] = np.array([lam,sp_data])
    lam_selc = np.where((lam > norm_0) &(lam<norm_1))
    
    good = np.where((lam > lim_d) & (lam < lim_u))
    norm =  np.where((lam > d_norm) & (lam < u_norm))
    filtered = sigma_clip(sp_data[norm], sigma=1.5, maxiters=None,cenfunc=mean, masked=False, copy=False)
    sp_data, lam  = sp_data[good], lam[good]
    sp_data = sp_data/np.mean(filtered)
    
    
    
    flat_flux = np.nanmean(sp_data[lam_selc])
    # print(flat_flux)
    # ax[i].set_ylim(flat_flux*0.5, flat_flux*2)
    ax[i,0].set_ylim(0, 2)
     # ax[i].plot(lam,sp_data, label = 'x=%s\ny=%s'%(cab['X_FULL']+1,cab['Y_FULL']+1))
    ax[i,0].plot(lam,sp_data, label = 'Y%s (Ks = %.2f)'%(i+1,ls_spec[i][-1]), lw = 1, color = colorines[i])
    ax[i,0].tick_params(left = False, right = False , labelleft = False , 
                labelbottom = False, bottom = False) 
    dic_yso['Y%s'%(i+1)] = np.array([lam,sp_data])
    if i == 3:
        ax[i,0].set_ylabel('Normalized flux', fontsize = 25)
    
    for l in lines:
        ax[i,0].axvline(l, color = 'grey', alpha = 0.5, ls = 'dashed') 
    ax[i,0].legend(loc = 2, fontsize = 15)
ax[i,0].tick_params(left = False, right = False , labelleft = False , 
              labelbottom = True, bottom = True, labelsize = 20) 
ax[i,0].set_xlabel('$\lambda (\mu m)$', fontsize = 25)    
ax2 = ax[0,0].twiny()
ax2.set_xticks(lines)
ax2.set_xlim(lim_d, lim_u)
tl = l_names
ax2.set_xticklabels(tl, fontsize = 12)
ax2.plot(lam,sp_data, alpha = 0)

# plt.savefig(pruebas + 'spectra_young.png', dpi =300, bbox_inches = 'tight')


# %%
# In this section we are going to plot the spectra of some Late type star for 
# showing a comparation with the YSO

# fig, ax = plt.subplots(len(ls_spec), 1, figsize = (15,18))
# fig.subplots_adjust(hspace=0)
late_t = np.loadtxt(pruebas + 'late_type_for_comparation.txt')
late_t[-1] = np.array([10,1,14.02350000])
# lim_d, lim_u = 2, COIII +0.01

dic_late = {}
for i,y in enumerate(late_t):
    ax[i,1].set_xlim(2, COIII +0.01)
    print(y)
    all_fol = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/cluster_spectra/ifu_%.0f/half_%.0f/'%(reductions[0],y[0],y[1])
    all_fol = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/cluster_spectra/ifu_%.0f/half_%.0f/'%('ABC',y[0],y[1])

    all_spec = np.loadtxt(all_fol + 'gns_lib_in_kmos_ifu%.0f_half%.0f.txt'%(y[0],y[1]))
    # print(all_spec[:,4])
    id_mag = np.where(all_spec[:,4] == y[2])
    x_l, y_l = all_spec[id_mag][0][-1],all_spec[id_mag][0][-2]
    print(x_l, y_l)
    spec_late = fits.open(all_fol + 'spec_ifu%.0f_half%.0f_%.0f_%.0f.fits'%(y[0],y[1],x_l,y_l))
    cab_late = spec_late[0].header
    sp_late = spec_late[0].data
    lam_late = np.array([cab_late['CRVAL1']+cab_late['CDELT1']*i for i in range(cab_late['NAXIS1'])])
    
    good = np.where((lam_late > lim_d) & (lam_late < lim_u))
    norm =  np.where((lam_late > d_norm) & (lam_late < u_norm))
    filtered_data = sigma_clip(sp_late[norm], sigma=2, maxiters=None,cenfunc=mean, masked=False, copy=False)
    sp_late, lam_late  = sp_late[good], lam_late[good]
    sp_late = sp_late/np.nanmean(filtered_data)
    
    lam_selc_late= np.where((lam_late > norm_0) &(lam_late<norm_1))
    flat_flux_late = np.nanmean(sp_data[lam_selc_late])
    print(flat_flux_late)
    ax[i,1].set_ylim(0, 2)
    
    ax[i,1].plot(lam_late, sp_late, label = 'L%s (Ks = %.2f)'%(i+1,late_t[i][-1]),lw =1 )
    ax[i,1].legend(loc = 2,fontsize = 15)
    ax[i,1].tick_params(left = False, right = False , labelleft = False , 
                labelbottom = False, bottom = False) 
    dic_late['L%s'%(i+1)] = np.array([lam_late,sp_late])
    for l in lines:
        ax[i,1].axvline(l, color = 'grey', alpha = 0.5, ls = 'dashed') 
        
ax[i,1].tick_params(left = False, right = False , labelleft = False , 
              labelbottom = True, bottom = True, labelsize = 20) 
ax[i,1].set_xlabel('$\lambda (\mu m)$', fontsize = 25)    
ax2 = ax[0,1].twiny()
ax2.set_xticks(lines)
ax2.set_xlim(lim_d, lim_u)
tl = l_names
ax2.set_xticklabels(tl, fontsize = 12)
ax2.plot(lam,sp_data, alpha = 0)
plt.subplots_adjust(wspace=0.04)

# plt.savefig(pruebas + 'spectra_young_and_late.png', dpi =300, bbox_inches = 'tight')
# sys.exit()

# %%

fig, ax = plt.subplots(len(ls_spec), 1, figsize = (10,18))
fig.subplots_adjust(hspace=0)
lam = dic_yso['Y1'][0]
for j in range(len(ls_spec)):
    ax[j].set_xlim(lim_d, lim_u)
    lam_selc = np.where((lam > norm_0) &(lam<norm_1))
    
    sp_data = dic_yso['Y%s'%(j+1)][1]
    sp_late = dic_late['L%s'%(j+1)][1]
    good = np.where((lam > lim_d) & (lam < lim_u))
    norm =  np.where((lam > d_norm) & (lam < u_norm))
    filtered = sigma_clip(sp_data[norm], sigma=1.5, maxiters=None,cenfunc=mean, masked=False, copy=False)
    sp_data, lam  = sp_data[good], lam[good]
    sp_data = sp_data/np.mean(filtered)
    
   
    flat_flux = np.nanmean(sp_data[lam_selc])
    # print(flat_flux)
    # ax[i].set_ylim(flat_flux*0.5, flat_flux*2)
    ax[j].set_ylim(0, 2.2)
    ax[j].plot(lam,sp_data +0.2, label = 'M%s (Ks = %.2f)'%(j+1,ls_spec[j][-1]), lw = 1, color = colorines[j])
    ax[j].plot(lam,sp_late -0.4, label = 'L%s (Ks = %.2f)'%(j+1,late_t[j][-1]), lw = 1, color = '#1f77b4', alpha = 0.5)

   
    ax[j].tick_params(left = False, right = False , labelleft = False , 
                labelbottom = False, bottom = False) 
    if j == 3:
        ax[j].set_ylabel('Normalized flux + constant', fontsize = 25)
    
    for l in lines:
        ax[j].axvline(l, color = 'grey', alpha = 0.5, ls = 'dashed') 
    
    
    lg = ax[j].legend(loc = 9, fontsize = 15)
    for l in lg.legend_handles:
        l.set_linewidth(5)
ax[j].tick_params(left = False, right = False , labelleft = False , 
              labelbottom = True, bottom = True, labelsize = 20) 
ax[j].set_xlabel('$\lambda (\mu m)$', fontsize = 25)    
ax2 = ax[0].twiny()
ax2.set_xticks(lines)
ax2.set_xlim(lim_d, lim_u)
tl = l_names
ax2.set_xticklabels(tl, fontsize = 20)
ax2.plot(lam,sp_data, alpha = 0)


# plt.savefig(pruebas + 'spectra_young_and_late_2.png', dpi =300, bbox_inches = 'tight')

# %%
norm_fold = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/tramos_reduction/young_candidates/normalized_spline/'

dic_yso_norm = {}
for j in range(7):
    sp_late = dic_late['L%s'%(j+1)][1]
    
    name = glob.glob(norm_fold + 'yso%s*'%(j))
    print(name)
    spec = np.loadtxt(name[0])
    dic_yso_norm['Y%s'%(j+1)] = spec
    sp_data_norm = dic_yso_norm['Y%s'%(j+1)][:,1]
    lam_norm = dic_yso_norm['Y%s'%(j+1)][:,0]
fig, ax = plt.subplots(len(ls_spec), 1, figsize = (10,18))
fig.subplots_adjust(hspace=0)
lam = dic_yso['Y1'][0]
for j in range(len(ls_spec)):
    ax[j].set_xlim(lim_d, lim_u)
    lam_selc = np.where((lam > norm_0) &(lam<norm_1))
    
    sp_data = dic_yso_norm['Y%s'%(j+1)][:,1]
    lam_norm = dic_yso_norm['Y%s'%(j+1)][:,0]
    sp_late = dic_late['L%s'%(j+1)][1]
    

    ax[j].set_ylim(0, 2.2)
    ax[j].plot(lam_norm,sp_data +0.2, label = 'M%s (Ks = %.2f)'%(j+1,ls_spec[j][-1]), lw = 1, color = colorines[j])
    ax[j].plot(lam,sp_late -0.4, label = 'L%s (Ks = %.2f)'%(j+1,late_t[j][-1]), lw = 1, color = '#1f77b4', alpha = 0.5)

   
    ax[j].tick_params(left = False, right = False , labelleft = False , 
                labelbottom = False, bottom = False) 
    if j == 3:
        ax[j].set_ylabel('Normalized flux + constant', fontsize = 25)
    
    for l in lines:
        ax[j].axvline(l, color = 'grey', alpha = 0.5, ls = 'dashed') 
    
    
    lg = ax[j].legend(loc = 9, fontsize = 15)
    for l in lg.legend_handles:
        l.set_linewidth(5)
ax[j].tick_params(left = False, right = False , labelleft = False , 
              labelbottom = True, bottom = True, labelsize = 20) 
ax[j].set_xlabel('$\lambda (\mu m)$', fontsize = 25)    
ax2 = ax[0].twiny()
ax2.set_xticks(lines)
ax2.set_xlim(lim_d, lim_u)
tl = l_names
ax2.set_xticklabels(tl, fontsize = 20)
ax2.plot(lam_norm,sp_data_norm, alpha = 0)






































