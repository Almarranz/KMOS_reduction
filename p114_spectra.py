#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 10:26:23 2024

@author: amartinez
"""

# Analysis of the P114 spectra

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
IPython.get_ipython().run_line_magic('matplotlib', 'auto')
# IPython.get_ipython().run_line_magic('matplotlib', 'inline')

cab_folder = '/Users/amartinez/Desktop/PhD/KMOS/P114_spectra/'
cab_p114 = np.loadtxt(cab_folder +  'P114-Arches-B1K-r2-2014-08-05_spec.dat')
cab_p114_2 = np.loadtxt(cab_folder +  'P114-B1K-2014-05-02_spec.dat')
cab_p114_3 = np.loadtxt(cab_folder +  'P114-Arches-B1K-r1-2014-08-05_spec.dat')

# fig, ax = plt.subplots(1,1)
# ax.plot(spec[:,0], spec[:,1])


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
He = 2.112
HeII = 2.189
N3 = 2.115
C4 = 2.078

# H2 = 2.12
l_names = ['HeI', '$^{12}$CO(2,0)', ' $^{12}$CO(3,1)\n','Br$\gamma$', 'He\n', 'HeII','$^{12}$CO(4,2)','NIII','CIV']
lines = [HeI, COI, COII,Brg, He, HeII,COIII, N3, C4]

# l_names = [ '$^{12}$CO(2,0)', ' $^{12}$CO(3,1)\n','Br$\gamma$', 'He','$^{12}$CO(4,2)']
# lines = [COI, COII,Brg, He,COIII]

age = 'young_candidates'
reductions = [reduction]
ls_spec = np.loadtxt(pruebas + 'ls_young.txt')
# ls_spec = np.delete(ls_spec,1, axis =0)
# colorines = ['#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2','#7f7f7f','#8c564b', '#e377c2','#7f7f7f']
colorines = ['#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2','#7f7f7f', '#bcbd22', '#17becf', '#67AFAD', '#B018CC', '#D24FBD', '#A12322', '#4CC351', '#54DF4F', '#7389D2', '#898EE0', '#289C88', '#18EAA4', '#9ECC27', '#71A317', '#421256', '#A23C97', '#44302F']

fig, ax = plt.subplots(1, 1, figsize = (30,18))
fig.suptitle('%s'%(reduction))
ax2 = ax.twiny()
# fig.subplots_adjust(hspace=0)

# lim_d, lim_u = 2.15, COIII +0.01
lim_d, lim_u = 2, COIII +0.12#!!!

norm_0, norm_1 = HeI, He
d_norm, u_norm =  2.2,2.28 #!!!
snr_d, snr_u = 2.2356, 2.2661#!!!
dic_yso = {}
dic_yso_snr = {}
# for i in range(len(ls_spec)):
ls_spec[0] = [7,2,13,50,10.84]



for i in range(1):

    ax.set_xlim(lim_d, lim_u)
    ifu_sel = ls_spec[i][0]
    half_ifu = ls_spec[i][1]
    x, y = ls_spec[i][2], ls_spec[i][3]
    spec_fol = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/%s/ifu_%.0f/half_%.0f/'%(reductions[0],age,ifu_sel,half_ifu)
    spec = fits.open(spec_fol + 'spec_ifu%.0f_half%.0f_%.0f_%.0f.fits'%(ifu_sel, half_ifu, x, y))
    cab = spec[0].header
    sp_data = spec[0].data
    lam = np.array([cab['CRVAL1']+cab['CDELT1']*i for i in range(cab['NAXIS1'])])
    # sys.exit()
    dic_yso_snr['Y%s'%(i+1)] = np.array([lam,sp_data])
    lam_selc = np.where((lam > norm_0) &(lam<norm_1))
    
    good = np.where((lam > lim_d) & (lam < lim_u))
    norm =  np.where((lam > d_norm) & (lam < u_norm))
    filtered = sigma_clip(sp_data[norm], sigma=1.5, maxiters=None,cenfunc=mean, masked=False, copy=False)
    # sp_data, lam  = sp_data[good], lam[good]
    sp_data = sp_data/np.mean(filtered)
    
    
    
    flat_flux = np.nanmean(sp_data[lam_selc])
    # print(flat_flux)
    # ax[i].set_ylim(flat_flux*0.5, flat_flux*2)
    # ax.set_ylim(0, 2)
     # ax[i].plot(lam,sp_data, label = 'x=%s\ny=%s'%(cab['X_FULL']+1,cab['Y_FULL']+1))
    ax.plot(lam,sp_data - i, label = 'Y%s (Ks = %.2f)'%(i+1,ls_spec[i][-1]), lw = 1, color = colorines[i])
    ax.plot(cab_p114[:,0], cab_p114[:,1]-0.5, label = 'CAB P114\n(Ago2)')
    ax.plot(cab_p114_2[:,0], cab_p114_2[:,1]-1.5, label = 'CAB P114\n(May)', color = 'g')
    ax.plot(cab_p114_3[:,0], cab_p114_3[:,1]-1, label = 'CAB P114\n(Ago1)', color = 'r')
    ax.tick_params(left = False, right = False , labelleft = False , 
                labelbottom = False, bottom = False) 
    dic_yso['Y%s'%(i+1)] = np.array([lam,sp_data])
    
    ax.set_ylabel('Normalized flux', fontsize = 25)
    
    for l in lines:
        ax.axvline(l, color = 'grey', alpha = 0.1, ls = 'dashed', zorder = 0) 
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 10 )
ax.tick_params(left = False, right = False , labelleft = False , 
              labelbottom = True, bottom = True, labelsize = 20) 
ax.set_xlabel('$\lambda (\mu m)$', fontsize = 25)    
# ax2 = ax.twiny()
ax2.set_xticks(lines)
ax2.set_xlim(lim_d, lim_u)
tl = l_names
ax2.set_xticklabels(tl, fontsize = 12)
ax2.plot(lam,sp_data, alpha = 0)
ax.set_ylim(-1*2,2)

chunk = 0     
no_lines = []   
continuum = []
continuum_inds = []
chunk_ind = []
delete_inds = []
deselect = []
def vertical(event):
    global x2, y2, chunk, no_lines, continuum_inds  # Use global variables
    if event.dblclick and event.button ==1:
        no_lines.append(event.xdata)
        ax2.axvline(event.xdata, ls = 'dashed',color = 'r' )
        if len(no_lines) == 2:
            ind1 = np.where(abs(lam-no_lines[0])== min(abs(lam-no_lines[0])))
            ind2 = np.where(abs(lam-no_lines[1])== min(abs(lam-no_lines[1])))
            
            ax2.axvline(lam[ind1], ls = 'dashed',color = 'green' )
            ax2.axvline(lam[ind2], ls = 'dashed',color = 'green' )
            # print(extract_spec)
            continuum.append(lam[ind1[0][0]:ind2[0][0]+1])
            continuum_inds.append(np.arange(ind1[0][0], ind2[0][0]+1,1))
            chunk_ind.append(np.full((ind2[0][0]+1-ind1[0][0]),chunk))
            ax2.axvspan(lam[ind1][0], lam[ind2][0],color = 'green', alpha = 0.2)
            # ax2.axvspan(no_lines[0], no_lines[1], color = 'r', alpha = 0.2)
            no_lines = []
            chunk +=1
            print('chunk_ind', chunk_ind)
            print('%.4f %.4f'%(lam[ind1][0],lam[ind2][0]))
            deselect.append((lam[ind1][0],lam[ind2][0]))
            np.savetxt(pruebas + 'lines_ifu%.0f_half%.0f_%.0f_%.0f.txt'%(ifu_sel, half_ifu, x, y),deselect, fmt = '%.4f')
    if event.dblclick and event.button ==3:
            del continuum[-1]
            del continuum_inds[-1]
            del chunk_ind[-1]
            no_lines = []
def snr_cal(event):
    props = dict(boxstyle='round', facecolor='white', alpha=1 )
    if event.key == 'c':
        cont_id = np.sort(np.concatenate(continuum_inds))   
        for i in range(len(dic_yso_snr)):
            # snr = np.nanmedian(dic_yso_snr['Y%s'%(i+1)][1][cont_id])/np.nanstd(dic_yso_snr['Y%s'%(i+1)][1][cont_id])
            snr = np.nanmedian(dic_yso['Y%s'%(i+1)][1][cont_id])/np.nanstd(dic_yso['Y%s'%(i+1)][1][cont_id])

            print('SNR, i = %.2f, %s'%(snr, i+1))
            # ax.text(2.3,1-i,'Y%s SNR = %.2f'%(i+1, snr))
            ax.text(2.01,1-i,'Y%s SNR = %.2f'%(i+1, snr), fontsize=8, verticalalignment='top', bbox=props)

def re_plot(event):
    if event.key == 'd':
        del continuum_inds[:]
        ax.cla()
        ax2.cla()
        for i in range(len(ls_spec)):
            ax.set_xlim(lim_d, lim_u)
            ifu_sel = ls_spec[i][0]
            half_ifu = ls_spec[i][1]
            x, y = ls_spec[i][2], ls_spec[i][3]
            spec_fol = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/%s/ifu_%.0f/half_%.0f/'%(reductions[0],age,ifu_sel,half_ifu)
            spec = fits.open(spec_fol + 'spec_ifu%.0f_half%.0f_%.0f_%.0f.fits'%(ifu_sel, half_ifu, x, y))
            cab = spec[0].header
            sp_data = spec[0].data
            lam = np.array([cab['CRVAL1']+cab['CDELT1']*i for i in range(cab['NAXIS1'])])
            lam_selc = np.where((lam > norm_0) &(lam<norm_1))
            
            good = np.where((lam > lim_d) & (lam < lim_u))
            norm =  np.where((lam > d_norm) & (lam < u_norm))
            filtered = sigma_clip(sp_data[norm], sigma=1.5, maxiters=None,cenfunc=mean, masked=False, copy=False)
            sp_data, lam  = sp_data[good], lam[good]
            sp_data = sp_data/np.mean(filtered)
            
            
            
            flat_flux = np.nanmean(sp_data[lam_selc])
            # print(flat_flux)
            # ax[i].set_ylim(flat_flux*0.5, flat_flux*2)
            # ax.set_ylim(0, 2)
             # ax[i].plot(lam,sp_data, label = 'x=%s\ny=%s'%(cab['X_FULL']+1,cab['Y_FULL']+1))
            ax.plot(lam,sp_data - i, label = 'Y%s (Ks = %.2f)'%(i+1,ls_spec[i][-1]), lw = 1, color = colorines[i])
            ax.tick_params(left = False, right = False , labelleft = False , 
                        labelbottom = False, bottom = False) 
            
            ax.set_ylabel('Normalized flux', fontsize = 25)
            
            for l in lines:
                ax.axvline(l, color = 'grey', alpha = 0.5, ls = 'dashed') 
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 10 )
        ax.tick_params(left = False, right = False , labelleft = False , 
                      labelbottom = True, bottom = True, labelsize = 20) 
        ax.set_xlabel('$\lambda (\mu m)$', fontsize = 25)    
        ax2.set_xticks(lines)
        ax2.set_xlim(lim_d, lim_u)
        tl = l_names
        ax2.set_xticklabels(tl, fontsize = 12)
        ax2.plot(lam,sp_data, alpha = 0)
        ax.set_ylim(-1*len(ls_spec),2)
        
        
cidv = fig.canvas.mpl_connect('button_press_event',vertical)
cid_snr = fig.canvas.mpl_connect('key_press_event',snr_cal)
cid_re = fig.canvas.mpl_connect('key_press_event',re_plot)
















