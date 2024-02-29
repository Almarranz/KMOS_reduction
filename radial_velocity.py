#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 13:30:02 2024

@author: amartinez
"""

# Estimation of the Radial velocity from Brg lines
import sys
import glob
from subprocess import call
import astropy.units as u
from astropy import constants as const
from astropy.utils.data import download_file
from astropy.io import fits  # We use fits to open the actual data file
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation
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
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit 
from astropy.time import Time
import copy
from astropy.constants import c
from astropy.coordinates import ICRS, LSR
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
rcParams.update({'font.size': 10})
rcParams.update({'figure.figsize':(10,10)})
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

pruebas = '/Users/amartinez/Desktop/PhD/KMOS/practice/'
yso_coor = np.loadtxt(pruebas + 'yso_coord.txt')
# reduction = 'ABC'
reduction = 'tramos'
colorines = ['#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2','#7f7f7f', '#bcbd22', '#17becf', '#67AFAD', '#B018CC', '#D24FBD', '#A12322', '#4CC351', '#54DF4F', '#7389D2', '#898EE0', '#289C88', '#18EAA4', '#9ECC27', '#71A317', '#421256', '#A23C97', '#44302F']
# yso = 4
# yso_ls = [1,2,4,5]
yso_ls = [1]
age = 'young_candidates'
reductions = [reduction]
ls_spec = np.loadtxt(pruebas + 'ls_young.txt')
dic_yso ={} 

norm_fold = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/tramos_reduction/young_candidates/normalized_spline/'
dic_yso_norm = {}

fig, ax = plt.subplots(2,1)
fig.subplots_adjust(hspace=0)
brg_line = 2.166120
delta = 0.005
# for i in range(yso,yso+1):
for i in yso_ls:
    #Raw spectrum (not normalized)
    yso_par = ls_spec[i]
    # name = fits.open('/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/tramos_reduction/young_candidates/ifu_6/half_1/spec_ifu6_half1_61_26.fits')
    name = fits.open('/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/tramos_reduction/young_candidates/ifu_%.0f/half_%.0f/spec_ifu%.0f_half%.0f_%.0f_%.0f.fits'%(yso_par[0],yso_par[1],yso_par[0],yso_par[1],yso_par[2],yso_par[3]))

    cab = name[0].header
    lam = np.array([cab['CRVAL1']+cab['CDELT1']*i for i in range(cab['NAXIS1'])])
    
    sp_data = np.c_[lam, name[0].data]
    
    #Already normalized spectrum
    # name = glob.glob(norm_fold + 'yso%s*'%(i))
    # sp_data = np.loadtxt(name[0])
    
    br_g = (sp_data[:,0] >brg_line-delta) & (sp_data[:,0] <brg_line+delta)
    sp_data_brg = sp_data[br_g]
    ax[0].plot(sp_data_brg[:,0], sp_data_brg[:,1], color = colorines[i], label = 'M%s'%(i+1))
    ax[0].axvline(brg_line, color = 'k', alpha = 0.5, label = 'Br$\gamma$ = %.4f'%(brg_line))
    ax[1].axvline(brg_line, color = 'k', alpha = 0.5)
    ax[0].scatter(sp_data_brg[:,0], sp_data_brg[:,1], color = colorines[i])
    ax[0].legend()
    # ax.set_ylim(0.3,1.2)

no_lines = []   
select_wave = [] 
smo = 5e-38
# smo = 0.1
k_spl = 3
rv =[]
def gaussian(x, amplitude, mean, stddev):
    return 1 + amplitude * np.exp(-((x - mean) / stddev) ** 2 / 2)


def vertical(event):
    global x2, y2, chunk, no_lines, continuum_inds, vel  # Use global variables
    if event.dblclick and event.button ==1:
        no_lines.append(event.xdata)
        ax[0].axvline(event.xdata, ls = 'dashed',color = 'r' )
        if len(no_lines) == 2:
            ind1 = np.where(abs(sp_data_brg[:,0]-no_lines[0])== min(abs(sp_data_brg[:,0]-no_lines[0])))
            ind2 = np.where(abs(sp_data_brg[:,0]-no_lines[1])== min(abs(sp_data_brg[:,0]-no_lines[1])))
            
            ax[0].axvline(sp_data_brg[:,0][ind1], ls = 'dashed',color = 'green' )
            ax[0].axvline(sp_data_brg[:,0][ind2], ls = 'dashed',color = 'green' )
            # print(extract_spec)
            ax[0].axvspan(sp_data_brg[:,0][ind1][0], sp_data_brg[:,0][ind2][0],color = 'green', alpha = 0.2)
            # ax2.axvspan(no_lines[0], no_lines[1], color = 'r', alpha = 0.2)
            no_lines = []
           
            select_wave.append((sp_data_brg[:,0][ind1][0],sp_data_brg[:,0][ind2][0]))
    if event.dblclick and event.button ==3:
        # print('yomamma')
        lam, flux = sp_data_brg[:,0],sp_data_brg[:,1]
        sel_x = np.concatenate([lam[(lam >= start) & (lam <= end)] for start, end in select_wave])
        sel_y = np.concatenate([flux[(lam >= start) & (lam <= end)] for start, end in select_wave])
        no_sel_x = (lam > select_wave[0][1]) & (lam<select_wave[1][0])
        no_fit = lam[no_sel_x]
        
        normfunc = UnivariateSpline(sel_x, sel_y, s=smo, k = k_spl)
        normfunction = normfunc(sp_data_brg)
        # Gaussian fit
        
        pi = [np.min(flux/normfunction[:,0]),
              sp_data_brg[np.where(flux/normfunction[:,0] == np.min(flux/normfunction[:,0]))][0][0],
              np.std(sp_data_brg[:,0])]  # Amplitude, mean, standard deviation initial guess
        params, covariance = curve_fit(gaussian, lam[no_sel_x], flux[no_sel_x]/normfunction[:,0][no_sel_x], p0=pi)  
        print('params',params)
        
        ax[1].plot(lam,flux/normfunction[:,0], label ='A =%.2f\n$\mu$=%.4f\n$\sigma$=%.2f'%(pi[0]-1,pi[1],pi[2]))
        ax[1].legend()
        ax[0].plot(sp_data_brg, normfunction[:,0],'k')
        
        ax[0].set_xlim(np.min(select_wave)-0.001,np.max(select_wave)+0.001)
        # ax.set_ylim(np.min(flux)-np.min(flux)*0.1,np.max(flux)+np.max(flux)*0.1)
        
        # ax[1].plot(sp_data_brg,flux/normfunction[:,0])
        # ax[1].plot(sp_data_brg[:,0], sp_data_brg[:,1]+1, color = colorines[i])        
        ax[1].set_xlim(np.min(select_wave)-0.001,np.max(select_wave)+0.001)
        # for ns in range(len(no_fit)):
        #     ax[0].axvline(no_fit[ns], color = 'k', ls ='dotted', alpha = 0.3)
        #     ax[1].axvline(no_fit[ns], color = 'k', ls ='dotted',alpha = 0.3)
        
        pi = [min(flux[no_sel_x]/normfunction[:,0][no_sel_x])-1,
              no_fit[np.where(flux[no_sel_x]==np.min(flux[no_sel_x]))[0][0]],
               np.std(no_fit)]  # Amplitude, mean, standard deviation initial guess
        print(pi)
        # params, covariance = curve_fit(gaussian, no_fit, flux[no_sel_x]/normfunction[:,0][no_sel_x], p0 = pi) 
        # ym = gaussian(no_fit,params[0],params[1],params[2])
        # ax[1].plot(no_fit,ym, 'k', label = 'A$_{fit}$ =%.2f\n$\mu_{fit}$ = %.4f\n$\sigma_{fit}$ = %.4f'%(params[0],params[1],params[2]))
 
        params, covariance = curve_fit(gaussian, lam,flux/normfunction[:,0], p0 = pi) 
        ym = gaussian(lam,params[0],params[1],params[2])
        ax[1].plot(lam,ym, 'k', label = 'A$_{fit}$ =%.2f\n$\mu_{fit}$ = %.4f\n$\sigma_{fit}$ = %.4f'%(params[0],params[1],params[2]))

        ax[1].legend()
        vel = (params[1]- brg_line)/brg_line*const.c.to(u.km/u.s)
        print(vel)
        rv.append(vel)
        
        
        ax[0].axvline(params[1],color = 'k', alpha = 0.5, ls = 'dotted')
        ax[0].axvline(pi[1],color = 'k', alpha = 0.5, ls = 'dotted')
        ax[1].axvline(params[1],color ='k', alpha = 0.5, ls = 'dotted')
        
        loc = EarthLocation.of_site('Cerro Paranal')  
        t = Time('2021-06-05T00:00:00', scale='utc')
        sc = SkyCoord(yso_coor[i][0]*u.deg, yso_coor[i][1]*u.deg)
        vcorr = sc.radial_velocity_correction(kind='barycentric', obstime=t, location=loc)
        # rv_corr = vel -  vcorr
        rv_corr = vel -  0

        my_observation = ICRS(ra=yso_coor[i][0]*u.deg, dec=yso_coor[i][1]*u.deg, \
             pm_ra_cosdec=-0*u.mas/u.yr, pm_dec=-0*u.mas/u.yr, \
             radial_velocity=rv_corr, distance = 8.2e3*u.pc)

        new_rv = my_observation.transform_to(LSR()).radial_velocity
        print('new_rv',new_rv)
        
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        texto = 'v = %.2f km/s [$(\Delta \lambda/\lambda)*c$]\nv$_{corr}$ = %.2f km/s'%(vel.value,new_rv.value)
        ax[1].text(0.6, 0.3, texto, transform=ax[1].transAxes, fontsize=14,
        verticalalignment='top', bbox=props)


cidv = fig.canvas.mpl_connect('button_press_event',vertical)

# %%
# Correction of the radial velocity


# loc = EarthLocation.of_site('Cerro Paranal')  
# t = Time('2021-06-05T00:00:00', scale='utc')
# sc = SkyCoord(yso_coor[0]*u.deg, yso_coor[1]*u.deg)
# vcorr = sc.radial_velocity_correction(kind='barycentric', obstime=t, location=loc)
# rv_corr = rv[0] -  vcorr

# my_observation = ICRS(ra=yso_coor[i][0]*u.deg, dec=yso_coor[i][1]*u.deg, \
#      pm_ra_cosdec=-0*u.mas/u.yr, pm_dec=-0*u.mas/u.yr, \
#      radial_velocity=rv_corr, distance = 8.2e3*u.pc)

# new_rv = my_observation.transform_to(LSR()).radial_velocity
# print(new_rv)



