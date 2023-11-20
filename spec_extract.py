
"""
Created on Tue Nov  7 09:57:30 2023

@author: amartinez
"""

# Extract the sky and the stars sky-subtracted
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



# %%
half_ifu = 2
ifu_sel =12
spec_folder = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/cluster_spectra/ifu_%s/half_%s/'%(reduction, ifu_sel, half_ifu)
aligned_cube = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/cubes/cube_ifu%s_half%s_mean.fits'%(reduction,ifu_sel,half_ifu)

delete_spec = 'no'#!!!

if delete_spec == 'yes':
    del_ques = input('Do you really want to delete files form ifu%s half%s?'%(ifu_sel, half_ifu))
    if del_ques == 'yes':
        shutil.rmtree(spec_folder + '/' )
        os.mkdir(spec_folder)
    elif del_ques == 'no':
        sys.exit(84)

dic_x = {}
dic_y = {}
for ifu in range(1,24):
    x_shifts = []
    y_shifts = []
    with open(log_5 + 'esorex.log', 'r') as f:
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

# %%

# ima = fits.open(esorex_cube_5)


ima = fits.open(aligned_cube)
mapa = WCS(ima[1].header).celestial

h0 = ima[0].header
h1 = ima[1].header

# This is the header for the spectrum fits file
# h_esp = h1.copy()
model = '/Users/amartinez/Desktop/PhD/KMOS/practice/spec_header_model.fits'
h_esp = fits.open(model)[0].header
# h_esp['CRVAL1'] = h_esp['CRVAL3']
# h_esp['CDELT1'] = h_esp['CDELT3']
# h_esp['CTYPE1'] = 'WAVE'
# h_esp['CUNIT1'] = 'um'
# del h_esp['CRVAL3'], h_esp['CRVAL2'], h_esp['CDELT2'], h_esp['CDELT3']

ima_1 = fits.open(esorex_ima_5)
wcs = WCS(ima[1].header)
h0_1 = ima_1[0].header
h1_1 = ima_1[1].header

cube = ima[1].data

# cube = cube[:,y_d:y_up, x_d:x_up]
cube_im = ima_1[1].data[y_d:y_up, x_d:x_up]

sky_cal = input('Choose: Sci or Sky?')#TODO
# sky_cal = 'Sky' #TODO
if sky_cal == 'Sci':
    if os.path.isfile(spec_folder  + 'sky_ifu%s_half%s.fits'%(ifu_sel,half_ifu)):
        print('There is difuse emission sprectrum!')
    else:
        print('There is no sky compted for ifu = %s, half = %s'%(ifu_sel,half_ifu))
        sys.exit('DO IT!')

fig, ax = plt.subplots(1,1,figsize =(10,10)) 
ax.set_title('%s, IFU %s, Half %s '%(sky_cal, ifu_sel, half_ifu))
#Uncomment these for the WCS projection 
# =============================================================================
# ax.axis("off")
# ax = plt.subplot(projection=wcs,slices=(1000, 'y', 'x'))
# =============================================================================
v_min =0.01e-18
v_max = 0.01e-16


# im = ax.imshow(cube_im, cmap='inferno', origin='lower',vmin=v_min,vmax = v_max,alpha =1)
cube_plot = np.nanmedian(cube, axis = 0)
im = ax.imshow(cube_plot, cmap='inferno',label='overlays',origin='lower',vmin=v_min,vmax = v_max,alpha =1)
fig.colorbar(im, ax=ax, orientation='vertical')

# x, y = 62,43
# yy,xx = np.indices(cube_im.shape)
# distances = np.sqrt((xx - x)**2 + (yy - y)**2)
# mask = np.where(distances <=3)
# values = cube[:,mask[0],mask[1]]
# spec = np.mean(values, axis = 1)
# ax.scatter(mask[1],mask[0],color = 'r')

x = None
y = None
wide = 4#TODO
w = int(wide/2)
def extract_spec(cube,y,x):
    # yy,xx = np.indices(cube_im.shape)
    xx,yy = np.indices(cube_im.shape)
    distances = np.sqrt((xx - x)**2 + (yy - y)**2)
    mask = np.where(distances <=w)
    spec = cube[:,mask[0],mask[1]]
    spec_mean = np.mean(spec, axis = 1)
    ax.scatter(mask[1],mask[0],color = 'r')
    if sky_cal == 'Sci':
        hdu_sky = fits.open(spec_folder  + 'sky_ifu%s_half%s.fits'%(ifu_sel,half_ifu))
        sky_data = hdu_sky[0].data
        spec_mean_sky = spec_mean - sky_data
        # spec = cube[:,x-w:x+w+1,y-w:y+w+1]
        # spec_mean = np.mean(spec, axis =(1,2)) 
        # print('SEPC',spec_mean[1200],spec_mean_1[1200]) 
        hdu_spec = fits.PrimaryHDU(data=spec_mean_sky, header=h_esp)
        hdu_spec.writeto(spec_folder + 'spec_ifu%s_half%s_%s_%s.fits'%(ifu_sel,half_ifu,x,y), overwrite= True)
    
        return spec_mean_sky
    
    
    else:
        return spec_mean

# Function to update the plot with the clicked point
def update_plot(x,y):
    x = int(np.rint(x))
    y = int(np.rint(y))
    plt.scatter(x, y, color='red', marker='x', s=100)  # Customize the marker style
    # square = Rectangle((x - wide/2-0.5, y - wide/2-0.5), wide, wide, linewidth=1, edgecolor='r', facecolor='none')
    circle = Circle((x, y),w, facecolor = 'none', edgecolor = 'r', linewidth=1)
    # ax.add_patch(square)
    ax.add_patch(circle)
    ax.text(x-w,y-w-1,'%.0f,%.0f'%(x,y), color = 'r',fontsize = 15)
    plt.draw()  # Redraw the plot with the updated point


cab = ima[1].header
lam = np.array([cab['CRVAL3']+cab['CDELT3']*i for i in range(cab['NAXIS3'])])
d_spec = lam[0]
u_spec = lam[-1]
good = np.where((lam > d_spec) & (lam < u_spec))
lam = lam[good]


fig2, ax2 = plt.subplots(1,1,figsize =(20,10))
HeI = 2.058
COI =  2.29322
COII = 2.32246
Brg = 2.165
He = 2.12
HeII = 2.189
# H2 = 2.12
l_names = ['HeI', 'COI', 'COII','Br$\gamma$', 'He', 'HeII']
lines = [HeI, COI, COII,Brg, He, HeII]
for l in lines:
    ax2.axvline(l, color = 'grey', alpha = 0.5, ls = 'dashed')
ax2b = ax2.twiny()
ax2b.set_xticks(lines)
tl = l_names
ax2b.set_xticklabels(tl)

def plot_spec(flux):
        
    flux = flux[good]
    sig = np.nanstd(flux)    
    # ax2.set_ylim(np.nanmean(flux)-sig, np.nanmean(flux) + sig)
    ax2.plot(lam,flux, label = '(%.0f,%.0f)'%(x,y))
    ax2.set_xlabel('$\lambda (\mu m)$')
    ax2.set_ylabel('flux')
    ax2.legend(fontsize=10)
    
    
    ax2b.plot(lam,flux, alpha = 0)
    
    
def bright(y,x):
    x = int(np.rint(x))
    y = int(np.rint(y))
    around = cube_im[y-w:y+w,x-w:x+w]
    bright_x, bright_y = np.unravel_index(around.argmax(), around.shape)
    bright_flux = np.max(around)
    if sky_cal == 'Sci':
        x = x-(w-bright_y)
        y = y-(w-bright_x)
        print( bright_flux,  bright_x,bright_y)
        print('max flus at:x = %s, y = %s'%(x-(w-bright_y),y-(w-bright_x)))
    elif sky_cal == 'Sky':
        x = int(np.rint(x))
        y = int(np.rint(y))
        
    return x,y 




def onclick(event):
    global x, y  # Use global variables
    if event.xdata is not None and event.ydata is not None:
        x = event.xdata
        y = event.ydata 
        x = int(np.rint(x))
        y = int(np.rint(y))
        x,y = bright(y,x)
        flux = extract_spec(cube,int(np.rint(x)), int(np.rint(y)))
        update_plot(x,y)
        plot_spec(flux)
        
        sig = np.nanstd(flux)/2
        print(f'Clicked at x={x}, y={y}')
np_spec = np.empty((5,3072))
clicks = 0 
def onclick_sky(event):  
        global x, y, clicks# Use global variables
        if event.xdata is not None and event.ydata is not None:
            x = event.xdata
            y = event.ydata 
            x = int(np.rint(x))
            y = int(np.rint(y))
            x,y = bright(y,x)
            flux = extract_spec(cube,int(np.rint(x)), int(np.rint(y)))
            update_plot(x,y)
            plot_spec(flux)            
            sig = np.nanstd(flux)/2
            clicks += 1
            np_spec[clicks-1] = flux 
            print(f'Clicked at x={x}, y={y}, clicks = {clicks}')
            # print(dic_spec)
            if clicks == 5:
                sky_mean = np.mean(np_spec,axis = 0)
                hdu_sky = fits.PrimaryHDU(data=sky_mean, header=h_esp)
                hdu_sky.writeto(spec_folder + 'sky_ifu%s_half%s.fits'%(ifu_sel,half_ifu), overwrite= True)
                print('Saved mean sky')
                plt.close()
                plt.close()
                # sys.exit(293)

def delete(event):
    if event.key == 'd':
        ax2.cla()  
        ax2b.cla()
        for l in lines:
            ax2.axvline(l, color = 'grey', alpha = 0.5, ls = 'dashed')
        ax2.set_xticks(lines)
        tl = l_names
        ax2.set_xticklabels(tl)   
    
    
    
if sky_cal == 'Sci':
    cid = fig.canvas.mpl_connect('button_press_event',onclick)
if sky_cal == 'Sky':
    cid = fig.canvas.mpl_connect('button_press_event',onclick_sky)
cid_del = fig2.canvas.mpl_connect('key_press_event',delete)

# %%
# import matplotlib.pyplot as plt
# import numpy as np

# # Create some sample data
# x = np.linspace(0, 10, 100)
# y = np.sin(x)

# # Create the first plot
# fig, ax1 = plt.subplots()
# ax1.plot(x, y, 'b-')
# ax1.set_xlabel('Primary X-axis')
# ax1.set_ylabel('Y-axis', color='b')
# ax1.tick_params('x', colors='b')

# # Create the secondary x-axis
# ax2 = ax1.twiny()
# ax2.set_xlabel('Secondary X-axis', color='r')
# ax2.tick_params('x', colors='r')

# plt.title('Primary and Secondary X-axes')
# plt.show()


    






