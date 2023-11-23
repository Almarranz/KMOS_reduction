
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
from matplotlib.widgets import Slider
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
ifu_sel = 7
spec_folder = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/cluster_spectra/ifu_%s/half_%s/'%(reduction, ifu_sel, half_ifu)
spec_young = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/young_candidates/ifu_%s/half_%s/'%(reduction, ifu_sel, half_ifu)
aligned_cube = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/cubes/cube_ifu%s_half%s_mean.fits'%(reduction,ifu_sel,half_ifu)

delete_spec = 'no'#!!!

if delete_spec == 'yes':
    del_ques = input('Do you really want to delete files form ifu%s half%s?'%(ifu_sel, half_ifu))
    if del_ques == 'yes':
        shutil.rmtree(spec_folder + '/' )
        os.mkdir(spec_folder)
    elif del_ques == 'no':
        # sys.exit(84)
        print('Not deleting')
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

# sky_cal = input('Choose: Sci or Sky?')#TODO
sky_cal = 'Sci' #TODO
if sky_cal == 'Sci':
    if os.path.isfile(spec_folder  + 'sky_ifu%s_half%s.fits'%(ifu_sel,half_ifu)):
        difuse_flux = fits.open(spec_folder  + 'sky_ifu%s_half%s.fits'%(ifu_sel,half_ifu))[0].data
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



# im = ax.imshow(cube_im, cmap='inferno', origin='lower',vmin=v_min,vmax = v_max,alpha =1)
cube_plot = np.nanmedian(cube, axis = 0)
v_min =np.nanmin(cube_plot)
v_max =np.nanmax(cube_plot)
im = ax.imshow(cube_plot, cmap='inferno',label='overlays',origin='lower',vmin=v_min,vmax = v_max,alpha =1)
fig.colorbar(im, ax=ax, orientation='vertical')

# Add a slider for vmin
ax_vmin = plt.axes([0.2, 0.1, 0.3, 0.03], facecolor='lightgoldenrodyellow')
vmin_slider = Slider(ax_vmin, 'vmin', v_min, v_max, valinit=v_min)

# Add a slider for vmax
ax_vmax = plt.axes([0.2, 0.05, 0.3, 0.03], facecolor='lightgoldenrodyellow')
vmax_slider = Slider(ax_vmax, 'vmax', v_min, v_max, valinit=v_max)

def update(val):
    vmin = vmin_slider.val
    vmax = vmax_slider.val
    im.set_clim(vmin, vmax)
    plt.draw()

vmin_slider.on_changed(update)
vmax_slider.on_changed(update)

plt.draw()


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

# event2 = 0
# def w_around(event):
#     if event.key is not None and event.key.isdigit() is True:
#         print('You pushed',event.key)
#         w2 = int(event.key)
#         return w2
# def ring(cube, y, x):

w1 = 2
w2 = 2

coor = []
def extract_spec(cube,y,x):
    global temp
    coor.append([y,x])
    if len(coor)>1 and coor[-2] == coor[-1]:
        temp.remove()
    # yy,xx = np.indices(cube_im.shape)
    xx,yy = np.indices(cube_plot.shape)
    distances = np.sqrt((xx - x)**2 + (yy - y)**2)
    mask = np.where(distances <=w)
    spec = cube[:,mask[0],mask[1]]
    spec_mean = np.mean(spec, axis = 1)
    ax.scatter(mask[1],mask[0],color = 'lime',s = 50)
    
    mask_around = np.where((distances>w+w1)&(distances<w+w1+w2))
    r_ind = np.random.choice(len(mask_around[0]), size=30, replace=False)
    mask_rand = tuple(array[r_ind] for array in mask_around)
    # rand_anillos.append(mask_rand)
    # ring = cube[:,mask_around[0],mask_around[1]]
    ring = cube[:, mask_rand[0], mask_rand[1]]
    print('mask_around',mask_around[0])
    print('mask_around',mask_around[1])
    ring_mean = np.mean(ring, axis =1)
    # ax.scatter(mask_around[1],mask_around[0], color = 'orange', s = 50)
    temp = ax.scatter(mask_rand[1],mask_rand[0], color = 'green', s = 20)
    # temp =ax.scatter(rand_anillos[-1][1],rand_anillos[-1][0], color = 'red', s = 20)
    
    plt.draw()
    # ax2.plot(lam,ring_mean, color = 'fuchsia')
    
    if sky_cal == 'Sci':
        # hdu_sky = fits.open(spec_folder  + 'sky_ifu%s_half%s.fits'%(ifu_sel,half_ifu))
        # sky_data = hdu_sky[0].data
        # spec_mean_sky = spec_mean - sky_data
        
        # h_esp = fits.open(model)[0].header
        # hdu_spec = fits.PrimaryHDU(data=spec_mean_sky, header=h_esp)
        # hdu_spec.writeto(spec_folder + 'spec_ifu%s_half%s_%s_%s.fits'%(ifu_sel,half_ifu,x,y), overwrite= True)
        spec_mean_sky = spec_mean-ring_mean
        # return spec_mean_sky
        return spec_mean
    
    
    else:
        return spec_mean
    
# rand_anillos = []
# def anillos():
#     xx,yy = np.indices(cube_plot.shape)
#     distances = np.sqrt((xx - x)**2 + (yy - y)**2)
    
#     mask_around = np.where((distances>w+w1)&(distances<w+w1+w2))
#     r_ind = np.random.choice(len(mask_around[0]), size=30, replace=False)
#     mask_rand = tuple(array[r_ind] for array in mask_around)
#     rand_anillos.append(mask_rand)
#     # ring = cube[:,mask_around[0],mask_around[1]]
#     ring = cube[:, mask_rand[0], mask_rand[1]]
#     print('mask_around',mask_around[0])
#     print('mask_around',mask_around[1])
#     ring_mean = np.mean(ring, axis =1)
#     # ax.scatter(mask_around[1],mask_around[0], color = 'orange', s = 50)
#     ax.scatter(mask_rand[1],mask_rand[0], color = 'green', s = 20)
#     for anillos in rand_anillos:
#         ax.scatter(anillos[1],anillos[0], color = 'green', s = 20)
        
#     plt.draw() 



# Function to update the plot with the clicked point
def update_plot(x,y):   
    color_upd = 'lime'
    x = int(np.rint(x))
    y = int(np.rint(y))
    plt.scatter(x, y, color= color_upd, marker='x', s=50)  # Customize the marker style
    circle = Circle((x, y),w, facecolor = 'none', edgecolor = color_upd, linewidth=1)
    # ax.add_patch(square)
    ax.add_patch(circle)
    ax.text(x-w,y-w-1,'%.0f,%.0f'%(x,y), color = color_upd,fontsize = 15)
     # Redraw the plot with the updated point


cab = ima[1].header
lam = np.array([cab['CRVAL3']+cab['CDELT3']*i for i in range(cab['NAXIS3'])])
d_spec = lam[0]
u_spec = lam[-1]
good = np.where((lam >= d_spec) & (lam <= u_spec))
lam = lam[good]

def onkey(event):
    if event.key == 'd':
        if rand_anillos:
            rand_anillos.pop()  # Remove the last clicked point
            update_plot()


fig2, ax2 = plt.subplots(1,1,figsize =(20,10))
props = dict(boxstyle='round', facecolor='grey', alpha=0.5)
ax2.text(0.05, 0.95, 'Sci:\n"d":clear plot\n"y":save young star\nSky:\n"d":delete last\n"g": save', transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
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



for l in lines:
    ax2.axvline(l, color = 'grey', alpha = 0.5, ls = 'dashed')
ax2b = ax2.twiny()
ax2b.set_xticks(lines)
tl = l_names
ax2b.set_xticklabels(tl)

# if sky_cal == 'Sci':
#     sig_w = 3
#     ind1 = np.where(abs(lam-He)== min(abs(lam-He)))
#     ind2 = np.where(abs(lam-Brg)== min(abs(lam-Brg)))
#     ax2b.plot(lam,difuse_flux[good], color ='k', alpha = 0.3,zorder = 1, label ='Difuse emission')
#     ax2b.axhline(np.nanmean(difuse_flux[ind1[0][0]:ind2[0][0]]), alpha = 0.3, color = 'k')
#     sig = np.nanstd(difuse_flux[ind1[0][0]:ind2[0][0]])
#     ax2.axhspan(np.nanmean(difuse_flux[ind1[0][0]:ind2[0][0]])- sig*sig_w,
#                 np.nanmean(difuse_flux[ind1[0][0]:ind2[0][0]])+ sig*sig_w, color = 'k', alpha = 0.1)

def plot_spec(flux):
    
    factor = 5e-18
    flux = flux[good]
    sig = np.nanstd(flux)    
    # ax2.set_ylim(np.nanmean(flux)-sig, np.nanmean(flux) + sig)
    if sky_cal == 'Sci':
        ax2.plot(lam,flux, label = '(%.0f,%.0f)'%(x,y))
            
    if sky_cal == 'Sky':
        ax2.plot(lam,flux+ clicks*factor, label = '(%.0f,%.0f)'%(x,y))
    ax2.set_xlabel('$\lambda (\mu m)$')
    ax2.set_ylabel('flux')
    ax2.legend(fontsize=10)
    
    
    ax2b.plot(lam,flux, alpha = 0)
    
    
def bright(y,x):
    x = int(np.rint(x))
    y = int(np.rint(y))
    around = cube_plot[y-w:y+w,x-w:x+w]
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
        # flux = extract_spec(cube,int(np.rint(x)), int(np.rint(y)))
        flux = extract_spec(cube,int(np.rint(x)), int(np.rint(y)))
        update_plot(x,y)
        plot_spec(flux)
        
        # sig = np.nanstd(flux)/2
        print(f'Clicked at x={x}, y={y}')
    if event.key:
        print('Your mamma',event.key)
        
        
np_spec = []
clicks = 0 
def onclick_sky(event):  
        global x, y, clicks, np_spec# Use global variables
        if event.xdata is not None and event.ydata is not None:
            x = event.xdata
            y = event.ydata 
            x = int(np.rint(x))
            y = int(np.rint(y))
            x,y = bright(y,x)
            flux = extract_spec(cube,int(np.rint(x)), int(np.rint(y)))
            update_plot(x,y)
            plot_spec(flux)            
            # sig = np.nanstd(flux)/2
            clicks += 1            
            print(f'Clicked at x={x}, y={y}, clicks = {clicks}')
            # print(dic_spec)
           
        
    
def delete(event):
    global np_spec
    if sky_cal == 'Sci':
        if event.key == 'd':
            ax2.cla()  
            ax2b.cla()
            for l in lines:
                ax2.axvline(l, color = 'grey', alpha = 0.5, ls = 'dashed')
            ax2.set_xticks(lines)
            tl = l_names
            ax2.set_xticklabels(tl)   
            ax2b.plot(lam,difuse_flux[good], color ='k', alpha = 0.3,zorder = 1, label ='Difuse emission')
            ax2b.axhline(np.nanmean(difuse_flux[ind1[0][0]:ind2[0][0]]), alpha = 0.3, color = 'k')
            sig = np.nanstd(difuse_flux[ind1[0][0]:ind2[0][0]])
            ax2.axhspan(np.nanmean(difuse_flux[ind1[0][0]:ind2[0][0]])- sig*sig_w,
                        np.nanmean(difuse_flux[ind1[0][0]:ind2[0][0]])+ sig*sig_w, color = 'k', alpha = 0.1)
        if event.key == 'y':
            flux_young = extract_spec(cube,int(np.rint(x)), int(np.rint(y)))
            h_esp = fits.open(model)[0].header
            h_esp.append(('X_full', x_d + x, '+1 on image'), end=True)
            h_esp.append(('Y_full', y_d + y, '+1 on image'), end=True)
            hdu_young = fits.PrimaryHDU(data = flux_young, header=h_esp)
            hdu_young.writeto(spec_young + 'spec_ifu%s_half%s_%s_%s.fits'%(ifu_sel,half_ifu,x,y), overwrite= True)
        
            
    if sky_cal == 'Sky':
        
        if event.key == 'g':
            flux = extract_spec(cube,int(np.rint(x)), int(np.rint(y)))
            print(f'Saved at x={x}, y={y}, clicks = {clicks}')
            np_spec.append(flux)
            print('Saved %s arrays'%(len(np_spec)))
        if event.key == 'd':
            ax2.cla()  
            ax2b.cla()
            for l in lines:
                ax2.axvline(l, color = 'grey', alpha = 0.5, ls = 'dashed')
            ax2.set_xticks(lines)
            tl = l_names
            ax2.set_xticklabels(tl)  
            if len(np_spec)>=1:
                del np_spec[-1]
        if len(np_spec)== 5:
            np_spec = np.array(np_spec)
            sky_mean = np.mean(np_spec,axis = 0)
            hdu_sky = fits.PrimaryHDU(data=sky_mean, header=h_esp)
            hdu_sky.writeto(spec_folder + 'sky_ifu%s_half%s.fits'%(ifu_sel,half_ifu), overwrite= True)
            print('Saved mean sky')
            plt.close()
            plt.close()  
    
    
    
if sky_cal == 'Sci':
    cid = fig.canvas.mpl_connect('button_press_event',onclick)
    cid_a = fig.canvas.mpl_connect('key_press_event', onkey)
if sky_cal == 'Sky':
    cid = fig.canvas.mpl_connect('button_press_event',onclick_sky)
cid_del = fig2.canvas.mpl_connect('key_press_event',delete)


# %%
# import numpy as np

# # Your original tuple of arrays
# array_tuple = (np.arange(1,61), np.arange(10,610,10))

# # Set a random seed for reproducibility
# np.random.seed(42)

# # Randomly select 10 indices
# selected_indices = np.random.choice(len(array_tuple[0]), size=10, replace=False)

# # Use the selected indices to extract elements from each array in the tuple
# selected_elements = tuple(array[selected_indices] for array in array_tuple)

# print(selected_elements)
# %%
# import matplotlib.pyplot as plt
# from matplotlib.widgets import Slider
# data = np.random.rand(10, 10)


# def w_around(event2):
#     if event2.key.isdigit():
#         print('You pushed',event2.key.isdigit())
#         print(int(event2.key)*5)
# # Create a figure and axis
# figA, ax = plt.subplots()
# plt.subplots_adjust(bottom=0.25)  # Adjust the bottom to make room for the slider

# # Create an initial image
# image = ax.imshow(data, cmap='viridis')
# cid_A = figA.canvas.mpl_connect('key_press_event',w_around)

# %%









