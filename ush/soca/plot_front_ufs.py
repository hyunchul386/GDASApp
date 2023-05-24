#!/usr/bin/env python3
import cartopy.crs as ccrs
#import cartopy.feature as cfeature
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import netCDF4 as nc
from netCDF4 import Dataset, date2num, num2date
import numpy as np
#import argparse
import glob
#import os
import yaml
from datetime import datetime

def IJ_find(imm,jmm,lons,lats,Xupper,Xlower,Yupper,Ylower):
    for j in range(jmm):
        if (abs(lats[j,ic]-Ylower) < little):
            jst = j 
        if (abs(lats[j,ic]-Yupper) < little):
            jed = j 
    for i in range(imm):
        if (abs(lons[jst,i]-Xlower) < little):
            ist = i 
        if (abs(lons[jed,i]-Xupper) < little):
            ied = i 
    return ist,ied+1,jst,jed+1

def Find_zfactor(kmm,z,ref_z):
    kmmm = kmm - 1 
    for k in range(kmmm):
        if ( z[k] <= ref_z ) and ( z[k+1] > ref_z ):
            dz = z[k+1] - z[k]
            dz1 = ref_z - z[k]
            z_fac = dz1/dz
            ref_k = k
    return z_fac, ref_k

def Find_fields(im,jm,T,ssh,ist,ied,jst,jed,z_fac, ref_k):
    TC=np.zeros((jm, im,))
    SSH=np.zeros((jm, im,))
   #SSH[jst:jed,ist:ied] = ssh[0,jst:jed,ist:ied]  
    for i in range(ist,ied+1):
        for j in range(jst,jed+1):
            TC[j,i] = T[0,ref_k,j,i]-(T[0,ref_k,j,i]-T[0,ref_k+1,j,i])*z_fac 
            SSH[j,i] = ssh[0,j,i]
    return TC, SSH

def plot_world_map(im,jm,km,ic,jc,little,ref_z,ref_T,lons, lats, z, T, ssh, ymdh, RG):
    imm = im - 1 ; jmm = jm - 1 ; kmm = km - 1
    if ( RG == 1 ):
    #-- Gulf Stream only
        Xupper = 320 ; Xlower = 278  ; Yupper = 45 ; Ylower = 26
        outfile='./front_output_GS_usf.png'
        RGname='Gulf Stream'
    elif (RG == 2):
    #-- Kuroshio only
        Xupper = 180 ; Xlower = 130 ; Yupper = 46 ; Ylower = 25
        outfile='./front_output_KS_ufs.png'
        RGname='Kuroshio'
    
#  #--- Xlower/Xupper [ lons : -280 ~ 80 ]
    if (float(Xlower) < -280):
        Xlower = float(Xlower) + 360
    elif (float(Xlower) > 80):
        Xlower = float(Xlower) - 360
    else :
        Xlower = float(Xlower)

    if (float(Xupper) < -280):
        Xupper = float(Xupper) + 360
    elif (float(Xupper) > 80):
        Xupper = float(Xupper) - 360
    else :
        Xupper = float(Xupper)

    fig = plt.figure(figsize=(12,8))
    cenlon = 0.0
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=cenlon))
    ax.set_extent([float(Xlower), float(Xupper), float(Ylower), float(Yupper)])
#   ax.coastlines(resolution='10m',zorder=2);
    ax.coastlines();
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.7, linestyle='--')
    gl.xlocator = mticker.MultipleLocator(5)
    gl.ylocator = mticker.MultipleLocator(5)
#   ax.add_feature(cfeature.COASTLINE)
#   ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
#------
    ist,ied,jst,jed = IJ_find(imm,jmm,lons,lats,Xupper,Xlower,Yupper,Ylower)
#   print("ist,ied,jst,jed",ist,ied,jst,jed)
    z_fac, ref_k = Find_zfactor(kmm,z,ref_z)
    TC,SSH = Find_fields(im,jm,T,ssh,ist,ied,jst,jed,z_fac, ref_k)
#------
  # cs = plt.contourf(lons[:,:],lats[:,:],ssh[0,:,:],cmap='jet')
    cs = plt.contourf(lons[:,:],lats[:,:],SSH[:,:],cmap='jet')
    cb = plt.colorbar(cs, orientation='horizontal', shrink=0.5, pad=.04)
    cs = ax.contour(lons[:,:],lats[:,:],TC[:,:], [ref_T], color='darkblue', linewidths=2.5)
   
    plttitle = '%s Front and SSH (m) on %s' % (RGname, ymdh)
 
    plt.title(plttitle)
    plt.savefig(outfile)
    plt.close()

    return()

def read_var(datapath):
#   obsfiles = glob.glob(datapath+'*')
    datanc = nc.Dataset(datapath)
    lats =datanc.variables['geolat']
    lons = datanc.variables['geolon']
    z = datanc.variables['z_l']
#   T = datanc.variables['Temp']
    T = datanc.variables['temp']
#   ssh = datanc.variables['ave_ssh']
    ssh = datanc.variables['SSH']
    time = datanc.variables['time']
    dates = num2date(time, datanc.variables['time'].units)
    ymdh = datetime.strptime(str(dates[0]),'%Y-%m-%d %H:%M:%S')
    return lons, lats, z, T, ssh, ymdh

def gen_figure(inpath,im,jm,km,ic,jc,little,ref_z,ref_T):
   #read the files to get the 2D array to plot
    lons, lats, z, T, ssh, ymdh = read_var(inpath)
    #-- for Gulf Steam
    plot_world_map(im,jm,km,ic,jc,little,ref_z,ref_T,lons, lats, z, T, ssh, ymdh, 1)
    #-- for Kuroshio
#   plot_world_map(im,jm,km,ic,jc,little,ref_z,ref_T,lons, lats, z, T, ssh, ymdh, 2)


if __name__ == "__main__":

   input = "./ufs_input.nc"

   inpyaml = open("./testinput/plot_front.yaml", 'r')
   inp = yaml.load(inpyaml, Loader=yaml.FullLoader)
   im = inp["im"]
   jm = inp["jm"]
   km = inp["km"]
   ic = inp["ic"]
   jc = inp["jc"]
   little = inp["little"]
   ref_z = inp["ref_z"]
   ref_T = inp["ref_T"]



   gen_figure(input,im,jm,km,ic,jc,little,ref_z,ref_T)

