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
import csv

def Line_lon_lat(inmrf):
    with open(inmrf, newline='') as csvfile:
        fread = csv.reader(csvfile, delimiter='/')
        icheck = 0
        line_lons = []
        line_lats = []
        line_lon = []
        line_lat = []
        ck_lists = ["TEXT","ARC","ENDAT"]
        line_lists = ["N. WALL GULF STREAM","S. WALL GULF STREAM","N. WALL KUROSHIO","S. WALL KUROSHIO"]
        #--- line reading
        ch_line = 0
        n = 0
        for row in fread:
            #- N Gulf
            n = n + 1
         #  print (n, icheck, row[0])
            if (row[0] == 'TEXT') and any(line_list in [row[6]] for line_list in line_lists):
                icheck = 1
            elif (icheck == 1):
                if not any(ck in [row[0],row[1]] for ck in ck_lists):
                    if (row[0] == "LINE"):
                        ch_line = ch_line + 1
                    if (row[0] == "LINE") and (ch_line >= 2):
                        line_lons.append(line_lon)
                        line_lats.append(line_lat)
                        line_lon = []
                        line_lat = []
                    for xy in row[:]:
                        if ( len(xy) == 7 ):
                            xx = float(xy[0:5])/100.
                            if ( xx < 180 ):
                                xx = -1.0 * xx
                            else:
                                xx = 360.0 - xx
                            line_lon.append(xx)
                        elif ( len(xy) == 6 ):
                            yy = float(xy[0:4])/100.
                            line_lat.append(yy)
                else:
                    icheck = 0
                    ch_line = 0
                    line_lons.append(line_lon)
                    line_lats.append(line_lat)
                    line_lon = []
                    line_lat = []
                    
      # print(len(line_lons),type(line_lons))
      # print(len(line_lats),type(line_lats))
    return line_lons, line_lats

def IJ_find(imm,jmm,lon,lat,Xupper,Xlower,Yupper,Ylower):
    
    jst = 0; jed =0
    for j in range(jmm-1):
        if (abs(lat[j,ic]-Ylower) < abs(lat[jst,ic]-Ylower)):
            jst = j 
        if (abs(lat[j,ic]-Yupper) < abs(lat[jed,ic]-Yupper)):
            jed = j 
    ist = 0; ied =0
    for i in range(imm-1):
        if (abs(lon[jc,i]-Xlower) < abs(lon[jc,ist]-Xlower)):
            ist = i 
        if (abs(lon[jc,i]-Xupper) < abs(lon[jc,ied]-Xupper)):
            ied = i 
    return ist,ied+1,jst,jed+1

def IJ_find_ori(lons,lats,Xupper,Xlower,Yupper,Ylower):
    jst = 0; jed =0
    for j in range(jmm):
        if (abs(lats[j,ic]-Ylower) < little):
            jst = j 
        if (abs(lats[j,ic]-Yupper) < little):
            jed = j 
    for i in range(imm):
        if (abs(lons[jc,i]-Xlower) < little):
            ist = i 
        if (abs(lons[jc,i]-Xupper) < little):
            ied = i 
    return ist,ied,jst,jed

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
    SST=np.zeros((jm, im,))
    for i in range(ist,ied+1):
        for j in range(jst,jed+1):
            TC[j,i] = T[0,ref_k,j,i]-(T[0,ref_k,j,i]-T[0,ref_k+1,j,i])*z_fac 
            SSH[j,i] = ssh[0,j,i]
            SST[j,i] = T[0,0,j,i]
    return TC, SSH, SST
 
def plot_world_map(im,jm,km,ic,jc,little,ref_z,ref_T,lons, lats, z, T, ssh, ymdh, RG, SF, inmrf):
    imm = im - 1 ; jmm = jm - 1 ; kmm = km - 1

    if ( SF == 1 ):
        FD = "SSH"
    elif ( SF == 2 ):
        FD = "SST"

    if ( RG == 1 ):
    #-- Gulf Stream only
        Xupper = 320 ; Xlower = 278  ; Yupper = 45 ; Ylower = 26
        outfile='./front_output_GS_'+FD+'.png'
        RGname='Gulf Stream'
    elif (RG == 2):
    #-- Kuroshio only
        Xupper = -182.5 ; Xlower = -232.5 ; Yupper = 46 ; Ylower = 25
        outfile='./front_output_KS_'+FD+'.png'
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

    l_lons, l_lats = Line_lon_lat(inmrf)

    fig = plt.figure(figsize=(12,8))
    cenlon = 0.0
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=cenlon))
    ax.set_extent([float(Xlower), float(Xupper), float(Ylower), float(Yupper)])
    ax.coastlines() 
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.7, linestyle='--')
    gl.xlocator = mticker.MultipleLocator(5)
    gl.ylocator = mticker.MultipleLocator(5)
#------
    ist,ied,jst,jed = IJ_find(imm,jmm,lons,lats,Xupper,Xlower,Yupper,Ylower)
#   print("ist,ied,jst,jed: ",ist,ied,jst,jed)
    z_fac, ref_k = Find_zfactor(kmm,z,ref_z)
    TC,SSH,SST = Find_fields(im,jm,T,ssh,ist,ied,jst,jed,z_fac, ref_k)
#------
    if (SF == 1):
        min_lev = -2.0 ;  max_lev = 2.0 ;  step_lev = 0.2
        levels = np.arange(min_lev,max_lev+step_lev,step_lev)
        cs = plt.contourf(lons[:,:],lats[:,:],SSH[:,:],levels,cmap='jet')
        cb = plt.colorbar(cs, orientation='horizontal', shrink=0.5, pad=.04)
        cs = ax.contour(lons[:,:],lats[:,:],TC[:,:], [ref_T], color='darkblue', linewidths=2.5)
#       print (len(l_lons),range(len(l_lons)))
        l_lon = []
        l_lat = []
        for i in range(len(l_lons)): 
            l_lon = l_lons[i]
            l_lat = l_lats[i]
#           if (i == 0):
#               print(l_lon,l_lat)
            cs = plt.plot(l_lon,l_lat, color='darkred', linewidth=2.5)
        plttitle = '%s Front and SSH (m) on %s' % (RGname, ymdh)
    elif (SF == 2):
        min_lev = 10.0 ;  max_lev = 34.0 ;  step_lev = 1.0
        levels = np.arange(min_lev,max_lev+step_lev,step_lev)
        cs = plt.contourf(lons[:,:],lats[:,:],SST[:,:],levels,cmap='jet')
        cb = plt.colorbar(cs, orientation='horizontal', shrink=0.5, pad=.04)
        cs = ax.contour(lons[:,:],lats[:,:],TC[:,:], [ref_T], color='darkblue', linewidths=2.5)
        l_lon = []
        l_lat = []
        for i in range(len(l_lons)): 
            l_lon = l_lons[i]
            l_lat = l_lats[i]
#           if (i == 0):
#               print(l_lon,l_lat)
            cs = plt.plot(l_lon,l_lat, color='darkred', linewidth=2.5)
        plttitle = '%s Front and SST ($^{o}C$) on %s' % (RGname, ymdh)
 
    plt.title(plttitle)
    plt.savefig(outfile)
    plt.close()

    return()

def read_var(datapath):
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

def gen_figure(inpath,inmrf_gs,inmrf_np,im,jm,km,ic,jc,little,ref_z,ref_T):
   #read the files to get the 2D array to plot
    lons, lats, z, T, ssh, ymdh = read_var(inpath)
    #-- for Gulf Steam
    plot_world_map(im,jm,km,ic,jc,little,ref_z,ref_T,lons, lats, z, T, ssh, ymdh, 1, 1,inmrf_gs)
    plot_world_map(im,jm,km,ic,jc,little,ref_z,ref_T,lons, lats, z, T, ssh, ymdh, 1, 2,inmrf_gs)
    #-- for Kuroshio 
   #plot_world_map(im,jm,km,ic,jc,little,ref_z,ref_T,lons, lats, z, T, ssh, ymdh, 2, 1,inmrf_np)
   #plot_world_map(im,jm,km,ic,jc,little,ref_z,ref_T,lons, lats, z, T, ssh, ymdh, 2, 2,inmrf_np)


if __name__ == "__main__":

   input = "./ufs_input.nc"
   inmrf_gs = "./gs.mrf"
   inmrf_np = "./np.mrf"

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

   gen_figure(input,inmrf_gs,inmrf_np,im,jm,km,ic,jc,little,ref_z,ref_T)

