#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import xarray as xr
import collections
import glob
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
from cartopy import feature as cf
import pandas as pd
#xr.set_options(display_style="text")
import sys
import numpy as np

print(f'Welcome to the AR data selection criteria 2')



'''Define the following 


arguement 1 : path to get the AR algorithm data from
arguement 2 : destination folder for saving files
arguement 3 : threshold to select the main landfall region 
arguement 4 : threshold for other regions 
arguement 5 : left lon
arguement 6 : right lon 
arguement 7 : bottom lat 
arguement 8 : top lat 

'''

#All the user input arguments have been put here.
data_path = sys.argv[1]   #where find  the AR data 
save_path = sys.argv[2] #where to store the data that is computed 
th1 = float(sys.argv[3]) ; th2 = float(sys.argv[4]) #thresholds 1 for the main box and 2 for the secondary boxes C and B 
llon = int(sys.argv[5]) ; rlon = int(sys.argv[6]) #longitudes and latitudes for the boxes 
blat = int(sys.argv[7]) ; tlat = int(sys.argv[8]) 

print(f'Setting Boundary conditions for the main land fall region {llon,rlon,blat,tlat}')

A_lon = [llon,rlon+1]; A_lats = [tlat-1,tlat+1]  #box A boundary 
B_lon = [llon,rlon]; B_lats = [blat,tlat]  #box B boundary 
C_lon = [llon,rlon+1]; C_lats = [blat-1,blat+1]  #box C boundary 

print('Done')

#path = '/global/project/projectdirs/m1517/cascade/external_datasets/ARTMIP/tier1/ftp_mirror'
flds = ['/cascade_bard_v1/','/guan_waliser/','/reid250/','/mundhenk/'] #algorithm names 
#fpath = [ path+pt for pt in flds] 
fpath = [ data_path+pt for pt in flds]   #create the combination of data paths and algorithm names to find data 


print('Loading Data') 
data = collections.defaultdict(list) #create a collection to hold data 

#open the AR datasets and store them before computing selecting the region 
for idx, dx in enumerate(fpath):
    file = [x for x in list(glob.glob(dx + "*.nc4", recursive = True))] #get the collection of data using a wildcard method 

    dat = xr.open_mfdataset(file,chunks={'time':1000},parallel = True).sel(lon = slice(-170,-30), lat = slice(10,60),time=slice('1979','2021')) #open data 
    data[flds[idx][1:-1]].append(dat) #store data 
    
print('Done')    
    

print(f'Selection based on criteria {th1} for main box and {th2} for upper and lower boxes')
#select the data based on the criteria 
for vbls in list(data.keys()):
    Tsp =[] #this is to store the datasets after the regions are selected 
    vble = 'ar_binary_tag'
    reg_data = data[vbls][0][vble]   #regional data for various algorithms 
    
    #select the boundary boxes for the AR classification 
    #A_lon = [-134,-119]; A_lats = [41,43]
    #B_lon = [-134,-120]; B_lats = [35,42]
    #C_lon = [-134,-119]; C_lats = [34,36]
    print(f'Check for {vbls} conditions that fit the required thresholds {th1} and {th2}')
    for tx,tm in enumerate(reg_data.time.values) :
        
        main_box = reg_data[tx].sel(lon=slice(B_lon[0],B_lon[1]), lat = slice(B_lats[0],B_lats[1])).mean()  #.sel(time = str(tm))  #main boundary box to check for data 
        print('Main box selection done')
        ab_box = reg_data[tx].sel(lon=slice(A_lon[0],A_lon[1]), lat = slice(A_lats[0],A_lats[1])).mean() #boundary box above the main box 
        print(' Box B selection done')
        be_box =reg_data[tx].sel(lon=slice(C_lon[0],C_lon[1]), lat = slice(C_lats[0],C_lats[1])).mean()    #boundary box below the main box 
        print(' Box C selection done')
        
        print(f'Begin with the thresholding ....   {len(Tsp)}')
        if (main_box.values >th1 and ab_box.values<=th2) or (main_box.values >th1 and  be_box.values<=th2) or (main_box.values >th1 and  ab_box.values<=th2 and be_box.values<=th2) : #use the thresholds to select the data and store the data in Tsp 
            Tsp.append(reg_data[tx])
            print(f'Done with time : {tm}')
    
    print(f'Begin concatenation for {vbls}')
    concat_data = xr.concat(Tsp,dim = 'time') #concatenate the data
    
    print(f'Saving netcdf data to {save_path}') 
    concat_data.to_netcdf(f'{save_path}{vbls}_landfall_ARs__{th1}-{th2}.nc') #save the data 
    print(f'Done Saving for {vbls}')

print('All processes completed successfully')