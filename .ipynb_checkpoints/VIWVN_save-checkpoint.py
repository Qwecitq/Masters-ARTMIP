#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 14:17:00 2022

@author: ktquagra
"""


'''This script is for selecting the consecutive days that ARs are captured in all the other datasets [MSLP,CLWC,TCVW,IVT and others]'''

#import packages and modules 
import os
import xarray as xr
import collections
import glob
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
from cartopy import feature as cf
import pandas as pd
xr.set_options(display_style="text")
import itertools
import numpy as np
import sys


'''
arg 1 : path to the data that has ARs at the specified region 
    
arg 2 : short name of the variable 
    
arg 3 : path to the new variable to be compared
    
arg 4 : directory to store files 
    
arg 5 : prefix for the filenames 

'''
path1 = sys.argv[1]  #path to the data that has ARs at the specified region 
vbl_nme=sys.argv[2]     #short name of the variable 
viwvn_path = sys.argv[3]  #path to the new variable to be compared 
direc = sys.argv[4]  #directory to store files 
prefix = sys.argv[5]  #prefix for the filenames 

###############################################################################################
#                                                                       AR composited dataset                                                                              #
###############################################################################################


path = '/global/project/projectdirs/m1517/cascade/external_datasets/ARTMIP/tier1/ftp_mirror'  #path to the main AR datasets 
flds = ['/cascade_bard_v1/','/guan_waliser/','/mundhenk/','/reid250/'] #the algorithms under cosideration 
fpath = [ path+pt for pt in flds]  #the file paths 
 
#call the AR california data 
data_12_14 = collections.defaultdict(list) #all stored AR data 

cwp = os.getcwd()  #get the present working directory 

#load the AR datasets and keep them all in one dictionary 
vns = [x[1:-1] for x in flds ] 
for ix,vx in enumerate(vns):
    data_12_14[vx].append(xr.open_mfdataset(f'{cwp}/{path1}{flds[ix]}'+vx+'_newlandfall_ARs__0.6-0.4_*.nc',engine='netcdf4'))
    
    
###############################################################################################
#                                                             Functions for Saving Data by Chunks                                                                 #
###############################################################################################


def split_by_chunks(dataset):
    chunk_slices = {}
    for dim, chunks in dataset.chunks.items(): #get the chunks from the data loaded by dask 
        slices = [] #a list of the slices 
        start = 0 #start to count the chucks 
        for chunk in chunks: #specific chunks 
            if start >= dataset.sizes[dim]:  #if the number of start >= size of the dataset which should not be, break
                break
            stop = start + chunk #set the end point for the chunk 
            slices.append(slice(start, stop)) #create a slice object from the chuck start to stop and append to slices 
            start = stop #reset the start point to the end point of the previous chunk's end point
        chunk_slices[dim] = slices #set the dictionary chunk_slices to a specific time chunk 
    for slices in itertools.product(*chunk_slices.values()): #use the slices created together with the itertools.product to call the chunck slices 
        selection = dict(zip(chunk_slices.keys(), slices)) #create a dict and a zip of the slices 
        yield dataset[selection]#apply the selection zipped dict on the data to select the chunks 


def create_filepath(ds, prefix='filename', root_path="."):
    """
    Generate a filepath when given an xarray dataset
    """
    start = str(ds.time.data[0])[:10]
    end = str(ds.time.data[0])[:10]
    filepath = f'{root_path}{prefix}_{start}_{end}.nc'
    return filepath



###############################################################################################
#                                                                            Load Atmospheric Varibles Data                                                           #
###############################################################################################

print('Begin loading data')

viwvn=[]  #list of appended data from the other variables 

first_date = list(os.listdir(viwvn_path)) #get the first time 
tx = [int(x) for x in first_date] ; tx.sort() #create a list of the dates and sort them 

#load the data 
for i in tx:
    viwvn.append(glob.glob(viwvn_path+str(tx)+'*/'+vbl_nme+'*.nc',recursive=True)) #append the data to a list 
print('Done loading data \n starting to concatenate data')


#concatenate the list of data loaded
print(viwvn[0][0]) #check the name of the first file in the list 
print(viwvn[1][0]) #check if they are repeated 
print(viwvn[0][0]==viwvn[1][0]) #check for the repitition. If yes, then continue 

viwvn = viwvn[0]
print('Beginning to load the multiple datasets')

viwvn = xr.open_mfdataset(viwvn,engine='h5netcdf',parallel = True)  #load the data for the specific variable 

print('Done loading data')
###############################################################################################
#                                        Fixing the longitudes from 0 to 360 to -180 to 180                                                                    #
###############################################################################################


print('Fixing Longitudes from -180 to 180')

viwvn_data = viwvn
viwvn_data.coords['longitude'] = (viwvn_data.coords['longitude']  + 180) % 360 - 180 #convert from 0-360 to -180 to 180
viwvn_data=viwvn_data.sortby(viwvn.longitude) #sort the lons 
viwvn_data=viwvn_data.sel(longitude = slice(-170,-30), latitude = slice(60,10)) #select the region under consideration 
print('Done \n Begin composites of AR days for various algorithms')



viwvn_new = collections.defaultdict(list) #create a dictionary of the AR composite data for the various variables 

###############################################################################################
#                                                                      Find Composites                                                                                         #
###############################################################################################
     
for ix, x in enumerate(flds):
    print(f'Working on data : {x[1:-1]}')
    data_in = data_12_14[x[1:-1]][0]  #resample data into 6H
    
    data_time = data_in.time.values #call out the time from the AR datasets 
    data_in = data_in.load()  #load the entire AR data into memory for the selection 
    print(f'Done loading data for composites')
    #composite times of clwc and ARDTs 
    print(data_time)
    print(viwvn_data)
    viwvn_new[x[1:-1]].append(viwvn_data.sel(time=np.array(data_time)))                                                        
    print(f'Working on data : {x[1:-1]}')

###############################################################################################
#                                                                       Saving Data                                                                                               #
###############################################################################################

print(f'Done \n Start Saving data into {direc}')
VIWVE_collection = collections.defaultdict(list) 

def ch_paths(data,loc,pfix):
    in_data = list(split_by_chunks(data))
    paths = [create_filepath(ds,prefix=pfix, root_path=loc) for ds in in_data]
    t=xr.save_mfdataset(in_data,paths)
    return t

for fx,f in enumerate(flds):
    f=f[1:-1]
    dirs = f'{direc}{f}_{prefix}/'
    if os.path.isdir(direc) == True:
        print(f'Found directory {direc}')
        os.mkdir(dirs)
        
        print(f'Started Saving Data : {f}')
        t=ch_paths(viwvn_new[f][0], dirs , prefix)
        print(f'Done Saving Data : {f}')
        
    elif os.path.isdir(direc) == False:
        print(f'Found no directory {direc} \n Creating directory {direc}')
        os.mkdir(direc)
        os.mkdir(dirs)
        
        print(f'Started Saving Data : {f}')
        t=ch_paths(viwvn_new[f][0], dirs , prefix)
        print(f'Done Saving Data : {f}')

