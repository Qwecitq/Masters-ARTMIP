#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 14:17:00 2022

@author: ktquagra
"""

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


path = '/global/project/projectdirs/m1517/cascade/external_datasets/ARTMIP/tier1/ftp_mirror'
flds = ['/cascade_bard_v1/','/guan_waliser/','/reid250/','/mundhenk/']
fpath = [ path+pt for pt in flds] 

#call the AR california data 
data_12_14 = collections.defaultdict(list)

cwp = os.getcwd()

#path = 'AR_California_new/'
vns = [x[1:-1] for x in flds ] 
for ix,vx in enumerate(vns):
    data_12_14[vx].append(xr.open_dataset(f'{cwp}/{path1}/'+vx+'1979_2015.nc',engine='hd5netcdf'))
    
    
###############################################################################################
#                                                             Functions for Saving Data by Chunks                                                                 #
###############################################################################################


def split_by_chunks(dataset):
    chunk_slices = {}
    for dim, chunks in dataset.chunks.items():
        slices = []
        start = 0
        for chunk in chunks:
            if start >= dataset.sizes[dim]:
                break
            stop = start + chunk
            slices.append(slice(start, stop))
            start = stop
        chunk_slices[dim] = slices
    for slices in itertools.product(*chunk_slices.values()):
        selection = dict(zip(chunk_slices.keys(), slices))
        yield dataset[selection]


def create_filepath(ds, prefix='filename', root_path="."):
    """
    Generate a filepath when given an xarray dataset
    """
    start = str(ds.time.data[0])[:10]
    end = str(ds.time.data[0])[:10]
    filepath = f'{root_path}{prefix}_{start}_{end}.nc'
    return filepath



###############################################################################################
#                                                                            Load Data                                                                                             #
###############################################################################################

print('Begin loading data')

viwvn=[]

first_date = list(os.listdir(viwvn_path[1:]))
tx = [int(x) for x in first_date] ; tx.sort()
#vbl_nm = 'e5.oper.an.vinteg.162_072_viwvn.ll025sc.'

for i in tx:
    viwvn.append(xr.open_dataset(glob.glob(viwvn_path[1:]+str(tx)+'*/'+vbl_nm+'*.nc',recursive=True)[0],engine='h5netcdf'))
print('Done loading data \n starting to concatenate data')


#concatenate the list of data loaded
viwvn = xr.concat(viwvn, dim='time') 

###############################################################################################
#                                        Fixing the longitudes from 0 to 360 to -180 to 180                                                                    #
###############################################################################################


print('Fixing Longitudes from -180 to 180')

viwvn_data = viwvn
viwvn_data.coords['longitude'] = (viwvn_data.coords['longitude']  + 180) % 360 - 180
viwvn_data=viwvn_data.sortby(viwvn.longitude)
viwvn_data=viwvn_data.sel(latitude=slice(49,15),longitude = slice(-145,-60))
print('Done \n Begin composites of AR days for various algorithms')



viwvn_new = collections.defaultdict(list)

###############################################################################################
#                                                                      Find Composites                                                                                         #
###############################################################################################
     
for ix, x in enumerate(flds):
    data_in = data_12_14[x[1:-1]][0]#.resample(time='6H').mean() #resample data into 6H
    
    data_time = data_in.time.values #call out the time from the AR datasets 
    data_in = data_in.load()
    #composite times of clwc and ARDTs 
    viwvn_new[x[1:-1]].append(viwvn_data.sel(time=np.in1d(viwvn_data.time.values,data_time)))                                                        


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
    if os.path.isdir(direc) == False:
        os.mkdir(direc)
        os.mkdir(dirs)
        t=ch_paths(viwvn_new[f][0], dirs , prefix)
































