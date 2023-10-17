#!/usr/bin/env python3
# -*- coding: utf-8 -*-

exec(open('imports.py').read())
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import svd
import saver #personal created module to save file in chunks

# Load and preprocess climate data (sample data)
# Replace this with your own climate data loading and preprocessing.
# Ensure that the data is organized in a grid, with time as rows and spatial locations as columns.

# Generate sample data (replace with your own data)
#data = np.random.rand(100, 50)  # 100 time steps, 50 spatial locations
#data = anom_data
# Normalize the data by subtracting the mean along each spatial location.
mpath = '/pscratch/sd/k/kquagra/ARs Work/Non-coinciding ARDTs/'
climpath = '/N/project/ktquagra_project/'
reg = 'west_europe'

fld = sys.argv[1] #ARDT under consideration
vb = sys.argv[2] #variable under consideration 
which_type = sys.argv[3] #type of transect
n_components = 5# sys.argv[4] 
#n_components = int(n_components)# Number of EOFs to calculate

#fld ='cascade_bard_v1' #ARDT under consideration
#vb = 'q' #variable under consideration 
#which_type = 'EW' #type of transect
#n_components = 4 #sys.argv[4] 
#n_components = int(n_components)# Number of EOFs to calculate

if vb == 'pot':
    VB = 'PT'
else:
    VB = vb.upper()


 #receive input for which type of transect is done 

#change from zonal to meridional (lon to lat) based on the type of transect, either NS or EW
if which_type == 'NS':
    tloc = 'latitude'
elif which_type =='EW':
    tloc = 'longitude'
    
#print('Loading climatology ...')
#clim_files = glob.glob(f'{mpath}Bash Runs/{VB}_climatologies_west_europe/*.nc',recursive=True)
#clim = [xr.open_dataset(i , chunks={'time':1,'level':37}) for i in clim_files]

#print('Concating climatology for processing ... ')
#clim = xr.concat(clim,dim='time')
#clim.coords['longitude'] = (clim['longitude']+180)%360 - 180
#clim = clim.sortby(clim.longitude) ;  clim = clim.sel(longitude=slice(-15,30),latitude=slice(80,35))
#clim_dec = clim.sel(time=clim.time.dt.month.isin([1,2,12]))#.mean(tloc).mean('time')
#print(clim_dec)

print('Loading Data')
if fld == 'cascade_bard_v1':
    #calculate the anomaly of the data based on the type of transect 
    file_loc = glob.glob(f'{reg}_{vb}/{fld}_{vb}/*.nc',recursive=True)
    print(file_loc)
    data = xr.open_mfdataset(file_loc)
    data.coords['longitude'] = (data['longitude']+180)%360 - 180
    data = data.sortby(data.longitude); data = data.sel(longitude=slice(-15,30),latitude=slice(80,35))
    
    '''
    if vb=='pv':
        _, index = np.unique(clim_dec['time'], return_index=True) #remove duplicate times 
        clim = clim_dec.isel(time=index)
    else:
        clim = clim_dec'''
        
    data = data #- clim.mean(tloc).mean('time') #anomaly 
    print(data)
    
else:
    #obtain file location 
    file_loc =  glob.glob(f'{reg}_{vb}/{fld}_{vb}/*.nc', recursive=True)
    print(file_loc)
    #laod data and calculate anomaly based on the type of transect
    data = [xr.open_dataset(ds,chunks={'time':1,'level':37}) for ds in file_loc]
    data = xr.concat(data,dim = 'time')
    data.coords['longitude'] = (data['longitude']+180)%360 - 180
    data = data.sortby(data.longitude); data = data.sel(longitude=slice(-15,30),latitude=slice(80,35))
    
    '''
    if vb=='pv':
        _, index = np.unique(clim_dec['time'], return_index=True) #remove duplicate times 
        clim = clim_dec.isel(time=index)
    else:
        clim = clim_dec'''
    
    data = data# - clim.mean(tloc).mean('time') #anomaly 
    print(data)
    

print(f'Selecting Transects for {which_type}')
if which_type == 'across' :

    #use the metpy option for this kind of transect for along 2 different x,y points 
    start_point = (67,-15) #start point for cs
    end_point = (67,28) #end point for cs

    #select cross-section by squeezing data and selecting the cross section
    trans_data = data[VB].metpy.parse_cf().squeeze()
    trans_data = cross_section(trans_data,start_point, end_point)

#this cross-section uses the selection based on latitude then, longitude
elif which_type=='EW':
    center_lat = 67 #65 #latitudinal center of the system
    center_lon = 8# longitudinal center of the system
    lon_span = 21.5 # degrees east and west of the center
    
    # extract the desired latitude
    trans_data = data[VB].sel(latitude=center_lat, method='nearest').sel(
        longitude = slice(center_lon - lon_span , center_lon + lon_span ))
           
elif which_type=='NS':
    center_lat = 55 #
    center_lon = 15# 12 # + 360 # Bloomington, IN
    lat_span = 20 # degrees
    
    # extract the desired latitude
    trans_data = data[VB].sel(longitude=center_lon, method='nearest').sel(
        latitude = slice( center_lat + lat_span,center_lat - lat_span ))
    
    
    
print('Processing PCA')
#fill nan values with the mean value
mean = trans_data.mean()
data_array = trans_data.fillna(mean).values


# Reshape the data to be 2D (time x space)
data_2d = data_array.reshape((len(trans_data['time']), -1))

# Perform EOF analysis (PCA)
from sklearn.decomposition import PCA

pca = PCA(n_components=n_components)
pca.fit(data_2d)

# Get the EOFs and explained variance ratios
eofs = pca.components_
explained_variance = pca.explained_variance_ratio_

EOF_dt = [xr.DataArray(data=eofs[i].reshape(trans_data['level'].size, trans_data[tloc].size), dims=['level',tloc],
                      coords={'level':trans_data['level'].values,tloc:trans_data[tloc].values}) for i in np.arange(0,n_components)]

EOF_dt = xr.concat(EOF_dt,dim='nodes')
EOF_Data = xr.Dataset()
EOF_Data['EOF'] = EOF_dt#.coords['var'] = explained_variance

EOF_Data['var'] = explained_variance

print('Saving ...')
########################################################################## SAVING DATA ########################################################

sv_path = f'EOF_data_{vb}'
#create paths if they do not exist!
if os.path.isdir(sv_path)==False:
    os.mkdir(sv_path)
    
    
#if os.path.isdir(f'{sv_path}/eof_{fld}_{vb}') == False:
#   os.mkdir(f'{sv_path}/eof_{fld}_{vb}')

enc={'EOF' : {'zlib':True,'complevel': 5, 'fletcher32': True }}
EOF_Data.to_netcdf(f'{sv_path}/new_e5.trans.{vb}.{fld}.eof.nc', format="NETCDF4",engine='netcdf4', encoding=enc)
#Save data in chunks using custom module saver
#saver.ch_paths(EOF_dt,f'{sv_path}/eof_{fld}_{vb}',f'e5.trans.{vb}')  'chunksizes': (37, len(trans_data[tloc].values),n_components,n_components)
print('Done')