#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#LOAD PACKAGES
from eofs.xarray import Eof
import xarray as xr
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import os 
import collections 
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
from cartopy import feature as cf
xr.set_options(display_style="text")
from eofs.multivariate.standard import MultivariateEof


#####################################################################################################################################################
################################                                                                      EOF FUNCTION.                                                                             #############################
#####################################################################################################################################################
def stdd(data):
    '''Computes the standard deviation of the data'''
    STDD=(data-data.mean())/np.std(data)
    return STDD

def pca(data,coords,nmodes):
    '''Computes the principal components of the data 
    
    Arguments
    --------------
    data : the data for the computation 
    coords: the specific coordinates needed for the computation (time)
    nmodes: number of modes that should be computed'''
    
    #working on the weights 
    coslat = np.cos(np.deg2rad(data.coords[coords].values)) #convert the longitudes into specific weighted longitudes 
    wgts = np.sqrt(coslat)[...,np.newaxis]
    
    #create the solver for the EOF
    solver = Eof(data,weights=None) #eof solver 
    
    
    #calculate the EOFs 
    eof1=solver.eofsAsCovariance(neofs=nmodes) #compute the covariance 
    pc1 = solver.pcs(npcs=nmodes,pcscaling=1) #compute the principal components
    eigenv2 = solver.eigenvalues(neigs=nmodes) #compute the eigen values 
    variance_fractions = solver.varianceFraction(neigs=nmodes) #compute variance fractions 
    return pc1,eof1,variance_fractions # returns pc results 

def cart_plot(axes,rows,cols):
    axes[rows,cols].add_feature(cf.BORDERS,linewidth=1)
    axes[rows,cols].add_feature(cf.COASTLINE,linewidth=1)
    gl=axes[rows,cols].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.01, linestyle='--')
    #axes[rows,cols].set_extent([-3.2,0.9,8.9,11.5])
    gl.top_labels = False; gl.left_labels = True
    gl.right_labels=False; gl.xlines = True
    axes[rows,cols].xformatter = LONGITUDE_FORMATTER;  axes[rows,cols].yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'color':'black','size':12}
    gl.ylabel_style = {'color':'black','size':12}

flds = ['/cascade_bard_v1/','/guan_waliser/','/mundhenk/','/reid250/']

#####################################################################################################################################################
################################                                                USE TO RECALIBERATE YOUR LONS AND LATS                                              #############################
#####################################################################################################################################################

path = f'/global/cscratch1/sd/kquagra/ARs Work/Non-coinciding ARDTs/TEMP_consensus_point/'
otpath='/new_new_algo/'
dlist = os.listdir(path)

lats = np.array(xr.open_dataset(path+dlist[1])['latitude'])
lons = np.array(xr.open_dataset(path+dlist[1])['longitude'])

#####################################################################################################################################################
################################                                                                      LOAD DATASETS                                                                            #############################
#####################################################################################################################################################

#IVT
path = f'/global/cscratch1/sd/kquagra/ARs Work/EOF_STDs'
vbl_nm = 'IVT'
fpath = '_IVT'

ivt_col = collections.defaultdict(list)

for f in flds:
    f = f[1:-1]
    
    ivt_col[f].append(xr.open_dataset(f'{path}{fpath}_new/std_eof_data{fpath}_{f}.nc')[vbl_nm])
print(f'Done Loading {vbl_nm}')
    
    
#GEOPTH
path = f'/global/cscratch1/sd/kquagra/ARs Work/EOF_STDs'
vbl_nm = 'Z'
fpath = '_GEOPTH'

gpt_col = collections.defaultdict(list)

for f in flds:
    f = f[1:-1]
    
    gpt_col[f].append(xr.open_dataset(f'{path}{fpath}_new/std_eof_data{fpath}_{f}.nc')[vbl_nm]/9.8) #open geopotential data and convert to geopotential heights 
print(f'Done Loading {vbl_nm}')


#PV   
path = f'/global/cscratch1/sd/kquagra/ARs Work/EOF_STDs'
vbl_nm = 'PV'
fpath = '_PV'

pv_col = collections.defaultdict(list)

for f in flds:
    f = f[1:-1]
    
    pv_col[f].append(xr.open_dataset(f'{path}{fpath}_new/std_eof_data{fpath}_{f}.nc')[vbl_nm])
print(f'Done Loading {vbl_nm}')

    
#MSL
path = f'/global/cscratch1/sd/kquagra/ARs Work/EOF_STDs'
vbl_nm = 'MSL'
fpath = '_'+vbl_nm

msl_col = collections.defaultdict(list)

for f in flds:
    f = f[1:-1]
    
    msl_col[f].append(xr.open_dataset(f'{path}{fpath}_new/std_eof_data{fpath}_{f}.nc')[vbl_nm])
print(f'Done Loading {vbl_nm}')

    
#TCWV
path = f'/global/cscratch1/sd/kquagra/ARs Work/EOF_STDs'
vbl_nm = 'TCWV'
fpath = '_'+vbl_nm

tcwv_col = collections.defaultdict(list)

for f in flds:
    f = f[1:-1]
    
    tcwv_col[f].append(xr.open_dataset(f'{path}{fpath}_new/std_eof_data{fpath}_{f}.nc')[vbl_nm])
print(f'Done Loading {vbl_nm}')


#TEMP
path = f'/global/cscratch1/sd/kquagra/ARs Work/EOF_STDs'
vbl_nm = 'T'
fpath = '_TEMP'

temp_col = collections.defaultdict(list)

for f in flds:
    f = f[1:-1]
    
    temp_col[f].append(xr.open_dataset(f'{path}{fpath}_new/std_eof_data{fpath}_{f}.nc')[vbl_nm])
print(f'Done Loading {vbl_nm}')

    
#STACK THE DATASET 
print(f'Begining to Stack Data')
big_data = collections.defaultdict(list)
vbls = ['gpt_col','ivt_col','tcwv_col','msl_col','temp_col','pv_col']
vb_nms = ['gpt','ivt','tcwv','msl','temp','pv']

for fx , f in enumerate(flds):
    f = f[1:-1]
    #load the dataset into a stack
    dds = xr.Dataset({vb_nms[0] : eval(vbls[0])[f][0], 
                      vb_nms[1]  : eval(vbls[1])[f][0],
                      vb_nms[2] :eval(vbls[2])[f][0],
                      vb_nms[3] : eval(vbls[3])[f][0],
                      vb_nms[4] : eval(vbls[4])[f][0],
                      vb_nms[5] : eval(vbls[5])[f][0]})
    print('Done creating dataset for the stack \n Begin stacking ... ')
    dds=dds.to_stacked_array('eof_vars', sample_dims=['longitude','latitude','time'] )   
    big_data[f].append(dds)#xr.open_dataset(f'{direc}{f}_stacked.nc'))#[vb_nms[0]])
    print('Done stacking the dataset \n Begining the EOF analysis on the stacked dataset ... ')
    
#Create a Collection to hold all the eof computations 

data_eof = collections.defaultdict(list) #dictionary for eofs of each ARDT
data_pc = collections.defaultdict(list) #dictionary for pcs of each ARDT
data_var = collections.defaultdict(list) #dictionary for variance of each ARDT
array_vbl = 'gpt'   # variable naarray_vbl the first variable in the stacked data. 

for fx, f in enumerate(flds):
    f = f[1:-1]
    dd = big_data[f][0]#.load() #lload the data into memory 
    print(f'Done Loading data [{f}] into memory \n Begining transpose and run EOF on {f}')
    dd = dd.transpose('time','eof_vars','longitude','latitude')  #tranpose data 
    data_pc[f],data_eof[f], data_var[f] = pca(dd,'latitude',6)  #calculate EOFs 
    
    print('Done computing the EOF \n Working on the longitudes and latitudes to be rearranged')
    data_eof[f]['longitude'] = lons
    data_eof[f]['latitude'] = lats
    print('Done rearranging longitudes and latitudes \n unstacking the data for plotting...')
    data_eof[f] = data_eof[f].to_unstacked_dataset(dim='eof_vars')  #unstack EOFs
    print(f'Done unstacking the data for {f} \n Saving computations')
    
    di1 = 'EOF_out/';
    di2 = 'PC_out/'
    di3 = 'VAR_out/'
    
    print('Checking for directories')
    if os.path.isdir(di1)==False:
            os.mkdir(di1)
        
    if os.path.isdir(di2)==False:
            os.mkdir(di2)
    
    if os.path.isdir(di3)==False:
            os.mkdir(di3)
    
    print(f'Done. \n Saving data from computation of eof,pc and var of {f}')
    data_eof[f].to_netcdf(f'{di1}EOF_output_{f}.nc')  #create the .nc eof files 
    data_pc[f].to_netcdf(f'{di2}PC_output_{f}.nc') #create the .nc pc files 
    data_var[f].to_netcdf(f'{di3}VAR_output_{f}.nc') #create the .nc variance files 
    
print('Done working on the EOF and Data Unstacking. \n Ready to Plot...')  
############################################################################################################################################
############################                                                                        Plot Section                                                                               ###########################
############################################################################################################################################
 
print('Begining the plotting sequence ...')
svb_nms = ['ivt','tcwv','mslp','geopth','pv','temp']
v_levs = [1000,950,850, 700, 500, 200]

for fx, f in enumerate(flds):
    
    f=f[1:-1]
    print('Creating Figure and Canvas for Images')
    fig,axes = plt.subplots(figsize=(16,8), ncols=3,nrows=2, subplot_kw={'projection':ccrs.PlateCarree()})
    plt.suptitle(f'{f}') 
    ct = 0 #counter 
    print('Making images from the data')
    
    mode_n = 0        #mode number 
    sign_lon = -135 ; sign_lat = 40        #lonlat for eof sign selection
    slp = data_eof[f]['msl'][mode_n]
    
    #change sign for easy interpretation
    if slp.sel(longitude=sign_lon, latitude = sign_lat) > 0 :
        eof_sign = -1
    else:
        eof_sign = 1 
    
    
    try:
        for i in range(2):
            for j in range(3):

                plt_data = data_eof[f][vb_nms[ct]][mode_n]*eof_sign #change EOF sign 
                #PLOT 
                if ct < 3:
                    plt_data.T.plot(ax=axes[i,j],vmin=-1,vmax=1,cmap='bwr' , cbar_kwargs = {'orientation':'horizontal','label':'EOFs'})
                    axes[i,j].set_title(f'{vb_nms[ct]}')
                else: 
                    lev = 0     #vertical levels for pv and temp 
                    (plt_data[lev]*6).T.plot(ax=axes[i,j],cmap='bwr' , cbar_kwargs = {'orientation':'horizontal','label':'EOFs'})                    
                    axes[i,j].set_title(f'{vb_nms[ct]} at {str(v_levs[lev])}')
                    
                cart_plot(axes,i,j)

                ct+=1
        
    except:
        pass
    
    lev=0
    print('Done with plot. \n Begining to save images ...')
    dirs = f'New_Images/new_levels_geopth_mean/M{mode_n+1}L{v_levs[lev]}/'
    print(f'The images would be saved in {dirs}')
    if os.path.isdir(dirs) == False:
        print('Path not found ... \n Creating path ... ')
        os.mkdir(dirs)
    print(f'Saving Images to {dirs}')
    plt.savefig(f'{dirs}eof_{f}.png',facecolor='white',dpi=300)

print('Done with all the plots and data analysis successfully!')