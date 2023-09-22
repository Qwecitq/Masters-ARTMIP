import itertools
import os
import xarray as xr
import collections
import glob
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
from cartopy import feature as cf
import pandas as pd
from functools import reduce
xr.set_options(display_style="text")
import sys
import re #used for splitting 
import warnings
warnings.filterwarnings('ignore')
import numpy as np


def cart_plot(axes):
    axes.add_feature(cf.BORDERS,linewidth=0.5) #add features (boarders )
    axes.add_feature(cf.COASTLINE,linewidth=0.5) #coastlines 
    gl=axes.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.01, linestyle='--') #create the map type with grids 
    axes.set_extent([-150,-60,15,42]) #set the extent of the map
    gl.top_labels = False; gl.left_labels = True #set which axis labels should show 
    gl.right_labels=False; gl.xlines = True
    axes.xformatter = LONGITUDE_FORMATTER;  axes.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'color':'black','size':12} #change font style for x-axis 
    gl.ylabel_style = {'color':'black','size':12}#change font style for y-axis 

def load_data(alg_names , vbl_nm, fpath , dname, wgt=False,eof =False, lvl=None, nod='single', main_path ='/global/cscratch1/sd/kquagra/ARs Work/'):
    
    '''
    This function loads multiple datasets if there are diferent datasets in the same base location that needs to be loaded into a dictionary. 
    Applying the "eof = True " initiates the computation of the standardised version of the data. 
    
    This function will save EOF datasets for EOF analysis to a new path that is consistent with the name 'EOF_STDs/std_eof_data_{vbl_nm}_{f}.nc'
    
    
    Arguments: 
    alg_names : names of the algorithms in a list 
    vblnm : the specific variable name as it appears in the dataset and the label name
    fpath : the specific path different from the main source path to the entire directory
    dname : common name of dataset
    eof : a boolean ; True or False. This will cause the standardisation of the data if True
    wgt : a boolean ; determine if the std data should be weighted by the levels or any number
    lvl :  list of the levels that are under consideration. 
    nod : number of datasets to open [takes 'multiple' or 'single']
    '''
    path = main_path+fpath
    print(path)
    dt = collections.defaultdict(list)
    
    
    #find the weight for the standardised data 
    if wgt == True:
        wgt = len(lvl) #set weight to the len of data levels 
    elif wgt == False:
        wgt = 1
    
    for fx,f in enumerate(alg_names):
        if f.startswith('/') and f.endswith('/'): #check if algorithms are the same as paths and remove the /
            f=f[1:-1]
        else:
            f = f 
        
        #####################################################################################################################################################
        ###################################################       If there are multiple levels in the datasets, use this       ################################################### 
        #####################################################################################################################################################
        # If we have multiple data levels, use this 
        if lvl != None:
            if nod =='multiple':
                ds = xr.open_mfdataset(f'{path}{f}_{vbl_nm}/*.nc',parallel = True).sel(level=lvl)  # open the datasets 
            elif nod == 'single':
                ds = xr.open_dataset(f'{path}/{f}{dname}_{vbl_nm}.nc').sel(level=lvl)  # open the datasets 
            print(f'Done loading the datasets : {f}')
            
            #####################################################################################################################################################
            ################################ Check if the data is for an EOF analysis. If Yes, Standardise the data and save it to a path that will be created #############################
            #####################################################################################################################################################
            if eof == True: #check if the operator wants just the data or the eof prepared data
                
                print(' Standardising data ')
                vnm = list(ds.keys())[0] ; print(f'In data varaible name : {vnm}')
                ds = ds[vnm]  # Get the variable name in the dataset so that we have just an array of that 
                stddiv = ((ds) - np.nanmean(ds))/ (np.nanstd(ds) * wgt)  #compute the standard deviations for the data 
                
                print(' Output standardised data ')
                dstd = stddiv.fillna(np.nanmean(stddiv)) #fill nan values with the mean stddiv
                 
                #store the calculated Standardised data files to a path std_path
                std_path = f'EOF_STDs_{vbl_nm}_new/' #path to store the std files 
                if os.path.isdir(std_path) == True:  #condition met? 
                    print('Path Found \n Storing File')
                    dstd.to_netcdf(f'{std_path}std_eof_data_{vbl_nm}_{f}.nc',engine='h5netcdf') #store file 
                    print('Done Storing FIle')
                    
                elif os.path.isdir(std_path) == False: #condition not met? 
                    print('Path not Found \n Creating the Directory')
                    os.mkdir(std_path) #create the directory path
                    print('Done creating Directory, \n Storing File')
                    dstd.to_netcdf(f'{std_path}std_eof_data_{vbl_nm}_{f}.nc',engine='h5netcdf') #store file 
                    print('Done Storing FIle')
                    
        
                dt[f].append(dstd) #append the standard deviations to the variable assigned
                print(' Done ')
                
            elif eof == False:
                print(' Output data ')
                dt[f].append(ds)
                print(' Done ')
                
                
       #####################################################################################################################################################
       ###################################                                          Check if the data has one level. If Yes, Use this                                   ######################################
       #####################################################################################################################################################
        elif lvl == None:
            
            if nod =='multiple':
                ds = xr.open_mfdataset(f'{path}{f}_{vbl_nm}/*.nc',parallel = True).sel(level=lvl)  # open the datasets 
            elif nod == 'single':
                ds = xr.open_dataset(f'{path}/{f}{dname}_{vbl_nm}.nc') # open the datasets 
            print(f'Done loading the datasets : {f}')
            
            #####################################################################################################################################################
            ################################ Check if the data is for an EOF analysis. If Yes, Standardise the data and save it to a path that will be created #############################
            #####################################################################################################################################################
            if eof == True:  #check if the operator wants just the data or the eof prepared data
               
                print( ' Standardising data ')
                vnm = list(ds.keys())[0] ; print(f'In data varaible name : {vnm}')
                ds = ds[vnm]  # Get the variable name in the dataset so that we have just an array of that 
                stddiv = ((ds) - np.nanmean(ds))/ (np.nanstd(ds) * wgt)  #compute the standard deviations for the data 
                
                print(' Output standardised data ')
                dstd = stddiv.fillna(np.nanmean(stddiv))
                
                #store the calculated Standardised data files to a path std_path
                std_path = f'EOF_STDs_{vbl_nm}_new/' #path to store the std files 
                if os.path.isdir(std_path) == True:  #condition met? 
                    print('Path Found \n Storing File')
                    dstd.to_netcdf(f'{std_path}std_eof_data_{vbl_nm}_{f}.nc',engine='h5netcdf') #store file 
                    print('Done Storing FIle')
                    
                elif os.path.isdir(std_path) == False: #condition not met? 
                    print('Path not Found \n Creating the Directory')
                    os.mkdir(std_path) #create the directory path
                    print('Done creating Directory, \n Storing File')
                    dstd.to_netcdf(f'{std_path}std_eof_data_{vbl_nm}_{f}.nc',engine='h5netcdf') #store file 
                    print('Done Storing FIle')
                    
                dt[f].append(dstd) #append the standard deviations to the variable assigned
                print(' Done ')
                
            elif eof == False:
                print(' Output data ')
                dt[f].append(ds)
            
    if eof == False:
        return dt
    else:
        print(f'This is for EOF analysis, so the data is being saved to {std_path} ')
        

##################################################################################################################
#                                                                                                 SINGLE DATA LOAD                                                                                             # 
# #################################################################################################################
def single_data(path,algo,common_fname=None):   
    '''This function returns a dictionary of data from different algorithms which are single files on their own
    
    Attributes
    -------------
    
    path : this is the path of the data to be loaded 
    algo: this is the path of the data of the various algorithms
    common_fname : Default is  : None. This can be replaced with a string of values that represent the common name in the input dataset naming  
    '''
    output = collections.defaultdict(list)
    for ax,al in enumerate(algo):
        if al.startswith('/'):
            al = al[1:-1]
        else:
            al = al 
        main_path = f'{path}{al}'   
        print(main_path)
        output[al].append(xr.open_dataset(f'{main_path}{common_fname}'))
    return output

##################################################################################################################
#                                                                                                 ARDTs DATA LOAD                                                                                             # 
# #################################################################################################################
def ARDTs_load(path,algo,common_fname=None):   
    '''This function returns a dictionary of data from different algorithms which are single files on their own
    
    Attributes
    -------------
    
    path : this is the path of the data to be loaded 
    algo: this is the path of the data of the various algorithms
    common_fname : Default is  : None. This can be replaced with a string of values that represent the common name in the input dataset naming  
    '''
    output = collections.defaultdict(list)
    for ax,al in enumerate(algo):
        if al.startswith('/'):
            al = al[1:-1]
        else:
            al = al 
        main_path = f'{path}{al}/{al}'   
        print(main_path)
        output[al].append(xr.open_mfdataset(f'{main_path}_{common_fname}*.nc'))
    return output


def cons_data(data_dict,consensus_days ,algorithms , time_var,vbl_name,sv_path,invert=False):
    
    """
    This fuction selects consensus / non-consensus days in a collection of datasets based on their algorithms.
    
    Variables
    ------------
    
    data_dict : this takes a dictionary of datasets that are present in the ref_data_dict. This is the input data which will be checked for consensus days.
    consensus_days : A list of timesteps that you want to check if they are in the data or not 
    algorithms : A list of the algorithms (keys) of the dictionary that holds the data. 
    time_var : the name of the time dimension as it appears in the data
    invert : takes a boolean [ True, False ]... to select either consensus days or an inverse of that 
    
    rerturns output: as a dictionary
    """
    if type(data_dict) == dict or type(data_dict ) == collections.defaultdict:  #check if the storage point of the output is a dictionary or not 
        pass 
    else:
        raise ValueError("The data_dict input should be a dictionary!") #raise a value error if it is not a dictionary 
        
    output = collections.defaultdict(list) #create a collection to output data 
    
    
    for fx,f in enumerate(algorithms): 
        if f.startswith('/'): #check if the algorithm name is a path, remove the / 
            f = f[1:-1]
        else:
             f = f 
        inp1 = data_dict[f][0] #get the first algorithm 
        out = inp1.sel(time=np.isin(inp1[time_var].values,consensus_days,invert=invert)) #results 
        
        if os.path.isdir(sv_path) == False:
            os.mkdir(sv_path)
        elif os.path.isdir(sv_path) == True:
            print(sv_path)
        
        out.to_netcdf(f'{sv_path}/{f}{not invert}_concensus_{vbl_name}.nc')
        
        
        output[f].append(out) #append the results to output 
        
        
    return output
