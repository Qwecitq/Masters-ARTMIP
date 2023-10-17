#!/usr/bin/env python3
# -*- coding: utf-8 -*-

exec(open('imports.py').read())
from metpy.interpolate import cross_section



#################################################################################################################
######################################## THIS IS FOR SAVING DATA IN CHUNKS #######################################
#################################################################################################################


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
    end = str(ds.time.data[-1])[:10]
    
    filepath = f'{root_path}{prefix}_{start}_{end}.nc'
    return filepath

def ch_paths(data,loc,pfix):
    in_data = list(split_by_chunks(data))
    paths = [create_filepath(ds,prefix=pfix, root_path=loc) for ds in in_data]
    t=xr.save_mfdataset(in_data,paths)
    return t

mpath = '/global/cfs/cdirs/m3522/cmip6/ERA5/e5.oper.an.pl/'
vb = sys.argv[1]

start_year = sys.argv[2]
start_year = int(start_year)
end_year = start_year+1

start_month = sys.argv[3]
start_month = int(start_month)
print('Starting loop')
#start a loop for calculating POT for each year and save 
for yy in np.arange(start_year,end_year):
    for mm in np.arange(start_month,start_month+1):
        
        print(f'Starting loop - year : {yy} month: {mm}')
        #load files using glob
        mon = str(mm).zfill(2)
        files = glob.glob(f'{mpath}{yy}{mon}/*_{vb}*.nc',recursive=True)
        print(files)
        open_ds_main = xr.open_mfdataset(files,chunks={'time':24})
        open_ds_main.coords['longitude'] = (open_ds_main.coords['longitude']  + 180) % 360 - 180 #convert from 0-360 to -180 to 180
        open_ds_main=open_ds_main.sortby(open_ds_main.longitude) #sort the lons
        
        #change the longitude and latitude extent to use this script for different regions.
        open_ds_main = open_ds_main.sel(longitude = slice(-20,30),latitude=slice(80,35))#,parallel=True) #open data 
        
        in_data = list(split_by_chunks(open_ds_main))
        #print(in_data)
        sv_data = []
        
        print('Starting the chuncked calculation loop')
        for idx,open_ds in enumerate(in_data):
            dataset=open_ds
    
            save_path = f'{vb.upper()}_climatologies_west_europe/'
            if os.path.isdir(save_path)== False:
                os.mkdir(save_path)
            print(f'Done with chunk calculation loop... \nSaving {vb} for {yy}')  
            enc={f'{vb.upper()}' : {'zlib':True,'complevel': 5, 'fletcher32': True,'chunksizes': (24, 37, 181, 201)}}
            dataset.to_netcdf(f'{save_path}e5.oper.an.pl_{vb}_{yy}{mm}_{idx+1}.nc', format="NETCDF4",engine='netcdf4', encoding=enc)

        #paths = [create_filepath(ds,prefix=pfix, root_path=loc) for ds in sv_data]
        #t=xr.save_mfdataset(sv_data,paths)
        
print('Completed process successfully')
        
        
