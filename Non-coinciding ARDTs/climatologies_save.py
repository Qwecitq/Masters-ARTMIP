#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#LOAD PACKAGES
exec(open('imports.py').read())




mpath = sys.argv[1]  #main path for data 
common_name = sys.argv[2]  #common name in files 
vbl = sys.argv[3] #name of variable 
direc = sys.argv[4] #path to store data 
t_start = sys.argv[5] #beginning time 
t_end = sys.argv[6] #ending time 
#lvl = sys.argv[7] #use this if the data has levels 
lvl = input('Type levels separated by comma (,) in the data or else type None foloowed by comma(,)')

if lvl == 'None':
    lvl = None
else:
    lvl = [int(i) for i in re.split(",", lvl)]
    
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
    end = str(ds.time.data[-1])[:10]
    
    filepath = f'{root_path}{prefix}_{start}_{end}.nc'
    return filepath

def ch_paths(data,loc,pfix):
    in_data = list(split_by_chunks(data))
    paths = [create_filepath(ds,prefix=pfix, root_path=loc) for ds in in_data]
    t=xr.save_mfdataset(in_data,paths)
    return t


###############################################################################################
#                                                             Functions for Computing  Anomalies                                                                 #
###############################################################################################

def save_clim(mpath,common_name,vbl,direc,t_start,t_end,lvl=None):
    
    """
    Compute the anomaly for the data selecting the longitudes and the latitudes (boundary box)
    
    Arguments 
    --------------
    
    mpath : paths for files 
    common_name : common name in file storage 
    vbl : variable name 
    direc : directory to save file 
    t_start : start time period
    t_end : end time period
    lvl : if the data has levels set this to a llist of the desired levels
    
    returns : d [climatological mean] 
    
    """
    
    tme = np.arange(int(t_start),int(t_end)+1,1)  #time 
    
    input_files = [] #store data files 
    for tx in tme:
        input_files.append(glob.glob(f'{mpath}*{common_name}*.nc*',recursive=True))  #open specific years
    input_files.sort() #sort files 
    
    
    #check for the storage location file and create it if its not there 
    if os.path.isdir(direc) == False:  #condition to check for the folder 
        print(f'Save folder :{direc}  - not found \nCreating the folder...')
        os.mkdir(f'{direc}') #create the main destination folder 
        
    if os.path.isdir(f'{direc}/{vbl}_climatologies/') == False: #specific variable storage folder 
        os.mkdir(f'{direc}/{vbl}_climatologies/') #create that folder 
    
    #loop through the input files and load for every year 
    for fx,f in enumerate(input_files):
        print(f'loading data \n{f[-1]}')
        
        if lvl!= None:
            d =xr.open_mfdataset(f,parallel=True,chunks={'time':1024,'latitude': 721, 'longitude': 1440}).sel(level=lvl).resample(time='1M').mean() #load data and regroup into monthly means over time 
            d.coords['longitude'] = (d.coords['longitude']  + 180) % 360 - 180 #convert from 0-360 to -180 to 180
            d=d.sortby(d.longitude) #sort the lons 
            d=d.sel(longitude = slice(-170,-30), latitude = slice(60,10)) #select the region under consideration 
            d=d.chunk(chunks={'time':6})
            
        elif lvl == None:
            d =xr.open_mfdataset(f,parallel=True,chunks={'time':1024,'level':37,'latitude': 721, 'longitude': 1440}).resample(time='1M').mean() #load data and regroup into monthly means over time 
            d.coords['longitude'] = (d.coords['longitude']  + 180) % 360 - 180 #convert from 0-360 to -180 to 180
            d=d.sortby(d.longitude) #sort the lons 
            d=d.sel(longitude = slice(-170,-30), latitude = slice(60,10)) #select the region under consideration 
            d=d.chunk(chunks={'time':6})
            
        #print(d)
        print(f'Done loading data \n{d} \nSaving Data')
        ch_paths(d,f'{direc}/{vbl}_climatologies/',f'climatology{fx+int(t_start)}')
        #d.to_netcdf(f'{direc}/{vbl}_climatologies/climatology{fx+int(t_start)}.nc')  #store the files 
        print(f'Done with climatologies for {fx+int(t_start)}')
    
    print('Process complete')

    
#use the multiprocessing tool to work on the data 
def mp_process(mpath,common_name,vbl,direc,t_start,t_end,lvl):

    p = Pool(processes=20)  #call number of processors for the job
    res  = p.map(save_clim(mpath,common_name,vbl,direc,t_start,t_end,lvl))  #call the function
    p.close()
    p.join()

if __name__ == "__main__":
    mp_process(mpath,common_name,vbl,direc,t_start,t_end,lvl) 