exec(open('imports.py').read())
from metpy.interpolate import cross_section

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


def create_filepath(ds, prefix='filename', root_path=".",sv_data_times=None):
    """
    Generate a filepath when given an xarray dataset
    """
    
    #do this if you want to use the stock file naming which adds the first 10 letters of the time variable in the data 
    #else, add your own texts that differentiate the saved data for different chunks as sv_data_times
    if sv_data_times==None:
        
        start = str(ds.time.data[0])[:14]
        end = str(ds.time.data[-1])[:14]

        filepath = f'{root_path}{prefix}_{start}_{end}.nc'
        
    else:
        filepath = f'{root_path}{prefix}.nc'
        
    return filepath

def ch_paths(data,loc,pfix,vb,sv_data_times=False):
    '''
    
    This automates the whole saving process and saves the files for the dataset
    
    Arguments
    --------------
    
    data : data to be saved in netcdf format (`xr.Dataset()`)
    loc :   location to store dataset
    pfix :  prefix for storing the data
    vb :   variable name for data encoding (this should be the same as the variable name found in the dataset)
    sv_data_times : default (`None`) change this to a list of the prefixes that you want to use.
    
    Return ch_paths
    
    '''
    
    in_data = list(split_by_chunks(data))
    #ind_data = sum(in_data[0],())
    print(in_data[0])
    if sv_data_times == 'Numbered':
        
        paths = [f'{loc}{pfix}_{str(sv_num)}' for sv_num in np.arange(1,len(in_data)+1)]
        #start_time = 
        #end_time = 
        #paths = [f'{loc}{pfix}_{ind.time.values[0][:10]}_{ind.time.values[-1][:10]}' for ind in in_data]
        
    elif sv_data_times == False:
        
        paths = [create_filepath(ds,prefix=pfix, root_path=loc, sv_data_times = False) for ds in in_data]
        
    else: 
        paths = [create_filepath(ds,prefix=pfix, root_path=loc, sv_data_times = None) for ds in in_data]
        
        
    enc_dict =  {'zlib':True,'complevel': 5, 'fletcher32': True}
    
    enc = { i: enc_dict for i in vb }
    t=xr.save_mfdataset(in_data,paths,encoding=enc)
    return t

