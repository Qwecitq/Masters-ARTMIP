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


def create_filepath(ds, prefix='filename', root_path="."):
    """
    Generate a filepath when given an xarray dataset
    """
    start = str(ds.time.data[0])[:10]
    end = str(ds.time.data[-1])[:10]
    
    filepath = f'{root_path}{prefix}_{start}.nc'
    return filepath

def ch_paths(data,loc,pfix,vb):
    in_data = list(split_by_chunks(data))
    paths = [create_filepath(ds,prefix=pfix, root_path=loc) for ds in in_data]
    enc={vb : {'zlib':True,'complevel': 5, 'fletcher32': True}}
    t=xr.save_mfdataset(in_data,paths,encoding=enc)
    return t
