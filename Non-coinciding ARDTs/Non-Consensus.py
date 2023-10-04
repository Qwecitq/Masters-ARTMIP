#!/usr/bin/env python3
# -*- coding: utf-8 -*-


exec(open('imports.py').read())
#change the region for other regions
reg = sys.argv[1] #'north_america'
import itertools

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
    start = str(ds.time.data[0])[:15]
    end = str(ds.time.data[0])[:15]
    filepath = f'{root_path}{prefix}_{start}_{end}.nc'
    return filepath


def ch_paths(data,loc,pfix):
    in_data = list(split_by_chunks(data))
    paths = [create_filepath(ds,prefix=pfix, root_path=loc) for ds in in_data]
    t=xr.save_mfdataset(in_data,paths,mode='w')
    return t

####################################################### MAIN SCRIPT #########################################################        

flds=[sys.argv[2]]
vbls = [sys.argv[3]]
typ = sys.argv[4]
print(typ)
print(f'{reg}, {flds}, {vbls}')
#flds = ['cascade_bard_v1']#['cascade_bard_v1','guan_waliser','reid250','mundhenk_v2']
mpath = '/global/cfs/cdirs/m3522/cmip6/ERA5/'
#vbls =['z']# ['msl','tcwv','pv','z','t']

###################################################################################################################################
##########################################Check if the variable being selected is in levels or a surface parameter##################################
###################################################################################################################################

if typ == 'levels':
    locs =['e5.oper.an.pl/']#,'e5.oper.an.pl/','e5.oper.an.pl/']# ['e5.oper.an.sfc/','e5.oper.an.sfc/','e5.oper.an.pl/','e5.oper.an.pl/','e5.oper.an.pl/']
    
elif typ == 'surface':
    locs = ['e5.oper.an.sfc/']

print(locs)
#####################################################################################################################################


files = [glob.glob(f'{mpath}{y}*/*{x}*',recursive=True) for x,y in zip(vbls,locs)]
ivt_path = '/pscratch/sd/k/kquagra/era5_IVT/'
ivt_dt = [glob.glob(f'{ivt_path}*.nc',recursive=True)]

files.insert(0,ivt_dt)  #insert IVT path into the files path

clim_data_path = [mpath+i for i in locs]

print(clim_data_path)

#iterate through dataset and obtain the climatology as a dictionary
climatology=collections.defaultdict(list)
###################################################################################################################################
###################################################### Start computing climatologies ####################################################
###################################################################################################################################
print('Beginning to compute climatologies ...')


for vb,pts in zip(vbls,clim_data_path):
    dd = []
    #select years to reduce memory usage
    for yy in range(1980,2018):
        
        if pts == '/pscratch/sd/k/kquagra/era5_IVT/' : 
            mpt = glob.glob(f'{pts}*{yy}*.nc',recursive=True) #new path to specific variable
            mpt.sort()
            
            #Load data
            ndt= xr.open_mfdataset(mpt , chunks={'time':10},parallel=True)
            dd.append(ndt)
            
        elif pts !=  '/pscratch/sd/k/kquagra/era5_IVT/':
            mpt = glob.glob(f'{pts}{yy}*/*_{vb}*.nc',recursive=True) #new path to specific variable
            mpt.sort()

            #print(glob.glob(f'{pts}{yy}*/*_{vb}*.nc',recursive=True))
            #Load data
            ndt= xr.open_mfdataset(mpt)
            
            dd.append(ndt)
            #dd.append(ndt2)
        print(f'{yy} loaded')
    ds = xr.concat(dd,dim='time')
    month_ds = ds.resample(time='1M').mean()
    nds = ds.sel(time=ds.time.dt.month.isin([1,2,12]))
    #climatology calculated 
    #mean_ds = nds.mean('time')
    
    
    #use ds to select the specific timesteps corresponding to the ARDTs
    
    print('starting ARDTS consensus')
    for fx,f in enumerate(flds):
        
        print(glob.glob(f'../larger_region_AR_capture/{reg}-{f}/*.nc',recursive=True))
        ar_data = xr.open_mfdataset(glob.glob(f'../larger_region_AR_capture/{reg}-{f}/*0.8-1.0*.nc',recursive=True),parallel=True)#.sel(time=slice('1980','2018'))
        #select AR composites 
        ar_composites = ds.sel(time=ar_data.time)
        print(len(ar_composites.time.values))
        print('saving process initiated ...')
        direc = f'{reg}_{vb}/'
        prefix = f'{vb}'
        
        #for fx,f in enumerate(flds):
        if os.path.isdir(direc)==False:
            os.mkdir(direc)
        dirs = f'{direc}{f}_{vb}/'
        if os.path.isdir(dirs) == False:
            os.mkdir(dirs)
            t=ch_paths(ar_composites, dirs , prefix)
            
    print('Done ...')     
    #save mean data 
    
    #print('Saving mean for variable ...')
    #mean_ds.to_netcdf(f'{direc}mean_1980_2017_{vb}.nc')
    #
    #
    
'''
print('Done loading files')

for px,pts in enumerate(clim_data_path):
    dd =[]
    #pts=pts[0]
    print(pts)
    #check if the file name is for msl,tcwv,z or temp since they have a different arrangement 
    if 0<px <=2 or 3<px<=5: 
        for yy in range(1980,2018): #loop through the files year by year 
            dim='month'  #set the concatenation dimension 
            tm=pd.date_range(f'{yy}-01',f'{yy+1}-01',freq='1M') #create the time axis 
            tm = [np.datetime64(x) for x in tm]  #convert to datetime64
            #print(glob.glob(f'{pts}*{yy}*.nc'))
            ds = xr.open_mfdataset(glob.glob(f'{pts}*{yy}*.nc',recursive=True),parallel=True)  #open the dataset

            ds['month'] =tm  #reassign the time 
            ds=ds.sel(month=ds.month.dt.month.isin([1,2,12]))  #select only DJF
            ds.coords['longitude'] = (ds.coords['longitude']  + 180) % 360 - 180
            ds=ds.sortby(ds.longitude)
            ds = ds.sel(longitude=slice(-170,-60),latitude=slice(55,10))
            dd.append(ds.mean('month')) #append to the dd list 
        new_dd = xr.concat(dd,dim=dim)  #concatenate the list of yearly datasets into one file
    elif px==0 or px ==3:
        for yy in range(1980,2018):
            dim='time'
            ds = xr.open_mfdataset(glob.glob(f'{pts}*{yy}*.nc',recursive=True),parallel=True)
            ds = ds.sel(time=ds.time.dt.month.isin([1,2,12]))
            ds.coords['longitude'] = (ds.coords['longitude']  + 180) % 360 - 180
            ds=ds.sortby(ds.longitude)
            ds = ds.sel(longitude=slice(-170,-60),latitude=slice(55,10))
            
            dd.append(ds.mean('time'))
        new_dd = xr.concat(dd,dim=dim)
    climatology[vbls[px]].append(new_dd.mean(dim)) #save the variable DJF climatology 
print('Done computing climatologies')'''
'''
print('Starting with climatologies ...')
reg='west-europe-'
non_con_paths = []
flds = ['cascade_bard_v1','guan_waliser','reid250','mundhenk']
vbls=['IVT','MSL','TCWV','PV','GEOPTH','TEMP']
cascade_bard_v1 = collections.defaultdict(list)
guan_waliser = collections.defaultdict(list)
reid250 = collections.defaultdict(list)
mundhenk = collections.defaultdict(list)

for cx,c in enumerate(list(climatology.keys())):
    for fx,f in enumerate(flds):
        clim = climatology[c][0]
        
        print('Started climatologies ...')
        #clim.coords['longitude'] = (clim.coords['longitude']  + 180) % 360 - 180
        #clim=clim.sortby(clim.longitude)
        clim = clim.load()
        if cx<=2:
            ds = xr.open_dataset(f'../Final_variable_data/{reg}{vbls[cx]}_consensus_point/{f}True_concensus_{vbls[cx]}.nc')#; ds.load()
        if cx>2: 
            ds = xr.open_dataset(f'../Final_variable_data/{reg}{vbls[cx]}_consensus_point/{f}True_concensus_{vbls[cx]}.nc').sel(level=500)#; ds.load()
        ds.coords['longitude'] = (ds.coords['longitude']  + 180) % 360 - 180
        ds=ds.sortby(ds.longitude)
        anom = ds- clim
        print(anom)
        eval(f)[c].append(anom)

print('Done with climatologies ... \nBegining plots ...')'''
'''
### PLOT THE CLIMATOLOGIES 
units = ['kgm$^{-1}$s$^{-1}$' ,'hPa' ,'kgm$^{-2}$' ,'m$^2$s$^{-1}$Kkg$^{-1}$' ,'m' ,'K']
vbls=['IVT','MSL','TCWV','500 hPa PV','500 hPa GEOPTH','500 hPa TEMP']
vb=['IVT','MSL','TCWV','PV','Z','T']
vmin = [-300, -10,-10,-3e-7,-650,-4]
vmax= [300,10,10,3e-7,650,4]
for fx,f in enumerate(flds):
    fig,axes=plt.subplots(ncols=3,nrows=2,figsize=(15,8),subplot_kw={'projection':ccrs.PlateCarree()})
    ct = 0
    for v,x in zip(vbls,axes.flat):
        if v == 'MSL':
            (eval(f)[vb[ct].lower()][0][vb[ct]]/100).mean('time').plot(cmap='bwr',ax=x,vmin=vmin[ct],vmax=vmax[ct],cbar_kwargs={'label' : f'{v} ({units[ct]})','orientation':'horizontal'})
        elif v == 'Z':
            (eval(f)[vb[ct].lower()][0][vb[ct]]/10).mean('time').plot(cmap='bwr',ax=x,vmin=vmin[ct],vmax=vmax[ct],cbar_kwargs={'label' : f'{v} ({units[ct]})','orientation':'horizontal'})
        else:
            eval(f)[vb[ct].lower()][0][vb[ct]].mean('time').plot(cmap='bwr',ax=x,vmin=vmin[ct],vmax=vmax[ct],cbar_kwargs={'label' : f'{v} ({units[ct]})','orientation':'horizontal'})
        cart_plot(x)
        x.set_title(f'{v}')
        ct+=1
    plt.savefig(f"anomaly-images/{reg}newDJF-anomaly{f}.png",facecolor='white')
print('Done with plots ...')

print('All processes done ...')'''