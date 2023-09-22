#!/usr/bin/env python3
# -*- coding: utf-8 -*-

exec(open('imports.py').read())

################### COMPUTE CLIMATOLOGIES FOR THE DIFFERENT VARIABLES IN DJF SEASON ########################

####################################################### MAIN SCRIPT #########################################################        

flds = ['mundhenk_v2']#['cascade_bard_v1','guan_waliser','reid250','mundhenk_v2']
mpath = '/global/cfs/cdirs/m3522/cmip6/ERA5/'
vbls =['u']# ['msl','tcwv','pv','z','t']
locs =['e5.oper.an.pl/']#,'e5.oper.an.pl/','e5.oper.an.pl/']# ['e5.oper.an.sfc/','e5.oper.an.sfc/','e5.oper.an.pl/','e5.oper.an.pl/','e5.oper.an.pl/'] #['e5.oper.an.sfc/','e5.oper.an.sfc/','e5.oper.an.pl/',]
files = [glob.glob(f'{mpath}{y}*/*{x}*',recursive=True) for x,y in zip(vbls,locs)]
ivt_path = '/pscratch/sd/k/kquagra/era5_IVT/'
ivt_dt = [glob.glob(f'{ivt_path}*.nc',recursive=True)]
files.insert(0,ivt_dt)  #insert IVT path into the files path

clim_data_path = [mpath+i for i in locs]
#clim_data_path.insert(0,ivt_path)
#sum(clim_data_path, [])

print(clim_data_path)

#iterate through dataset and obtain the climatology as a dictionary
climatology=collections.defaultdict(list)
vbls =['u']# ['ivt','msl','tcwv','pv','z','t']
print('Beginning to compute climatologies ...')


for vb,pts in zip(vbls,clim_data_path):
    dd = []
    #select years to reduce memory usage
    for yy in range(1980,2000):
        
        if pts == '/pscratch/sd/k/kquagra/era5_IVT/' : 
            mpt = glob.glob(f'{pts}*{yy}*.nc',recursive=True) #new path to specific variable
            mpt.sort()
            
            #Load data
            ndt= xr.open_mfdataset(mpt , chunks=dict(time=1,level=1),parallel=True)
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
    nds = ds.sel(time=ds.time.dt.month.isin([1,2,12]))
    #climatology calculated 
    
    #create climatology path
    save_path = 'DJF_Climatologies/'
    if os.path.isdir(save_path)== False:
        os.mkdir(save_path)
        
    print('started saving climatoloical data')
    mean_ds = nds.mean('time')
    mean_ds.to_netcdf(f'{save_path}{vb}_1980_2000_climatology_001.nc')
    print('Done saving successfully')
    #climatology[fd] = mean_ds#.load()
