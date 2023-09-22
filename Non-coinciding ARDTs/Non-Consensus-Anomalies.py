#!/usr/bin/env python3
# -*- coding: utf-8 -*-


exec(open('imports.py').read())
reg = 'west-europe-'
#data paths 
clim_data_path = ['/global/cscratch1/sd/kquagra/era5_IVT/',
                  '/global/cscratch1/sd/kquagra/ARs Work/Non-coinciding ARDTs/Climatologies/msl_climatologies/',
                  '/global/cscratch1/sd/kquagra/ARs Work/Non-coinciding ARDTs/Climatologies/tcwv_climatologies/',
                  '/global/cscratch1/sd/kquagra/new-ERA5_PV/',
                  '/global/cscratch1/sd/kquagra/ARs Work/Non-coinciding ARDTs/Climatologies/z_climatologies/',
                 '/global/cscratch1/sd/kquagra/ARs Work/Non-coinciding ARDTs/Climatologies/temp_climatologies/',
                 ]
#variable names 
vbls = ['ivt','msl','tcwv','pv','z','t']

#files paths 
files = [glob.glob(f'{clim_data_path[vx]}*.nc',recursive=True) for vx,v in enumerate(vbls)]
#files.sort()


#iterate through dataset and obtain the climatology as a dictionary
climatology=collections.defaultdict(list)
vbls = ['ivt','msl','tcwv','pv','z','t']
for px,pts in enumerate(clim_data_path):
    dd =[]
    
    #check if the file name is for msl,tcwv,z or temp since they have a different arrangement 
    if 0<px <=2 or 3<px<=5: 
        for yy in range(1980,2018): #loop through the files year by year 
            dim='month'  #set the concatenation dimension 
            tm=pd.date_range(f'{yy}-01',f'{yy+1}-01',freq='1M') #create the time axis 
            tm = [np.datetime64(x) for x in tm]  #convert to datetime64
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
