#!/usr/bin/env python3
# -*- coding: utf-8 -*-


exec(open('imports.py').read())

#data paths 
clim_data_path = ['/global/cscratch1/sd/kquagra/era5_IVT/']#,
'''/global/cscratch1/sd/kquagra/ARs Work/Non-coinciding ARDTs/Climatologies/msl_climatologies/',
  '/global/cscratch1/sd/kquagra/ARs Work/Non-coinciding ARDTs/Climatologies/tcwv_climatologies/',
  '/global/cscratch1/sd/kquagra/new-ERA5_PV/',
  '/global/cscratch1/sd/kquagra/ARs Work/Non-coinciding ARDTs/Climatologies/z_climatologies/',
 '/global/cscratch1/sd/kquagra/ARs Work/Non-coinciding ARDTs/Climatologies/temp_climatologies/',
 ]'''
#variable names 
vbls = ['ivt']#,'msl','tcwv','pv','z','t']

#files paths 
files = [glob.glob(f'{clim_data_path[vx]}*.nc',recursive=True) for vx,v in enumerate(vbls)]
#files.sort()


#iterate through dataset and obtain the climatology as a dictionary
climatology=collections.defaultdict(list)
vbls = ['ivt']#,'msl','tcwv','pv','z','t']
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
    elif px==0 or px ==3: #working on data that is arranged differently from the others 
        for yy in range(1980,2018):
            dim='time'
            ds = xr.open_mfdataset(glob.glob(f'{pts}*{yy}*.nc',recursive=True),parallel=True) #open data 
            ds = ds.sel(time=ds.time.dt.month.isin([1,2,12])) #select the months 
            ds.coords['longitude'] = (ds.coords['longitude']  + 180) % 360 - 180 #fix lons 
            ds=ds.sortby(ds.longitude)
            ds = ds.sel(longitude=slice(-170,-60),latitude=slice(55,10)) 
            
            dd.append(ds.mean('time')) #append the averages 
        new_dd = xr.concat(dd,dim=dim) #concat the yearly averages for DJF
    climatology[vbls[px]].append(new_dd.mean(dim)) #save the variable DJF climatology 
    
    

#######################################################################################################
#                                                                                CALCULATE ANOMALIES                                                                               #
#######################################################################################################
    
non_con_paths = []
flds =['guan_waliser']#'guan_waliser','reid250','mundhenk']
vbls=['IVT']#,'MSL','TCWV','PV','GEOPTH','TEMP']
cascade_bard_v1 = collections.defaultdict(list)
guan_waliser = collections.defaultdict(list)
reid250 = collections.defaultdict(list)
mundhenk = collections.defaultdict(list)

for cx,c in enumerate(list(climatology.keys())):
    for fx,f in enumerate(flds):
        clim = climatology[c][0]
        clim.coords['longitude'] = (clim.coords['longitude']  + 180) % 360 - 180
        clim=clim.sortby(clim.longitude)
        
        clim = clim.load()
        if cx<=2:
            ds = xr.open_dataset(f'../Final_variable_data/{vbls[cx]}_consensus_point/{f}True_concensus_{vbls[cx]}.nc'); ds.load()
            
        if cx>2: 
            ds = xr.open_dataset(f'../Final_variable_data/{vbls[cx]}_consensus_point/{f}True_concensus_{vbls[cx]}.nc').sel(level=500); ds.load()
            
        ds.coords['longitude'] = (ds.coords['longitude']  + 180) % 360 - 180
        ds=ds.sortby(ds.longitude)
        anom = ds- clim
        print(anom)
        eval(f)[c].append(anom)

#######################################################################################################
#                                                                                BOOTSTRAPPING                                                                                           #
#######################################################################################################

#create dictionary to hold the data 
cascade_bard_v1_boot = collections.defaultdict(list)
guan_waliser_boot = collections.defaultdict(list)
mundhenk_boot = collections.defaultdict(list)
reid250_boot = collections.defaultdict(list)
vb=['IVT']#,'MSL','TCWV','PV','Z','T']
flds = ['guan_waliser']#,'guan_waliser','reid250','mundhenk']

print('Beginning Bootstrap computations')
for fx,f in enumerate(flds):  #loop through algorithms 
    
    for vx,v in enumerate(vb):  #loop through the variables 
        data = [] #hold data temporarily 
        n=1000
        for num_times in range(n):  #loop for bootstrap sample n times 
        
            dtime = eval(f)[v.lower()][0].time.values #get the time values 
            dtim = np.random.choice(dtime,size=len(dtime))  #randomly sample the times 
            df = eval(f)[v.lower()][0][v]  #data 
            df = df.sel(time=dtim).mean('time') #select the randomly sampled time 
            df.coords['tstep'] = num_times #create timestamp for the data 
            data.append(df)   #
        concated = xr.concat(data,dim='tstep') #concatenate all the bootstrapped data 
        #eval(f'{f}_boot')[v].append(concated)
        concated.to_netcdf(f'Bootstrapped Data/{f}_consensus-boot-{v}.nc') #store as netcdf 
print('Done with all procedures')

'''
print('Done with computations, plotting data ...')
### PLOT THE CLIMATOLOGIES 
units = ['kgm$^{-1}$s$^{-1}$' ,'hPa' ,'kgm$^{-2}$' ,'m$^2$s$^{-1}$Kkg$^{-1}$' ,'m' ,'K']
vbls=['IVT','MSL','TCWV','PV','GEOPTH','TEMP']
flds = ['cascade_bard_v1','guan_waliser','reid250','mundhenk']
vb=['IVT','MSL','TCWV','PV','Z','T']
vmin = [-300, -10,-10,-3e-7,-150,-4]
vmax= [300,10,10,3e-7,150,4]
for fx,f in enumerate(flds):
    fig,axes=plt.subplots(ncols=3,nrows=2,figsize=(15,8),subplot_kw={'projection':ccrs.PlateCarree()})
    ct = 0
    for v,x in zip(vb,axes.flat):
        if v == 'MSL':
            ((eval(f'{f}_boot')[v][0]/100).mean('tstep')- (eval(f)[vb[ct].lower()][0][vb[ct]]/100).mean('time')).plot(cmap='bwr',ax=x,cbar_kwargs={'label' : f'{v} ({units[ct]})','orientation':'horizontal'})
        elif v == 'Z':
            ((eval(f'{f}_boot')[v][0]/10).mean('tstep')- (eval(f)[vb[ct].lower()][0][vb[ct]]/10).mean('time')).plot(cmap='bwr',ax=x,cbar_kwargs={'label' : f'{v} ({units[ct]})','orientation':'horizontal'})
        else:
            (eval(f'{f}_boot')[v][0].mean('tstep')-eval(f)[vb[ct].lower()][0][vb[ct]].mean('time')).plot(cmap='bwr',ax=x,cbar_kwargs={'label' : f'{v} ({units[ct]})','orientation':'horizontal'})
        cart_plot(x)
        x.set_title(f'{v}')
        ct+=1
    plt.savefig(f"anomaly-images/diff-Bootstrap-DJF-anomaly{f}.png",facecolor='white')
'''
