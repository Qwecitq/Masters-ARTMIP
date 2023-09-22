#!/usr/bin/env python3
# -*- coding: utf-8 -*-


exec(open('imports.py').read())

#data paths 
'''clim_data_path = ['/global/cscratch1/sd/kquagra/era5_IVT/',
                  '/global/cscratch1/sd/kquagra/ARs Work/Non-coinciding ARDTs/Climatologies/msl_climatologies/',
                  '/global/cscratch1/sd/kquagra/ARs Work/Non-coinciding ARDTs/Climatologies/tcwv_climatologies/',
                  '/global/cscratch1/sd/kquagra/new-ERA5_PV/',
                  '/global/cscratch1/sd/kquagra/ARs Work/Non-coinciding ARDTs/Climatologies/z_climatologies/',
                 '/global/cscratch1/sd/kquagra/ARs Work/Non-coinciding ARDTs/Climatologies/temp_climatologies/',
                 ]
'''
mpath = '/global/cfs/cdirs/m3522/cmip6/ERA5/'
vbls = ['msl']#,'tcwv']#,'pv','z','t']
locs = ['e5.oper.an.sfc/']#,'e5.oper.an.sfc/']#,'e5.oper.an.pl/','e5.oper.an.pl/','e5.oper.an.pl/'] #['e5.oper.an.sfc/','e5.oper.an.sfc/','e5.oper.an.pl/',]
files = [glob.glob(f'{mpath}{y}*/*{x}*',recursive=True) for x,y in zip(vbls,locs)]

print(f'{mpath}{locs[0]}')
print(glob.glob(f'{mpath}{locs[0]}*/*{vbls[0]}*',recursive=True) )

reg='west-europe-'
#variable names 
#vbls = ['ivt','msl','tcwv','pv','z','t']

#files paths 
#files = [glob.glob(f'{clim_data_path[vx]}*.nc',recursive=True) for vx,v in enumerate(vbls)]
#files.sort()


#iterate through dataset and obtain the climatology as a dictionary
climatology = collections.defaultdict(list)

cntr = 0
for fx,fi in enumerate(files):
    print(f'Starting with climatologies for {vbls[fx]} ...')
    if fx <=1:
        ds = xr.open_mfdataset(fi,parallel=True,engine='h5netcdf')
        ds=ds.sel(time=ds.time.dt.month.isin([1,2,12])).mean('time')
        climatology[f'{vbls[cntr]}'].append(ds)
    elif fx>1:
        ds = xr.open_mfdataset(fi,parallel = True).sel(level=500)
        ds=ds.sel(time=ds.time.dt.month.isin([1,2,12])).mean('time')
        climatology[f'{vbls[cntr]}'].append(ds)
    cntr+=1
print('Done ...')

#######################################################################################################
#                                                                                CALCULATE ANOMALIES                                                                               #
#######################################################################################################
    
non_con_paths = []
flds =['mundhenk']#,'mundhenk']# ['cascade_bard_v1','guan_waliser',
vbls=['MSL']#,'TCWV']#,'PV','GEOPTH','TEMP']
vab = ['msl']#,'tcwv']#,'pv','z','t']

cascade_bard_v1 = collections.defaultdict(list)
guan_waliser = collections.defaultdict(list)
reid250 = collections.defaultdict(list)
mundhenk = collections.defaultdict(list)



for cx,c in enumerate(vab):
    print(f'Calculate anomalies for {c} ...')
    for fx,f in enumerate(flds):
        clim = climatology[c][0]
        clim.coords['longitude'] = (clim.coords['longitude']  + 180) % 360 - 180
        clim=clim.sortby(clim.longitude)
        
        #clim = clim.load()
        if cx<=2:
            ds = xr.open_dataset(f'../Final_variable_data/{reg}{vbls[cx]}_non_consensus_point/{f}False_concensus_{vbls[cx]}.nc'); #ds.load()
            
        if cx>2: 
            ds = xr.open_dataset(f'../Final_variable_data/{reg}{vbls[cx]}_non_consensus_point/{f}False_concensus_{vbls[cx]}.nc').sel(level=500); #ds.load()
            
        ds.coords['longitude'] = (ds.coords['longitude']  + 180) % 360 - 180
        ds=ds.sortby(ds.longitude)
        anom = ds- clim
        #print(anom)
        eval(f)[c].append(anom)
print('Done ...')
#######################################################################################################
#                                                                                BOOTSTRAPPING                                                                                           #
#######################################################################################################

#create dictionary to hold the data 
cascade_bard_v1_boot = collections.defaultdict(list)
guan_waliser_boot = collections.defaultdict(list)
mundhenk_boot = collections.defaultdict(list)
reid250_boot = collections.defaultdict(list)
vb=['MSL']#,'TCWV']#,'PV','Z','T']
vab = 'msl'#,'tcwv']#,'pv','z','t']
flds =['mundhenk']#,'mundhenk']# ['cascade_bard_v1','guan_waliser',

print('Beginning Bootstrap computations')
for fx,f in enumerate(flds):  #loop through algorithms 
    print(f'Beginning Bootstrap for {f}')
    for vx,v in enumerate(vb):  #loop through the variables 
        data = [] #hold data temporarily 
        n=1000
        for num_times in range(n):  #loop for bootstrap sample n times 
            
            df = eval(f)[vab][0][v]#data 
            #dtime = eval(f)[vab][0].time.values #get the time values 
            dtime = df.time.values
            dtim = np.random.choice(dtime,size=len(dtime))  #randomly sample the times 
              
            df = df.sel(time=dtim).mean('time') #select the randomly sampled time 
            df.coords['tstep'] = num_times #create timestamp for the data 
            data.append(df)   #
        concated = xr.concat(data,dim='tstep') #concatenate all the bootstrapped data 
        #eval(f'{f}_boot')[v].append(concated)
        concated.to_netcdf(f'Bootstrapped Data/{reg}{f}_non_consensus-boot-{v}.nc') #store as netcdf 
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
