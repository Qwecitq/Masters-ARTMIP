#!/usr/bin/env python3
# -*- coding: utf-8 -*-

exec(open('imports.py').read())

flds=['cascade_bard_v1','guan_waliser','mundhenk','reid250']
#calcualte anomalies for the non-consensus AR days for the variables

#calculate all the climatologies 
mpath = '/global/cscratch1/sd/kquagra/'  #main path 
vbls = ['IVT','TCWV','MSL','PV','Z','T'] #variable names 
file_locs = [mpath+f'era5_{x}/' for x in vbls ] #file locations 

climatologies = collections.defaultdict(list) #dict to hold climatologies 
for ix,x in enumerate(file_locs):
    print(f'Beginning to compute climatology for {vbls[ix]}')
    ds=xr.open_mfdataset(f'{x}*.nc',parallel=True)
    dmn = ds.sel(time=ds.time.dt.month.isin([1,2,12])).mean('time') ; print('Done Loading')
    print('Converting Longitudes to -180 to 180')
    dmn.coords['longitude'] = (dmn.coords['longitude']  + 180) % 360 - 180 #convert from 0-360 to -180 to 180
    dmn=dmn.sortby(dmn.longitude) #sort the lons 
    dmn=dmn.sel(longitude = slice(-170,-30), latitude = slice(60,10))
    
    print('Computing anomalies')
    for fx,f in enumerate(flds):
            
        dd = xr.open_mfdataset(f'{vbls[ix]}_non_consensus_point/{f}*{vbls[ix]}*.nc',parallel=True, chunks={'time':1}) ; print('Done loading data')
        anom = dd - dmn ; print('Done calculating anomaly')
        anom=anom.mean('time') ; print('Done calculating time mean')
        anom.to_netcdf(f'{mpath}ARs Work/Non-coinciding ARDTs/Climatologies/ERA5-Climatology/era5_clim-{f}-{vbls[ix]}.nc') ;print('Done Saving Data')
        #climatologies[vbls[ix]].append(ds) #compute the climatology and append to the didctionary 
    print(f'Done computing climatology for {vbls[ix]}')