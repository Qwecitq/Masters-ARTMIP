#!/usr/bin/env python
# -*- coding: utf-8 -*-

exec(open('imports.py').read())

#open the data paths 
vbl_nme='z'
geopth_path = '/global/cfs/projectdirs/m3522/cmip6/ERA5/e5.oper.an.pl/' #main path
tme = np.arange(1990,2020) #data span 


geopth=[] #list appended 

#loop to get the paths for the data 
for i in tme:
    geopth.append(glob.glob(geopth_path+str(i)+'*/*'+vbl_nme+'*.nc', recursive=True))
print('Done loading paths ...')
#join the lists into one list 
geopth=sum(geopth,[])

#select the level of interest from the data 
lvl = 700 #level 

#open the list comprehension 
print('Begin loading data')
geopth=xr.open_mfdataset(geopth,parallel=True,engine='h5netcdf').sel(level=lvl)
geopth = geopth.sel(time=geopth.time.dt.month==8)
print('Done')

#geopth=xr.concat(geopth,dim='time').sel(time=slice('1990','2020'),level=lvl) 

print('Fix the lons and lats')
geopth_data = geopth
geopth_data.coords['longitude'] = (geopth_data.coords['longitude']  + 180) % 360 - 180
geopth_data=geopth_data.sortby(geopth.longitude)
gpt=geopth_data.sel(longitude= slice(-5,5), latitude=slice(8,4))  #select the region

gptheight = gpt/9.8

print('Beginning to resample data to 1D')
gptheight = gptheight.resample(time='1D').mean()

print('Begin to plot')
fig,axes = plt.subplots(ncols = 1, figsize = (7,4.5),subplot_kw={'projection':ccrs.PlateCarree()})
p1=gptheight['Z'].mean('time').plot.contour(ax=axes,add_colorbar=True,cbar_kwargs = {'orientation':'horizontal','label':'Geopotential Heights'})
cart_plot(axes)
plt.savefig('Ghana_gpt_1990-2020.png',facecolor = 'white')
print('Done plotting')