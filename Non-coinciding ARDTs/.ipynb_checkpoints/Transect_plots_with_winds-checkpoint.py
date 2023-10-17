#!/usr/bin/env python3
# -*- coding: utf-8 -*-

exec(open('imports.py').read())
from metpy.interpolate import cross_section

reg = 'west_europe_'
vbls = ['pv','q','pot']#,'u','v']
cl_vbls = ['pv','q','pt']
#################################################################################################################
################################ TAKE INPUT FOR FIELDS AND ARDT TRANSECT NAMES ################################
#################################################################################################################
flds=[sys.argv[1]]
tr_ardts=[sys.argv[2]]

which_type = sys.argv[3] #receive input for which type of transect is done 

#change from zonal to meridional (lon to lat) based on the type of transect, either NS or EW
if which_type == 'NS':
    tloc = 'latitude'
elif which_type =='EW':
    tloc = 'longitude'


cascade_bard_v1_data=collections.defaultdict(list)
guan_waliser_data = collections.defaultdict(list)
mundhenk_v2_data = collections.defaultdict(list)
reid250_data = collections.defaultdict(list)

#################################################################################################################
################################ LOAD TEMPERATURE DATASET AND CALCULATE ANOMALY  ############################
#################################################################################################################
#load temperature dataset and calculate anomaly
mpath = '/global/cfs/cdirs/m3522/cmip6/ERA5/'
locs ='e5.oper.an.pl/'
dl = glob.glob(f'{mpath}{locs}*/*pv*.nc',recursive=True)

#with xr.open_mfdataset( glob.glob(f'{mpath}{locs}*/*pv*.nc',recursive=True) ,chunks={'time':1}) as nds:
#    clim = nds.mean('time')
#    clim.to_netcdf('PV_climatology_1980_2020.nc')

#open the surface pressure dataset and select the region 
surf_pressure = xr.open_mfdataset(f'{mpath}e5.oper.an.sfc/200001/*sp*.nc')
surf_pressure.coords['longitude'] = (surf_pressure.coords['longitude']  + 180) % 360 - 180 #convert from 0-360 to -180 to 180
surf_pressure=surf_pressure.sortby(surf_pressure.longitude) #sort the lons
surf_pressure = surf_pressure.sel(longitude=slice(-15,30),latitude=slice(80,35))


for f in flds:
    for vb,cvb in zip(vbls,cl_vbls):
        print(f' {f} {vb} {cvb}')
        
        #To obtain this, use the Climatologies.sh script in the Bash Runs folder to get the climatologies from the ERA5 hub
        clim = xr.open_mfdataset(f'Bash Runs/{cvb.upper()}_climatologies_west_europe/*.nc')
        clim = clim.sel(time=clim.time.dt.month.isin([1,2,12])).mean(tloc)
        
        if f != 'cascade_bard_v1':
                
            files = glob.glob(f'{reg}{vb}/{f}_{vb}/*.nc', recursive = True)
            #print(files[0]) #check for files in path 
            #one = xr.open_mfdataset(files,parallel=True,chunks={'time':1,'level':37}).sel(time=slice('1980-01','2000-01'))
            
            
            #files = glob.glob(f'{reg}{vb}/{f}_{vb}/*20*.nc', recursive = True)
            #print(files[0])#check for files in path 
            #two = xr.open_mfdataset(files,parallel=True,chunks={'time':1,'level':37}).sel(time=slice('2000-01','2018-01'))
            #eval(f'{f}_data')[vb]  = xr.concat([two,one],dim='time') 
            one  =  [xr.open_dataset(f,chunks={'time':1,'level':37}) for f in files]
            print('Done Loading Files \nConcatenating FIles')
            eval(f'{f}_data')[vb]  = xr.concat(one,dim='time')
            
            
            #subtract the longitudinal climatology from the data
            if vb == 'pv' or vb =='pot': #add variables that you want the anomaly to be calculated
                
                eval(f'{f}_data')[vb] = eval(f'{f}_data')[vb] - clim               
                
            #elif vb == 'pot':
               
                open_ds_main = eval(f'{f}_data')[vb]
                open_ds_main.coords['longitude'] = (open_ds_main.coords['longitude']  + 180) % 360 - 180 #convert from 0-360 to -180 to 180
                open_ds_main=open_ds_main.sortby(open_ds_main.longitude) #sort the lons
                eval(f'{f}_data')[vb] = open_ds_main.sel(longitude=slice(-15,30),latitude=slice(80,35))
                eval(f'{f}_data')[vb] = eval(f'{f}_data')[vb] #- clim
                
            else:
                
                #commenting this for the case of specific humidity.
                eval(f'{f}_data')[vb] = eval(f'{f}_data')[vb] #- clim    
                open_ds_main = eval(f'{f}_data')[vb]
                open_ds_main.coords['longitude'] = (open_ds_main.coords['longitude']  + 180) % 360 - 180 #convert from 0-360 to -180 to 180
                open_ds_main=open_ds_main.sortby(open_ds_main.longitude) #sort the lons
                eval(f'{f}_data')[vb] = open_ds_main.sel(longitude=slice(-15,30),latitude=slice(80,35))
                print(eval(f'{f}_data')[vb])
                
        elif f == 'cascade_bard_v1':
            files = glob.glob(f'{reg}{vb}/{f}_{vb}/*.nc', recursive = True)
            print(files[0])#check for files in path
            eval(f'{f}_data')[vb]=xr.open_mfdataset(files,chunks={'time':1,'level':37},parallel=True)
            
            
            if vb == 'pv' or vb=='pot':
                                
                eval(f'{f}_data')[vb] = eval(f'{f}_data')[vb] - clim
                
            #elif vb == 'pot':

                open_ds_main = eval(f'{f}_data')[vb]
                open_ds_main.coords['longitude'] = (open_ds_main.coords['longitude']  + 180) % 360 - 180 #convert from 0-360 to -180 to 180
                open_ds_main=open_ds_main.sortby(open_ds_main.longitude) #sort the lons
                eval(f'{f}_data')[vb] = open_ds_main.sel(longitude=slice(-15,30),latitude=slice(80,35))
                eval(f'{f}_data')[vb] = eval(f'{f}_data')[vb] #- clim
                
            else:
                eval(f'{f}_data')[vb] = eval(f'{f}_data')[vb] #- clim
            #if vb!='pot':
                open_ds_main = eval(f'{f}_data')[vb]
                print(open_ds_main)
                open_ds_main.coords['longitude'] = (open_ds_main.coords['longitude']  + 180) % 360 - 180 #convert from 0-360 to -180 to 180
                open_ds_main=open_ds_main.sortby(open_ds_main.longitude) #sort the lons
                eval(f'{f}_data')[vb] = open_ds_main.sel(longitude=slice(-15,30),latitude=slice(80,35))
                print(eval(f'{f}_data')[vb])
            
vbls = ['u','w']

for f in flds:
    for vb in vbls:
        print(f' {f} {vb}')
        files = glob.glob(f'{reg}{vb}/{f}_{vb}/*.nc', recursive = True)
        print(files[0])#check for files in path
        eval(f'{f}_data')[vb]=xr.open_mfdataset(files,chunks={'time':1,'level':37}) #parallel=True
        open_ds_main=eval(f'{f}_data')[vb]
        open_ds_main.coords['longitude'] = (open_ds_main.coords['longitude']  + 180) % 360 - 180 #convert from 0-360 to -180 to 180
        open_ds_main=open_ds_main.sortby(open_ds_main.longitude) #sort the lons
        eval(f'{f}_data')[vb] = open_ds_main.sel(longitude=slice(-15,30),latitude=slice(80,35))
        print(eval(f'{f}_data')[vb])
        
#################################################################################################################
############################################ BEGIN TRANSECT SECTION ###########################################
#################################################################################################################

from metpy.interpolate import cross_section

#start_point = (0,58)
#end_point = (28,75)
#collection to hold data
teca_tr = collections.defaultdict(list)
gw_tr = collections.defaultdict(list)
mdk_tr = collections.defaultdict(list)
rd_tr = collections.defaultdict(list)

#tr_ardts = ['mdk_tr']# ['teca_tr']# ['gw_tr']# ['rd_tr']#
#flds =['mundhenk_v2']#['cascade_bard_v1']#['guan_waliser']#['reid250']#
#flds = flds[1:]

#receive input from shell for the field names (ARDTS) and the transect ARDT. variable name
flds=[sys.argv[1]]
tr_ardts=[sys.argv[2]]
vbls = ['pot','pv','q','u','w']
vnms = ['PT','PV','Q','U','W']

#set the points for the line (lat, lon)
#start_point = (45,0) #start point for cs
#end_point = (78,22) #end point for cs


if which_type == 'across' :

    #use the metpy option for this kind of transect for along 2 different x,y points 
    start_point = (67,-15) #start point for cs
    end_point = (67,28) #end point for cs

    #select cross-section by squeezing data and selecting the cross section
    for tx,trs in enumerate(tr_ardts):

        for vx,vb in enumerate(vbls):

            dt = eval(f'{flds[tx]}_data')[vb].metpy.parse_cf().squeeze()
            eval(trs)[vb] =  cross_section(dt[vnms[vx]],start_point,end_point)#.set_coords(('lat', 'lon'))
            
            sp_dt = surf_pressure['SP'].metpy.parse_cf().squeeze()
            sp_dt = cross_section(sp_dt,start_point,end_point)

#this cross-section uses the selection based on latitude then, longitude
elif which_type=='EW':
    center_lat = 67 #65 #latitudinal center of the system
    center_lon = 8# longitudinal center of the system
    lon_span = 21.5 # degrees east and west of the center
    for tx,trs in enumerate(tr_ardts):

        for vx,vb in enumerate(vbls):
            # extract the desired latitude
             eval(trs)[vb] =eval(f'{flds[tx]}_data')[vb][vnms[vx]].sel(latitude = center_lat, method='nearest').sel(
                 longitude = slice(center_lon - lon_span , center_lon + lon_span ))
    
    sp_dt = surf_pressure.sel(latitude = center_lat, method='nearest').sel(
                longitude = slice(center_lon - lon_span , center_lon + lon_span ))
    
    cent = center_lat  #set the center for conversion to degrees/s
    
#this cross-section uses the selection based on longitude then, latitude            
elif which_type=='NS':
    center_lat = 55 #
    center_lon = 15# 12 # + 360 # Bloomington, IN
    lat_span = 20 # degrees
    for tx,trs in enumerate(tr_ardts):

            for vx,vb in enumerate(vbls):
                # extract the desired latitude
                eval(trs)[vb] =eval(f'{flds[tx]}_data')[vb][vnms[vx]].sel(longitude = center_lon, method='nearest').sel(
                    latitude = slice( center_lat + lat_span,center_lat - lat_span ))
                
    sp_dt = surf_pressure.sel(longitude = center_lon, method='nearest').sel(
                    latitude = slice( center_lat + lat_span,center_lat - lat_span ))
    
    cent = center_lon #set the center for conversion to degrees/s
vbls = ['pot','q','pv']#,'u','v']
vnms = ['PT','Q','PV']#,'U','V']
#tr_ardts = ['mdk_tr']# ['teca_tr']# ['gw_tr']# ['rd_tr']#
#flds =['mundhenk_v2']#['cascade_bard_v1']#['guan_waliser']#['reid250']#
flds=[sys.argv[1]]
tr_ardts=[sys.argv[2]]

#intervals = [np.arange(-4,4.5,0.5),  np.arange(-0.4,0.41,0.02), np.arange(-12,12,1) ] #set intervals for plot contours
from matplotlib.colors import LinearSegmentedColormap

colors = [(0, 'blue'), (0.2, 'cornflowerblue'), (0.5, 'bisque'), (0.8, 'orange'), (1, 'red')]
pot_cmap = LinearSegmentedColormap.from_list("custom_colormap", colors)
intervals = [np.arange(-7,7.1,.5), np.arange(0,7.1,0.5), np.arange(-0.3,0.31,0.04)]
cmps = ['bwr','YlOrRd','bwr']

for trs in tr_ardts:
    ct = 0
    
    fig,axes = plt.subplots(ncols=3, figsize=(16, 4))

    for vb,ax,ivs in zip(vbls,axes.flat,intervals):
        
        #check for the type of transect input and assign index or leave as it is. 
        if which_type!= 'EW' or which_type!='NS':
            
            eval(f'{trs}')[vb]['index']= eval(f'{trs}')[vb].longitude #reassign the index to longitude
            u = eval(f'{trs}')['u'].mean('time') 
            u['index'] = u.longitude
            #v = eval(f'{trs}')['v'].mean('time') 
            # convert wind from m/s to degrees/s
            Rearth = 6371000 # radius of the earth in meters
            u = u * 2*np.pi / (Rearth * np.cos(np.deg2rad(cent))) * 180/np.pi

            # convert omega from Pa/s to mb/s
            w_mbps = eval(f'{trs}')['w'] / 100
            print(w_mbps)
        else:
            
            #obtain the horizontal windspeed and convert from m/s to degrees/s
            u = eval(f'{trs}')['u'].mean('time') 
            Rearth = 6371000 # radius of the earth in meters
            u = u * 2*np.pi / (Rearth * np.cos(np.deg2rad(cent))) * 180/np.pi

            # obtain vertical windspeed and convert omega from Pa/s to mb/s
            w_mbps = eval(f'{trs}')['w'] / 100
            print(w_mbps)
            
        
        if vb == 'z':    
            
            ds =  (eval(f'{trs}')[vb]/(10*100)).mean('time') 
            ds.plot.contourf(ax=ax, cmap= cmps[ct], levels=ivs, cbar_kwargs=dict(orientation='horizontal'))
           
        if vb == 'q':    
            
            ds =  (eval(f'{trs}')[vb]*1000).mean('time')
            ds.plot.contourf(ax=ax, cmap= cmps[ct], levels=ivs, cbar_kwargs=dict(orientation='horizontal'))
                             
        elif vb == 'pv':
            
            ds =  (eval(f'{trs}')[vb]/(10e-6)).mean('time') #- (eval(f'{tr_ardts[0]}')[vb]/(10e-6)).mean('time')
            ds.plot.contourf(ax=ax, cmap= cmps[ct], levels=ivs, cbar_kwargs=dict(orientation='horizontal'))
          
            
        elif vb == 'pot':
            
            ds =  (eval(f'{trs}'))[vb].mean('time') #- 273.15 #- (eval(f'{tr_ardts[0]}')[vb]/(10e-6)).mean('time')
            
            print(ds)
            ds.plot.contourf(ax=ax, cmap= cmps[ct], levels=ivs, cbar_kwargs=dict(orientation='horizontal'))
           
        else:
            
            #print(eval(f'{trs}')[vb].mean('time'))
            print(eval(tr_ardts[0])[vb])
            ds =  eval(f'{trs}')[vb].mean('time') #- eval(tr_ardts[0])[vb].mean('time')
            ds.plot.contourf(ax=ax, cmap= cmps[ct], levels=ivs, cbar_kwargs=dict(orientation='horizontal'))
           
        
        
        

        # plot the wind vectors
        
        #add the quiver plots to all axes 
        ax.quiver(
            u[tloc][::3], # x-coordinates of the vectors
            u.level[::3], # y-coordinates of the vectors
            u[::3,::3], # u-component of the wind (in units of the transect)
            -w_mbps.mean('time')[::3,::3], # vertical velocity (in units of the transect [add the negative sign b/c the axes are reversed])
        )
        
        ax.plot(sp_dt['SP'][tloc].values,(sp_dt['SP']/100).mean('time').values,color='grey')
        ax.invert_yaxis()
        x=sp_dt['SP'][tloc].values
        y=(sp_dt['SP']/100).mean('time').values
        ax.fill_between(x, y, 1000, interpolate=True, color='black', alpha=1)
        
        
        ax.set_yticklabels(np.arange(1000, 50, -100))
        ax.set_ylim(1000, 100)
        ax.set_yticks(np.arange(1000, 50, -100))
        ax.set_title(f' Composite Average {vb.upper()}')
        ct+=1
    #plt.show()
        plt.savefig(f'transect images/{trs}_uv_{which_type}.png')

print('All processes completed')