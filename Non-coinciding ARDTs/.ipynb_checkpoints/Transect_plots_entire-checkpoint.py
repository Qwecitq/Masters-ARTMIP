#!/usr/bin/env python3
# -*- coding: utf-8 -*-

exec(open('imports.py').read())
from metpy.interpolate import cross_section

reg = 'west_europe_'
vbls = ['pv','z','pot']
flds = ['cascade_bard_v1','guan_waliser','mundhenk_v2','reid250']
#flds = ['guan_waliser']#,'reid250']

cascade_bard_v1_data=collections.defaultdict(list)
guan_waliser_data = collections.defaultdict(list)
mundhenk_v2_data = collections.defaultdict(list)
reid250_data = collections.defaultdict(list)
#load temperature dataset and calculate anomaly
for f in flds:
    for vb in vbls:
        print(f' {f} {vb}')

        if f == 'guan_waliser':
            
            files = glob.glob(f'{reg}{vb}/{f}_{vb}/*20*.nc', recursive = True)
            print(files[0]) #check for files in path 
            one = xr.open_mfdataset(files,parallel=True).sel(time=slice('2000-01','2018-01'))
            
            files = glob.glob(f'{reg}{vb}/{f}_{vb}/*19*.nc', recursive = True)
            print(files[0])#check for files in path 
            two = xr.open_mfdataset(files,parallel=True).sel(time=slice('1980-01','2000-01'))
            eval(f'{f}_data')[vb] = xr.concat([two,one],dim='time')
            
        else:
            files = glob.glob(f'{reg}{vb}/{f}_{vb}/*.nc', recursive = True)
            print(files[0])#check for files in path
            eval(f'{f}_data')[vb]=xr.open_mfdataset(files,chunks={'time':1},parallel=True)

        

#start_point = (0,58)
#end_point = (28,75)
#collection to hold data
teca_tr = collections.defaultdict(list)
gw_tr = collections.defaultdict(list)
mdk_tr = collections.defaultdict(list)
rd_tr = collections.defaultdict(list)

tr_ardts = ['teca_tr','gw_tr','mdk_tr','rd_tr']
flds = ['cascade_bard_v1','guan_waliser','mundhenk_v2','reid250']
#flds = flds[1:]
vbls = ['pot','pv','z']
vnms = ['PT','PV','Z']

#set the points for the line
start_point = (45,15) #start point for cs
end_point = (78,15) #end point for cs

#select cross-section by squeezing data and selecting the cross section
for tx,trs in enumerate(tr_ardts):

    for vx,vb in enumerate(vbls):

        dt = eval(f'{flds[tx]}_data')[vb].metpy.parse_cf().squeeze()
        eval(trs)[vb] =  cross_section(dt[vnms[vx]],start_point,end_point)#.set_coords(('lat', 'lon'))

vbls = ['pot','pv','z']
vnms = ['PT','PV','Z']
tr_ardts = ['mdk_tr','gw_tr','rd_tr','teca_tr']
flds = ['mundhenk_v2','guan_waliser','reid250','cascade_bard_v1']
flds = flds[1:]


intervals = [np.arange(-8,8.1,0.5), np.arange(-0.4,0.41,0.1), np.arange(-150,151,10) ] #set intervals for plot contours
cmps = ['coolwarm','bwr','bwr']

#ct=0 #set a counter for the plots
for trs in tr_ardts[1:]:
    ct = 0
    
    fig,axes = plt.subplots(ncols=3, figsize=(16, 4))

    for vb,ax,iv in zip(vbls,axes.flat,intervals):

        eval(f'{trs}')[vb]['index']= eval(f'{trs}')[vb].latitude #reassign the index to longitude
        
        if vb == 'z':
            ds =  (eval(f'{trs}')[vb]/(10)).mean('time') - (eval(f'{tr_ardts[0]}')[vb]/(10)).mean('time')
            ds.plot.contourf(ax=ax, cmap= cmps[ct], levels = iv, cbar_kwargs=dict(orientation='horizontal'))
            ax.set_ylim(1000,100)

        elif vb == 'pv':
            ds =  (eval(f'{trs}')[vb]/(10e-6)).mean('time') - (eval(f'{tr_ardts[0]}')[vb]/(10e-6)).mean('time')
            ds.plot.contourf(ax=ax, cmap= cmps[ct], levels = iv, cbar_kwargs=dict(orientation='horizontal'))
            ax.set_ylim(1000,100)
        
        else:
            
            #print(eval(f'{trs}')[vb].mean('time'))
            print(eval(tr_ardts[0])[vb])
            ds =  eval(f'{trs}')[vb].mean('time') - eval(tr_ardts[0])[vb].mean('time')
            ds.plot.contourf(ax=ax, cmap= cmps[ct], levels = iv, cbar_kwargs=dict(orientation='horizontal'))

        ax.set_ylim(1000,100)
        ax.set_title(f' Composite Average {vb.upper()}')
        ct+=1
    #plt.show()
    plt.savefig(f'transect images/{trs}_gw_crt.png')
