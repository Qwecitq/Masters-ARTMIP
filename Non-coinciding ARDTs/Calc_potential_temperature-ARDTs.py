#!/usr/bin/env python3
# -*- coding: utf-8 -*-

exec(open('imports.py').read())
from metpy.interpolate import cross_section


print('Welcome to this program \nThis program computes the  potential temperature over a specific transect')
def pot(T, P, P_0=1016, Rd=0.286):
    """This function is used to compute the potential temperature of an air parcel
    
    Argument
    -------------
    T : array of temperature values for different levels 
    P : data/ air parcel's pressure levels 
    P_0: reference pressure level (default `1016`)
    Rd : constant for R/Cp
    """
    res = np.zeros_like(T)
    
    for ip, p  in enumerate(P):
        res[:,ip,:,:]  = T[:,ip,:,:] * ( (P_0/ p )**Rd )
        
    #res = T * ((P_0/P) ** Rd)
    
    return res 


reg = 'west_europe_'
vbls = ['t']
#flds = ['cascade_bard_v1','guan_waliser','mundhenk_v2','reid250']
flds = ['guan_waliser']#,'reid250']

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
            
            
#collection to hold data 
teca_tr = collections.defaultdict(list)
gw_tr = collections.defaultdict(list)
mdk_tr = collections.defaultdict(list)
rd_tr = collections.defaultdict(list)

#tr_ardts = ['teca_tr','gw_tr','mdk_tr','rd_tr']
tr_ardts = ['gw_tr']#,'rd_tr']

#flds = ['cascade_bard_v1','guan_waliser','mundhenk_v2','reid250']
flds = ['guan_waliser']#,'reid250']

vbls = ['t']

#divide the data into 3 and compute the potential temperature and save

tdiv = ["1980",'1996','2005']
tdiv2 = ['1995','2004','2017']


print('Calculating Potential Temperature, Data Manipulation and Saving')
for rx, r in enumerate(flds):
    
    #for ix,x in enumerate(tdiv):
    for ix in np.arange(1980,2018,2):
        iy = ix +1
    
    
        print(f'Calculating potential temperature for {ix} to {iy}.')
        
        trim_data = eval(f'{flds[rx]}_data')['t'].sel(time=slice(str(ix),str(iy)))
        print(trim_data.time.values)
        ds = pot(trim_data['T'].values, trim_data.level.values)

        print('Creating Data Array')
        newdt = xr.DataArray(ds,
                            dims= ['time','level','latitude','longitude'],
                            coords = [trim_data.time,
                                      trim_data.level,
                                      trim_data.latitude,
                                      trim_data.longitude])
        print('Converting to Dataset')
        dataset = xr.Dataset()
        dataset['PT'] = newdt

        print('Saving Data')
        dataset.to_netcdf(f'west_europe_pot/{r}_pot-temp_{ix}-{iy}.nc')
