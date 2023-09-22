#!/usr/bin/env python3
# -*- coding: utf-8 -*-


exec(open('imports.py').read())
main_path = '/global/pscratch/sd/k/kquagra/ARs Work/larger_region_AR_capture/west-europe-' #main path
flds = ['/cascade_bard_v1/','/guan_waliser/','/mundhenk/','/reid250/'] #algorithms 
#ARDTs = single_data(main_path,flds,'landfall_ARs__0.9-0.4.nc')
ARDTs = ARDTs_load(main_path,flds,'west-europe-landfall_ARs__0.8-1.0') #load main AR data c
#SELECT DAYS THAT APPEAR IN ALL THE ALGORITHMS 

consensus_days = reduce(np.intersect1d,(ARDTs['guan_waliser'][0].time.values,ARDTs['cascade_bard_v1'][0].time.values,ARDTs['reid250'][0].time.values,ARDTs['mundhenk'][0].time.values))

print('Loading ARTMIP data ')
cons_ARDTs =   { } # specific consensus days when ARs are captured 
non_cons_ARDTs ={ }# select the specific non_consensus days for AR captured 

for fx,f in enumerate(flds):
    f = f[1:-1]
    cons_ARDTs[f]=[] ;non_cons_ARDTs[f]=[] #empty list to hold data 
    
    print('Loading Consensus and Non-Consesus times')
    cons_ARDTs[f].append(ARDTs[f][0].sel(time=np.isin(ARDTs[f][0].time.values,consensus_days))) #compute the consensus days 
    non_cons_ARDTs[f].append(ARDTs[f][0].sel(time = np.isin(ARDTs[f][0].time.values,consensus_days,invert=True))) #invert to get non-consensus days 

print('Done loading ... \nBegin loading atmospheric data')    
mpath = '/global/cfs/cdirs/m3522/cmip6/ERA5/'
vbls = ['w']#,'t'] #['msl','tcwv','pv',
locs = ['e5.oper.an.pl/']#,'e5.oper.an.pl/'] #['e5.oper.an.sfc/','e5.oper.an.sfc/','e5.oper.an.pl/',]
dpaths = [glob.glob(f'{mpath}{y}*/*{x}*',recursive=True) for x,y in zip(vbls,locs)]

#Load DATA

for dx,dp in enumerate(dpaths):
    vb=vbls[dx].upper()
    
    if vb == 'Z' or vb =='T' or vb == 'PV' or vb == 'W':
        data = xr.open_mfdataset(dp,parallel=True).sel(level=500)
    else :
        data = xr.open_mfdataset(dp,parallel=True)

    non_cascade_bard_v1 = collections.defaultdict(list)
    non_guan_waliser = collections.defaultdict(list)
    non_mundhenk = collections.defaultdict(list)
    non_reid250 = collections.defaultdict(list)
    flds = ['/cascade_bard_v1/','/guan_waliser/','/mundhenk/','/reid250/'] #algorithms 
    

    for fx,f in enumerate(flds):
        f = f[1:-1]
        print(f'Working on {vb}')
        newd = data.sel(time=non_cons_ARDTs[f][0].time.values)
        eval(f'non_{f}')[vb].append(newd)

        newd = newd.chunk({'time':1})
        print('Saving data')
        SAVE_DIR=f'/global/pscratch/sd/k/kquagra/ARs Work/Final_variable_data/west-europe-{vb}_non_consensus_point/'
        if os.path.isdir(SAVE_DIR) == False:
            os.mkdir(SAVE_DIR)
        newd.to_netcdf(f'{SAVE_DIR}{f}False_concensus_{vb}.nc')
    print(f'Done working on {vb}')