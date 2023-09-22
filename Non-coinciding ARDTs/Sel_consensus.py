#!/usr/bin/env python3
# -*- coding: utf-8 -*-

exec(open('imports.py').read())
reg = 'west-europe-'
main_path = '/global/cscratch1/sd/kquagra/ARs Work/larger_region_AR_capture/west-europe-' #main path
flds = ['/cascade_bard_v1/','/guan_waliser/','/mundhenk/','/reid250/'] #algorithms 
#ARDTs = single_data(main_path,flds,'landfall_ARs__0.9-0.4.nc')
ARDTs = ARDTs_load(main_path,flds,'west-europe-landfall_ARs__0.8-1.0') #load main AR data c


print(ARDTs['mundhenk'][0])#.time.values)
#SELECT DAYS THAT APPEAR IN ALL THE ALGORITHMS 
consensus_days = reduce(np.intersect1d,(ARDTs['guan_waliser'][0].time.values,ARDTs['cascade_bard_v1'][0].time.values,ARDTs['reid250'][0].time.values,ARDTs['mundhenk'][0].time.values))

cons_ARDTs =   { } # specific consensus days when ARs are captured 
non_cons_ARDTs = { }# select the specific non_consensus days for AR captured 

for fx,f in enumerate(flds):
    f = f[1:-1]
    cons_ARDTs[f]=[] ;non_cons_ARDTs[f]=[] #empty list to hold data 
    
    print(ARDTs[f])
    ARDT = ARDTs[f][0]
    print(f'{f} \n {ARDT.time.values}')
    cons_ARDTs[f].append(ARDT.sel(time=np.isin(ARDT.time.values,consensus_days))) #compute the consensus days 
    non_cons_ARDTs[f].append(ARDT.sel(time = np.isin(ARDT.time.values,consensus_days,invert=True))) #invert to get non-consensus days 
    print(f'Done with {f}')
"""main_path = '/global/cscratch1/sd/kquagra/ARs Work/Final_variable_data/' #main path 
fpath = f'{reg}TCWV_ARS_capture_2_point/'
flds = ['/cascade_bard_v1/','/guan_waliser/','/mundhenk/','/reid250/'] #algorithms 
TCWV = load_data(flds,'TCWV',fpath, main_path=main_path) #load the data

TCWV_non_con = cons_data(TCWV,consensus_days,flds,'time',vbl_name='TCWV',sv_path=f'{reg}TCWV_non_consensus_point',invert=True)  #check for non-consensus days 
TCWV_con = cons_data(TCWV,consensus_days,flds,'time',vbl_name='TCWV',sv_path=f'{reg}TCWV_consensus_point',invert=False) #check for consensus days""" 

"""main_path = f'/global/cscratch1/sd/kquagra/ARs Work/Final_variable_data/{reg}PV_ARS_capture_2_point/' #main path 
fpath = ''
flds = ['/cascade_bard_v1/','/guan_waliser/','/mundhenk/','/reid250/'] #algorithms 
PV = load_data(flds,'PV',fpath, main_path=main_path) #load the data 

PV_non_con = cons_data(PV,consensus_days,flds,'time',vbl_name='PV',sv_path=f'{reg}PV_non_consensus_point',invert=True) #check for non-consensus days 
PV_con = cons_data(PV,consensus_days,flds,'time',vbl_name='PV',sv_path=f'{reg}PV_consensus_point',invert=False) #check for consensus days """

main_path = f'/global/cscratch1/sd/kquagra/ARs Work/Final_variable_data/{reg}GEOPTH_ARS_capture_2_point/' #main path 
fpath = ''
flds = ['/cascade_bard_v1/','/guan_waliser/','/mundhenk/','/reid250/'] #algorithms 
GEOPTH = load_data(flds,'GEOPTH',fpath, main_path=main_path)#load the data

GEOPTH_non_con = cons_data(GEOPTH,consensus_days,flds,'time',vbl_name='GEOPTH',sv_path=f'{reg}GEOPTH_non_consensus_point',invert=True) #check for non-consensus days
GEOPTH_con = cons_data(GEOPTH,consensus_days,flds,'time',vbl_name='GEOPTH',sv_path=f'{reg}GEOPTH_consensus_point',invert=False) #check for consensus days 

"""main_path = '/global/cscratch1/sd/kquagra/ARs Work/IVT_ARS_capture_2_point/' #main path 
fpath = ''
flds = ['/cascade_bard_v1/','/guan_waliser/','/mundhenk/','/reid250/'] #algorithms 
IVT = single_data(main_path,flds,'_IVT.nc') #load the data

IVT_non_con = cons_data(IVT,consensus_days,flds,'time',vbl_name='IVT',sv_path='IVT_non_consensus_point',invert=True)#check for non-consensus days
IVT_con = cons_data(IVT,consensus_days,flds,'time',vbl_name='IVT',sv_path='IVT_consensus_point',invert=False) #check for consensus days 
"""
"""main_path = f'/global/cscratch1/sd/kquagra/ARs Work/Final_variable_data/{reg}MSL_ARS_capture_2_point/' #main path 
fpath = ''
flds = ['/cascade_bard_v1/','/guan_waliser/','/mundhenk/','/reid250/'] #algorithms 
MSL = load_data(flds,'MSL',fpath, main_path=main_path)  #load the data

MSL_non_con = cons_data(MSL,consensus_days,flds,'time',vbl_name='MSL',sv_path=f'{reg}MSL_non_consensus_point',invert=True) #check for non-consensus days
MSL_con = cons_data(MSL,consensus_days,flds,'time',vbl_name='MSL',sv_path=f'{reg}MSL_consensus_point',invert=False)  #check for consensus days """

"""main_path = f'/global/cscratch1/sd/kquagra/ARs Work/Final_variable_data/{reg}TEMP_ARS_capture_2_point/' #main path 
fpath = ''
flds = ['/cascade_bard_v1/','/guan_waliser/','/mundhenk/','/reid250/']  #algorithms 
TEMP = load_data(flds,'TEMP',fpath, main_path=main_path) #load the data

TEMP_non_con = cons_data(TEMP,consensus_days,flds,'time',vbl_name='TEMP',sv_path=f'{reg}TEMP_non_consensus_point',invert=True) #check for non-consensus days
TEMP_con = cons_data(TEMP,consensus_days,flds,'time',vbl_name='TEMP',sv_path=f'{reg}TEMP_consensus_point',invert=False)  #check for consensus days """