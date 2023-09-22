#!/usr/bin/env python3
# -*- coding: utf-8 -*-

exec(open('imports.py').read())


main_path = '/global/pscratch/sd/k/kquagra/ARs Work/larger_region_AR_capture/west-europe-' #main path
flds = ['/cascade_bard_v1/','/guan_waliser/','/mundhenk/','/reid250/'] #algorithms 
#ARDTs = single_data(main_path,flds,'landfall_ARs__0.9-0.4.nc')
ARDTs = ARDTs_load(main_path,flds,'west-europe-landfall_ARs__0.8-1.0') #load main AR data 

flds =['guan_waliser','reid250','mundhenk_v2', 'cascade_bard_v1']
for ff in flds:
    print(f'{ff} : {len(ARDTs[ff][0].time)}')