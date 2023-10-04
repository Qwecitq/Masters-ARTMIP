#!/usr/bin/env python3
# -*- coding: utf-8 -*-


exec(open('imports.py').read())


#################################################################################################################
########################################  COLLECT THE INPUT VARIABLES TO RUN ###################################
#################################################################################################################

reg=sys.argv[1] #REGION
mpath = ''
vb = sys.argv[2] # VARIABLE NAME
fld = sys.argv[3] #ARDT NAME 
ff = sys.argv[4] # FILE TO BE REWRITTEN

print(f'Processing {ff} file')

#################################################################################################################
######################################## OPENNING AND RESAVING THE FILES #######################################
#################################################################################################################
dataset = xr.open_dataset(f'{ff}') 

#COLLECT DATA INFORMATION [TIME, LEVELS, LONGITUDE AND LATITUDE]
ln = len(dataset.longitude)
lt = len(dataset.latitude)
tm = len(dataset.time)
levs = len(dataset.level)
print(tm , levs, lt ,ln)

#SET ENCODING TO COMPRESS OR DEFLATE THE FILE SIZE
enc={f'{vb.upper()}' : {'zlib':True,'complevel': 5, 'fletcher32': True,'chunksizes': (1, levs, lt, ln)}}

#CREATE THE SAVING PATHS 
save_path = f'{reg}_{vb}_resaved/'

if os.path.isdir(save_path)== False:
    os.mkdir(save_path)

fld_save_path = f'{fld}_{vb}/'

if os.path.isdir(f'{save_path}{fld_save_path}') == False:
    os.mkdir(f'{save_path}{fld_save_path}')

pth = ff[len(f'{reg}_{vb}/{fld}_{vb}/'):]
print(pth)

#SAVE THE DATASET
print(f'Saving {ff}')
dataset.to_netcdf(f'{save_path}{fld_save_path}{pth}', format="NETCDF4",engine='netcdf4', encoding=enc)