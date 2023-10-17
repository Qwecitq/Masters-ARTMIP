#!/usr/bin/env python3
# -*- coding: utf-8 -*-


exec(open('imports.py').read())
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import saver

reg = 'west_europe'

fld = sys.argv[1] #ARDT under consideration
vb = sys.argv[2] #variable under consideration 
which_type = sys.argv[3] #type of transect
k = 4

if vb == 'pot':
    VB = 'PT'
else:
    VB = vb.upper()


 #receive input for which type of transect is done 

#change from zonal to meridional (lon to lat) based on the type of transect, either NS or EW
if which_type == 'NS':
    tloc = 'latitude'
elif which_type =='EW':
    tloc = 'longitude'
    

print('Loading Data')
if fld == 'cascade_bard_v1':
    #calculate the anomaly of the data based on the type of transect 
    file_loc = glob.glob(f'{reg}_{vb}/{fld}_{vb}/*.nc',recursive=True)
    print(file_loc)
    data = xr.open_mfdataset(file_loc)
    data.coords['longitude'] = (data['longitude']+180)%360 - 180
    data = data.sortby(data.longitude); data = data.sel(longitude=slice(-15,30),latitude=slice(80,35))
        
    data = data #- clim.mean(tloc).mean('time') #anomaly 
    print(data)
    
else:
    #obtain file location 
    file_loc =  glob.glob(f'{reg}_{vb}/{fld}_{vb}/*.nc', recursive=True)
    print(file_loc)
    #laod data and calculate anomaly based on the type of transect
    data = [xr.open_dataset(ds,chunks={'time':1,'level':37}) for ds in file_loc]
    data = xr.concat(data,dim = 'time')
    data.coords['longitude'] = (data['longitude']+180)%360 - 180
    data = data.sortby(data.longitude); data = data.sel(longitude=slice(-15,30),latitude=slice(80,35))
    
    data = data# - clim.mean(tloc).mean('time') #anomaly 
    print(data)
    

print(f'Selecting Transects for {which_type}')
if which_type == 'across' :

    #use the metpy option for this kind of transect for along 2 different x,y points 
    start_point = (67,-15) #start point for cs
    end_point = (67,28) #end point for cs

    #select cross-section by squeezing data and selecting the cross section
    trans_data = data[VB].metpy.parse_cf().squeeze()
    trans_data = cross_section(trans_data,start_point, end_point)

#this cross-section uses the selection based on latitude then, longitude
elif which_type=='EW':
    center_lat = 67 #65 #latitudinal center of the system
    center_lon = 8# longitudinal center of the system
    lon_span = 21.5 # degrees east and west of the center
    
    # extract the desired latitude
    trans_data = data[VB].sel(latitude=center_lat, method='nearest').sel(
        longitude = slice(center_lon - lon_span , center_lon + lon_span ))
           
elif which_type=='NS':
    center_lat = 55 #
    center_lon = 15# 12 # + 360 # Bloomington, IN
    lat_span = 20 # degrees
    
    # extract the desired latitude
    trans_data = data[VB].sel(longitude=center_lon, method='nearest').sel(
        latitude = slice( center_lat + lat_span,center_lat - lat_span ))
    

    
scaler = StandardScaler()
new_data = np.reshape(trans_data.values,(len(trans_data.time.values),len(trans_data.level)*len(trans_data[tloc].values)))
scaled_data = scaler.fit_transform(new_data)

#k = 3 # The number of clusters
kmeans = KMeans(n_clusters=k, random_state=0)
clustered_data = kmeans.fit_predict(scaled_data)



trans=xr.Dataset({VB:trans_data})

trans['cluster'] = kmeans.fit_predict(scaled_data)

trans=trans[VB].expand_dims({'cluster' :trans.cluster.values})
trans = trans.assign_coords({'cluster':trans.cluster.values}) 

trans = xr.Dataset({f'{VB}_trans':trans})
print(trans)

sv_path = f'k_means_{vb}'

if os.path.isdir(sv_path) == False:
    os.mkdir(sv_path)
    
ds_sv_path = f'{sv_path}/{fld}_{vb}/'
if os.path.isdir(ds_sv_path) == False:
    os.mkdir(ds_sv_path)

saver.ch_paths(trans,ds_sv_path,f'e5.{vb}.{fld}.kmeans')