#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''THIS SCRIPT SELECTS THE LANDFALLING REGIONS FOR ATMOSPHERIC RIVERS FOR THE INPUT ALGORITHMS'''

exec(open('imports.py').read())

#### CHANGE ME ######
reg='south_america'
#input details from user
ardt = str(sys.argv[1])                                               #name of ardt 

blat = sys.argv[2] ; blat = int(blat)                             #bottom latitude
llon = sys.argv[3] ; llon = int(llon)                              #left longitude 

threshold = sys.argv[4] ; threshold = float(threshold)                             #set the threshold for the selection of landfall (usually between 0 and 1)

seas = sys.argv[5]                                         #this selects the seasons or the period the user wants to select for landfalling ARs
seas = seas.split(',')


#select the landfall region
lfreg = sys.argv[6]
lfreg = lfreg.split('/')                                       #split the string input by ','
lfreg = [float(x) for x in lfreg]                        #convert the string numbers into floats

sv_path=sys.argv[7]                                    #location to save data


#############################################################################################################
#############################################  RECEIVED INPUTS  ##############################################
#############################################################################################################
#check if all inputs are right 

print(f'ARDT name: {ardt} \nBottom Latitude: {blat} \nLeft Longitude: {llon} \nThreshold: {threshold} \nSeasons: {seas} \nSel Region Lon/Lat: {lfreg} ')



#############################################################################################################
########################################## Tunable Parameters #################################################
#############################################################################################################

#path to the dataset for selection
mpath = '/global/cfs/projectdirs/m1517/cascade/external_datasets/ARTMIP/tier1/ftp_mirror/' 

#You can change this to suit what you want 

lat_width = 40                                    #region latitude width
lon_width = 60                                   #region longitude width 

selection_box_width = 2                    #landfall region box width in both lon and lat directions

#set the latitude and longitude range
tlat = blat + lat_width
rlon = llon + lon_width

#############################################################################################################
############################################ Select Season ####################################################
#############################################################################################################


#assign the calendar months into numbers
months_in_year = [calendar.month_abbr[m] for m in np.arange(1,13)] 
months_idx = {mon:idx+1 for idx,mon in enumerate(months_in_year)}

#select the period the user wants to select landfall for
period = [months_idx.get(se) for se in seas]
    

#############################################################################################################
############################################## Load data ######################################################
#############################################################################################################


#load files 
files= glob.glob(f'{mpath}{ardt}/*.nc*',recursive=True)
print(files)
data = xr.open_mfdataset(files,parallel=True).sel(lon=slice(llon,rlon), lat=slice(blat,tlat))
data = data.isel(time=data.time.dt.month.isin(period))
if ardt=='guan_waliser':
    data = data['ar_binary_tag']
    ds = xr.Dataset()
    
    #for guan and waliser, you need to select the time slices in 10 years interval manually or you can just attach a for loop to do that!
    ds['ar_binary_tag'] = data.sel(time=slice('1980','1989')) #rechunk data time axis to be optimal
    data = ds.chunk({'time':728})
 #print out data to inspect the selection is correct 
print(data)

#############################################################################################################
######################################## Set Threshold for selection  ##############################################
#############################################################################################################

#select the region for the landfall test 

selection_box = data.sel(lon=slice(lfreg[0], lfreg[0]+selection_box_width),
                        lat=slice( lfreg[1],  lfreg[1]+selection_box_width))

print(selection_box)
#check if the threshold criteria is met within the selected region 

sel_box=selection_box.chunk({'time':1})
selected_times= sel_box.time.where(sel_box.mean(['lon','lat']) >= threshold,drop=True)['time'].values

print(len(selected_times))
selected_data = data.sel(time=selected_times)

#############################################################################################################
########################################### Begin Saving Data  #################################################
#############################################################################################################

print(f'Started saving for {ardt} ...')
if os.path.isdir(sv_path)==False:
    os.mkdir(sv_path)

final_path = f'{sv_path}/{reg}_{ardt}_phd/'

if os.path.isdir(final_path)==False:
    os.mkdir(final_path)
    
saver.ch_paths(selected_data, f'{final_path}',f'{ardt}-{reg}-ARs.gt.8','ar_binary_tag')
print('Done ...')