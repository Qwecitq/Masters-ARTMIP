#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''THIS SCRIPT SELECTS THE LANDFALLING REGIONS FOR ATMOSPHERIC RIVERS FOR THE INPUT ALGORITHMS'''

exec(open('imports.py').read())
import dask

#### CHANGE ME ######
regs=['south_africa', 'south_america','north_america']

for reg in regs: 
    #input details from user
    #ardt = str(sys.argv[1])                                               #name of ardt 

    blat = sys.argv[1] ; blat = int(blat)                             #bottom latitude
    llon = sys.argv[2] ; llon = int(llon)                              #left longitude 

    threshold = sys.argv[3] ; threshold = float(threshold)                             #set the threshold for the selection of landfall (usually between 0 and 1)

    seas = sys.argv[4]                                         #this selects the seasons or the period the user wants to select for landfalling ARs
    seas = seas.split(',')


    #select the landfall region
    lfreg = sys.argv[5]
    lfreg = lfreg.split('/')                                       #split the string input by ','
    lfreg = [float(x) for x in lfreg]                        #convert the string numbers into floats

    sv_path=sys.argv[6]                                    #location to save data


    #############################################################################################################
    #############################################  RECEIVED INPUTS  ##############################################
    #############################################################################################################
    #check if all inputs are right 

    print(f'Variable name: MCS \nBottom Latitude: {blat} \nLeft Longitude: {llon} \nThreshold: {threshold} \nSeasons: {seas} \nSel Region Lon/Lat: {lfreg} ')



    #############################################################################################################
    ########################################## Tunable Parameters #################################################
    #############################################################################################################

    #path to the dataset for selection
    mpath = '/global/cfs/projectdirs/m1517/cascade/external_datasets/ARTMIP/tier1/ftp_mirror/' 

    #You can change this to suit what you want 

    lat_width = 60                                    #region latitude width
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
    #from dask.distributed import Client
    #c = Client(n_workers=25, threads_per_worker=1)

    mcs_path = '/global/cfs/cdirs/m4374/catalogues/raw_catalogue_files/observations/MCS/mcs_global/mcstracking_orig/'

    varbs = ['cloudtracknumber','pcptracknumber','tb','precipitation']
    dpts = []

    for yx in range(2000,2007):

        if yx == 2000:
            suff1 = '0601.0000'
            suff2 = '0101.0000'
        else: 
            suff1 = '0101.0000'
            suff2 = '0101.0000'

        mcs_ps = glob.glob(f'{mcs_path}{yx}{suff1}_{yx+1}{suff2}/*.nc', recursive=True)
        print('Done tracing paths')

        datasets = [dask.delayed(xr.open_mfdataset)(mcs_ps)]

        print(f'Loading datasets for {yx} ...')
        dataset = datasets[0].compute()

        data = dataset.chunk({'time':50}).sel(time=dataset.time.dt.month.isin(period))

        print(f'Done loading and chunking data for {yx}')

        #############################################################################################################
        ######################################## Set Threshold for selection  ##############################################
        #############################################################################################################

        #select the region for the landfall test 

        selection_box = data['cloudtracknumber'].sel(lon=slice(lfreg[0], lfreg[0]+selection_box_width),
                                lat=slice( lfreg[1],  lfreg[1]+selection_box_width))

        print(selection_box)
        #check if the threshold criteria is met within the selected region 

        sel_box=selection_box.chunk({'time':1})
        selected_times= sel_box.time.where(sel_box.mean(['lon','lat']) > threshold,drop=True)['time'].values


        print(f'Done with selection for {yx} \n{len(selected_times)}')
        selected_data = data.sel(time=selected_times)

        #############################################################################################################
        ########################################### Begin Saving Data  #################################################
        #############################################################################################################

        print(f'Started saving for {yx} ...')
        if os.path.isdir(sv_path)==False:
            os.mkdir(sv_path)

        final_path = f'{sv_path}/{reg}_MCS_phd/'

        if os.path.isdir(final_path)==False:
            os.mkdir(final_path)

        vbs = list(data.variables)[3:-3]
        saver.ch_paths(selected_data, f'{final_path}',f'MCS.gt.0.',vbs,sv_data_times=None)    
        print(f'Done with {yx} ...')
    print('Done with all years')