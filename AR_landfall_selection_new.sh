#!/bin/bash

#ardt=( 'cascade_bard_v1' 'mundhenk' 'guan_waliser' 'reid250' ) #name of ardt 

ardt=( "guan_waliser" )

blat=-60  #bottom latitude
llon=-115  #left longitude 

threshold=0.8 # sys.argv[4] #set the threshold for the selection of landfall (usually between 0 and 1)


#This takes a cordinate and expands to the right and top of it (i.e, if (60,15), then lon 60 to 61 and lat 15 to 16 will be selected as the box
#echo 'What is the landfall region? (input : lon/lat)'
#read lfreg
lfreg='-73.5/-44'
#echo "What months do you like to select? (Dec,Jan,Feb, or Feb,Mar,Apr or Apr,Jul,Oct...)"

#read seas
seas='Dec,Jan,Feb'
sv_path='larger_region_AR_capture'
for ff in "${!ardt[@]}"; do
    ./AR_landfall_selection_new.py ${ardt[ff]} $blat $llon $threshold $seas $lfreg $sv_path
    
    done