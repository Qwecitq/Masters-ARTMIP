#!/bin/bash

#tr_ardts=( "mdk_tr" "teca_tr" "gw_tr" "rd_tr" ) 
#flds=( 'mundhenk_v2' 'cascade_bard_v1' 'guan_waliser' 'reid250' )

tr_ardts=( "teca_tr" ) 
flds=('cascade_bard_v1' )

for ix in "${!flds[@]}"; do
    #printf "${flds[ix]} is for ${tr_ardts[ix]}"
    echo 'What cross section would you want? (options are NS, EW and across)'
    read whichtype
    ./Transect_plots_with_winds.py ${flds[ix]} ${tr_ardts[ix]} $whichtype&
    done