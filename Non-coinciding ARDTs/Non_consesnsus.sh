#!/bin/bash

#tr_ardts=( "mdk_tr" "teca_tr" "gw_tr" "rd_tr" ) 
#flds=( 'mundhenk_v2' 'cascade_bard_v1' 'guan_waliser' 'reid250' )

echo 'What is the variable name? (z, pv, t, etc )'
read varb 
echo 'What region are you looking at? (north_america, west_europe, south_america, etc)'
read reg

echo 'Which ARDT are you looking at? (cascade_bard_v1, mundhenk_v2, guan_waliser, reid250, etc)'
read fld

echo 'What is the type of data? (levels or surface data)'
read typ
#tr_ardts=( "mdk_tr" ) 
flds=($fld)

for ix in "${!flds[@]}"; do
    ./Non-Consensus.py ${reg} ${flds[ix]} ${varb} ${typ}&
    done