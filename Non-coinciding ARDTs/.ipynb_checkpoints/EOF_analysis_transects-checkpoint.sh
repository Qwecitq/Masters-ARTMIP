#!/bin/bash

#flds=('cascade_bard_v1' )
#flds=( 'mundhenk_v2' 'cascade_bard_v1' )
flds=( 'guan_waliser' 'reid250' )

#echo 'What cross section would you want? (options are NS, EW and across)'
#read whichtype
whichtype='NS' 
#echo 'What is the variable name or symbol? (q,pot,z,cape)'
#read vb
vb=( 'q' 'pv' 'pot' )

#echo 'How many PCAs? '
#read n_components
#n_components=4


for ix in "${!flds[@]}"; do
    f=${flds[ix]}
    echo $f
        for kx in "${!vb[@]}"; do
        echo ${vb[kx]}
        echo $whichtype
        ./EOF_analysis_transects.py $f ${vb[kx]} $whichtype 
        done
    done

