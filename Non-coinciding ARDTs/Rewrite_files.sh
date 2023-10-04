#!/bin/bash



#################################################################################################################
######################################## RECEIVE INPUT FROM THE USER ###########################################
#################################################################################################################
echo 'What region are you looking at? (north_america, west_europe, south_america, etc)'
read reg

echo 'What is the variable name? (z, pv, t, etc )'
read vb 

echo 'Which ARDT are you looking at? (cascade_bard_v1, mundhenk_v2, guan_waliser, reid250, etc)'
read fld


#files=$(ls ${reg}_${vb}/${fld}_${vb}) #LIST THE FILES IN THIS DIRECTORY INTO A VARIABLE
FILES=$(ls ${reg}_${vb}/${fld}_${vb}/*.nc)
for  f in $FILES
do
    #echo "Process $f file"
    ./Rewrite_files.py ${reg} ${vb} ${fld} $f& #APPLY INPUTS 
    #ADD OR REMOVE & TO WRITE DATA FROM FILES SIMULTANEOUSLY
    done