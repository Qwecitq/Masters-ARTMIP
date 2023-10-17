#!/bin/bash


FILES="Global_Potential_Temperature_sb/"
dir='zipped_Global_PT'

mkdir -p $dir

#Loop through the years
echo 'What is the variable name? (z, pv, t, etc )'
read varb 
for yy in {1980..2019}

do
    #Run the python script after
    
    
    for mm in {1..12}
    do 
        
                #written in the order of the sys.argv inputs 
            ./Climatologies.py $varb $yy $mm &
       done
    sleep 25m
done
            
            #dur=$(ls ${FILES}*${yy}*${mm}*)
            #num=1
#        else
#           ./Calc_potential_temperature.py $yy $mm 
        
 #           dur=$(ls ${FILES}*${yy}*${mm}*)
  #          num=1
   #             for ff in $dur #${FILES}*${yy}${mm}*${num}.nc
    #            do
     #               infile=$ff
      #              outfile=${FILES}e5.oper.an.pl_pt_${yy}${mm}


  #                  nccopy -c time/1024,level/37,longitude/1440,latitude/721 -d1 $infile ${outfile}_comp${num}.nc &
#

    #                num=$((num+1))

 

#for yy in {1980..1981}
#do
#for mm in {01..12}
#do 

    
#        infile=${FILES}*${yy}${mm}*.nc
 #       outfile=${dir}/e5.oper.an.pl_pt_${yy}${mm}
 #       echo "$yy Started processing ${FILES}*${yy}${mm}*.nc"
 #       nccopy -c time/1,level/1,longitude/561,latitude/201 -d1 $infile ${outfile}_comp.nc
 #       cdo -f nc mergetime $infile ${outfile}.nc
 #       
  #      #rm $outfile
        
#done
#done    
    
#module load python 
#conda activate -i climate_py38 
#-f nc4c -z zip_6

#./Calc_potential_temperature.py