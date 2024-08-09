#!/bin/bash
#module load python3/recommended
MP=0.12
BP=0.04
xbox=51
seed=1
Time=1509600
t_stretch=1e4
OGMol=30 # LIST : 
for seed in 1 2 3 #400 300 200 100 #500
do
for XStretch in 51 76
do
foldername='MP'${MP}'_BP'${BP}'/Stretch/run_Simu_Dense_tstretch'${t_stretch}'_TimeFrom'${Time}'_xystretch'${XStretch}'_seed'${seed}
input='data.poststretchRel'
#foldername='Restart_PostMakeBonds/DENSE_BD1.3_BP0.0_MD1.3_MP1.0/'
#input='data.postreset'
#SYNTAX = $ python AddMonomersToDataFile.py {initial mon #} {initial xbox} {stretch amount}
cd ${foldername}
python /nfs/scistore15/saricgrp/bmeadowc/Scratch/Collagen/RiccyProject2/AddMonomersToDataFile2.py ${input} ${OGMol} ${xbox} ${XStretch}
cd .. 
cd ..
cd ..
done
done
