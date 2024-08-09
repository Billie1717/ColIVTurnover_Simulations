#!/bin/bash
#module load python3/recommended
MP=0.12
BP=0.04
k_angle=1.0
Blen=0.5
seed=1
t_stretch=1e4
Time=1509600
#Equilibrate/run_tbonds'${t_bonds}'_tsample'${t_sample}'_N_evSamp'${NevSamp}'_N_totSamp'${NtotSamp}'_MP'${MakeProb}'_BP'${BreakProb}'_MD'${MakeDist}'_BD'${BreakDist}'_seed'${seed}
#input='MP'${MP}'_BP'${BP}'/Equilibrate/run_T1.45000.restart'
#input='runsEquibSlab/run_tbonds1e6_tsample1e3_N_evSamp20_N_totSamp100_XStretch108_MP'${MP}'_BP'${BP}'_MD1.5_BD2.4_seed1/restart/run_T1.683400.restart' # {MP0.15 : 1917600}, {MP0.12 : 1642200}, {MP0.06 : 1815600}, {MP0.02 : 683400}
for seed in 1 3 #400 300 200 100 #500
do
for XStretch in 51 76
do
input='MP'${MP}'_BP'${BP}'/Equilibrate/runMMons_tbonds1e6_tsample1e3_N_evSamp20_N_totSamp100_MP'${MP}'_BP'${BP}'_MD1.3_BD2.4_seed'${seed}'/restart/run_T1.'${Time}'.restart'
#MP0.03_BP0.01/Equilibrate/runMMons_tbonds1e6_tsample1e3_N_evSamp20_N_totSamp100_MP0.03_BP0.01_MD1.3_BD2.4_seed1/restart/run_T1.1509600.restart
foldername1='MP'${MP}'_BP'${BP}'/Stretch/'
mkdir ${foldername1}
foldernameadd='MP'${MP}'_BP'${BP}'/Stretch/run_Simu_Dense_tstretch'${t_stretch}'_TimeFrom'${Time}'_xystretch'${XStretch}'_seed'${seed}
mkdir ${foldernameadd}
#input='Restart_PostMakeBonds/BD'${BreakDist}'_BP'${BreakProb}'_MD'${MakeDist}'_MP'${MakeProb}'/restart.postmakebonds'
cp ${input} ${foldernameadd}'/data'
cp -r 'templates' ${foldernameadd}
python3 build_StretchSimu.py ${foldernameadd} ${t_stretch} ${XStretch} ${k_angle} ${Blen} ${seed}
cd ${foldernameadd}
runscriptfile='runscript.sh'
sbatch $runscriptfile
cd .. 
cd ..
cd ..
done
done
