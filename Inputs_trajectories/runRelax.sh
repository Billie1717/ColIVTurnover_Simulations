#!/bin/bash
#module load python3/recommended
MakeDist=1.5
BreakDist=2.4
k_angle=1.0
XStretch=72
BreakProb=0.1
MakeProb=0.01
Blen=0.5
t_mix=1e4
t_stretch=1e4
#t_relax=5e6
t_sample=1e3
t_bonds=1e5
NBM_Cyc=100
Nevery=500
NevSamp=20
NtotSamp=100
XStretch=108
Time=1509600
for XStretch in 51 76 #400 300 200 100 #500
do
for BreakProb in 0.04 #0.03 #0.1 0.03 0.01
do
for MakeProb in 0.12 #0.1 0.03 0.01
do
for seed in 1 2 3 #4 5
do
# 80, 190, 330
# 100, 250, 440
# 140 320 540
# 60,90,120
input='MP'${MakeProb}'_BP'${BreakProb}'/Stretch/run_Simu_Dense_tstretch'${t_stretch}'_TimeFrom'${Time}'_xystretch'${XStretch}'_seed'${seed}'/data.poststretchReladdedMons'
foldername1='MP'${MakeProb}'_BP'${BreakProb}'/Relax'
foldernameadd='MP'${MakeProb}'_BP'${BreakProb}'/Relax/run_Simu_Dense_TimeFrom'${Time}'_tbonds'${t_bonds}'_tsample'${t_sample}'_Xstretch_'${XStretch}'_N_evSamp'${NevSamp}'_N_totSamp'${NtotSamp}'_MP'${MakeProb}'_BP'${BreakProb}'_MD'${MakeDist}'_BD'${BreakDist}'_seed'${seed}
mkdir ${foldername1}
mkdir ${foldernameadd}
#input='Restart_PostMakeBonds/BD'${BreakDist}'_BP'${BreakProb}'_MD'${MakeDist}'_MP'${MakeProb}'/restart.postmakebonds'
cp ${input} ${foldernameadd}'/data'
cp -r 'templates' ${foldernameadd}
seed=2
python3 build_Relax.py ${foldernameadd} ${t_mix} ${t_bonds} ${t_stretch} ${t_sample} ${NevSamp} ${NtotSamp} ${NBM_Cyc} ${XStretch} ${k_angle} ${Nevery} ${Blen} ${MakeDist} ${BreakDist} ${MakeProb} ${BreakProb} ${seed}
cd ${foldernameadd}
runscriptfile='runscript.sh'
sbatch $runscriptfile
cd .. 
cd ..
cd ..
done
done
done
done
