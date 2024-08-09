#!/bin/bash
#module load python3/recommended
MakeDist=1.3
BreakDist=1.3
k_angle=1.0
XStretch=72
BreakProb=0.0
MakeProb=1.0
Blen=0.5
t_mix=1e4
t_stretch=1e4
#t_relax=5e6
t_sample=1e3
t_bonds=1e6
NBM_Cyc=100
Nevery=500
NevSamp=20
NtotSamp=100
#input='Restart_PostMakeBonds/DENSE_BD'${BreakDist}'_BP'${BreakProb}'_MD'${MakeDist}'_MP'${MakeProb}'/restart.postmakebonds'
#input='Restart_PostMakeBonds/DENSE_BD'${BreakDist}'_BP'${BreakProb}'_MD'${MakeDist}'_MP'${MakeProb}'/run_T1.45000.restart'
input='Restart_PostMakeBonds/DENSE_BD'${BreakDist}'_BP'${BreakProb}'_MD'${MakeDist}'_MP'${MakeProb}'/data.added40mols'
for BreakDist in 2.4 #400 300 200 100 #500
do
for BreakProb in 0.02 #0.03 #0.1 0.03 0.01
do
for MakeProb in 0.06 #0.1 0.03 0.01
do
for seed in 1 #2 3 4 5
do
foldername1='MP'${MakeProb}'_BP'${BreakProb}
mkdir ${foldername1}
foldername2='MP'${MakeProb}'_BP'${BreakProb}'/Equilibrate/'
mkdir ${foldername2}
foldernameadd='MP'${MakeProb}'_BP'${BreakProb}'/Equilibrate/runMMons_tbonds'${t_bonds}'_tsample'${t_sample}'_N_evSamp'${NevSamp}'_N_totSamp'${NtotSamp}'_MP'${MakeProb}'_BP'${BreakProb}'_MD'${MakeDist}'_BD'${BreakDist}'_seed'${seed}
#foldernameadd='MakingDatafile'
mkdir ${foldernameadd}
#input='Restart_PostMakeBonds/BD'${BreakDist}'_BP'${BreakProb}'_MD'${MakeDist}'_MP'${MakeProb}'/restart.postmakebonds'
cp ${input} ${foldernameadd}'/data'
cp -r 'templates' ${foldernameadd}
python3 build_Equilibrate.py ${foldernameadd} ${t_mix} ${t_bonds} ${t_stretch} ${t_sample} ${NevSamp} ${NtotSamp} ${NBM_Cyc} ${XStretch} ${k_angle} ${Nevery} ${Blen} ${MakeDist} ${BreakDist} ${MakeProb} ${BreakProb} ${seed}
#python3 build_StretchSimu.py ${foldernameadd} ${t_stretch} ${XStretch} ${k_angle} ${Blen} ${seed}
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
