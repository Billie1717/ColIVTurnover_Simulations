# The inputs required and typical parameters used for each simulation protocol:

Python files run by eg:

python3 build_Equilibrate.py ${foldernameadd} ${t_mix} ${t_bonds} ${t_stretch} ${t_sample} ${NevSamp} ${NtotSamp} ${NBM_Cyc} ${XStretch} ${k_angle} ${Nevery} ${Blen} ${MakeDist} ${BreakDist} ${MakeProb} ${BreakProb} ${seed}

## /slab_stretchXY/


Running this python script to create lammps input files requires the following inputs:


 - build_Equilibrate.py:
 foldernameadd----- folder_name------ folder name for simulation 
 t_mix------------- 1e4-------------- time to mix protomers in simulation steps  
 t_bonds----------- 1e6-------------- time to make bonds in simulation steps 
 t_stretch--------- 1e4-------------- time to stretch in simulation steps 
 t_sample---------- 1e3-------------- The length of time to sample much more frequently the bonds (In this time, bonds are dumped every NevSample timesteps (see below)
 NevSamp----------- 20--------------- The rate of sampling (every NevSamp steps) for the frequent bond sampling stages
 NtotSamp---------- 100-------------- The total number of samples taken in the simulation
 NBM_Cyc----------- 100-------------- number of make and break bond cycles ##(actual number > as in relax = 5x NBM_cyc and make bonds 1x NBM_cyc)???##
 XStretch---------- 51--------------- amount of stretch in x and y-dirction in simulation length units, vary between 51 and 76
 k_angle----------- 1.0-------------- Angular bond constant in kBT
 Nevery------------ 500-------------- How often we check for breaking a bond, must be quite high to allow bonds to react (type change) in between. Vary between 200 and 500.
 Blen-------------- 0.5-------------- Length of harmonic bond between protomers (in simulation units)
 MakeDist---------- 1.3-------------- The factor of the bond length at which a bond is considered to be made (it will be considered at MakeDist * Blen)
 BreakDist--------- 2.4-------------- The factor of the bond length at which a bond is considered to be broken (it will be considered at BreakDist * Blen) varies between {2.4,4}
 MakeProb---------- 0.06------------- The probability of making a bond between two protomers which are < MakeDist * Blen away. Kept at 3x BreakProb for constant ratio
 BreakProb--------- 0.02------------- The probability of breaking a bond when it is > BreakDist * Blen away. Varies between {0.007,0.04}
 seed-------------- 1---------------- Random seed for - integer, for creating random seeds for initial velocity and lengevin thermostat

 - build_StretchSimu.py
 foldernameadd----- folder_name------ folder name for simulation 
 t_stretch--------- 1e4-------------- time to stretch in simulation steps 
 XStretch---------- 51--------------- amount of stretch in x and y-dirction in simulation length units, vary between 51 and 76
 k_angle----------- 1.0-------------- Angular bond constant in kBT
 Blen-------------- 0.5-------------- Length of harmonic bond between protomers (in simulation units)
 seed-------------- 1---------------- Random seed for - integer, for creating random seeds for initial velocity and lengevin thermostat
 
 - build_Relax.py
 Same inputs as build_Equilibrate.py
 
 - AddMonomersToDataFile.py
 input-------------'data.poststretchRel' - The input data which will be reproduced with some added monomers
 OGMol------------- 30--------------- Number of free monomers in the input simulation
 xbox-------------- 51--------------- The size of the simulation box in x (& y)
 XStretch---------- 51--------------- amount of stretch that happened in that simulation in x- and y- directions
