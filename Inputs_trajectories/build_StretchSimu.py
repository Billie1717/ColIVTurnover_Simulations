from __future__ import division, print_function
import numpy as np
#import matplotlib.pyplot as plt
import math
import time
import random
import sys, getopt
import os


def main():

    ################ Parameters ################
    workingdir = str(sys.argv[1])
    tstretch = str(sys.argv[2])
    XStretch = float(sys.argv[3])
    K_angle = str(sys.argv[4])
    Blen=float(sys.argv[5])
    seed  = str(sys.argv[6])
    Nstretch = 20
    trelax = float(tstretch)*10
    
    #cutoff = np.sqrt((ChemBond-np.log(0.01))/6)+Blen #Cutoff the making probabbilty at Prob = 0.01, with bond constant 6 
    #print("cutoff"+str(cutoff))
    ################ write in.local ################
    print('Writing collagen_dumbbell_wall_bond_create.in\n')
    outfile2 = workingdir + '/collagen_dumbbell_wall_bond_create.in'
    fo=open(outfile2,'w')
    fo.write(
        '''#####################################################
#PARAMETERS

#NUMBER OF STEPS OF RUNS (stages)
''')
    fo.write('variable tstretch			equal 	'+str(tstretch)+'	#NVT steps for stretching\n')
    fo.write('variable trelax			equal 	'+str(trelax)+'	#NVT steps for stretching\n')
    fo.write('variable tstretchcycl			equal 	'+str(float(tstretch)/Nstretch)+'	#NVT steps for stretching\n')
    fo.write(
        '''

#THERMODYNAMIC PARAMETERS
variable mytemp         equal 	1.0	    	#Temperature

#INFORMATION PRINTING PERIODS
variable ttotal     equal   $(4*v_tstretch)		#total timesteps #v_trelax
variable thermodump     equal   $(floor(v_ttotal/1e3))		#Thermo info dump period
variable trajdump       equal   $(floor(v_ttotal/20))    	#Config info dump period''')

    fo.write(
        '''
#OTHER SIMULATION PARAMETERS
variable tstep        	equal   0.01   	#Integration time step for first part of simulation
variable damp   		equal   0.1	   	#Langevin thermostat damping coefficient [F_friction = -(m/damp)*v]
variable neigh_cutoff	equal	0.4  		#Neighbor list cutoff
variable lj_cutoff      equal   1.122462	#For WCA: 1.122462048 
#variable comm_cutoff	equal   5. 			#Ghost cutoff\n\n
''')
    fo.write('variable k_nc1			equal	6.0 		#Elastic constant A type bond\n')
    fo.write('variable r0_nc1			equal	'+str(Blen)+' 		#Rest length A type bond\n\n')

    fo.write('variable k_7s			equal	6.0 		#Elastic constant B type bond\n')
    fo.write('variable r0_7s			equal	'+str(Blen)+' 		#Rest length B type bond\n\n')
    #fo.write('variable ChemBond			equal	'+str(ChemBond)+' 		#Rest length B type bond\n\n')


    fo.write(
        '''variable zero_cutoff    equal   6.0				#Cutoff for zero potential; must be larger that max bonding distance 

#SEEDS
variable lseed			file 	lseed.dat	#Seed for fix langevin
variable vseed			file 	vseed.dat	#Seed for creation of velocities

#STRINGS
variable fname 			string 	data		#name of the data file to be read to start the simulation
variable thermofile     string  thermo.dat	#name of the file with thermodynamic infos
variable simname 		string 	run_T${mytemp}	#name of this simulation

#####################################################
#INITALIZATION
units           lj          	#use lj units
boundary        p p f       	#periodic boundary conditions
atom_style     	molecular

#BOND STYLE: 
bond_style      harmonic

#ANGLE STYLE: 
angle_style      harmonic

#TEMPORARY PAIR STYLE TO START SIMULATION
pair_style		lj/cut 1.0

#Starting from an initial configuration WITHOUT BONDS (format lammpsdata):
#read_data       ${fname} #extra/bond/types 2 extra/bond/per/atom 10 extra/special/per/atom 100 extra/angle/types 2 extra/angle/per/atom 20
read_restart ${fname}
#reset_timestep  0

#GROUPS

group  			nc1 type 1
group  			7s type 2 3 4 5 6 7
group  			7sBreaking  type 2 3 4 5
group  			7s_trimer  type 4

#BONDS COEFFICIENTS

bond_coeff     	1 	100.0	3.0
bond_coeff     	2 ${k_nc1} ${r0_nc1}
bond_coeff		3 ${k_7s} ${r0_7s}

#ANGLE COEFFICIENTS\n
''')

    fo.write('angle_coeff     1 	'+str(K_angle)+'	180.0\n')
    fo.write('angle_coeff     2 	'+str(K_angle)+'	120.0\n')
    fo.write(
        '''

#NEIGHBOR LISTS
neighbor		${neigh_cutoff} bin 		#value = skin = extra distance beyond cutoff
neigh_modify	every 10 delay 0 check yes

#PAIR STYLE (WCA)
pair_style     	hybrid/overlay lj/cut 1.0 zero 1.0
pair_coeff      * * lj/cut 0.0 0.0 ${lj_cutoff}
#pair_modify     pair lj/cut shift yes
pair_coeff      * * zero ${zero_cutoff}

#comm_modify     cutoff ${comm_cutoff} #Extend ghost cutoff

#SET TIME STEP
timestep	${tstep}

#COMPUTE BOND PROPERTIES
compute 		1 all property/local batom1 batom2 btype
compute 		4 all property/local atype aatom1 aatom2 aatom3
compute 		2 all pe bond
compute 		3 all pe angle
compute perAtomStress all stress/atom NULL
compute totalStressX all reduce sum c_perAtomStress[1]
compute totalStressY all reduce sum c_perAtomStress[2]
compute totalStressZ all reduce sum c_perAtomStress[3]
compute forceBond all bond/local force
compute forceBondAll all reduce ave c_forceBond

####-PRINTING STUFF-##############################

#DEFINITION OF VARIABLES FOR USE IN FIX PRINT
variable step	equal step
variable temp	equal temp
variable etot	equal etotal
variable ke		equal ke
variable peV	equal pe
variable press	equal press
variable peBond	equal c_2
variable peAngle	equal c_3
variable StressX equal c_totalStressX
variable StressY equal c_totalStressY
variable StressZ equal c_totalStressZ
variable Fbond equal c_forceBondAll

#THERMODYNAMIC INFO
thermo_style	custom step etotal ke pe temp press c_totalStressX c_totalStressY c_totalStressZ
thermo		${thermodump}
thermo_modify	norm no

##################################################
#NVT (NVE+LANGEVIN) RUN

fix 		fix_wall_lo 	all wall/lj126 zlo EDGE 1.0 1.0 ${lj_cutoff} #Lower wall
fix 		fix_wall_hi 	all wall/lj126 zhi EDGE 1.0 1.0 ${lj_cutoff} #Upper wall
fix 		fix_lengevin	all langevin ${mytemp} ${mytemp} ${damp} ${lseed}
fix	fNVE	all	nve 

#CREATE VELOCITIES USING MAXWELL-BOLTZMANN DISTRIBUTION
velocity 	all create ${mytemp} ${vseed} dist gaussian mom yes rot yes 
run 0
velocity	all scale ${mytemp}

fix		fix_print all print ${thermodump} "${step} ${etot} ${ke} ${peBond} ${peAngle} ${temp} ${press} ${StressX} ${StressY} ${StressZ} ${Fbond}" file ${thermofile} screen no title "step etot ke peBond peAngle temp press stressX stressY stressZ AvBondForce"

#Save configurations in linear time
dump		2 all custom ${trajdump} dumplin/dump.${simname}.*.lammpstrj type mol xu yu zu
dump_modify	2 sort id
dump_modify	2 pad 10
dump_modify	2 format line "%d %d %.8f %.8f %.8f"

#Save bond information
dump 		3 all local ${trajdump} dumplin_bonds/bonds.${simname}.* index c_1[1] c_1[2] c_1[3] #bond index, id of atom 1 and 2 in the bond, bond type
dump_modify	3 format line "%d %0.0f %0.0f %0.0f"
dump_modify	3 pad 10

neigh_modify exclude molecule/inter 7s_trimer
reset_atoms mol 7s

special_bonds lj 0.0 1.0 1.0
run                    ${tstretch}
''')
    for i in range(Nstretch):
        #fo.write('run                    ${tstretch}\n')
        fo.write('fix                    fDeformX all deform 1 x delta '+str(-XStretch*0.5/Nstretch)+' '+str(XStretch*0.5/Nstretch)+ ' remap x\n')
        fo.write('run                    ${tstretchcycl}\n')
        fo.write('unfix fDeformX\n')
        fo.write('fix                    fDeformY all deform 1 y delta '+str(-XStretch*0.5/Nstretch)+' '+str(XStretch*0.5/Nstretch)+ ' remap x\n')
        fo.write('run                    ${tstretchcycl}\n')
        fo.write('unfix fDeformY\n')
    
    fo.write('''
write_data	data.poststretch
write_restart	restart.poststretch
run                    ${trelax}
write_data	data.poststretchRel
write_restart	restart.poststretchRel

\n''')






#\n''')
    
    fo.close()

    ################ write runscript ################
    runscriptfilename =  workingdir +"/runscript.sh"
    print('Writing runscript.sh\n')
    f = open(runscriptfilename, "w")
    f.write(
        '''#!/bin/bash
#
#SBATCH -J MakingBreakingOutputs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=400M 
#SBATCH --time=72:00:00 
#SBATCH --mail-user=ucapbbm@ucl.ac.uk
#SBATCH --mail-type=END
#SBATCH --export=NONE

module load python
module load gcc
module load openmpi
module load lammps/20220623b\n
mkdir dumplin dumplin_bonds restart dumplin_bondsSample
printf $RANDOM > lseed.dat
printf $RANDOM > vseed.dat\n
mpirun -np 1 /nfs/scistore15/saricgrp/bmeadowc/Scratch/lammps-15Jun2023/src/lmp_mpi -in collagen_dumbbell_wall_bond_create.in
\n''')
    f.close()

##SBATCH --constraint="epsilon|delta|beta|leonid|serbyn|gamma"     
if __name__ == "__main__":
    main()
