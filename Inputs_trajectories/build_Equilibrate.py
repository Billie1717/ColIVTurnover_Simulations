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
    NumFixCycl = int(sys.argv[6])
    NumStretch = 1 #int(sys.argv[6])
    workingdir = str(sys.argv[1])
    tmix = str(sys.argv[2])
    tbonds = str(sys.argv[3])
    tstretch = str(sys.argv[4])
    tSample = str(sys.argv[5])
    NevSample = str(sys.argv[6])
    NtotSample = str(sys.argv[7])
   # trelax = str(sys.argv[5])
    NumFixCycl = int(sys.argv[8])
    XStretch = float(sys.argv[9])
    K_angle = str(sys.argv[10])
    Nevery = str(sys.argv[11])
    Nevery_Cr2 = 1 
    MakeDist2 = 1.5
    Blen=float(sys.argv[12])
    MakeDist = str(sys.argv[13])
    BreakDist = str(sys.argv[14])
    MakeProb = str(sys.argv[15])
    BreakProb = str(sys.argv[16])
    
    seed  = str(sys.argv[17])
    
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
    fo.write('variable tmix			equal 	'+str(tmix)+'	#NVT steps for molecule mixing\n')
    fo.write('variable tbonds			equal 	'+str(tbonds)+'	#NVT steps for bond equilibration\n')
    fo.write('variable tstretch			equal 	'+str(tstretch)+'	#NVT steps for stretching\n')
    fo.write('variable tsampleTop			equal 	'+str(tSample)+'	#NVT steps for sampling bond topology\n')
    fo.write('variable NevSample			equal 	'+str(NevSample)+'	#NVT every this many steps during sampling\n')
    fo.write('variable NevSample			equal 	'+str(NevSample)+'	#NVT every this many steps during sampling\n')
    fo.write('variable NtotSample			equal 	'+str(NtotSample)+'	#This number of total samples taken in the whole run\n')
    #fo.write('variable total			equal 	'+str(tbonds+tmix+tstretch)+'	#NVT steps for post stretch relaxation\n')
    fo.write('variable Tbetween			equal 	$(floor(((v_tbonds)/(v_NtotSample)-v_tsampleTop)/2))	#NVT steps between each sample\n')
    fo.write('variable trelax			equal 	'+str(0)+'	#NVT steps for post stretch relaxation\n')
    
    fo.write(
        '''

#THERMODYNAMIC PARAMETERS
variable mytemp         equal 	1.0	    	#Temperature

#INFORMATION PRINTING PERIODS
variable ttotal     equal   $(v_tmix+v_tbonds+v_tstretch)		#total timesteps #v_trelax
variable thermodump     equal   $(floor(v_ttotal/5e4))		#Thermo info dump period
variable trajdump       equal   $(floor(v_ttotal/500))    	#Config info dump period
variable trestart       equal   $(floor(v_ttotal/100))    	#Restart configuration saving period\n''')

    #fo.write('variable fixstretch   equal   $(floor(v_nsteps/20))    	#Restart configuration saving period\n')
    fo.write('variable tcycleMB   equal   $(floor((v_tbonds)/'+str(2*NumFixCycl)+'))    	#Restart configuration saving period\n')
    fo.write('variable tcycleRel   equal   $(floor((v_trelax)/'+str(2*NumFixCycl)+'))    	#Restart configuration saving period\n')
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

#BOND CREATION/BREAKING PARAMETERS 
#NC1 bond (1 per molecule)\n''')
    
    fo.write('variable break_prob_nc1	equal   '+str(BreakProb)+'				#Probability to create a bond\n')
    fo.write('variable make_prob_nc1	equal   '+str(MakeProb)+'				#Probability to create a bond\n')
    fo.write('variable make_dist_nc1  	equal 	$(v_r0_nc1*'+str(MakeDist)+') 	#Bonding distance\n')
    fo.write('variable break_dist_nc1  	equal 	$(v_r0_nc1*'+str(BreakDist)+') 	#Breaking distance\n')
    fo.write('variable max_bonds_nc1	equal   1 				#Max. number of bonds of type A that each particle can form\n')
    fo.write('variable nevery_Cr		equal	'+str(Nevery)+' 				#Attempt bond creation every "nevery" time steps\n')
    fo.write('variable nevery_Cr2		equal	'+str(Nevery_Cr2)+' 				#Attempt bond creation every "nevery" time steps\n')
    fo.write('variable nevery_Br		equal	'+str(Nevery)+' 				#Attempt bond breakage every "nevery" time steps\n')

    fo.write('#7S bonds (3 per molecule)\n')
    fo.write('variable break_prob_7s	equal   '+str(BreakProb)+'				#Probability to create a bond\n')
    fo.write('variable make_prob_7s	equal   '+str(MakeProb)+'				#Probability to create a bond\n')
    fo.write('variable make_dist_7s  	equal 	$(v_r0_nc1*'+str(MakeDist)+') 	#Bonding distance\n')
    fo.write('variable break_dist_7s  	equal 	$(v_r0_nc1*'+str(BreakDist)+') 	#Breaking distance\n')
    fo.write('variable max_bonds_7s	equal   3				#Max. number of bonds of type A that each particle can form\n')
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
read_data       ${fname} extra/bond/per/atom 10 extra/special/per/atom 100 extra/angle/per/atom 20
#read_restart ${fname}
reset_timestep  0

#GROUPS
group  			nc1 type 1
group  			7s  type 2 3 4 5 6 7
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

#MOLECULE TEMPLATES
molecule            mPredimer templates/dimer/pre_dimer.txt
molecule            mPostdimer templates/dimer/post_dimer.txt
molecule            mPredimer_1 templates/dimer/pre_dimer_1.txt
molecule            mPostdimer_1 templates/dimer/post_dimer_1.txt
molecule            mPredimer_2 templates/dimer/pre_dimer_2.txt
molecule            mPostdimer_2 templates/dimer/post_dimer_2.txt
molecule            mPredimer_3 templates/dimer/pre_dimer_3.txt
molecule            mPostdimer_3 templates/dimer/post_dimer_3.txt

molecule            mPretrimer templates/trimer/pre_trimer.txt
molecule            mPosttrimer templates/trimer/post_trimer.txt
molecule            mPretrimer_1 templates/trimer/pre_trimer_1.txt
molecule            mPosttrimer_1 templates/trimer/post_trimer_1.txt
molecule            mPretrimer_2 templates/trimer/pre_trimer_2.txt
molecule            mPosttrimer_2 templates/trimer/post_trimer_2.txt
molecule            mPretrimer_3 templates/trimer/pre_trimer_3.txt
molecule            mPosttrimer_3 templates/trimer/post_trimer_3.txt
molecule            mPretrimer_4 templates/trimer/pre_trimer_4.txt
molecule            mPosttrimer_4 templates/trimer/post_trimer_4.txt
molecule            mPretrimer_5 templates/trimer/pre_trimer_5.txt
molecule            mPosttrimer_5 templates/trimer/post_trimer_5.txt
molecule            mPretrimer_6 templates/trimer/pre_trimer_6.txt
molecule            mPosttrimer_6 templates/trimer/post_trimer_6.txt
molecule            mPretrimer_7 templates/trimer/pre_trimer_7.txt
molecule            mPosttrimer_7 templates/trimer/post_trimer_7.txt

molecule            mPretrimer_1_1 templates/trimer/pre_trimer_1_1.txt
molecule            mPosttrimer_1_1 templates/trimer/post_trimer_1_1.txt
molecule            mPretrimer_1_2 templates/trimer/pre_trimer_1_2.txt
molecule            mPosttrimer_1_2 templates/trimer/post_trimer_1_2.txt
molecule            mPretrimer_1_3 templates/trimer/pre_trimer_1_3.txt
molecule            mPosttrimer_1_3 templates/trimer/post_trimer_1_3.txt
molecule            mPretrimer_1_4 templates/trimer/pre_trimer_1_4.txt
molecule            mPosttrimer_1_4 templates/trimer/post_trimer_1_4.txt
molecule            mPretrimer_1_5 templates/trimer/pre_trimer_1_5.txt
molecule            mPosttrimer_1_5 templates/trimer/post_trimer_1_5.txt
molecule            mPretrimer_1_6 templates/trimer/pre_trimer_1_6.txt
molecule            mPosttrimer_1_6 templates/trimer/post_trimer_1_6.txt

molecule            mPretrimer_F templates/trimerFive/pre_trimerF.txt
molecule            mPosttrimer_F templates/trimerFive/post_trimerF.txt
molecule            mPretrimer_F2 templates/trimerFive/pre_trimerF_2.txt
molecule            mPosttrimer_F2 templates/trimerFive/post_trimerF_2.txt
molecule            mPretrimer_F3 templates/trimerFive/pre_trimerF_3.txt
molecule            mPosttrimer_F3 templates/trimerFive/post_trimerF_3.txt
molecule            mPretrimer_F4 templates/trimerFive/pre_trimerF_4.txt
molecule            mPosttrimer_F4 templates/trimerFive/post_trimerF_4.txt
molecule            mPretrimer_F_1 templates/trimerFive/pre_trimerF_1.txt
molecule            mPosttrimer_F_1 templates/trimerFive/post_trimerF_1.txt
molecule            mPretrimer_F_12 templates/trimerFive/pre_trimerF_12.txt
molecule            mPosttrimer_F_12 templates/trimerFive/post_trimerF_12.txt

molecule            mPretetramer_zero templates/tetramer/pre_tetramer_zero.txt
molecule            mPosttetramer_zero templates/tetramer/post_tetramer_zero.txt
molecule            mPretetramer_one_1 templates/tetramer/pre_tetramer_one_1.txt
molecule            mPosttetramer_one_1 templates/tetramer/post_tetramer_one_1.txt
molecule            mPretetramer_one_2 templates/tetramer/pre_tetramer_one_2.txt
molecule            mPosttetramer_one_2 templates/tetramer/post_tetramer_one_2.txt
molecule            mPretetramer_one_3 templates/tetramer/pre_tetramer_one_3.txt
molecule            mPosttetramer_one_3 templates/tetramer/post_tetramer_one_3.txt
molecule            mPretetramer_two_1 templates/tetramer/pre_tetramer_two_1.txt
molecule            mPosttetramer_two_1 templates/tetramer/post_tetramer_two_1.txt
molecule            mPretetramer_two_2 templates/tetramer/pre_tetramer_two_2.txt
molecule            mPosttetramer_two_2 templates/tetramer/post_tetramer_two_2.txt
molecule            mPretetramer_two_3 templates/tetramer/pre_tetramer_two_3.txt
molecule            mPosttetramer_two_3 templates/tetramer/post_tetramer_two_3.txt
molecule            mPretetramer_two_4 templates/tetramer/pre_tetramer_two_4.txt
molecule            mPosttetramer_two_4 templates/tetramer/post_tetramer_two_4.txt
molecule            mPretetramer_three_1 templates/tetramer/pre_tetramer_three_1.txt
molecule            mPosttetramer_three_1 templates/tetramer/post_tetramer_three_1.txt
molecule            mPretetramer_three_2 templates/tetramer/pre_tetramer_three_2.txt
molecule            mPosttetramer_three_2 templates/tetramer/post_tetramer_three_2.txt
molecule            mPretetramer_three_3 templates/tetramer/pre_tetramer_three_3.txt
molecule            mPosttetramer_three_3 templates/tetramer/post_tetramer_three_3.txt
molecule            mPretetramer_four templates/tetramer/pre_tetramer_four.txt
molecule            mPosttetramer_four templates/tetramer/post_tetramer_four.txt

molecule            mPretetramer_1bond_two_1 templates/tetramer/pre_tetramer_1bond_two_1.txt
molecule            mPosttetramer_1bond_two_1 templates/tetramer/post_tetramer_1bond_two_1.txt
molecule            mPretetramer_1bond_two_2 templates/tetramer/pre_tetramer_1bond_two_2.txt
molecule            mPosttetramer_1bond_two_2 templates/tetramer/post_tetramer_1bond_two_2.txt
molecule            mPretetramer_1bond_two_3 templates/tetramer/pre_tetramer_1bond_two_3.txt
molecule            mPosttetramer_1bond_two_3 templates/tetramer/post_tetramer_1bond_two_3.txt

molecule            mPretetramer_1bond_three_1 templates/tetramer/pre_tetramer_1bond_three_1.txt
molecule            mPosttetramer_1bond_three_1 templates/tetramer/post_tetramer_1bond_three_1.txt
molecule            mPretetramer_1bond_three_2 templates/tetramer/pre_tetramer_1bond_three_2.txt
molecule            mPosttetramer_1bond_three_2 templates/tetramer/post_tetramer_1bond_three_2.txt
molecule            mPretetramer_1bond_three_3 templates/tetramer/pre_tetramer_1bond_three_3.txt
molecule            mPosttetramer_1bond_three_3 templates/tetramer/post_tetramer_1bond_three_3.txt
molecule            mPretetramer_1bond_three_4 templates/tetramer/pre_tetramer_1bond_three_4.txt
molecule            mPosttetramer_1bond_three_4 templates/tetramer/post_tetramer_1bond_three_4.txt
molecule            mPretetramer_1bond_three_5 templates/tetramer/pre_tetramer_1bond_three_5.txt
molecule            mPosttetramer_1bond_three_5 templates/tetramer/post_tetramer_1bond_three_5.txt

molecule            mPretetramer_1bond_four_1 templates/tetramer/pre_tetramer_1bond_four_1.txt
molecule            mPosttetramer_1bond_four_1 templates/tetramer/post_tetramer_1bond_four_1.txt
molecule            mPretetramer_1bond_four_2 templates/tetramer/pre_tetramer_1bond_four_2.txt
molecule            mPosttetramer_1bond_four_2 templates/tetramer/post_tetramer_1bond_four_2.txt
molecule            mPretetramer_1bond_four_3 templates/tetramer/pre_tetramer_1bond_four_3.txt
molecule            mPosttetramer_1bond_four_3 templates/tetramer/post_tetramer_1bond_four_3.txt

molecule            mPretetramer_2bond_1 templates/tetramer/pre_tetramer_2bond_1.txt
molecule            mPosttetramer_2bond_1 templates/tetramer/post_tetramer_2bond_1.txt
molecule            mPretetramer_2bond_2 templates/tetramer/pre_tetramer_2bond_2.txt
molecule            mPosttetramer_2bond_2 templates/tetramer/post_tetramer_2bond_2.txt

#BREAKING MOLECULE TEMPLATES

molecule            mPredimer1Br templates/7sbreaking/monomers/pre_dimer1.txt
molecule            mPostdimer1Br templates/7sbreaking/monomers/post_dimer1.txt
molecule            mPredimer2Br templates/7sbreaking/monomers/pre_dimer2.txt
molecule            mPostdimer2Br templates/7sbreaking/monomers/post_dimer2.txt
molecule            mPredimer1Br_6 templates/7sbreaking/monomers/pre_dimer1_6.txt
molecule            mPredimer2Br_6 templates/7sbreaking/monomers/pre_dimer2_6.txt
molecule            mPredimer1Br_4 templates/7sbreaking/monomers/pre_dimer1_4.txt
molecule            mPredimer2Br_4 templates/7sbreaking/monomers/pre_dimer2_4.txt
molecule            mPredimer1Br_5 templates/7sbreaking/monomers/pre_dimer1_5.txt #Shouldn't really be here but just in case?
molecule            mPredimer2Br_5 templates/7sbreaking/monomers/pre_dimer2_5.txt
molecule            mPredimer1Br_7 templates/7sbreaking/monomers/pre_dimer1_7.txt #Shouldn't really be here but just in case?
molecule            mPredimer2Br_7 templates/7sbreaking/monomers/pre_dimer2_7.txt

molecule            mPretrimer1Br templates/7sbreaking/pre_trimerF_1.txt
molecule            mPosttrimer1Br templates/7sbreaking/post_trimerF_1.txt
molecule            mPretrimer2Br templates/7sbreaking/pre_trimerF_2.txt
molecule            mPosttrimer2Br templates/7sbreaking/post_trimerF_2.txt
molecule            mPretrimer3Br templates/7sbreaking/pre_trimerF_3.txt
molecule            mPosttrimer3Br templates/7sbreaking/post_trimerF_3.txt
molecule            mPretrimer4Br templates/7sbreaking/pre_trimerF_4.txt
molecule            mPosttrimer4Br templates/7sbreaking/post_trimerF_4.txt
molecule            mPretrimer5Br templates/7sbreaking/pre_trimerF_5.txt
molecule            mPosttrimer5Br templates/7sbreaking/post_trimerF_5.txt
molecule            mPretrimer6Br templates/7sbreaking/pre_trimerF_6.txt
molecule            mPosttrimer6Br templates/7sbreaking/post_trimerF_6.txt

molecule            mPredimerb templates/7sbreaking/dimer/pre_dimer.txt
molecule            mPredimerb_1 templates/7sbreaking/dimer/pre_dimer_1.txt
molecule            mPredimerb_2 templates/7sbreaking/dimer/pre_dimer_2.txt
molecule            mPredimerb_3 templates/7sbreaking/dimer/pre_dimer_3.txt

molecule            mPretrimerb_2 templates/7sbreaking/trimer/pre_trimer_2.txt
molecule            mPretrimerb_3 templates/7sbreaking/trimer/pre_trimer_3.txt
molecule            mPretrimerb_5 templates/7sbreaking/trimer/pre_trimer_5.txt
molecule            mPretrimerb_4 templates/7sbreaking/trimer/pre_trimer_4.txt
molecule            mPretrimerb_7 templates/7sbreaking/trimer/pre_trimer_7.txt
molecule            mPretrimerb_1_1 templates/7sbreaking/trimer/pre_trimer_1_1.txt
molecule            mPretrimerb_1_3 templates/7sbreaking/trimer/pre_trimer_1_3.txt
molecule            mPretrimerb_1_4 templates/7sbreaking/trimer/pre_trimer_1_4.txt
molecule            mPretrimerb_1_6 templates/7sbreaking/trimer/pre_trimer_1_6.txt
molecule            mPretrimerb templates/7sbreaking/trimer/pre_trimer.txt

molecule            mPretetramerb templates/7sbreaking/tetramer/pre_tetramer_zero.txt
molecule            mPosttetramerb templates/7sbreaking/tetramer/post_tetramer_zero.txt
molecule            mPretetramerb1 templates/7sbreaking/tetramer/pre_tetramer_one.txt
molecule            mPosttetramerb1 templates/7sbreaking/tetramer/post_tetramer_one.txt
molecule            mPretetramerb2 templates/7sbreaking/tetramer/pre_tetramer_two.txt
molecule            mPretetramerb22 templates/7sbreaking/tetramer/pre_tetramer_two2.txt
molecule            mPosttetramerb2 templates/7sbreaking/tetramer/post_tetramer_two.txt
molecule            mPosttetramerb22 templates/7sbreaking/tetramer/post_tetramer_two2.txt

molecule            mPretetramer1Bb1 templates/7sbreaking/tetramer/pre_tetramer_1bond_two.txt
molecule            mPosttetramer1Bb1 templates/7sbreaking/tetramer/post_tetramer_1bond_two.txt
molecule            mPretetramer1Bb12 templates/7sbreaking/tetramer/pre_tetramer_1bond_two2.txt
molecule            mPosttetramer1Bb12 templates/7sbreaking/tetramer/post_tetramer_1bond_two2.txt
molecule            mPretetramerb3 templates/7sbreaking/tetramer/pre_tetramer_three.txt
molecule            mPosttetramerb3 templates/7sbreaking/tetramer/post_tetramer_three.txt
molecule            mPretetramer1Bb2 templates/7sbreaking/tetramer/pre_tetramer_1bond_three.txt
molecule            mPosttetramer1Bb2 templates/7sbreaking/tetramer/post_tetramer_1bond_three.txt
molecule            mPretetramer1Bb22 templates/7sbreaking/tetramer/pre_tetramer_1bond_three2.txt
molecule            mPosttetramer1Bb22 templates/7sbreaking/tetramer/post_tetramer_1bond_three2.txt
molecule            mPretetramer1Bb4 templates/7sbreaking/tetramer/pre_tetramer_1bond_four.txt
molecule            mPosttetramer1Bb4 templates/7sbreaking/tetramer/post_tetramer_1bond_four.txt
molecule            mPretetramer1Bb42 templates/7sbreaking/tetramer/pre_tetramer_1bond_four2.txt
molecule            mPosttetramer1Bb42 templates/7sbreaking/tetramer/post_tetramer_1bond_four2.txt
molecule            mPretetramer2Bb1 templates/7sbreaking/tetramer/pre_tetramer_2bond.txt
molecule            mPosttetramer2Bb1 templates/7sbreaking/tetramer/post_tetramer_2bond.txt
molecule            mPretetramer2Bb2 templates/7sbreaking/tetramer/pre_tetramer_2bond2.txt
molecule            mPosttetramer2Bb2 templates/7sbreaking/tetramer/post_tetramer_2bond2.txt
molecule            mPretetramerb4 templates/7sbreaking/tetramer/pre_tetramer_four.txt
molecule            mPosttetramerb4 templates/7sbreaking/tetramer/post_tetramer_four.txt

# Nucleation reactions molecular templates
molecule            mPreNucleation templates/Nucleation/pre_Nucleation.txt
molecule            mPostNucleation templates/Nucleation/post_Nucleation.txt


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

variable C_bonds_nc1 equal 0
variable C_bonds_nc1_Cu equal 0
variable C_bonds_7s1 equal 0
variable C_bonds_7s2 equal 0
variable C_bonds_7s3 equal 0
variable C_bonds_7s4 equal 0
variable C_bonds_7s5 equal 0
variable C_bonds_7s_Cu1 equal 0
variable C_bonds_7s_Cu2 equal 0
variable C_bonds_7s_Cu3 equal 0
variable C_bonds_7s_Cu4 equal 0
variable C_bonds_7s_Cu5 equal 0
variable B_bonds_nc1 equal 0
variable B_bonds_nc1_Cu equal 0
variable B_bonds_7s equal 0
variable B_bonds_7s_Cu equal 0

#CREATE VELOCITIES USING MAXWELL-BOLTZMANN DISTRIBUTION
velocity 	all create ${mytemp} ${vseed} dist gaussian mom yes rot yes 
run 0
velocity	all scale ${mytemp}

fix		fix_print all print ${thermodump} "${step} ${etot} ${ke} ${peBond} ${peAngle} ${temp} ${press} ${StressX} ${StressY} ${StressZ} ${Fbond} ${C_bonds_nc1} ${C_bonds_nc1_Cu} ${C_bonds_7s1} ${C_bonds_7s2} ${C_bonds_7s3} ${C_bonds_7s4} ${C_bonds_7s5} ${C_bonds_7s_Cu1} ${C_bonds_7s_Cu2} ${C_bonds_7s_Cu3} ${C_bonds_7s_Cu4} ${C_bonds_7s_Cu5} ${B_bonds_nc1} ${B_bonds_nc1_Cu} ${B_bonds_7s} ${B_bonds_7s_Cu}" file ${thermofile} screen no title "step etot ke peBond peAngle temp press stressX stressY stressZ AvBondForce CreateNC1 CumuCreateNC1 Create7s1 Create7s2 Create7s3 Create7s4 Create7s5 CumuCreate7s1 CumuCreate7s2 CumuCreate7s3 CumuCreate7s4 CumuCreate7s5 BreakNC1 CumuBreakNC1 Break7s CumuBreak7s"

#Save configurations in linear time
dump		2 all custom ${trajdump} dumplin/dump.${simname}.*.lammpstrj type mol xu yu zu
dump_modify	2 sort id
dump_modify	2 pad 10
dump_modify	2 format line "%d %d %.8f %.8f %.8f"

#Save bond information
dump 		3 all local ${trajdump} dumplin_bonds/bonds.${simname}.* index c_1[1] c_1[2] c_1[3] #bond index, id of atom 1 and 2 in the bond, bond type
dump_modify	3 format line "%d %0.0f %0.0f %0.0f"
dump_modify	3 pad 10

#Save angle information
#dump 		4 all local ${trajdump} dumplin_angles/angles.${simname}.* index c_4[1] c_4[2] c_4[3] c_4[4] #id of atom 1 and 2 in the bond 
#dump_modify	4 format line "%d %0.0f %0.0f %0.0f %0.0f"
#dump_modify	4 pad 10

restart         ${trestart} restart/${simname}.*.restart

#run ${tmix}

#write_restart	restart.postmix

neigh_modify exclude molecule/inter 7s_trimer
reset_atoms mol 7s

fix             freact all bond/react stabilization no reset_mol_ids no &
                react rdimer all 1 0.900000 3.300000 mPredimer mPostdimer templates/dimer/map_dimer.txt prob 1.0 ${lseed} &
                react rdimer_1 all 1 0.900000 3.300000 mPredimer_1 mPostdimer_1 templates/dimer/map_dimer_1.txt prob 1.0 ${lseed} &
                react rdimer_2 all 1 0.900000 3.300000 mPredimer_2 mPostdimer_2 templates/dimer/map_dimer_2.txt prob 1.0 ${lseed} &
                react rdimer_3 all 1 0.900000 3.300000 mPredimer_3 mPostdimer_3 templates/dimer/map_dimer_3.txt prob 1.0 ${lseed} &
                react rtrimer all 1 0.900000 3.300000 mPretrimer mPosttrimer templates/trimer/map_trimer.txt prob 1.0 ${lseed} &
                react rtrimer_1 all 1 0.900000 3.300000 mPretrimer_1 mPosttrimer_1 templates/trimer/map_trimer_1.txt prob 1.0 ${lseed} &
                react rtrimer_2 all 1 0.900000 3.300000 mPretrimer_2 mPosttrimer_2 templates/trimer/map_trimer_1.txt prob 1.0 ${lseed} &
                react rtrimer_3 all 1 0.900000 3.300000 mPretrimer_3 mPosttrimer_3 templates/trimer/map_trimer_1.txt prob 1.0 ${lseed} &
                react rtrimer_4 all 1 0.900000 3.300000 mPretrimer_4 mPosttrimer_4 templates/trimer/map_trimer_4.txt prob 1.0 ${lseed} &
                react rtrimer_5 all 1 0.900000 3.300000 mPretrimer_5 mPosttrimer_5 templates/trimer/map_trimer_4.txt prob 1.0 ${lseed} &
                react rtrimer_6 all 1 0.900000 3.300000 mPretrimer_6 mPosttrimer_6 templates/trimer/map_trimer_4.txt prob 1.0 ${lseed} &
                react rtrimer_7 all 1 0.900000 3.300000 mPretrimer_7 mPosttrimer_7 templates/trimer/map_trimer_7.txt prob 1.0 ${lseed} &
                react rtrimer_1_1 all 1 0.900000 3.300000 mPretrimer_1_1 mPosttrimer_1_1 templates/trimer/map_trimer.txt prob 1.0 ${lseed} &
                react rtrimer_1_2 all 1 0.900000 3.300000 mPretrimer_1_2 mPosttrimer_1_2 templates/trimer/map_trimer_1.txt prob 1.0 ${lseed} &
                react rtrimer_1_3 all 1 0.900000 3.300000 mPretrimer_1_3 mPosttrimer_1_3 templates/trimer/map_trimer.txt prob 1.0 ${lseed} &
                react rtrimer_1_4 all 1 0.900000 3.300000 mPretrimer_1_4 mPosttrimer_1_4 templates/trimer/map_trimer_1.txt prob 1.0 ${lseed} &
                react rtrimer_1_5 all 1 0.900000 3.300000 mPretrimer_1_5 mPosttrimer_1_5 templates/trimer/map_trimer.txt prob 1.0 ${lseed} &
                react rtrimer_1_6 all 1 0.900000 3.300000 mPretrimer_1_6 mPosttrimer_1_6 templates/trimer/map_trimer_1.txt prob 1.0 ${lseed}&
                react rtrimer_F all 1 0.900000 3.300000 mPretrimer_F mPosttrimer_F templates/trimerFive/map_trimerF.txt prob 1.0 ${lseed} &
                react rtrimer_F2 all 1 0.900000 3.300000 mPretrimer_F2 mPosttrimer_F2 templates/trimerFive/map_trimerF_2.txt prob 1.0 ${lseed} &
                react rtrimer_F3 all 1 0.900000 3.300000 mPretrimer_F3 mPosttrimer_F3 templates/trimerFive/map_trimerF_3.txt prob 1.0 ${lseed} &
                react rtrimer_F4 all 1 0.900000 3.300000 mPretrimer_F4 mPosttrimer_F4 templates/trimerFive/map_trimerF_4.txt prob 1.0 ${lseed} &
                react rtrimer_F_1 all 1 0.900000 3.300000 mPretrimer_F_1 mPosttrimer_F_1 templates/trimerFive/map_trimerF.txt prob 1.0 ${lseed} &
                react rtrimer_F_12 all 1 0.900000 3.300000 mPretrimer_F_12 mPosttrimer_F_12 templates/trimerFive/map_trimerF_2.txt prob 1.0 ${lseed}&
                react rtetramer_zero all 1 0.900000 3.300000 mPretetramer_zero mPosttetramer_zero templates/tetramer/map_tetramer_zero.txt prob 1.0 ${lseed} &
                react rtetramer_one1 all 1 0.900000 3.300000 mPretetramer_one_1 mPosttetramer_one_1 templates/tetramer/map_tetramer_one.txt prob 1.0 ${lseed} &
                react rtetramer_one2 all 1 0.900000 3.300000 mPretetramer_one_2 mPosttetramer_one_2 templates/tetramer/map_tetramer_one.txt prob 1.0 ${lseed} &
                react rtetramer_one3 all 1 0.900000 3.300000 mPretetramer_one_3 mPosttetramer_one_3 templates/tetramer/map_tetramer_one.txt prob 1.0 ${lseed} &
                react rtetramer_two1 all 1 0.900000 3.300000 mPretetramer_two_1 mPosttetramer_two_1 templates/tetramer/map_tetramer_two.txt prob 1.0 ${lseed} &
                react rtetramer_two2 all 1 0.900000 3.300000 mPretetramer_two_2 mPosttetramer_two_2 templates/tetramer/map_tetramer_two.txt prob 1.0 ${lseed} &
                react rtetramer_two3 all 1 0.900000 3.300000 mPretetramer_two_3 mPosttetramer_two_3 templates/tetramer/map_tetramer_two.txt prob 1.0 ${lseed} &
                react rtetramer_two4 all 1 0.900000 3.300000 mPretetramer_two_4 mPosttetramer_two_4 templates/tetramer/map_tetramer_two.txt prob 1.0 ${lseed} &
                react rtetramer_three1 all 1 0.900000 3.300000 mPretetramer_three_1 mPosttetramer_three_1 templates/tetramer/map_tetramer_three.txt prob 1.0 ${lseed} &
                react rtetramer_three2 all 1 0.900000 3.300000 mPretetramer_three_2 mPosttetramer_three_2 templates/tetramer/map_tetramer_three.txt prob 1.0 ${lseed} &
                react rtetramer_three3 all 1 0.900000 3.300000 mPretetramer_three_3 mPosttetramer_three_3 templates/tetramer/map_tetramer_three.txt prob 1.0 ${lseed} &
                react rtetramer_four all 1 0.900000 3.300000 mPretetramer_four mPosttetramer_four templates/tetramer/map_tetramer_four.txt prob 1.0 ${lseed} &
                react rtetramer_1bond_two1 all 1 0.900000 3.300000 mPretetramer_1bond_two_1 mPosttetramer_1bond_two_1 templates/tetramer/map_tetramer_1bond_two.txt prob 1.0 ${lseed} &
                react rtetramer_1bond_two2 all 1 0.900000 3.300000 mPretetramer_1bond_two_2 mPosttetramer_1bond_two_2 templates/tetramer/map_tetramer_1bond_two.txt prob 1.0 ${lseed} &
                react rtetramer_1bond_two3 all 1 0.900000 3.300000 mPretetramer_1bond_two_3 mPosttetramer_1bond_two_3 templates/tetramer/map_tetramer_1bond_two.txt prob 1.0 ${lseed} &
                react rtetramer_1bond_three1 all 1 0.900000 3.300000 mPretetramer_1bond_three_1 mPosttetramer_1bond_three_1 templates/tetramer/map_tetramer_1bond_three.txt prob 1.0 ${lseed} &
                react rtetramer_1bond_three2 all 1 0.900000 3.300000 mPretetramer_1bond_three_2 mPosttetramer_1bond_three_2 templates/tetramer/map_tetramer_1bond_three.txt prob 1.0 ${lseed} &
                react rtetramer_1bond_three3 all 1 0.900000 3.300000 mPretetramer_1bond_three_3 mPosttetramer_1bond_three_3 templates/tetramer/map_tetramer_1bond_three.txt prob 1.0 ${lseed} &
                react rtetramer_1bond_three4 all 1 0.900000 3.300000 mPretetramer_1bond_three_4 mPosttetramer_1bond_three_4 templates/tetramer/map_tetramer_1bond_three.txt prob 1.0 ${lseed} &
                react rtetramer_1bond_three5 all 1 0.900000 3.300000 mPretetramer_1bond_three_5 mPosttetramer_1bond_three_5 templates/tetramer/map_tetramer_1bond_three.txt prob 1.0 ${lseed} &
                react rtetramer_1bond_four1 all 1 0.900000 3.300000 mPretetramer_1bond_four_1 mPosttetramer_1bond_four_1 templates/tetramer/map_tetramer_1bond_four.txt prob 1.0 ${lseed} &
                react rtetramer_1bond_four2 all 1 0.900000 3.300000 mPretetramer_1bond_four_2 mPosttetramer_1bond_four_2 templates/tetramer/map_tetramer_1bond_four.txt prob 1.0 ${lseed} &
                react rtetramer_1bond_four3 all 1 0.900000 3.300000 mPretetramer_1bond_four_3 mPosttetramer_1bond_four_3 templates/tetramer/map_tetramer_1bond_four.txt prob 1.0 ${lseed} &
                react rtetramer_2bond_1 all 1 0.900000 3.300000 mPretetramer_2bond_1 mPosttetramer_2bond_1 templates/tetramer/map_tetramer_2bond.txt prob 1.0 ${lseed} &
                react rtetramer_2bond_2 all 1 0.900000 3.300000 mPretetramer_2bond_2 mPosttetramer_2bond_2 templates/tetramer/map_tetramer_2bond.txt prob 1.0 ${lseed} &
                react rdimer1 all 1 0.100000 3.300000 mPredimer1Br mPostdimer1Br templates/7sbreaking/map_dimer1.txt prob 1.0 ${lseed} & 
                react rdimer2 all 1 0.100000 3.300000 mPredimer2Br mPostdimer2Br templates/7sbreaking/map_dimer2.txt prob 1.0 ${lseed} &
                react rdimer1_6 all 1 0.100000 4.300000 mPredimer1Br_6 mPostdimer1Br templates/7sbreaking/map_dimer1.txt prob 1.0 ${lseed} &
                react rdimer2_6 all 1 0.100000 4.300000 mPredimer2Br_6 mPostdimer2Br templates/7sbreaking/map_dimer2.txt prob 1.0 ${lseed} &
                react rdimer1_4 all 1 0.100000 4.300000 mPredimer1Br_4 mPostdimer1Br templates/7sbreaking/map_dimer1.txt prob 1.0 ${lseed} &
                react rdimer2_4 all 1 0.100000 4.300000 mPredimer2Br_4 mPostdimer2Br templates/7sbreaking/map_dimer2.txt prob 1.0 ${lseed} &
                react rdimer1_5 all 1 0.100000 4.300000 mPredimer1Br_5 mPostdimer1Br templates/7sbreaking/map_dimer1.txt prob 1.0 ${lseed} &
                react rdimer2_5 all 1 0.100000 4.300000 mPredimer2Br_5 mPostdimer2Br templates/7sbreaking/map_dimer2.txt prob 1.0 ${lseed} &
                react rdimer1_7 all 1 0.100000 4.300000 mPredimer1Br_7 mPostdimer1Br templates/7sbreaking/map_dimer1.txt prob 1.0 ${lseed} &
                react rdimer2_7 all 1 0.100000 4.300000 mPredimer2Br_7 mPostdimer2Br templates/7sbreaking/map_dimer2.txt prob 1.0 ${lseed} &
                react rdimerb all 1 0.900000 3.300000 mPredimerb mPostdimer templates/dimer/map_dimer.txt prob 1.0 ${lseed} &
                react rdimerb_1 all 1 0.900000 3.300000 mPredimerb_1 mPostdimer_1 templates/dimer/map_dimer_1.txt prob 1.0 ${lseed} &
                react rdimerb_2 all 1 0.900000 3.300000 mPredimerb_2 mPostdimer_2 templates/dimer/map_dimer_2.txt prob 1.0 ${lseed} &
                react rdimerb_3 all 1 0.900000 3.300000 mPredimerb_3 mPostdimer_3 templates/dimer/map_dimer_3.txt prob 1.0 ${lseed} &
                react rtrimer1 all 1 0.900000 3.300000 mPretrimer1Br mPosttrimer1Br templates/7sbreaking/map_trimerF.txt prob 1.0 ${lseed} &
                react rtrimer2 all 1 0.900000 3.300000 mPretrimer2Br mPosttrimer2Br templates/7sbreaking/map_trimerF_2.txt prob 1.0 ${lseed} &
                react rtrimer3 all 1 0.900000 3.300000 mPretrimer3Br mPosttrimer3Br templates/7sbreaking/map_trimerF_3.txt prob 1.0 ${lseed} &
                react rtrimer4 all 1 0.900000 3.300000 mPretrimer4Br mPosttrimer4Br templates/7sbreaking/map_trimerF_4.txt prob 1.0 ${lseed} &
                react rtrimer5 all 1 0.900000 3.300000 mPretrimer5Br mPosttrimer5Br templates/7sbreaking/map_trimerF_2.txt prob 1.0 ${lseed} &
                react rtrimer6 all 1 0.900000 3.300000 mPretrimer6Br mPosttrimer6Br templates/7sbreaking/map_trimerF.txt prob 1.0 ${lseed} &
                react rtrimerb_2 all 1 0.900000 3.300000 mPretrimerb_2 mPosttrimer_2 templates/trimer/map_trimer_1.txt prob 1.0 ${lseed} &
                react rtrimerb_3 all 1 0.900000 3.300000 mPretrimerb_3 mPosttrimer_3 templates/trimer/map_trimer_1.txt prob 1.0 ${lseed} &
                react rtrimerb_4 all 1 0.900000 3.300000 mPretrimerb_4 mPosttrimer_4 templates/trimer/map_trimer_4.txt prob 1.0 ${lseed} &
                react rtrimerb_5 all 1 0.900000 3.300000 mPretrimerb_5 mPosttrimer_5 templates/trimer/map_trimer_4.txt prob 1.0 ${lseed} &
                react rtrimerb_7 all 1 0.900000 3.300000 mPretrimerb_7 mPosttrimer_7 templates/trimer/map_trimer_7.txt prob 1.0 ${lseed} &
                react rtrimerb all 1 0.900000 3.300000 mPretrimerb mPosttrimer templates/trimer/map_trimer.txt prob 1.0 ${lseed} &
                react rtrimerb_1_1 all 1 0.900000 3.300000 mPretrimerb_1_1 mPosttrimer_1_1 templates/trimer/map_trimer.txt prob 1.0 ${lseed} &
                react rtrimerb_1_3 all 1 0.900000 3.300000 mPretrimerb_1_3 mPosttrimer_1_3 templates/trimer/map_trimer.txt prob 1.0 ${lseed} &
                react rtrimerb_1_6 all 1 0.900000 3.300000 mPretrimerb_1_6 mPosttrimer_1_6 templates/trimer/map_trimer_1.txt prob 1.0 ${lseed} &
                react rtrimerb_1_4 all 1 0.900000 3.300000 mPretrimerb_1_4 mPosttrimer_1_4 templates/trimer/map_trimer_1.txt prob 1.0 ${lseed} &
                react rtetramerB_zero all 1 0.900000 3.300000 mPretetramerb mPosttetramerb templates/tetramer/map_tetramer_zero.txt prob 1.0 ${lseed} &
                react rtetramerB_one all 1 0.900000 3.300000 mPretetramerb1 mPosttetramerb1 templates/tetramer/map_tetramer_one.txt prob 1.0 ${lseed} &
                react rtetramerB_two all 1 0.900000 3.300000 mPretetramerb2 mPosttetramerb2 templates/tetramer/map_tetramer_two.txt prob 1.0 ${lseed} &
                react rtetramerB_two2 all 1 0.900000 3.300000 mPretetramerb22 mPosttetramerb22 templates/tetramer/map_tetramer_two.txt prob 1.0 ${lseed} &
                react rtetramer_1bondB_one all 1 0.900000 3.300000 mPretetramer1Bb1 mPosttetramer1Bb1 templates/tetramer/map_tetramer_1bond_two.txt prob 1.0 ${lseed} &
                react rtetramer_1bondB_one2 all 1 0.900000 3.300000 mPretetramer1Bb12 mPosttetramer1Bb12 templates/tetramer/map_tetramer_1bond_two.txt prob 1.0 ${lseed} &
                react rtetramerB_three all 1 0.900000 3.300000 mPretetramerb3 mPosttetramerb3 templates/tetramer/map_tetramer_three.txt prob 1.0 ${lseed} &
                react rtetramer_1bondB_two all 1 0.900000 3.300000 mPretetramer1Bb2 mPosttetramer1Bb2 templates/tetramer/map_tetramer_1bond_three.txt prob 1.0 ${lseed} &
                react rtetramer_1bondB_two2 all 1 0.900000 3.300000 mPretetramer1Bb22 mPosttetramer1Bb22 templates/tetramer/map_tetramer_1bond_three.txt prob 1.0 ${lseed} &
                react rtetramerB_four all 1 0.900000 3.300000 mPretetramer1Bb4 mPretetramer1Bb4 templates/tetramer/map_tetramer_1bond_four.txt prob 1.0 ${lseed} &
                react rtetramerB_four2 all 1 0.900000 3.300000 mPretetramer1Bb42 mPretetramer1Bb42 templates/tetramer/map_tetramer_1bond_four.txt prob 1.0 ${lseed} &
                react rtetramer_2bondB all 1 0.900000 3.300000 mPretetramer2Bb1 mPosttetramer2Bb1 templates/tetramer/map_tetramer_2bond.txt prob 1.0 ${lseed} &
                react rtetramer_2bondB2 all 1 0.900000 3.300000 mPretetramer2Bb2 mPosttetramer2Bb2 templates/tetramer/map_tetramer_2bond.txt prob 1.0 ${lseed} &
                react rtetramerB_four all 1 0.900000 3.300000 mPretetramerb4 mPosttetramerb4 templates/tetramer/map_tetramer_four.txt prob 1.0 ${lseed} &


special_bonds lj 0.0 1.0 1.0
''')
    for n in range(int(NumFixCycl)):
        fo.write(
        '''
special_bonds lj 0.0 1.0 1.0
reset_atoms mol 7s
group 7s_trimer  type 4
group 7s_tetramer  type 6 7
neigh_modify exclude molecule/inter 7s_trimer 
neigh_modify exclude molecule/inter 7s_tetramer 
fix 		fix_bonds_nc1 	nc1 bond/create  ${nevery_Cr} 1 1 ${make_dist_nc1} 2 iparam ${max_bonds_nc1} 1 jparam ${max_bonds_nc1} 1 prob ${make_prob_nc1} 49829 atype 1 
fix 		fix_bonds_7s	7s 	bond/create  ${nevery_Cr} 2 2 ${make_dist_7s} 3 iparam ${max_bonds_7s} 7 jparam ${max_bonds_7s} 7 prob ${make_prob_7s} 29829 atype 2 
fix 		fix_bonds_7s_2	7s bond/create  ${nevery_Cr} 2 3 ${make_dist_7s} 3 iparam ${max_bonds_7s} 7 jparam ${max_bonds_7s} 7 prob ${make_prob_7s} 29829 atype 2 
fix 		fix_bonds_7s_3	7s bond/create ${nevery_Cr2} 4 4 ${make_dist_7s} 3 iparam ${max_bonds_7s} 7 jparam ${max_bonds_7s} 7 prob ${make_prob_7s} 29829 atype 2 
fix 		fix_bonds_7s_4	7s bond/create  ${nevery_Cr} 2 5 ${make_dist_7s} 3 iparam ${max_bonds_7s} 2 jparam ${max_bonds_7s} 5 prob ${make_prob_7s} 29829 atype 2 
fix 		fix_bonds_7s_5	7s bond/create  ${nevery_Cr2} 6 6 ${make_dist_7s} 3 iparam ${max_bonds_7s} 6 jparam ${max_bonds_7s} 6 prob ${make_prob_7s} 29829 atype 2 


variable C_bonds_nc1 equal f_fix_bonds_nc1[1]
variable C_bonds_nc1_Cu equal f_fix_bonds_nc1[2]
variable C_bonds_7s1 equal f_fix_bonds_7s[1]
variable C_bonds_7s2 equal f_fix_bonds_7s_2[1]
variable C_bonds_7s3 equal f_fix_bonds_7s_3[1]
variable C_bonds_7s4 equal f_fix_bonds_7s_4[1]
variable C_bonds_7s5 equal f_fix_bonds_7s_5[1]
variable C_bonds_7s_Cu1 equal f_fix_bonds_7s[2] 
variable C_bonds_7s_Cu2 equal f_fix_bonds_7s_2[2] 
variable C_bonds_7s_Cu3 equal f_fix_bonds_7s_3[2] 
variable C_bonds_7s_Cu4 equal f_fix_bonds_7s_4[2] 
variable C_bonds_7s_Cu5 equal f_fix_bonds_7s_5[2] 
variable B_bonds_nc1 equal 0
variable B_bonds_nc1_Cu equal 0
variable B_bonds_7s equal 0
variable B_bonds_7s_Cu equal 0


run		${Tbetween}

#Sample with fast dumping for bond topology changes
dump 		4 all local ${NevSample} dumplin_bondsSample/bonds.${simname}.* index c_1[1] c_1[2] c_1[3] #bond index, id of atom 1 and 2 in the bond, bond type
dump_modify	4 format line "%d %0.0f %0.0f %0.0f"
dump_modify	4 pad 10

run		${tsampleTop}

undump 4

run		${Tbetween}

unfix fix_bonds_nc1
unfix fix_bonds_7s
unfix fix_bonds_7s_2
unfix fix_bonds_7s_3
unfix fix_bonds_7s_4
unfix fix_bonds_7s_5

fix 		fix_bonds_break_nc1 	nc1 bond/break  ${nevery_Br} 2 ${break_dist_nc1} prob ${break_prob_nc1} 49829
fix 		fix_bonds_break_7s 	7s bond/break  ${nevery_Br} 3 ${break_dist_nc1} prob ${break_prob_7s} 49825

variable C_bonds_nc1 equal 0
variable C_bonds_nc1_Cu equal 0
variable C_bonds_7s1 equal 0
variable C_bonds_7s2 equal 0
variable C_bonds_7s3 equal 0
variable C_bonds_7s4 equal 0
variable C_bonds_7s5 equal 0
variable C_bonds_7s_Cu1 equal 0
variable C_bonds_7s_Cu2 equal 0
variable C_bonds_7s_Cu3 equal 0
variable C_bonds_7s_Cu4 equal 0
variable C_bonds_7s_Cu5 equal 0
variable B_bonds_nc1 equal f_fix_bonds_break_nc1[1]
variable B_bonds_nc1_Cu equal f_fix_bonds_break_nc1[2]
variable B_bonds_7s equal f_fix_bonds_break_7s[1]
variable B_bonds_7s_Cu equal f_fix_bonds_break_7s[2]

run		${Tbetween}

#Sample with fast dumping for bond topology changes
dump 		4 all local ${NevSample} dumplin_bondsSample/bonds.${simname}.* index c_1[1] c_1[2] c_1[3] #bond index, id of atom 1 and 2 in the bond, bond type
dump_modify	4 format line "%d %0.0f %0.0f %0.0f"
dump_modify	4 pad 10

run		${tsampleTop}

undump 4

run		${Tbetween}
    
unfix fix_bonds_break_nc1
unfix fix_bonds_break_7s

write_data	data.*
    ''')
    fo.write('''
write_restart	restart.postmakebonds\n''')
    fo.write('fix                    fDeform all deform 1 x delta '+str(-XStretch*0.5)+' '+str(XStretch*0.5)+ ' remap x\n')
    fo.write('''
    
variable C_bonds_nc1 equal 0
variable C_bonds_nc1_Cu equal 0
variable C_bonds_7s1 equal 0
variable C_bonds_7s2 equal 0
variable C_bonds_7s3 equal 0
variable C_bonds_7s4 equal 0
variable C_bonds_7s5 equal 0
variable C_bonds_7s_Cu1 equal 0
variable C_bonds_7s_Cu2 equal 0
variable C_bonds_7s_Cu3 equal 0
variable C_bonds_7s_Cu4 equal 0
variable C_bonds_7s_Cu5 equal 0
variable B_bonds_nc1 equal 0
variable B_bonds_nc1_Cu equal 0
variable B_bonds_7s equal 0
variable B_bonds_7s_Cu equal 0
    
run                    ${tstretch}
unfix                  fDeform
write_restart	restart.poststretch\n''')

    for n in range(int(NumFixCycl*5)):
        fo.write(
        '''
special_bonds lj 0.0 1.0 1.0
reset_atoms mol 7s
group 7s_trimer  type 4
group 7s_tetramer  type 6
neigh_modify exclude molecule/inter 7s_trimer 
neigh_modify exclude molecule/inter 7s_tetramer 
fix 		fix_bonds_nc1 	nc1 bond/create  ${nevery_Cr} 1 1 ${make_dist_nc1} 2 iparam ${max_bonds_nc1} 1 jparam ${max_bonds_nc1} 1 prob ${make_prob_nc1} 49829 atype 1 
fix 		fix_bonds_7s	7s 	bond/create  ${nevery_Cr} 2 2 ${make_dist_7s} 3 iparam ${max_bonds_7s} 7 jparam ${max_bonds_7s} 7 prob ${make_prob_7s} 29829 atype 2 
fix 		fix_bonds_7s_2	7s bond/create  ${nevery_Cr} 2 3 ${make_dist_7s} 3 iparam ${max_bonds_7s} 7 jparam ${max_bonds_7s} 7 prob ${make_prob_7s} 29829 atype 2 
fix 		fix_bonds_7s_3	7s bond/create ${nevery_Cr2} 4 4 ${make_dist_7s} 3 iparam ${max_bonds_7s} 7 jparam ${max_bonds_7s} 7 prob ${make_prob_7s} 29829 atype 2 
fix 		fix_bonds_7s_4	7s bond/create  ${nevery_Cr} 2 5 ${make_dist_7s} 3 iparam ${max_bonds_7s} 2 jparam ${max_bonds_7s} 5 prob ${make_prob_7s} 29829 atype 2 
fix 		fix_bonds_7s_5	7s bond/create  ${nevery_Cr2} 6 6 ${make_dist_7s} 3 iparam ${max_bonds_7s} 6 jparam ${max_bonds_7s} 6 prob ${make_prob_7s} 29829 atype 2 

variable C_bonds_nc1 equal f_fix_bonds_nc1[1]
variable C_bonds_nc1_Cu equal f_fix_bonds_nc1[2]
variable C_bonds_7s1 equal f_fix_bonds_7s[1]
variable C_bonds_7s2 equal f_fix_bonds_7s_2[1]
variable C_bonds_7s3 equal f_fix_bonds_7s_3[1]
variable C_bonds_7s4 equal f_fix_bonds_7s_4[1]
variable C_bonds_7s5 equal f_fix_bonds_7s_5[1]
variable C_bonds_7s_Cu1 equal f_fix_bonds_7s[2] 
variable C_bonds_7s_Cu2 equal f_fix_bonds_7s_2[2] 
variable C_bonds_7s_Cu3 equal f_fix_bonds_7s_3[2] 
variable C_bonds_7s_Cu4 equal f_fix_bonds_7s_4[2] 
variable C_bonds_7s_Cu5 equal f_fix_bonds_7s_5[2] 
variable B_bonds_nc1 equal 0
variable B_bonds_nc1_Cu equal 0
variable B_bonds_7s equal 0
variable B_bonds_7s_Cu equal 0


run		${tcycleRel}
unfix fix_bonds_nc1
unfix fix_bonds_7s
unfix fix_bonds_7s_2
unfix fix_bonds_7s_3
unfix fix_bonds_7s_4
unfix fix_bonds_7s_5
special_bonds lj 0.0 1.0 1.0
reset_atoms mol 7s
group 7s_trimer  type 4
group 7s_tetramer  type 6
group  			7sBreaking  type 2 3 4 5

neigh_modify exclude molecule/inter 7s_trimer 
neigh_modify exclude molecule/inter 7s_tetramer 

fix 		fix_bonds_break_nc1 	nc1 bond/break  ${nevery_Br} 2 ${break_dist_nc1} prob ${break_prob_nc1} 49829
fix 		fix_bonds_break_7s 	7s bond/break  ${nevery_Br} 3 ${break_dist_nc1} prob ${break_prob_7s} 49825

variable C_bonds_nc1 equal 0
variable C_bonds_nc1_Cu equal 0
variable C_bonds_7s1 equal 0
variable C_bonds_7s2 equal 0
variable C_bonds_7s3 equal 0
variable C_bonds_7s4 equal 0
variable C_bonds_7s5 equal 0
variable C_bonds_7s_Cu1 equal 0
variable C_bonds_7s_Cu2 equal 0
variable C_bonds_7s_Cu3 equal 0
variable C_bonds_7s_Cu4 equal 0
variable C_bonds_7s_Cu5 equal 0
variable B_bonds_nc1 equal f_fix_bonds_break_nc1[1]
variable B_bonds_nc1_Cu equal f_fix_bonds_break_nc1[2]
variable B_bonds_7s equal f_fix_bonds_break_7s[1]
variable B_bonds_7s_Cu equal f_fix_bonds_break_7s[2]

run ${tcycleRel}
    
unfix fix_bonds_break_nc1
unfix fix_bonds_break_7s
    
    ''')
    fo.write('''
write_restart	restart.postrelax\n''')
    




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
#SBATCH -N 4
#SBATCH -n 4
#SBATCH --mem=400M 
#SBATCH --time=72:00:00 
#SBATCH --mail-user=ucapbbm@ucl.ac.uk
#SBATCH --mail-type=END
#SBATCH --export=NONE
#SBATCH --exclude zeta247

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
