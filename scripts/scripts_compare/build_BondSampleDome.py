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
    #${foldernameadd} ${Time_Eq} ${NumFrames} ${Nevery} ${MD} ${BD} ${MP} ${BP} ${seed}
    tEq = str(sys.argv[2])
    tStr = str(sys.argv[3])
    tRel = str(sys.argv[4])
    Frames =  int(sys.argv[5])
    Nevery = str(sys.argv[6])
    tstep = str(sys.argv[7])
    MakeDist = str(sys.argv[8])
    BreakDist = str(sys.argv[9])
    MakeDist_7S = str(sys.argv[10])
    BreakDist_7S = str(sys.argv[11])
    
    MakeProb = str(sys.argv[12])
    MakeProb_7S = str(sys.argv[13])
    FactorMult = float(sys.argv[14])
    BreakProb = str(FactorMult*float(MakeProb))
    BreakProb_7S = str(FactorMult*float(MakeProb_7S))
    
    KangNC1 =str(sys.argv[15])
    PrefMonNum = str(sys.argv[16])
    XStretch = float(sys.argv[17])
    OvaR = str(sys.argv[18])
    
    kappa_density = str(sys.argv[19])
    VolInit = str(sys.argv[20])
    seed  = str(sys.argv[21])
    Nstretch = 20
    tmix = 0 #1e5
    tAssemble = 1e5
    
    #cutoff = np.sqrt((ChemBond-np.log(0.01))/6)+Blen #Cutoff the making probabbilty at Prob = 0.01, with bond constant 6 
    #print("cutoff"+str(cutoff))
    ################ write in.local ################
    print('Writing collagen.in\n')
    outfile2 = workingdir + '/collagen.in'
    fo=open(outfile2,'w')
    fo.write(
        '''#####################################################
#PARAMETERS

#NUMBER OF STEPS OF RUNS (stages)
''')
    fo.write('variable tmix			equal 	'+str(tmix)+'	#NVT steps for molecule mixing\n')
    fo.write('variable tbonds			equal 	'+str(tEq)+'	#NVT steps for bond equilibration\n')
    fo.write('variable tassemble       equal 	'+str(tAssemble)+'	#NVT steps for bond equilibration\n')
    fo.write('variable tstretch			equal 	'+str(tStr)+'	#NVT steps for bond equilibration\n')
    fo.write('variable tstretchcycl			equal 	'+str(float(tStr)/Nstretch)+'	#NVT steps for stretching\n')
    fo.write('variable trelax			equal 	'+str(tRel)+'	#NVT steps for bond equilibration\n')
    fo.write('variable Frames			equal 	'+str(Frames)+'	#NVT steps for bond equilibration\n')
    fo.write('variable FramesMany			equal 	'+str(float(Frames)*10)+'	#NVT steps for bond equilibration\n')
    fo.write(
        '''

#THERMODYNAMIC PARAMETERS
variable mytemp         equal 	1.0	    	#Temperature

#INFORMATION PRINTING PERIODS
variable ttotal     equal   $(v_tmix+v_tbonds+v_tstretch+v_trelax)		#total timesteps
variable thermodump     equal   $(floor(v_ttotal/5e4))		#Thermo info dump period
variable trajdump       equal   $(floor(v_ttotal/v_Frames))    	#Config info dump period
variable trajdumpFast       equal   $(floor(v_ttotal/v_FramesMany))    	#Config info dump period
variable trestart       equal   $(floor(v_ttotal/20))    	#Restart configuration saving period\n''')

    fo.write(
        '''
#OTHER SIMULATION PARAMETERS
''')
    fo.write('variable tstep        	equal   '+str(tstep)+'   	#Integration time step for first part of simulation')
    fo.write(
        '''
variable damp   		equal   0.1	   	#Langevin thermostat damping coefficient [F_friction = -(m/damp)*v]
variable neigh_cutoff	equal	0.4  		#Neighbor list cutoff
variable lj_cutoff      equal   1.122462	#For WCA: 1.122462048 

#new variables Fernanda
variable RminF  equal  0.
variable RmaxFB  equal  30

variable NeveryFast equal 1

''')
    fo.write('variable RmaxF  equal '+str(MakeDist)+' #Distance to make a bond\n')
    fo.write('variable probF  equal '+str(MakeProb)+' #Distance to make a bond\n')
    fo.write('variable RminFB  equal '+str(BreakDist)+' #Distance to make a bond\n')
    fo.write('variable probFB  equal '+str(BreakProb)+' #Distance to make a bond\n')
    fo.write('variable RmaxF_7S  equal '+str(MakeDist_7S)+' #Distance to make a bond\n')
    fo.write('variable probF_7S  equal '+str(MakeProb_7S)+' #Distance to make a bond\n')
    fo.write('variable RminFB_7S  equal '+str(BreakDist_7S)+' #Distance to make a bond\n')
    fo.write('variable probFB_7S  equal '+str(BreakProb_7S)+' #Distance to make a bond\n')
    fo.write('variable N_preferred equal '+str(PrefMonNum)+' #Preferred number of monomres\n')
    fo.write('variable OvaR equal '+str(OvaR)+' #Space needed to add monomer\n')
    fo.write('variable KappaD equal '+str(kappa_density)+' #Strength of density conservation\n')
    fo.write('variable volInit equal '+str(VolInit)+' #Initial Volume\n')
    
    fo.write(
        '''

variable lseed_base equal 13579
variable lseed1 equal ${lseed_base}+1 #NC1 creation
variable lseed2 equal ${lseed_base}+2 #dimer creation
variable lseed3 equal ${lseed_base}+2 #line-trimer creation
variable lseed4 equal ${lseed_base}+2 #triangle-trimer creation
variable lseed5 equal ${lseed_base}+2 #open-tetramer creation
variable lseed6 equal ${lseed_base}+2 #rhomboid-tetramer creation
variable lseed7 equal ${lseed_base}+2 #squared-tetramer creation
variable lseed8 equal ${lseed_base}+2 #squared-tetramer break
variable lseed9 equal ${lseed_base}+2 #rhomboid-tetramer break
variable lseed10 equal ${lseed_base}+2 #open-tetramer break
variable lseed11 equal ${lseed_base}+2 #triangle-trimer break
variable lseed12 equal ${lseed_base}+2 #line-trimer break
variable lseed13 equal ${lseed_base}+2 #dimer break
variable lseed14 equal ${lseed_base}+2 #NC1 break
''')
    fo.write('variable NeverySlow equal '+str(Nevery)+' #how often we attempt reactions\n')
    fo.write('variable k_nc1			equal	6.0 		#Elastic constant A type bond\n')
    fo.write('variable r0_nc1			equal	0.5 		#Rest length A type bond\n\n')

    fo.write('variable k_7s			equal	6.0 		#Elastic constant B type bond\n')
    fo.write('variable r0_7s			equal	0.5 		#Rest length B type bond\n\n')

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
boundary        p p p      	#periodic boundary conditions
atom_style     	molecular

#BOND STYLE: 
bond_style      harmonic

#ANGLE STYLE: 
angle_style      harmonic

#TEMPORARY PAIR STYLE TO START SIMULATION
pair_style		lj/cut 1.0

#Starting from an initial configuration WITHOUT BONDS (format lammpsdata):
#read_data       ${fname} extra/bond/per/atom 10 extra/special/per/atom 100 extra/angle/per/atom 20
read_restart ${fname}
#reset_timestep  0


#GROUPS
#group  			nc1 type 1
#group  			7s  type 2 3 4 5 6 7
#group  			7sBreaking  type 2 3 4 5
#group  			7s_trimer  type 4
group  		    mobile type 1 2 3 4 5 6 7 8 9
group PureMon type 9
#BONDS COEFFICIENTS

bond_coeff     	1 	100.0	3.0
bond_coeff     	2 ${k_nc1} ${r0_nc1}
bond_coeff		3 ${k_7s} ${r0_7s}
bond_coeff		4 1000.0 2.0

#ANGLE COEFFICIENTS\n
''')
    fo.write('angle_coeff     1 	'+KangNC1+'	180.0\n')
    fo.write('angle_coeff     2 	1.0	155.0\n')
    fo.write('angle_coeff     3 	1.0	60.0\n')
    fo.write(
        '''

#NEIGHBOR LISTS
neighbor		${neigh_cutoff} bin 		#value = skin = extra distance beyond cutoff
neigh_modify	every 10 delay 0 check yes

#PAIR STYLE (WCA)
pair_style     	hybrid/overlay lj/cut 1.0 zero 1.0
pair_coeff      * * lj/cut 0.0 0.0 ${lj_cutoff}
pair_coeff      * * zero ${zero_cutoff}
pair_coeff      11 1 lj/cut 2.0 1.0 1.22
pair_coeff      11 2 lj/cut 2.0 1.0 1.22
pair_coeff      11 3 lj/cut 2.0 1.0 1.22
pair_coeff      11 4 lj/cut 2.0 1.0 1.22
pair_coeff      11 5 lj/cut 2.0 1.0 1.22
pair_coeff      11 6 lj/cut 2.0 1.0 1.22
pair_coeff      11 7 lj/cut 2.0 1.0 1.22
pair_coeff      11 8 lj/cut 2.0 1.0 1.22
pair_coeff      11 9 lj/cut 2.0 1.0 1.22
pair_coeff      11 10 lj/cut 2.0 1.0 1.22

#SET TIME STEP
timestep	${tstep}
reset_atoms id sort yes

####-PRINTING STUFF-##############################


############# MOLECULE TEMPLATES ################

#MOLECULE TEMPLATES: Nc1-2N
molecule            Nc1_2N_01_pre ../../templates/NC1_2N/01_pre.txt
molecule            Nc1_2N_01_post ../../templates/NC1_2N/01_post.txt
molecule            Nc1_2N_02_pre ../../templates/NC1_2N/02_pre.txt
molecule            Nc1_2N_02_post ../../templates/NC1_2N/02_post.txt
molecule            Nc1_2N_03_pre ../../templates/NC1_2N/03_pre.txt
molecule            Nc1_2N_03_post ../../templates/NC1_2N/03_post.txt
molecule            Nc1_2N_04_pre ../../templates/NC1_2N/04_pre.txt
molecule            Nc1_2N_04_post ../../templates/NC1_2N/04_post.txt
molecule            Nc1_2N_05_pre ../../templates/NC1_2N/05_pre.txt
molecule            Nc1_2N_05_post ../../templates/NC1_2N/05_post.txt
molecule            Nc1_2N_06_pre ../../templates/NC1_2N/06_pre.txt
molecule            Nc1_2N_06_post ../../templates/NC1_2N/06_post.txt
molecule            Nc1_2N_07_pre ../../templates/NC1_2N/07_pre.txt
molecule            Nc1_2N_07_post ../../templates/NC1_2N/07_post.txt
molecule            Nc1_2N_08_pre ../../templates/NC1_2N/08_pre.txt
molecule            Nc1_2N_08_post ../../templates/NC1_2N/08_post.txt
molecule            Nc1_2N_09_pre ../../templates/NC1_2N/09_pre.txt
molecule            Nc1_2N_09_post ../../templates/NC1_2N/09_post.txt
molecule            Nc1_2N_10_pre ../../templates/NC1_2N/10_pre.txt
molecule            Nc1_2N_10_post ../../templates/NC1_2N/10_post.txt
molecule            Nc1_2N_11_pre ../../templates/NC1_2N/11_pre.txt
molecule            Nc1_2N_11_post ../../templates/NC1_2N/11_post.txt
#MOLECULE TEMPLATES: Nc1-3N
molecule            Nc1_3N_01_pre ../../templates/NC1_3N/01_pre.txt
molecule            Nc1_3N_01_post ../../templates/NC1_3N/01_post.txt
molecule            Nc1_3N_02_pre ../../templates/NC1_3N/02_pre.txt
molecule            Nc1_3N_02_post ../../templates/NC1_3N/02_post.txt
molecule            Nc1_3N_03_pre ../../templates/NC1_3N/03_pre.txt
molecule            Nc1_3N_03_post ../../templates/NC1_3N/03_post.txt
molecule            Nc1_3N_04_pre ../../templates/NC1_3N/04_pre.txt
molecule            Nc1_3N_04_post ../../templates/NC1_3N/04_post.txt
molecule            Nc1_3N_05_pre ../../templates/NC1_3N/05_pre.txt
molecule            Nc1_3N_05_post ../../templates/NC1_3N/05_post.txt
molecule            Nc1_3N_06_pre ../../templates/NC1_3N/06_pre.txt
molecule            Nc1_3N_06_post ../../templates/NC1_3N/06_post.txt
molecule            Nc1_3N_07_pre ../../templates/NC1_3N/07_pre.txt
molecule            Nc1_3N_07_post ../../templates/NC1_3N/07_post.txt
molecule            Nc1_3N_08_pre ../../templates/NC1_3N/08_pre.txt
molecule            Nc1_3N_08_post ../../templates/NC1_3N/08_post.txt
molecule            Nc1_3N_09_pre ../../templates/NC1_3N/09_pre.txt
molecule            Nc1_3N_09_post ../../templates/NC1_3N/09_post.txt
molecule            Nc1_3N_10_pre ../../templates/NC1_3N/10_pre.txt
molecule            Nc1_3N_10_post ../../templates/NC1_3N/10_post.txt
#MOLECULE TEMPLATES: Nc1-4N
molecule            Nc1_4N_01_pre ../../templates/NC1_4N/01_pre.txt
molecule            Nc1_4N_01_post ../../templates/NC1_4N/01_post.txt
molecule            Nc1_4N_02_pre ../../templates/NC1_4N/02_pre.txt
molecule            Nc1_4N_02_post ../../templates/NC1_4N/02_post.txt
molecule            Nc1_4N_03_pre ../../templates/NC1_4N/03_pre.txt
molecule            Nc1_4N_03_post ../../templates/NC1_4N/03_post.txt
molecule            Nc1_4N_04_pre ../../templates/NC1_4N/04_pre.txt
molecule            Nc1_4N_04_post ../../templates/NC1_4N/04_post.txt
molecule            Nc1_4N_05_pre ../../templates/NC1_4N/05_pre.txt
molecule            Nc1_4N_05_post ../../templates/NC1_4N/05_post.txt
molecule            Nc1_4N_06_pre ../../templates/NC1_4N/06_pre.txt
molecule            Nc1_4N_06_post ../../templates/NC1_4N/06_post.txt
molecule            Nc1_4N_07_pre ../../templates/NC1_4N/07_pre.txt
molecule            Nc1_4N_07_post ../../templates/NC1_4N/07_post.txt
molecule            Nc1_4N_08_pre ../../templates/NC1_4N/08_pre.txt
molecule            Nc1_4N_08_post ../../templates/NC1_4N/08_post.txt
molecule            Nc1_4N_09_pre ../../templates/NC1_4N/09_pre.txt
molecule            Nc1_4N_09_post ../../templates/NC1_4N/09_post.txt
molecule            Nc1_4N_10_pre ../../templates/NC1_4N/10_pre.txt
molecule            Nc1_4N_10_post ../../templates/NC1_4N/10_post.txt
molecule            Nc1_4N_11_pre ../../templates/NC1_4N/11_pre.txt
molecule            Nc1_4N_11_post ../../templates/NC1_4N/11_post.txt
molecule            Nc1_4N_12_pre ../../templates/NC1_4N/12_pre.txt
molecule            Nc1_4N_12_post ../../templates/NC1_4N/12_post.txt
molecule            Nc1_4N_13_pre ../../templates/NC1_4N/13_pre.txt
molecule            Nc1_4N_13_post ../../templates/NC1_4N/13_post.txt
molecule            Nc1_4N_14_pre ../../templates/NC1_4N/14_pre.txt
molecule            Nc1_4N_14_post ../../templates/NC1_4N/14_post.txt
molecule            Nc1_4N_15_pre ../../templates/NC1_4N/15_pre.txt
molecule            Nc1_4N_15_post ../../templates/NC1_4N/15_post.txt
molecule            Nc1_4N_16_pre ../../templates/NC1_4N/16_pre.txt
molecule            Nc1_4N_16_post ../../templates/NC1_4N/16_post.txt
molecule            Nc1_4N_17_pre ../../templates/NC1_4N/17_pre.txt
molecule            Nc1_4N_17_post ../../templates/NC1_4N/17_post.txt
#MOLECULE TEMPLATES: Nc1-5N
molecule            Nc1_5N_01_pre ../../templates/NC1_5N/01_pre.txt
molecule            Nc1_5N_01_post ../../templates/NC1_5N/01_post.txt
molecule            Nc1_5N_02_pre ../../templates/NC1_5N/02_pre.txt
molecule            Nc1_5N_02_post ../../templates/NC1_5N/02_post.txt
molecule            Nc1_5N_03_pre ../../templates/NC1_5N/03_pre.txt
molecule            Nc1_5N_03_post ../../templates/NC1_5N/03_post.txt
molecule            Nc1_5N_04_pre ../../templates/NC1_5N/04_pre.txt
molecule            Nc1_5N_04_post ../../templates/NC1_5N/04_post.txt
molecule            Nc1_5N_05_pre ../../templates/NC1_5N/05_pre.txt
molecule            Nc1_5N_05_post ../../templates/NC1_5N/05_post.txt
molecule            Nc1_5N_06_pre ../../templates/NC1_5N/06_pre.txt
molecule            Nc1_5N_06_post ../../templates/NC1_5N/06_post.txt
molecule            Nc1_5N_07_pre ../../templates/NC1_5N/07_pre.txt
molecule            Nc1_5N_07_post ../../templates/NC1_5N/07_post.txt
#MOLECULE TEMPLATES: Nc1-6N
molecule            Nc1_6N_01_pre ../../templates/NC1_6N/01_pre.txt
molecule            Nc1_6N_01_post ../../templates/NC1_6N/01_post.txt
molecule            Nc1_6N_02_pre ../../templates/NC1_6N/02_pre.txt
molecule            Nc1_6N_02_post ../../templates/NC1_6N/02_post.txt
molecule            Nc1_6N_03_pre ../../templates/NC1_6N/03_pre.txt
molecule            Nc1_6N_03_post ../../templates/NC1_6N/03_post.txt
molecule            Nc1_6N_04_pre ../../templates/NC1_6N/04_pre.txt
molecule            Nc1_6N_04_post ../../templates/NC1_6N/04_post.txt
molecule            Nc1_6N_05_pre ../../templates/NC1_6N/05_pre.txt
molecule            Nc1_6N_05_post ../../templates/NC1_6N/05_post.txt
molecule            Nc1_6N_06_pre ../../templates/NC1_6N/06_pre.txt
molecule            Nc1_6N_06_post ../../templates/NC1_6N/06_post.txt
molecule            Nc1_6N_07_pre ../../templates/NC1_6N/07_pre.txt
molecule            Nc1_6N_07_post ../../templates/NC1_6N/07_post.txt
molecule            Nc1_6N_08_pre ../../templates/NC1_6N/08_pre.txt
molecule            Nc1_6N_08_post ../../templates/NC1_6N/08_post.txt
molecule            Nc1_6N_09_pre ../../templates/NC1_6N/09_pre.txt
molecule            Nc1_6N_09_post ../../templates/NC1_6N/09_post.txt
molecule            Nc1_6N_10_pre ../../templates/NC1_6N/10_pre.txt
molecule            Nc1_6N_10_post ../../templates/NC1_6N/10_post.txt
molecule            Nc1_6N_11_pre ../../templates/NC1_6N/11_pre.txt
molecule            Nc1_6N_11_post ../../templates/NC1_6N/11_post.txt
molecule            Nc1_6N_12_pre ../../templates/NC1_6N/12_pre.txt
molecule            Nc1_6N_12_post ../../templates/NC1_6N/12_post.txt
molecule            Nc1_6N_13_pre ../../templates/NC1_6N/13_pre.txt
molecule            Nc1_6N_13_post ../../templates/NC1_6N/13_post.txt
molecule            Nc1_6N_14_pre ../../templates/NC1_6N/14_pre.txt
molecule            Nc1_6N_14_post ../../templates/NC1_6N/14_post.txt
molecule            Nc1_6N_15_pre ../../templates/NC1_6N/15_pre.txt
molecule            Nc1_6N_15_post ../../templates/NC1_6N/15_post.txt
#MOLECULE TEMPLATES: Nc1-7N
molecule            Nc1_7N_01_pre ../../templates/NC1_7N/01_pre.txt
molecule            Nc1_7N_01_post ../../templates/NC1_7N/01_post.txt
molecule            Nc1_7N_02_pre ../../templates/NC1_7N/02_pre.txt
molecule            Nc1_7N_02_post ../../templates/NC1_7N/02_post.txt
molecule            Nc1_7N_03_pre ../../templates/NC1_7N/03_pre.txt
molecule            Nc1_7N_03_post ../../templates/NC1_7N/03_post.txt
molecule            Nc1_7N_04_pre ../../templates/NC1_7N/04_pre.txt
molecule            Nc1_7N_04_post ../../templates/NC1_7N/04_post.txt
molecule            Nc1_7N_05_pre ../../templates/NC1_7N/05_pre.txt
molecule            Nc1_7N_05_post ../../templates/NC1_7N/05_post.txt
#MOLECULE TEMPLATES: Nc1-8N
molecule            Nc1_8N_01_pre ../../templates/NC1_8N/01_pre.txt
molecule            Nc1_8N_01_post ../../templates/NC1_8N/01_post.txt
#MOLECULE TEMPLATES: loop
molecule            Nc1_loop_01_pre ../../templates/NC1_loop/01_pre.txt
molecule            Nc1_loop_01_post ../../templates/NC1_loop/01_post.txt
molecule            Nc1_loop_02_pre ../../templates/NC1_loop/02_pre.txt
molecule            Nc1_loop_02_post ../../templates/NC1_loop/02_post.txt
molecule            Nc1_loop_03_pre ../../templates/NC1_loop/03_pre.txt
molecule            Nc1_loop_03_post ../../templates/NC1_loop/03_post.txt
molecule            Nc1_loop_04_pre ../../templates/NC1_loop/04_pre.txt
molecule            Nc1_loop_04_post ../../templates/NC1_loop/04_post.txt
molecule            Nc1_loop_05_pre ../../templates/NC1_loop/05_pre.txt
molecule            Nc1_loop_05_post ../../templates/NC1_loop/05_post.txt
molecule            Nc1_loop_06_pre ../../templates/NC1_loop/06_pre.txt
molecule            Nc1_loop_06_post ../../templates/NC1_loop/06_post.txt
molecule            Nc1_loop_07_pre ../../templates/NC1_loop/07_pre.txt
molecule            Nc1_loop_07_post ../../templates/NC1_loop/07_post.txt
molecule            Nc1_loop_08_pre ../../templates/NC1_loop/08_pre.txt
molecule            Nc1_loop_08_post ../../templates/NC1_loop/08_post.txt
molecule            Nc1_loop_09_pre ../../templates/NC1_loop/09_pre.txt
molecule            Nc1_loop_09_post ../../templates/NC1_loop/09_post.txt
molecule            Nc1_loop_10_pre ../../templates/NC1_loop/10_pre.txt
molecule            Nc1_loop_10_post ../../templates/NC1_loop/10_post.txt
molecule            Nc1_loop_11_pre ../../templates/NC1_loop/11_pre.txt
molecule            Nc1_loop_11_post ../../templates/NC1_loop/11_post.txt


#MOLECULE TEMPLATES: MonoToDimer
molecule            MonoToDimer_01_pre ../../templates/MonomerToDimer/01_pre.txt
molecule            MonoToDimer_01_post ../../templates/MonomerToDimer/01_post.txt
molecule            MonoToDimer_02_pre ../../templates/MonomerToDimer/02_pre.txt
molecule            MonoToDimer_02_post ../../templates/MonomerToDimer/02_post.txt
molecule            MonoToDimer_03_pre ../../templates/MonomerToDimer/03_pre.txt
molecule            MonoToDimer_03_post ../../templates/MonomerToDimer/03_post.txt
molecule            MonoToDimer_04_pre ../../templates/MonomerToDimer/04_pre.txt
molecule            MonoToDimer_04_post ../../templates/MonomerToDimer/04_post.txt
#MOLECULE TEMPLATES: DimerToLineTrimer
molecule            DimerToLineTrimer_01_pre ../../templates/DimerToLineTrimer/01_pre.txt
molecule            DimerToLineTrimer_01_post ../../templates/DimerToLineTrimer/01_post.txt
molecule            DimerToLineTrimer_02_pre ../../templates/DimerToLineTrimer/02_pre.txt
molecule            DimerToLineTrimer_02_post ../../templates/DimerToLineTrimer/02_post.txt
molecule            DimerToLineTrimer_03_pre ../../templates/DimerToLineTrimer/03_pre.txt
molecule            DimerToLineTrimer_03_post ../../templates/DimerToLineTrimer/03_post.txt
molecule            DimerToLineTrimer_04_pre ../../templates/DimerToLineTrimer/04_pre.txt
molecule            DimerToLineTrimer_04_post ../../templates/DimerToLineTrimer/04_post.txt
molecule            DimerToLineTrimer_05_pre ../../templates/DimerToLineTrimer/05_pre.txt
molecule            DimerToLineTrimer_05_post ../../templates/DimerToLineTrimer/05_post.txt
molecule            DimerToLineTrimer_06_pre ../../templates/DimerToLineTrimer/06_pre.txt
molecule            DimerToLineTrimer_06_post ../../templates/DimerToLineTrimer/06_post.txt
molecule            DimerToLineTrimer_07_pre ../../templates/DimerToLineTrimer/07_pre.txt
molecule            DimerToLineTrimer_07_post ../../templates/DimerToLineTrimer/07_post.txt
molecule            DimerToLineTrimer_08_pre ../../templates/DimerToLineTrimer/08_pre.txt
molecule            DimerToLineTrimer_08_post ../../templates/DimerToLineTrimer/08_post.txt
molecule            DimerToLineTrimer_09_pre ../../templates/DimerToLineTrimer/09_pre.txt
molecule            DimerToLineTrimer_09_post ../../templates/DimerToLineTrimer/09_post.txt
molecule            DimerToLineTrimer_10_pre ../../templates/DimerToLineTrimer/10_pre.txt
molecule            DimerToLineTrimer_10_post ../../templates/DimerToLineTrimer/10_post.txt
molecule            DimerToLineTrimer_11_pre ../../templates/DimerToLineTrimer/11_pre.txt
molecule            DimerToLineTrimer_11_post ../../templates/DimerToLineTrimer/11_post.txt
molecule            DimerToLineTrimer_12_pre ../../templates/DimerToLineTrimer/12_pre.txt
molecule            DimerToLineTrimer_12_post ../../templates/DimerToLineTrimer/12_post.txt
molecule            DimerToLineTrimer_13_pre ../../templates/DimerToLineTrimer/13_pre.txt
molecule            DimerToLineTrimer_13_post ../../templates/DimerToLineTrimer/13_post.txt
molecule            DimerToLineTrimer_14_pre ../../templates/DimerToLineTrimer/14_pre.txt
molecule            DimerToLineTrimer_14_post ../../templates/DimerToLineTrimer/14_post.txt
#MOLECULE TEMPLATES: LineTrimerToTriangleTrimer
molecule            LineTrimerToTriangleTrimer_01_pre ../../templates/LineTrimerToTriangleTrimer/01_pre.txt
molecule            LineTrimerToTriangleTrimer_01_post ../../templates/LineTrimerToTriangleTrimer/01_post.txt
molecule            LineTrimerToTriangleTrimer_02_pre ../../templates/LineTrimerToTriangleTrimer/02_pre.txt
molecule            LineTrimerToTriangleTrimer_02_post ../../templates/LineTrimerToTriangleTrimer/02_post.txt
molecule            LineTrimerToTriangleTrimer_03_pre ../../templates/LineTrimerToTriangleTrimer/03_pre.txt
molecule            LineTrimerToTriangleTrimer_03_post ../../templates/LineTrimerToTriangleTrimer/03_post.txt
molecule            LineTrimerToTriangleTrimer_04_pre ../../templates/LineTrimerToTriangleTrimer/04_pre.txt
molecule            LineTrimerToTriangleTrimer_04_post ../../templates/LineTrimerToTriangleTrimer/04_post.txt
molecule            LineTrimerToTriangleTrimer_05_pre ../../templates/LineTrimerToTriangleTrimer/05_pre.txt
molecule            LineTrimerToTriangleTrimer_05_post ../../templates/LineTrimerToTriangleTrimer/05_post.txt
molecule            LineTrimerToTriangleTrimer_06_pre ../../templates/LineTrimerToTriangleTrimer/06_pre.txt
molecule            LineTrimerToTriangleTrimer_06_post ../../templates/LineTrimerToTriangleTrimer/06_post.txt
molecule            LineTrimerToTriangleTrimer_07_pre ../../templates/LineTrimerToTriangleTrimer/07_pre.txt
molecule            LineTrimerToTriangleTrimer_07_post ../../templates/LineTrimerToTriangleTrimer/07_post.txt
molecule            LineTrimerToTriangleTrimer_08_pre ../../templates/LineTrimerToTriangleTrimer/08_pre.txt
molecule            LineTrimerToTriangleTrimer_08_post ../../templates/LineTrimerToTriangleTrimer/08_post.txt
molecule            LineTrimerToTriangleTrimer_09_pre ../../templates/LineTrimerToTriangleTrimer/09_pre.txt
molecule            LineTrimerToTriangleTrimer_09_post ../../templates/LineTrimerToTriangleTrimer/09_post.txt
molecule            LineTrimerToTriangleTrimer_10_pre ../../templates/LineTrimerToTriangleTrimer/10_pre.txt
molecule            LineTrimerToTriangleTrimer_10_post ../../templates/LineTrimerToTriangleTrimer/10_post.txt
#MOLECULE TEMPLATES: TriangleTrimerToOpenTetramer
molecule            TriangleTrimerToOpenTetramer_01_pre ../../templates/TriangleTrimerToOpenTetramer/01_pre.txt
molecule            TriangleTrimerToOpenTetramer_01_post ../../templates/TriangleTrimerToOpenTetramer/01_post.txt
molecule            TriangleTrimerToOpenTetramer_02_pre ../../templates/TriangleTrimerToOpenTetramer/02_pre.txt
molecule            TriangleTrimerToOpenTetramer_02_post ../../templates/TriangleTrimerToOpenTetramer/02_post.txt
molecule            TriangleTrimerToOpenTetramer_03_pre ../../templates/TriangleTrimerToOpenTetramer/03_pre.txt
molecule            TriangleTrimerToOpenTetramer_03_post ../../templates/TriangleTrimerToOpenTetramer/03_post.txt
molecule            TriangleTrimerToOpenTetramer_04_pre ../../templates/TriangleTrimerToOpenTetramer/04_pre.txt
molecule            TriangleTrimerToOpenTetramer_04_post ../../templates/TriangleTrimerToOpenTetramer/04_post.txt
molecule            TriangleTrimerToOpenTetramer_05_pre ../../templates/TriangleTrimerToOpenTetramer/05_pre.txt
molecule            TriangleTrimerToOpenTetramer_05_post ../../templates/TriangleTrimerToOpenTetramer/05_post.txt
molecule            TriangleTrimerToOpenTetramer_06_pre ../../templates/TriangleTrimerToOpenTetramer/06_pre.txt
molecule            TriangleTrimerToOpenTetramer_06_post ../../templates/TriangleTrimerToOpenTetramer/06_post.txt
molecule            TriangleTrimerToOpenTetramer_07_pre ../../templates/TriangleTrimerToOpenTetramer/07_pre.txt
molecule            TriangleTrimerToOpenTetramer_07_post ../../templates/TriangleTrimerToOpenTetramer/07_post.txt
molecule            TriangleTrimerToOpenTetramer_08_pre ../../templates/TriangleTrimerToOpenTetramer/08_pre.txt
molecule            TriangleTrimerToOpenTetramer_08_post ../../templates/TriangleTrimerToOpenTetramer/08_post.txt
molecule            TriangleTrimerToOpenTetramer_09_pre ../../templates/TriangleTrimerToOpenTetramer/09_pre.txt
molecule            TriangleTrimerToOpenTetramer_09_post ../../templates/TriangleTrimerToOpenTetramer/09_post.txt
molecule            TriangleTrimerToOpenTetramer_10_pre ../../templates/TriangleTrimerToOpenTetramer/10_pre.txt
molecule            TriangleTrimerToOpenTetramer_10_post ../../templates/TriangleTrimerToOpenTetramer/10_post.txt
molecule            TriangleTrimerToOpenTetramer_11_pre ../../templates/TriangleTrimerToOpenTetramer/11_pre.txt
molecule            TriangleTrimerToOpenTetramer_11_post ../../templates/TriangleTrimerToOpenTetramer/11_post.txt
molecule            TriangleTrimerToOpenTetramer_12_pre ../../templates/TriangleTrimerToOpenTetramer/12_pre.txt
molecule            TriangleTrimerToOpenTetramer_12_post ../../templates/TriangleTrimerToOpenTetramer/12_post.txt
molecule            TriangleTrimerToOpenTetramer_13_pre ../../templates/TriangleTrimerToOpenTetramer/13_pre.txt
molecule            TriangleTrimerToOpenTetramer_13_post ../../templates/TriangleTrimerToOpenTetramer/13_post.txt
molecule            TriangleTrimerToOpenTetramer_14_pre ../../templates/TriangleTrimerToOpenTetramer/14_pre.txt
molecule            TriangleTrimerToOpenTetramer_14_post ../../templates/TriangleTrimerToOpenTetramer/14_post.txt
molecule            TriangleTrimerToOpenTetramer_15_pre ../../templates/TriangleTrimerToOpenTetramer/15_pre.txt
molecule            TriangleTrimerToOpenTetramer_15_post ../../templates/TriangleTrimerToOpenTetramer/15_post.txt
molecule            TriangleTrimerToOpenTetramer_16_pre ../../templates/TriangleTrimerToOpenTetramer/16_pre.txt
molecule            TriangleTrimerToOpenTetramer_16_post ../../templates/TriangleTrimerToOpenTetramer/16_post.txt
molecule            TriangleTrimerToOpenTetramer_17_pre ../../templates/TriangleTrimerToOpenTetramer/17_pre.txt
molecule            TriangleTrimerToOpenTetramer_17_post ../../templates/TriangleTrimerToOpenTetramer/17_post.txt
molecule            TriangleTrimerToOpenTetramer_18_pre ../../templates/TriangleTrimerToOpenTetramer/18_pre.txt
molecule            TriangleTrimerToOpenTetramer_18_post ../../templates/TriangleTrimerToOpenTetramer/18_post.txt
molecule            TriangleTrimerToOpenTetramer_19_pre ../../templates/TriangleTrimerToOpenTetramer/19_pre.txt
molecule            TriangleTrimerToOpenTetramer_19_post ../../templates/TriangleTrimerToOpenTetramer/19_post.txt
molecule            TriangleTrimerToOpenTetramer_20_pre ../../templates/TriangleTrimerToOpenTetramer/20_pre.txt
molecule            TriangleTrimerToOpenTetramer_20_post ../../templates/TriangleTrimerToOpenTetramer/20_post.txt
molecule            TriangleTrimerToOpenTetramer_21_pre ../../templates/TriangleTrimerToOpenTetramer/21_pre.txt
molecule            TriangleTrimerToOpenTetramer_21_post ../../templates/TriangleTrimerToOpenTetramer/21_post.txt
molecule            TriangleTrimerToOpenTetramer_22_pre ../../templates/TriangleTrimerToOpenTetramer/22_pre.txt
molecule            TriangleTrimerToOpenTetramer_22_post ../../templates/TriangleTrimerToOpenTetramer/22_post.txt
molecule            TriangleTrimerToOpenTetramer_23_pre ../../templates/TriangleTrimerToOpenTetramer/23_pre.txt
molecule            TriangleTrimerToOpenTetramer_23_post ../../templates/TriangleTrimerToOpenTetramer/23_post.txt
molecule            TriangleTrimerToOpenTetramer_24_pre ../../templates/TriangleTrimerToOpenTetramer/24_pre.txt
molecule            TriangleTrimerToOpenTetramer_24_post ../../templates/TriangleTrimerToOpenTetramer/24_post.txt
molecule            TriangleTrimerToOpenTetramer_25_pre ../../templates/TriangleTrimerToOpenTetramer/25_pre.txt
molecule            TriangleTrimerToOpenTetramer_25_post ../../templates/TriangleTrimerToOpenTetramer/25_post.txt
molecule            TriangleTrimerToOpenTetramer_26_pre ../../templates/TriangleTrimerToOpenTetramer/26_pre.txt
molecule            TriangleTrimerToOpenTetramer_26_post ../../templates/TriangleTrimerToOpenTetramer/26_post.txt
molecule            TriangleTrimerToOpenTetramer_27_pre ../../templates/TriangleTrimerToOpenTetramer/27_pre.txt
molecule            TriangleTrimerToOpenTetramer_27_post ../../templates/TriangleTrimerToOpenTetramer/27_post.txt
molecule            TriangleTrimerToOpenTetramer_28_pre ../../templates/TriangleTrimerToOpenTetramer/28_pre.txt
molecule            TriangleTrimerToOpenTetramer_28_post ../../templates/TriangleTrimerToOpenTetramer/28_post.txt
molecule            TriangleTrimerToOpenTetramer_29_pre ../../templates/TriangleTrimerToOpenTetramer/29_pre.txt
molecule            TriangleTrimerToOpenTetramer_29_post ../../templates/TriangleTrimerToOpenTetramer/29_post.txt
#MOLECULE TEMPLATES: OpenTetramerToRhomboidTetramer
molecule            OpenTetramerToRhomboidTetramer_01_pre ../../templates/TriangleTrimerToOpenTetramer/01_post.txt
molecule            OpenTetramerToRhomboidTetramer_01_post ../../templates/OpenTetramerToRhomboidTetramer/01_post.txt
molecule            OpenTetramerToRhomboidTetramer_02_pre ../../templates/TriangleTrimerToOpenTetramer/02_post.txt
molecule            OpenTetramerToRhomboidTetramer_02_post ../../templates/OpenTetramerToRhomboidTetramer/02_post.txt
molecule            OpenTetramerToRhomboidTetramer_03_pre ../../templates/TriangleTrimerToOpenTetramer/03_post.txt
molecule            OpenTetramerToRhomboidTetramer_03_post ../../templates/OpenTetramerToRhomboidTetramer/03_post.txt
molecule            OpenTetramerToRhomboidTetramer_04_pre ../../templates/TriangleTrimerToOpenTetramer/04_post.txt
molecule            OpenTetramerToRhomboidTetramer_04_post ../../templates/OpenTetramerToRhomboidTetramer/04_post.txt
molecule            OpenTetramerToRhomboidTetramer_05_pre ../../templates/TriangleTrimerToOpenTetramer/05_post.txt
molecule            OpenTetramerToRhomboidTetramer_05_post ../../templates/OpenTetramerToRhomboidTetramer/05_post.txt
molecule            OpenTetramerToRhomboidTetramer_06_pre ../../templates/TriangleTrimerToOpenTetramer/06_post.txt
molecule            OpenTetramerToRhomboidTetramer_06_post ../../templates/OpenTetramerToRhomboidTetramer/06_post.txt
molecule            OpenTetramerToRhomboidTetramer_07_pre ../../templates/TriangleTrimerToOpenTetramer/07_post.txt
molecule            OpenTetramerToRhomboidTetramer_07_post ../../templates/OpenTetramerToRhomboidTetramer/07_post.txt
molecule            OpenTetramerToRhomboidTetramer_08_pre ../../templates/TriangleTrimerToOpenTetramer/08_post.txt
molecule            OpenTetramerToRhomboidTetramer_08_post ../../templates/OpenTetramerToRhomboidTetramer/08_post.txt
molecule            OpenTetramerToRhomboidTetramer_09_pre ../../templates/TriangleTrimerToOpenTetramer/09_post.txt
molecule            OpenTetramerToRhomboidTetramer_09_post ../../templates/OpenTetramerToRhomboidTetramer/09_post.txt
molecule            OpenTetramerToRhomboidTetramer_10_pre ../../templates/TriangleTrimerToOpenTetramer/10_post.txt
molecule            OpenTetramerToRhomboidTetramer_10_post ../../templates/OpenTetramerToRhomboidTetramer/10_post.txt
molecule            OpenTetramerToRhomboidTetramer_11_pre ../../templates/TriangleTrimerToOpenTetramer/11_post.txt
molecule            OpenTetramerToRhomboidTetramer_11_post ../../templates/OpenTetramerToRhomboidTetramer/11_post.txt
molecule            OpenTetramerToRhomboidTetramer_12_pre ../../templates/TriangleTrimerToOpenTetramer/12_post.txt
molecule            OpenTetramerToRhomboidTetramer_12_post ../../templates/OpenTetramerToRhomboidTetramer/12_post.txt
molecule            OpenTetramerToRhomboidTetramer_13_pre ../../templates/TriangleTrimerToOpenTetramer/13_post.txt
molecule            OpenTetramerToRhomboidTetramer_13_post ../../templates/OpenTetramerToRhomboidTetramer/13_post.txt
molecule            OpenTetramerToRhomboidTetramer_14_pre ../../templates/TriangleTrimerToOpenTetramer/14_post.txt
molecule            OpenTetramerToRhomboidTetramer_14_post ../../templates/OpenTetramerToRhomboidTetramer/14_post.txt
molecule            OpenTetramerToRhomboidTetramer_15_pre ../../templates/TriangleTrimerToOpenTetramer/15_post.txt
molecule            OpenTetramerToRhomboidTetramer_15_post ../../templates/OpenTetramerToRhomboidTetramer/15_post.txt
molecule            OpenTetramerToRhomboidTetramer_16_pre ../../templates/TriangleTrimerToOpenTetramer/16_post.txt
molecule            OpenTetramerToRhomboidTetramer_16_post ../../templates/OpenTetramerToRhomboidTetramer/16_post.txt
molecule            OpenTetramerToRhomboidTetramer_17_pre ../../templates/TriangleTrimerToOpenTetramer/17_post.txt
molecule            OpenTetramerToRhomboidTetramer_17_post ../../templates/OpenTetramerToRhomboidTetramer/17_post.txt
molecule            OpenTetramerToRhomboidTetramer_18_pre ../../templates/TriangleTrimerToOpenTetramer/18_post.txt
molecule            OpenTetramerToRhomboidTetramer_18_post ../../templates/OpenTetramerToRhomboidTetramer/18_post.txt
molecule            OpenTetramerToRhomboidTetramer_19_pre ../../templates/TriangleTrimerToOpenTetramer/19_post.txt
molecule            OpenTetramerToRhomboidTetramer_19_post ../../templates/OpenTetramerToRhomboidTetramer/19_post.txt
molecule            OpenTetramerToRhomboidTetramer_20_pre ../../templates/TriangleTrimerToOpenTetramer/20_post.txt
molecule            OpenTetramerToRhomboidTetramer_20_post ../../templates/OpenTetramerToRhomboidTetramer/20_post.txt
molecule            OpenTetramerToRhomboidTetramer_21_pre ../../templates/TriangleTrimerToOpenTetramer/21_post.txt
molecule            OpenTetramerToRhomboidTetramer_21_post ../../templates/OpenTetramerToRhomboidTetramer/21_post.txt
molecule            OpenTetramerToRhomboidTetramer_22_pre ../../templates/TriangleTrimerToOpenTetramer/22_post.txt
molecule            OpenTetramerToRhomboidTetramer_22_post ../../templates/OpenTetramerToRhomboidTetramer/22_post.txt
molecule            OpenTetramerToRhomboidTetramer_23_pre ../../templates/TriangleTrimerToOpenTetramer/23_post.txt
molecule            OpenTetramerToRhomboidTetramer_23_post ../../templates/OpenTetramerToRhomboidTetramer/23_post.txt
molecule            OpenTetramerToRhomboidTetramer_24_pre ../../templates/TriangleTrimerToOpenTetramer/24_post.txt
molecule            OpenTetramerToRhomboidTetramer_24_post ../../templates/OpenTetramerToRhomboidTetramer/24_post.txt
molecule            OpenTetramerToRhomboidTetramer_25_pre ../../templates/TriangleTrimerToOpenTetramer/25_post.txt
molecule            OpenTetramerToRhomboidTetramer_25_post ../../templates/OpenTetramerToRhomboidTetramer/25_post.txt
molecule            OpenTetramerToRhomboidTetramer_26_pre ../../templates/TriangleTrimerToOpenTetramer/26_post.txt
molecule            OpenTetramerToRhomboidTetramer_26_post ../../templates/OpenTetramerToRhomboidTetramer/26_post.txt
molecule            OpenTetramerToRhomboidTetramer_27_pre ../../templates/TriangleTrimerToOpenTetramer/27_post.txt
molecule            OpenTetramerToRhomboidTetramer_27_post ../../templates/OpenTetramerToRhomboidTetramer/27_post.txt
molecule            OpenTetramerToRhomboidTetramer_28_pre ../../templates/TriangleTrimerToOpenTetramer/28_post.txt
molecule            OpenTetramerToRhomboidTetramer_28_post ../../templates/OpenTetramerToRhomboidTetramer/28_post.txt
molecule            OpenTetramerToRhomboidTetramer_29_pre ../../templates/TriangleTrimerToOpenTetramer/29_post.txt
molecule            OpenTetramerToRhomboidTetramer_29_post ../../templates/OpenTetramerToRhomboidTetramer/29_post.txt
molecule            OpenTetramerToRhomboidTetramer_30_pre ../../templates/TriangleTrimerToOpenTetramer/03_post.txt
molecule            OpenTetramerToRhomboidTetramer_30_post ../../templates/OpenTetramerToRhomboidTetramer/30_post.txt
molecule            OpenTetramerToRhomboidTetramer_31_pre ../../templates/TriangleTrimerToOpenTetramer/05_post.txt
molecule            OpenTetramerToRhomboidTetramer_31_post ../../templates/OpenTetramerToRhomboidTetramer/31_post.txt
molecule            OpenTetramerToRhomboidTetramer_32_pre ../../templates/TriangleTrimerToOpenTetramer/08_post.txt
molecule            OpenTetramerToRhomboidTetramer_32_post ../../templates/OpenTetramerToRhomboidTetramer/32_post.txt
molecule            OpenTetramerToRhomboidTetramer_33_pre ../../templates/TriangleTrimerToOpenTetramer/09_post.txt
molecule            OpenTetramerToRhomboidTetramer_33_post ../../templates/OpenTetramerToRhomboidTetramer/33_post.txt
molecule            OpenTetramerToRhomboidTetramer_34_pre ../../templates/TriangleTrimerToOpenTetramer/12_post.txt
molecule            OpenTetramerToRhomboidTetramer_34_post ../../templates/OpenTetramerToRhomboidTetramer/34_post.txt
molecule            OpenTetramerToRhomboidTetramer_35_pre ../../templates/TriangleTrimerToOpenTetramer/14_post.txt
molecule            OpenTetramerToRhomboidTetramer_35_post ../../templates/OpenTetramerToRhomboidTetramer/35_post.txt
molecule            OpenTetramerToRhomboidTetramer_36_pre ../../templates/TriangleTrimerToOpenTetramer/16_post.txt
molecule            OpenTetramerToRhomboidTetramer_36_post ../../templates/OpenTetramerToRhomboidTetramer/36_post.txt
molecule            OpenTetramerToRhomboidTetramer_37_pre ../../templates/TriangleTrimerToOpenTetramer/17_post.txt
molecule            OpenTetramerToRhomboidTetramer_37_post ../../templates/OpenTetramerToRhomboidTetramer/37_post.txt
molecule            OpenTetramerToRhomboidTetramer_38_pre ../../templates/TriangleTrimerToOpenTetramer/18_post.txt
molecule            OpenTetramerToRhomboidTetramer_38_post ../../templates/OpenTetramerToRhomboidTetramer/38_post.txt
molecule            OpenTetramerToRhomboidTetramer_39_pre ../../templates/TriangleTrimerToOpenTetramer/21_post.txt
molecule            OpenTetramerToRhomboidTetramer_39_post ../../templates/OpenTetramerToRhomboidTetramer/39_post.txt
molecule            OpenTetramerToRhomboidTetramer_40_pre ../../templates/TriangleTrimerToOpenTetramer/22_post.txt
molecule            OpenTetramerToRhomboidTetramer_40_post ../../templates/OpenTetramerToRhomboidTetramer/40_post.txt
molecule            OpenTetramerToRhomboidTetramer_41_pre ../../templates/TriangleTrimerToOpenTetramer/24_post.txt
molecule            OpenTetramerToRhomboidTetramer_41_post ../../templates/OpenTetramerToRhomboidTetramer/41_post.txt
#MOLECULE TEMPLATES: RhomboidTetramerToSquareTetramer
molecule            RhomboidTetramerToSquareTetramer_01_pre ../../templates/OpenTetramerToRhomboidTetramer/01_post.txt
molecule            RhomboidTetramerToSquareTetramer_01_post ../../templates/RhomboidTetramerToSquareTetramer/01_post.txt
molecule            RhomboidTetramerToSquareTetramer_02_pre ../../templates/OpenTetramerToRhomboidTetramer/02_post.txt
molecule            RhomboidTetramerToSquareTetramer_02_post ../../templates/RhomboidTetramerToSquareTetramer/02_post.txt
molecule            RhomboidTetramerToSquareTetramer_03_pre ../../templates/OpenTetramerToRhomboidTetramer/03_post.txt
molecule            RhomboidTetramerToSquareTetramer_03_post ../../templates/RhomboidTetramerToSquareTetramer/03_post.txt
molecule            RhomboidTetramerToSquareTetramer_04_pre ../../templates/OpenTetramerToRhomboidTetramer/04_post.txt
molecule            RhomboidTetramerToSquareTetramer_04_post ../../templates/RhomboidTetramerToSquareTetramer/04_post.txt
molecule            RhomboidTetramerToSquareTetramer_05_pre ../../templates/OpenTetramerToRhomboidTetramer/05_post.txt
molecule            RhomboidTetramerToSquareTetramer_05_post ../../templates/RhomboidTetramerToSquareTetramer/05_post.txt
molecule            RhomboidTetramerToSquareTetramer_06_pre ../../templates/OpenTetramerToRhomboidTetramer/06_post.txt
molecule            RhomboidTetramerToSquareTetramer_06_post ../../templates/RhomboidTetramerToSquareTetramer/06_post.txt
molecule            RhomboidTetramerToSquareTetramer_07_pre ../../templates/OpenTetramerToRhomboidTetramer/07_post.txt
molecule            RhomboidTetramerToSquareTetramer_07_post ../../templates/RhomboidTetramerToSquareTetramer/07_post.txt
molecule            RhomboidTetramerToSquareTetramer_08_pre ../../templates/OpenTetramerToRhomboidTetramer/08_post.txt
molecule            RhomboidTetramerToSquareTetramer_08_post ../../templates/RhomboidTetramerToSquareTetramer/08_post.txt
molecule            RhomboidTetramerToSquareTetramer_09_pre ../../templates/OpenTetramerToRhomboidTetramer/09_post.txt
molecule            RhomboidTetramerToSquareTetramer_09_post ../../templates/RhomboidTetramerToSquareTetramer/09_post.txt
molecule            RhomboidTetramerToSquareTetramer_10_pre ../../templates/OpenTetramerToRhomboidTetramer/10_post.txt
molecule            RhomboidTetramerToSquareTetramer_10_post ../../templates/RhomboidTetramerToSquareTetramer/10_post.txt
molecule            RhomboidTetramerToSquareTetramer_11_pre ../../templates/OpenTetramerToRhomboidTetramer/11_post.txt
molecule            RhomboidTetramerToSquareTetramer_11_post ../../templates/RhomboidTetramerToSquareTetramer/11_post.txt
molecule            RhomboidTetramerToSquareTetramer_12_pre ../../templates/OpenTetramerToRhomboidTetramer/12_post.txt
molecule            RhomboidTetramerToSquareTetramer_12_post ../../templates/RhomboidTetramerToSquareTetramer/12_post.txt
molecule            RhomboidTetramerToSquareTetramer_13_pre ../../templates/OpenTetramerToRhomboidTetramer/13_post.txt
molecule            RhomboidTetramerToSquareTetramer_13_post ../../templates/RhomboidTetramerToSquareTetramer/13_post.txt
molecule            RhomboidTetramerToSquareTetramer_14_pre ../../templates/OpenTetramerToRhomboidTetramer/14_post.txt
molecule            RhomboidTetramerToSquareTetramer_14_post ../../templates/RhomboidTetramerToSquareTetramer/14_post.txt
molecule            RhomboidTetramerToSquareTetramer_15_pre ../../templates/OpenTetramerToRhomboidTetramer/15_post.txt
molecule            RhomboidTetramerToSquareTetramer_15_post ../../templates/RhomboidTetramerToSquareTetramer/15_post.txt
molecule            RhomboidTetramerToSquareTetramer_16_pre ../../templates/OpenTetramerToRhomboidTetramer/16_post.txt
molecule            RhomboidTetramerToSquareTetramer_16_post ../../templates/RhomboidTetramerToSquareTetramer/16_post.txt
molecule            RhomboidTetramerToSquareTetramer_17_pre ../../templates/OpenTetramerToRhomboidTetramer/17_post.txt
molecule            RhomboidTetramerToSquareTetramer_17_post ../../templates/RhomboidTetramerToSquareTetramer/17_post.txt
molecule            RhomboidTetramerToSquareTetramer_18_pre ../../templates/OpenTetramerToRhomboidTetramer/18_post.txt
molecule            RhomboidTetramerToSquareTetramer_18_post ../../templates/RhomboidTetramerToSquareTetramer/18_post.txt
molecule            RhomboidTetramerToSquareTetramer_19_pre ../../templates/OpenTetramerToRhomboidTetramer/19_post.txt
molecule            RhomboidTetramerToSquareTetramer_19_post ../../templates/RhomboidTetramerToSquareTetramer/19_post.txt
molecule            RhomboidTetramerToSquareTetramer_20_pre ../../templates/OpenTetramerToRhomboidTetramer/20_post.txt
molecule            RhomboidTetramerToSquareTetramer_20_post ../../templates/RhomboidTetramerToSquareTetramer/20_post.txt
molecule            RhomboidTetramerToSquareTetramer_21_pre ../../templates/OpenTetramerToRhomboidTetramer/21_post.txt
molecule            RhomboidTetramerToSquareTetramer_21_post ../../templates/RhomboidTetramerToSquareTetramer/21_post.txt
molecule            RhomboidTetramerToSquareTetramer_22_pre ../../templates/OpenTetramerToRhomboidTetramer/22_post.txt
molecule            RhomboidTetramerToSquareTetramer_22_post ../../templates/RhomboidTetramerToSquareTetramer/22_post.txt
molecule            RhomboidTetramerToSquareTetramer_23_pre ../../templates/OpenTetramerToRhomboidTetramer/23_post.txt
molecule            RhomboidTetramerToSquareTetramer_23_post ../../templates/RhomboidTetramerToSquareTetramer/23_post.txt
molecule            RhomboidTetramerToSquareTetramer_24_pre ../../templates/OpenTetramerToRhomboidTetramer/24_post.txt
molecule            RhomboidTetramerToSquareTetramer_24_post ../../templates/RhomboidTetramerToSquareTetramer/24_post.txt
molecule            RhomboidTetramerToSquareTetramer_25_pre ../../templates/OpenTetramerToRhomboidTetramer/25_post.txt
molecule            RhomboidTetramerToSquareTetramer_25_post ../../templates/RhomboidTetramerToSquareTetramer/25_post.txt
molecule            RhomboidTetramerToSquareTetramer_26_pre ../../templates/OpenTetramerToRhomboidTetramer/26_post.txt
molecule            RhomboidTetramerToSquareTetramer_26_post ../../templates/RhomboidTetramerToSquareTetramer/26_post.txt
molecule            RhomboidTetramerToSquareTetramer_27_pre ../../templates/OpenTetramerToRhomboidTetramer/27_post.txt
molecule            RhomboidTetramerToSquareTetramer_27_post ../../templates/RhomboidTetramerToSquareTetramer/27_post.txt
molecule            RhomboidTetramerToSquareTetramer_28_pre ../../templates/OpenTetramerToRhomboidTetramer/28_post.txt
molecule            RhomboidTetramerToSquareTetramer_28_post ../../templates/RhomboidTetramerToSquareTetramer/28_post.txt
molecule            RhomboidTetramerToSquareTetramer_29_pre ../../templates/OpenTetramerToRhomboidTetramer/29_post.txt
molecule            RhomboidTetramerToSquareTetramer_29_post ../../templates/RhomboidTetramerToSquareTetramer/29_post.txt
molecule            RhomboidTetramerToSquareTetramer_30_pre ../../templates/OpenTetramerToRhomboidTetramer/30_post.txt
molecule            RhomboidTetramerToSquareTetramer_30_post ../../templates/RhomboidTetramerToSquareTetramer/30_post.txt
molecule            RhomboidTetramerToSquareTetramer_31_pre ../../templates/OpenTetramerToRhomboidTetramer/31_post.txt
molecule            RhomboidTetramerToSquareTetramer_31_post ../../templates/RhomboidTetramerToSquareTetramer/31_post.txt
molecule            RhomboidTetramerToSquareTetramer_32_pre ../../templates/OpenTetramerToRhomboidTetramer/32_post.txt
molecule            RhomboidTetramerToSquareTetramer_32_post ../../templates/RhomboidTetramerToSquareTetramer/32_post.txt
molecule            RhomboidTetramerToSquareTetramer_33_pre ../../templates/OpenTetramerToRhomboidTetramer/33_post.txt
molecule            RhomboidTetramerToSquareTetramer_33_post ../../templates/RhomboidTetramerToSquareTetramer/33_post.txt
molecule            RhomboidTetramerToSquareTetramer_34_pre ../../templates/OpenTetramerToRhomboidTetramer/34_post.txt
molecule            RhomboidTetramerToSquareTetramer_34_post ../../templates/RhomboidTetramerToSquareTetramer/34_post.txt
molecule            RhomboidTetramerToSquareTetramer_35_pre ../../templates/OpenTetramerToRhomboidTetramer/35_post.txt
molecule            RhomboidTetramerToSquareTetramer_35_post ../../templates/RhomboidTetramerToSquareTetramer/35_post.txt
molecule            RhomboidTetramerToSquareTetramer_36_pre ../../templates/OpenTetramerToRhomboidTetramer/36_post.txt
molecule            RhomboidTetramerToSquareTetramer_36_post ../../templates/RhomboidTetramerToSquareTetramer/36_post.txt
molecule            RhomboidTetramerToSquareTetramer_37_pre ../../templates/OpenTetramerToRhomboidTetramer/37_post.txt
molecule            RhomboidTetramerToSquareTetramer_37_post ../../templates/RhomboidTetramerToSquareTetramer/37_post.txt
molecule            RhomboidTetramerToSquareTetramer_38_pre ../../templates/OpenTetramerToRhomboidTetramer/38_post.txt
molecule            RhomboidTetramerToSquareTetramer_38_post ../../templates/RhomboidTetramerToSquareTetramer/38_post.txt
molecule            RhomboidTetramerToSquareTetramer_39_pre ../../templates/OpenTetramerToRhomboidTetramer/39_post.txt
molecule            RhomboidTetramerToSquareTetramer_39_post ../../templates/RhomboidTetramerToSquareTetramer/39_post.txt
molecule            RhomboidTetramerToSquareTetramer_40_pre ../../templates/OpenTetramerToRhomboidTetramer/40_post.txt
molecule            RhomboidTetramerToSquareTetramer_40_post ../../templates/RhomboidTetramerToSquareTetramer/40_post.txt
molecule            RhomboidTetramerToSquareTetramer_41_pre ../../templates/OpenTetramerToRhomboidTetramer/41_post.txt
molecule            RhomboidTetramerToSquareTetramer_41_post ../../templates/RhomboidTetramerToSquareTetramer/41_post.txt

##NUCLEATION TEMPLATES##

molecule nucleation_pre ../../templates/Nucleation/pre_Nucleation.txt
molecule nucleation_post ../../templates/Nucleation/post_Nucleation.txt
molecule death_pre ../../templates/Nucleation/pre_Death.txt
molecule death_post ../../templates/Nucleation/post_Death.txt


#COMPUTE BOND PROPERTIES
compute 		1 all property/local batom1 batom2 btype
compute 		4 all property/local atype aatom1 aatom2 aatom3
compute 		2 all pe bond
compute 		3 all pe angle
compute perAtomStress mobile stress/atom NULL
compute totalStressX mobile reduce sum c_perAtomStress[1]
compute totalStressY mobile reduce sum c_perAtomStress[2]
compute totalStressZ mobile reduce sum c_perAtomStress[3]
compute forceBond all bond/local force
compute forceBondAll all reduce ave c_forceBond

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

#fix 		fix_wall_lo 	all wall/lj126 zlo EDGE 1.0 1.0 ${lj_cutoff} #Lower wall
#fix 		fix_wall_hi 	all wall/lj126 zhi EDGE 1.0 1.0 ${lj_cutoff} #Upper wall

# Mask: mobile if type is 1..9
variable vMobileMask atom "type>=1 && type<=9"

# Dynamic group updated every step (or every N steps)
group mobile_dyn dynamic all var vMobileMask every 1

# Thermostat + integrator applied to dynamic group
fix fix_langevin mobile_dyn langevin ${mytemp} ${mytemp} ${damp} ${lseed}
fix fNVE         mobile_dyn nve



#CREATE VELOCITIES USING MAXWELL-BOLTZMANN DISTRIBUTION
velocity 	mobile_dyn create ${mytemp} ${vseed} dist gaussian mom yes rot yes 
run 0
velocity	mobile_dyn scale ${mytemp}

variable box_lengthX equal xhi-xlo
variable box_lengthY equal yhi-ylo
variable box_lengthZ equal zhi-zlo
variable box_Vol equal v_box_lengthX*v_box_lengthY*v_box_lengthZ
variable conc_preferred equal v_N_preferred/v_volInit #This is specific to this starting box volume, Change for different intitial data volumes
variable conc_pref_tot equal 1404/v_volInit 

variable            vMaskType9 atom "type==9"
group               type9 dynamic all var vMaskType9 every 1
variable            NumType9 equal count(type9)

variable            vMaskType1 atom "type==1"
group               type1 dynamic all var vMaskType1 every 1
variable            NumType1 equal count(type1)
variable            NumTot equal v_NumType1+v_NumType9
variable            NumTotPref equal v_conc_pref_tot*v_box_Vol

variable CurrentConc equal v_NumType9/v_box_Vol
variable CurrentConc_a equal (v_NumType9+1)/(v_box_Vol)

variable TermId_add equal (v_conc_preferred/(v_CurrentConc_a*2)) #divide by two because it recognises twice the initialisation bonded molecules (two configs 1-2 2-1)
variable NumType9Boot equal abs(0.5*(v_NumType9+1.0-abs(v_NumType9-1.0))-1.0)+v_NumType9 # weird way to avoid v_NumType9  = 0
variable TermId_remov equal (v_CurrentConc/(v_conc_preferred*v_NumType9Boot))
variable ConstRate equal 1/100

# variable En_o_add equal (v_KappaD/v_box_Vol)*(v_NumTot-v_conc_pref_tot*v_box_Vol)*(v_NumTot-v_conc_pref_tot*v_box_Vol)
# variable En_n_add equal (v_KappaD/v_box_Vol)*((v_NumTot+1)-v_conc_pref_tot*v_box_Vol)*((v_NumTot+1)-v_conc_pref_tot*v_box_Vol)

variable En_o_add equal (v_KappaD/1)*(v_NumTot-v_conc_pref_tot*v_box_Vol)*(v_NumTot-v_conc_pref_tot*v_box_Vol)
variable En_n_add equal (v_KappaD/1)*((v_NumTot+1)-v_conc_pref_tot*v_box_Vol)*((v_NumTot+1)-v_conc_pref_tot*v_box_Vol)
variable TermDen_add equal exp(-(v_En_n_add-v_En_o_add)) #Acc = np.min([1,C*t1*np.exp(-(En_n-En_o))])

# variable En_o_remov equal (v_KappaD/v_box_Vol)*(v_NumTot-v_conc_pref_tot*v_box_Vol)*(v_NumTot-v_conc_pref_tot*v_box_Vol)
# variable En_n_remov equal (v_KappaD/v_box_Vol)*((v_NumTot-1)-v_conc_pref_tot*v_box_Vol)*((v_NumTot-1)-v_conc_pref_tot*v_box_Vol)

variable En_o_remov equal (v_KappaD/1)*(v_NumTot-v_conc_pref_tot*v_box_Vol)*(v_NumTot-v_conc_pref_tot*v_box_Vol)
variable En_n_remov equal (v_KappaD/1)*((v_NumTot-1)-v_conc_pref_tot*v_box_Vol)*((v_NumTot-1)-v_conc_pref_tot*v_box_Vol)
variable TermDen_remov equal exp(-(v_En_n_remov-v_En_o_remov)) #Acc = np.min([1,C*t1*np.exp(-(En_n-En_o))])

variable P_add equal v_ConstRate*v_TermDen_add*v_TermId_add
variable P_remov equal v_ConstRate*v_TermDen_remov*v_TermId_remov

# Per-atom variable that just uses the scalar each step
variable P_add_atom atom v_P_add
variable P_remov_atom atom v_P_remov


fix		fix_print all print ${thermodump} "${step} ${etot} ${ke} ${peBond} ${peAngle} ${temp} ${press} ${StressX} ${StressY} ${StressZ} ${Fbond} ${NumType9} ${box_lengthX} ${box_Vol} ${CurrentConc} ${conc_preferred}" file ${thermofile} screen no title "step etot ke peBond peAngle temp press stressX stressY stressZ AvBondForce NumMons BoxLength BoxVol CurrentConc ConcPref"

#Save configurations in linear time
dump		2 all custom ${trajdump} dumplin/dump.${simname}.*.lammpstrj type id xu yu zu
dump_modify	2 sort id
dump_modify	2 pad 10
dump_modify	2 format line "%d %d %.8f %.8f %.8f"

#Save bond information
dump 		3 all local ${trajdump} dumplin_bonds/bonds.${simname}.* index c_1[1] c_1[2] c_1[3] #bond index, id of atom 1 and 2 in the bond, bond type
dump_modify	3 format line "%d %0.0f %0.0f %0.0f"
dump_modify	3 pad 10

compute binfo all bond/local dist
dump 33 all local ${trajdump} dumplin_bonds/bondlen.${simname}.* index c_1[1] c_1[2] c_1[3] c_binfo
dump_modify 33 pad 10

#Save bond information
dump 		4 all local ${trajdumpFast} dumplin_bonds_Sample/bonds.${simname}.* index c_1[1] c_1[2] c_1[3] #bond index, id of atom 1 and 2 in the bond, bond type
dump_modify	4 format line "%d %0.0f %0.0f %0.0f"
dump_modify	4 pad 10


restart         ${trestart} restart/${simname}.*.restart

run ${tmix}

write_restart	restart.postmix

special_bonds lj 0.0 1.0 1.0

fix             freact_creation all bond/react stabilization no reset_mol_ids no &
                react Nucleation all ${NeveryFast} ${RminF} ${RmaxFB} nucleation_pre nucleation_post ../../templates/Nucleation/map_Nucleation.txt prob 1 ${lseed1} modify_create nuc xwalls -2 &
                react Death all ${NeveryFast} ${RminF} ${RmaxFB} death_pre death_post ../../templates/Nucleation/map_Death.txt prob 1 ${lseed1} &
                react Rc1_2N_01 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_2N_01_pre Nc1_2N_01_post ../../templates/NC1_2N/01_map.txt prob v_probF ${lseed1} &
                react Rc1_2N_02 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_2N_02_pre Nc1_2N_02_post ../../templates/NC1_2N/02_map.txt prob v_probF ${lseed1} &
                react Rc1_2N_03 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_2N_03_pre Nc1_2N_03_post ../../templates/NC1_2N/03_map.txt prob v_probF ${lseed1} &
                react Rc1_2N_04 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_2N_04_pre Nc1_2N_04_post ../../templates/NC1_2N/04_map.txt prob v_probF ${lseed1} &
                react Rc1_2N_05 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_2N_05_pre Nc1_2N_05_post ../../templates/NC1_2N/05_map.txt prob v_probF ${lseed1} &
                react Rc1_2N_06 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_2N_06_pre Nc1_2N_06_post ../../templates/NC1_2N/06_map.txt prob v_probF ${lseed1} &
                react Rc1_2N_07 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_2N_07_pre Nc1_2N_07_post ../../templates/NC1_2N/07_map.txt prob v_probF ${lseed1} &
                react Rc1_2N_08 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_2N_08_pre Nc1_2N_08_post ../../templates/NC1_2N/08_map.txt prob v_probF ${lseed1} &
                react Rc1_2N_09 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_2N_09_pre Nc1_2N_09_post ../../templates/NC1_2N/09_map.txt prob v_probF ${lseed1} &
                react Rc1_2N_10 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_2N_10_pre Nc1_2N_10_post ../../templates/NC1_2N/10_map.txt prob v_probF ${lseed1} &
                react Rc1_2N_11 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_2N_11_pre Nc1_2N_11_post ../../templates/NC1_2N/11_map.txt prob v_probF ${lseed1} &
                react Rc1_3N_01 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_3N_01_pre Nc1_3N_01_post ../../templates/NC1_3N/01_map.txt prob v_probF ${lseed1} &
                react Rc1_3N_02 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_3N_02_pre Nc1_3N_02_post ../../templates/NC1_3N/02_map.txt prob v_probF ${lseed1} &
                react Rc1_3N_03 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_3N_03_pre Nc1_3N_03_post ../../templates/NC1_3N/03_map.txt prob v_probF ${lseed1} &
                react Rc1_3N_04 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_3N_04_pre Nc1_3N_04_post ../../templates/NC1_3N/04_map.txt prob v_probF ${lseed1} &
                react Rc1_3N_05 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_3N_05_pre Nc1_3N_05_post ../../templates/NC1_3N/05_map.txt prob v_probF ${lseed1} &
                react Rc1_3N_06 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_3N_06_pre Nc1_3N_06_post ../../templates/NC1_3N/06_map.txt prob v_probF ${lseed1} &
                react Rc1_3N_07 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_3N_07_pre Nc1_3N_07_post ../../templates/NC1_3N/07_map.txt prob v_probF ${lseed1} &
                react Rc1_3N_08 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_3N_08_pre Nc1_3N_08_post ../../templates/NC1_3N/08_map.txt prob v_probF ${lseed1} &
                react Rc1_3N_09 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_3N_09_pre Nc1_3N_09_post ../../templates/NC1_3N/09_map.txt prob v_probF ${lseed1} &
                react Rc1_3N_10 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_3N_10_pre Nc1_3N_10_post ../../templates/NC1_3N/10_map.txt prob v_probF ${lseed1} &
                react Rc1_4N_01 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_4N_01_pre Nc1_4N_01_post ../../templates/NC1_4N/01_map.txt prob v_probF ${lseed1} &
                react Rc1_4N_02 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_4N_02_pre Nc1_4N_02_post ../../templates/NC1_4N/02_map.txt prob v_probF ${lseed1} &
                react Rc1_4N_03 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_4N_03_pre Nc1_4N_03_post ../../templates/NC1_4N/03_map.txt prob v_probF ${lseed1} &
                react Rc1_4N_04 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_4N_04_pre Nc1_4N_04_post ../../templates/NC1_4N/04_map.txt prob v_probF ${lseed1} &
                react Rc1_4N_05 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_4N_05_pre Nc1_4N_05_post ../../templates/NC1_4N/05_map.txt prob v_probF ${lseed1} &
                react Rc1_4N_06 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_4N_06_pre Nc1_4N_06_post ../../templates/NC1_4N/06_map.txt prob v_probF ${lseed1} &
                react Rc1_4N_07 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_4N_07_pre Nc1_4N_07_post ../../templates/NC1_4N/07_map.txt prob v_probF ${lseed1} &
                react Rc1_4N_08 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_4N_08_pre Nc1_4N_08_post ../../templates/NC1_4N/08_map.txt prob v_probF ${lseed1} &
                react Rc1_4N_09 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_4N_09_pre Nc1_4N_09_post ../../templates/NC1_4N/09_map.txt prob v_probF ${lseed1} &
                react Rc1_4N_10 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_4N_10_pre Nc1_4N_10_post ../../templates/NC1_4N/10_map.txt prob v_probF ${lseed1} &
                react Rc1_4N_11 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_4N_11_pre Nc1_4N_11_post ../../templates/NC1_4N/11_map.txt prob v_probF ${lseed1} &
                react Rc1_4N_12 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_4N_12_pre Nc1_4N_12_post ../../templates/NC1_4N/12_map.txt prob v_probF ${lseed1} &
                react Rc1_4N_13 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_4N_13_pre Nc1_4N_13_post ../../templates/NC1_4N/13_map.txt prob v_probF ${lseed1} &
                react Rc1_4N_14 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_4N_14_pre Nc1_4N_14_post ../../templates/NC1_4N/14_map.txt prob v_probF ${lseed1} &
                react Rc1_4N_15 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_4N_15_pre Nc1_4N_15_post ../../templates/NC1_4N/15_map.txt prob v_probF ${lseed1} &
                react Rc1_4N_16 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_4N_16_pre Nc1_4N_16_post ../../templates/NC1_4N/16_map.txt prob v_probF ${lseed1} &
                react Rc1_4N_17 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_4N_17_pre Nc1_4N_17_post ../../templates/NC1_4N/17_map.txt prob v_probF ${lseed1} &
                react Rc1_5N_01 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_5N_01_pre Nc1_5N_01_post ../../templates/NC1_5N/01_map.txt prob v_probF ${lseed1} &
                react Rc1_5N_02 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_5N_02_pre Nc1_5N_02_post ../../templates/NC1_5N/02_map.txt prob v_probF ${lseed1} &
                react Rc1_5N_03 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_5N_03_pre Nc1_5N_03_post ../../templates/NC1_5N/03_map.txt prob v_probF ${lseed1} &
                react Rc1_5N_04 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_5N_04_pre Nc1_5N_04_post ../../templates/NC1_5N/04_map.txt prob v_probF ${lseed1} &
                react Rc1_5N_05 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_5N_05_pre Nc1_5N_05_post ../../templates/NC1_5N/05_map.txt prob v_probF ${lseed1} &
                react Rc1_5N_06 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_5N_06_pre Nc1_5N_06_post ../../templates/NC1_5N/06_map.txt prob v_probF ${lseed1} &
                react Rc1_5N_07 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_5N_07_pre Nc1_5N_07_post ../../templates/NC1_5N/07_map.txt prob v_probF ${lseed1} &
                react Rc1_6N_01 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_6N_01_pre Nc1_6N_01_post ../../templates/NC1_6N/01_map.txt prob v_probF ${lseed1} &
                react Rc1_6N_02 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_6N_02_pre Nc1_6N_02_post ../../templates/NC1_6N/02_map.txt prob v_probF ${lseed1} &
                react Rc1_6N_03 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_6N_03_pre Nc1_6N_03_post ../../templates/NC1_6N/03_map.txt prob v_probF ${lseed1} &
                react Rc1_6N_04 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_6N_04_pre Nc1_6N_04_post ../../templates/NC1_6N/04_map.txt prob v_probF ${lseed1} &
                react Rc1_6N_05 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_6N_05_pre Nc1_6N_05_post ../../templates/NC1_6N/05_map.txt prob v_probF ${lseed1} &
                react Rc1_6N_06 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_6N_06_pre Nc1_6N_06_post ../../templates/NC1_6N/06_map.txt prob v_probF ${lseed1} &
                react Rc1_6N_07 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_6N_07_pre Nc1_6N_07_post ../../templates/NC1_6N/07_map.txt prob v_probF ${lseed1} &
                react Rc1_6N_08 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_6N_08_pre Nc1_6N_08_post ../../templates/NC1_6N/08_map.txt prob v_probF ${lseed1} &
                react Rc1_6N_09 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_6N_09_pre Nc1_6N_09_post ../../templates/NC1_6N/09_map.txt prob v_probF ${lseed1} &
                react Rc1_6N_10 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_6N_10_pre Nc1_6N_10_post ../../templates/NC1_6N/10_map.txt prob v_probF ${lseed1} &
                react Rc1_6N_11 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_6N_11_pre Nc1_6N_11_post ../../templates/NC1_6N/11_map.txt prob v_probF ${lseed1} &
                react Rc1_6N_12 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_6N_12_pre Nc1_6N_12_post ../../templates/NC1_6N/12_map.txt prob v_probF ${lseed1} &
                react Rc1_6N_13 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_6N_13_pre Nc1_6N_13_post ../../templates/NC1_6N/13_map.txt prob v_probF ${lseed1} &
                react Rc1_6N_14 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_6N_14_pre Nc1_6N_14_post ../../templates/NC1_6N/14_map.txt prob v_probF ${lseed1} &
                react Rc1_6N_15 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_6N_15_pre Nc1_6N_15_post ../../templates/NC1_6N/15_map.txt prob v_probF ${lseed1} &
                react Rc1_7N_01 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_7N_01_pre Nc1_7N_01_post ../../templates/NC1_7N/01_map.txt prob v_probF ${lseed1} &
                react Rc1_7N_02 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_7N_02_pre Nc1_7N_02_post ../../templates/NC1_7N/02_map.txt prob v_probF ${lseed1} &
                react Rc1_7N_03 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_7N_03_pre Nc1_7N_03_post ../../templates/NC1_7N/03_map.txt prob v_probF ${lseed1} &
                react Rc1_7N_04 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_7N_04_pre Nc1_7N_04_post ../../templates/NC1_7N/04_map.txt prob v_probF ${lseed1} &
                react Rc1_7N_05 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_7N_05_pre Nc1_7N_05_post ../../templates/NC1_7N/05_map.txt prob v_probF ${lseed1} &
                react Rc1_8N_01 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_8N_01_pre Nc1_8N_01_post ../../templates/NC1_8N/01_map.txt prob v_probF ${lseed1} &
                react Rc1_loop_01 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_loop_01_pre Nc1_loop_01_post ../../templates/NC1_loop/01_map.txt prob v_probF ${lseed1} &
                react Rc1_loop_02 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_loop_02_pre Nc1_loop_02_post ../../templates/NC1_loop/02_map.txt prob v_probF ${lseed1} &
                react Rc1_loop_03 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_loop_03_pre Nc1_loop_03_post ../../templates/NC1_loop/03_map.txt prob v_probF ${lseed1} &
                react Rc1_loop_04 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_loop_04_pre Nc1_loop_04_post ../../templates/NC1_loop/04_map.txt prob v_probF ${lseed1} &
                react Rc1_loop_05 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_loop_05_pre Nc1_loop_05_post ../../templates/NC1_loop/05_map.txt prob v_probF ${lseed1} &
                react Rc1_loop_06 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_loop_06_pre Nc1_loop_06_post ../../templates/NC1_loop/06_map.txt prob v_probF ${lseed1} &
                react Rc1_loop_07 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_loop_07_pre Nc1_loop_07_post ../../templates/NC1_loop/07_map.txt prob v_probF ${lseed1} &
                react Rc1_loop_08 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_loop_08_pre Nc1_loop_08_post ../../templates/NC1_loop/08_map.txt prob v_probF ${lseed1} &
                react Rc1_loop_09 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_loop_09_pre Nc1_loop_09_post ../../templates/NC1_loop/09_map.txt prob v_probF ${lseed1} &
                react Rc1_loop_10 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_loop_10_pre Nc1_loop_10_post ../../templates/NC1_loop/10_map.txt prob v_probF ${lseed1} &
                react Rc1_loop_11 all ${NeverySlow} ${RminF} ${RmaxF} Nc1_loop_11_pre Nc1_loop_11_post ../../templates/NC1_loop/11_map.txt prob v_probF ${lseed1} &
                react MonoToDimer_01 all ${NeverySlow} ${RminF} ${RmaxF_7S} MonoToDimer_01_pre MonoToDimer_01_post ../../templates/MonomerToDimer/01_map.txt prob v_probF_7S ${lseed2} &
                react MonoToDimer_02 all ${NeverySlow} ${RminF} ${RmaxF_7S} MonoToDimer_02_pre MonoToDimer_02_post ../../templates/MonomerToDimer/02_map.txt prob v_probF_7S ${lseed2} &
                react MonoToDimer_03 all ${NeverySlow} ${RminF} ${RmaxF_7S} MonoToDimer_03_pre MonoToDimer_03_post ../../templates/MonomerToDimer/03_map.txt prob v_probF_7S ${lseed2} &
                react MonoToDimer_04 all ${NeverySlow} ${RminF} ${RmaxF_7S} MonoToDimer_04_pre MonoToDimer_04_post ../../templates/MonomerToDimer/04_map.txt prob v_probF_7S ${lseed2} &
                react DimerToLineTrimer_01 all ${NeverySlow} ${RminF} ${RmaxF_7S} DimerToLineTrimer_01_pre DimerToLineTrimer_01_post ../../templates/DimerToLineTrimer/01_map.txt prob v_probF_7S ${lseed3} &
                react DimerToLineTrimer_02 all ${NeverySlow} ${RminF} ${RmaxF_7S} DimerToLineTrimer_02_pre DimerToLineTrimer_02_post ../../templates/DimerToLineTrimer/02_map.txt prob v_probF_7S ${lseed3} &
                react DimerToLineTrimer_03 all ${NeverySlow} ${RminF} ${RmaxF_7S} DimerToLineTrimer_03_pre DimerToLineTrimer_03_post ../../templates/DimerToLineTrimer/03_map.txt prob v_probF_7S ${lseed3} &
                react DimerToLineTrimer_04 all ${NeverySlow} ${RminF} ${RmaxF_7S} DimerToLineTrimer_04_pre DimerToLineTrimer_04_post ../../templates/DimerToLineTrimer/04_map.txt prob v_probF_7S ${lseed3} &
                react DimerToLineTrimer_05 all ${NeverySlow} ${RminF} ${RmaxF_7S} DimerToLineTrimer_05_pre DimerToLineTrimer_05_post ../../templates/DimerToLineTrimer/05_map.txt prob v_probF_7S ${lseed3} &
                react DimerToLineTrimer_06 all ${NeverySlow} ${RminF} ${RmaxF_7S} DimerToLineTrimer_06_pre DimerToLineTrimer_06_post ../../templates/DimerToLineTrimer/06_map.txt prob v_probF_7S ${lseed3} &
                react DimerToLineTrimer_07 all ${NeverySlow} ${RminF} ${RmaxF_7S} DimerToLineTrimer_07_pre DimerToLineTrimer_07_post ../../templates/DimerToLineTrimer/07_map.txt prob v_probF_7S ${lseed3} &
                react DimerToLineTrimer_08 all ${NeverySlow} ${RminF} ${RmaxF_7S} DimerToLineTrimer_08_pre DimerToLineTrimer_08_post ../../templates/DimerToLineTrimer/08_map.txt prob v_probF_7S ${lseed3} &
                react DimerToLineTrimer_09 all ${NeverySlow} ${RminF} ${RmaxF_7S} DimerToLineTrimer_09_pre DimerToLineTrimer_09_post ../../templates/DimerToLineTrimer/09_map.txt prob v_probF_7S ${lseed3} &
                react DimerToLineTrimer_10 all ${NeverySlow} ${RminF} ${RmaxF_7S} DimerToLineTrimer_10_pre DimerToLineTrimer_10_post ../../templates/DimerToLineTrimer/10_map.txt prob v_probF_7S ${lseed3} &
                react DimerToLineTrimer_11 all ${NeverySlow} ${RminF} ${RmaxF_7S} DimerToLineTrimer_11_pre DimerToLineTrimer_11_post ../../templates/DimerToLineTrimer/11_map.txt prob v_probF_7S ${lseed3} &
                react DimerToLineTrimer_12 all ${NeverySlow} ${RminF} ${RmaxF_7S} DimerToLineTrimer_12_pre DimerToLineTrimer_12_post ../../templates/DimerToLineTrimer/12_map.txt prob v_probF_7S ${lseed3} &
                react DimerToLineTrimer_13 all ${NeverySlow} ${RminF} ${RmaxF_7S} DimerToLineTrimer_13_pre DimerToLineTrimer_13_post ../../templates/DimerToLineTrimer/13_map.txt prob v_probF_7S ${lseed3} &
                react DimerToLineTrimer_14 all ${NeverySlow} ${RminF} ${RmaxF_7S} DimerToLineTrimer_14_pre DimerToLineTrimer_14_post ../../templates/DimerToLineTrimer/14_map.txt prob v_probF_7S ${lseed3} &
                react LineTrimerToTriangleTrimer_01 all ${NeveryFast} ${RminF} ${RmaxF_7S} LineTrimerToTriangleTrimer_01_pre LineTrimerToTriangleTrimer_01_post ../../templates/LineTrimerToTriangleTrimer/01_map.txt prob v_probF_7S ${lseed4} &
                react LineTrimerToTriangleTrimer_02 all ${NeveryFast} ${RminF} ${RmaxF_7S} LineTrimerToTriangleTrimer_02_pre LineTrimerToTriangleTrimer_02_post ../../templates/LineTrimerToTriangleTrimer/02_map.txt prob v_probF_7S ${lseed4} &
                react LineTrimerToTriangleTrimer_03 all ${NeveryFast} ${RminF} ${RmaxF_7S} LineTrimerToTriangleTrimer_03_pre LineTrimerToTriangleTrimer_03_post ../../templates/LineTrimerToTriangleTrimer/03_map.txt prob v_probF_7S ${lseed4} &
                react LineTrimerToTriangleTrimer_04 all ${NeveryFast} ${RminF} ${RmaxF_7S} LineTrimerToTriangleTrimer_04_pre LineTrimerToTriangleTrimer_04_post ../../templates/LineTrimerToTriangleTrimer/04_map.txt prob v_probF_7S ${lseed4} &
                react LineTrimerToTriangleTrimer_05 all ${NeveryFast} ${RminF} ${RmaxF_7S} LineTrimerToTriangleTrimer_05_pre LineTrimerToTriangleTrimer_05_post ../../templates/LineTrimerToTriangleTrimer/05_map.txt prob v_probF_7S ${lseed4} &
                react LineTrimerToTriangleTrimer_06 all ${NeveryFast} ${RminF} ${RmaxF_7S} LineTrimerToTriangleTrimer_06_pre LineTrimerToTriangleTrimer_06_post ../../templates/LineTrimerToTriangleTrimer/06_map.txt prob v_probF_7S ${lseed4} &
                react LineTrimerToTriangleTrimer_07 all ${NeveryFast} ${RminF} ${RmaxF_7S} LineTrimerToTriangleTrimer_07_pre LineTrimerToTriangleTrimer_07_post ../../templates/LineTrimerToTriangleTrimer/07_map.txt prob v_probF_7S ${lseed4} &
                react LineTrimerToTriangleTrimer_08 all ${NeveryFast} ${RminF} ${RmaxF_7S} LineTrimerToTriangleTrimer_08_pre LineTrimerToTriangleTrimer_08_post ../../templates/LineTrimerToTriangleTrimer/08_map.txt prob v_probF_7S ${lseed4} &
                react LineTrimerToTriangleTrimer_09 all ${NeveryFast} ${RminF} ${RmaxF_7S} LineTrimerToTriangleTrimer_09_pre LineTrimerToTriangleTrimer_09_post ../../templates/LineTrimerToTriangleTrimer/09_map.txt prob v_probF_7S ${lseed4} &
                react LineTrimerToTriangleTrimer_10 all ${NeveryFast} ${RminF} ${RmaxF_7S} LineTrimerToTriangleTrimer_10_pre LineTrimerToTriangleTrimer_10_post ../../templates/LineTrimerToTriangleTrimer/10_map.txt prob v_probF_7S ${lseed4} &
                react TriangleTrimerToOpenTetramer_01 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_01_pre TriangleTrimerToOpenTetramer_01_post ../../templates/TriangleTrimerToOpenTetramer/01_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_02 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_02_pre TriangleTrimerToOpenTetramer_02_post ../../templates/TriangleTrimerToOpenTetramer/02_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_03 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_03_pre TriangleTrimerToOpenTetramer_03_post ../../templates/TriangleTrimerToOpenTetramer/03_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_04 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_04_pre TriangleTrimerToOpenTetramer_04_post ../../templates/TriangleTrimerToOpenTetramer/04_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_05 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_05_pre TriangleTrimerToOpenTetramer_05_post ../../templates/TriangleTrimerToOpenTetramer/05_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_06 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_06_pre TriangleTrimerToOpenTetramer_06_post ../../templates/TriangleTrimerToOpenTetramer/06_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_07 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_07_pre TriangleTrimerToOpenTetramer_07_post ../../templates/TriangleTrimerToOpenTetramer/07_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_08 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_08_pre TriangleTrimerToOpenTetramer_08_post ../../templates/TriangleTrimerToOpenTetramer/08_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_09 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_09_pre TriangleTrimerToOpenTetramer_09_post ../../templates/TriangleTrimerToOpenTetramer/09_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_10 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_10_pre TriangleTrimerToOpenTetramer_10_post ../../templates/TriangleTrimerToOpenTetramer/10_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_11 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_11_pre TriangleTrimerToOpenTetramer_11_post ../../templates/TriangleTrimerToOpenTetramer/11_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_12 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_12_pre TriangleTrimerToOpenTetramer_12_post ../../templates/TriangleTrimerToOpenTetramer/12_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_13 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_13_pre TriangleTrimerToOpenTetramer_13_post ../../templates/TriangleTrimerToOpenTetramer/13_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_14 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_14_pre TriangleTrimerToOpenTetramer_14_post ../../templates/TriangleTrimerToOpenTetramer/14_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_15 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_15_pre TriangleTrimerToOpenTetramer_15_post ../../templates/TriangleTrimerToOpenTetramer/15_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_16 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_16_pre TriangleTrimerToOpenTetramer_16_post ../../templates/TriangleTrimerToOpenTetramer/16_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_17 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_17_pre TriangleTrimerToOpenTetramer_17_post ../../templates/TriangleTrimerToOpenTetramer/17_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_18 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_18_pre TriangleTrimerToOpenTetramer_18_post ../../templates/TriangleTrimerToOpenTetramer/18_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_19 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_19_pre TriangleTrimerToOpenTetramer_19_post ../../templates/TriangleTrimerToOpenTetramer/19_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_20 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_20_pre TriangleTrimerToOpenTetramer_20_post ../../templates/TriangleTrimerToOpenTetramer/20_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_21 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_21_pre TriangleTrimerToOpenTetramer_21_post ../../templates/TriangleTrimerToOpenTetramer/21_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_22 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_22_pre TriangleTrimerToOpenTetramer_22_post ../../templates/TriangleTrimerToOpenTetramer/22_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_23 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_23_pre TriangleTrimerToOpenTetramer_23_post ../../templates/TriangleTrimerToOpenTetramer/23_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_24 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_24_pre TriangleTrimerToOpenTetramer_24_post ../../templates/TriangleTrimerToOpenTetramer/24_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_25 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_25_pre TriangleTrimerToOpenTetramer_25_post ../../templates/TriangleTrimerToOpenTetramer/25_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_26 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_26_pre TriangleTrimerToOpenTetramer_26_post ../../templates/TriangleTrimerToOpenTetramer/26_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_27 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_27_pre TriangleTrimerToOpenTetramer_27_post ../../templates/TriangleTrimerToOpenTetramer/27_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_28 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_28_pre TriangleTrimerToOpenTetramer_28_post ../../templates/TriangleTrimerToOpenTetramer/28_map.txt prob v_probF_7S ${lseed5} &
                react TriangleTrimerToOpenTetramer_29 all ${NeverySlow} ${RminF} ${RmaxF_7S} TriangleTrimerToOpenTetramer_29_pre TriangleTrimerToOpenTetramer_29_post ../../templates/TriangleTrimerToOpenTetramer/29_map.txt prob v_probF_7S ${lseed5} &
                react OpenTetramerToRhomboidTetramer_01 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_01_pre OpenTetramerToRhomboidTetramer_01_post ../../templates/OpenTetramerToRhomboidTetramer/01_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_02 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_02_pre OpenTetramerToRhomboidTetramer_02_post ../../templates/OpenTetramerToRhomboidTetramer/02_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_03 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_03_pre OpenTetramerToRhomboidTetramer_03_post ../../templates/OpenTetramerToRhomboidTetramer/03_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_04 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_04_pre OpenTetramerToRhomboidTetramer_04_post ../../templates/OpenTetramerToRhomboidTetramer/04_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_05 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_05_pre OpenTetramerToRhomboidTetramer_05_post ../../templates/OpenTetramerToRhomboidTetramer/05_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_06 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_06_pre OpenTetramerToRhomboidTetramer_06_post ../../templates/OpenTetramerToRhomboidTetramer/06_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_07 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_07_pre OpenTetramerToRhomboidTetramer_07_post ../../templates/OpenTetramerToRhomboidTetramer/07_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_08 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_08_pre OpenTetramerToRhomboidTetramer_08_post ../../templates/OpenTetramerToRhomboidTetramer/08_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_09 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_09_pre OpenTetramerToRhomboidTetramer_09_post ../../templates/OpenTetramerToRhomboidTetramer/09_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_10 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_10_pre OpenTetramerToRhomboidTetramer_10_post ../../templates/OpenTetramerToRhomboidTetramer/10_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_11 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_11_pre OpenTetramerToRhomboidTetramer_11_post ../../templates/OpenTetramerToRhomboidTetramer/11_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_12 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_12_pre OpenTetramerToRhomboidTetramer_12_post ../../templates/OpenTetramerToRhomboidTetramer/12_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_13 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_13_pre OpenTetramerToRhomboidTetramer_13_post ../../templates/OpenTetramerToRhomboidTetramer/13_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_14 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_14_pre OpenTetramerToRhomboidTetramer_14_post ../../templates/OpenTetramerToRhomboidTetramer/14_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_15 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_15_pre OpenTetramerToRhomboidTetramer_15_post ../../templates/OpenTetramerToRhomboidTetramer/15_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_16 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_16_pre OpenTetramerToRhomboidTetramer_16_post ../../templates/OpenTetramerToRhomboidTetramer/16_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_17 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_17_pre OpenTetramerToRhomboidTetramer_17_post ../../templates/OpenTetramerToRhomboidTetramer/17_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_18 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_18_pre OpenTetramerToRhomboidTetramer_18_post ../../templates/OpenTetramerToRhomboidTetramer/18_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_19 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_19_pre OpenTetramerToRhomboidTetramer_19_post ../../templates/OpenTetramerToRhomboidTetramer/19_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_20 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_20_pre OpenTetramerToRhomboidTetramer_20_post ../../templates/OpenTetramerToRhomboidTetramer/20_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_21 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_21_pre OpenTetramerToRhomboidTetramer_21_post ../../templates/OpenTetramerToRhomboidTetramer/21_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_22 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_22_pre OpenTetramerToRhomboidTetramer_22_post ../../templates/OpenTetramerToRhomboidTetramer/22_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_23 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_23_pre OpenTetramerToRhomboidTetramer_23_post ../../templates/OpenTetramerToRhomboidTetramer/23_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_24 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_24_pre OpenTetramerToRhomboidTetramer_24_post ../../templates/OpenTetramerToRhomboidTetramer/24_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_25 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_25_pre OpenTetramerToRhomboidTetramer_25_post ../../templates/OpenTetramerToRhomboidTetramer/25_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_26 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_26_pre OpenTetramerToRhomboidTetramer_26_post ../../templates/OpenTetramerToRhomboidTetramer/26_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_27 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_27_pre OpenTetramerToRhomboidTetramer_27_post ../../templates/OpenTetramerToRhomboidTetramer/27_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_28 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_28_pre OpenTetramerToRhomboidTetramer_28_post ../../templates/OpenTetramerToRhomboidTetramer/28_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_29 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_29_pre OpenTetramerToRhomboidTetramer_29_post ../../templates/OpenTetramerToRhomboidTetramer/29_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_30 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_30_pre OpenTetramerToRhomboidTetramer_30_post ../../templates/OpenTetramerToRhomboidTetramer/30_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_31 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_31_pre OpenTetramerToRhomboidTetramer_31_post ../../templates/OpenTetramerToRhomboidTetramer/31_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_32 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_32_pre OpenTetramerToRhomboidTetramer_32_post ../../templates/OpenTetramerToRhomboidTetramer/32_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_33 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_33_pre OpenTetramerToRhomboidTetramer_33_post ../../templates/OpenTetramerToRhomboidTetramer/33_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_34 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_34_pre OpenTetramerToRhomboidTetramer_34_post ../../templates/OpenTetramerToRhomboidTetramer/34_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_35 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_35_pre OpenTetramerToRhomboidTetramer_35_post ../../templates/OpenTetramerToRhomboidTetramer/35_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_36 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_36_pre OpenTetramerToRhomboidTetramer_36_post ../../templates/OpenTetramerToRhomboidTetramer/36_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_37 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_37_pre OpenTetramerToRhomboidTetramer_37_post ../../templates/OpenTetramerToRhomboidTetramer/37_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_38 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_38_pre OpenTetramerToRhomboidTetramer_38_post ../../templates/OpenTetramerToRhomboidTetramer/38_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_39 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_39_pre OpenTetramerToRhomboidTetramer_39_post ../../templates/OpenTetramerToRhomboidTetramer/39_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_40 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_40_pre OpenTetramerToRhomboidTetramer_40_post ../../templates/OpenTetramerToRhomboidTetramer/40_map.txt prob v_probF_7S ${lseed6} &
                react OpenTetramerToRhomboidTetramer_41 all ${NeveryFast} ${RminF} ${RmaxF_7S} OpenTetramerToRhomboidTetramer_41_pre OpenTetramerToRhomboidTetramer_41_post ../../templates/OpenTetramerToRhomboidTetramer/41_map.txt prob v_probF_7S ${lseed6} &
                react RhomboidTetramerToSquareTetramer_01 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_01_pre RhomboidTetramerToSquareTetramer_01_post ../../templates/RhomboidTetramerToSquareTetramer/01_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_02 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_02_pre RhomboidTetramerToSquareTetramer_02_post ../../templates/RhomboidTetramerToSquareTetramer/02_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_03 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_03_pre RhomboidTetramerToSquareTetramer_03_post ../../templates/RhomboidTetramerToSquareTetramer/03_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_04 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_04_pre RhomboidTetramerToSquareTetramer_04_post ../../templates/RhomboidTetramerToSquareTetramer/04_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_05 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_05_pre RhomboidTetramerToSquareTetramer_05_post ../../templates/RhomboidTetramerToSquareTetramer/05_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_06 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_06_pre RhomboidTetramerToSquareTetramer_06_post ../../templates/RhomboidTetramerToSquareTetramer/06_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_07 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_07_pre RhomboidTetramerToSquareTetramer_07_post ../../templates/RhomboidTetramerToSquareTetramer/07_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_08 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_08_pre RhomboidTetramerToSquareTetramer_08_post ../../templates/RhomboidTetramerToSquareTetramer/08_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_09 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_09_pre RhomboidTetramerToSquareTetramer_09_post ../../templates/RhomboidTetramerToSquareTetramer/09_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_10 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_10_pre RhomboidTetramerToSquareTetramer_10_post ../../templates/RhomboidTetramerToSquareTetramer/10_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_11 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_11_pre RhomboidTetramerToSquareTetramer_11_post ../../templates/RhomboidTetramerToSquareTetramer/11_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_12 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_12_pre RhomboidTetramerToSquareTetramer_12_post ../../templates/RhomboidTetramerToSquareTetramer/12_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_13 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_13_pre RhomboidTetramerToSquareTetramer_13_post ../../templates/RhomboidTetramerToSquareTetramer/13_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_14 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_14_pre RhomboidTetramerToSquareTetramer_14_post ../../templates/RhomboidTetramerToSquareTetramer/14_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_15 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_15_pre RhomboidTetramerToSquareTetramer_15_post ../../templates/RhomboidTetramerToSquareTetramer/15_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_16 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_16_pre RhomboidTetramerToSquareTetramer_16_post ../../templates/RhomboidTetramerToSquareTetramer/16_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_17 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_17_pre RhomboidTetramerToSquareTetramer_17_post ../../templates/RhomboidTetramerToSquareTetramer/17_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_18 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_18_pre RhomboidTetramerToSquareTetramer_18_post ../../templates/RhomboidTetramerToSquareTetramer/18_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_19 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_19_pre RhomboidTetramerToSquareTetramer_19_post ../../templates/RhomboidTetramerToSquareTetramer/19_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_20 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_20_pre RhomboidTetramerToSquareTetramer_20_post ../../templates/RhomboidTetramerToSquareTetramer/20_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_21 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_21_pre RhomboidTetramerToSquareTetramer_21_post ../../templates/RhomboidTetramerToSquareTetramer/21_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_22 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_22_pre RhomboidTetramerToSquareTetramer_22_post ../../templates/RhomboidTetramerToSquareTetramer/22_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_23 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_23_pre RhomboidTetramerToSquareTetramer_23_post ../../templates/RhomboidTetramerToSquareTetramer/23_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_24 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_24_pre RhomboidTetramerToSquareTetramer_24_post ../../templates/RhomboidTetramerToSquareTetramer/24_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_25 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_25_pre RhomboidTetramerToSquareTetramer_25_post ../../templates/RhomboidTetramerToSquareTetramer/25_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_26 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_26_pre RhomboidTetramerToSquareTetramer_26_post ../../templates/RhomboidTetramerToSquareTetramer/26_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_27 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_27_pre RhomboidTetramerToSquareTetramer_27_post ../../templates/RhomboidTetramerToSquareTetramer/27_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_28 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_28_pre RhomboidTetramerToSquareTetramer_28_post ../../templates/RhomboidTetramerToSquareTetramer/28_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_29 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_29_pre RhomboidTetramerToSquareTetramer_29_post ../../templates/RhomboidTetramerToSquareTetramer/29_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_30 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_30_pre RhomboidTetramerToSquareTetramer_30_post ../../templates/RhomboidTetramerToSquareTetramer/30_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_31 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_31_pre RhomboidTetramerToSquareTetramer_31_post ../../templates/RhomboidTetramerToSquareTetramer/31_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_32 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_32_pre RhomboidTetramerToSquareTetramer_32_post ../../templates/RhomboidTetramerToSquareTetramer/32_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_33 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_33_pre RhomboidTetramerToSquareTetramer_33_post ../../templates/RhomboidTetramerToSquareTetramer/33_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_34 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_34_pre RhomboidTetramerToSquareTetramer_34_post ../../templates/RhomboidTetramerToSquareTetramer/34_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_35 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_35_pre RhomboidTetramerToSquareTetramer_35_post ../../templates/RhomboidTetramerToSquareTetramer/35_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_36 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_36_pre RhomboidTetramerToSquareTetramer_36_post ../../templates/RhomboidTetramerToSquareTetramer/36_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_37 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_37_pre RhomboidTetramerToSquareTetramer_37_post ../../templates/RhomboidTetramerToSquareTetramer/37_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_38 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_38_pre RhomboidTetramerToSquareTetramer_38_post ../../templates/RhomboidTetramerToSquareTetramer/38_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_39 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_39_pre RhomboidTetramerToSquareTetramer_39_post ../../templates/RhomboidTetramerToSquareTetramer/39_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_40 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_40_pre RhomboidTetramerToSquareTetramer_40_post ../../templates/RhomboidTetramerToSquareTetramer/40_map.txt prob v_probF_7S ${lseed7} &
                react RhomboidTetramerToSquareTetramer_41 all ${NeveryFast} ${RminF} ${RmaxF_7S} RhomboidTetramerToSquareTetramer_41_pre RhomboidTetramerToSquareTetramer_41_post ../../templates/RhomboidTetramerToSquareTetramer/41_map.txt prob v_probF_7S ${lseed7} &      
                react BRhomboidTetramerToSquareTetramer_01 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_01_post RhomboidTetramerToSquareTetramer_01_pre ../../templates/RhomboidTetramerToSquareTetramer/01_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_02 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_02_post RhomboidTetramerToSquareTetramer_02_pre ../../templates/RhomboidTetramerToSquareTetramer/02_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_03 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_03_post RhomboidTetramerToSquareTetramer_03_pre ../../templates/RhomboidTetramerToSquareTetramer/03_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_04 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_04_post RhomboidTetramerToSquareTetramer_04_pre ../../templates/RhomboidTetramerToSquareTetramer/04_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_05 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_05_post RhomboidTetramerToSquareTetramer_05_pre ../../templates/RhomboidTetramerToSquareTetramer/05_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_06 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_06_post RhomboidTetramerToSquareTetramer_06_pre ../../templates/RhomboidTetramerToSquareTetramer/06_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_07 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_07_post RhomboidTetramerToSquareTetramer_07_pre ../../templates/RhomboidTetramerToSquareTetramer/07_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_08 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_08_post RhomboidTetramerToSquareTetramer_08_pre ../../templates/RhomboidTetramerToSquareTetramer/08_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_09 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_09_post RhomboidTetramerToSquareTetramer_09_pre ../../templates/RhomboidTetramerToSquareTetramer/09_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_10 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_10_post RhomboidTetramerToSquareTetramer_10_pre ../../templates/RhomboidTetramerToSquareTetramer/10_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_11 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_11_post RhomboidTetramerToSquareTetramer_11_pre ../../templates/RhomboidTetramerToSquareTetramer/11_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_12 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_12_post RhomboidTetramerToSquareTetramer_12_pre ../../templates/RhomboidTetramerToSquareTetramer/12_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_13 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_13_post RhomboidTetramerToSquareTetramer_13_pre ../../templates/RhomboidTetramerToSquareTetramer/13_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_14 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_14_post RhomboidTetramerToSquareTetramer_14_pre ../../templates/RhomboidTetramerToSquareTetramer/14_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_15 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_15_post RhomboidTetramerToSquareTetramer_15_pre ../../templates/RhomboidTetramerToSquareTetramer/15_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_16 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_16_post RhomboidTetramerToSquareTetramer_16_pre ../../templates/RhomboidTetramerToSquareTetramer/16_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_17 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_17_post RhomboidTetramerToSquareTetramer_17_pre ../../templates/RhomboidTetramerToSquareTetramer/17_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_18 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_18_post RhomboidTetramerToSquareTetramer_18_pre ../../templates/RhomboidTetramerToSquareTetramer/18_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_19 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_19_post RhomboidTetramerToSquareTetramer_19_pre ../../templates/RhomboidTetramerToSquareTetramer/19_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_20 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_20_post RhomboidTetramerToSquareTetramer_20_pre ../../templates/RhomboidTetramerToSquareTetramer/20_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_21 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_21_post RhomboidTetramerToSquareTetramer_21_pre ../../templates/RhomboidTetramerToSquareTetramer/21_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_22 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_22_post RhomboidTetramerToSquareTetramer_22_pre ../../templates/RhomboidTetramerToSquareTetramer/22_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_23 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_23_post RhomboidTetramerToSquareTetramer_23_pre ../../templates/RhomboidTetramerToSquareTetramer/23_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_24 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_24_post RhomboidTetramerToSquareTetramer_24_pre ../../templates/RhomboidTetramerToSquareTetramer/24_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_25 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_25_post RhomboidTetramerToSquareTetramer_25_pre ../../templates/RhomboidTetramerToSquareTetramer/25_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_26 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_26_post RhomboidTetramerToSquareTetramer_26_pre ../../templates/RhomboidTetramerToSquareTetramer/26_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_27 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_27_post RhomboidTetramerToSquareTetramer_27_pre ../../templates/RhomboidTetramerToSquareTetramer/27_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_28 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_28_post RhomboidTetramerToSquareTetramer_28_pre ../../templates/RhomboidTetramerToSquareTetramer/28_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_29 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_29_post RhomboidTetramerToSquareTetramer_29_pre ../../templates/RhomboidTetramerToSquareTetramer/29_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_30 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_30_post RhomboidTetramerToSquareTetramer_30_pre ../../templates/RhomboidTetramerToSquareTetramer/30_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_31 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_31_post RhomboidTetramerToSquareTetramer_31_pre ../../templates/RhomboidTetramerToSquareTetramer/31_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_32 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_32_post RhomboidTetramerToSquareTetramer_32_pre ../../templates/RhomboidTetramerToSquareTetramer/32_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_33 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_33_post RhomboidTetramerToSquareTetramer_33_pre ../../templates/RhomboidTetramerToSquareTetramer/33_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_34 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_34_post RhomboidTetramerToSquareTetramer_34_pre ../../templates/RhomboidTetramerToSquareTetramer/34_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_35 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_35_post RhomboidTetramerToSquareTetramer_35_pre ../../templates/RhomboidTetramerToSquareTetramer/35_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_36 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_36_post RhomboidTetramerToSquareTetramer_36_pre ../../templates/RhomboidTetramerToSquareTetramer/36_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_37 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_37_post RhomboidTetramerToSquareTetramer_37_pre ../../templates/RhomboidTetramerToSquareTetramer/37_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_38 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_38_post RhomboidTetramerToSquareTetramer_38_pre ../../templates/RhomboidTetramerToSquareTetramer/38_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_39 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_39_post RhomboidTetramerToSquareTetramer_39_pre ../../templates/RhomboidTetramerToSquareTetramer/39_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_40 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_40_post RhomboidTetramerToSquareTetramer_40_pre ../../templates/RhomboidTetramerToSquareTetramer/40_map.txt prob v_probFB_7S ${lseed8} &
                react BRhomboidTetramerToSquareTetramer_41 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} RhomboidTetramerToSquareTetramer_41_post RhomboidTetramerToSquareTetramer_41_pre ../../templates/RhomboidTetramerToSquareTetramer/41_map.txt prob v_probFB_7S ${lseed8} &
                react BOpenTetramerToRhomboidTetramer_01 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_01_post OpenTetramerToRhomboidTetramer_01_pre ../../templates/OpenTetramerToRhomboidTetramer/01_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_02 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_02_post OpenTetramerToRhomboidTetramer_02_pre ../../templates/OpenTetramerToRhomboidTetramer/02_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_03 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_03_post OpenTetramerToRhomboidTetramer_03_pre ../../templates/OpenTetramerToRhomboidTetramer/03_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_04 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_04_post OpenTetramerToRhomboidTetramer_04_pre ../../templates/OpenTetramerToRhomboidTetramer/04_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_05 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_05_post OpenTetramerToRhomboidTetramer_05_pre ../../templates/OpenTetramerToRhomboidTetramer/05_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_06 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_06_post OpenTetramerToRhomboidTetramer_06_pre ../../templates/OpenTetramerToRhomboidTetramer/06_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_07 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_07_post OpenTetramerToRhomboidTetramer_07_pre ../../templates/OpenTetramerToRhomboidTetramer/07_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_08 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_08_post OpenTetramerToRhomboidTetramer_08_pre ../../templates/OpenTetramerToRhomboidTetramer/08_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_09 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_09_post OpenTetramerToRhomboidTetramer_09_pre ../../templates/OpenTetramerToRhomboidTetramer/09_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_10 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_10_post OpenTetramerToRhomboidTetramer_10_pre ../../templates/OpenTetramerToRhomboidTetramer/10_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_11 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_11_post OpenTetramerToRhomboidTetramer_11_pre ../../templates/OpenTetramerToRhomboidTetramer/11_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_12 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_12_post OpenTetramerToRhomboidTetramer_12_pre ../../templates/OpenTetramerToRhomboidTetramer/12_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_13 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_13_post OpenTetramerToRhomboidTetramer_13_pre ../../templates/OpenTetramerToRhomboidTetramer/13_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_14 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_14_post OpenTetramerToRhomboidTetramer_14_pre ../../templates/OpenTetramerToRhomboidTetramer/14_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_15 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_15_post OpenTetramerToRhomboidTetramer_15_pre ../../templates/OpenTetramerToRhomboidTetramer/15_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_16 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_16_post OpenTetramerToRhomboidTetramer_16_pre ../../templates/OpenTetramerToRhomboidTetramer/16_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_17 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_17_post OpenTetramerToRhomboidTetramer_17_pre ../../templates/OpenTetramerToRhomboidTetramer/17_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_18 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_18_post OpenTetramerToRhomboidTetramer_18_pre ../../templates/OpenTetramerToRhomboidTetramer/18_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_19 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_19_post OpenTetramerToRhomboidTetramer_19_pre ../../templates/OpenTetramerToRhomboidTetramer/19_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_20 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_20_post OpenTetramerToRhomboidTetramer_20_pre ../../templates/OpenTetramerToRhomboidTetramer/20_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_21 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_21_post OpenTetramerToRhomboidTetramer_21_pre ../../templates/OpenTetramerToRhomboidTetramer/21_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_22 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_22_post OpenTetramerToRhomboidTetramer_22_pre ../../templates/OpenTetramerToRhomboidTetramer/22_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_23 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_23_post OpenTetramerToRhomboidTetramer_23_pre ../../templates/OpenTetramerToRhomboidTetramer/23_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_24 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_24_post OpenTetramerToRhomboidTetramer_24_pre ../../templates/OpenTetramerToRhomboidTetramer/24_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_25 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_25_post OpenTetramerToRhomboidTetramer_25_pre ../../templates/OpenTetramerToRhomboidTetramer/25_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_26 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_26_post OpenTetramerToRhomboidTetramer_26_pre ../../templates/OpenTetramerToRhomboidTetramer/26_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_27 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_27_post OpenTetramerToRhomboidTetramer_27_pre ../../templates/OpenTetramerToRhomboidTetramer/27_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_28 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_28_post OpenTetramerToRhomboidTetramer_28_pre ../../templates/OpenTetramerToRhomboidTetramer/28_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_29 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_29_post OpenTetramerToRhomboidTetramer_29_pre ../../templates/OpenTetramerToRhomboidTetramer/29_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_30 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_30_post OpenTetramerToRhomboidTetramer_30_pre ../../templates/OpenTetramerToRhomboidTetramer/30_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_31 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_31_post OpenTetramerToRhomboidTetramer_31_pre ../../templates/OpenTetramerToRhomboidTetramer/31_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_32 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_32_post OpenTetramerToRhomboidTetramer_32_pre ../../templates/OpenTetramerToRhomboidTetramer/32_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_33 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_33_post OpenTetramerToRhomboidTetramer_33_pre ../../templates/OpenTetramerToRhomboidTetramer/33_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_34 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_34_post OpenTetramerToRhomboidTetramer_34_pre ../../templates/OpenTetramerToRhomboidTetramer/34_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_35 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_35_post OpenTetramerToRhomboidTetramer_35_pre ../../templates/OpenTetramerToRhomboidTetramer/35_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_36 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_36_post OpenTetramerToRhomboidTetramer_36_pre ../../templates/OpenTetramerToRhomboidTetramer/36_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_37 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_37_post OpenTetramerToRhomboidTetramer_37_pre ../../templates/OpenTetramerToRhomboidTetramer/37_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_38 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_38_post OpenTetramerToRhomboidTetramer_38_pre ../../templates/OpenTetramerToRhomboidTetramer/38_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_39 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_39_post OpenTetramerToRhomboidTetramer_39_pre ../../templates/OpenTetramerToRhomboidTetramer/39_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_40 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_40_post OpenTetramerToRhomboidTetramer_40_pre ../../templates/OpenTetramerToRhomboidTetramer/40_map.txt prob v_probFB_7S ${lseed9} &
                react BOpenTetramerToRhomboidTetramer_41 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} OpenTetramerToRhomboidTetramer_41_post OpenTetramerToRhomboidTetramer_41_pre ../../templates/OpenTetramerToRhomboidTetramer/41_map.txt prob v_probFB_7S ${lseed9} &
                react BTriangleTrimerToOpenTetramer_01 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_01_post TriangleTrimerToOpenTetramer_01_pre ../../templates/TriangleTrimerToOpenTetramer/01_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_02 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_02_post TriangleTrimerToOpenTetramer_02_pre ../../templates/TriangleTrimerToOpenTetramer/02_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_03 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_03_post TriangleTrimerToOpenTetramer_03_pre ../../templates/TriangleTrimerToOpenTetramer/03_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_04 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_04_post TriangleTrimerToOpenTetramer_04_pre ../../templates/TriangleTrimerToOpenTetramer/04_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_05 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_05_post TriangleTrimerToOpenTetramer_05_pre ../../templates/TriangleTrimerToOpenTetramer/05_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_06 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_06_post TriangleTrimerToOpenTetramer_06_pre ../../templates/TriangleTrimerToOpenTetramer/06_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_07 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_07_post TriangleTrimerToOpenTetramer_07_pre ../../templates/TriangleTrimerToOpenTetramer/07_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_08 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_08_post TriangleTrimerToOpenTetramer_08_pre ../../templates/TriangleTrimerToOpenTetramer/08_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_09 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_09_post TriangleTrimerToOpenTetramer_09_pre ../../templates/TriangleTrimerToOpenTetramer/09_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_10 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_10_post TriangleTrimerToOpenTetramer_10_pre ../../templates/TriangleTrimerToOpenTetramer/10_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_11 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_11_post TriangleTrimerToOpenTetramer_11_pre ../../templates/TriangleTrimerToOpenTetramer/11_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_12 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_12_post TriangleTrimerToOpenTetramer_12_pre ../../templates/TriangleTrimerToOpenTetramer/12_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_13 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_13_post TriangleTrimerToOpenTetramer_13_pre ../../templates/TriangleTrimerToOpenTetramer/13_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_14 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_14_post TriangleTrimerToOpenTetramer_14_pre ../../templates/TriangleTrimerToOpenTetramer/14_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_15 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_15_post TriangleTrimerToOpenTetramer_15_pre ../../templates/TriangleTrimerToOpenTetramer/15_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_16 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_16_post TriangleTrimerToOpenTetramer_16_pre ../../templates/TriangleTrimerToOpenTetramer/16_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_17 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_17_post TriangleTrimerToOpenTetramer_17_pre ../../templates/TriangleTrimerToOpenTetramer/17_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_18 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_18_post TriangleTrimerToOpenTetramer_18_pre ../../templates/TriangleTrimerToOpenTetramer/18_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_19 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_19_post TriangleTrimerToOpenTetramer_19_pre ../../templates/TriangleTrimerToOpenTetramer/19_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_20 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_20_post TriangleTrimerToOpenTetramer_20_pre ../../templates/TriangleTrimerToOpenTetramer/20_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_21 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_21_post TriangleTrimerToOpenTetramer_21_pre ../../templates/TriangleTrimerToOpenTetramer/21_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_22 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_22_post TriangleTrimerToOpenTetramer_22_pre ../../templates/TriangleTrimerToOpenTetramer/22_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_23 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_23_post TriangleTrimerToOpenTetramer_23_pre ../../templates/TriangleTrimerToOpenTetramer/23_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_24 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_24_post TriangleTrimerToOpenTetramer_24_pre ../../templates/TriangleTrimerToOpenTetramer/24_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_25 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_25_post TriangleTrimerToOpenTetramer_25_pre ../../templates/TriangleTrimerToOpenTetramer/25_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_26 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_26_post TriangleTrimerToOpenTetramer_26_pre ../../templates/TriangleTrimerToOpenTetramer/26_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_27 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_27_post TriangleTrimerToOpenTetramer_27_pre ../../templates/TriangleTrimerToOpenTetramer/27_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_28 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_28_post TriangleTrimerToOpenTetramer_28_pre ../../templates/TriangleTrimerToOpenTetramer/28_map.txt prob v_probFB_7S ${lseed10} &
                react BTriangleTrimerToOpenTetramer_29 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} TriangleTrimerToOpenTetramer_29_post TriangleTrimerToOpenTetramer_29_pre ../../templates/TriangleTrimerToOpenTetramer/29_map.txt prob v_probFB_7S ${lseed10} &
                react BLineTrimerToTriangleTrimer_01 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} LineTrimerToTriangleTrimer_01_post LineTrimerToTriangleTrimer_01_pre ../../templates/LineTrimerToTriangleTrimer/01_map.txt prob v_probFB_7S ${lseed11} &
                react BLineTrimerToTriangleTrimer_02 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} LineTrimerToTriangleTrimer_02_post LineTrimerToTriangleTrimer_02_pre ../../templates/LineTrimerToTriangleTrimer/02_map.txt prob v_probFB_7S ${lseed11} &
                react BLineTrimerToTriangleTrimer_03 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} LineTrimerToTriangleTrimer_03_post LineTrimerToTriangleTrimer_03_pre ../../templates/LineTrimerToTriangleTrimer/03_map.txt prob v_probFB_7S ${lseed11} &
                react BLineTrimerToTriangleTrimer_04 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} LineTrimerToTriangleTrimer_04_post LineTrimerToTriangleTrimer_04_pre ../../templates/LineTrimerToTriangleTrimer/04_map.txt prob v_probFB_7S ${lseed11} &
                react BLineTrimerToTriangleTrimer_05 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} LineTrimerToTriangleTrimer_05_post LineTrimerToTriangleTrimer_05_pre ../../templates/LineTrimerToTriangleTrimer/05_map.txt prob v_probFB_7S ${lseed11} &
                react BLineTrimerToTriangleTrimer_06 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} LineTrimerToTriangleTrimer_06_post LineTrimerToTriangleTrimer_06_pre ../../templates/LineTrimerToTriangleTrimer/06_map.txt prob v_probFB_7S ${lseed11} &
                react BLineTrimerToTriangleTrimer_07 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} LineTrimerToTriangleTrimer_07_post LineTrimerToTriangleTrimer_07_pre ../../templates/LineTrimerToTriangleTrimer/07_map.txt prob v_probFB_7S ${lseed11} &
                react BLineTrimerToTriangleTrimer_08 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} LineTrimerToTriangleTrimer_08_post LineTrimerToTriangleTrimer_08_pre ../../templates/LineTrimerToTriangleTrimer/08_map.txt prob v_probFB_7S ${lseed11} &
                react BLineTrimerToTriangleTrimer_09 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} LineTrimerToTriangleTrimer_09_post LineTrimerToTriangleTrimer_09_pre ../../templates/LineTrimerToTriangleTrimer/09_map.txt prob v_probFB_7S ${lseed11} &
                react BLineTrimerToTriangleTrimer_10 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} LineTrimerToTriangleTrimer_10_post LineTrimerToTriangleTrimer_10_pre ../../templates/LineTrimerToTriangleTrimer/10_map.txt prob v_probFB_7S ${lseed11} &
                react BDimerToLineTrimer_01 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} DimerToLineTrimer_01_post DimerToLineTrimer_01_pre ../../templates/DimerToLineTrimer/01_map.txt prob v_probFB_7S ${lseed12} &
                react BDimerToLineTrimer_02 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} DimerToLineTrimer_02_post DimerToLineTrimer_02_pre ../../templates/DimerToLineTrimer/02_map.txt prob v_probFB_7S ${lseed12} &
                react BDimerToLineTrimer_03 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} DimerToLineTrimer_03_post DimerToLineTrimer_03_pre ../../templates/DimerToLineTrimer/03_map.txt prob v_probFB_7S ${lseed12} &
                react BDimerToLineTrimer_04 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} DimerToLineTrimer_04_post DimerToLineTrimer_04_pre ../../templates/DimerToLineTrimer/04_map.txt prob v_probFB_7S ${lseed12} &
                react BDimerToLineTrimer_05 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} DimerToLineTrimer_05_post DimerToLineTrimer_05_pre ../../templates/DimerToLineTrimer/05_map.txt prob v_probFB_7S ${lseed12} &
                react BDimerToLineTrimer_06 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} DimerToLineTrimer_06_post DimerToLineTrimer_06_pre ../../templates/DimerToLineTrimer/06_map.txt prob v_probFB_7S ${lseed12} &
                react BDimerToLineTrimer_07 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} DimerToLineTrimer_07_post DimerToLineTrimer_07_pre ../../templates/DimerToLineTrimer/07_map.txt prob v_probFB_7S ${lseed12} &
                react BDimerToLineTrimer_08 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} DimerToLineTrimer_08_post DimerToLineTrimer_08_pre ../../templates/DimerToLineTrimer/08_map.txt prob v_probFB_7S ${lseed12} &
                react BDimerToLineTrimer_09 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} DimerToLineTrimer_09_post DimerToLineTrimer_09_pre ../../templates/DimerToLineTrimer/09_map.txt prob v_probFB_7S ${lseed12} &
                react BDimerToLineTrimer_10 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} DimerToLineTrimer_10_post DimerToLineTrimer_10_pre ../../templates/DimerToLineTrimer/10_map.txt prob v_probFB_7S ${lseed12} &
                react BDimerToLineTrimer_11 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} DimerToLineTrimer_11_post DimerToLineTrimer_11_pre ../../templates/DimerToLineTrimer/11_map.txt prob v_probFB_7S ${lseed12} &
                react BDimerToLineTrimer_12 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} DimerToLineTrimer_12_post DimerToLineTrimer_12_pre ../../templates/DimerToLineTrimer/12_map.txt prob v_probFB_7S ${lseed12} &
                react BDimerToLineTrimer_13 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} DimerToLineTrimer_13_post DimerToLineTrimer_13_pre ../../templates/DimerToLineTrimer/13_map.txt prob v_probFB_7S ${lseed12} &
                react BDimerToLineTrimer_14 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} DimerToLineTrimer_14_post DimerToLineTrimer_14_pre ../../templates/DimerToLineTrimer/14_map.txt prob v_probFB_7S ${lseed12} &
                react BMonoToDimer_01 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} MonoToDimer_01_post MonoToDimer_01_pre ../../templates/MonomerToDimer/01_map.txt prob v_probFB_7S ${lseed13} &
                react BMonoToDimer_02 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} MonoToDimer_02_post MonoToDimer_02_pre ../../templates/MonomerToDimer/02_map.txt prob v_probFB_7S ${lseed13} &
                react BMonoToDimer_03 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} MonoToDimer_03_post MonoToDimer_03_pre ../../templates/MonomerToDimer/03_map.txt prob v_probFB_7S ${lseed13} &
                react BMonoToDimer_04 all ${NeverySlow} ${RminFB_7S} ${RmaxFB} MonoToDimer_04_post MonoToDimer_04_pre ../../templates/MonomerToDimer/04_map.txt prob v_probFB_7S ${lseed13} &
                react BRc1_2N_01 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_2N_01_post Nc1_2N_01_pre ../../templates/NC1_2N/01_map.txt prob v_probFB ${lseed14} &
                react BRc1_2N_02 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_2N_02_post Nc1_2N_02_pre ../../templates/NC1_2N/02_map.txt prob v_probFB ${lseed14} &
                react BRc1_2N_03 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_2N_03_post Nc1_2N_03_pre ../../templates/NC1_2N/03_map.txt prob v_probFB ${lseed14} &
                react BRc1_2N_04 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_2N_04_post Nc1_2N_04_pre ../../templates/NC1_2N/04_map.txt prob v_probFB ${lseed14} &
                react BRc1_2N_05 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_2N_05_post Nc1_2N_05_pre ../../templates/NC1_2N/05_map.txt prob v_probFB ${lseed14} &
                react BRc1_2N_06 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_2N_06_post Nc1_2N_06_pre ../../templates/NC1_2N/06_map.txt prob v_probFB ${lseed14} &
                react BRc1_2N_07 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_2N_07_post Nc1_2N_07_pre ../../templates/NC1_2N/07_map.txt prob v_probFB ${lseed14} &
                react BRc1_2N_08 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_2N_08_post Nc1_2N_08_pre ../../templates/NC1_2N/08_map.txt prob v_probFB ${lseed14} &
                react BRc1_2N_09 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_2N_09_post Nc1_2N_09_pre ../../templates/NC1_2N/09_map.txt prob v_probFB ${lseed14} &
                react BRc1_2N_10 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_2N_10_post Nc1_2N_10_pre ../../templates/NC1_2N/10_map.txt prob v_probFB ${lseed14} &
                react BRc1_2N_11 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_2N_11_post Nc1_2N_11_pre ../../templates/NC1_2N/11_map.txt prob v_probFB ${lseed14} &
                react BRc1_3N_01 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_3N_01_post Nc1_3N_01_pre ../../templates/NC1_3N/01_map.txt prob v_probFB ${lseed14} &
                react BRc1_3N_02 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_3N_02_post Nc1_3N_02_pre ../../templates/NC1_3N/02_map.txt prob v_probFB ${lseed14} &
                react BRc1_3N_03 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_3N_03_post Nc1_3N_03_pre ../../templates/NC1_3N/03_map.txt prob v_probFB ${lseed14} &
                react BRc1_3N_04 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_3N_04_post Nc1_3N_04_pre ../../templates/NC1_3N/04_map.txt prob v_probFB ${lseed14} &
                react BRc1_3N_05 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_3N_05_post Nc1_3N_05_pre ../../templates/NC1_3N/05_map.txt prob v_probFB ${lseed14} &
                react BRc1_3N_06 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_3N_06_post Nc1_3N_06_pre ../../templates/NC1_3N/06_map.txt prob v_probFB ${lseed14} &
                react BRc1_3N_07 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_3N_07_post Nc1_3N_07_pre ../../templates/NC1_3N/07_map.txt prob v_probFB ${lseed14} &
                react BRc1_3N_08 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_3N_08_post Nc1_3N_08_pre ../../templates/NC1_3N/08_map.txt prob v_probFB ${lseed14} &
                react BRc1_3N_09 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_3N_09_post Nc1_3N_09_pre ../../templates/NC1_3N/09_map.txt prob v_probFB ${lseed14} &
                react BRc1_3N_10 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_3N_10_post Nc1_3N_10_pre ../../templates/NC1_3N/10_map.txt prob v_probFB ${lseed14} &
                react BRc1_4N_01 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_4N_01_post Nc1_4N_01_pre ../../templates/NC1_4N/01_map.txt prob v_probFB ${lseed14} &
                react BRc1_4N_02 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_4N_02_post Nc1_4N_02_pre ../../templates/NC1_4N/02_map.txt prob v_probFB ${lseed14} &
                react BRc1_4N_03 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_4N_03_post Nc1_4N_03_pre ../../templates/NC1_4N/03_map.txt prob v_probFB ${lseed14} &
                react BRc1_4N_04 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_4N_04_post Nc1_4N_04_pre ../../templates/NC1_4N/04_map.txt prob v_probFB ${lseed14} &
                react BRc1_4N_05 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_4N_05_post Nc1_4N_05_pre ../../templates/NC1_4N/05_map.txt prob v_probFB ${lseed14} &
                react BRc1_4N_06 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_4N_06_post Nc1_4N_06_pre ../../templates/NC1_4N/06_map.txt prob v_probFB ${lseed14} &
                react BRc1_4N_07 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_4N_07_post Nc1_4N_07_pre ../../templates/NC1_4N/07_map.txt prob v_probFB ${lseed14} &
                react BRc1_4N_08 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_4N_08_post Nc1_4N_08_pre ../../templates/NC1_4N/08_map.txt prob v_probFB ${lseed14} &
                react BRc1_4N_09 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_4N_09_post Nc1_4N_09_pre ../../templates/NC1_4N/09_map.txt prob v_probFB ${lseed14} &
                react BRc1_4N_10 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_4N_10_post Nc1_4N_10_pre ../../templates/NC1_4N/10_map.txt prob v_probFB ${lseed14} &
                react BRc1_4N_11 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_4N_11_post Nc1_4N_11_pre ../../templates/NC1_4N/11_map.txt prob v_probFB ${lseed14} &
                react BRc1_4N_12 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_4N_12_post Nc1_4N_12_pre ../../templates/NC1_4N/12_map.txt prob v_probFB ${lseed14} &
                react BRc1_4N_13 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_4N_13_post Nc1_4N_13_pre ../../templates/NC1_4N/13_map.txt prob v_probFB ${lseed14} &
                react BRc1_4N_14 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_4N_14_post Nc1_4N_14_pre ../../templates/NC1_4N/14_map.txt prob v_probFB ${lseed14} &
                react BRc1_4N_15 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_4N_15_post Nc1_4N_15_pre ../../templates/NC1_4N/15_map.txt prob v_probFB ${lseed14} &
                react BRc1_4N_16 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_4N_16_post Nc1_4N_16_pre ../../templates/NC1_4N/16_map.txt prob v_probFB ${lseed14} &
                react BRc1_4N_17 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_4N_17_post Nc1_4N_17_pre ../../templates/NC1_4N/17_map.txt prob v_probFB ${lseed14} &
                react BRc1_5N_01 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_5N_01_post Nc1_5N_01_pre ../../templates/NC1_5N/01_map.txt prob v_probFB ${lseed14} &
                react BRc1_5N_02 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_5N_02_post Nc1_5N_02_pre ../../templates/NC1_5N/02_map.txt prob v_probFB ${lseed14} &
                react BRc1_5N_03 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_5N_03_post Nc1_5N_03_pre ../../templates/NC1_5N/03_map.txt prob v_probFB ${lseed14} &
                react BRc1_5N_04 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_5N_04_post Nc1_5N_04_pre ../../templates/NC1_5N/04_map.txt prob v_probFB ${lseed14} &
                react BRc1_5N_05 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_5N_05_post Nc1_5N_05_pre ../../templates/NC1_5N/05_map.txt prob v_probFB ${lseed14} &
                react BRc1_5N_06 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_5N_06_post Nc1_5N_06_pre ../../templates/NC1_5N/06_map.txt prob v_probFB ${lseed14} &
                react BRc1_5N_07 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_5N_07_post Nc1_5N_07_pre ../../templates/NC1_5N/07_map.txt prob v_probFB ${lseed14} &
                react BRc1_6N_01 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_6N_01_post Nc1_6N_01_pre ../../templates/NC1_6N/01_map.txt prob v_probFB ${lseed14} &
                react BRc1_6N_02 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_6N_02_post Nc1_6N_02_pre ../../templates/NC1_6N/02_map.txt prob v_probFB ${lseed14} &
                react BRc1_6N_03 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_6N_03_post Nc1_6N_03_pre ../../templates/NC1_6N/03_map.txt prob v_probFB ${lseed14} &
                react BRc1_6N_04 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_6N_04_post Nc1_6N_04_pre ../../templates/NC1_6N/04_map.txt prob v_probFB ${lseed14} &
                react BRc1_6N_05 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_6N_05_post Nc1_6N_05_pre ../../templates/NC1_6N/05_map.txt prob v_probFB ${lseed14} &
                react BRc1_6N_06 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_6N_06_post Nc1_6N_06_pre ../../templates/NC1_6N/06_map.txt prob v_probFB ${lseed14} &
                react BRc1_6N_07 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_6N_07_post Nc1_6N_07_pre ../../templates/NC1_6N/07_map.txt prob v_probFB ${lseed14} &
                react BRc1_6N_08 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_6N_08_post Nc1_6N_08_pre ../../templates/NC1_6N/08_map.txt prob v_probFB ${lseed14} &
                react BRc1_6N_09 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_6N_09_post Nc1_6N_09_pre ../../templates/NC1_6N/09_map.txt prob v_probFB ${lseed14} &
                react BRc1_6N_10 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_6N_10_post Nc1_6N_10_pre ../../templates/NC1_6N/10_map.txt prob v_probFB ${lseed14} &
                react BRc1_6N_11 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_6N_11_post Nc1_6N_11_pre ../../templates/NC1_6N/11_map.txt prob v_probFB ${lseed14} &
                react BRc1_6N_12 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_6N_12_post Nc1_6N_12_pre ../../templates/NC1_6N/12_map.txt prob v_probFB ${lseed14} &
                react BRc1_6N_13 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_6N_13_post Nc1_6N_13_pre ../../templates/NC1_6N/13_map.txt prob v_probFB ${lseed14} &
                react BRc1_6N_14 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_6N_14_post Nc1_6N_14_pre ../../templates/NC1_6N/14_map.txt prob v_probFB ${lseed14} &
                react BRc1_6N_15 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_6N_15_post Nc1_6N_15_pre ../../templates/NC1_6N/15_map.txt prob v_probFB ${lseed14} &
                react BRc1_7N_01 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_7N_01_post Nc1_7N_01_pre ../../templates/NC1_7N/01_map.txt prob v_probFB ${lseed14} &
                react BRc1_7N_02 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_7N_02_post Nc1_7N_02_pre ../../templates/NC1_7N/02_map.txt prob v_probFB ${lseed14} &
                react BRc1_7N_03 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_7N_03_post Nc1_7N_03_pre ../../templates/NC1_7N/03_map.txt prob v_probFB ${lseed14} &
                react BRc1_7N_04 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_7N_04_post Nc1_7N_04_pre ../../templates/NC1_7N/04_map.txt prob v_probFB ${lseed14} &
                react BRc1_7N_05 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_7N_05_post Nc1_7N_05_pre ../../templates/NC1_7N/05_map.txt prob v_probFB ${lseed14} &
                react BRc1_8N_01 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_8N_01_post Nc1_8N_01_pre ../../templates/NC1_8N/01_map.txt prob v_probFB ${lseed14} &
                react BRc1_loop_01 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_loop_01_post Nc1_loop_01_pre ../../templates/NC1_loop/01_map.txt prob v_probFB ${lseed14} &
                react BRc1_loop_02 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_loop_02_post Nc1_loop_02_pre ../../templates/NC1_loop/02_map.txt prob v_probFB ${lseed14} &
                react BRc1_loop_03 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_loop_03_post Nc1_loop_03_pre ../../templates/NC1_loop/03_map.txt prob v_probFB ${lseed14} &
                react BRc1_loop_04 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_loop_04_post Nc1_loop_04_pre ../../templates/NC1_loop/04_map.txt prob v_probFB ${lseed14} &
                react BRc1_loop_05 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_loop_05_post Nc1_loop_05_pre ../../templates/NC1_loop/05_map.txt prob v_probFB ${lseed14} &
                react BRc1_loop_06 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_loop_06_post Nc1_loop_06_pre ../../templates/NC1_loop/06_map.txt prob v_probFB ${lseed14} &
                react BRc1_loop_07 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_loop_07_post Nc1_loop_07_pre ../../templates/NC1_loop/07_map.txt prob v_probFB ${lseed14} &
                react BRc1_loop_08 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_loop_08_post Nc1_loop_08_pre ../../templates/NC1_loop/08_map.txt prob v_probFB ${lseed14} &
                react BRc1_loop_09 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_loop_09_post Nc1_loop_09_pre ../../templates/NC1_loop/09_map.txt prob v_probFB ${lseed14} &
                react BRc1_loop_10 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_loop_10_post Nc1_loop_10_pre ../../templates/NC1_loop/10_map.txt prob v_probFB ${lseed14} &
                react BRc1_loop_11 all ${NeverySlow} ${RminFB} ${RmaxFB} Nc1_loop_11_post Nc1_loop_11_pre ../../templates/NC1_loop/11_map.txt prob v_probFB ${lseed14} &
                
variable            NumNuc equal f_freact_creation[1]
variable            NumDeath equal f_freact_creation[2]

fix		fix_print all print ${thermodump} "${step} ${etot} ${ke} ${peBond} ${peAngle} ${temp} ${press} ${StressX} ${StressY} ${StressZ} ${Fbond} ${NumType9} ${box_lengthX} ${box_Vol} ${CurrentConc} ${conc_preferred} ${P_add} ${P_remov} ${NumNuc} ${NumDeath} ${TermDen_add} ${TermDen_remov} ${TermId_add} ${TermId_remov} ${NumTot}" file ${thermofile} screen no title "step etot ke peBond peAngle temp press stressX stressY stressZ AvBondForce NumMons BoxLength BoxVol CurrentConc ConcPref P_add P_remov NreactNuc NreactDeath T_den_add T_den_remov T_id_add T_id_remov NumTot"

reset_atoms id sort no

run		${tbonds}
write_data	data.postmakebonds
write_restart	restart.postmakebonds

''')
    
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
#SBATCH --mem=500M 
#SBATCH --time=100:00:00 
#SBATCH --mail-user=ucapbbm@ucl.ac.uk
#SBATCH --mail-type=END
#SBATCH --export=NONE
#SBATCH --exclude=eta331

module load gcc
module load openmpi

mkdir dumplin dumplin_bonds dumplin_bonds_Sample restart
printf $RANDOM > lseed.dat
printf $RANDOM > vseed.dat

mpirun -np 1 /nfs/scistore26/saricgrp/bmeadowc/Scratch/lammps-15Jun2023/src/lmp_mpi -in collagen.in
\n''')
    f.close()

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
#SBATCH --mem=500M 
#SBATCH --time=100:00:00 
#SBATCH --mail-user=ucapbbm@ucl.ac.uk
#SBATCH --mail-type=END
#SBATCH --export=NONE
#SBATCH --exclude=eta331

module load gcc
module load openmpi

mkdir dumplin dumplin_bonds dumplin_bonds_Sample restart
printf $RANDOM > lseed.dat
printf $RANDOM > vseed.dat

mpirun -np 1 /nfs/scistore26/saricgrp/bmeadowc/Scratch/lammps-15Jun2023/src/lmp_mpi -in collagen.in
\n''')
    f.close()

##SBATCH --constraint="epsilon|delta|beta|leonid|serbyn|gamma"     
if __name__ == "__main__":
    main()

