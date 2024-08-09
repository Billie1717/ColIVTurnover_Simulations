# Explaining the model

The model is run in lammps and uses the REACTION package. The reactions allow the changing of atom types depending on their bonding, allowing us to choose the sequence of bonding.

## atom-type changing templates

To change types depending on bonding, it is required to specify many atom templates which singles out atoms bonded in a certain way and changes their type accordingly. These are contained in the folder /templates. This folder is essential and must be kept in the same directory as the lammps input file (and lammps initial data file) for the simulations to run.

## initial conditions

The files lseed.dat and vseed.dat contain random numbers which set the initial random lengevin thermostat and initial random velocity, respectivly.

## inital data file

The file 'data' contains the lammps data file which for these simulations is an equilibrated bonded network in a simulation box 51x51x13. This example has 3165 protomers, 6330 atoms. 

## input file

The file input_example.in is not complete for running a simulation but contains all the ingredients of the input file including comments, it skips the many 'cycles' of runs explained below.

The file input.in is complete and should be able to be run with the correct lammps version + src code edits.

Because the command fix create and fix break for bonds in lammps shouldn't be used simultaneously, the main point to note is that the code requires bonds to be make and broken separately, with type-IDs of atoms being updated in between. These cycles have to be short enough so that atoms types can be updated sufficiently frequently without a single atom bonding more than once. This makes the lammps input script very long as it may have 100s of cycles of breaking and making bonds. 


## Outputs

The turnover is always measured by sampling the network topoligy from the bonds. Therefore, the main other useful output from lammps is the stress in the x- and y- directions, so that the stress-relaxation can be measured. The following commands can be found in the lammps input script:

    compute totalStressX all reduce sum c_perAtomStress[1]
    compute totalStressY all reduce sum c_perAtomStress[2]

    variable StressX equal c_totalStressX
    variable StressY equal c_totalStressY

    fix		fix_print all print ${thermodump} "${step} ${etot} ${ke} ${peBond} ${peAngle} ${temp} ${press} ${StressX} ${StressY} ${StressZ} ${Fbond}" file ${thermofile} screen no title "step etot ke peBond peAngle temp press stressX stressY stressZ AvBondForce"

So that the xx and yy components of the stress tensor are outputted in 'thermo.dat' which we read out for each time point in each simulation.

