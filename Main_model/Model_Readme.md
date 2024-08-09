# Explaining the model

The model is run in lammps and uses the REACTION package. The reactions allow the changing of atom types depending on their bonding, allowing us to choose the sequence of bonding.

## atom-type changing templates

To change types depending on bonding, it is required to specify many atom templates which singles out atoms bonded in a certain way and changes their type accordingly. These are contained in the folder /templates. This folder is essential and must be kept in the same directory as the lammps input file for the simulations to run.

## initial conditions

The files lseed.dat and vseed.dat contain random numbers which set the initial random lengevin thermostat and initial random velocity, respectivly

## inital data file

The file 'data' contains a grid of non-bonded protomers all of type 1 and 2 in a simulation box 72x72x12. This example has 3125 protomers, 6250 atoms. 

## input file

the file input_example.in is not complete for running a simulation but contains all the ingredients of the input file included comments, it skips the many 'cycles' of runs explained below.

the file input.in is complete and should be able to be run with the correct lammps version + src code edits (Note : it is a very long lammps input file!)

Because the command fix create and fix break for bonds in lammps shouldn't be used simultaneously, the main point to note is that the code requires bonds to be make and broken separately, with type-IDs of atoms being updated in between. These cycles have to be short enough so that atoms types can be updated sufficiently frequently without a single atom bonding more than once. This makes the lammps input script very long as it may have 100s of cycles of breaking and making bonds. 

