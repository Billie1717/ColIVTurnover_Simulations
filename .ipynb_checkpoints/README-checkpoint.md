## CollagenIV turnover simulations

This repository contains all the code required to create simulations of collagenIV turnover and relaxation. 


<ins>Simulation overview:</ins>

CollagenIV molecules are modelled as two beads connected by a stiff spring where the beads are reactive. One bead is modelling the 7s end of the collagenIV molecule which can react with three other 7s beads to form a tetramer bond. The other is the NC1 end of the molecule which can react with one other NC1 bead to form a dimer bond. The bonding reactions of the 7S are via a set sequence of dimer-trimer-tetramer where the end result is 4 beads all connected (by a total of 6 bonds). 

Many protomers are simulated in a molecular dynamics simulation using LAMMPS (https://www.lammps.org/). There are 3 stages of simulating the network turnover and relaxation. Firstly the molecules are left to form a bonded network. This bonded network is dynamically remodelling with the breaking and making of bonds. In dynamic equilibrium, the turnover (average amount of time a molecule resides within the network) is measured by frequent sampling of the network topology. The network is then stretched in the x- and y- directions. Finally the network relaxes at this fixed strain while the stress in the x- and y- directions is measured. 


## Folders contained in this repository:

<ins>Lammps:</ins>

Building the correct Lammps version and replacing certain cpp source codes with the edited versions

<ins>Main_model:</ins>

The files required to run the simulations and a simplified commented input script. There is also a description in 'readme' of the useful outputs that this lammps script will produce.

<ins>Inputs_trajectories:</ins>

The python and bash scripts for creating the input scripts for: the equilibium, stretching, monomer addition and relaxation stages of the simulation protocol.
A full list of parameters used in each stage is documented here in the file 'Parameters.md'

<ins>Analysis:</ins>

Analysis scripts used in this work (graph analysis and molecule-angle analysis).

<ins>Data_and_Figures:</ins>

The plotting scripts and basic data used in the main figures of this work

