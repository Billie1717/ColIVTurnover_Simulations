# This folder contains a set of different simulation protocols (with the same underlying model) which are used to run the simulations

The simulations are then analysed using the scripts in /Analysis folder and finally plotting using the plotting in /Figs_and_Data folder

## There are 4 main simulation protocols used in this work:

 - Assembly (with various bond energies)
 - Stretching in the x-direction and monitoring relaxation (different seeds from initially unbonded atoms)
 - Non-equilibrium bond remodelling
 - Stretching in x- and y- directions in the non-equilibrium remodelling regime, allowing some relaxation without bond remodelling and then full relaxation (see /Data_and_Figures/RicardoInVivoTurnover/Fig_Simprot/). Monomers are added after stretch to keep free monomer concentration constant.
 
Below, the python script for creating input files for the four protocols are outlined.


## build_Assembly.py:

Created simulations in the remote directory :
Scratch/Collagen/model_strict_tetramer/EquilibriumSeededWorkflow/runsTestStress/

Called Scratch/Collagen/model_strict_tetramer/EquilibriumSeededWorkflow/build_SeededPrintingBondTurnover.py in remote directory

Analysed directly using jupyter notebook with graph analysis on data /dump_bonds/ (remote notebook Scratch/Collagen/model_strict_tetramer/EquilibriumSeededWorkflow/Plotting/AdjacencyMatrix.ipynb)
 

## build_StretchEquilib.py:

Created simulations in the remote directory :
Scratch/Collagen/model_strict_tetramer/EquilibriumSeededWorkflow/runsTauPhiSweep/

Called Scratch/Collagen/model_strict_tetramer/EquilibriumSeededWorkflow/build_SeededEquilib.py in remote directory

Analysed using a mixture of graph analysis from jupyter notebooks and MoleculeAlignment.py which outputs the alignment of molecules with the x-axis
The two important plotting scripts for this data are :
1) Scratch/Collagen/model_strict_tetramer/EquilibriumSeededWorkflow/Plotting/MasterMechanicsCB5.0_seeds.ipynb
2) Scratch/Collagen/model_strict_tetramer/EquilibriumSeededWorkflow/Plotting/AlignmentAndFittingNEWEQ.ipynb


## build_StretchNonEquilib.py

Created simulations in the remote directory :
Scratch/Collagen/model_strict_tetramer/NonEquilibriumSeededWorkflow/

Called Scratch/Collagen/model_strict_tetramer/NonEquilibriumSeededWorkflow/build_SeededNonEquilib.py

Analysed using a mixture of graph analysis from jupyter notebooks and MoleculeAlignment.py which outputs the alignment of molecules with the x-axis
The important plotting scripts for this data are :
Scratch/Collagen/model_strict_tetramer/NonEquilibriumSeededWorkflow/Plotting/AlignmentAndFitting.ipynb


## /slab_stretchXY/

The workflow for this protocol is as follows:

1. Equilibrate a slab with a given monomer turnover/ Break & making bonds probability (runMakeScripts.sh --> build_Equilibrate.py)
2. Stretch this network some amount (runStretch.sh --> build_Stretch.py)
3. Add a number of monomers st freely floating monomer number is kept constant --> AddMonomersToDataFile.py
4. Allow network to relax with given break & making bonds probability from before (runRelax.sh --> build_Relax.py)
5. Measure stress relaxation

Where all the mentioned scripts are in the folder /slab_stretchXY/


