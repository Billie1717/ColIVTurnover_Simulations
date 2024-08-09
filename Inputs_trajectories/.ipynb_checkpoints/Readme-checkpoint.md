# This folder contains a set of different simulation protocols (with the same underlying model) which are used to run the simulations

The simulations are then analysed using the scripts in /Analysis folder and finally plotting using the plotting in /Figs_and_Data folder


## /slab_stretchXY/

The workflow for this protocol is as follows:

1. Equilibrate a slab with a given monomer turnover/ Break & making bonds probability (runMakeScripts.sh --> build_Equilibrate.py)
2. Stretch this network some amount (runStretch.sh --> build_Stretch.py)
3. Add a number of monomers st freely floating monomer number is kept constant --> AddMonomersToDataFile.py
4. Allow network to relax with given break & making bonds probability from before (runRelax.sh --> build_Relax.py)
5. Measure stress relaxation

Where all the mentioned scripts are in the folder /slab_stretchXY/


