From the lammps downloads:

https://download.lammps.org/tars/

The version used in this work is: lammps-15Jun2023

With some packages edited, namely the files in the REACTION package and some files in the MC package (see below)

all simulations were run on Lammps built with the following packages :

make yes-dipole
make yes-extra-pair
make yes-rigid
make yes-molecule
make yes-reaction
make yes-MC

# Files edited in folder /src

For making and breaking bonds:

MC/fix_bond_break.cpp
MC/fix_bond_break.h
MC/fix_bond_create.cpp
MC/fix_bond_create.h

For bug fixes and other things in reaction package:

REACTION/fix_bond_react.cpp
REACTION/fix_bond_react.h

The model only works with these edited files, they have to replace the original files in the lammps src folder and lammps has to be built with these files.
