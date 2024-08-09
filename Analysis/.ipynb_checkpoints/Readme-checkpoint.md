# Analysis scripts

All the simulations are analysed either using graph analysis, as in the python notebook, or using 'MoleculeAngle.py' or 'Histrogram_Types.py'.

The stress output and number of bonds broken/made is read directly from the lammps output thermo.dat.

AdjacencyMatrix.ipynb

<ins>Description</ins>

Contains an example of the use of the graph analysis. The second cell contains the graph analysis definitions, which take as inputs the bond topology files
    
<ins>Usage</ins>

As shown in the notebook, the inputs are the files within /dump_bonds/ folders from simulation runs. Various quantities are calculated using the functions in cell 2 like fraction of protomers in largest cluster, average number of bonds per monomer (average degree) etc.

## MoleculeAngle.py

<ins>Description</ins>

Takes as input the positions of protomers (in /dumplin/ folder output from simulations) and returns the dot protuct of the vector of that protomer with the x-axis
    
<ins>Usage</ins>

python MoleculeAngle.py ${filename} ${fileOut} ${NumMols}
Where filename is .../dumplin/lammpsT.000
fileout is name of the output file
NumMols is the number of protomers in your simulation
  
    
## Histrogram_Types.py 
    
<ins>Description</ins> 

Outputs a histogram of particle types from an input with multiple timesteps (uses the Ovito module, https://www.ovito.org/docs/current/python/)

<ins>Usage</ins> 

python Histrogram_Types.py ${filepattern0} ${fileOut}

filepattern0 is whole folder in which multiple timeteps are contained 