from __future__ import division, print_function
import numpy as np
#import matplotlib.pyplot as plt
import math
import time
import random
import sys, getopt
import os
import pandas as pd

 #str(sys.argv[1])


def MakePositions(xStretch,OrigBox,zBox,nmols):
    x=[]
    y=[]
    z=[]
    HalfStretch = xStretch/2
    FullNewX = OrigBox+xStretch
    #print(rng1,rng2)
    vx,vy,vz = [],[],[]
    for n in range(nmols):
        r1,r2,r3,phiR,thetaR = np.random.default_rng().random(),np.random.default_rng().random(),np.random.default_rng().random(),np.random.default_rng().random(),np.random.default_rng().random()
        vxR,vyR,vzR =np.random.default_rng().random(),np.random.default_rng().random(),np.random.default_rng().random()
        XX = r1*FullNewX-HalfStretch
        x.append(XX)
        x.append(XX+3.*np.sin(thetaR*2*np.pi)*np.cos(phiR*2*np.pi))
        YY = r2*FullNewX-HalfStretch
        y.append(YY)
        y.append(YY+3.*np.sin(thetaR*2*np.pi)*np.sin(phiR*2*np.pi))
        ZZ = r3*zBox
        z.append(ZZ)
        z.append(ZZ+3.*np.cos(thetaR*2*np.pi))

        vx.append(np.add(vxR,-0.5))
        vx.append(np.add(vxR,-0.5))
        vy.append(np.add(vyR,-0.5))
        vy.append(np.add(vyR,-0.5))
        vz.append(np.add(vzR,-0.5))
        vz.append(np.add(vzR,-0.5))
    return x, y, z, vx, vy, vz

# Define the file paths
input_file_path = str(sys.argv[1])
OriginalNumMols = int(sys.argv[2])
OrigXbox=float(sys.argv[3])
Xstretch=float(sys.argv[4])

OldArea = OrigXbox*OrigXbox
NewArea=(OrigXbox+Xstretch)*(OrigXbox+Xstretch)
Ratio = NewArea/OldArea
NewNumMols = Ratio*OriginalNumMols
NumMolsAdded = int(NewNumMols-OriginalNumMols)
NumMonsAdded = int(NumMolsAdded*2)

output_file_path = input_file_path+'addedMons'
OrigNumMons=6330
OrigNumMols=int(OrigNumMons/2)
TotalMons = str(NumMonsAdded+OrigNumMons)

x, y, z, vx, vy, vz = MakePositions(Xstretch,OrigXbox,12,NumMolsAdded)

# Step 1: Read the file line by line and store the lines
with open(input_file_path, 'r') as file:
    lines = file.readlines()

# Step 2: Initialize an empty list to store modified lines
modified_lines = []
counter = 0
# Step 3: Process each line
for line in lines:
    counter+=1
    # Check if the line starts with "XXX"
    
    if counter==3:
        # Modify the line as needed, e.g., append "_edited" to the end of the line
        modified_line = TotalMons + ' atoms\n'
        modified_lines.append(modified_line)
    elif counter==5:
        numbonds = line.split()[0]
        print(line.split()[0])
        modified_line = str(int(numbonds)+NumMolsAdded) + ' bonds\n'
        modified_lines.append(modified_line)
    if counter!=3 and counter!=5:
        modified_line = line
        modified_lines.append(modified_line)
    if line.startswith("Velocities"): #Appending more atoms on the end of atom data
        print(modified_lines[counter-1])
        nmon=0
        for nmol in range(NumMolsAdded):
            extra_line = str(OrigNumMons+1+nmon)+" "+str(OrigNumMols+1+nmol) +" 1 "+str(x[nmon])+" "+str(y[nmon])+" "+str(z[nmon]) + " 0 0 0\n"
            if nmol==0:
                modified_lines[counter-2] =extra_line
            else:
                modified_lines.append(extra_line)
            nmon+=1
            extra_line = str(OrigNumMons+1+nmon)+" "+str(OrigNumMols+1+nmol) +" 2 "+str(x[nmon]+3)+" "+str(y[nmon])+" "+str(z[nmon]) + " 0 0 0\n" #for now, all moleules in the x-direction?
            if nmol==0:
                modified_lines[counter-1] =extra_line
            else:
                modified_lines.append(extra_line)
            counter+=2
            nmon+=1
        modified_lines.append("\n")
        modified_lines.append("Velocities\n")
        
    if line.startswith("Bonds"): #Appending velocities
        for nmon in range(NumMonsAdded):
            extra_line = str(OrigNumMons+1+nmon)+" "+str(vx[nmon])+" "+str(vy[nmon])+" "+str(vz[nmon]) + "\n"
            if nmon==0 or nmon==1:
                modified_lines[counter-2] =extra_line
            else:
                modified_lines.append(extra_line)
            counter+=1
        modified_lines.append("\n")
        modified_lines.append("Bonds\n")
        
    if line.startswith("Angles"): #Appending bonds
        StartNum=modified_lines[counter-3].split()[0]
        print(StartNum)
        nmon=0
        for nmol in range(NumMolsAdded):
            extra_line = str(int(StartNum)+1+nmol)+" 1 "+str(OrigNumMons+1+nmon)+" "+str(OrigNumMons+2+nmon) + "\n"               
            nmon+=2
            if nmol==0 or nmol==1:
                modified_lines[counter-2] =extra_line
            else:
                modified_lines.append(extra_line)
            counter+=1
        modified_lines.append("\n")
        modified_lines.append("Angles\n")
    if line.startswith("bond_react_props_internal"):
        break
#for nmon in range(NumMonsAdded):
#    extra_line = str(6251+nmon)+" 0 0\n"
#    modified_lines.append(extra_line)
    


# Step 4: Combine the original lines and the modified lines
#all_lines = lines + modified_lines

# Step 5: Write all lines to a new file
with open(output_file_path, 'w') as file:
    file.writelines(modified_lines)

### ADD BOTTOM LINES XXX 0 0
### ADD POSITIONS, VELOCITIES, BONDS