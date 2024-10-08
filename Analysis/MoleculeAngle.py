# Boilerplate code generated by OVITO Pro 3.8.2
import sys
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *
import numpy as np

def angle(x1,y1,z1,x2,y2,z2):
    #dx, dy, dz = x1-x2, y1-y2, z1-z2
    mod1 = np.sqrt(x1*x1+y1*y1+z1*z1)
    mod2 = np.sqrt(x2*x2+y2*y2+z2*z2)
    dot = x1*x2 + y1*y2 + z1*z2
    angle = np.arccos(dot/(mod1*mod2))
    
    return angle*(180/np.pi)

def thetaphis(x1,y1,z1,x2,y2,z2):
    X, Y, Z = x1-x2, y1-y2, z1-z2
    R = np.sqrt(X*X+Y*Y+Z*Z)
    R2 = np.sqrt(X*X+Y*Y)
    phi = np.arccos(Z/R)*180/np.pi
    theta = np.arctan(Y/X)*180/np.pi
    #if Y <0:
    #    theta=360-theta
    return phi,theta

def DotX(x1,y1,z1,x2,y2,z2):
    X, Y, Z = x1-x2, y1-y2, z1-z2
    R = np.sqrt(X*X+Y*Y+Z*Z)
    dotX = X/R
    #if Y <0:
    #    theta=360-theta
    return dotX

numMolecules = int(sys.argv[3])
f = open(sys.argv[2],'a')    
filename = sys.argv[1]
pipeline = open(filename, 'r')
line = pipeline.readline()
nextline = pipeline.readline()
timestep = float(nextline.strip('\n'))
#print('timestep',timestep)
f.write('timestep '+str(timestep)+'\n')
for i in range(7):
    nextline = pipeline.readline()

for j in range(numMolecules):
    nextline = pipeline.readline()
    Particle1=nextline.split()
    X1,Y1,Z1 = float(Particle1[2]),float(Particle1[3]),float(Particle1[4])
    nextline = pipeline.readline()
    Particle2=nextline.split()
    X2,Y2,Z2 = float(Particle2[2]),float(Particle2[3]),float(Particle2[4])
    theta,phi = thetaphis(X1,Y1,Z1,X2,Y2,Z2)
    dotX =DotX(X1,Y1,Z1,X2,Y2,Z2)
    f.write(str(dotX)+'\n')
        #f.write(str(theta)+' '+str(phi)+'\n')
f.close()



