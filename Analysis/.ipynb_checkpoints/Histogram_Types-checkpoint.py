import sys
from ovito.io import *
from ovito.modifiers import *
from ovito.io import export_file
from ovito.data import *
import numpy as np


def main(): 
    bins  = 500
    pipeline = import_file(sys.argv[1]+'/dumplin/dump.run_T1.*.lammpstrj', multiple_frames = True, columns = ["Particle Type","Molecule Identifier", "Position.X","Position.Y","Position.Z"])
    
    # Wrap at periodic boundaries:
    pipeline.modifiers.append(WrapPeriodicImagesModifier())

    modifier2 = HistogramModifier(bin_count=7, property = 'Particle Type',fix_xrange = True, xrange_start = 0.5, xrange_end = 7.5) 
    pipeline.modifiers.append(modifier2)
    data = pipeline.compute()
    export_file(pipeline, str(sys.argv[2]), "txt/table", key="histogram[Particle Type]",multiple_frames = True)

if __name__ == "__main__":
    main()
