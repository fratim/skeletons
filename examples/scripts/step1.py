from skeletons.connectome import wiring
import sys
from skeletons.utilities import dataIO
import os
import time

prefix = 'Zebrafinch'

output_directories = [  'synapses_projected/',
'anchorpoints_computed/',
'walls/',
'segment_IDs/',
'skeleton/',
'distances/',
'widths/',
'running_times/',
'running_times/refinement/',
'running_times/skeleton/']

start_blocks = dataIO.StartBlocks(prefix)
n_blocks = dataIO.NBlocks(prefix)

# pass arguments
if(len(sys.argv))!=2:
    raise ValueError(" Scripts needs exactley 1 input argument (bz)")
else:
    bz = int(sys.argv[1])

for by in range(start_blocks[1], start_blocks[1] + n_blocks[1]):
    for bx in range(start_blocks[2], start_blocks[2] + n_blocks[2]):

        output_folder = dataIO.OutputDirectory(prefix)+'output-'+str(bz).zfill(4)+'z-'+str(by).zfill(4)+'y-'+str(bx).zfill(4)+'x'+'/'
        wiring.SaveWalls(prefix, output_folder, bz, by, bx)
