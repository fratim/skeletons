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
if(len(sys.argv))!=4:
    raise ValueError(" Scripts needs exactley 3 input arguments (bz by bx)")
else:
    bz = int(sys.argv[1])
    by = int(sys.argv[2])
    bx = int(sys.argv[3])

start_time_thinning = time.time()

output_folder = dataIO.OutputDirectory(prefix)
wiring.GenerateSkeleton(prefix, output_folder, bz, by, bx)

print("total time thinning:" + str(time.time()-start_time_thinning))
