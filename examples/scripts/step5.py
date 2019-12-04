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

start_time_refinement = time.time()
wiring.RefineSkeleton(  prefix, dataIO.OutputDirectory(prefix),start_blocks[0],start_blocks[1],start_blocks[2],
                        start_blocks[0]+n_blocks[0]-1,start_blocks[1]+n_blocks[1]-1,start_blocks[2]+n_blocks[2]-1)
print("total time refinement:" + str(time.time()-start_time_refinement))
