from skeletons.connectome import wiring
import sys
from skeletons.utilities import dataIO
import os
import time

prefix = 'Zebrafinch'

start_blocks = dataIO.StartBlocks(prefix)
n_blocks = dataIO.NBlocks(prefix)

# pass arguments
if(len(sys.argv))!=3:
    raise ValueError(" Scripts needs exactley 2 input arguments (ID_start ID_end))")
else:
    ID_start = int(sys.argv[1])
    ID_end = int(sys.argv[2])

start_time_refinement = time.time()
wiring.RefineSkeleton(  prefix, dataIO.OutputDirectory(prefix),start_blocks[0],start_blocks[1],start_blocks[2],
                        start_blocks[0]+n_blocks[0]-1,start_blocks[1]+n_blocks[1]-1,start_blocks[2]+n_blocks[2]-1, ID_start, ID_end)

print("total time refinement:" + str(time.time()-start_time_refinement))
