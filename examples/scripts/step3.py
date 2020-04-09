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
if(len(sys.argv)!=3):
    raise ValueError(" Scripts needs exactley 2 inputs (bz output_file_step2_bz) ")
else:
    bz = int(sys.argv[1])
    out_file_S1 = sys.argv[2]

inp_file = open(out_file_S1)
inp_text = inp_file.read()
inp_file.close()

if inp_text[0]!="0":
    print(inp_text)
    raise ValueError("Execution Stopped: Wrong Error Code (!=0)")

for by in range(start_blocks[1], start_blocks[1] + n_blocks[1]):
    for bx in range(start_blocks[2], start_blocks[2] + n_blocks[2]):

    #start_time_thinning = time.time()

    output_folder = dataIO.OutputDirectory(prefix)
    wiring.GenerateSkeleton(prefix, output_folder, bz, by, bx)

    #print("total time thinning:" + str(time.time()-start_time_thinning))

print("0")
