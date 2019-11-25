from skeletons.connectome import wiring
import sys
from skeletons.utilities import dataIO
import os


if(len(sys.argv))!=4:
    raise ValueError(" Scripts needs exactley 3 input arguments (bz by bx)")
else:
    block_z = int(sys.argv[1])
    block_y = int(sys.argv[2])
    block_x = int(sys.argv[3])

prefix = 'Zebrafinch'

output_directories = [  'synapses_projected/',
'anchorpoints_computed/',
'anchorpoints_seeded/',
'connectomes/',
'IDs_processed/',
'IDs_present/',
'skeleton/',
'distances/',
'widths/',
'running_times/',
'running_times/refinement/',
'running_times/skeleton/']

# for bz in range(block_z+1):
#     for by in range(block_y+1):
#         for bx in range(block_x+1):
#
#             #create output directories for this block
#             output_folder = dataIO.OutputDirectory(prefix)+'output-'+str(bz).zfill(4)+'z-'+str(by).zfill(4)+'y-'+str(bx).zfill(4)+'x'+'/'
#             if not os.path.exists(output_folder): os.mkdir(output_folder)
#
#             for output_directory in output_directories:
#                 directory_path = output_folder+output_directory
#                 if not os.path.exists(directory_path): os.mkdir(directory_path)
#                 if not os.path.exists('{}/{}'.format(directory_path, prefix)): os.mkdir('{}/{}'.format(directory_path, prefix))
#
#             wiring.GenerateSkeleton(prefix, output_folder, bz, by, bx)

wiring.RefineSkeleton(prefix, dataIO.OutputDirectory(prefix),149,0,0,0,1,1,1)
