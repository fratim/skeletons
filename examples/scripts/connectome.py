from skeletons.connectome import wiring
import sys
from skeletons.utilities import dataIO
import os
import time


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
'walls/',
'segment_IDs/',
'skeleton/',
'distances/',
'widths/',
'running_times/',
'running_times/refinement/',
'running_times/skeleton/']

# for bz in range(block_z+1+2):
#     for by in range(block_y+1+2):
#         for bx in range(block_x+1+2):
#             #create output directories for this block
#             output_folder = dataIO.OutputDirectory(prefix)+'output-'+str(bz).zfill(4)+'z-'+str(by).zfill(4)+'y-'+str(bx).zfill(4)+'x'+'/'
#             if not os.path.exists(output_folder): os.mkdir(output_folder)
#             for output_directory in output_directories:
#                 directory_path = output_folder+output_directory
#                 if not os.path.exists(directory_path): os.mkdir(directory_path)
#                 if not os.path.exists('{}/{}'.format(directory_path, prefix)): os.mkdir('{}/{}'.format(directory_path, prefix))
#
# start_time_savewalls = time.time()
# for bz in range(block_z+1+1):
#     for by in range(block_y+1+1):
#         for bx in range(block_x+1+1):
#             output_folder = dataIO.OutputDirectory(prefix)+'output-'+str(bz).zfill(4)+'z-'+str(by).zfill(4)+'y-'+str(bx).zfill(4)+'x'+'/'
#             wiring.SaveWalls(prefix, output_folder, bz, by, bx)
# print("total time save walls:" + str(time.time()-start_time_savewalls))
#
# for bz in range(block_z+1):
#     for by in range(block_y+1):
#         for bx in range(block_x+1):
#             output_folder = dataIO.OutputDirectory(prefix)
#             wiring.MakeAnchorpoints(prefix, output_folder, bz, by, bx)
#
# start_time_thinning = time.time()
#
# for bz in range(block_z+1):
#     for by in range(block_y+1):
#         for bx in range(block_x+1):
#             output_folder = dataIO.OutputDirectory(prefix)+'output-'+str(bz).zfill(4)+'z-'+str(by).zfill(4)+'y-'+str(bx).zfill(4)+'x'+'/'
#             wiring.GenerateSkeleton(prefix, output_folder, bz, by, bx)
#
# print("total time thinning:" + str(time.time()-start_time_thinning))


start_time_refinement = time.time()

wiring.RefineSkeleton(prefix, dataIO.OutputDirectory(prefix),0,0,0,1,1,1)

print("total time refinement:" + str(time.time()-start_time_refinement))
