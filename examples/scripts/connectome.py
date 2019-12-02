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

# for bz in range(start_blocks[0], start_blocks[0]+n_blocks[0]):
#     for by in range(start_blocks[2], start_blocks[1]+n_blocks[1]):
#         for bx in range(start_blocks[2], start_blocks[2]+n_blocks[2]):
#             #create output directories for this block
#             output_folder = dataIO.OutputDirectory(prefix)+'output-'+str(bz).zfill(4)+'z-'+str(by).zfill(4)+'y-'+str(bx).zfill(4)+'x'+'/'
#             if not os.path.exists(output_folder): os.mkdir(output_folder)
#             for output_directory in output_directories:
#                 directory_path = output_folder+output_directory
#                 if not os.path.exists(directory_path): os.mkdir(directory_path)
#                 if not os.path.exists('{}/{}'.format(directory_path, prefix)): os.mkdir('{}/{}'.format(directory_path, prefix))
#
# start_time_savewalls = time.time()
# for bz in range(start_blocks[0], start_blocks[0]+n_blocks[0]):
#     for by in range(start_blocks[2], start_blocks[1]+n_blocks[1]):
#         for bx in range(start_blocks[2], start_blocks[2]+n_blocks[2]):
#             output_folder = dataIO.OutputDirectory(prefix)+'output-'+str(bz).zfill(4)+'z-'+str(by).zfill(4)+'y-'+str(bx).zfill(4)+'x'+'/'
#             wiring.SaveWalls(prefix, output_folder, bz, by, bx)
# print("total time save walls:" + str(time.time()-start_time_savewalls))
#
# start_time_anchors = time.time()
# for bz in range(start_blocks[0], start_blocks[0]+n_blocks[0]):
#     for by in range(start_blocks[2], start_blocks[1]+n_blocks[1]):
#         for bx in range(start_blocks[2], start_blocks[2]+n_blocks[2]):
#             output_folder = dataIO.OutputDirectory(prefix)
#             wiring.MakeAnchorpoints(prefix, output_folder, bz, by, bx)
# print("total time make anchorpoints:" + str(time.time()-start_time_anchors))
#
# start_time_thinning = time.time()
# for bz in range(start_blocks[0], start_blocks[0]+n_blocks[0]):
#     for by in range(start_blocks[2], start_blocks[1]+n_blocks[1]):
#         for bx in range(start_blocks[2], start_blocks[2]+n_blocks[2]):
#             output_folder = dataIO.OutputDirectory(prefix)
#             wiring.GenerateSkeleton(prefix, output_folder, bz, by, bx)
# print("total time thinning:" + str(time.time()-start_time_thinning))

start_time_refinement = time.time()
wiring.RefineSkeleton(  prefix, dataIO.OutputDirectory(prefix),start_blocks[0],start_blocks[1],start_blocks[2],
                        start_blocks[0]+n_blocks[0]-1,start_blocks[1]+n_blocks[1]-1,start_blocks[2]+n_blocks[2]-1)
print("total time refinement:" + str(time.time()-start_time_refinement))
