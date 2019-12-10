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
'somae_surface',
'distances/',
'widths/',
'running_times/',
'running_times/refinement/',
'running_times/skeleton/']

start_blocks = dataIO.StartBlocks(prefix)
n_blocks = dataIO.NBlocks(prefix)

output_folder = dataIO.OutputDirectory(prefix)
if not os.path.exists(output_folder): os.mkdir(output_folder)

output_folder_skeletons = dataIO.OutputDirectory(prefix)+'skeletons/'
if not os.path.exists(output_folder_skeletons): os.mkdir(output_folder_skeletons)
if not os.path.exists('{}/{}'.format(output_folder_skeletons, prefix)): os.mkdir('{}/{}'.format(output_folder_skeletons, prefix))

output_folder_times = dataIO.OutputDirectory(prefix)+'running_times/'
if not os.path.exists(output_folder_times): os.mkdir(output_folder_times)
if not os.path.exists('{}/{}'.format(output_folder_times, prefix)): os.mkdir('{}/{}'.format(output_folder_times, prefix))

for bz in range(start_blocks[0], start_blocks[0]+n_blocks[0]):
    for by in range(start_blocks[2], start_blocks[1]+n_blocks[1]):
        for bx in range(start_blocks[2], start_blocks[2]+n_blocks[2]):
            #create output directories for this block
            output_folder = dataIO.OutputDirectory(prefix)+'output-'+str(bz).zfill(4)+'z-'+str(by).zfill(4)+'y-'+str(bx).zfill(4)+'x'+'/'
            if not os.path.exists(output_folder): os.mkdir(output_folder)
            for output_directory in output_directories:
                directory_path = output_folder+output_directory
                if not os.path.exists(directory_path): os.mkdir(directory_path)
                if not os.path.exists('{}/{}'.format(directory_path, prefix)): os.mkdir('{}/{}'.format(directory_path, prefix))
