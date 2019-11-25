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

output_folder = dataIO.OutputDirectory(prefix)+'output-'+str(block_z).zfill(4)+'z-'+str(block_y).zfill(4)+'y-'+str(block_x).zfill(4)+'x'+'/'

#create output directories for this block
if not os.path.exists(output_folder): os.mkdir(output_folder)

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

for output_directory in output_directories:
    directory_path = output_folder+output_directory
    if not os.path.exists(directory_path): os.mkdir(directory_path)
    if not os.path.exists('{}/{}'.format(directory_path, prefix)): os.mkdir('{}/{}'.format(directory_path, prefix))

wiring.GenerateSkeleton(prefix, output_folder, block_z, block_y, block_x)
# wiring.RefineSkeleton(prefix, 149, block_z, block_y, block_x)
