import os
import sys


prefix = sys.argv[1]

if not os.path.exists('output'): os.mkdir('output')

output_directories = [  'synapses_projected',
                        'anchorpoints_computed',
                        'anchorpoints_seeded',
                        'connectomes',
                        'IDs_processed',
                        'IDs_present',
                        'skeleton',
                        'distances',
                        'widths',
                        'running_times',
                        'running_times/refinement',
                        'running_times/skeleton']


for output_directory in output_directories:
    directory_path = 'output/'+output_directory
    if not os.path.exists(directory_path): os.mkdir(directory_path)
    if not os.path.exists('{}/{}'.format(directory_path, prefix)): os.mkdir('{}/{}'.format(directory_path, prefix))
