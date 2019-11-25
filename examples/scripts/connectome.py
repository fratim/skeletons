# from synapseaware.isthmus import topological_thinning
# from synapseaware.teaser import teaser
from skeletons.connectome import wiring
import sys


if(len(sys.argv))!=4:
    raise ValueError(" Scripts needs exactley 3 input arguments (bz by bx)")
else:
    block_z = int(sys.argv[1])
    block_y = int(sys.argv[2])
    block_x = int(sys.argv[3])

prefix = 'Zebrafinch'

# topological_thinning.TopologicalThinning(prefix, label)
# teaser.TEASER(prefix, label)
# wiring.GenerateSkeleton(prefix, block_z, block_y, block_x)
wiring.RefineSkeleton(prefix, 149, block_z, block_y, block_x)
