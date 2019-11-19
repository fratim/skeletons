# from synapseaware.isthmus import topological_thinning
# from synapseaware.teaser import teaser
from skeletons.connectome import wiring


prefix = 'Zebrafinch'
block_z = 1
block_y = 1
block_x = 1

# topological_thinning.TopologicalThinning(prefix, label)
# teaser.TEASER(prefix, label)
wiring.GenerateSkeleton(prefix, block_z, block_y, block_x)
# wiring.RefineSkeleton(prefix, label)
