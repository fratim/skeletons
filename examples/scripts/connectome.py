# from synapseaware.isthmus import topological_thinning
# from synapseaware.teaser import teaser
# from synapseaware.connectome 
import wiring




prefix = 'Zebrafinch'
label = 55

# topological_thinning.TopologicalThinning(prefix, label)
# teaser.TEASER(prefix, label)
wiring.GenerateSkeleton(prefix, label)
# wiring.RefineSkeleton(prefix, label)
