This repository enables topological thinning of connectomics data, guarateeing a one-to-open correspondence of skeleton endpoints and synapses. 
Input:
 1. Segmentation volume saved in blocks (h5 files)
 2. Synapse locations as 3D coordinates (point list)
(3.) somae segmentation in blocks (h5 file)

It operates in a MapReduce style, composed of four different steps:
1. At first, output folders are created and the boundary slices in z,y,x direction are saved for each block contained in the volume
2. Anchorpoints are computed for each boundary slice, then saved to the corresponding output folders
3. Skeletons are generated for each block, using 3D topological thinning. Synapses and anchor-points are set to be fix points, hence are part of the output skeleton.
4. On a global scale, skeletons for each neuron are refined.
