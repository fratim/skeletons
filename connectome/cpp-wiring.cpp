/* c++ file for mutual functions and variables */

#include <string.h>
#include "cpp-wiring.h"
#include <iostream>
#include <stdexcept>



// set default values

long padded_blocksize[3] = { -1, -1, -1 };
long input_blocksize[3] = { -1, -1, -1 };
long volumesize[3] = { -1, -1, -1 };
float resolution[3] = { -1, -1, -1 };
long block_z = -1;
long block_y = -1;
long block_x = -1;
long nentries = -1;
long sheet_size = -1;
long row_size = -1;
long infinity = -1;

// create new default variables (will be overwritten) TODO> do I really need these declarations?
std::unordered_map<long, char> segment = std::unordered_map<long, char>();
std::unordered_set<long> synapses = std::unordered_set<long>();
std::unordered_set<long> IDs_in_block = std::unordered_set<long>();


void CppUpdateResolution(float input_resolution[3])
{
    resolution[OR_Z] = input_resolution[OR_Z];
    resolution[OR_Y] = input_resolution[OR_Y];
    resolution[OR_X] = input_resolution[OR_X];
}

void CppUpdateBlocksize(long inp_blocksize[3])
{
    input_blocksize[OR_Z] = inp_blocksize[OR_Z];
    input_blocksize[OR_Y] = inp_blocksize[OR_Y];
    input_blocksize[OR_X] = inp_blocksize[OR_X];

    padded_blocksize[OR_Z] = inp_blocksize[OR_Z]+2;
    padded_blocksize[OR_Y] = inp_blocksize[OR_Y]+2;
    padded_blocksize[OR_X] = inp_blocksize[OR_X]+2;

}

void CppUpdateVolumesize(long volume_size[3])
{
    volumesize[OR_Z] = volume_size[OR_Z];
    volumesize[OR_Y] = volume_size[OR_Y];
    volumesize[OR_X] = volume_size[OR_X];

}

void CppUpdateBlockindices(long inp_block_z, long inp_block_y, long inp_block_x)
{
    block_z = inp_block_z;
    block_y = inp_block_y;
    block_x = inp_block_x;
}

///////////////////////////////////////
//// POINT CLOUD UTILITY FUNCTIONS ////
///////////////////////////////////////

/* conventient I/O function */
void CppPopulatePointCloudFromH5(long label, long *labels) {

    std::unordered_map<long, std::unordered_map<long,char>> Pointclouds;
    std::unordered_set<long> synapses;
    std::unordered_set<long> IDs_in_block;

    // indexing parameters for indexing within current block
    nentries = padded_blocksize[OR_Z] * padded_blocksize[OR_Y] * padded_blocksize[OR_X];
    sheet_size = padded_blocksize[OR_Y] * padded_blocksize[OR_X];
    row_size = padded_blocksize[OR_X];
    infinity = padded_blocksize[OR_Z] * padded_blocksize[OR_Z] + padded_blocksize[OR_Y] * padded_blocksize[OR_Y] + padded_blocksize[OR_X] * padded_blocksize[OR_X];

    long n_points = input_blocksize[0]*input_blocksize[1]*input_blocksize[2];

    std::cout << "n_points is: " << n_points << std::endl << std::flush;

    for (long voxel_index = 0; voxel_index < n_points; voxel_index++){

      // get label of current index and skip if is zero
      long curr_label = labels[voxel_index];
      if (!curr_label) continue;

      // find coordinates of this voxel_index
      long iz = voxel_index / (input_blocksize[OR_Y] * input_blocksize[OR_X]);
      long iy = (voxel_index - iz * (input_blocksize[OR_Y] * input_blocksize[OR_X])) / input_blocksize[OR_X];
      long ix = voxel_index % input_blocksize[OR_X];

      //  pad the location by one
      iz += 1; iy += 1; ix += 1;

      // find the new voxel index
      long iv = IndicesToIndex(ix, iy, iz);

      // check if pointcloud of this label already exists, otherwise add new pointcloud
      if (Pointclouds.find(curr_label) == Pointclouds.end()) {
        Pointclouds[curr_label] = std::unordered_map<long,char>();
        // std::cout << "New label detected: " << curr_label << std::endl << std::flush;
        IDs_in_block.insert(curr_label);
      }

      // add index to pointcloud
      Pointclouds[curr_label][iv] = 1;

    }

    // save segment as the pointcloud for a specific label TODO: remove this and return array of labels to be iterated over
    segment = Pointclouds[label];

    std::cout << "Printing all IDs in this block" <<std::endl<<std::flush;

    // iterate over set to print all components
    std::unordered_set<long>::iterator it = IDs_in_block.begin();
    while (it != IDs_in_block.end())
    {
      std::cout << *it << ",";
    	it++;
    }

    std::cout << std::endl << "Printed all IDS and created map of pointclouds" << std::endl << std::flush;


}

/* conventient I/O function */
void CppPopulatePointCloud(const char *prefix, const char *dataset, long label) {

    std::cout << "Next: Populate Point Cloud" << std::endl << std::flush;
    // read in the point cloud for this label
    char filename[4096];
    sprintf(filename, "%s/%s/%06ld.pts", dataset, prefix, label);
    // sprintf(filename, "/home/frtim/Documents/Code/skeletons/examples/synapses/Zebrafinch/syn_0055.txt", dataset, prefix, label);

    FILE *fp = fopen(filename, "rb");
    if (!fp) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }

    long read_volumesize[3];
    long npoints;
    if (fread(&(read_volumesize[OR_Z]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }
    if (fread(&(read_volumesize[OR_Y]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }
    if (fread(&(read_volumesize[OR_X]), sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }
    if (fread(&npoints, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }

    //
    std::cout<<"Volume size read in: "<<read_volumesize[0]<<","<<read_volumesize[1]<<","<<read_volumesize[2]<<","<<std::endl<<std::flush;
    std::cout<<"npoints: "<<npoints<<std::endl<<std::flush;

    if (read_volumesize[0]!=volumesize[0] || read_volumesize[1]!=volumesize[1] || read_volumesize[2]!=volumesize[2]) {
        throw std::invalid_argument("read_volumesize not equal to volumesize");
    }

    // set block indexing parameters
    nentries = padded_blocksize[OR_Z] * padded_blocksize[OR_Y] * padded_blocksize[OR_X];
    sheet_size = padded_blocksize[OR_Y] * padded_blocksize[OR_X];
    row_size = padded_blocksize[OR_X];
    infinity = padded_blocksize[OR_Z] * padded_blocksize[OR_Z] + padded_blocksize[OR_Y] * padded_blocksize[OR_Y] + padded_blocksize[OR_X] * padded_blocksize[OR_X];

    for (long ip = 0; ip < npoints; ++ip) {
        long voxel_index;
        if (fread(&voxel_index, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", filename); exit(-1); }

        long iz = voxel_index / (volumesize[OR_Y] * volumesize[OR_X]);
        long iy = (voxel_index - iz * (volumesize[OR_Y] * volumesize[OR_X])) / volumesize[OR_X];
        long ix = voxel_index % volumesize[OR_X];

        // offset indices according to current block position
        iz -= block_z * input_blocksize[OR_Z];
        iy -= block_y * input_blocksize[OR_Y];
        ix -= block_x * input_blocksize[OR_X];

        // skip point if not within block
        if (iz<0 || iy<0 || ix<0 || iz>=input_blocksize[OR_Z] ||  iy>=input_blocksize[OR_Y] || ix>=input_blocksize[OR_X]) continue;

        //  pad the location by one
        iz += 1; iy += 1; ix += 1;

        std::cout << "synapse added at: "<<iz<<","<<iy<<","<<ix<< std::endl << std::flush;

        // find the new voxel index
        long iv = IndicesToIndex(ix, iy, iz);

        if (!strcmp(dataset, "segmentations")) {
            segment[iv] = 1;
        }
        else if (!strcmp(dataset, "skeletons")) {
            segment[iv] = 1;
        }
        else if (!strcmp(dataset, "synapses")) {
            segment[iv] = 3;
            synapses.insert(iv);
        }
        else if (!strcmp(dataset, "somae")) {
            segment[iv] = 4;
        }
        else if (!strcmp(dataset, "volumetric_somae/surfaces")) {
            segment[iv] = 4;
        }
        else { fprintf(stderr, "Unrecognized point cloud: %s.\n", dataset); exit(-1); }
    }

    // close file
    fclose(fp);
}
