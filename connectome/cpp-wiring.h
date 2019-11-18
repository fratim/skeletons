#ifndef __CPP_WIRING__
#define __CPP_WIRING__
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>



// function calls across cpp files

void CppUpdateResolution(float resolution[3]);
void CppUpdateBlocksize(long inp_blocksize[3]);
void CppUpdateVolumesize(long volumesize[3]);
void CppUpdateBlockindices(long block_z, long block_y, long block_x);
void CppUpdateDirectories(char* synapses_dir, char* somae_dir, char* skeleton_dir);

void CppSkeletonGeneration(const char *prefix, const char *lookup_table_directory, long *inp_labels);
// void CppSkeletonRefinement(const char *prefix, long segment_ID, double resolution[3]);
void CppPopulatePointCloud(const char *prefix, const char *dataset, long segment_ID);
void CppPopulatePointCloudFromH5(long *inp_labels);
void CppUpdateDirectories(const char* synapses_dir, const char* somae_dir, const char* skeleton_dir);
void CPPcreateDataBlock(const char *prefix, const char *lookup_table_directory, long *inp_labels, float input_resolution[3], long inp_blocksize[3], long volume_size[3], long block_ind[3], const char* synapses_dir, const char* somae_dir, const char* skeleton_dir);



// universal variables and functions

#define OR_Z 0
#define OR_Y 1
#define OR_X 2

// global variables

// extern long padded_blocksize[3];
// extern long input_blocksize[3];
// extern long volumesize[3];
// extern float resolution[3];
// extern long block_z;
// extern long block_y;
// extern long block_x;
// extern long nentries;
// extern long sheet_size;
// extern long row_size;
// extern long infinity;
// extern std::unordered_map<long, short> segment;
// // extern std::unordered_map<long, std::unordered_set<long> > synapses_cloud;
// extern std::unordered_set<long> IDs_in_block;
// extern std::unordered_map<long, std::unordered_map<long,short>> Pointclouds;
// extern std::unordered_map<long, std::unordered_map<long,std::unordered_set<long>>> borderpoints;
// extern std::unordered_map<long,std::unordered_set<long>> borderpoints_segment;
// extern const char *synapses_directory;
// extern const char *somae_directory;
// extern const char *skeleton_directory;



//////////////////////////////////////
//// COORDINATE UTILITY FUNCTIONS ////
//////////////////////////////////////

inline void IndexToIndices(long iv, long &ix, long &iy, long &iz, long sheet_size, long row_size)
{
    iz = iv / sheet_size;
    iy = (iv - iz * sheet_size) / row_size;
    ix = iv % row_size;
}



inline long IndicesToIndex(long ix, long iy, long iz, long sheet_size, long row_size)
{
    return iz * sheet_size + iy * row_size + ix;
}



#endif
