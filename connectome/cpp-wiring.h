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

void CppSkeletonGeneration(const char *prefix, long label, const char *lookup_table_directory, long *labels);
void CppSkeletonRefinement(const char *prefix, long label, double resolution[3]);
void CppPopulatePointCloud(const char *prefix, const char *dataset, long label);
void CppPopulatePointCloudFromH5(long label, long *labels);




// universal variables and functions

#define OR_Z 0
#define OR_Y 1
#define OR_X 2

// global variables

extern long padded_blocksize[3];
extern long input_blocksize[3];
extern long volumesize[3];
extern float resolution[3];
extern long block_z;
extern long block_y;
extern long block_x;
extern long nentries;
extern long sheet_size;
extern long row_size;
extern long infinity;
extern std::unordered_map<long, char> segment;
extern std::unordered_set<long> synapses;



//////////////////////////////////////
//// COORDINATE UTILITY FUNCTIONS ////
//////////////////////////////////////

inline void IndexToIndices(long iv, long &ix, long &iy, long &iz)
{
    iz = iv / sheet_size;
    iy = (iv - iz * sheet_size) / row_size;
    ix = iv % row_size;
}



inline long IndicesToIndex(long ix, long iy, long iz)
{
    return iz * sheet_size + iy * row_size + ix;
}



#endif
