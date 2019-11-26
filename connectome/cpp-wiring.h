#ifndef __CPP_WIRING__
#define __CPP_WIRING__
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <stdio.h>




// function calls across cpp files
void CPPcreateDataBlock(const char *prefix, const char *lookup_table_directory, long *inp_labels, float input_resolution[3], long inp_blocksize[3], long volume_size[3], long block_ind[3],
        const char* synapses_dir, const char* somae_dir, const char* output_dir, const char* output_dir_z_inp, const char* output_dir_y_inp, const char* output_dir_x_inp);
void CppSkeletonRefinement(const char *prefix, float input_resolution[3], long inp_blocksize[3], long volume_size[3], long block_ind_begin[3], long block_ind_end[3], const char* output_dir);
void ComputeAnchorPoints(const char *prefix, const char* output_dir, long input_blocksize[3], long *z_min_wall, long *z_max_wall, long *y_min_wall, long *y_max_wall, long *x_min_wall, long *x_max_wall);
// universal variables and functions
#define OR_Z 0
#define OR_Y 1
#define OR_X 2


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

inline long PadIndex(long index_unpadded, long unpadded_sheet_size, long unpadded_row_size, long padded_sheet_size, long padded_row_size)
{

    long iz = index_unpadded / unpadded_sheet_size;
    long iy = (index_unpadded - iz * unpadded_sheet_size) / unpadded_row_size;
    long ix = index_unpadded % unpadded_row_size;

    //  pad the location by one
    long iz_padded = iz + 1;
    long iy_padded = iy + 1;
    long ix_padded = ix + 1;

    // find the new voxel index
    long index_padded = iz_padded * padded_sheet_size + iy_padded * padded_row_size + ix_padded;

    return index_padded;
}



#endif
