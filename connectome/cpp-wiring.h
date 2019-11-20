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
void CPPcreateDataBlock(const char *prefix, const char *lookup_table_directory, long *inp_labels, float input_resolution[3], long inp_blocksize[3], long volume_size[3], long block_ind[3], const char* synapses_dir, const char* somae_dir, const char* skeleton_dir);

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



#endif
