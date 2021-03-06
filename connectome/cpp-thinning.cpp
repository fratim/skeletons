/* c++ file to extract wiring diagram */
#include <limits>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <vector>
#include <time.h>
#include <stdio.h>
#include <string>
#include <algorithm>
#include "cpp-MinBinaryHeap.h"
#include "cpp-wiring.h"


// constant variables
static const int lookup_table_size = 1 << 23;

// DO NOT CHANGE THIS ORDERING
static const int NTHINNING_DIRECTIONS = 6;
static const int UP = 0;
static const int DOWN = 1;
static const int NORTH = 2;
static const int SOUTH = 3;
static const int EAST = 4;
static const int WEST = 5;

// lookup tables
static unsigned char *lut_simple;

// mask variables for bitwise operations
static long long_mask[26];
static unsigned char char_mask[8];
static long n26_offsets[26];
static long n6_offsets[6];

// TODO remove this at some point
double time_projSynapses = 0;
clock_t time_beforeprojSynapses = clock();

static void set_long_mask(void)
{
  long_mask[ 0] = 0x00000001;
  long_mask[ 1] = 0x00000002;
  long_mask[ 2] = 0x00000004;
  long_mask[ 3] = 0x00000008;
  long_mask[ 4] = 0x00000010;
  long_mask[ 5] = 0x00000020;
  long_mask[ 6] = 0x00000040;
  long_mask[ 7] = 0x00000080;
  long_mask[ 8] = 0x00000100;
  long_mask[ 9] = 0x00000200;
  long_mask[10] = 0x00000400;
  long_mask[11] = 0x00000800;
  long_mask[12] = 0x00001000;
  long_mask[13] = 0x00002000;
  long_mask[14] = 0x00004000;
  long_mask[15] = 0x00008000;
  long_mask[16] = 0x00010000;
  long_mask[17] = 0x00020000;
  long_mask[18] = 0x00040000;
  long_mask[19] = 0x00080000;
  long_mask[20] = 0x00100000;
  long_mask[21] = 0x00200000;
  long_mask[22] = 0x00400000;
  long_mask[23] = 0x00800000;
  long_mask[24] = 0x01000000;
  long_mask[25] = 0x02000000;
}

static void set_char_mask(void)
{
  char_mask[0] = 0x01;
  char_mask[1] = 0x02;
  char_mask[2] = 0x04;
  char_mask[3] = 0x08;
  char_mask[4] = 0x10;
  char_mask[5] = 0x20;
  char_mask[6] = 0x40;
  char_mask[7] = 0x80;
}

static void PopulateOffsets(long padded_blocksize[3])
{
  n26_offsets[0] = -1 * padded_blocksize[OR_Y] * padded_blocksize[OR_X] - padded_blocksize[OR_X] - 1;
  n26_offsets[1] = -1 * padded_blocksize[OR_Y] * padded_blocksize[OR_X] - padded_blocksize[OR_X];
  n26_offsets[2] = -1 * padded_blocksize[OR_Y] * padded_blocksize[OR_X] - padded_blocksize[OR_X] + 1;
  n26_offsets[3] = -1 * padded_blocksize[OR_Y] * padded_blocksize[OR_X] - 1;
  n26_offsets[4] = -1 * padded_blocksize[OR_Y] * padded_blocksize[OR_X];
  n26_offsets[5] = -1 * padded_blocksize[OR_Y] * padded_blocksize[OR_X] + 1;
  n26_offsets[6] = -1 * padded_blocksize[OR_Y] * padded_blocksize[OR_X] + padded_blocksize[OR_X] - 1;
  n26_offsets[7] = -1 * padded_blocksize[OR_Y] * padded_blocksize[OR_X] + padded_blocksize[OR_X];
  n26_offsets[8] = -1 * padded_blocksize[OR_Y] * padded_blocksize[OR_X] + padded_blocksize[OR_X] + 1;

  n26_offsets[9] = -1 * padded_blocksize[OR_X] - 1;
  n26_offsets[10] = -1 * padded_blocksize[OR_X];
  n26_offsets[11] = -1 * padded_blocksize[OR_X] + 1;
  n26_offsets[12] = -1;
  n26_offsets[13] = +1;
  n26_offsets[14] = padded_blocksize[OR_X] - 1;
  n26_offsets[15] = padded_blocksize[OR_X];
  n26_offsets[16] = padded_blocksize[OR_X] + 1;

  n26_offsets[17] = padded_blocksize[OR_Y] * padded_blocksize[OR_X] - padded_blocksize[OR_X] - 1;
  n26_offsets[18] = padded_blocksize[OR_Y] * padded_blocksize[OR_X] - padded_blocksize[OR_X];
  n26_offsets[19] = padded_blocksize[OR_Y] * padded_blocksize[OR_X] - padded_blocksize[OR_X] + 1;
  n26_offsets[20] = padded_blocksize[OR_Y] * padded_blocksize[OR_X] - 1;
  n26_offsets[21] = padded_blocksize[OR_Y] * padded_blocksize[OR_X];
  n26_offsets[22] = padded_blocksize[OR_Y] * padded_blocksize[OR_X] + 1;
  n26_offsets[23] = padded_blocksize[OR_Y] * padded_blocksize[OR_X] + padded_blocksize[OR_X] - 1;
  n26_offsets[24] = padded_blocksize[OR_Y] * padded_blocksize[OR_X] + padded_blocksize[OR_X];
  n26_offsets[25] = padded_blocksize[OR_Y] * padded_blocksize[OR_X] + padded_blocksize[OR_X] + 1;

  // use this order to go UP, DOWN, NORTH, SOUTH, EAST, WEST
  // DO NOT CHANGE THIS ORDERING
  n6_offsets[0] = -1 * padded_blocksize[OR_X]; //neg y direction
  n6_offsets[1] = padded_blocksize[OR_X]; // pos y direction
  n6_offsets[2] = -1 * padded_blocksize[OR_Y] * padded_blocksize[OR_X]; // neg z direction
  n6_offsets[3] = padded_blocksize[OR_Y] * padded_blocksize[OR_X]; // pos z direction
  n6_offsets[4] = +1; // pos x direction
  n6_offsets[5] = -1; // neg x direction
}

// very simple double linked list data structure
typedef struct {
  long iv, ix, iy, iz;
  void *next;
  void *prev;
} ListElement;

typedef struct {
  void *first;
  void *last;
} List;

typedef struct {
  long iv, ix, iy, iz;
} Voxel;

typedef struct {
  Voxel v;
  ListElement *ptr;
  void *next;
} Cell;

typedef struct {
  Cell *head;
  Cell *tail;
  int length;
} PointList;

static void NewSurfaceVoxel(long iv, long ix, long iy, long iz, List &surface_voxels);
static void RemoveSurfaceVoxel(ListElement *LE, List &surface_voxels);
static void CreatePointList(PointList *s);
static void AddToList(PointList *s, Voxel e, ListElement *ptr);
static Voxel GetFromList(PointList *s, ListElement **ptr);
static void DestroyPointList(PointList *s);
static void InitializeLookupTables(const char *lookup_table_directory);
static bool Simple26_6(unsigned int neighbors);

class DataBlock{
protected:
  long input_blocksize[3] = {-1,-1,-1};
  long volumesize[3] = {-1,-1,-1};
  float resolution[3] = {-1,-1,-1};
  long block_ind[3] = {-1,-1,-1};
  long block_ind_min[3] = {-1,-1,-1};
  long block_ind_max[3] = {-1,-1,-1};
  long padded_row_size_block = -1;
  long padded_sheet_size_block = -1;
  long input_row_size_block = -1;
  long input_sheet_size_block = -1;
  const char *prefix;
  const char *synapses_directory;
  const char *output_directory;

public:
  long padded_blocksize[3];
  uoSet IDs_to_process = uoSet();
  uoSet IDs_in_block = uoSet();
  pointclouds_obj Pointclouds = pointclouds_obj();
  borderpoints_obj borderpoints = borderpoints_obj();
  map_numTovector anchors_comp = map_numTovector();
  map_numTovector synapses = map_numTovector();
  map_numTovector synapses_off = map_numTovector();
  map_numToset somae_interiorpoints = map_numToset();
  map_numToset somae_surfacepoints = map_numToset();

  DataBlock(const char* prefix_inp, float input_resolution[3], long inp_blocksize[3], long volume_size[3], long block_ind_inp[3], long block_ind_start_inp[3], long block_ind_end_inp[3], const char* synapses_dir, const char* output_dir){

    resolution[OR_Z] = input_resolution[OR_Z];
    resolution[OR_Y] = input_resolution[OR_Y];
    resolution[OR_X] = input_resolution[OR_X];

    // std::cout << "Resolution set to: " << resolution[OR_Z] << "," << resolution[OR_Y] << "," << resolution[OR_X] << "," << std::endl;

    input_blocksize[OR_Z] = inp_blocksize[OR_Z];
    input_blocksize[OR_Y] = inp_blocksize[OR_Y];
    input_blocksize[OR_X] = inp_blocksize[OR_X];

    padded_blocksize[OR_Z] = inp_blocksize[OR_Z]+2;
    padded_blocksize[OR_Y] = inp_blocksize[OR_Y]+2;
    padded_blocksize[OR_X] = inp_blocksize[OR_X]+2;

    padded_sheet_size_block = padded_blocksize[OR_Y] * padded_blocksize[OR_X];
    padded_row_size_block = padded_blocksize[OR_X];

    input_sheet_size_block = input_blocksize[OR_Y] * input_blocksize[OR_X];
    input_row_size_block = input_blocksize[OR_X];

    // std::cout << "Blocksize input set to: " << input_blocksize[OR_Z] << "," << input_blocksize[OR_Y] << "," << input_blocksize[OR_X] << "," << std::endl;
    // std::cout << "Blocksize padded set to: " << padded_blocksize[OR_Z] << "," << padded_blocksize[OR_Y] << "," << padded_blocksize[OR_X] << "," << std::endl;
    // std::cout << "Sheetsize padded set to: " << padded_sheet_size_block << std::endl;
    // std::cout << "Rowsize padded set to: " << padded_row_size_block << std::endl;
    // std::cout << "Sheetsize input set to: " << input_sheet_size_block << std::endl;
    // std::cout << "Rowsize input set to: " << input_row_size_block << std::endl;

    volumesize[OR_Z] = volume_size[OR_Z];
    volumesize[OR_Y] = volume_size[OR_Y];
    volumesize[OR_X] = volume_size[OR_X];

    //std::cout << "Volumesize set to: " << volumesize[OR_Z] << "," << volumesize[OR_Y] << "," << volumesize[OR_X] << "," << std::endl;

    block_ind[OR_Z]= block_ind_inp[OR_Z];
    block_ind[OR_Y] = block_ind_inp[OR_Y];
    block_ind[OR_X] = block_ind_inp[OR_X];

    block_ind_min[OR_Z]= block_ind_start_inp[OR_Z];
    block_ind_min[OR_Y] = block_ind_start_inp[OR_Y];
    block_ind_min[OR_X] = block_ind_start_inp[OR_X];

    block_ind_max[OR_Z]= block_ind_end_inp[OR_Z];
    block_ind_max[OR_Y] = block_ind_end_inp[OR_Y];
    block_ind_max[OR_X] = block_ind_end_inp[OR_X];

    //std::cout << "Block indices set to: " << block_ind[OR_Z] << "," << block_ind[OR_Y] << "," << block_ind[OR_X] << "," << std::endl;

    prefix = prefix_inp;

    synapses_directory = synapses_dir;
    output_directory = output_dir;

    // std::cout << "Directories set. " << std::endl;

  }

  DataBlock(DataBlock &Block)
  {
    std::copy(std::begin(Block.input_blocksize), std::end(Block.input_blocksize), std::begin(input_blocksize));
    std::copy(std::begin(Block.padded_blocksize), std::end(Block.padded_blocksize), std::begin(padded_blocksize));
    std::copy(std::begin(Block.volumesize), std::end(Block.volumesize), std::begin(volumesize));
    std::copy(std::begin(Block.resolution), std::end(Block.resolution), std::begin(resolution));
    std::copy(std::begin(Block.block_ind), std::end(Block.block_ind), std::begin(block_ind));
    std::copy(std::begin(Block.block_ind_min), std::end(Block.block_ind_min), std::begin(block_ind_min));
    std::copy(std::begin(Block.block_ind_max), std::end(Block.block_ind_max), std::begin(block_ind_max));
    output_directory = Block.output_directory;
    prefix = Block.prefix;

    padded_row_size_block = Block.padded_row_size_block;
    padded_sheet_size_block = Block.padded_sheet_size_block;

    input_row_size_block = Block.input_row_size_block;
    input_sheet_size_block = Block.input_sheet_size_block;
  }

  void CppPopulatePointCloudFromH5(long *inp_labels)
  {

    long n_points = input_blocksize[OR_Z]*input_blocksize[OR_Y]*input_blocksize[OR_X];

    //std::cout << "n_points is: " << n_points << std::endl;

    long z_max =input_blocksize[OR_Z]-1;
    long y_max =input_blocksize[OR_Y]-1;
    long x_max =input_blocksize[OR_X]-1;

    for (long up_iv_local = 0; up_iv_local < n_points; up_iv_local++){

      // get segment_ID of current index and skip if is zero
      long curr_label = inp_labels[up_iv_local];
      if (!curr_label) continue;

      // find the new voxel index
      long p_iv_local = PadIndex(up_iv_local, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);

      // check if pointcloud of this segment_ID already exists, otherwise add new pointcloud
      if (Pointclouds.find(curr_label) == Pointclouds.end()) {
        Pointclouds[curr_label] = map_numTochar();
        // std::cout << "New segment_ID detected: " << curr_label << std::endl << std::flush;
        IDs_in_block.insert(curr_label);
      }


      if (somae_interiorpoints[curr_label].find(p_iv_local)==somae_interiorpoints[curr_label].end()){
        Pointclouds[curr_label][p_iv_local] = 1;
      }


      long ix, iy, iz;
      IndexToIndices(up_iv_local, ix, iy, iz, input_sheet_size_block, input_row_size_block);



      // add points to borderpoints component, if on any of the 6 walls
      // max key = 1
      // min key = 0
      // do not add points that are on edges or corners, as those otherwise could never be removed

      if (iz == z_max) {
          if (not (iy == 0 or iy == y_max or ix == 0 or ix == x_max)) borderpoints[curr_label][OR_Z][1].insert(p_iv_local);
      }
      if (iy == y_max) {
          if (not (iz == 0 or iz == z_max or ix == 0 or ix == x_max)) borderpoints[curr_label][OR_Y][1].insert(p_iv_local);
      }
      if (ix == x_max) {
          if (not (iz == 0 or iz == z_max or iy == 0 or iy == y_max)) borderpoints[curr_label][OR_X][1].insert(p_iv_local);
      }
      if (iz == 0) {
          if (not (iy == 0 or iy == y_max or ix == 0 or ix == x_max)) borderpoints[curr_label][OR_Z][0].insert(p_iv_local);
      }
      if (iy == 0) {
          if (not (iz == 0 or iz == z_max or ix == 0 or ix == x_max)) borderpoints[curr_label][OR_Y][0].insert(p_iv_local);
      }
      if (ix == 0) {
          if (not (iz == 0 or iz == z_max or iy == 0 or iy == y_max)) borderpoints[curr_label][OR_X][0].insert(p_iv_local);
      }
    }

    for (map_numToset::iterator itr = somae_surfacepoints.begin(); itr!=somae_surfacepoints.end(); ++itr){
      long label = itr->first;
      //std::cout << "Seg ID: " << label << std::endl;
      //std::cout << "points: " << somae_surfacepoints[label].size() << std::endl;
      long p_iv_local;
      for (uoSet::iterator itr2 = somae_surfacepoints[label].begin(); itr2!=somae_surfacepoints[label].end(); ++itr2){
        p_iv_local = PadIndex(*itr2, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);
        Pointclouds[label][p_iv_local] = 4;
      }
    }


  }

  void CppPopulateSomaeFromH5(long *inp_somae)
  {
    // save characteristics of the downsampled volume, which are only needed and therfore defined within the function

    long dsp = SOMA_DSP;
    long input_blocksize_dsp[3] = {input_blocksize[OR_Z]/dsp, input_blocksize[OR_Y]/dsp, input_blocksize[OR_X]/dsp};
    long padded_blocksize_dsp[3] = {input_blocksize[OR_Z]/dsp+2, input_blocksize[OR_Y]/dsp+2, input_blocksize[OR_X]/dsp+2};
    long input_sheet_size_dsp = input_blocksize_dsp[OR_Y]*input_blocksize_dsp[OR_X];
    long input_row_size_dsp = input_blocksize_dsp[OR_X];
    long padded_sheet_size_dsp = padded_blocksize_dsp[OR_Y]*padded_blocksize_dsp[OR_X];
    long padded_row_size_dsp = padded_blocksize_dsp[OR_X];

    long n_points = input_blocksize_dsp[OR_Z]*input_blocksize_dsp[OR_Y]*input_blocksize_dsp[OR_X];

    //std::cout << "n_points somae is: " << n_points << std::endl;

    for (long up_iv_local_dsp = 0; up_iv_local_dsp < n_points; up_iv_local_dsp++){

      // get segment_ID of current index and skip if is zero
      long curr_label = inp_somae[up_iv_local_dsp];
      if (!curr_label) continue;

      // get downsampled coordinates
      long ix_dsp, iy_dsp, iz_dsp;
      IndexToIndices(up_iv_local_dsp, ix_dsp, iy_dsp, iz_dsp, input_sheet_size_dsp, input_row_size_dsp);

      long n6_offsets_dsp[6];
      n6_offsets_dsp[0] = -1 * padded_blocksize_dsp[OR_X]; //neg y direction
      n6_offsets_dsp[1] = padded_blocksize_dsp[OR_X]; // neg y direction
      n6_offsets_dsp[2] = -1 * padded_blocksize_dsp[OR_Y] * padded_blocksize_dsp[OR_X]; // neg z direction
      n6_offsets_dsp[3] = padded_blocksize_dsp[OR_Y] * padded_blocksize_dsp[OR_X]; // pos z direction
      n6_offsets_dsp[4] = +1; // pos x direction
      n6_offsets_dsp[5] = -1; // neg x direction

      // check if this is a surface voxel in either of the 6 direction
      bool isSurface_z_pos = 0;
      bool isSurface_y_pos = 0;
      bool isSurface_x_pos = 0;
      bool isSurface_z_neg = 0;
      bool isSurface_y_neg = 0;
      bool isSurface_x_neg = 0;

      long p_iv_local_dsp = PadIndex(up_iv_local_dsp, input_sheet_size_dsp, input_row_size_dsp, padded_sheet_size_dsp, padded_row_size_dsp);

      // check all 6-connected neighboring voxels, determining if this is a soma surface point
      for (long dir = 0; dir < NTHINNING_DIRECTIONS; dir++) {

        long p_neighbor_index = p_iv_local_dsp + n6_offsets_dsp[dir];

        // skip the boundary elements
        long ii, ij, ik;
        IndexToIndices(p_neighbor_index, ii, ij, ik, padded_sheet_size_dsp, padded_row_size_dsp);
        if ((ii == 0) or (ii == padded_blocksize_dsp[OR_X] - 1)) continue;
        if ((ij == 0) or (ij == padded_blocksize_dsp[OR_Y] - 1)) continue;
        if ((ik == 0) or (ik == padded_blocksize_dsp[OR_Z] - 1)) continue;

        long up_neighbor_index = UnpadIndex(p_neighbor_index, input_sheet_size_dsp, input_row_size_dsp, padded_sheet_size_dsp, padded_row_size_dsp);

        if (inp_somae[up_neighbor_index] != curr_label) {
          if (dir==0) isSurface_y_neg = 1;
          if (dir==1) isSurface_y_pos = 1;
          if (dir==2) isSurface_z_neg = 1;
          if (dir==3) isSurface_z_pos = 1;
          if (dir==4) isSurface_x_pos = 1;
          if (dir==5) isSurface_x_neg = 1;
        }
      }

      // get normal coordinates
      long ix = ix_dsp*dsp;
      long iy = iy_dsp*dsp;
      long iz = iz_dsp*dsp;

      long up_iv_local_add;
      long p_iv_local_add;

      // general idea: each somae voxel represents a volume of (dsp*dsp*dsp) in the segmentation space.
      // if the voxel is a surface voxel, the respective face of this volume is added to the segmentation as a soma surface.
      // otherwise, respective face of the volume is added to the somae interior n_points
      // the inside voxels of the volume (a total of dsp*dsp*dsp-dsp*dsp*6) are always added to the somae interior points
      // when loading the segmentation, somae interior points can the be omitted

      if (isSurface_z_neg){
        for (int iv = iy; iv<iy+dsp; iv++){
          for (int iu = ix; iu<ix+dsp; iu++){
            int iw = iz;
            up_iv_local_add = IndicesToIndex(iu, iv, iw, input_sheet_size_block, input_row_size_block);
            p_iv_local_add = PadIndex(up_iv_local_add, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);
            // Pointclouds[curr_label][p_iv_local_add] = 4;
            somae_surfacepoints[curr_label].insert(up_iv_local_add);
          }
        }
      }

      else{
        for (int iv = iy; iv<iy+dsp; iv++){
          for (int iu = ix; iu<ix+dsp; iu++){
            int iw = iz;
            up_iv_local_add = IndicesToIndex(iu, iv, iw, input_sheet_size_block, input_row_size_block);
            p_iv_local_add = PadIndex(up_iv_local_add, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);
            somae_interiorpoints[curr_label].insert(p_iv_local_add);
          }
        }
      }

      if (isSurface_y_neg){
        for (int iw = iz; iw<iz+dsp; iw++){
          for (int iu = ix; iu<ix+dsp; iu++){
            int iv = iy;
            up_iv_local_add = IndicesToIndex(iu, iv, iw, input_sheet_size_block, input_row_size_block);
            p_iv_local_add = PadIndex(up_iv_local_add, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);
            // Pointclouds[curr_label][p_iv_local_add] = 4;
            somae_surfacepoints[curr_label].insert(up_iv_local_add);
          }
        }
      }

      else{
        for (int iw = iz; iw<iz+dsp; iw++){
          for (int iu = ix; iu<ix+dsp; iu++){
            int iv = iy;
            up_iv_local_add = IndicesToIndex(iu, iv, iw, input_sheet_size_block, input_row_size_block);
            p_iv_local_add = PadIndex(up_iv_local_add, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);
            somae_interiorpoints[curr_label].insert(p_iv_local_add);
          }
        }
      }

      if (isSurface_x_neg){
        for (int iw = iz; iw<iz+dsp; iw++){
          for (int iv = iy; iv<iy+dsp; iv++){
            int iu = ix;
            up_iv_local_add = IndicesToIndex(iu, iv, iw, input_sheet_size_block, input_row_size_block);
            p_iv_local_add = PadIndex(up_iv_local_add, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);
            // Pointclouds[curr_label][p_iv_local_add] = 4;
            somae_surfacepoints[curr_label].insert(up_iv_local_add);
          }
        }
      }

      else{
        for (int iw = iz; iw<iz+dsp; iw++){
          for (int iv = iy; iv<iy+dsp; iv++){
            int iu = ix;
            up_iv_local_add = IndicesToIndex(iu, iv, iw, input_sheet_size_block, input_row_size_block);
            p_iv_local_add = PadIndex(up_iv_local_add, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);
            somae_interiorpoints[curr_label].insert(p_iv_local_add);
          }
        }
      }

      if (isSurface_z_pos){
        for (int iv = iy; iv<iy+dsp; iv++){
          for (int iu = ix; iu<ix+dsp; iu++){
            int iw = iz+dsp-1;
            up_iv_local_add = IndicesToIndex(iu, iv, iw, input_sheet_size_block, input_row_size_block);
            p_iv_local_add = PadIndex(up_iv_local_add, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);
            // Pointclouds[curr_label][p_iv_local_add] = 4;
            somae_surfacepoints[curr_label].insert(up_iv_local_add);
          }
        }
      }

      else{
        for (int iv = iy; iv<iy+dsp; iv++){
          for (int iu = ix; iu<ix+dsp; iu++){
            int iw = iz+dsp-1;
            up_iv_local_add = IndicesToIndex(iu, iv, iw, input_sheet_size_block, input_row_size_block);
            p_iv_local_add = PadIndex(up_iv_local_add, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);
            somae_interiorpoints[curr_label].insert(p_iv_local_add);
          }
        }
      }

      if (isSurface_y_pos){
        for (int iw = iz; iw<iz+dsp; iw++){
          for (int iu = ix; iu<ix+dsp; iu++){
            int iv = iy+dsp-1;
            up_iv_local_add = IndicesToIndex(iu, iv, iw, input_sheet_size_block, input_row_size_block);
            p_iv_local_add = PadIndex(up_iv_local_add, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);
            // Pointclouds[curr_label][p_iv_local_add] = 4;
            somae_surfacepoints[curr_label].insert(up_iv_local_add);
          }
        }
      }

      else{
        for (int iw = iz; iw<iz+dsp; iw++){
          for (int iu = ix; iu<ix+dsp; iu++){
            int iv = iy+dsp-1;
            up_iv_local_add = IndicesToIndex(iu, iv, iw, input_sheet_size_block, input_row_size_block);
            p_iv_local_add = PadIndex(up_iv_local_add, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);
            somae_interiorpoints[curr_label].insert(p_iv_local_add);
          }
        }
      }

      if (isSurface_x_pos){
        for (int iw = iz; iw<iz+dsp; iw++){
          for (int iv = iy; iv<iy+dsp; iv++){
            int iu = ix+dsp-1;
            up_iv_local_add = IndicesToIndex(iu, iv, iw, input_sheet_size_block, input_row_size_block);
            p_iv_local_add = PadIndex(up_iv_local_add, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);
            // Pointclouds[curr_label][p_iv_local_add] = 4;
            somae_surfacepoints[curr_label].insert(up_iv_local_add);
          }
        }
      }

      else{
        for (int iw = iz; iw<iz+dsp; iw++){
          for (int iv = iy; iv<iy+dsp; iv++){
            int iu = ix+dsp-1;
            up_iv_local_add = IndicesToIndex(iu, iv, iw, input_sheet_size_block, input_row_size_block);
            p_iv_local_add = PadIndex(up_iv_local_add, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);
            somae_interiorpoints[curr_label].insert(p_iv_local_add);
          }
        }
      }

      for (int iw = iz+1; iw<iz+dsp-1; iw++){
        for (int iv = iy+1; iv<iy+dsp-1; iv++){
          for (int iu = ix+1; iu<ix+dsp-1; iu++){
            up_iv_local_add = IndicesToIndex(iu, iv, iw, input_sheet_size_block, input_row_size_block);
            p_iv_local_add = PadIndex(up_iv_local_add, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);
            somae_interiorpoints[curr_label].insert(p_iv_local_add);
          }
        }
      }

    }

    // somae surfaces are written to one file per block per neuron
    WriteSomaeSurface();

  }

  int ReadSynapses(void)
  {

    // read the synapses
    char synapse_filename[4096];
    snprintf(synapse_filename, 4096, "%s/%s/%s-synapses-%04ldz-%04ldy-%04ldx.pts", synapses_directory, prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

    FILE *fp = fopen(synapse_filename, "rb");
    if (!fp) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

    long nneurons;
    ReadHeader(fp, nneurons);

    for (long n = 0; n < nneurons; ++n) {
      // get the label and number of synapses
      long segment_ID;
      long nsynapses;

      if (fread(&segment_ID, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
      if (fread(&nsynapses, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

      // ignore the global coordinates
      for (long is = 0; is < nsynapses; ++is) {
        long dummy_index;
        if (fread(&dummy_index, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
      }

      // add the local coordinates (offset by one in each direction)
      for (long is = 0; is < nsynapses; ++is) {

        long up_iv_local;
        if (fread(&up_iv_local, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
        // add them to the synapses vector (unpadded)

        // find the new voxel index
        long p_iv_local = PadIndex(up_iv_local, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);

        // set ID of this index to 3 (==Synapse)
        // only add synapse if it is on the body
        // if synapse is not on the body, it is added to synapses_off, which are projected onto the surface in a later step
        if (Pointclouds[segment_ID].find(p_iv_local) == Pointclouds[segment_ID].end()) {
            synapses_off[segment_ID].push_back(up_iv_local);
        }
        else {
            Pointclouds[segment_ID][p_iv_local] = 3;
            synapses[segment_ID].push_back(up_iv_local);
        }

      }

    }


    // close file
    fclose(fp);

    return 1;

  }

  int ReadAnchorpoints(void)
  {
    // Read in precomputed anchorpoints for the surfaces needed

    // read in z anchors (Zmin)
    if (block_ind[OR_Z]>block_ind_min[OR_Z]){
      char anchor_ID[] = "Anchors_Comp_ZMin";
      ReadAnchorFile(anchor_ID);
    }

    // read in z anchors (Zmax)
    if (block_ind[OR_Z]<block_ind_max[OR_Z]){
      char anchor_ID[] = "Anchors_Comp_ZMax";
      ReadAnchorFile(anchor_ID);
    }

    // read in z anchors (Ymin)
    if (block_ind[OR_Y]>block_ind_min[OR_Y]){
      char anchor_ID[] = "Anchors_Comp_YMin";
      ReadAnchorFile(anchor_ID);
    }

    // read in z anchors (Ymax)
    if (block_ind[OR_Y]<block_ind_max[OR_Y]){
      char anchor_ID[] = "Anchors_Comp_YMax";
      ReadAnchorFile(anchor_ID);
    }

    // read in z anchors (Xmin)
    if (block_ind[OR_X]>block_ind_min[OR_X]){
      char anchor_ID[] = "Anchors_Comp_XMin";
      ReadAnchorFile(anchor_ID);
    }

    // read in z anchors (Xmax)
    if (block_ind[OR_X]<block_ind_max[OR_X]){
      char anchor_ID[] = "Anchors_Comp_XMax";
      ReadAnchorFile(anchor_ID);
    }

    return 1;

  }

  void ReadAnchorFile(const char* identifier)
  {
    // make filename of adjacent z lock (in negative direction)
    char output_filename[4096];
    sprintf(output_filename, "%s/output-%04ldz-%04ldy-%04ldx/anchorpoints_computed/%s/%s-%s-%04ldz-%04ldy-%04ldx.pts",
    output_directory,  block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X], prefix, prefix, identifier,
    block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

    FILE *fp = fopen(output_filename, "rb");
    if (!fp) { fprintf(stderr, "Failed to read %s.\n", output_filename); exit(-1); }

    long nsegments;
    ReadHeader(fp, nsegments);
    long checksum = 0;

    for (long i=0; i<nsegments; i++) {

      long seg_ID;
      long n_anchors;

      if (fread(&seg_ID, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s\n", output_filename); exit(-1); }
      if (fread(&n_anchors, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s\n", output_filename); exit(-1); }

      // ignore global indices
      for (long pos=0; pos<n_anchors; pos++) {

        long dummy;
        if (fread(&dummy, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s\n", output_filename); exit(-1); }
        checksum += dummy;
      }

      for (long pos=0; pos<n_anchors; pos++) {

        long up_iv_local;
        if (fread(&up_iv_local, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s\n", output_filename); exit(-1); }
        checksum += up_iv_local;

        long ix,iy,iz;
        IndexToIndices(up_iv_local, ix, iy, iz, input_sheet_size_block, input_row_size_block);

        long p_iv_local = PadIndex(up_iv_local, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);

        Pointclouds[seg_ID][p_iv_local] = 3;
        anchors_comp[seg_ID].push_back(up_iv_local);
      }
    }

    long checksum_read;
    if (fread(&checksum_read, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s\n", output_filename); exit(-1);}
    if (checksum!=checksum_read) { fprintf(stderr, "Read ERROR. Checksum incorrect for %s \n", output_filename); exit(-1);}

    // close file
    fclose(fp);


  }

  void writeIDs(void)
  {
    // write IDs present in a block and Ids being processed to output files
    {
      // write IDs that are processed
      char output_filename_IDs[4096];
      sprintf(output_filename_IDs, "%s/output-%04ldz-%04ldy-%04ldx/segment_IDs/%s/%s-IDs_processed-%04ldz-%04ldy-%04ldx.pts",
      output_directory,  block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X],
      prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

      FILE *fpid = fopen(output_filename_IDs, "wb");
      if (!fpid) { fprintf(stderr, "Failed to open %s\n", output_filename_IDs); exit(-1); }

      long num = IDs_to_process.size();

      WriteHeader(fpid, num);

      for (uoSet::iterator iter = IDs_to_process.begin(); iter != IDs_to_process.end(); ++iter) {

        long seg_ID = *iter;

        // IDs
        if (fwrite(&seg_ID, sizeof(long), 1, fpid) != 1) { fprintf(stderr, "Failed to write to IDs in Block to %s \n", output_filename_IDs); exit(-1); }
      }
      fclose(fpid);
    }
    {
      // write all IDs present in this Block
      char output_filename_IDs[4096];
      sprintf(output_filename_IDs, "%s/output-%04ldz-%04ldy-%04ldx/segment_IDs/%s/%s-IDs_present-%04ldz-%04ldy-%04ldx.pts",
      output_directory, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X],
      prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

      FILE *fpid = fopen(output_filename_IDs, "wb");
      if (!fpid) { fprintf(stderr, "Failed to open %s\n", output_filename_IDs); exit(-1); }

      long num = IDs_in_block.size();

      WriteHeader(fpid, num);

      for (uoSet::iterator iter = IDs_in_block.begin(); iter != IDs_in_block.end(); ++iter) {

        long seg_ID = *iter;

        // IDs
        if (fwrite(&seg_ID, sizeof(long), 1, fpid) != 1) { fprintf(stderr, "Failed to write to IDs in Block to %s \n", output_filename_IDs); exit(-1); }
      }
      fclose(fpid);
    }
  }

  void WriteHeader(FILE *fp, long &num)
  {

    // write the header parameters to the top of the output file
    int check = 0;
    int size_l = sizeof(long);

    check += fwrite(&(volumesize[OR_Z]), size_l, 1, fp);
    check += fwrite(&(volumesize[OR_Y]), size_l, 1, fp);
    check += fwrite(&(volumesize[OR_X]), size_l, 1, fp);
    check += fwrite(&(input_blocksize[OR_Z]), size_l, 1, fp);
    check += fwrite(&(input_blocksize[OR_Y]), size_l, 1, fp);
    check += fwrite(&(input_blocksize[OR_X]), size_l, 1, fp);
    check += fwrite(&num, size_l, 1, fp);

    if (check != 7) { fprintf(stderr, "Failed to write file in writeheader\n"); exit(-1); }
  }

  void ReadHeader(FILE *fp, long &num)
  {
    // read the header parameters from the top of the file and check if they agree

    long z_input_volume_size, y_input_volume_size, x_input_volume_size;
    long z_input_block_size, y_input_block_size, x_input_block_size;

    if (fread(&z_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read rty.\n"); exit(-1); }
    if (fread(&y_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read rty.\n"); exit(-1); }
    if (fread(&x_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read rty.\n"); exit(-1); }

    if (z_input_volume_size != volumesize[OR_Z]) { fprintf(stderr, "Volume Size not equal to input volume size.\n"); exit(-1); }
    if (y_input_volume_size != volumesize[OR_Y]) { fprintf(stderr, "Volume Size not equal to input volume size.\n"); exit(-1); }
    if (x_input_volume_size != volumesize[OR_X]) { fprintf(stderr, "Volume Size not equal to input volume size.\n"); exit(-1); }


    if (fread(&z_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read rty.\n"); exit(-1); }
    if (fread(&y_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read rty.\n"); exit(-1); }
    if (fread(&x_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read rty.\n"); exit(-1); }

    if (z_input_block_size != input_blocksize[OR_Z]) { fprintf(stderr, "Block Size not equal to input block size.\n"); exit(-1); }
    if (y_input_block_size != input_blocksize[OR_Y]) { fprintf(stderr, "Block Size not equal to input block size.\n"); exit(-1); }
    if (x_input_block_size != input_blocksize[OR_X]) { fprintf(stderr, "Block Size not equal to input block size.\n"); exit(-1); }

    if (fread(&num, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read rty.\n"); exit(-1); }

  }

  void WriteProjectedSynapses(void)
  {
    //get number of anchor points
    long n_neurons = IDs_to_process.size();

    // create an output file for the points
    char output_filename[4096];
    sprintf(output_filename, "%s/output-%04ldz-%04ldy-%04ldx/synapses_projected/%s/%s-synapses_projected-%04ldz-%04ldy-%04ldx.pts",
    output_directory,  block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X],
    prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

    FILE *wfp = fopen(output_filename, "wb");
    if (!wfp) { fprintf(stderr, "Failed to open %s\n", output_filename); exit(-1); }

    // write the characteristics header
    WriteHeader(wfp, n_neurons);
    long checksum = 0;

    for (map_numTovector::iterator itr = synapses.begin(); itr!=synapses.end(); ++itr){
      long seg_id = itr->first;

      if (IDs_to_process.count(seg_id)==0) continue;

      long n_synapses = synapses[seg_id].size();

      if (fwrite(&seg_id, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
      if (fwrite(&n_synapses, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
      long index_local[n_synapses];
      long pos = 0;

      for (vec::iterator itr2 = synapses[seg_id].begin(); itr2!=synapses[seg_id].end(); ++itr2, ++pos){

        long up_iv_local = *itr2;
        long up_iv_global = IndexLocalToGlobal(up_iv_local, block_ind, input_blocksize, volumesize);
        if (fwrite(&up_iv_global, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
        checksum += up_iv_global;

        index_local[pos] = up_iv_local;
        checksum += up_iv_local;
      }

      for (int j=0; j<n_synapses; j++){
        if (fwrite(&index_local[j], sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
      }

    }

    if (fwrite(&checksum, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
    fclose(wfp);
  }

  void WriteSomaeSurface(void)
    {

      for (map_numToset::iterator itr = somae_surfacepoints.begin(); itr!=somae_surfacepoints.end(); ++itr){

        long seg_id = itr->first;
        long n_points = somae_surfacepoints[seg_id].size();

        // create an output file for the points
        char output_filename[4096];
        sprintf(output_filename, "%s/output-%04ldz-%04ldy-%04ldx/somae_surface/%s/%s-somae_surfaces-%04ldz-%04ldy-%04ldx-ID-%012ld.pts",
        output_directory,  block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X],
        prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X], seg_id);

        FILE *wfp = fopen(output_filename, "wb");
        if (!wfp) { fprintf(stderr, "Failed to open %s\n", output_filename); exit(-1); }

        // write the characteristics header
        WriteHeader(wfp, seg_id);
        long checksum = 0;

        if (fwrite(&n_points, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }

        vec index_local = vec();

        long pos = 0;

        for (uoSet::iterator itr2 = somae_surfacepoints[seg_id].begin(); itr2!=somae_surfacepoints[seg_id].end(); ++itr2, ++pos){

          long up_iv_local = *itr2;
          long up_iv_global = IndexLocalToGlobal(up_iv_local, block_ind, input_blocksize, volumesize);
          if (fwrite(&up_iv_global, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
          checksum += up_iv_global;

          index_local.push_back(up_iv_local);
          checksum += up_iv_local;
        }

        for (int j=0; j<n_points; j++){
          if (fwrite(&index_local[j], sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
        }

        if (fwrite(&checksum, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
        fclose(wfp);

      }

    }

  };

    class BlockSegment : public DataBlock{
    private:
      long initial_points_input = -1;
      long segment_ID = -1;
      List surface_voxels;
      map_numTochar segment = map_numTochar();
      map_numTonumToset borderpoints_segment = map_numTonumToset();
      std::map<long, float> widths = std::map<long, float>();

    public:
      BlockSegment(long segment_ID_inp, DataBlock &Block):DataBlock(Block)
      {

        segment_ID = segment_ID_inp;
        segment = Block.Pointclouds[segment_ID];
        borderpoints_segment = Block.borderpoints[segment_ID];

        //std::cout << "-----------------------------------" << std::endl;
        //std::cout << "Processing segment_ID " << segment_ID << std::endl;

        initial_points_input = segment.size();
        //printf("segment_ID %ld initial points: %ld\n", segment_ID, initial_points_input);

        surface_voxels.first = NULL;
        surface_voxels.last = NULL;

      }

      void SequentialThinning(DataBlock &Block)
      {
        // create a vector of surface voxels
        CollectSurfaceVoxels();
        // project the synapses on the surface (surface voxels must be indentified BEFORE
        Projectsynapses(Block);

        // clear the list of surface voxels
        while (surface_voxels.first != NULL) {
            // get the surface voxels
            ListElement *LE = (ListElement *) surface_voxels.first;
            long iv = LE->iv;
            if (segment[iv] == 2) segment[iv] = 1;
            // remove this voxel
             RemoveSurfaceVoxel(LE, surface_voxels);
        }
        // create a vector of surface voxels
        CollectSurfaceVoxels();

        time_beforeprojSynapses = clock();
        time_projSynapses += (double) (clock()-time_beforeprojSynapses) / CLOCKS_PER_SEC;

        int iteration = 0;
        long changed = 0;
        do {
          changed = ThinningIterationStep(Block);
          iteration++;

          //printf("  Iteration %d: %ld points removed\n", iteration, changed);
        } while (changed);

        //printf("Needed %d iterations\n", iteration);
      }

      void CollectSurfaceVoxels(void)
      {

        long n_surface_voxels = 0;

        // go through all voxels and check their six neighbors
        for (map_numTochar::iterator it = segment.begin(); it != segment.end(); ++it) {

          // all of these elements are either 1 or 3 and in the segment
          long index = it->first;

          // initialize widths to maximum float value
          widths[index] = std::numeric_limits<float>::max();

          long ix, iy, iz;
          IndexToIndices(index, ix, iy, iz, padded_sheet_size_block, padded_row_size_block);
          // check the 6 neighbors
          for (long dir = 0; dir < NTHINNING_DIRECTIONS; ++dir) {
            long neighbor_index = index + n6_offsets[dir];

            long ii, ij, ik;

            IndexToIndices(neighbor_index, ii, ij, ik, padded_sheet_size_block, padded_row_size_block);

            // skip the fake boundary elements
            if ((ii == 0) or (ii == padded_blocksize[OR_X] - 1)) continue;
            if ((ij == 0) or (ij == padded_blocksize[OR_Y] - 1)) continue;
            if ((ik == 0) or (ik == padded_blocksize[OR_Z] - 1)) continue;

            if (segment.find(neighbor_index) == segment.end()) {
              // this location is a boundary so create a surface voxel and break
              // cannot update it->second if it is synapse so need this test!!
              if (it->second == 1) {
                it->second = 2;
                NewSurfaceVoxel(index, ix, iy, iz, surface_voxels);
                n_surface_voxels ++;
              }

              // note this location as surface, if not on a block surface
              if (  (borderpoints_segment[OR_Y][1].find(index)==borderpoints_segment[OR_Y][1].end()) &&
                    (borderpoints_segment[OR_Z][1].find(index)==borderpoints_segment[OR_Z][1].end()) &&
                    (borderpoints_segment[OR_X][1].find(index)==borderpoints_segment[OR_X][1].end()) &&
                    (borderpoints_segment[OR_Y][0].find(index)==borderpoints_segment[OR_Y][0].end()) &&
                    (borderpoints_segment[OR_Z][0].find(index)==borderpoints_segment[OR_Z][0].end()) &&
                    (borderpoints_segment[OR_X][0].find(index)==borderpoints_segment[OR_X][0].end())    ){

                      widths[index] = 0;

              }

              break;
            }
          }
        }

        //std::cout << "Initial Surface Voxels: " << n_surface_voxels << std::endl << std::flush;

      }

      long ThinningIterationStep(DataBlock &Block)
      {
        long changed = 0;

        // iterate through every direction
        for (int direction = 0; direction < NTHINNING_DIRECTIONS; ++direction) {
          PointList deletable_points;
          ListElement *ptr;

          CreatePointList(&deletable_points);

          DetectSimpleBorderPoints(&deletable_points, direction);

          while (deletable_points.length) {
            Voxel voxel = GetFromList(&deletable_points, &ptr);

            long index = voxel.iv;
            long ix = voxel.ix;
            long iy = voxel.iy;
            long iz = voxel.iz;

            // do not remove voxel in the dirction of the outer facing surface
            if ((borderpoints_segment[OR_Y][1].find(index)!=borderpoints_segment[OR_Y][1].end()) && direction==0) continue;
            if ((borderpoints_segment[OR_Y][0].find(index)!=borderpoints_segment[OR_Y][0].end()) && direction==1) continue;
            if ((borderpoints_segment[OR_Z][1].find(index)!=borderpoints_segment[OR_Z][1].end()) && direction==2) continue;
            if ((borderpoints_segment[OR_Z][0].find(index)!=borderpoints_segment[OR_Z][0].end()) && direction==3) continue;
            if ((borderpoints_segment[OR_X][0].find(index)!=borderpoints_segment[OR_X][0].end()) && direction==4) continue;
            if ((borderpoints_segment[OR_X][1].find(index)!=borderpoints_segment[OR_X][1].end()) && direction==5) continue;

            // check if simple, if so, delete it
            unsigned int neighbors = Collect26Neighbors(ix, iy, iz);
            if (Simple26_6(neighbors)) {
              // delete the simple point
              segment[index]=0;
              // add the new surface voxels
              for (long ip = 0; ip < NTHINNING_DIRECTIONS; ++ip) {
                long neighbor_index = index + n6_offsets[ip];

                // previously not on the surface but is in the object
                // widths of voxels start at maximum and first updated when put on surface
                if (segment[neighbor_index] == 1) {
                  long iu, iv, iw;
                  IndexToIndices(neighbor_index, iu, iv, iw, padded_sheet_size_block, padded_row_size_block);
                  NewSurfaceVoxel(neighbor_index, iu, iv, iw, surface_voxels);

                  // convert to a surface point
                  segment[neighbor_index] = 2;
                }
              }

              // check all 26 neighbors to see if width is better going through this voxel
              for (long ip = 0; ip < 26; ++ip) {
                long neighbor_index = index + n26_offsets[ip];
                if (!segment[neighbor_index]) continue;

                // get this index in (x, y, z)
                long iu, iv, iw;
                IndexToIndices(neighbor_index, iu, iv, iw, padded_sheet_size_block, padded_row_size_block);

                // get the distance from the voxel to be deleted
                float diffx = resolution[OR_X] * (ix - iu);
                float diffy = resolution[OR_Y] * (iy - iv);
                float diffz = resolution[OR_Z] * (iz - iw);

                float distance = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
                float current_width = widths[neighbor_index];

                if (widths[index] + distance < current_width) {
                  widths[neighbor_index] = widths[index] + distance;
                }
              }

              // remove this from the surface voxels
              RemoveSurfaceVoxel(ptr, surface_voxels);
              changed += 1;
            }

          }

          DestroyPointList(&deletable_points);

        }

        // return the number of changes
        return changed;
      }

      unsigned int Collect26Neighbors(long ix, long iy, long iz)
      {
        unsigned int neighbors = 0;

        long index = IndicesToIndex(ix, iy, iz, padded_sheet_size_block, padded_row_size_block);

        // some of these lookups will create a new entry but the region is
        // shrinking so memory overhead is minimal
        for (long iv = 0; iv < 26; ++iv) {
          if (segment[index + n26_offsets[iv]]) neighbors |= long_mask[iv];
        }

        return neighbors;
      }

      void DetectSimpleBorderPoints(PointList *deletable_points, int direction)
      {
        ListElement *LE = (ListElement *)surface_voxels.first;
        while (LE != NULL) {

          long iv = LE->iv;
          long ix = LE->ix;
          long iy = LE->iy;
          long iz = LE->iz;

          // not a synapse endpoint (need this here since endpoints are on the list of surfaces)
          // this will only be called on things on the surface already so already in map
          if (segment[iv] == 2) {

            long value = 0;
            // is the neighbor in the corresponding direction not in the segment
            // some of these keys will not exist but will default to 0 value
            // the search region retracts in from the boundary so limited memory overhead
            // the n6_offsets are in the order UP, DOWN, NORTH, SOUTH, EAST, WEST
            value = segment[iv + n6_offsets[direction]];

            // see if the required point belongs to a different segment
            if (!value) {
              unsigned int neighbors = Collect26Neighbors(ix, iy, iz);

              // deletable point
              if (Simple26_6(neighbors)) {
                Voxel voxel;
                voxel.iv = iv;
                voxel.ix = ix;
                voxel.iy = iy;
                voxel.iz = iz;
                AddToList(deletable_points, voxel, LE);
              }

            }

          }
          LE = (ListElement *) LE->next;
        }
      }

      void WriteOutputfiles(DataBlock &Block)
      {
        // count the number of remaining points
        long num = 0;
        ListElement *LE = (ListElement *) surface_voxels.first;
        while (LE != NULL) {
          num++;
          LE = (ListElement *)LE->next;
        }

        //get number of anchor points
        long n_anchors_comp =  Block.anchors_comp[segment_ID].size();
        long n_synapses =  Block.synapses[segment_ID].size();


        // create an output file for the points
        char output_filename[4096];
        sprintf(output_filename, "%s/output-%04ldz-%04ldy-%04ldx/skeleton/%s/%s-skeleton-%04ldz-%04ldy-%04ldx-ID-%012ld.pts",
        output_directory,  block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X],
        prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X] ,segment_ID);

        // write the widths to file
        char widths_filename[4096];
        sprintf(widths_filename, "%s/output-%04ldz-%04ldy-%04ldx/widths/%s/%06ld.pts",
        output_directory, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X], prefix, segment_ID);

        FILE *wfp = fopen(output_filename, "wb");
        if (!wfp) { fprintf(stderr, "Failed to open %s\n", output_filename); exit(-1); }

        FILE *width_fp = fopen(widths_filename, "wb");
        if (!width_fp) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }

        long total_points = num+n_anchors_comp+n_synapses;

        // write the characteristics header
        WriteHeaderSegID(wfp, total_points);

        // characteristics
        WriteHeaderSegID(width_fp, total_points);

        //printf("Remaining voxels: %ld\n", num);
        //printf("Anchors: %ld\n", n_anchors_comp);
        //printf("Synapses: %ld\n", n_synapses);

        long index_local[num];
        long counter_it = 0;

        long checksum_outp = 0;

        // write surface voxels (global index)
        while (surface_voxels.first != NULL) {
          // get the surface voxels
          ListElement *LE = (ListElement *) surface_voxels.first;

          long p_iv_local = LE->iv;
          float width = widths[p_iv_local];

          long up_iv_local = UnpadIndex(p_iv_local, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);
          long up_iv_global = IndexLocalToGlobal(up_iv_local, block_ind, input_blocksize, volumesize);

          if (fwrite(&up_iv_global, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
          if (fwrite(&up_iv_global, sizeof(long), 1, width_fp) != 1) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
          if (fwrite(&width, sizeof(float), 1, width_fp) != 1)  { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }

          index_local[counter_it]=up_iv_local;

          checksum_outp += up_iv_global;
          checksum_outp += up_iv_local;

          // remove this voxel
          RemoveSurfaceVoxel(LE, surface_voxels);
          counter_it++;
        }

        // write anchor points (global index)
        long index_local_anchors_comp[n_anchors_comp];
        for (long pos=0; pos<n_anchors_comp; pos++) {

          long up_iv_local = Block.anchors_comp[segment_ID][pos];
          long up_iv_global = IndexLocalToGlobal(up_iv_local, block_ind, input_blocksize, volumesize);

          if (fwrite(&up_iv_global, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }

          long p_iv_local = PadIndex(up_iv_local, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);
          float width = widths[p_iv_local];

          if (fwrite(&up_iv_global, sizeof(long), 1, width_fp) != 1) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
          if (fwrite(&width, sizeof(float), 1, width_fp) != 1)  { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }

          index_local_anchors_comp[pos] = up_iv_local;

          checksum_outp += up_iv_global;
          checksum_outp += up_iv_local;
        }

        // write the synapses (global index)
        long index_local_synapses[n_synapses];
        for (long pos=0; pos<n_synapses; pos++) {

          long up_iv_local = Block.synapses[segment_ID][pos];
          long up_iv_global = IndexLocalToGlobal(up_iv_local, block_ind, input_blocksize, volumesize);

          if (fwrite(&up_iv_global, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }

          long p_iv_local = PadIndex(up_iv_local, input_sheet_size_block, input_row_size_block, padded_sheet_size_block, padded_row_size_block);
          float width = widths[p_iv_local];

          if (fwrite(&up_iv_global, sizeof(long), 1, width_fp) != 1) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
          if (fwrite(&width, sizeof(float), 1, width_fp) != 1)  { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }

          index_local_synapses[pos] = up_iv_local;

          checksum_outp += up_iv_global;
          checksum_outp += up_iv_local;
        }

        // write surface voxels (local index)
        for (int j=0; j<num; j++){
          if (fwrite(&index_local[j], sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
        }

        // write z anchor points computed in this step (local index)
        for (int j=0; j<n_anchors_comp; j++){
          if (fwrite(&index_local_anchors_comp[j], sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
        }

        // write the synapses (local index)
        for (int j=0; j<n_synapses; j++){
          if (fwrite(&index_local_synapses[j], sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
        }

        if (fwrite(&checksum_outp, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }

        // close the I/O files
        fclose(wfp);
        fclose(width_fp);

      }

      void WriteTimeFile(clock_t start_time_seg)
      {
        char output_filename[4096];
        sprintf(output_filename, "%s/output-%04ldz-%04ldy-%04ldx/running_times/skeleton/%s/%s-total_time-%04ldz-%04ldy-%04ldx.pts",
        output_directory,  block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X], prefix, prefix,
        block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);
        FILE * fptime = fopen (output_filename,"a");

        double seg_time = (double) (clock() - start_time_seg) / CLOCKS_PER_SEC;

        fprintf(fptime,"%12ld, %12.2f \n", initial_points_input, seg_time);

        fclose(fptime);
      }

      void WriteHeaderSegID(FILE *fp, long &num)
      {
        int check = 0;
        int size_l = sizeof(long);

        check += fwrite(&(volumesize[OR_Z]), size_l, 1, fp);
        check += fwrite(&(volumesize[OR_Y]), size_l, 1, fp);
        check += fwrite(&(volumesize[OR_X]), size_l, 1, fp);
        check += fwrite(&(input_blocksize[OR_Z]), size_l, 1, fp);
        check += fwrite(&(input_blocksize[OR_Y]), size_l, 1, fp);
        check += fwrite(&(input_blocksize[OR_X]), size_l, 1, fp);
        check += fwrite(&segment_ID, size_l, 1, fp);
        check += fwrite(&num, size_l, 1, fp);

        if (check != 8) { fprintf(stderr, "Failed to write file in writeheader\n"); exit(-1); }
      }

      void Projectsynapses(DataBlock &Block)
      {

        //std::cout << "Projecting n synapses: " << Block.synapses_off[segment_ID].size() << std::endl;

        for (vec::iterator it = Block.synapses_off[segment_ID].begin(); it != Block.synapses_off[segment_ID].end(); ++it){
          long linear_index = *it;
          long iz_unpadded, iy_unpadded, ix_unpadded;
          IndexToIndices(linear_index, ix_unpadded, iy_unpadded, iz_unpadded, input_sheet_size_block, input_row_size_block);
          long iz_padded = iz_unpadded + 1;
          long iy_padded = iy_unpadded + 1;
          long ix_padded = ix_unpadded + 1;
          long cubesize = 3;
          bool projection_found = 0;
          bool unable_to_project = 0;
          long closest_index_padded = -1;

          float* i1 = std::min_element(resolution, resolution+2);
          float resolution_min = *i1;

          while (projection_found==0&&unable_to_project==0){

            double min_distance = (double)(cubesize*resolution_min);

            for (long iw = iz_padded - cubesize; iw <= iz_padded + cubesize; ++iw) {
              for (long iv = iy_padded - cubesize; iv <= iy_padded + cubesize; ++iv) {
                for (long iu = ix_padded - cubesize; iu <= ix_padded + cubesize; ++iu) {


                  long index_padded = IndicesToIndex(iu, iv, iw, padded_sheet_size_block, padded_row_size_block);
                  if (segment.find(index_padded)!=segment.end()){
                    if (segment[index_padded]==2){
                      // get the distance
                      double deltaz = resolution[OR_Z] * (double)(iw - iz_padded);
                      double deltay = resolution[OR_Y] * (double)(iv - iy_padded);
                      double deltax = resolution[OR_X] * (double)(iu - ix_padded);
                      double distance = sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);

                      // check if distance is smaller than smallest distance and smaller than radius of cube
                      if (distance<min_distance){
                        projection_found = 1;
                        min_distance = distance;
                        closest_index_padded = index_padded;
                      }
                    }
                  }
                }
              }
            }

            cubesize = (long)(cubesize*1.5);

            if (min_distance>(resolution_min*20)){
              // std::cout << "WARNING: Projection Synapses distance very large: "<< min_distance << std::endl;
            }

            if (min_distance>(resolution_min*30)){
              unable_to_project = 1;
              //std::cout << "Failed to find projection point, search distance too large" << std::endl;
            }

          }

          if (projection_found==1){
            // std::cout << "Projection found with cubesize: " << cubesize << std::endl;
            Block.Pointclouds[segment_ID][closest_index_padded] = 3;
            segment[closest_index_padded] = 3;

            IndexToIndices(closest_index_padded, ix_padded, iy_padded, iz_padded, padded_sheet_size_block, padded_row_size_block);

            iz_unpadded = iz_padded - 1;
            iy_unpadded = iy_padded - 1;
            ix_unpadded = ix_padded - 1;

            long projected_index = IndicesToIndex(ix_unpadded, iy_unpadded, iz_unpadded, input_sheet_size_block, input_row_size_block);
            Block.synapses[segment_ID].push_back(projected_index);
          }
        }
      }

    };


    static void NewSurfaceVoxel(long iv, long ix, long iy, long iz, List &surface_voxels)
    {
      ListElement *LE = new ListElement();
      LE->iv = iv;
      LE->ix = ix;
      LE->iy = iy;
      LE->iz = iz;

      LE->next = NULL;
      LE->prev = surface_voxels.last;

      if (surface_voxels.last != NULL) ((ListElement *) surface_voxels.last)->next = LE;
      surface_voxels.last = LE;
      if (surface_voxels.first == NULL) surface_voxels.first = LE;
    }

    static void RemoveSurfaceVoxel(ListElement *LE, List &surface_voxels)
    {
      ListElement *LE2;
      if (surface_voxels.first == LE) surface_voxels.first = LE->next;
      if (surface_voxels.last == LE) surface_voxels.last = LE->prev;

      if (LE->next != NULL) {
        LE2 = (ListElement *)(LE->next);
        LE2->prev = LE->prev;
      }
      if (LE->prev != NULL) {
        LE2 = (ListElement *)(LE->prev);
        LE2->next = LE->next;
      }
      delete LE;

    }

    static void CreatePointList(PointList *s)
    {
      s->head = NULL;
      s->tail = NULL;
      s->length = 0;
    }

    static void AddToList(PointList *s, Voxel e, ListElement *ptr)
    {
      Cell *newcell = new Cell();
      newcell->v = e;
      newcell->ptr = ptr;
      newcell->next = NULL;

      if (s->head == NULL) {
        s->head = newcell;
        s->tail = newcell;
        s->length = 1;
      }
      else {
        s->tail->next = newcell;
        s->tail = newcell;
        s->length++;
      }
    }

    static Voxel GetFromList(PointList *s, ListElement **ptr)
    {
      Voxel V;
      Cell *tmp;
      V.iv = -1;
      V.ix = -1;
      V.iy = -1;
      V.iz = -1;
      (*ptr) = NULL;
      if (s->length == 0) return V;
      else {
        V = s->head->v;
        (*ptr) = s->head->ptr;
        tmp = (Cell *) s->head->next;
        delete s->head;
        s->head = tmp;
        s->length--;
        if (s->length == 0) {
          s->head = NULL;
          s->tail = NULL;
        }
        return V;
      }
    }

    static void DestroyPointList(PointList *s)
    {
      ListElement *ptr;
      while (s->length) GetFromList(s, &ptr);
    }

    static void InitializeLookupTables(const char *lookup_table_directory)
    {
      char lut_filename[4096];
      FILE *lut_file;


      // read the simple lookup table
      sprintf(lut_filename, "%s/lut_simple.dat", lookup_table_directory);
      lut_simple = new unsigned char[lookup_table_size];
      lut_file = fopen(lut_filename, "rb");
      if (!lut_file) {
        fprintf(stderr, "Failed to read %s\n", lut_filename);
        exit(-1);
      }
      if (fread(lut_simple, 1, lookup_table_size, lut_file) != lookup_table_size) {
        fprintf(stderr, "Failed to read %s\n", lut_filename);
        exit(-1);
      }
      fclose(lut_file);

      // set the mask variables
      set_char_mask();
      set_long_mask();
    }

    static bool Simple26_6(unsigned int neighbors)
    {
      return lut_simple[(neighbors >> 3)] & char_mask[neighbors % 8];
    }

    void CPPcreateDataBlock(const char *prefix, const char *lookup_table_directory, long *inp_labels, long *inp_somae, float input_resolution[3],
      long inp_blocksize[3], long volume_size[3], long block_ind_inp[3], long block_ind_start_inp[3], long block_ind_end_inp[3],
      const char* synapses_dir, const char* output_dir){

        clock_t start_time_total = clock();

        // create new Datablock and set the input variables
        clock_t start_time_createDataBlock = clock();
        DataBlock*  BlockA = new DataBlock (prefix, input_resolution, inp_blocksize, volume_size, block_ind_inp, block_ind_start_inp, block_ind_end_inp, synapses_dir, output_dir);
        double time_createDataBlock = (double) (clock() - start_time_createDataBlock) / CLOCKS_PER_SEC;

        // process Somae
        clock_t start_time_popSomae = clock();
        BlockA->CppPopulateSomaeFromH5(inp_somae);
        double time_popSomae = (double) (clock() - start_time_popSomae) / CLOCKS_PER_SEC;

        // process the input labels
        clock_t start_time_popPointcloud = clock();
        BlockA->CppPopulatePointCloudFromH5(inp_labels);
        double time_popPointcloud = (double) (clock() - start_time_popPointcloud) / CLOCKS_PER_SEC;

        // read Synapses
        clock_t start_time_readSynapses = clock();
        if (!BlockA->ReadSynapses()) exit(-1);
        double time_readSynapses = (double) (clock() - start_time_readSynapses) / CLOCKS_PER_SEC;

        // read Anchor points (if block before does exist)
        clock_t start_time_readAnchors = clock();
        if (!BlockA->ReadAnchorpoints()) exit(-1);
        double time_readAnchors = (double) (clock() - start_time_readAnchors) / CLOCKS_PER_SEC;


        double time_setup = 0;
        double time_thinning = 0;
        double time_WriteOutput = 0;

        clock_t time_beforesetup = clock();
        clock_t time_beforethinning = clock();
        clock_t time_beforeWriteOutput = clock();


        // initialize lookup tables
        time_beforesetup = clock();
        InitializeLookupTables(lookup_table_directory);

        // needs to happen after PopulatePointCloud()
        PopulateOffsets(BlockA->padded_blocksize);

        // insert IDs that should be processed (IDs_in_Block if all)
        BlockA->IDs_to_process = BlockA->IDs_in_block;
        BlockA->writeIDs();
        time_setup += (double) (clock() - time_beforesetup) / CLOCKS_PER_SEC;

        // iterate over all IDs that should be processed
        uoSet::iterator itr = BlockA->IDs_to_process.begin();
        while (itr != BlockA->IDs_to_process.end())
        {
          // start timer for this segment
          clock_t start_time_seg = clock();

          time_beforesetup = clock();
          // initialize segment using the Block object and the segment ID to process
          BlockSegment* segA = new BlockSegment(*itr, *BlockA);
          time_setup += (double) (clock()-time_beforesetup) / CLOCKS_PER_SEC;

          // call the sequential thinning algorithm
          time_beforethinning = clock();
          segA->SequentialThinning(*BlockA);
          time_thinning += (double) (clock() - time_beforethinning) / CLOCKS_PER_SEC;

          // write skeletons and widths, add anchor points to
          time_beforeWriteOutput = clock();
          segA->WriteOutputfiles(*BlockA);
          // write timing to file
          segA->WriteTimeFile(start_time_seg);
          time_WriteOutput += (double) (clock() - time_beforeWriteOutput) / CLOCKS_PER_SEC;

          time_beforesetup = clock();
          delete segA;
          time_setup += (double) (clock()-time_beforesetup) / CLOCKS_PER_SEC;

          itr++;

        }

        time_beforeWriteOutput = clock();
        // write border anchor points to file
        BlockA->WriteProjectedSynapses();
        time_WriteOutput += (double) (clock() - time_beforeWriteOutput);

        double time_total = (double) (clock() - start_time_total) / CLOCKS_PER_SEC;

        time_beforesetup = clock();
        delete BlockA;
        time_setup += (double) (clock()-time_beforesetup) / CLOCKS_PER_SEC;


        {
          char output_filename[4096];
          sprintf(output_filename, "%s/running_times/%s/%s-total_time_thinning.pts", output_dir, prefix, prefix);
          //std::cout << "Writing time for block to : " << output_filename << std::endl;
          FILE * fptime = fopen (output_filename,"a");

          fprintf(fptime,"bz,by,bx\n");
          fprintf(fptime,"%04ld, %04ld, %04ld\n",
                  block_ind_inp[OR_Z], block_ind_inp[OR_Y], block_ind_inp[OR_X]);

          fprintf(fptime,"time_total, time_createDataBlock, time_popSomae, time_popPointcloud, time_readSynapses\n");
          fprintf(fptime,"%12.2f, %12.2f, %12.2f, %12.2f, %12.2f\n",
                  time_total, time_createDataBlock, time_popSomae, time_popPointcloud, time_readSynapses);

          fprintf(fptime,"time_readAnchors, time_setup, time_thinning, time_WriteOutput, time_projSynapses\n");
          fprintf(fptime,"%12.2f, %12.2f, %12.2f, %12.2f, %12.2f\n",
                  time_readAnchors, time_setup, time_thinning, time_WriteOutput, time_projSynapses);

          fclose(fptime);
        }

      }
