/* c++ file to extract wiring diagram */
#include <limits>
#include "cpp-wiring.h"
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <vector>
#include <time.h>
#include <stdio.h>
#include <string>
#include <unistd.h>
#include "cpp-MinBinaryHeap.h"
#include "cpp-wiring.h"
#include <algorithm>

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
static std::unordered_map<long, float> widths;

// mask variables for bitwise operations
static long long_mask[26];
static unsigned char char_mask[8];
static long n26_offsets[26];
static long n6_offsets[6];

double time_added = 0;

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
    n6_offsets[0] = -1 * padded_blocksize[OR_X];
    n6_offsets[1] = padded_blocksize[OR_X];
    n6_offsets[2] = -1 * padded_blocksize[OR_Y] * padded_blocksize[OR_X];
    n6_offsets[3] = padded_blocksize[OR_Y] * padded_blocksize[OR_X];
    n6_offsets[4] = +1;
    n6_offsets[5] = -1;
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

typedef struct {
    ListElement *first;
    ListElement *last;
} DoubleList;

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
    long padded_row_size = -1;
    long padded_sheet_size = -1;
    long input_row_size = -1;
    long input_sheet_size = -1;
    long volume_row_size = -1;
    long volume_sheet_size = -1;
    const char *synapses_directory;
    const char *somae_directory;
    const char *skeleton_directory;

  public:
    long padded_blocksize[3];
    std::unordered_set<long> IDs_to_process = std::unordered_set<long>();
    std::unordered_set<long> IDs_in_block = std::unordered_set<long>();
    std::unordered_map<long, std::unordered_map<long, char>> Pointclouds = std::unordered_map<long, std::unordered_map<long, char>>();
    std::unordered_map<long, std::unordered_map<long, std::unordered_set<long>>> borderpoints = std::unordered_map<long, std::unordered_map<long, std::unordered_set<long>>>();
    std::unordered_map<long, std::unordered_map<short, std::vector<long>>> max_anchors_comp = std::unordered_map<long, std::unordered_map<short, std::vector<long>>>();
    std::unordered_map<long, std::unordered_map<short, std::vector<long>>> min_anchors_seeded = std::unordered_map<long, std::unordered_map<short, std::vector<long>>>();
    std::unordered_map<long, std::vector<long>> synapses = std::unordered_map<long, std::vector<long>>();
    std::unordered_map<long, std::vector<long>> synapses_off = std::unordered_map<long, std::vector<long>>();

    DataBlock(float input_resolution[3], long inp_blocksize[3], long volume_size[3], long block_ind_inp[3], const char* synapses_dir, const char* somae_dir, const char* skeleton_dir){

      resolution[OR_Z] = input_resolution[OR_Z];
      resolution[OR_Y] = input_resolution[OR_Y];
      resolution[OR_X] = input_resolution[OR_X];

      std::cout << "Resolution set to: " << resolution[OR_Z] << "," << resolution[OR_Y] << "," << resolution[OR_X] << "," << std::endl;

      input_blocksize[OR_Z] = inp_blocksize[OR_Z];
      input_blocksize[OR_Y] = inp_blocksize[OR_Y];
      input_blocksize[OR_X] = inp_blocksize[OR_X];

      padded_blocksize[OR_Z] = inp_blocksize[OR_Z]+2;
      padded_blocksize[OR_Y] = inp_blocksize[OR_Y]+2;
      padded_blocksize[OR_X] = inp_blocksize[OR_X]+2;

      padded_sheet_size = padded_blocksize[OR_Y] * padded_blocksize[OR_X];
      padded_row_size = padded_blocksize[OR_X];

      input_sheet_size = input_blocksize[OR_Y] * input_blocksize[OR_X];
      input_row_size = input_blocksize[OR_X];

      std::cout << "Blocksize input set to: " << input_blocksize[OR_Z] << "," << input_blocksize[OR_Y] << "," << input_blocksize[OR_X] << "," << std::endl;
      std::cout << "Blocksize padded set to: " << padded_blocksize[OR_Z] << "," << padded_blocksize[OR_Y] << "," << padded_blocksize[OR_X] << "," << std::endl;
      std::cout << "Sheetsize padded set to: " << padded_sheet_size << std::endl;
      std::cout << "Rowsize padded set to: " << padded_row_size << std::endl;
      std::cout << "Sheetsize input set to: " << input_sheet_size << std::endl;
      std::cout << "Rowsize input set to: " << input_row_size << std::endl;

      volumesize[OR_Z] = volume_size[OR_Z];
      volumesize[OR_Y] = volume_size[OR_Y];
      volumesize[OR_X] = volume_size[OR_X];

      std::cout << "Volumesize set to: " << volumesize[OR_Z] << "," << volumesize[OR_Y] << "," << volumesize[OR_X] << "," << std::endl;

      block_ind[OR_Z]= block_ind_inp[OR_Z];
      block_ind[OR_Y] = block_ind_inp[OR_Y];
      block_ind[OR_X] = block_ind_inp[OR_X];

      std::cout << "Block indices set to: " << block_ind[OR_Z] << "," << block_ind[OR_Y] << "," << block_ind[OR_X] << "," << std::endl;

      synapses_directory = synapses_dir;
      somae_directory = somae_dir;
      skeleton_directory = skeleton_dir;

      std::cout << "Directories set. " << std::endl;

    }

    DataBlock(DataBlock &Block){
      std::copy(std::begin(Block.input_blocksize), std::end(Block.input_blocksize), std::begin(input_blocksize));
      std::copy(std::begin(Block.padded_blocksize), std::end(Block.padded_blocksize), std::begin(padded_blocksize));
      std::copy(std::begin(Block.volumesize), std::end(Block.volumesize), std::begin(volumesize));
      std::copy(std::begin(Block.resolution), std::end(Block.resolution), std::begin(resolution));
      std::copy(std::begin(Block.block_ind), std::end(Block.block_ind), std::begin(block_ind));
      skeleton_directory = Block.skeleton_directory;

      padded_row_size = Block.padded_row_size;
      padded_sheet_size = Block.padded_sheet_size;

      input_row_size = Block.input_row_size;
      input_sheet_size = Block.input_sheet_size;

    }

    void CppPopulatePointCloudFromH5(long *inp_labels)
    {

        long n_points = input_blocksize[0]*input_blocksize[1]*input_blocksize[2];

        std::cout << "n_points is: " << n_points << std::endl << std::flush;

        long z_max =input_blocksize[OR_Z]-1;
        long y_max =input_blocksize[OR_Y]-1;
        long x_max =input_blocksize[OR_X]-1;

        for (long voxel_index = 0; voxel_index < n_points; voxel_index++){

          // get segment_ID of current index and skip if is zero
          long curr_label = inp_labels[voxel_index];
          if (!curr_label) continue;

          // find the new voxel index
          long iv = PadIndex(voxel_index, input_sheet_size, input_row_size, padded_sheet_size, padded_row_size);

          // check if pointcloud of this segment_ID already exists, otherwise add new pointcloud
          if (Pointclouds.find(curr_label) == Pointclouds.end()) {
            Pointclouds[curr_label] = std::unordered_map<long,char>();
            // std::cout << "New segment_ID detected: " << curr_label << std::endl << std::flush;
            IDs_in_block.insert(curr_label);
          }

          Pointclouds[curr_label][iv] = 1;

          // add index to borderpoints unordered_map (only for max walls, as these anchor points are then copied to next level)
          long ix, iy, iz;
          IndexToIndices(voxel_index, ix, iy, iz, input_sheet_size, input_row_size);

          if ((iz==z_max && iy!=y_max) && ix!=x_max) borderpoints[curr_label][OR_Z].insert(iv);
          else if ((iz!=z_max && iy==y_max) && ix!=x_max) borderpoints[curr_label][OR_Y].insert(iv);
          else if ((iz!=z_max && iy!=y_max) && ix==x_max) borderpoints[curr_label][OR_X].insert(iv);

        }

    }

    int ReadSynapses(const char *prefix)
    {

      // read the synapses
      char synapse_filename[4096];
      snprintf(synapse_filename, 4096, "%s/%s/%s-synapses-%04ldz-%04ldy-%04ldx.pts", synapses_directory, prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

      FILE *fp = fopen(synapse_filename, "rb");
      if (!fp) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

      long nneurons;
      ReadHeader(fp, nneurons);

      for (long iv = 0; iv < nneurons; ++iv) {
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

              long linear_index;
              if (fread(&linear_index, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
              // add them to the synapses vector (unpadded)

              // find the new voxel index
              long iv = PadIndex(linear_index, input_sheet_size, input_row_size, padded_sheet_size, padded_row_size);

              // set ID of this index to 3 (==Synapse)
              // only add synapse if it is on the body
              if (Pointclouds[segment_ID][iv] == 1){
                Pointclouds[segment_ID][iv] = 3;
                synapses[segment_ID].push_back(linear_index);
              }
              else {
                synapses_off[segment_ID].push_back(linear_index);
              }

          }

      }


      // close file
      fclose(fp);

      return 1;

    }

    int ReadAnchorpoints(const char *prefix)
    {

      // read in z anchors
      if (block_ind[OR_Z]>0)
      {
        // make filename of adjacent z lock (in negative direction)
        char output_filename_zmax[4096];
        sprintf(output_filename_zmax, "%s/%s/%s-Anchors_Comp_Z-%04ldz-%04ldy-%04ldx.pts", skeleton_directory, prefix, prefix, block_ind[OR_Z]-1, block_ind[OR_Y], block_ind[OR_X]);

        FILE *fpzmax = fopen(output_filename_zmax, "rb");
        if (!fpzmax) { fprintf(stderr, "Failed to read %s.\n", output_filename_zmax); return 0; }

        long nsegments;
        ReadHeader(fpzmax, nsegments);

        for (long i=0; i<nsegments; i++) {

          long seg_ID;
          long n_anchors;

          if (fread(&seg_ID, sizeof(long), 1, fpzmax) != 1) { fprintf(stderr, "Failed to read %s\n", output_filename_zmax); exit(-1); }
          if (fread(&n_anchors, sizeof(long), 1, fpzmax) != 1) { fprintf(stderr, "Failed to read %s\n", output_filename_zmax); exit(-1); }

          for (long pos=0; pos<n_anchors; pos++) {

            long iv_local;
            if (fread(&iv_local, sizeof(long), 1, fpzmax) != 1) { fprintf(stderr, "Failed to read %s\n", output_filename_zmax); exit(-1); }

            long iz, iy, ix;
            IndexToIndices(iv_local, ix, iy, iz, input_sheet_size, input_row_size);

            long iy_padded = iy + 1;
            long ix_padded = ix + 1;
            long iz_padded;

            long iv_padded;
            long iv_unpadded;

            iz_padded = 1;
            iz = 0;
            iv_padded = IndicesToIndex(ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);
            iv_unpadded = IndicesToIndex(ix, iy, iz, input_sheet_size, input_row_size);
            Pointclouds[seg_ID][iv_padded] = 3;
            min_anchors_seeded[seg_ID][OR_Z].push_back(iv_unpadded);

            iz_padded = 2;
            iz = 1;
            iv_padded = IndicesToIndex(ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);
            iv_unpadded = IndicesToIndex(ix, iy, iz, input_sheet_size, input_row_size);
            Pointclouds[seg_ID][iv_padded] = 3;
            min_anchors_seeded[seg_ID][OR_Z].push_back(iv_unpadded);

            iz_padded = 3;
            iz = 2;
            iv_padded = IndicesToIndex(ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);
            iv_unpadded = IndicesToIndex(ix, iy, iz, input_sheet_size, input_row_size);
            Pointclouds[seg_ID][iv_padded] = 3;
            min_anchors_seeded[seg_ID][OR_Z].push_back(iv_unpadded);
          }
        }

        // close file
        fclose(fpzmax);
      }

      // read in y anchors
      if (block_ind[OR_Y]>0)
      {
        // make filename of adjacent z lock (in negative direction)
        char output_filename_zmax[4096];
        sprintf(output_filename_zmax, "%s/%s/%s-Anchors_Comp_Y-%04ldz-%04ldy-%04ldx.pts", skeleton_directory, prefix, prefix, block_ind[OR_Z], block_ind[OR_Y]-1, block_ind[OR_X]);

        FILE *fpzmax = fopen(output_filename_zmax, "rb");
        if (!fpzmax) { fprintf(stderr, "Failed to read %s.\n", output_filename_zmax); return 0; }

        long nsegments;
        ReadHeader(fpzmax, nsegments);

        for (long i=0; i<nsegments; i++) {

          long seg_ID;
          long n_anchors;

          if (fread(&seg_ID, sizeof(long), 1, fpzmax) != 1) { fprintf(stderr, "Failed to read %s\n", output_filename_zmax); exit(-1); }
          if (fread(&n_anchors, sizeof(long), 1, fpzmax) != 1) { fprintf(stderr, "Failed to read %s\n", output_filename_zmax); exit(-1); }

          for (long pos=0; pos<n_anchors; pos++) {

            long iv_local;
            if (fread(&iv_local, sizeof(long), 1, fpzmax) != 1) { fprintf(stderr, "Failed to read %s\n", output_filename_zmax); exit(-1); }

            long iz, iy, ix;
            IndexToIndices(iv_local, ix, iy, iz, input_sheet_size, input_row_size);

            long iz_padded = iz + 1;
            long ix_padded = ix + 1;
            long iy_padded;

            long iv_padded;
            long iv_unpadded;

            iy_padded = 1;
            iy = 0;
            iv_padded = IndicesToIndex(ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);
            iv_unpadded = IndicesToIndex(ix, iy, iz, input_sheet_size, input_row_size);
            Pointclouds[seg_ID][iv_padded] = 3;
            min_anchors_seeded[seg_ID][OR_Y].push_back(iv_unpadded);

            iy_padded = 2;
            iy = 1;
            iv_padded = IndicesToIndex(ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);
            iv_unpadded = IndicesToIndex(ix, iy, iz, input_sheet_size, input_row_size);
            Pointclouds[seg_ID][iv_padded] = 3;
            min_anchors_seeded[seg_ID][OR_Y].push_back(iv_unpadded);

            iy_padded = 3;
            iy = 2;
            iv_padded = IndicesToIndex(ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);
            iv_unpadded = IndicesToIndex(ix, iy, iz, input_sheet_size, input_row_size);
            Pointclouds[seg_ID][iv_padded] = 3;
            min_anchors_seeded[seg_ID][OR_Y].push_back(iv_unpadded);

          }
        }
        // close file
        fclose(fpzmax);
      }

      // read in y anchors
      if (block_ind[OR_X]>0)
      {
        // make filename of adjacent z lock (in negative direction)
        char output_filename_zmax[4096];
        sprintf(output_filename_zmax, "%s/%s/%s-Anchors_Comp_X-%04ldz-%04ldy-%04ldx.pts", skeleton_directory, prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]-1);

        FILE *fpzmax = fopen(output_filename_zmax, "rb");
        if (!fpzmax) { fprintf(stderr, "Failed to read %s.\n", output_filename_zmax); return 0; }

        long nsegments;
        ReadHeader(fpzmax, nsegments);

        for (long i=0; i<nsegments; i++) {

          long seg_ID;
          long n_anchors;

          if (fread(&seg_ID, sizeof(long), 1, fpzmax) != 1) { fprintf(stderr, "Failed to read %s\n", output_filename_zmax); exit(-1); }
          if (fread(&n_anchors, sizeof(long), 1, fpzmax) != 1) { fprintf(stderr, "Failed to read %s\n", output_filename_zmax); exit(-1); }

          for (long pos=0; pos<n_anchors; pos++) {

            long iv_local;
            if (fread(&iv_local, sizeof(long), 1, fpzmax) != 1) { fprintf(stderr, "Failed to read %s\n", output_filename_zmax); exit(-1); }

            long iz, iy, ix;
            IndexToIndices(iv_local, ix, iy, iz, input_sheet_size, input_row_size);

            long iy_padded = iy + 1;
            long iz_padded = iz + 1;

            long iv_padded;
            long iv_unpadded;
            long ix_padded;

            ix_padded = 1;
            ix = 0;
            iv_padded = IndicesToIndex(ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);
            iv_unpadded = IndicesToIndex(ix, iy, iz, input_sheet_size, input_row_size);
            Pointclouds[seg_ID][iv_padded] = 3;
            min_anchors_seeded[seg_ID][OR_X].push_back(iv_unpadded);

            ix_padded = 2;
            ix = 1;
            iv_padded = IndicesToIndex(ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);
            iv_unpadded = IndicesToIndex(ix, iy, iz, input_sheet_size, input_row_size);
            Pointclouds[seg_ID][iv_padded] = 3;
            min_anchors_seeded[seg_ID][OR_X].push_back(iv_unpadded);

            ix_padded = 3;
            ix = 2;
            iv_padded = IndicesToIndex(ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);
            iv_unpadded = IndicesToIndex(ix, iy, iz, input_sheet_size, input_row_size);
            Pointclouds[seg_ID][iv_padded] = 3;
            min_anchors_seeded[seg_ID][OR_X].push_back(iv_unpadded);
          }
        }

      // close file
      fclose(fpzmax);
    }

      return 1;

    }

    void writeAnchorsComputed (const char *prefix)
    {
      // write the zmax anchor points to a file (so far stored in vectors)
      {
        char output_filename_zmax[4096];
        sprintf(output_filename_zmax, "%s/%s/%s-Anchors_Comp_Z-%04ldz-%04ldy-%04ldx.pts", skeleton_directory, prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

        FILE *zmaxfp = fopen(output_filename_zmax, "wb");
        if (!zmaxfp) { fprintf(stderr, "Failed to open %s\n", output_filename_zmax); exit(-1); }

        long nsegments = max_anchors_comp.size();
        WriteHeader(zmaxfp, nsegments);

        for (std::unordered_map<long, std::unordered_map<short, std::vector<long>>>::iterator iter = max_anchors_comp.begin(); iter != max_anchors_comp.end(); ++iter) {

          long seg_ID = iter->first;

          long n_anchors;

          // write z anchors
          n_anchors =  max_anchors_comp[seg_ID][OR_Z].size();
          if (fwrite(&seg_ID, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
          if (fwrite(&n_anchors, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
          for (long pos=0; pos<n_anchors; pos++) {
            long iv_local = max_anchors_comp[seg_ID][OR_Z][pos];
            if (fwrite(&iv_local, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
          }
        }
        fclose(zmaxfp);

      }
      // write the ymax anchor points to a file (so far stored in vectors)
      {
        char output_filename_zmax[4096];
        sprintf(output_filename_zmax, "%s/%s/%s-Anchors_Comp_Y-%04ldz-%04ldy-%04ldx.pts", skeleton_directory, prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

        FILE *zmaxfp = fopen(output_filename_zmax, "wb");
        if (!zmaxfp) { fprintf(stderr, "Failed to open %s\n", output_filename_zmax); exit(-1); }

        long nsegments = max_anchors_comp.size();
        WriteHeader(zmaxfp, nsegments);

        for (std::unordered_map<long, std::unordered_map<short, std::vector<long>>>::iterator iter = max_anchors_comp.begin(); iter != max_anchors_comp.end(); ++iter) {

          long seg_ID = iter->first;

          long n_anchors;

          // write z anchors
          n_anchors =  max_anchors_comp[seg_ID][OR_Y].size();
          if (fwrite(&seg_ID, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
          if (fwrite(&n_anchors, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
          for (long pos=0; pos<n_anchors; pos++) {
            long iv_local = max_anchors_comp[seg_ID][OR_Y][pos];
            if (fwrite(&iv_local, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
          }
        }
        fclose(zmaxfp);
      }
      // write the xmax anchor points to a file (so far stored in vectors)
      {
        char output_filename_zmax[4096];
        sprintf(output_filename_zmax, "%s/%s/%s-Anchors_Comp_X-%04ldz-%04ldy-%04ldx.pts", skeleton_directory, prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

        FILE *zmaxfp = fopen(output_filename_zmax, "wb");
        if (!zmaxfp) { fprintf(stderr, "Failed to open %s\n", output_filename_zmax); exit(-1); }

        long nsegments = max_anchors_comp.size();
        WriteHeader(zmaxfp, nsegments);

        for (std::unordered_map<long, std::unordered_map<short, std::vector<long>>>::iterator iter = max_anchors_comp.begin(); iter != max_anchors_comp.end(); ++iter) {

          long seg_ID = iter->first;

          long n_anchors;

          // write z anchors
          n_anchors =  max_anchors_comp[seg_ID][OR_X].size();
          if (fwrite(&seg_ID, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
          if (fwrite(&n_anchors, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
          for (long pos=0; pos<n_anchors; pos++) {
            long iv_local = max_anchors_comp[seg_ID][OR_X][pos];
            if (fwrite(&iv_local, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
          }

        }
        fclose(zmaxfp);
      }
    }

    void writeIDsToProcess(const char *prefix)
    {
      char output_filename_IDs[4096];
      sprintf(output_filename_IDs, "%s/%s/%s-IDsProcessed-%04ldz-%04ldy-%04ldx.pts", skeleton_directory, prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

      FILE *fpid = fopen(output_filename_IDs, "wb");
      if (!fpid) { fprintf(stderr, "Failed to open %s\n", output_filename_IDs); exit(-1); }

      long num = IDs_to_process.size();
      WriteHeader(fpid, num);

      for (std::unordered_set<long>::iterator iter = IDs_to_process.begin(); iter != IDs_to_process.end(); ++iter) {

        long seg_ID = *iter;

        // IDs
        if (fwrite(&seg_ID, sizeof(long), 1, fpid) != 1) { fprintf(stderr, "Failed to write to IDs in Block to %s \n", output_filename_IDs); exit(-1); }
      }
      fclose(fpid);
    }

    void writeanchorsSeeded (const char *prefix)
    {
      // write the zmin anchors that were seeded in this block
      {
        char output_filename_zmax[4096];
        sprintf(output_filename_zmax, "%s/%s/%s-Anchors_Seeded_Z-%04ldz-%04ldy-%04ldx.pts", skeleton_directory, prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

        FILE *zmaxfp = fopen(output_filename_zmax, "wb");
        if (!zmaxfp) { fprintf(stderr, "Failed to open %s\n", output_filename_zmax); exit(-1); }

        long nsegments = min_anchors_seeded.size();
        WriteHeader(zmaxfp, nsegments);

        for (std::unordered_map<long, std::unordered_map<short, std::vector<long>>>::iterator iter = min_anchors_seeded.begin(); iter != min_anchors_seeded.end(); ++iter) {

        long seg_ID = iter->first;

        long n_anchors;

        // write z anchors
        n_anchors =  min_anchors_seeded[seg_ID][OR_Z].size();
        if (fwrite(&seg_ID, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
        if (fwrite(&n_anchors, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
        for (long pos=0; pos<n_anchors; pos++) {
          long iv_local = min_anchors_seeded[seg_ID][OR_Z][pos];
          if (fwrite(&iv_local, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
        }

        }
        fclose(zmaxfp);

      }
      // write the ymin anchors that were seeded in this block
      {
        char output_filename_zmax[4096];
        sprintf(output_filename_zmax, "%s/%s/%s-Anchors_Seeded_Y-%04ldz-%04ldy-%04ldx.pts", skeleton_directory, prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

        FILE *zmaxfp = fopen(output_filename_zmax, "wb");
        if (!zmaxfp) { fprintf(stderr, "Failed to open %s\n", output_filename_zmax); exit(-1); }

        long nsegments = min_anchors_seeded.size();
        WriteHeader(zmaxfp, nsegments);

        for (std::unordered_map<long, std::unordered_map<short, std::vector<long>>>::iterator iter = min_anchors_seeded.begin(); iter != min_anchors_seeded.end(); ++iter) {

          long seg_ID = iter->first;

          long n_anchors;

          // write z anchors
          n_anchors =  min_anchors_seeded[seg_ID][OR_Y].size();
          if (fwrite(&seg_ID, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
          if (fwrite(&n_anchors, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
          for (long pos=0; pos<n_anchors; pos++) {
            long iv_local = min_anchors_seeded[seg_ID][OR_Y][pos];
            if (fwrite(&iv_local, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
          }
        }
        fclose(zmaxfp);
      }
      // write the xmin anchors that were seeded in this block
      {
        char output_filename_zmax[4096];
        sprintf(output_filename_zmax, "%s/%s/%s-Anchors_Seeded_X-%04ldz-%04ldy-%04ldx.pts", skeleton_directory, prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

        FILE *zmaxfp = fopen(output_filename_zmax, "wb");
        if (!zmaxfp) { fprintf(stderr, "Failed to open %s\n", output_filename_zmax); exit(-1); }

        long nsegments = min_anchors_seeded.size();
        WriteHeader(zmaxfp, nsegments);

        for (std::unordered_map<long, std::unordered_map<short, std::vector<long>>>::iterator iter = min_anchors_seeded.begin(); iter != min_anchors_seeded.end(); ++iter) {

        long seg_ID = iter->first;

        long n_anchors;

        // write z anchors
        n_anchors =  min_anchors_seeded[seg_ID][OR_X].size();
        if (fwrite(&seg_ID, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
        if (fwrite(&n_anchors, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
        for (long pos=0; pos<n_anchors; pos++) {
          long iv_local = min_anchors_seeded[seg_ID][OR_X][pos];
          if (fwrite(&iv_local, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
        }


        }
        fclose(zmaxfp);
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

    void WriteProjectedSynapses(const char *prefix)
    {
        //get number of anchor points
        long n_neurons = IDs_to_process.size();

        // create an output file for the points
        char output_filename[4096];
        sprintf(output_filename, "%s/%s/%s-projectedSynapses-%04ldz-%04ldy-%04ldx.pts", synapses_directory, prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

        FILE *wfp = fopen(output_filename, "wb");
        if (!wfp) { fprintf(stderr, "Failed to open %s\n", output_filename); exit(-1); }

        // write the characteristics header
        WriteHeader(wfp, n_neurons);

        for (std::unordered_map<long,std::vector<long>>::iterator itr = synapses.begin(); itr!=synapses.end(); ++itr){
            long seg_id = itr->first;
            long n_synapses = synapses[seg_id].size();

            if (fwrite(&seg_id, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
            if (fwrite(&n_synapses, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }

            long index_local[n_synapses];
            long pos = 0;

            for (std::vector<long>::iterator itr2 = synapses[seg_id].begin(); itr2!=synapses[seg_id].end(); ++itr2, ++pos){

              long iv_local = *itr2;

              long iz_local, iy_local, ix_local;
              IndexToIndices(iv_local, ix_local, iy_local, iz_local, input_sheet_size, input_row_size);

              long iz_global = iz_local + block_ind[OR_Z]*input_blocksize[OR_Z];
              long iy_global = iy_local + block_ind[OR_Y]*input_blocksize[OR_Y];
              long ix_global = ix_local + block_ind[OR_X]*input_blocksize[OR_X];

              long iv_global = IndicesToIndex(ix_global, iy_global, iz_global, volume_sheet_size, volume_row_size);
              if (fwrite(&iv_global, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }

              index_local[pos] = iv_local;
            }

            for (int j=0; j<n_synapses; j++){
              if (fwrite(&index_local[j], sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
            }

        }
    }

};

class BlockSegment : public DataBlock{
  private:
    long initial_points = -1;
    long segment_ID = -1;
    List surface_voxels;
    long padded_infinity = -1;
    std::unordered_map<long, char> segment = std::unordered_map<long, char>();
    std::unordered_map<long, std::unordered_set<long>> borderpoints_segment = std::unordered_map<long, std::unordered_set<long>>();
    std::unordered_map<long, float> widths = std::unordered_map<long, float>();
    std::unordered_map<long, char> skeleton = std::unordered_map<long, char>();

  public:
    BlockSegment(long segment_ID_inp, DataBlock &Block):DataBlock(Block){

      segment_ID = segment_ID_inp;
      segment = Block.Pointclouds[segment_ID];
      borderpoints_segment = Block.borderpoints[segment_ID];

      std::cout << "-----------------------------------" << std::endl;
      std::cout << "Processing segment_ID " << segment_ID << std::endl;

      initial_points = segment.size();
      printf("segment_ID %ld initial points: %ld\n", segment_ID, initial_points);

      surface_voxels.first = NULL;
      surface_voxels.last = NULL;

    }

    void SequentialThinning(const char *prefix, DataBlock &Block)
    {
        // create a vector of surface voxels
        CollectSurfaceVoxels();
        Projectsynapses(Block);
        int iteration = 0;
        long changed = 0;
        do {
            changed = ThinningIterationStep(Block);
            iteration++;
            // printf("  Iteration %d deleted %ld points\n", iteration, changed);
        } while (changed);
        printf("Needed %d iterations\n", iteration);
    }

    void CollectSurfaceVoxels(void)
    {

        // go through all voxels and check their six neighbors
        for (std::unordered_map<long, char>::iterator it = segment.begin(); it != segment.end(); ++it) {

            // all of these elements are either 1 or 3 and in the segment
            long index = it->first;

            // initialize widths to maximum float value
            widths[index] = std::numeric_limits<float>::max();

            long ix, iy, iz;
            IndexToIndices(index, ix, iy, iz, padded_sheet_size, padded_row_size);
            // check the 6 neighbors
            for (long iv = 0; iv < NTHINNING_DIRECTIONS; ++iv) {
                long neighbor_index = index + n6_offsets[iv];

                long ii, ij, ik;

                IndexToIndices(neighbor_index, ii, ij, ik, padded_sheet_size, padded_row_size);

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
                    }

                    // note this location as surface
                    widths[index] = 0;

                    break;
                }
            }
        }
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
                if ((borderpoints_segment[OR_Y].find(index)!=borderpoints_segment[OR_Y].end()) && direction==0) continue;
                if ((borderpoints_segment[OR_Z].find(index)!=borderpoints_segment[OR_Z].end()) && direction==2) continue;
                if ((borderpoints_segment[OR_X].find(index)!=borderpoints_segment[OR_X].end()) && direction==5) continue;

                bool isAnchorPoint = 0;

                if (borderpoints_segment[OR_Y].find(index)!=borderpoints_segment[OR_Y].end()){

                    long sum_of_neighbors = 0; //collect voxels of neighbors on x-z plane
                    sum_of_neighbors += segment[index+ n26_offsets[3]];
                    sum_of_neighbors += segment[index+ n26_offsets[4]];
                    sum_of_neighbors += segment[index+ n26_offsets[5]];
                    sum_of_neighbors += segment[index+ n26_offsets[12]];
                    sum_of_neighbors += segment[index+ n26_offsets[13]];
                    sum_of_neighbors += segment[index+ n26_offsets[20]];
                    sum_of_neighbors += segment[index+ n26_offsets[21]];
                    sum_of_neighbors += segment[index+ n26_offsets[22]];

                    if (sum_of_neighbors==0) {
                      isAnchorPoint = 1;

                      // fix anchor point (as a fake synapse) TODO: might need another label here to identify anchor points seperately
                      long ix_local, iy_local, iz_local;
                      IndexToIndices(index, ix_local, iy_local, iz_local, padded_sheet_size, padded_row_size);
                      ix_local-=1; iy_local-=1; iz_local-=1;

                      long iv_local = IndicesToIndex(ix_local,iy_local,iz_local,input_sheet_size, input_row_size);

                      // write to zmax anchor file
                      Block.max_anchors_comp[segment_ID][OR_Y].push_back(iv_local);

                    }
                }
                else if (borderpoints_segment[OR_Z].find(index)!=borderpoints_segment[OR_Z].end()){

                    long sum_of_neighbors = 0; //collect voxels of neighbors on x-y plane
                    sum_of_neighbors += segment[index+ n26_offsets[9]];
                    sum_of_neighbors += segment[index+ n26_offsets[10]];
                    sum_of_neighbors += segment[index+ n26_offsets[11]];
                    sum_of_neighbors += segment[index+ n26_offsets[12]];
                    sum_of_neighbors += segment[index+ n26_offsets[13]];
                    sum_of_neighbors += segment[index+ n26_offsets[14]];
                    sum_of_neighbors += segment[index+ n26_offsets[15]];
                    sum_of_neighbors += segment[index+ n26_offsets[16]];

                    if (sum_of_neighbors==0) {
                      isAnchorPoint = 1;

                      // fix anchor point (as a fake synapse) TODO: might need another label here to identify anchor points seperately
                      long ix_local, iy_local, iz_local;
                      IndexToIndices(index, ix_local, iy_local, iz_local, padded_sheet_size, padded_row_size);
                      ix_local-=1; iy_local-=1; iz_local-=1;

                      long iv_local = IndicesToIndex(ix_local,iy_local,iz_local,input_sheet_size, input_row_size);

                      // write to zmax anchor file
                      Block.max_anchors_comp[segment_ID][OR_Z].push_back(iv_local);

                    }
                }
                else if (borderpoints_segment[OR_X].find(index)!=borderpoints_segment[OR_X].end()){

                    long sum_of_neighbors = 0; //collect voxels of neighbors on y-z plane
                    sum_of_neighbors += segment[index+ n26_offsets[7]];
                    sum_of_neighbors += segment[index+ n26_offsets[4]];
                    sum_of_neighbors += segment[index+ n26_offsets[1]];
                    sum_of_neighbors += segment[index+ n26_offsets[10]];
                    sum_of_neighbors += segment[index+ n26_offsets[15]];
                    sum_of_neighbors += segment[index+ n26_offsets[24]];
                    sum_of_neighbors += segment[index+ n26_offsets[21]];
                    sum_of_neighbors += segment[index+ n26_offsets[18]];

                    if (sum_of_neighbors==0) {
                      isAnchorPoint = 1;

                      // fix anchor point (as a fake synapse) TODO: might need another label here to identify anchor points seperately
                      long ix_local, iy_local, iz_local;
                      IndexToIndices(index, ix_local, iy_local, iz_local, padded_sheet_size, padded_row_size);
                      ix_local-=1; iy_local-=1; iz_local-=1;

                      long iv_local = IndicesToIndex(ix_local,iy_local,iz_local,input_sheet_size, input_row_size);

                      // write to zmax anchor file
                      Block.max_anchors_comp[segment_ID][OR_X].push_back(iv_local);

                    }
                }

                // if anchor point detected, fix it, do not check if simple and remove from surface voxels
                if (isAnchorPoint){
                  segment[index]=3;
                  RemoveSurfaceVoxel(ptr, surface_voxels);
                  continue;
                }

                // otherise do normal procedure - ceck if simple, if so, delete it
                else{
                  unsigned int neighbors = Collect26Neighbors(ix, iy, iz);
                  if (Simple26_6(neighbors)) {
                      // delete the simple point
                      segment[index] = 0;

                      // add the new surface voxels
                      for (long ip = 0; ip < NTHINNING_DIRECTIONS; ++ip) {
                          long neighbor_index = index + n6_offsets[ip];

                          // previously not on the surface but is in the object
                          // widths of voxels start at maximum and first updated when put on surface
                          if (segment[neighbor_index] == 1) {
                              long iu, iv, iw;
                              IndexToIndices(neighbor_index, iu, iv, iw, padded_sheet_size, padded_row_size);
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
                          IndexToIndices(neighbor_index, iu, iv, iw, padded_sheet_size, padded_row_size);

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
            }

            DestroyPointList(&deletable_points);

        }

        // return the number of changes
        return changed;
    }

    unsigned int Collect26Neighbors(long ix, long iy, long iz)
    {
        unsigned int neighbors = 0;

        long index = IndicesToIndex(ix, iy, iz, padded_sheet_size, padded_row_size);

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
            // this will only be called on things on the surface already so already in unordered_map
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

    void WriteOutputfiles(const char *prefix, DataBlock &Block)
    {

        // count the number of remaining points
        long num = 0;
        ListElement *LE = (ListElement *) surface_voxels.first;
        while (LE != NULL) {
            num++;
            LE = (ListElement *)LE->next;
        }

        //get number of anchor points
        long n_anchors_comp_z =  Block.max_anchors_comp[segment_ID][OR_Z].size();
        long n_anchors_seeded_z =  Block.min_anchors_seeded[segment_ID][OR_Z].size();
        long n_anchors_comp_y =  Block.max_anchors_comp[segment_ID][OR_Y].size();
        long n_anchors_seeded_y =  Block.min_anchors_seeded[segment_ID][OR_Y].size();
        long n_anchors_comp_x =  Block.max_anchors_comp[segment_ID][OR_X].size();
        long n_anchors_seeded_x =  Block.min_anchors_seeded[segment_ID][OR_X].size();
        long n_synapses =  Block.synapses[segment_ID].size();


        // create an output file for the points
        char output_filename[4096];
        sprintf(output_filename, "%s/%s/%s-skeleton-%04ldz-%04ldy-%04ldx-ID-%012ld.pts", skeleton_directory, prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X] ,segment_ID);

        // write the widths to file
        char widths_filename[4096];
        sprintf(widths_filename, "widths/%s/%06ld.pts", prefix, segment_ID);

        FILE *wfp = fopen(output_filename, "wb");
        if (!wfp) { fprintf(stderr, "Failed to open %s\n", output_filename); exit(-1); }

        FILE *width_fp = fopen(widths_filename, "wb");
        if (!width_fp) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }

        long total_points = num+n_anchors_comp_z+n_anchors_seeded_z+n_anchors_comp_y+n_anchors_seeded_y+n_anchors_comp_x+n_anchors_seeded_x+n_synapses;

        // write the characteristics header
        WriteHeaderSegID(wfp, total_points);

        // characteristics
        WriteHeaderSegID(width_fp, total_points);

        printf("Remaining voxels: %ld\n", num);

        long index_local[num];
        long counter_it = 0;

        // write surface voxels (global index)
        while (surface_voxels.first != NULL) {
            // get the surface voxels
            ListElement *LE = (ListElement *) surface_voxels.first;

            long index = LE->iv;
            float width = widths[index];

            // get the coordinates for this skeleton point in the global frame
            long iz_local = LE->iz - 1;
            long iy_local = LE->iy - 1;
            long ix_local = LE->ix - 1;

            long iz_global = iz_local + block_ind[OR_Z]*input_blocksize[OR_Z];
            long iy_global = iy_local + block_ind[OR_Y]*input_blocksize[OR_Y];
            long ix_global = ix_local + block_ind[OR_X]*input_blocksize[OR_X];

            long iv_local = IndicesToIndex(ix_local, iy_local, iz_local, input_sheet_size, input_row_size);
            long iv_global = IndicesToIndex(ix_global, iy_global, iz_global, volume_sheet_size, volume_row_size);

            if (fwrite(&iv_global, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
            if (fwrite(&iv_global, sizeof(long), 1, width_fp) != 1) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
            if (fwrite(&width, sizeof(float), 1, width_fp) != 1)  { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }

            index_local[counter_it]=iv_local;

            // remove this voxel
            RemoveSurfaceVoxel(LE, surface_voxels);
            counter_it++;
        }

        // write z anchor points computed in this step (global index)
        long index_local_anchors_comp_z[n_anchors_comp_z];
        for (long pos=0; pos<n_anchors_comp_z; pos++) {

          long iv_local = Block.max_anchors_comp[segment_ID][OR_Z][pos];

          long iz_local, iy_local, ix_local;
          IndexToIndices(iv_local, ix_local, iy_local, iz_local, input_sheet_size, input_row_size);


          long iz_padded = iz_local + 1;
          long iy_padded = iy_local + 1;
          long ix_padded = ix_local + 1;

          long iv_padded = IndicesToIndex(ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);
          float width = widths[iv_padded];

          long iz_global = iz_local + block_ind[OR_Z]*input_blocksize[OR_Z];
          long iy_global = iy_local + block_ind[OR_Y]*input_blocksize[OR_Y];
          long ix_global = ix_local + block_ind[OR_X]*input_blocksize[OR_X];

          long iv_global = IndicesToIndex(ix_global, iy_global, iz_global, volume_sheet_size, volume_row_size);

          if (fwrite(&iv_global, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
          if (fwrite(&iv_global, sizeof(long), 1, width_fp) != 1) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
          if (fwrite(&width, sizeof(float), 1, width_fp) != 1)  { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }

          index_local_anchors_comp_z[pos] = iv_local;

        }

        // write z anchor points seeded in this step (global index)
        long index_local_anchors_seeded_z[n_anchors_seeded_z];
        for (long pos=0; pos<n_anchors_seeded_z; pos++) {

          long iv_local = Block.min_anchors_seeded[segment_ID][OR_Z][pos];



          long iz_local, iy_local, ix_local;
          IndexToIndices(iv_local, ix_local, iy_local, iz_local, input_sheet_size, input_row_size);


          long iz_padded = iz_local + 1;
          long iy_padded = iy_local + 1;
          long ix_padded = ix_local + 1;

          long iv_padded = IndicesToIndex(ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);
          float width = widths[iv_padded];

          long iz_global = iz_local + block_ind[OR_Z]*input_blocksize[OR_Z];
          long iy_global = iy_local + block_ind[OR_Y]*input_blocksize[OR_Y];
          long ix_global = ix_local + block_ind[OR_X]*input_blocksize[OR_X];

          long iv_global = IndicesToIndex(ix_global, iy_global, iz_global, volume_sheet_size, volume_row_size);

          if (fwrite(&iv_global, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
          if (fwrite(&iv_global, sizeof(long), 1, width_fp) != 1) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
          if (fwrite(&width, sizeof(float), 1, width_fp) != 1)  { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }

          index_local_anchors_seeded_z[pos] = iv_local;

        }

        // write y anchor points computed in this step (global index)
        long index_local_anchors_comp_y[n_anchors_comp_y];
        for (long pos=0; pos<n_anchors_comp_y; pos++) {

          long iv_local = Block.max_anchors_comp[segment_ID][OR_Y][pos];

          long iz_local, iy_local, ix_local;
          IndexToIndices(iv_local, ix_local, iy_local, iz_local, input_sheet_size, input_row_size);

          long iz_padded = iz_local + 1;
          long iy_padded = iy_local + 1;
          long ix_padded = ix_local + 1;

          long iv_padded = IndicesToIndex(ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);
          float width = widths[iv_padded];

          long iz_global = iz_local + block_ind[OR_Z]*input_blocksize[OR_Z];
          long iy_global = iy_local + block_ind[OR_Y]*input_blocksize[OR_Y];
          long ix_global = ix_local + block_ind[OR_X]*input_blocksize[OR_X];

          long iv_global = IndicesToIndex(ix_global, iy_global, iz_global, volume_sheet_size, volume_row_size);

          if (fwrite(&iv_global, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
          if (fwrite(&iv_global, sizeof(long), 1, width_fp) != 1) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
          if (fwrite(&width, sizeof(float), 1, width_fp) != 1)  { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }

          index_local_anchors_comp_y[pos] = iv_local;

        }

        // write y anchor points seeded in this step (global index)
        long index_local_anchors_seeded_y[n_anchors_seeded_y];
        for (long pos=0; pos<n_anchors_seeded_y; pos++) {

          long iv_local = Block.min_anchors_seeded[segment_ID][OR_Y][pos];

          long iz_local, iy_local, ix_local;
          IndexToIndices(iv_local, ix_local, iy_local, iz_local, input_sheet_size, input_row_size);

          long iz_padded = iz_local + 1;
          long iy_padded = iy_local + 1;
          long ix_padded = ix_local + 1;

          long iv_padded = IndicesToIndex(ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);
          float width = widths[iv_padded];

          long iz_global = iz_local + block_ind[OR_Z]*input_blocksize[OR_Z];
          long iy_global = iy_local + block_ind[OR_Y]*input_blocksize[OR_Y];
          long ix_global = ix_local + block_ind[OR_X]*input_blocksize[OR_X];

          long iv_global = IndicesToIndex(ix_global, iy_global, iz_global, volume_sheet_size, volume_row_size);

          if (fwrite(&iv_global, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
          if (fwrite(&iv_global, sizeof(long), 1, width_fp) != 1) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
          if (fwrite(&width, sizeof(float), 1, width_fp) != 1)  { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }

          index_local_anchors_seeded_y[pos] = iv_local;

        }

        // write x anchor points computed in this step (global index)
        long index_local_anchors_comp_x[n_anchors_comp_x];
        for (long pos=0; pos<n_anchors_comp_x; pos++) {

          long iv_local = Block.max_anchors_comp[segment_ID][OR_X][pos];

          long iz_local, iy_local, ix_local;
          IndexToIndices(iv_local, ix_local, iy_local, iz_local, input_sheet_size, input_row_size);


          long iz_padded = iz_local + 1;
          long iy_padded = iy_local + 1;
          long ix_padded = ix_local + 1;

          long iv_padded = IndicesToIndex(ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);
          float width = widths[iv_padded];

          long iz_global = iz_local + block_ind[OR_Z]*input_blocksize[OR_Z];
          long iy_global = iy_local + block_ind[OR_Y]*input_blocksize[OR_Y];
          long ix_global = ix_local + block_ind[OR_X]*input_blocksize[OR_X];

          long iv_global = IndicesToIndex(ix_global, iy_global, iz_global, volume_sheet_size, volume_row_size);

          if (fwrite(&iv_global, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
          if (fwrite(&iv_global, sizeof(long), 1, width_fp) != 1) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
          if (fwrite(&width, sizeof(float), 1, width_fp) != 1)  { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }

          index_local_anchors_comp_x[pos] = iv_local;

        }

        // write x anchor points seeded in this step (global index)
        long index_local_anchors_seeded_x[n_anchors_seeded_x];
        for (long pos=0; pos<n_anchors_seeded_x; pos++) {

          long iv_local = Block.min_anchors_seeded[segment_ID][OR_X][pos];

          long iz_local, iy_local, ix_local;
          IndexToIndices(iv_local, ix_local, iy_local, iz_local, input_sheet_size, input_row_size);


          long iz_padded = iz_local + 1;
          long iy_padded = iy_local + 1;
          long ix_padded = ix_local + 1;

          long iv_padded = IndicesToIndex(ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);
          float width = widths[iv_padded];

          long iz_global = iz_local + block_ind[OR_Z]*input_blocksize[OR_Z];
          long iy_global = iy_local + block_ind[OR_Y]*input_blocksize[OR_Y];
          long ix_global = ix_local + block_ind[OR_X]*input_blocksize[OR_X];

          long iv_global = IndicesToIndex(ix_global, iy_global, iz_global, volume_sheet_size, volume_row_size);

          if (fwrite(&iv_global, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
          if (fwrite(&iv_global, sizeof(long), 1, width_fp) != 1) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
          if (fwrite(&width, sizeof(float), 1, width_fp) != 1)  { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }

          index_local_anchors_seeded_x[pos] = iv_local;

        }

        // write the synapses (global index)
        long index_local_synapses[n_synapses];
        for (long pos=0; pos<n_synapses; pos++) {

          long iv_local = Block.synapses[segment_ID][pos];

          long iz_local, iy_local, ix_local;
          IndexToIndices(iv_local, ix_local, iy_local, iz_local, input_sheet_size, input_row_size);

          long iz_global = iz_local + block_ind[OR_Z]*input_blocksize[OR_Z];
          long iy_global = iy_local + block_ind[OR_Y]*input_blocksize[OR_Y];
          long ix_global = ix_local + block_ind[OR_X]*input_blocksize[OR_X];

          long iv_global = IndicesToIndex(ix_global, iy_global, iz_global, volume_sheet_size, volume_row_size);

          if (fwrite(&iv_global, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }

          index_local_synapses[pos] = iv_local;

        }

        // write surface voxels (local index)
        for (int j=0; j<num; j++){
          if (fwrite(&index_local[j], sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
          skeleton[index_local[j]]=1;
        }

        // write z anchor points computed in this step (local index)
        for (int j=0; j<n_anchors_comp_z; j++){
          if (fwrite(&index_local_anchors_comp_z[j], sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
          skeleton[index_local_anchors_comp_z[j]]=1;
        }

        // write z anchor points seeded in this step (local index)
        for (int j=0; j<n_anchors_seeded_z; j++){
          if (fwrite(&index_local_anchors_seeded_z[j], sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
          skeleton[index_local_anchors_seeded_z[j]]=1;
        }

        // write y anchor points computed in this step (local index)
        for (int j=0; j<n_anchors_comp_y; j++){
          if (fwrite(&index_local_anchors_comp_y[j], sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
          skeleton[index_local_anchors_comp_y[j]]=1;
        }

        // write y anchor points seeded in this step (local index)
        for (int j=0; j<n_anchors_seeded_y; j++){
          if (fwrite(&index_local_anchors_seeded_y[j], sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
          skeleton[index_local_anchors_seeded_y[j]]=1;
        }

        // write x anchor points computed in this step (local index)
        for (int j=0; j<n_anchors_comp_x; j++){
          if (fwrite(&index_local_anchors_comp_x[j], sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
          skeleton[index_local_anchors_comp_x[j]]=1;
        }

        // write x anchor points seeded in this step (local index)
        for (int j=0; j<n_anchors_seeded_x; j++){
          if (fwrite(&index_local_anchors_seeded_x[j], sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
          skeleton[index_local_anchors_seeded_x[j]]=1;
        }

        // write the synapses (local index)
        for (int j=0; j<n_synapses; j++){
          if (fwrite(&index_local_synapses[j], sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
          skeleton[index_local_synapses[j]]=1;
        }

        // write checkvalue at end
        long checkvalue = 2147483647;
        if (fwrite(&checkvalue, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }

        // add randomly the first surface voxel as the soma
        skeleton[index_local[0]]=4;

        // close the I/O files
        fclose(wfp);
        fclose(width_fp);


        // double total_time = (double) (clock() - start_time) / CLOCKS_PER_SEC;
        //
        // char time_filename[4096];
        // sprintf(time_filename, "running_times/skeletons/%s-%06ld.time", prefix, segment_ID);
        //
        // FILE *tfp = fopen(time_filename, "wb");
        // if (!tfp) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }
        //
        // // write the number of points and the total time to file
        // if (fwrite(&initial_points, sizeof(long), 1, tfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }
        // if (fwrite(&total_time, sizeof(double), 1, tfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }


        // // close file
        // fclose(tfp);


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

      std::cout << "Projecting n synapses: " << Block.synapses_off[segment_ID].size() << std::endl;

      for (std::vector<long>::iterator it = Block.synapses_off[segment_ID].begin(); it != Block.synapses_off[segment_ID].end(); ++it){

        long linear_index = *it;

        long iz_unpadded, iy_unpadded, ix_unpadded;
        IndexToIndices(linear_index, ix_unpadded, iy_unpadded, iz_unpadded, input_sheet_size, input_row_size);
        // std::cout << "Original: " << iz_unpadded <<","<< iy_unpadded<<","<<ix_unpadded << std::endl;
        long iz_padded = iz_unpadded + 1;
        long iy_padded = iy_unpadded + 1;
        long ix_padded = ix_unpadded + 1;

        long cubesize = 3;
        bool projection_found = 0;
        long closest_index_padded = -1;

        float* i1 = std::min_element(resolution, resolution+2);
        float resolution_min = *i1;

        // std::cout <<"min resolution is: " << resolution_min << std::endl;

        while (projection_found==0){

          cubesize = (long)(cubesize*1.5);
          double min_distance = (double)(cubesize*resolution_min);
          if (min_distance>(resolution_min*30)){ fprintf(stderr, "Fail in finding projection point, search distance too large.\n"); exit(-1); }

          for (long iw = iz_padded - cubesize; iw <= iz_padded + cubesize; ++iw) {
              for (long iv = iy_padded - cubesize; iv <= iy_padded + cubesize; ++iv) {
                  for (long iu = ix_padded - cubesize; iu <= ix_padded + cubesize; ++iu) {


                    long index_padded = IndicesToIndex(iu, iv, iw, padded_sheet_size, padded_row_size);
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

        if (closest_index_padded == -1){ fprintf(stderr, "Fail in finding projection point.\n"); exit(-1); }
        // std::cout << "Projection found with cubesize: " << cubesize << std::endl;
        Block.Pointclouds[segment_ID][closest_index_padded] = 3;
        segment[closest_index_padded] = 3;

        IndexToIndices(closest_index_padded, ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);

        iz_unpadded = iz_padded - 1;
        iy_unpadded = iy_padded - 1;
        ix_unpadded = ix_padded - 1;

        // std::cout << "Projected to: " << iz_unpadded<<","<<iy_unpadded<<","<<ix_unpadded << std::endl;

        long projected_index = IndicesToIndex(ix_unpadded, iy_unpadded, iz_unpadded, input_sheet_size, input_row_size);
        Block.synapses[segment_ID].push_back(projected_index);

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

void CPPcreateDataBlock(const char *prefix, const char *lookup_table_directory, long *inp_labels, float input_resolution[3], long inp_blocksize[3], long volume_size[3], long block_ind[3], const char* synapses_dir, const char* somae_dir, const char* skeleton_dir){
  // create new Datablock and set the input variables

  clock_t start_time_total = clock();

  DataBlock BlockA(input_resolution, inp_blocksize, volume_size, block_ind, synapses_dir, somae_dir, skeleton_dir);

  BlockA.CppPopulatePointCloudFromH5(inp_labels);

  // read Synapses
  if (!BlockA.ReadSynapses(prefix)) exit(-1);

  // read Anchor points (if block before does exist)
  if (!BlockA.ReadAnchorpoints(prefix)) exit(-1);

  // initialize lookup tables
  InitializeLookupTables(lookup_table_directory);

  // needs to happen after PopulatePointCloud()
  PopulateOffsets(BlockA.padded_blocksize);

  // BlockA.IDs_to_process.insert({149});
  // insert IDs that should be processed
  BlockA.IDs_to_process = BlockA.IDs_in_block;
  BlockA.IDs_to_process.erase(55);
  BlockA.IDs_to_process.erase(81);
  BlockA.IDs_to_process.erase(301);

  BlockA.writeIDsToProcess(prefix);

  std::unordered_set<long>::iterator itr = BlockA.IDs_to_process.begin();

  while (itr != BlockA.IDs_to_process.end())
  {

      // initialize segment using the Block object and the segment ID to process
      BlockSegment* segA = new BlockSegment(*itr, BlockA);

      // call the sequential thinning algorithm
      segA->SequentialThinning(prefix, BlockA);

      // write skeletons and widths, add anchor points to
      segA->WriteOutputfiles(prefix, BlockA);

      delete segA;

      itr++;

  }

  // write border anchor points to file
  BlockA.writeAnchorsComputed(prefix);
  BlockA.writeanchorsSeeded(prefix);
  BlockA.WriteProjectedSynapses(prefix);
  double time_total = (double) (clock() - start_time_total) / CLOCKS_PER_SEC;
  std::cout << "time added summed: " << time_total << std::endl;

}

// clock_t start_time = clock();
// time_added += (double) (clock() - start_time) / CLOCKS_PER_SEC;
// std::cout << "time added summed: " << time_added << std::endl;
