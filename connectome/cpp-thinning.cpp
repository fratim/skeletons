/* c++ file to extract wiring diagram */
#include <limits>
#include "cpp-wiring.h"
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <vector>
#include <time.h>
#include <stdio.h>

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
    long input_blocksize[3];
    long volumesize[3];
    float resolution[3];
    long block_ind[3];
    long padded_row_size;
    long padded_sheet_size;
    long input_row_size;
    long input_sheet_size;
    const char *synapses_directory;
    const char *somae_directory;
    const char *skeleton_directory;

  public:
    long padded_blocksize[3];
    std::unordered_set<long> IDs_to_process;
    std::unordered_set<long> IDs_in_block;
    std::unordered_map<long, std::unordered_map<long, char>> Pointclouds;
    std::unordered_map<long, std::unordered_map<long, std::unordered_set<long>>> borderpoints;
    std::vector<long> zmax_iy_local = std::vector<long>();
    std::vector<long> zmax_ix_local = std::vector<long>();
    std::vector<long> zmax_segment_ID = std::vector<long>();

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

        for (long voxel_index = 0; voxel_index < n_points; voxel_index++){

          // get segment_ID of current index and skip if is zero
          long curr_label = inp_labels[voxel_index];
          if (!curr_label) continue;

          long ix, iy, iz;
          IndexToIndices(voxel_index, ix, iy, iz, input_sheet_size, input_row_size);

          long iz_padded; long iy_padded; long ix_padded;

          //  pad the location by one
          iz_padded = iz + 1;
          iy_padded = iy + 1;
          ix_padded = ix + 1;

          // find the new voxel index
          long iv = IndicesToIndex(ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);

          // check if pointcloud of this segment_ID already exists, otherwise add new pointcloud
          if (Pointclouds.find(curr_label) == Pointclouds.end()) {
            Pointclouds[curr_label] = std::unordered_map<long,char>();
            // std::cout << "New segment_ID detected: " << curr_label << std::endl << std::flush;
            IDs_in_block.insert(curr_label);
          }

          Pointclouds[curr_label][iv] = 1;

          // add index to borderpoints unordered_map (only for max walls, as these anchor points are then copied to next level)
          if (iz==(input_blocksize[OR_Z]-1)) borderpoints[curr_label][OR_Z].insert(iv);
          else if (iy==(input_blocksize[OR_Y]-1)) borderpoints[curr_label][OR_Y].insert(iv);
          else if (ix==(input_blocksize[OR_X]-1)) borderpoints[curr_label][OR_X].insert(iv);

        }

    }

    int ReadSynapses(const char *prefix)
    {

      // read the synapses
      char synapse_filename[4096];
      snprintf(synapse_filename, 4096, "%s/%s/%s-synapses-%04ldz-%04ldy-%04ldx.pts", synapses_directory, prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

      FILE *fp = fopen(synapse_filename, "rb");
      if (!fp) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

      long z_input_volume_size, y_input_volume_size, x_input_volume_size;
      long z_input_block_size, y_input_block_size, x_input_block_size;

      if (fread(&z_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
      if (fread(&y_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
      if (fread(&x_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

      if (z_input_volume_size != volumesize[OR_Z]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
      if (y_input_volume_size != volumesize[OR_Y]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
      if (x_input_volume_size != volumesize[OR_X]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

      if (fread(&z_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
      if (fread(&y_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
      if (fread(&x_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

      if (z_input_block_size != input_blocksize[OR_Z]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
      if (y_input_block_size != input_blocksize[OR_Y]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
      if (x_input_block_size != input_blocksize[OR_X]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

      long nneurons;
      if (fread(&nneurons, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

      for (long iv = 0; iv < nneurons; ++iv) {
          // get the label and number of synapses
          long segment_ID;
          long nsynapses;

          if (fread(&segment_ID, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
          if (fread(&nsynapses, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

          // synapses[segment_ID] = std::vector<long>();

          // ignore the global coordinates
          for (long is = 0; is < nsynapses; ++is) {
              long dummy_index;
              if (fread(&dummy_index, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
          }


          // add the local coordinates (offset by one in each direction)
          for (long is = 0; is < nsynapses; ++is) {
              long linear_index;
              if (fread(&linear_index, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

              long ix, iy, iz;
              IndexToIndices(linear_index, ix, iy, iz, input_sheet_size, input_row_size);

              //  pad the location by one
              iz += 1; iy += 1; ix += 1;

              // find the new voxel index
              long iv = IndicesToIndex(ix, iy, iz, padded_sheet_size, padded_row_size);

              Pointclouds[segment_ID][iv] = 3;
          }

      }

      // close file
      fclose(fp);

      return 1;

    }

    int ReadAnchorpoints(const char *prefix)
    {
      // make filename of adjacent z lock (in negative direction)
      char output_filename_zmax[4096];
      sprintf(output_filename_zmax, "%s/%s/%s-borderZMax-%04ldz-%04ldy-%04ldx.pts", skeleton_directory, prefix, prefix, block_ind[OR_Z]-1, block_ind[OR_Y], block_ind[OR_X]);

      std::cout << "Reading Anchor points from: " << output_filename_zmax << std::endl;

      FILE *fpzmax = fopen(output_filename_zmax, "rb");
      if (!fpzmax) { fprintf(stderr, "Failed to read %s.\n", output_filename_zmax); return 0; }

      long npoints;
      if (fread(&npoints, sizeof(long), 1, fpzmax) != 1) { fprintf(stderr, "Failed to read %s.\n", output_filename_zmax); return 0; }

      std::cout << "number of anchor points: " << npoints << std::endl;

      for (long pos = 0; pos < npoints; pos++) {

          // get the label and number of synapses
          long iy_local;
          long ix_local;
          long segment_ID;

          if (fread(&iy_local, sizeof(long), 1, fpzmax) != 1) { fprintf(stderr, "Failed to read %s.\n", output_filename_zmax); return 0; }
          if (fread(&ix_local, sizeof(long), 1, fpzmax) != 1) { fprintf(stderr, "Failed to read %s.\n", output_filename_zmax); return 0; }
          if (fread(&segment_ID, sizeof(long), 1, fpzmax) != 1) { fprintf(stderr, "Failed to read %s.\n", output_filename_zmax); return 0; }

          // populate offset
          long iz_padded = 1;
          long iy_padded = iy_local + 1;
          long ix_padded = ix_local + 1;


          //  find the new voxel index
          long iv = IndicesToIndex(ix_padded, iy_padded, iz_padded, padded_sheet_size, padded_row_size);

          Pointclouds[segment_ID][iv] = 3;

          // std::cout << "added anchor point at: " << iy_local << "," << ix_local << "segment ID: " << segment_ID << std::endl;

      }

      // close file
      fclose(fpzmax);

      return 1;

    }

    void writeZmaxBlock (const char *prefix)
    {

      // write the zmax border points to a file (so far stored in vectors)
      // create an output file that saves th epoints on the positive z bounday, as (global index, local index, segment_ID)*n_points , npoints
      char output_filename_zmax[4096];
      sprintf(output_filename_zmax, "%s/%s/%s-borderZMax-%04ldz-%04ldy-%04ldx.pts", skeleton_directory, prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

      FILE *zmaxfp = fopen(output_filename_zmax, "wb");
      if (!zmaxfp) { fprintf(stderr, "Failed to open %s\n", output_filename_zmax); exit(-1); }

      if (!((zmax_iy_local.size()==zmax_ix_local.size())&&(zmax_ix_local.size()==zmax_segment_ID.size()))){
        throw std::invalid_argument("z vectors not of same size!");
      }

      long n_points_zmax = zmax_iy_local.size();

      std::cout << "Anchor points found: " << n_points_zmax << std::endl;

      if (fwrite(&n_points_zmax, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }

      for (int pos = 0; pos < n_points_zmax; pos++){
        if (fwrite(&zmax_iy_local[pos], sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
        if (fwrite(&zmax_ix_local[pos], sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
        if (fwrite(&zmax_segment_ID[pos], sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
      }

      fclose(zmaxfp);
    }

};

class BlockSegment : public DataBlock{
  private:
    long initial_points;
    long segment_ID;
    List surface_voxels;
    std::unordered_map<long, char> segment;
    std::unordered_map<long, std::unordered_set<long>> borderpoints_segment = std::unordered_map<long, std::unordered_set<long>>();
    std::unordered_map<long, float> widths = std::unordered_map<long, float>();

  public:
    BlockSegment(long segment_ID_inp, DataBlock &Blockx):DataBlock(Blockx){

      segment_ID = segment_ID_inp;
      segment = Blockx.Pointclouds[segment_ID];
      borderpoints_segment = Blockx.borderpoints[segment_ID];

      std::cout << "-----------------------------------" << std::endl;
      std::cout << "Processing segment_ID " << segment_ID << std::endl;

      initial_points = segment.size();
      printf("segment_ID %ld initial points: %ld\n", segment_ID, initial_points);

    }

    void SequentialThinning(const char *prefix)
    {
        // create a vector of surface voxels
        CollectSurfaceVoxels();
        int iteration = 0;
        long changed = 0;
        do {
            changed = ThinningIterationStep();
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

    long ThinningIterationStep(void)
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

                bool isAnchorPoint = 0;

                if (borderpoints_segment[OR_Y].count(index)){

                    long sum_of_neighbors = 0; //collect voxels of neighbors on x-z plane
                    sum_of_neighbors += segment[index+ n26_offsets[3]];
                    sum_of_neighbors += segment[index+ n26_offsets[4]];
                    sum_of_neighbors += segment[index+ n26_offsets[5]];
                    sum_of_neighbors += segment[index+ n26_offsets[12]];
                    sum_of_neighbors += segment[index+ n26_offsets[13]];
                    sum_of_neighbors += segment[index+ n26_offsets[20]];
                    sum_of_neighbors += segment[index+ n26_offsets[21]];
                    sum_of_neighbors += segment[index+ n26_offsets[22]];

                    if (sum_of_neighbors==0) isAnchorPoint = 1;
                }
                else if (borderpoints_segment[OR_Z].count(index)){

                    long sum_of_neighbors = 0; //collect voxels of neighbors on x-y plane
                    sum_of_neighbors += segment[index+ n26_offsets[9]];
                    sum_of_neighbors += segment[index+ n26_offsets[10]];
                    sum_of_neighbors += segment[index+ n26_offsets[11]];
                    sum_of_neighbors += segment[index+ n26_offsets[12]];
                    sum_of_neighbors += segment[index+ n26_offsets[13]];
                    sum_of_neighbors += segment[index+ n26_offsets[14]];
                    sum_of_neighbors += segment[index+ n26_offsets[15]];
                    sum_of_neighbors += segment[index+ n26_offsets[16]];

                    if (sum_of_neighbors==0) isAnchorPoint = 1;

                }
                else if (borderpoints_segment[OR_X].count(index)){

                    long sum_of_neighbors = 0; //collect voxels of neighbors on y-z plane
                    sum_of_neighbors += segment[index+ n26_offsets[7]];
                    sum_of_neighbors += segment[index+ n26_offsets[4]];
                    sum_of_neighbors += segment[index+ n26_offsets[1]];
                    sum_of_neighbors += segment[index+ n26_offsets[10]];
                    sum_of_neighbors += segment[index+ n26_offsets[15]];
                    sum_of_neighbors += segment[index+ n26_offsets[24]];
                    sum_of_neighbors += segment[index+ n26_offsets[21]];
                    sum_of_neighbors += segment[index+ n26_offsets[18]];

                    if (sum_of_neighbors==0) isAnchorPoint = 1;
                }

                // if anchor point detected, fix it and do not check if simple
                if (isAnchorPoint){
                  // fix anchor point (as a fake synapse) TODO: might need another label here to identify anchor points seperately
                  segment[index]=3;
                }

                // otherise do normal procedure - hceck if simple, if so, delete it
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

    void WriteOutputfiles(const char *prefix, DataBlock &Blockx)
    {

        // count the number of remaining points
        long num = 0;
        ListElement *LE = (ListElement *) surface_voxels.first;
        while (LE != NULL) {
            num++;
            LE = (ListElement *)LE->next;
        }

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

        // write the characteristics header
        WriteHeader(wfp, num);

        // characteristics
        WriteHeader(width_fp, num);

        printf("Remaining voxels: %ld\n", num);

        long index_local[num];
        long counter_it = 0;

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

            long iv_local = iz_local * input_blocksize[OR_X] * input_blocksize[OR_Y] + iy_local * input_blocksize[OR_X] + ix_local;
            long iv_global = iz_global * volumesize[OR_X] * volumesize[OR_Y] + iy_global * volumesize[OR_X] + ix_global;

            // write to zmax border file, if on border
            if (iz_local==(input_blocksize[OR_Z]-1)){
                Blockx.zmax_iy_local.push_back(iy_local);
                Blockx.zmax_ix_local.push_back(ix_local);
                Blockx.zmax_segment_ID.push_back(segment_ID);
            }

            if (fwrite(&iv_global, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
            if (fwrite(&iv_global, sizeof(long), 1, width_fp) != 1) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
            if (fwrite(&width, sizeof(float), 1, width_fp) != 1)  { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }

            index_local[counter_it]=iv_local;

            // remove this voxel
            RemoveSurfaceVoxel(LE, surface_voxels);
            counter_it++;
        }

        for (int j=0; j<num; j++){
            if (fwrite(&index_local[j], sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
        }

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

    void WriteHeader(FILE *fp, long &num){
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

  // insert IDs that should be processed
  // BlockA.IDs_to_process.insert({91});

  BlockA.CppPopulatePointCloudFromH5(inp_labels);

  // read Synapses
  if (!BlockA.ReadSynapses(prefix)) exit(-1);

  // read Anchor points (if block before does exist)
  // if (!BlockA.ReadAnchorpoints(prefix)) exit(-1);

  // initialize lookup tables
  InitializeLookupTables(lookup_table_directory);

  // needs to happen after PopulatePointCloud()
  PopulateOffsets(BlockA.padded_blocksize);

  // create iterator over set
  std::unordered_set<long>::iterator itr = BlockA.IDs_in_block.begin();

  while (itr != BlockA.IDs_in_block.end())
  {
      // initialize segment using the Block object and the segment ID to process
      BlockSegment segA(*itr, BlockA);

      // call the sequential thinning algorithm
      segA.SequentialThinning(prefix);

      // write skeletons and widths, add anchor points to
      segA.WriteOutputfiles(prefix, BlockA);

      // increment iterator
      itr++;
  }

  // write border anchor points to file
  BlockA.writeZmaxBlock(prefix);
}

// clock_t start_time = clock();
// time_added += (double) (clock() - start_time) / CLOCKS_PER_SEC;
// std::cout << "time added summed: " << time_added << std::endl;
