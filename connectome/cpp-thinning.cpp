/* c++ file to extract wiring diagram */
#include <limits>
#include "cpp-wiring.h"
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <vector>

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

class DataBlock{
  public:
    float resolution_test[3];
    long input_blocksize[3];
    long padded_blocksize[3];
    long volumesize[3];
    float resolution[3];
    long block_z;
    long block_y;
    long block_x;
    const char *synapses_directory;
    const char *somae_directory;
    const char *skeleton_directory;
    std::unordered_map<long, std::unordered_map<long, short>> Pointclouds;
    std::unordered_map<long, std::unordered_map<long, std::unordered_set<long>>> borderpoints;

    std::unordered_set<long> IDs_to_process;
    std::unordered_set<long> IDs_in_block;
    std::vector<long> zmax_iy_local = std::vector<long>();
    std::vector<long> zmax_ix_local = std::vector<long>();
    std::vector<long> zmax_segment_ID = std::vector<long>();

    void CppUpdateResolution(float input_resolution[3])
    {
        resolution_test[OR_Z] = input_resolution[OR_Z];
        resolution_test[OR_Y] = input_resolution[OR_Y];
        resolution_test[OR_X] = input_resolution[OR_X];

        std::cout << "Resolution set to: " << resolution_test[OR_Z] << "," << resolution_test[OR_Y] << "," << resolution_test[OR_X] << "," << std::endl;

    }

    void CppUpdateBlocksize(long inp_blocksize[3])
    {
        input_blocksize[OR_Z] = inp_blocksize[OR_Z];
        input_blocksize[OR_Y] = inp_blocksize[OR_Y];
        input_blocksize[OR_X] = inp_blocksize[OR_X];

        padded_blocksize[OR_Z] = inp_blocksize[OR_Z]+2;
        padded_blocksize[OR_Y] = inp_blocksize[OR_Y]+2;
        padded_blocksize[OR_X] = inp_blocksize[OR_X]+2;

        std::cout << "Blocksize input set to: " << input_blocksize[OR_Z] << "," << input_blocksize[OR_Y] << "," << input_blocksize[OR_X] << "," << std::endl;
        std::cout << "Blocksize padded set to: " << padded_blocksize[OR_Z] << "," << padded_blocksize[OR_Y] << "," << padded_blocksize[OR_X] << "," << std::endl;

    }

    void CppUpdateVolumesize(long volume_size[3])
    {
        volumesize[OR_Z] = volume_size[OR_Z];
        volumesize[OR_Y] = volume_size[OR_Y];
        volumesize[OR_X] = volume_size[OR_X];

        std::cout << "Volumesize set to: " << volumesize[OR_Z] << "," << volumesize[OR_Y] << "," << volumesize[OR_X] << "," << std::endl;

    }

    void CppUpdateBlockindices(long block_ind_inp[3])
    {
        block_z = block_ind_inp[OR_Z];
        block_y = block_ind_inp[OR_Y];
        block_x = block_ind_inp[OR_X];

        std::cout << "Block indices set to: " << block_z << "," << block_y << "," << block_x << "," << std::endl;
    }

    void CppUpdateDirectories(const char* synapses_dir, const char* somae_dir, const char* skeleton_dir)
    {
        synapses_directory = synapses_dir;
        somae_directory = somae_dir;
        skeleton_directory = skeleton_dir;

        std::cout << "Directories set. " << std::endl;
    }

    void CppPopulatePointCloudFromH5(long *inp_labels)
    {

        // indexing parameters for indexing within current block
        long nentries = padded_blocksize[OR_Z] * padded_blocksize[OR_Y] * padded_blocksize[OR_X];
        long sheet_size = padded_blocksize[OR_Y] * padded_blocksize[OR_X];
        long row_size = padded_blocksize[OR_X];
        long infinity = padded_blocksize[OR_Z] * padded_blocksize[OR_Z] + padded_blocksize[OR_Y] * padded_blocksize[OR_Y] + padded_blocksize[OR_X] * padded_blocksize[OR_X];

        long n_points = input_blocksize[0]*input_blocksize[1]*input_blocksize[2];

        std::cout << "n_points is: " << n_points << std::endl << std::flush;

        for (long voxel_index = 0; voxel_index < n_points; voxel_index++){

          // get segment_ID of current index and skip if is zero
          long curr_label = inp_labels[voxel_index];
          if (!curr_label) continue;

          // find coordinates of this voxel_index
          long iz = voxel_index / (input_blocksize[OR_Y] * input_blocksize[OR_X]);
          long iy = (voxel_index - iz * (input_blocksize[OR_Y] * input_blocksize[OR_X])) / input_blocksize[OR_X];
          long ix = voxel_index % input_blocksize[OR_X];

          long iz_padded; long iy_padded; long ix_padded;

          //  pad the location by one
          iz_padded = iz + 1;
          iy_padded = iy + 1;
          ix_padded = ix + 1;

          // find the new voxel index
          long iv = IndicesToIndex(ix_padded, iy_padded, iz_padded, sheet_size, row_size);

          // check if pointcloud of this segment_ID already exists, otherwise add new pointcloud
          if (Pointclouds.find(curr_label) == Pointclouds.end()) {
            Pointclouds[curr_label] = std::unordered_map<long,short>();
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
      snprintf(synapse_filename, 4096, "%s/%s/%s-synapses-%04ldz-%04ldy-%04ldx.pts", synapses_directory, prefix, prefix, block_z, block_y, block_x);

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

          // set block indexing parameters
          long nentries = padded_blocksize[OR_Z] * padded_blocksize[OR_Y] * padded_blocksize[OR_X];
          long sheet_size = padded_blocksize[OR_Y] * padded_blocksize[OR_X];
          long row_size = padded_blocksize[OR_X];
          long infinity = padded_blocksize[OR_Z] * padded_blocksize[OR_Z] + padded_blocksize[OR_Y] * padded_blocksize[OR_Y] + padded_blocksize[OR_X] * padded_blocksize[OR_X];

          // add the local coordinates (offset by one in each direction)
          for (long is = 0; is < nsynapses; ++is) {
              long linear_index;
              if (fread(&linear_index, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }

              long iz = linear_index / (input_blocksize[OR_Y] * input_blocksize[OR_X]);
              long iy = (linear_index - iz * (input_blocksize[OR_Y] * input_blocksize[OR_X])) / input_blocksize[OR_X];
              long ix = linear_index % input_blocksize[OR_X];

              //  pad the location by one
              iz += 1; iy += 1; ix += 1;

              // find the new voxel index
              long iv = IndicesToIndex(ix, iy, iz, sheet_size, row_size);

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
      sprintf(output_filename_zmax, "%s/%s/%s-borderZMax-%04ldz-%04ldy-%04ldx.pts", skeleton_directory, prefix, prefix, block_z-1, block_y, block_x);

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

          long nentries = padded_blocksize[OR_Z] * padded_blocksize[OR_Y] * padded_blocksize[OR_X];
          long sheet_size = padded_blocksize[OR_Y] * padded_blocksize[OR_X];
          long row_size = padded_blocksize[OR_X];
          long infinity = padded_blocksize[OR_Z] * padded_blocksize[OR_Z] + padded_blocksize[OR_Y] * padded_blocksize[OR_Y] + padded_blocksize[OR_X] * padded_blocksize[OR_X];

          //  find the new voxel index
          long iv = IndicesToIndex(ix_padded, iy_padded, iz_padded, sheet_size, row_size);

          Pointclouds[segment_ID][iv] = 3;

          // std::cout << "added anchor point at: " << iy_local << "," << ix_local << "segment ID: " << segment_ID << std::endl;

      }

      // close file
      fclose(fpzmax);

      return 1;

    }

};

class BlockSegment : public DataBlock{
  private:
    std::unordered_map<long, short> segment;
    long segment_ID;
    std::unordered_map<long, std::unordered_set<long>> borderpoints_segment;
    std::unordered_map<long, float> widths;
    List surface_voxels;

  public:
    BlockSegment(long segment_ID_inp, DataBlock Blockx){
      std::copy(std::begin(Blockx.resolution_test), std::end(Blockx.resolution_test), std::begin(resolution_test));
      std::copy(std::begin(Blockx.input_blocksize), std::end(Blockx.input_blocksize), std::begin(input_blocksize));
      std::copy(std::begin(Blockx.padded_blocksize), std::end(Blockx.padded_blocksize), std::begin(padded_blocksize));
      std::copy(std::begin(Blockx.volumesize), std::end(Blockx.volumesize), std::begin(volumesize));
      std::copy(std::begin(Blockx.resolution), std::end(Blockx.resolution), std::begin(resolution));
      block_z = Blockx.block_z;
      block_y = Blockx.block_y;
      block_x = Blockx.block_x;
      synapses_directory = Blockx.synapses_directory;
      somae_directory = Blockx.somae_directory;
      skeleton_directory = Blockx.skeleton_directory;

      segment_ID = segment_ID_inp;
      segment = Blockx.Pointclouds[segment_ID];
      borderpoints_segment = Blockx.borderpoints[segment_ID];

      std::cout << "-----------------------------------" << std::endl;
      std::cout << "Processing segment_ID " << segment_ID << std::endl;

      long initial_points = segment.size();
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

};

// static void NewSurfaceVoxel(long iv, long ix, long iy, long iz)
// {
//     ListElement *LE = new ListElement();
//     LE->iv = iv;
//     LE->ix = ix;
//     LE->iy = iy;
//     LE->iz = iz;
//
//     LE->next = NULL;
//     LE->prev = surface_voxels.last;
//
//     if (surface_voxels.last != NULL) ((ListElement *) surface_voxels.last)->next = LE;
//     surface_voxels.last = LE;
//     if (surface_voxels.first == NULL) surface_voxels.first = LE;
// }
//
// static void RemoveSurfaceVoxel(ListElement *LE)
// {
//     ListElement *LE2;
//     if (surface_voxels.first == LE) surface_voxels.first = LE->next;
//     if (surface_voxels.last == LE) surface_voxels.last = LE->prev;
//
//     if (LE->next != NULL) {
//         LE2 = (ListElement *)(LE->next);
//         LE2->prev = LE->prev;
//     }
//     if (LE->prev != NULL) {
//         LE2 = (ListElement *)(LE->prev);
//         LE2->next = LE->next;
//     }
//     delete LE;
// }
//
// static void CreatePointList(PointList *s)
// {
//     s->head = NULL;
//     s->tail = NULL;
//     s->length = 0;
// }
//
// static void AddToList(PointList *s, Voxel e, ListElement *ptr)
// {
//     Cell *newcell = new Cell();
//     newcell->v = e;
//     newcell->ptr = ptr;
//     newcell->next = NULL;
//
//     if (s->head == NULL) {
//         s->head = newcell;
//         s->tail = newcell;
//         s->length = 1;
//     }
//     else {
//         s->tail->next = newcell;
//         s->tail = newcell;
//         s->length++;
//     }
// }
//
// static Voxel GetFromList(PointList *s, ListElement **ptr)
// {
//     Voxel V;
//     Cell *tmp;
//     V.iv = -1;
//     V.ix = -1;
//     V.iy = -1;
//     V.iz = -1;
//     (*ptr) = NULL;
//     if (s->length == 0) return V;
//     else {
//         V = s->head->v;
//         (*ptr) = s->head->ptr;
//         tmp = (Cell *) s->head->next;
//         delete s->head;
//         s->head = tmp;
//         s->length--;
//         if (s->length == 0) {
//             s->head = NULL;
//             s->tail = NULL;
//         }
//         return V;
//     }
// }
//
// static void DestroyPointList(PointList *s)
// {
//     ListElement *ptr;
//     while (s->length) GetFromList(s, &ptr);
// }

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

// static void CollectSurfaceVoxels(void)
// {
//     // go through all voxels and check their six neighbors
//     for (std::unordered_map<long, short>::iterator it = segment.begin(); it != segment.end(); ++it) {
//         // all of these elements are either 1 or 3 and in the segment
//         long index = it->first;
//
//         // initialize widths to maximum float value
//         widths[index] = std::numeric_limits<float>::max();
//
//         long ix, iy, iz;
//         IndexToIndices(index, ix, iy, iz);
//         // check the 6 neighbors
//         for (long iv = 0; iv < NTHINNING_DIRECTIONS; ++iv) {
//             long neighbor_index = index + n6_offsets[iv];
//
//             long ii, ij, ik;
//             IndexToIndices(neighbor_index, ii, ij, ik);
//
//             // skip the fake boundary elements
//             if ((ii == 0) or (ii == padded_blocksize[OR_X] - 1)) continue;
//             if ((ij == 0) or (ij == padded_blocksize[OR_Y] - 1)) continue;
//             if ((ik == 0) or (ik == padded_blocksize[OR_Z] - 1)) continue;
//
//             if (segment.find(neighbor_index) == segment.end()) {
//                 // this location is a boundary so create a surface voxel and break
//                 // cannot update it->second if it is synapse so need this test!!
//                 if (it->second == 1) {
//                     it->second = 2;
//                     NewSurfaceVoxel(index, ix, iy, iz);
//                 }
//
//                 // note this location as surface
//                 widths[index] = 0;
//
//                 break;
//             }
//         }
//     }
// }
//
// static unsigned int Collect26Neighbors(long ix, long iy, long iz)
// {
//     unsigned int neighbors = 0;
//     long index = IndicesToIndex(ix, iy, iz);
//
//     // some of these lookups will create a new entry but the region is
//     // shrinking so memory overhead is minimal
//     for (long iv = 0; iv < 26; ++iv) {
//         if (segment[index + n26_offsets[iv]]) neighbors |= long_mask[iv];
//     }
//
//     return neighbors;
// }
//
// static bool Simple26_6(unsigned int neighbors)
// {
//     return lut_simple[(neighbors >> 3)] & char_mask[neighbors % 8];
// }
//
// static void DetectSimpleBorderPoints(PointList *deletable_points, int direction)
// {
//     ListElement *LE = (ListElement *)surface_voxels.first;
//     while (LE != NULL) {
//         long iv = LE->iv;
//         long ix = LE->ix;
//         long iy = LE->iy;
//         long iz = LE->iz;
//
//         // not a synapse endpoint (need this here since endpoints are on the list of surfaces)
//         // this will only be called on things on the surface already so already in unordered_map
//         if (segment[iv] == 2) {
//             long value = 0;
//             // is the neighbor in the corresponding direction not in the segment
//             // some of these keys will not exist but will default to 0 value
//             // the search region retracts in from the boundary so limited memory overhead
//             // the n6_offsets are in the order UP, DOWN, NORTH, SOUTH, EAST, WEST
//             value = segment[iv + n6_offsets[direction]];
//
//             // see if the required point belongs to a different segment
//             if (!value) {
//                 unsigned int neighbors = Collect26Neighbors(ix, iy, iz);
//
//                 // deletable point
//                 if (Simple26_6(neighbors)) {
//                     Voxel voxel;
//                     voxel.iv = iv;
//                     voxel.ix = ix;
//                     voxel.iy = iy;
//                     voxel.iz = iz;
//                     AddToList(deletable_points, voxel, LE);
//                 }
//             }
//         }
//         LE = (ListElement *) LE->next;
//     }
// }
//
// static long ThinningIterationStep(void)
// {
//     long changed = 0;
//
//     // iterate through every direction
//     for (int direction = 0; direction < NTHINNING_DIRECTIONS; ++direction) {
//         PointList deletable_points;
//         ListElement *ptr;
//
//         CreatePointList(&deletable_points);
//         DetectSimpleBorderPoints(&deletable_points, direction);
//
//         while (deletable_points.length) {
//             Voxel voxel = GetFromList(&deletable_points, &ptr);
//
//             long index = voxel.iv;
//             long ix = voxel.ix;
//             long iy = voxel.iy;
//             long iz = voxel.iz;
//
//             bool isAnchorPoint = 0;
//
//             if (borderpoints_segment[OR_Y].count(index)){
//
//                 long sum_of_neighbors = 0; //collect voxels of neighbors on x-z plane
//                 sum_of_neighbors += segment[index+ n26_offsets[3]];
//                 sum_of_neighbors += segment[index+ n26_offsets[4]];
//                 sum_of_neighbors += segment[index+ n26_offsets[5]];
//                 sum_of_neighbors += segment[index+ n26_offsets[12]];
//                 sum_of_neighbors += segment[index+ n26_offsets[13]];
//                 sum_of_neighbors += segment[index+ n26_offsets[20]];
//                 sum_of_neighbors += segment[index+ n26_offsets[21]];
//                 sum_of_neighbors += segment[index+ n26_offsets[22]];
//
//                 if (sum_of_neighbors==0) isAnchorPoint = 1;
//             }
//             else if (borderpoints_segment[OR_Z].count(index)){
//
//                 long sum_of_neighbors = 0; //collect voxels of neighbors on x-y plane
//                 sum_of_neighbors += segment[index+ n26_offsets[9]];
//                 sum_of_neighbors += segment[index+ n26_offsets[10]];
//                 sum_of_neighbors += segment[index+ n26_offsets[11]];
//                 sum_of_neighbors += segment[index+ n26_offsets[12]];
//                 sum_of_neighbors += segment[index+ n26_offsets[13]];
//                 sum_of_neighbors += segment[index+ n26_offsets[14]];
//                 sum_of_neighbors += segment[index+ n26_offsets[15]];
//                 sum_of_neighbors += segment[index+ n26_offsets[16]];
//
//                 if (sum_of_neighbors==0) isAnchorPoint = 1;
//
//             }
//             else if (borderpoints_segment[OR_X].count(index)){
//
//                 long sum_of_neighbors = 0; //collect voxels of neighbors on y-z plane
//                 sum_of_neighbors += segment[index+ n26_offsets[7]];
//                 sum_of_neighbors += segment[index+ n26_offsets[4]];
//                 sum_of_neighbors += segment[index+ n26_offsets[1]];
//                 sum_of_neighbors += segment[index+ n26_offsets[10]];
//                 sum_of_neighbors += segment[index+ n26_offsets[15]];
//                 sum_of_neighbors += segment[index+ n26_offsets[24]];
//                 sum_of_neighbors += segment[index+ n26_offsets[21]];
//                 sum_of_neighbors += segment[index+ n26_offsets[18]];
//
//                 if (sum_of_neighbors==0) isAnchorPoint = 1;
//             }
//
//             // if anchor point detected, fix it and do not check if simple
//             if (isAnchorPoint){
//               // fix anchor point (as a fake synapse) TODO: might need another label here to identify anchor points seperately
//               segment[index]=3;
//             }
//
//             // otherise do normal procedure - hceck if simple, if so, delete it
//             else{
//
//               unsigned int neighbors = Collect26Neighbors(ix, iy, iz);
//               if (Simple26_6(neighbors)) {
//                   // delete the simple point
//                   segment[index] = 0;
//
//                   // add the new surface voxels
//                   for (long ip = 0; ip < NTHINNING_DIRECTIONS; ++ip) {
//                       long neighbor_index = index + n6_offsets[ip];
//
//                       // previously not on the surface but is in the object
//                       // widths of voxels start at maximum and first updated when put on surface
//                       if (segment[neighbor_index] == 1) {
//                           long iu, iv, iw;
//                           IndexToIndices(neighbor_index, iu, iv, iw);
//                           NewSurfaceVoxel(neighbor_index, iu, iv, iw);
//
//                           // convert to a surface point
//                           segment[neighbor_index] = 2;
//                       }
//                   }
//
//                   // check all 26 neighbors to see if width is better going through this voxel
//                   for (long ip = 0; ip < 26; ++ip) {
//                       long neighbor_index = index + n26_offsets[ip];
//                       if (!segment[neighbor_index]) continue;
//
//                       // get this index in (x, y, z)
//                       long iu, iv, iw;
//                       IndexToIndices(neighbor_index, iu, iv, iw);
//
//                       // get the distance from the voxel to be deleted
//                       float diffx = resolution[OR_X] * (ix - iu);
//                       float diffy = resolution[OR_Y] * (iy - iv);
//                       float diffz = resolution[OR_Z] * (iz - iw);
//
//                       float distance = sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
//                       float current_width = widths[neighbor_index];
//
//                       if (widths[index] + distance < current_width) {
//                           widths[neighbor_index] = widths[index] + distance;
//                       }
//                   }
//
//                   // remove this from the surface voxels
//                   RemoveSurfaceVoxel(ptr);
//                   changed += 1;
//                 }
//             }
//         }
//         DestroyPointList(&deletable_points);
//     }
//
//
//     // return the number of changes
//     return changed;
// }
//
// static void SequentialThinning(const char *prefix, long segment_ID)
// {
//     // create a vector of surface voxels
//     CollectSurfaceVoxels();
//     int iteration = 0;
//     long changed = 0;
//     do {
//         changed = ThinningIterationStep();
//         iteration++;
//         // printf("  Iteration %d deleted %ld points\n", iteration, changed);
//     } while (changed);
//     printf("Needed %d iterations\n", iteration);
//
// }
//
// static void WriteOutputfiles(const char *prefix, long segment_ID, clock_t start_time, long initial_points, std::vector<long> &zmax_iy_local, std::vector<long> &zmax_ix_local, std::vector<long> &zmax_segment_ID)
// {
//
//     // count the number of remaining points
//     long num = 0;
//     ListElement *LE = (ListElement *) surface_voxels.first;
//     while (LE != NULL) {
//         num++;
//         LE = (ListElement *)LE->next;
//     }
//
//     // create an output file for the points
//     char output_filename[4096];
//     sprintf(output_filename, "%s/%s/%s-skeleton-%04ldz-%04ldy-%04ldx-ID-%012ld.pts", skeleton_directory, prefix, prefix, block_z, block_y, block_x ,segment_ID);
//
//     FILE *wfp = fopen(output_filename, "wb");
//     if (!wfp) { fprintf(stderr, "Failed to open %s\n", output_filename); exit(-1); }
//
//     // write the number of elements
//     if (fwrite(&(volumesize[OR_Z]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
//     if (fwrite(&(volumesize[OR_Y]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
//     if (fwrite(&(volumesize[OR_X]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
//     if (fwrite(&(input_blocksize[OR_Z]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
//     if (fwrite(&(input_blocksize[OR_Y]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
//     if (fwrite(&(input_blocksize[OR_X]), sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
//     if (fwrite(&segment_ID, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
//     if (fwrite(&num, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
//
//     // write the widths to file
//     char widths_filename[4096];
//     sprintf(widths_filename, "widths/%s/%06ld.pts", prefix, segment_ID);
//
//     FILE *width_fp = fopen(widths_filename, "wb");
//     if (!width_fp) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
//
//     // write the number of elements
//     if (fwrite(&(volumesize[OR_Z]), sizeof(long), 1, width_fp) != 1) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
//     if (fwrite(&(volumesize[OR_Y]), sizeof(long), 1, width_fp) != 1) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
//     if (fwrite(&(volumesize[OR_X]), sizeof(long), 1, width_fp) != 1) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
//     if (fwrite(&num, sizeof(long), 1, width_fp) != 1) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
//
//     printf("Remaining voxels: %ld\n", num);
//
//     long index_local[num];
//     long counter_it = 0;
//
//     while (surface_voxels.first != NULL) {
//         // get the surface voxels
//         ListElement *LE = (ListElement *) surface_voxels.first;
//
//         long index = LE->iv;
//         float width = widths[index];
//
//         // get the coordinates for this skeleton point in the global frame
//         long iz_local = LE->iz - 1;
//         long iy_local = LE->iy - 1;
//         long ix_local = LE->ix - 1;
//
//         long iz_global = iz_local + block_z*input_blocksize[OR_Z];
//         long iy_global = iy_local + block_y*input_blocksize[OR_Y];
//         long ix_global = ix_local + block_x*input_blocksize[OR_X];
//
//         // check that indices are not out of volume size
//         if (iz_global>=volumesize[OR_Z] ||  iy_global>=volumesize[OR_Y] || ix_global>=volumesize[OR_X]){
//           throw std::invalid_argument("Output global indices outside of volumesize!");
//         }
//
//         if (iz_local>=input_blocksize[OR_Z] ||  iy_local>=input_blocksize[OR_Y] || ix_local>=input_blocksize[OR_X]){
//           throw std::invalid_argument("Output local indices outside of blocksize_input!");
//         }
//
//         // TODO: transform to global frame here
//         long iv_local = iz_local * input_blocksize[OR_X] * input_blocksize[OR_Y] + iy_local * input_blocksize[OR_X] + ix_local;
//         long iv_global = iz_global * volumesize[OR_X] * volumesize[OR_Y] + iy_global * volumesize[OR_X] + ix_global;
//
//         // write to zmax border file, if on border
//         if (iz_local==(input_blocksize[OR_Z]-1)){
//
//             zmax_iy_local.push_back(iy_local);
//             zmax_ix_local.push_back(ix_local);
//             zmax_segment_ID.push_back(segment_ID);
//
//             // std::cout << "added anchor point - iz_local, inp_blocksize[OR_Z]-1: " << iz_local << "," << (input_blocksize[OR_Z]-1) << std::endl;
//
//         }
//
//         if (fwrite(&iv_global, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
//         if (fwrite(&iv_global, sizeof(long), 1, width_fp) != 1) { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
//         if (fwrite(&width, sizeof(float), 1, width_fp) != 1)  { fprintf(stderr, "Failed to write to %s\n", widths_filename); exit(-1); }
//
//         index_local[counter_it]=iv_local;
//
//         // remove this voxel
//         RemoveSurfaceVoxel(LE);
//         counter_it++;
//     }
//
//     for (int j=0; j<num; j++){
//         if (fwrite(&index_local[j], sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
//     }
//
//     // close the I/O files
//     fclose(wfp);
//     fclose(width_fp);
//
//     double total_time = (double) (clock() - start_time) / CLOCKS_PER_SEC;
//
//     char time_filename[4096];
//     sprintf(time_filename, "running_times/skeletons/%s-%06ld.time", prefix, segment_ID);
//
//     FILE *tfp = fopen(time_filename, "wb");
//     if (!tfp) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }
//
//     // write the number of points and the total time to file
//     if (fwrite(&initial_points, sizeof(long), 1, tfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }
//     if (fwrite(&total_time, sizeof(double), 1, tfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }
//
//
//     // close file
//     fclose(tfp);
//
//
// }
//
// static void writeZmaxBlock (const char *prefix, std::vector<long> &zmax_iy_local, std::vector<long> &zmax_ix_local, std::vector<long> &zmax_segment_ID)
// {
//
//   // write the zmax border points to a file (so far stored in vectors)
//   // create an output file that saves th epoints on the positive z bounday, as (global index, local index, segment_ID)*n_points , npoints
//   char output_filename_zmax[4096];
//   sprintf(output_filename_zmax, "%s/%s/%s-borderZMax-%04ldz-%04ldy-%04ldx.pts", skeleton_directory, prefix, prefix, block_z, block_y, block_x);
//
//   FILE *zmaxfp = fopen(output_filename_zmax, "wb");
//   if (!zmaxfp) { fprintf(stderr, "Failed to open %s\n", output_filename_zmax); exit(-1); }
//
//   if (!((zmax_iy_local.size()==zmax_ix_local.size())&&(zmax_ix_local.size()==zmax_segment_ID.size()))){
//     throw std::invalid_argument("z vectors not of same size!");
//   }
//
//   long n_points_zmax = zmax_iy_local.size();
//
//   std::cout << "Anchor points found: " << n_points_zmax << std::endl;
//
//   if (fwrite(&n_points_zmax, sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
//
//   for (int pos = 0; pos < n_points_zmax; pos++){
//     if (fwrite(&zmax_iy_local[pos], sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
//     if (fwrite(&zmax_ix_local[pos], sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
//     if (fwrite(&zmax_segment_ID[pos], sizeof(long), 1, zmaxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_zmax); exit(-1); }
//   }
//
//   fclose(zmaxfp);
// }
//
// static int ReadSynapses(const char *prefix)
// {
//
//   // read the synapses
//   char synapse_filename[4096];
//   snprintf(synapse_filename, 4096, "%s/%s/%s-synapses-%04ldz-%04ldy-%04ldx.pts", synapses_directory, prefix, prefix, block_z, block_y, block_x);
//
//   FILE *fp = fopen(synapse_filename, "rb");
//   if (!fp) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//
//   long z_input_volume_size, y_input_volume_size, x_input_volume_size;
//   long z_input_block_size, y_input_block_size, x_input_block_size;
//
//   if (fread(&z_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//   if (fread(&y_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//   if (fread(&x_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//
//   if (z_input_volume_size != volumesize[OR_Z]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//   if (y_input_volume_size != volumesize[OR_Y]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//   if (x_input_volume_size != volumesize[OR_X]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//
//   if (fread(&z_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//   if (fread(&y_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//   if (fread(&x_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//
//   if (z_input_block_size != input_blocksize[OR_Z]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//   if (y_input_block_size != input_blocksize[OR_Y]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//   if (x_input_block_size != input_blocksize[OR_X]) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//
//   long nneurons;
//   if (fread(&nneurons, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//
//   for (long iv = 0; iv < nneurons; ++iv) {
//       // get the label and number of synapses
//       long segment_ID;
//       long nsynapses;
//
//       if (fread(&segment_ID, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//       if (fread(&nsynapses, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//
//       // synapses[segment_ID] = std::vector<long>();
//
//       // ignore the global coordinates
//       for (long is = 0; is < nsynapses; ++is) {
//           long dummy_index;
//           if (fread(&dummy_index, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//       }
//
//       // set block indexing parameters
//       long nentries = padded_blocksize[OR_Z] * padded_blocksize[OR_Y] * padded_blocksize[OR_X];
//       long sheet_size = padded_blocksize[OR_Y] * padded_blocksize[OR_X];
//       long row_size = padded_blocksize[OR_X];
//       long infinity = padded_blocksize[OR_Z] * padded_blocksize[OR_Z] + padded_blocksize[OR_Y] * padded_blocksize[OR_Y] + padded_blocksize[OR_X] * padded_blocksize[OR_X];
//
//       // add the local coordinates (offset by one in each direction)
//       for (long is = 0; is < nsynapses; ++is) {
//           long linear_index;
//           if (fread(&linear_index, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", synapse_filename); return 0; }
//
//           long iz = linear_index / (input_blocksize[OR_Y] * input_blocksize[OR_X]);
//           long iy = (linear_index - iz * (input_blocksize[OR_Y] * input_blocksize[OR_X])) / input_blocksize[OR_X];
//           long ix = linear_index % input_blocksize[OR_X];
//
//           //  pad the location by one
//           iz += 1; iy += 1; ix += 1;
//
//           // find the new voxel index
//           long iv = IndicesToIndex(ix, iy, iz, sheet_size, row_size);
//
//           Pointclouds[segment_ID][iv] = 3;
//       }
//
//   }
//
//   // close file
//   fclose(fp);
//
//   return 1;
//
// }
//
// static int ReadAnchorpoints(const char *prefix)
// {
//   // make filename of adjacent z lock (in negative direction)
//   char output_filename_zmax[4096];
//   sprintf(output_filename_zmax, "%s/%s/%s-borderZMax-%04ldz-%04ldy-%04ldx.pts", skeleton_directory, prefix, prefix, block_z-1, block_y, block_x);
//
//   std::cout << "Reading Anchor points from: " << output_filename_zmax << std::endl;
//
//   FILE *fpzmax = fopen(output_filename_zmax, "rb");
//   if (!fpzmax) { fprintf(stderr, "Failed to read %s.\n", output_filename_zmax); return 0; }
//
//   long npoints;
//   if (fread(&npoints, sizeof(long), 1, fpzmax) != 1) { fprintf(stderr, "Failed to read %s.\n", output_filename_zmax); return 0; }
//
//   std::cout << "number of anchor points: " << npoints << std::endl;
//
//   for (long pos = 0; pos < npoints; pos++) {
//
//       // get the label and number of synapses
//       long iy_local;
//       long ix_local;
//       long segment_ID;
//
//       if (fread(&iy_local, sizeof(long), 1, fpzmax) != 1) { fprintf(stderr, "Failed to read %s.\n", output_filename_zmax); return 0; }
//       if (fread(&ix_local, sizeof(long), 1, fpzmax) != 1) { fprintf(stderr, "Failed to read %s.\n", output_filename_zmax); return 0; }
//       if (fread(&segment_ID, sizeof(long), 1, fpzmax) != 1) { fprintf(stderr, "Failed to read %s.\n", output_filename_zmax); return 0; }
//
//       // populate offset
//       long iz_padded = 1;
//       long iy_padded = iy_local + 1;
//       long ix_padded = ix_local + 1;
//
//       nentries = padded_blocksize[OR_Z] * padded_blocksize[OR_Y] * padded_blocksize[OR_X];
//       sheet_size = padded_blocksize[OR_Y] * padded_blocksize[OR_X];
//       row_size = padded_blocksize[OR_X];
//       infinity = padded_blocksize[OR_Z] * padded_blocksize[OR_Z] + padded_blocksize[OR_Y] * padded_blocksize[OR_Y] + padded_blocksize[OR_X] * padded_blocksize[OR_X];
//
//       //  find the new voxel index
//       long iv = IndicesToIndex(ix_padded, iy_padded, iz_padded, sheet_size, row_size);
//
//       Pointclouds[segment_ID][iv] = 3;
//
//       std::cout << "added anchor point at: " << iy_local << "," << ix_local << "segment ID: " << segment_ID << std::endl;
//
//   }
//
//   // close file
//   fclose(fpzmax);
//
//   return 1;
//
// }

// void CppSkeletonGeneration(const char *prefix, const char *lookup_table_directory, long *inp_labels)
// {
//     // initialize and clear set to hold all IDs that are present in this block
//     IDs_in_block = std::unordered_set<long>();
//     std::unordered_set<long> IDs_to_process;
//
//     // insert IDs that should be processed
//     IDs_to_process.insert({48,71,72,73,91,111,124,137,158,209});
//
//     // create vectors to store border points at zmax direction (size unknown so far)
//     std::vector<long> zmax_iy_local = std::vector<long>();
//     std::vector<long> zmax_ix_local = std::vector<long>();
//     std::vector<long> zmax_segment_ID = std::vector<long>();
//
//     // retrive points clouds from h5 file
//     CppPopulatePointCloudFromH5(inp_labels);
//
//     // read all synapses for this block
//     if (!ReadSynapses(prefix)) exit(-1);
//
//     // read all anchor points for this block
//     if (!ReadAnchorpoints(prefix)) exit(-1);
//
//     // initialize all of the lookup tables
//     InitializeLookupTables(lookup_table_directory);
//
//     // initialize variable segment_ID, that holds the segment_ID which is currently processed
//     long segment_ID;
//
//     // create iterator over set
//     std::unordered_set<long>::iterator itr = IDs_to_process.begin();
//
//     // iterate over all elements in this set and compute and save their skeletons
//     long loop_executions = 0;
//
//     while (itr != IDs_to_process.end())
//     {
//       // create (and clear) the global variables
//       segment = std::unordered_map<long, short>();
//       borderpoints_segment = std::unordered_map<long,std::unordered_set<long>>();
//       widths = std::unordered_map<long, float>();
//
//       // Reset surface voxels list (TODO: is this correct?)
//       surface_voxels.first = NULL;
//       surface_voxels.last = NULL;
//
//       // start timing statistics
//       clock_t start_time = clock();
//
//       // set segment_ID to current ID
//       segment_ID = *itr;
//
//       std::cout << "-----------------------------------" << std::endl;
//       std::cout << "Processing segment_ID " << segment_ID << std::endl;
//
//       // set segment to the pointcloud of current segment_ID
//       segment = Pointclouds[segment_ID];
//       borderpoints_segment = borderpoints[segment_ID];
//
//       // std::cout << "value of segment 95 at 4218037: " << segment[4218037] << std::endl;
//
//       // print number of synapses in this pointcloud
//       // std::cout << "Synapses: " << synapses.size() << std::endl;
//
//       // get the number of points
//       long initial_points = segment.size();
//       printf("segment_ID %ld initial points: %ld\n", segment_ID, initial_points);
//
//       // needs to happen after PopulatePointCloud()
//       PopulateOffsets();
//
//       // call the sequential thinning algorithm
//       SequentialThinning(prefix, segment_ID);
//
//       // Write time, skeleton and width output files
//       WriteOutputfiles(prefix, segment_ID, start_time, initial_points, zmax_iy_local, zmax_ix_local, zmax_segment_ID);
//
//       // increment iterator
//       itr++;
//       loop_executions += 1;
//
//     }
//
//     writeZmaxBlock(prefix, zmax_iy_local, zmax_ix_local, zmax_segment_ID);
//
//     delete[] lut_simple;
//
// }

void CPPcreateDataBlock(const char *prefix, const char *lookup_table_directory, long *inp_labels, float input_resolution[3], long inp_blocksize[3], long volume_size[3], long block_ind[3], const char* synapses_dir, const char* somae_dir, const char* skeleton_dir){
  // create new Datablock and set the input variables

  DataBlock BlockA;
  BlockA.CppUpdateResolution(input_resolution);
  BlockA.CppUpdateBlocksize(inp_blocksize);
  BlockA.CppUpdateVolumesize(volume_size);
  BlockA.CppUpdateBlockindices(block_ind);
  BlockA.CppUpdateDirectories(synapses_dir, somae_dir, skeleton_dir);

  // insert IDs that should be processed
  BlockA.IDs_to_process.insert({91});
  BlockA.CppPopulatePointCloudFromH5(inp_labels);
  if (!BlockA.ReadSynapses(prefix)) exit(-1);
  if (!BlockA.ReadAnchorpoints(prefix)) exit(-1);
  InitializeLookupTables(lookup_table_directory);

  long initial_points = BlockA.Pointclouds[91].size();
  printf("segment_ID 91 initial points: %ld\n", initial_points);

  // create iterator over set
  std::unordered_set<long>::iterator itr = BlockA.IDs_to_process.begin();

  // iterate over all elements in this set and compute and save their skeletons
  long loop_executions = 0;

  while (itr != BlockA.IDs_to_process.end())
  {
      BlockSegment segA(*itr, BlockA);

      // needs to happen after PopulatePointCloud()
      PopulateOffsets(segA.padded_blocksize);

      // call the sequential thinning algorithm
      SequentialThinning(prefix, segment_ID);

      itr++;
  }


}
