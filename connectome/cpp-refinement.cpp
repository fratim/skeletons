/* c++ file for running skeleton refinement algorithm */

#include "cpp-MinBinaryHeap.h"
#include "cpp-wiring.h"
#include <algorithm>
#include <unistd.h>
#include <time.h>

struct DijkstraData {
    long iv;
    DijkstraData *prev;
    double distance;
    bool visited;
};

void WriteHeader(FILE *fp, long &num);
void ReadHeader(FILE *fp, long &num);
void WriteHeaderSegID(FILE *fp, long &num, long&segID);
void ReadHeaderSegID(FILE *fp, long &num, long&segID);
void CppSkeletonRefinement(const char *prefix, float input_resolution[3], long inp_blocksize[3], long volume_size[3], long block_ind_begin[3], long block_ind_end[3], const char* output_dir);
void ReadSynapses(const char *prefix, std::unordered_map<long, std::unordered_map<long, char>> &segment, std::unordered_map<long, std::unordered_set<long>> &synapses, std::unordered_set<long> IDsToProcess, long (&block_ind)[3]);
void ReadSkeleton(const char *prefix, std::unordered_map<long, std::unordered_map<long, char>> &segment, std::unordered_set<long> IDsToProcess, long (&block_ind)[3]);
void setParameters(float input_resolution[3], long inp_blocksize[3], long volume_size[3], long block_ind_inp[3], const char* output_dir);
void ReadSomaeSurface(const char *prefix, std::unordered_map<long, std::unordered_map<long, char>> &segment, std::unordered_set<long> IDsToProcess, long (&block_ind)[3]);
void WriteProjectedSynapses(const char *prefix, std::unordered_map<long, std::unordered_set<long>> &synapses);

float resolution[3] = {-1,-1,-1};
long input_blocksize[3] = {-1,-1,-1};
long padded_blocksize[3] = {-1,-1,-1};
long volumesize[3]  = {-1,-1,-1};
long padded_volumesize[3]  = {-1,-1,-1};

long input_sheet_size_block = -1;
long input_row_size_block = -1;
long padded_sheet_size_block = -1;
long padded_row_size_block = -1;

long input_sheet_size_volume = -1;
long input_row_size_volume = -1;
long padded_sheet_size_volume = -1;
long padded_row_size_volume = -1;

long infinity = -1;

long block_search_start[3] = {-1,-1,-1};
long block_search_end[3] = {-1,-1,-1};

const char *output_directory;

void setParameters(float input_resolution[3], long inp_blocksize[3], long volume_size[3], long block_ind_begin[3], long block_ind_end[3], const char* output_dir)
{
  resolution[OR_Z] = input_resolution[OR_Z];
  resolution[OR_Y] = input_resolution[OR_Y];
  resolution[OR_X] = input_resolution[OR_X];

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

  volumesize[OR_Z] = volume_size[OR_Z];
  volumesize[OR_Y] = volume_size[OR_Y];
  volumesize[OR_X] = volume_size[OR_X];

  std::cout << "Volumesize is: " << volumesize[OR_Z] << ","<< volumesize[OR_Y] << ","<< volumesize[OR_X] << std::endl;

  padded_volumesize[OR_Z] = volume_size[OR_Z]+2;
  padded_volumesize[OR_Y] = volume_size[OR_Y]+2;
  padded_volumesize[OR_X] = volume_size[OR_X]+2;

  input_sheet_size_volume = volumesize[OR_Y] * volumesize[OR_X];
  input_row_size_volume = volumesize[OR_X];
  padded_sheet_size_volume = padded_volumesize[OR_Y] * padded_volumesize[OR_X];
  padded_row_size_volume = padded_volumesize[OR_X];

  block_search_start[OR_Z] = block_ind_begin[OR_Z];
  block_search_start[OR_Y] = block_ind_begin[OR_Y];
  block_search_start[OR_X] = block_ind_begin[OR_X];

  block_search_end[OR_Z] = block_ind_end[OR_Z];
  block_search_end[OR_Y] = block_ind_end[OR_Y];
  block_search_end[OR_X] = block_ind_end[OR_X];

  infinity = padded_volumesize[OR_Z] * padded_volumesize[OR_Z] + padded_volumesize[OR_Y] * padded_volumesize[OR_Y] + padded_volumesize[OR_X] * padded_volumesize[OR_X];

  output_directory = output_dir;

  std::cout << "Block indices start set to: " << block_search_start[OR_Z] << "," << block_search_start[OR_Y] << "," << block_search_start[OR_X] << std::endl;


}

void CppSkeletonRefinement(const char *prefix, float input_resolution[3], long inp_blocksize[3], long volume_size[3], long block_ind_begin[3], long block_ind_end[3], const char* output_dir)
{
  // start timing statistics
  double time_total = 0;
  double time_read_IDstoProcess = 0;
  double time_read_Synapses = 0;
  double time_read_Skeleton = 0;
  double time_read_SomaeSurface = 0;
  double time_dijkstra = 0;
  double time_write = 0;

  clock_t start_time_total = 0;
  clock_t start_time_read_IDstoProcess = 0;
  clock_t start_time_read_Synapses = 0;
  clock_t start_time_read_Skeleton = 0;
  clock_t start_time_read_SomaeSurface = 0;
  clock_t start_time_dijkstra = 0;
  clock_t start_time_write = 0;

  start_time_total = clock();


  setParameters(input_resolution, inp_blocksize, volume_size, block_ind_begin, block_ind_end, output_dir);

  std::unordered_set<long> IDsToProcess = std::unordered_set<long>();

  start_time_read_IDstoProcess = clock();
  for (long bz = block_search_start[OR_Z]; bz<=block_search_end[OR_Z]; bz++){
    for (long by =  block_search_start[OR_Y]; by<=block_search_end[OR_Y]; by++){
      for (long bx =  block_search_start[OR_X]; bx<=block_search_end[OR_X]; bx++){


        char output_filename_IDs[4096];
        sprintf(output_filename_IDs, "%s/output-%04ldz-%04ldy-%04ldx/segment_IDs/%s/%s-IDs_processed-%04ldz-%04ldy-%04ldx.pts", output_directory, bz, by, bx, prefix, prefix, bz,by,bx);

        FILE *fpid = fopen(output_filename_IDs, "rb");
        if (!fpid) { fprintf(stderr, "Failed to open %s\n", output_filename_IDs); exit(-1); }

        long n_ids;
        ReadHeader(fpid, n_ids);

        for (long pos = 0; pos<n_ids; pos++) {

          long seg_ID;
          if (fread(&seg_ID, sizeof(long), 1, fpid) != 1) { fprintf(stderr, "Failed to read to IDs in Block to %s \n", output_filename_IDs); exit(-1); }
          IDsToProcess.insert(seg_ID);
        }
        fclose(fpid);

      }
    }
  }

  time_read_IDstoProcess += (double) (clock() - start_time_read_IDstoProcess) / CLOCKS_PER_SEC;

  std::unordered_map<long, std::unordered_map<long, char>> segment = std::unordered_map<long, std::unordered_map<long, char>>();
  std::unordered_map<long, std::unordered_set<long>> synapses = std::unordered_map<long, std::unordered_set<long>>();

  for (long bz = block_search_start[OR_Z]; bz<=block_search_end[OR_Z]; bz++){
    for (long by =  block_search_start[OR_Y]; by<=block_search_end[OR_Y]; by++){
      for (long bx =  block_search_start[OR_X]; bx<=block_search_end[OR_X]; bx++){

        long block_ind[3];

        block_ind[0]=bz;
        block_ind[1]=by;
        block_ind[2]=bx;

        start_time_read_Synapses = clock();
        ReadSynapses(prefix, segment, synapses, IDsToProcess, block_ind);
        time_read_Synapses += (double) (clock() - start_time_read_Synapses) / CLOCKS_PER_SEC;
        start_time_read_Skeleton = clock();
        ReadSkeleton(prefix, segment, IDsToProcess, block_ind);
        time_read_Skeleton += (double) (clock() - start_time_read_Skeleton) / CLOCKS_PER_SEC;
        start_time_read_SomaeSurface = clock();
        ReadSomaeSurface(prefix, segment, IDsToProcess, block_ind);
        time_read_SomaeSurface += (double) (clock() - start_time_read_SomaeSurface) / CLOCKS_PER_SEC;

      }
    }
  }

  for (std::unordered_set<long>::iterator iter = IDsToProcess.begin(); iter != IDsToProcess.end(); ++iter){

    long ID_query = *iter;
    std::cout << "----------------------------------"<<std::endl;
    std::cout << "processing: "<<ID_query<<std::endl;

    std::unordered_map<long, long> dijkstra_map = std::unordered_map<long, long>();

    start_time_dijkstra = clock();

    long nelements = segment[ID_query].size();

    std::cout << "points before refinement: "<<nelements<<std::endl;
    std::cout << "Synapses: "<<synapses[ID_query].size()<<std::endl;

    DijkstraData *voxel_data = new DijkstraData[nelements];
    if (!voxel_data) exit(-1);

    // initialize the priority queue
    DijkstraData tmp;
    MinBinaryHeap<DijkstraData *> voxel_heap(&tmp, (&tmp.distance), segment[ID_query].size());

    // initialize all data
    long index = 0;
    for (std::unordered_map<long, char>::iterator it = segment[ID_query].begin(); it != segment[ID_query].end(); ++it, ++index) {
      voxel_data[index].iv = it->first;
      voxel_data[index].prev = NULL;
      voxel_data[index].distance = infinity;
      voxel_data[index].visited = false;
      dijkstra_map[it->first] = index;

      // this is the soma
      if (it->second == 4) {
        // insert the source into the heap
        voxel_data[index].distance = 0.0;
        voxel_data[index].visited = true;
        voxel_heap.Insert(index, &(voxel_data[index]));
      }
    }

    // visit all vertices
    long voxel_index;
    while (!voxel_heap.IsEmpty()) {
      DijkstraData *current = voxel_heap.DeleteMin();
      voxel_index = current->iv;

      // visit all 26 neighbors of this index
      long ix, iy, iz;
      IndexToIndices(voxel_index, ix, iy, iz, padded_sheet_size_volume, padded_row_size_volume);

      for (long iw = iz - 1; iw <= iz + 1; ++iw) {
        for (long iv = iy - 1; iv <= iy + 1; ++iv) {
          for (long iu = ix - 1; iu <= ix + 1; ++iu) {
            // get the linear index for this voxel
            long neighbor_index = IndicesToIndex(iu, iv, iw, padded_sheet_size_volume, padded_row_size_volume);

            // skip if background
            if (!segment[ID_query][neighbor_index])continue;

            // get the corresponding neighbor data
            long dijkstra_index = dijkstra_map[neighbor_index];
            DijkstraData *neighbor_data = &(voxel_data[dijkstra_index]);

            // find the distance between these voxels
            long deltaz = resolution[OR_Z] * (iw - iz);
            long deltay = resolution[OR_Y] * (iv - iy);
            long deltax = resolution[OR_X] * (iu - ix);

            // get the distance between (ix, iy, iz) and (iu, iv, iw)
            double distance = sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);

            // get the distance to get to this voxel through the current voxel (requires a penalty for visiting this voxel)
            double distance_through_current = current->distance + distance;
            double distance_without_current = neighbor_data->distance;

            if (!neighbor_data->visited) {
              neighbor_data->prev = current;
              neighbor_data->distance = distance_through_current;
              neighbor_data->visited = true;
              voxel_heap.Insert(dijkstra_index, neighbor_data);
            }
            else if (distance_through_current < distance_without_current) {
              neighbor_data->prev = current;
              neighbor_data->distance = distance_through_current;
              voxel_heap.DecreaseKey(dijkstra_index, neighbor_data);
            }
          }
        }
      }
    }

    std::unordered_set<long> wiring_diagram = std::unordered_set<long>();

    // go through all of the synapses and add all of the skeleton points to the source
    for (std::unordered_set<long>::iterator it = synapses[ID_query].begin(); it != synapses[ID_query].end(); ++it) {
      // get the voxel and corresponding entry in the dijkstra data frame
      long voxel_index = *it;
      long dijkstra_index = dijkstra_map[voxel_index];

      DijkstraData *data = &(voxel_data[dijkstra_index]);

      while (data != NULL) {
        // add to the list of skeleton points
        long iv = data->iv;

        // convert to unpadded coordinates
        long ix, iy, iz;
        IndexToIndices(iv, ix, iy, iz, padded_sheet_size_volume, padded_row_size_volume);
        // unpad x, y, and z
        ix -= 1; iy -= 1; iz -= 1;
        // reconvert to linear coordinates
        long iv_unpadded = IndicesToIndex(ix, iy, iz, input_sheet_size_volume, input_row_size_volume);

        wiring_diagram.insert(iv_unpadded);

        data = data->prev;
      }
    }

    long nskeleton_points = wiring_diagram.size();
    std::cout << "points after refinement: " << nskeleton_points << std::endl;

    time_dijkstra += (double) (clock() - start_time_dijkstra) / CLOCKS_PER_SEC;
    start_time_write = clock();

    char wiring_filename[4096];
    sprintf(wiring_filename, "%s/skeletons/%s/%s-connectomes-ID-%012ld.pts", output_directory, prefix, prefix, ID_query);
    char distance_filename[4096];
    sprintf(distance_filename, "%s/skeletons/%s/%s-distances-ID-%012ld.pts", output_directory, prefix, prefix, ID_query);

    FILE *wfp = fopen(wiring_filename, "wb");
    if (!wfp) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }
    FILE *dfp = fopen(distance_filename, "wb");
    if (!dfp) { fprintf(stderr, "Failed to write to %s.\n", distance_filename); exit(-1); }

    WriteHeaderSegID(wfp, nskeleton_points, ID_query);
    WriteHeaderSegID(dfp, nskeleton_points, ID_query);

    long pos = 0;
    long checksum_con = 0;
    long checksum_dis = 0;

    for (std::unordered_set<long>::iterator it = wiring_diagram.begin(); it != wiring_diagram.end(); ++it, ++pos) {
      long index_global_up = *it;
      if (fwrite(&index_global_up, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }
      checksum_con += index_global_up;
      // get the upsampled location
      long ix_global_up, iy_global_up, iz_global_up;
      IndexToIndices(index_global_up, ix_global_up, iy_global_up, iz_global_up, input_sheet_size_volume, input_row_size_volume);

      // need to repad everything
      long iz_padded_global = iz_global_up +1;
      long iy_padded_global = iy_global_up +1;
      long ix_padded_global = ix_global_up +1;

      long padded_index_global = IndicesToIndex(ix_padded_global, iy_padded_global, iz_padded_global, padded_sheet_size_volume, padded_row_size_volume);

      // get the corresponding neighbor data
      long dijkstra_index = dijkstra_map[padded_index_global];
      DijkstraData *dijkstra_data = &(voxel_data[dijkstra_index]);
      double distance = dijkstra_data->distance;

      if (fwrite(&voxel_index, sizeof(long), 1, dfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", distance_filename); exit(-1); }
      if (fwrite(&distance, sizeof(double), 1, dfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", distance_filename); exit(-1); }
      checksum_dis += voxel_index;
      checksum_dis += distance;

    }

    // write checkvalue at end
    if (fwrite(&checksum_con, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", wiring_filename); exit(-1); }
    if (fwrite(&checksum_dis, sizeof(long), 1, dfp) != 1) { fprintf(stderr, "Failed to write to %s\n", distance_filename); exit(-1); }

    fclose(wfp);
    fclose(dfp);

    delete[] voxel_data;

    time_write += (double) (clock() - start_time_write) / CLOCKS_PER_SEC;

  }

  time_total += (double) (clock() - start_time_total) / CLOCKS_PER_SEC;

  {
    char output_filename[4096];
    sprintf(output_filename, "%s/running_times/%s/%s-refinement_times.pts", output_dir, prefix, prefix);
    std::cout << "Writing time for block to : " << output_filename << std::endl;
    FILE * fptime = fopen (output_filename,"a");
    fprintf(fptime,"time_total, time_read_IDstoProcess, time_read_Synapses, time_read_Skeleton, time_read_SomaeSurface, time_dijkstra, time_write \n");
    fprintf(fptime,"%8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, \n",time_total, time_read_IDstoProcess, time_read_Synapses, time_read_Skeleton, time_read_SomaeSurface, time_dijkstra, time_write);
    fclose(fptime);
  }

  WriteProjectedSynapses(prefix, synapses);

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

  if (fread(&z_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read rty1.\n"); exit(-1); }
  if (fread(&y_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read rty2.\n"); exit(-1); }
  if (fread(&x_input_volume_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read rty3.\n"); exit(-1); }

  if (z_input_volume_size != volumesize[OR_Z]) { fprintf(stderr, "Volume Size not equal to input volume size.\n"); exit(-1); }
  if (y_input_volume_size != volumesize[OR_Y]) { fprintf(stderr, "Volume Size not equal to input volume size.\n"); exit(-1); }
  if (x_input_volume_size != volumesize[OR_X]) { fprintf(stderr, "Volume Size not equal to input volume size.\n"); exit(-1); }

  if (fread(&z_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read rty4.\n"); exit(-1); }
  if (fread(&y_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read rty5.\n"); exit(-1); }
  if (fread(&x_input_block_size, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read rty6.\n"); exit(-1); }

  if (z_input_block_size != input_blocksize[OR_Z]) { fprintf(stderr, "Block Size not equal to input block size.\n"); exit(-1); }
  if (y_input_block_size != input_blocksize[OR_Y]) { fprintf(stderr, "Block Size not equal to input block size.\n"); exit(-1); }
  if (x_input_block_size != input_blocksize[OR_X]) { fprintf(stderr, "Block Size not equal to input block size.\n"); exit(-1); }

  if (fread(&num, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read rty7.\n"); exit(-1); }

}

void ReadHeaderSegID(FILE *fp, long &num, long&segID)
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

  if (fread(&segID, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read rty.\n"); exit(-1); }
  if (fread(&num, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read rty.\n"); exit(-1); }

}

void WriteHeaderSegID(FILE *fp, long &num, long&segID)
{
  int check = 0;
  int size_l = sizeof(long);

  check += fwrite(&(volumesize[OR_Z]), size_l, 1, fp);
  check += fwrite(&(volumesize[OR_Y]), size_l, 1, fp);
  check += fwrite(&(volumesize[OR_X]), size_l, 1, fp);
  check += fwrite(&(input_blocksize[OR_Z]), size_l, 1, fp);
  check += fwrite(&(input_blocksize[OR_Y]), size_l, 1, fp);
  check += fwrite(&(input_blocksize[OR_X]), size_l, 1, fp);
  check += fwrite(&segID, size_l, 1, fp);
  check += fwrite(&num, size_l, 1, fp);

  if (check != 8) { fprintf(stderr, "Failed to write file in writeheader\n"); exit(-1); }
}

void ReadSynapses(const char *prefix, std::unordered_map<long, std::unordered_map<long, char>> &segment, std::unordered_map<long, std::unordered_set<long>> &synapses, std::unordered_set<long> IDsToProcess, long (&block_ind)[3])
{

  // read the synapses
  char synapse_filename[4096];
  snprintf(synapse_filename, 4096, "%s/output-%04ldz-%04ldy-%04ldx/synapses_projected/%s/%s-synapses_projected-%04ldz-%04ldy-%04ldx.pts", output_directory, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X], prefix,prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

  FILE *fp = fopen(synapse_filename, "rb");
  if (fp) {

    long nneurons;
    ReadHeader(fp, nneurons);
    long checksum = 0;

    for (long iv = 0; iv < nneurons; ++iv) {
        // get the label and number of synapses
        long segment_ID;
        long nsynapses;

        if (fread(&segment_ID, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); exit(-1); }
        if (fread(&nsynapses, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename);  exit(-1); }

        bool isQuery = 0;
        if (IDsToProcess.find(segment_ID) != IDsToProcess.end()) isQuery = 1;

        // read in global indices
        for (long is = 0; is < nsynapses; ++is) {
            long up_iv_global;
            if (fread(&up_iv_global, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", synapse_filename);  exit(-1);}
            checksum += up_iv_global;

            if (isQuery){
              long p_iv_global = PadIndex(up_iv_global, input_sheet_size_volume, input_row_size_volume, padded_sheet_size_volume, padded_row_size_volume);
              segment[segment_ID][p_iv_global]=3;
              synapses[segment_ID].insert(p_iv_global);
            }

        }

        // ignore the local index
        for (long is = 0; is < nsynapses; ++is) {
            long iv_local;
            if (fread(&iv_local, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", synapse_filename);  exit(-1); }
            checksum += iv_local;
        }

    }

    long checksum_read;
    if (fread(&checksum_read, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to write to %s\n", synapse_filename);  exit(-1); }
    if (checksum != checksum_read) { fprintf(stderr, "Checksum for %s incorrect\n", synapse_filename); exit(-1); }

  }


  // close file
  fclose(fp);

}

void ReadSomaeSurface(const char *prefix, std::unordered_map<long, std::unordered_map<long, char>> &segment, std::unordered_set<long> IDsToProcess, long (&block_ind)[3])
{

  // read the synapses
  char somaefilename[4096];
  snprintf(somaefilename, 4096, "%s/output-%04ldz-%04ldy-%04ldx/somae_surface/%s/%s-somae_surfaces-%04ldz-%04ldy-%04ldx.pts", output_directory, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X], prefix,prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

  FILE *fp = fopen(somaefilename, "rb");
  if (fp) {

    long nneurons;
    ReadHeader(fp, nneurons);
    long checksum = 0;

    for (long iv = 0; iv < nneurons; ++iv) {
        // get the label and number of synapses
        long segment_ID;
        long n_points;

        if (fread(&segment_ID, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", somaefilename); exit(-1); }
        if (fread(&n_points, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", somaefilename);  exit(-1); }

        bool isQuery = 0;
        if (IDsToProcess.find(segment_ID) != IDsToProcess.end()) isQuery = 1;

        // read in global indices
        for (long is = 0; is < n_points; ++is) {
            long up_iv_global;
            if (fread(&up_iv_global, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", somaefilename);  exit(-1);}
            checksum += up_iv_global;

            if (isQuery){
              long p_iv_global = PadIndex(up_iv_global, input_sheet_size_volume, input_row_size_volume, padded_sheet_size_volume, padded_row_size_volume);
              segment[segment_ID][p_iv_global]=4;
            }

        }

        // ignore the local index
        for (long is = 0; is < n_points; ++is) {
            long iv_local;
            if (fread(&iv_local, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", somaefilename);  exit(-1); }
            checksum += iv_local;
        }

    }

    long checksum_read;
    if (fread(&checksum_read, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to write to %s\n", somaefilename);  exit(-1); }
    if (checksum != checksum_read) { fprintf(stderr, "Checksum for %s incorrect\n", somaefilename); exit(-1); }

  }


  // close file
  fclose(fp);

}

void ReadSkeleton(const char *prefix, std::unordered_map<long, std::unordered_map<long, char>> &segment, std::unordered_set<long> IDsToProcess, long (&block_ind)[3])
{

  for (std::unordered_set<long>::iterator iter = IDsToProcess.begin(); iter != IDsToProcess.end(); ++iter){

    long ID_query = *iter;

    // read the synapses
    char skeleton_filename[4096];
    snprintf(skeleton_filename, 4096, "%s/output-%04ldz-%04ldy-%04ldx/skeleton/%s/%s-skeleton-%04ldz-%04ldy-%04ldx-ID-%012ld.pts", output_directory, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X], prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X], ID_query);

    FILE *fp = fopen(skeleton_filename, "rb");
    if (fp) {

      // get the label and number of synapses
      long segment_ID;
      long nskeleton_points;

      ReadHeaderSegID(fp, nskeleton_points, segment_ID);
      long checksum = 0;

      if (ID_query != segment_ID) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename);  exit(-1); }

      // read in the global coordinates
      for (long is = 0; is < nskeleton_points; ++is) {
          long up_iv_global;
          if (fread(&up_iv_global, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", skeleton_filename);  exit(-1); }
          long p_iv_global = PadIndex(up_iv_global, input_sheet_size_volume, input_row_size_volume, padded_sheet_size_volume, padded_row_size_volume);
          segment[segment_ID][p_iv_global]=1;
          checksum += up_iv_global;

      }

      // ignore the local ones
      for (long is = 0; is < nskeleton_points; ++is) {
          long iv_local;
          if (fread(&iv_local, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", skeleton_filename);  exit(-1); }
          checksum += iv_local;

      }

      long checksum_read;
      if (fread(&checksum_read, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to write to %s\n", skeleton_filename);  exit(-1); }
      if (checksum != checksum_read) { fprintf(stderr, "Checksum for %s incorrect\n", skeleton_filename); exit(-1); }

      // close file
      fclose(fp);
    }
  }
}

void WriteProjectedSynapses(const char *prefix, std::unordered_map<long, std::unordered_set<long>> &synapses)
{
  //get number of anchor points
  long n_neurons = synapses.size();

  // create an output file for the points
  char output_filename[4096];
  sprintf(output_filename, "%s/synapses_projected/%s/%s-synapses_projected.pts", output_directory, prefix, prefix);

  FILE *wfp = fopen(output_filename, "wb");
  if (!wfp) { fprintf(stderr, "Failed to open %s\n", output_filename); exit(-1); }

  // write the characteristics header
  WriteHeader(wfp, n_neurons);
  long checksum = 0;

  for (std::unordered_map<long,std::unordered_set<long>>::iterator itr = synapses.begin(); itr!=synapses.end(); ++itr){
    long seg_id = itr->first;

    long n_synapses = synapses[seg_id].size();

    if (fwrite(&seg_id, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
    if (fwrite(&n_synapses, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }

    for (std::unordered_set<long>::iterator itr2 = synapses[seg_id].begin(); itr2!=synapses[seg_id].end(); ++itr2){

      long p_iv_global = *itr2;
      long up_iv_global = UnpadIndex(p_iv_global, input_sheet_size_volume, input_row_size_volume, padded_sheet_size_volume, padded_row_size_volume);
      if (fwrite(&up_iv_global, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
      checksum += up_iv_global;
    }

  }

  if (fwrite(&checksum, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename); exit(-1); }
  fclose(wfp);
}
