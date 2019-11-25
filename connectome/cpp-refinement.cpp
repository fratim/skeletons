/* c++ file for running skeleton refinement algorithm */

#include "cpp-MinBinaryHeap.h"
#include "cpp-wiring.h"
#include <algorithm>
#include <unistd.h>

struct DijkstraData {
    long iv;
    DijkstraData *prev;
    double distance;
    bool visited;
};

void WriteHeader(FILE *fp, long &num);
void ReadHeader(FILE *fp, long &num);
void CppSkeletonRefinement(const char *prefix, long segment_ID_query, long block_ind_inp[3]);
void ReadSynapses(const char *prefix, std::unordered_map<long, char> &segment, std::unordered_set<long> &synapses, long segment_ID_query);
void ReadSkeleton(const char *prefix, std::unordered_map<long, char> &segment, long segment_ID_query);


float resolution[3] = {18,20,20};
long input_blocksize[3] = {512,512,512};
long padded_blocksize[3] = {514,514,514};
long volumesize[3]  = {6144,5632,5632};
long padded_volumesize[3]  = {6146,5634,5634};
long block_ind[3]  = {-1,-1,-1};
long infinity = padded_volumesize[OR_Z] * padded_volumesize[OR_Z] + padded_volumesize[OR_Y] * padded_volumesize[OR_Y] + padded_volumesize[OR_X] * padded_volumesize[OR_X];

long padded_sheet_size_block = padded_blocksize[OR_Y] * padded_blocksize[OR_X];
long padded_row_size_block = padded_blocksize[OR_X];
long input_sheet_size_block = input_blocksize[OR_Y] * input_blocksize[OR_X];
long input_row_size_block = input_blocksize[OR_X];
long padded_sheet_size_volume = padded_volumesize[OR_Y] * padded_volumesize[OR_X];
long padded_row_size_volume = padded_volumesize[OR_X];
long input_sheet_size_volume = volumesize[OR_Y] * volumesize[OR_X];
long input_row_size_volume = volumesize[OR_X];

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

void CppSkeletonRefinement(const char *prefix, long segment_ID_query, long block_ind_inp[3])
{
    // start timing statistics
    // clock_t start_time = clock();
    std::cout << "----------------------------------"<<std::endl;

    // clear the global variables
    std::unordered_map<long, char> segment = std::unordered_map<long, char>();
    std::unordered_set<long> synapses = std::unordered_set<long>();
    std::unordered_map<long, long> dijkstra_map = std::unordered_map<long, long>();

    for (long bz = 0; bz<2; bz++){
      for (long by = 0; by<2; by++){
        for (long bx = 0; bx<2; bx++){

          block_ind[0]=block_ind_inp[0]+bz;
          block_ind[1]=block_ind_inp[1]+by;
          block_ind[2]=block_ind_inp[2]+bx;

          ReadSynapses(prefix, segment, synapses, segment_ID_query);
          ReadSkeleton(prefix, segment, segment_ID_query);

          std::cout << "processing block: "<<block_ind[0]<<","<< block_ind[1]<<","<< block_ind[2]<<std::endl;
        }
      }
    }

    // set random somae
    // get the number of elements in the skeleton
    long nelements = segment.size();
    long nsynapses = synapses.size();

    std::cout << "points before refinement: "<<nelements<<std::endl;
    std::cout << "Synapses: "<<synapses.size()<<std::endl;

    if (nsynapses>0) {
      std::unordered_set<long>::iterator it2 = synapses.begin();
      segment[*it2]=4;
    }

    DijkstraData *voxel_data = new DijkstraData[nelements];
    if (!voxel_data) exit(-1);

    // initialize the priority queue
    DijkstraData tmp;
    MinBinaryHeap<DijkstraData *> voxel_heap(&tmp, (&tmp.distance), segment.size());

    // initialize all data
    long index = 0;
    for (std::unordered_map<long, char>::iterator it = segment.begin(); it != segment.end(); ++it, ++index) {
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

            std::cout << "detected soma at: " << index << std::endl;
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
                    if (!segment[neighbor_index])continue;

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
    for (std::unordered_set<long>::iterator it = synapses.begin(); it != synapses.end(); ++it) {
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

    char wiring_filename[4096];
    sprintf(wiring_filename, "skeletons/%s/%s-connectomes-ID-%012ld.pts", prefix, prefix, segment_ID_query);
    char distance_filename[4096];
    sprintf(distance_filename, "skeletons/%s/%s-distances-ID-%012ld.pts", prefix, prefix, segment_ID_query);

    FILE *wfp = fopen(wiring_filename, "wb");
    if (!wfp) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }
    FILE *dfp = fopen(distance_filename, "wb");
    if (!dfp) { fprintf(stderr, "Failed to write to %s.\n", distance_filename); exit(-1); }


    WriteHeaderSegID(wfp, nskeleton_points, segment_ID_query);
    WriteHeaderSegID(dfp, nskeleton_points, segment_ID_query);

    long pos = 0;
    long local_index[nskeleton_points];

    for (std::unordered_set<long>::iterator it = wiring_diagram.begin(); it != wiring_diagram.end(); ++it, ++pos) {
        long index_global_up = *it;
        if (fwrite(&index_global_up, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }

        // get the upsampled location
        long ix_global_up, iy_global_up, iz_global_up;
        IndexToIndices(index_global_up, ix_global_up, iy_global_up, iz_global_up, input_sheet_size_volume, input_row_size_volume);

        long iz_local_up = iz_global_up - block_ind[OR_Z]*input_blocksize[OR_Z];
        long iy_local_up = iy_global_up - block_ind[OR_Y]*input_blocksize[OR_Y];
        long ix_local_up = ix_global_up - block_ind[OR_X]*input_blocksize[OR_X];

        // std::cout<<iz_local_up<<","<<iy_local_up<<","<<ix_local_up<<std::endl;

        long index_local_up = IndicesToIndex(ix_local_up, iy_local_up, iz_local_up, input_sheet_size_block, input_row_size_block);
        local_index[pos]=index_local_up;

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
    }

    for(long i=0; i<nskeleton_points;i++){
        if (fwrite(&local_index[i], sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", wiring_filename); exit(-1); }
    }

    // write checkvalue at end
    long checkvalue = 2147483647;
    if (fwrite(&checkvalue, sizeof(long), 1, wfp) != 1) { fprintf(stderr, "Failed to write to %s\n", wiring_filename); exit(-1); }

    fclose(wfp);
    fclose(dfp);

    delete[] voxel_data;

    // double total_time = (double) (clock() - start_time) / CLOCKS_PER_SEC;

    // char time_filename[4096];
    // sprintf(time_filename, "running_times/refinement/%s-%06ld.time", prefix, segment_ID_query);
    //
    // FILE *tfp = fopen(time_filename, "wb");
    // if (!tfp) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }
    //
    // // write the number of points and the total time to file
    // if (fwrite(&nelements, sizeof(long), 1, tfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }
    // if (fwrite(&total_time, sizeof(double), 1, tfp) != 1) { fprintf(stderr, "Failed to write to %s.\n", time_filename); exit(-1); }
    //
    // // close file
    // fclose(tfp);
}

void ReadSynapses(const char *prefix, std::unordered_map<long, char> &segment, std::unordered_set<long> &synapses, long segment_ID_query)
{

  // read the synapses
  char synapse_filename[4096];
  snprintf(synapse_filename, 4096, "synapses_projected/%s/%s-synapses_projected-%04ldz-%04ldy-%04ldx.pts", prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

  FILE *fp = fopen(synapse_filename, "rb");
  if (!fp) { fprintf(stderr, "Failed to read %s.\n", synapse_filename);  exit(-1);}

  long nneurons;
  ReadHeader(fp, nneurons);

  for (long iv = 0; iv < nneurons; ++iv) {
      // get the label and number of synapses
      long segment_ID;
      long nsynapses;

      if (fread(&segment_ID, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename); exit(-1); }
      if (fread(&nsynapses, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to read %s.\n", synapse_filename);  exit(-1); }

      // read in global indices
      for (long is = 0; is < nsynapses; ++is) {
          long iv_global;
          if (fread(&iv_global, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", synapse_filename);  exit(-1);}
      }

      // ignore the local index
      for (long is = 0; is < nsynapses; ++is) {

          long iv_local;
          if (fread(&iv_local, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", synapse_filename);  exit(-1); }

          if (segment_ID==segment_ID_query){
            long iz_local, iy_local, ix_local;
            IndexToIndices(iv_local, ix_local, iy_local, iz_local, input_sheet_size_block, input_row_size_block);

            long iz_global = iz_local + block_ind[OR_Z]*input_blocksize[OR_Z];
            long iy_global = iy_local + block_ind[OR_Y]*input_blocksize[OR_Y];
            long ix_global = ix_local + block_ind[OR_X]*input_blocksize[OR_X];

            long iv_global = IndicesToIndex(ix_global, iy_global, iz_global, input_sheet_size_volume, input_row_size_volume);

            long iv_global_padded = PadIndex(iv_global, input_sheet_size_volume, input_row_size_volume, padded_sheet_size_volume, padded_row_size_volume);
            segment[iv_global_padded]=3;
            synapses.insert(iv_global_padded);
          }
      }

  }

  // close file
  fclose(fp);

}

void ReadSkeleton(const char *prefix, std::unordered_map<long, char> &segment, long segment_ID_query)
{
  // read the synapses
  char skeleton_filename[4096];
  snprintf(skeleton_filename, 4096, "skeletons/%s/%s-skeleton-%04ldz-%04ldy-%04ldx-ID-%012ld.pts", prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X], segment_ID_query);

  FILE *fp = fopen(skeleton_filename, "rb");
  if (!fp) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename);  exit(-1); }

  // get the label and number of synapses
  long input_neuron_id;
  long nskeleton_points;

  ReadHeaderSegID(fp, nskeleton_points, input_neuron_id);

  if (segment_ID_query != input_neuron_id) { fprintf(stderr, "Failed to read %s.\n", skeleton_filename);  exit(-1); }

  // read in the global coordinates
  for (long is = 0; is < nskeleton_points; ++is) {
      long iv_global;
      if (fread(&iv_global, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", skeleton_filename);  exit(-1); }

  }

  // ignore the local ones
  for (long is = 0; is < nskeleton_points; ++is) {
      long iv_local;
      if (fread(&iv_local, sizeof(long), 1, fp) != 1)  { fprintf(stderr, "Failed to read %s.\n", skeleton_filename);  exit(-1); }

      long iz_local, iy_local, ix_local;
      IndexToIndices(iv_local, ix_local, iy_local, iz_local, input_sheet_size_block, input_row_size_block);

      long iz_global = iz_local + block_ind[OR_Z]*input_blocksize[OR_Z];
      long iy_global = iy_local + block_ind[OR_Y]*input_blocksize[OR_Y];
      long ix_global = ix_local + block_ind[OR_X]*input_blocksize[OR_X];

      long iv_global = IndicesToIndex(ix_global, iy_global, iz_global, input_sheet_size_volume, input_row_size_volume);

      long iv_global_padded = PadIndex(iv_global, input_sheet_size_volume, input_row_size_volume, padded_sheet_size_volume, padded_row_size_volume);
      segment[iv_global_padded]=1;
  }

  long checkvalue;
  if (fread(&checkvalue, sizeof(long), 1, fp) != 1) { fprintf(stderr, "Failed to write to %s\n", skeleton_filename);  exit(-1); }
  if (checkvalue != 2147483647) { fprintf(stderr, "Checkvalue for %s incorrect\n", skeleton_filename); exit(-1); }

  // close file
  fclose(fp);
}
