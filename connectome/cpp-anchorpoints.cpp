/* c++ file to extract wiring diagram */
#include <limits>
#include "cpp-wiring.h"
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <vector>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include "powersof2.h"
#include <malloc.h>
#include <string.h>
#include <queue>

#define CONTOUR_MASK 3
#define OBJECT 1

using namespace std;

class PGMImage {
 public:
   unsigned char* data;
   int width, height, depth;
   PGMImage( int width, int height, int depth );
   ~PGMImage();
   void createBorder( int bwidth, unsigned char color = 0 );
};

PGMImage::PGMImage( int width, int height, int depth ) {
 this->width = width;
 this->height = height;
 this->depth = depth;

 this->data = NULL;
 this->data = ( unsigned char* ) malloc ( width * height * sizeof( unsigned char ) );
 if ( this->data == NULL ) {
   fprintf(stderr, "PGMImage: Allocation error. Not enough memory. \n"); exit(-1);;

   return;
 }
 int size = width * height;
 for ( int i = 0; i < size; i++ ) {
   this->data[ i ] = 0;
 }


}

PGMImage::~PGMImage() {
   delete[] this->data;
}

void PGMImage::createBorder( int bwidth, unsigned char color ) {
 if ( this->data == NULL ) {
   fprintf(stderr, "ERROR: PGMImage::createBorder: data does not exist\n "); exit(-1);;
 }
 int new_width, new_height;
 new_width = this->width + 2 * bwidth;
 new_height = this->height + 2 * bwidth;
 unsigned long new_length = new_width * new_height;
 unsigned char *temp;
 temp = ( unsigned char* ) malloc ( new_length * sizeof( unsigned char ) );
 if ( !temp ) {
   fprintf(stderr, "ERROR: PGMImage::createBorder: not enough memory for temp array\n "); exit(-1);;
 }
 unsigned long orig_index, new_index;

 for ( unsigned long i = 0; i < new_length; i++ ) temp[ i ] = color;

 if ( bwidth > 0 ) {
   for ( int y = 0; y < height; y++ ) {
     for ( int x = 0; x < width; x++ ) {
       orig_index = y * width + x;
       new_index = ( y + bwidth ) * (width + 2*bwidth ) + (x + bwidth);
       temp[ new_index ] = this->data[orig_index];
     }
   }
 } else {
   for ( int y = -bwidth; y < height + bwidth; y++ ) {
     for ( int x = -bwidth; x < width + bwidth; x++ ) {
       orig_index = y * width + x;
       new_index = ( y + bwidth ) * (width + 2*bwidth ) + (x + bwidth);
       temp[ new_index ] = this->data[orig_index];
     }
   }
 }

 free ( data );
 data = temp;
 this->width += 2*bwidth;
 this->height += 2 * bwidth;

}

unsigned long palagyi_fpta( PGMImage* img, queue<unsigned long> &contour, unsigned char* lut );
void fpta_thinning( PGMImage* img, unsigned char *lut, unsigned char *lut2 );
void ThinImage(PGMImage* img, std::vector<long> &iu_centers, std::vector<long> &iv_centers);
void Writeheader(FILE *fp, long &num);
void ProcessZAnchors(const char *prefix, const char* output_dir, long *z_min_wall, long *z_max_wall);
void ProcessYAnchors(const char *prefix, const char* output_dir, long *y_min_wall, long *y_max_wall);
void ProcessXAnchors(const char *prefix, const char* output_dir, long *x_min_wall, long *x_max_wall);

long inp_volumesize[3] = {-1,-1,-1};
long block_ind[3] = {-1,-1,-1};
long inp_blocksize[3] = {-1,-1,-1};
long inp_sheetsize_block = -1;
long inp_rowsize_block = -1;

void ComputeAnchorPoints(const char *prefix, const char* output_dir, long inp_blocksize_inp[3], long blockind_inp[3], long block_ind_start_inp[3], long block_ind_end_inp[3],
        long inp_volumesize_inp[3], long *z_min_wall, long *z_max_wall, long *y_min_wall, long *y_max_wall, long *x_min_wall, long *x_max_wall){

  inp_blocksize[OR_Z] = inp_blocksize_inp[OR_Z];
  inp_blocksize[OR_Y] = inp_blocksize_inp[OR_Y];
  inp_blocksize[OR_X] = inp_blocksize_inp[OR_X];
  inp_sheetsize_block = inp_blocksize_inp[OR_Y]*inp_blocksize_inp[OR_X];
  inp_rowsize_block = inp_blocksize_inp[OR_X];
  block_ind[OR_Z] = blockind_inp[OR_Z];
  block_ind[OR_Y] = blockind_inp[OR_Y];
  block_ind[OR_X] = blockind_inp[OR_X];
  inp_volumesize[OR_Z] = inp_volumesize_inp[OR_Z];
  inp_volumesize[OR_Y] = inp_volumesize_inp[OR_Y];
  inp_volumesize[OR_X] = inp_volumesize_inp[OR_X];

  if (block_ind[OR_Z]<block_ind_end_inp[OR_Z]) ProcessZAnchors(prefix, output_dir, z_min_wall,z_max_wall);
  if (block_ind[OR_Y]<block_ind_end_inp[OR_Y]) ProcessYAnchors(prefix, output_dir, y_min_wall,y_max_wall);
  if (block_ind[OR_X]<block_ind_end_inp[OR_X]) ProcessXAnchors(prefix, output_dir, x_min_wall,x_max_wall);

}

void ProcessZAnchors(const char *prefix, const char* output_dir, long *z_min_wall, long *z_max_wall){

  long entries = (inp_blocksize[OR_Y]*inp_blocksize[OR_X]);
  std::unordered_set<long> IDs_present = std::unordered_set<long>();
  std::unordered_map <long, PGMImage*> images = std::unordered_map <long, PGMImage*>();
  std::unordered_map<long, std::vector<long>> iu_centers = std::unordered_map<long, std::vector<long>>();
  std::unordered_map<long, std::vector<long>> iv_centers = std::unordered_map<long, std::vector<long>>();

  long current_ID;
  int u_size = (int)inp_blocksize[OR_X]; // width
  int v_size = (int)inp_blocksize[OR_Y]; // height
  int depth = 1;


  // compute anchors on the zmax plane and write them to the folder and the z max adjacent folder
  for (long pos =0; pos<entries; pos++){
    if (z_max_wall[pos]==z_min_wall[pos]) current_ID = z_max_wall[pos];
    else current_ID = 0;

    if (current_ID!=0){
      if (IDs_present.find(current_ID)==IDs_present.end()){
        IDs_present.insert(current_ID);
        images[current_ID] = new PGMImage(u_size,v_size,depth);
        images[current_ID]->data[pos]=current_ID;
      }
      else{
        images[current_ID]->data[pos]=current_ID;
      }
    }
  }


  for (std::unordered_map <long, PGMImage*>::iterator iter = images.begin(); iter != images.end(); iter++){
    ThinImage(iter->second, iu_centers[iter->first], iv_centers[iter->first]);

  }
  // write to file
  {
    char output_filename_max[4096];
    sprintf(output_filename_max, "%s/output-%04ldz-%04ldy-%04ldx/anchorpoints_computed/%s/%s-Anchors_Comp_ZMax-%04ldz-%04ldy-%04ldx.pts", output_dir, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X], prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

    FILE *maxfp = fopen(output_filename_max, "wb");
    if (!maxfp) { fprintf(stderr, "Failed to open %s\n", output_filename_max); exit(-1); }

    char output_filename_min[4096];
    sprintf(output_filename_min, "%s/output-%04ldz-%04ldy-%04ldx/anchorpoints_computed/%s/%s-Anchors_Comp_ZMin-%04ldz-%04ldy-%04ldx.pts", output_dir, block_ind[OR_Z]+1, block_ind[OR_Y], block_ind[OR_X], prefix, prefix, block_ind[OR_Z]+1, block_ind[OR_Y], block_ind[OR_X]);

    FILE *minfp = fopen(output_filename_min, "wb");
    if (!minfp) { fprintf(stderr, "Failed to open %s\n", output_filename_min); exit(-1); }

    // std::cout << "printing to " << output_filename_max << std::endl;
    // std::cout << "printing to " << output_filename_min << std::endl;


    long nsegments = iu_centers.size();
    Writeheader(maxfp, nsegments);
    Writeheader(minfp, nsegments);
    long checksum_max = 0;
    long checksum_min = 0;

    for (std::unordered_map<long, std::vector<long>>::iterator iter = iu_centers.begin(); iter != iu_centers.end(); ++iter) {

      long seg_ID = iter->first;

      long n_anchors =  iu_centers[seg_ID].size();
      if (n_anchors != (long)iv_centers[seg_ID].size()) { fprintf(stderr, "Different number of entries in iu iv segment: %ld\n", seg_ID); exit(-1); }

      if (fwrite(&seg_ID, sizeof(long), 1, maxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_max); exit(-1); }
      if (fwrite(&n_anchors, sizeof(long), 1, maxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_max); exit(-1); }

      if (fwrite(&seg_ID, sizeof(long), 1, minfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_min); exit(-1); }
      if (fwrite(&n_anchors, sizeof(long), 1, minfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_min); exit(-1); }

      long local_index_max[n_anchors];
      long local_index_min[n_anchors];

      for (long pos=0; pos<n_anchors; pos++) {
        long iy = iv_centers[seg_ID][pos];
        long ix = iu_centers[seg_ID][pos];

        long iz = inp_blocksize[OR_Z]-1;
        long up_iv_local_max = IndicesToIndex(ix, iy, iz, inp_sheetsize_block, inp_rowsize_block);
        local_index_max[pos] = up_iv_local_max;
        long up_iv_global_max = IndexLocalToGlobal(up_iv_local_max, block_ind, inp_blocksize, inp_volumesize);
        if (fwrite(&up_iv_global_max, sizeof(long), 1, maxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_max); exit(-1); }
        checksum_max += (up_iv_local_max+up_iv_global_max);

        iz = 0;
        long block_ind_zmax[3] = {block_ind[OR_Z]+1,block_ind[OR_Y],block_ind[OR_X]};
        long up_iv_local_min = IndicesToIndex(ix, iy, iz, inp_sheetsize_block, inp_rowsize_block);
        local_index_min[pos] = up_iv_local_min;
        long up_iv_global_min = IndexLocalToGlobal(up_iv_local_min, block_ind_zmax, inp_blocksize, inp_volumesize);
        if (fwrite(&up_iv_global_min, sizeof(long), 1, minfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_min); exit(-1); }
        checksum_min += (up_iv_local_min+up_iv_global_min);

      }

      for (long pos=0; pos<n_anchors; pos++) {
        if (fwrite(&local_index_max[pos], sizeof(long), 1, maxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_max); exit(-1); }
        if (fwrite(&local_index_min[pos], sizeof(long), 1, minfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_min); exit(-1); }
      }



    }

    if (fwrite(&checksum_max, sizeof(long), 1, maxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_max); exit(-1); }
    if (fwrite(&checksum_min, sizeof(long), 1, minfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_min); exit(-1); }

    fclose(maxfp);
    fclose(minfp);

  }

}

void ProcessYAnchors(const char *prefix, const char* output_dir, long *y_min_wall, long *y_max_wall){

  long entries = (inp_blocksize[OR_Z]*inp_blocksize[OR_X]);
  std::unordered_set<long> IDs_present = std::unordered_set<long>();
  std::unordered_map <long, PGMImage*> images = std::unordered_map <long, PGMImage*>();
  std::unordered_map<long, std::vector<long>> iu_centers = std::unordered_map<long, std::vector<long>>();
  std::unordered_map<long, std::vector<long>> iv_centers = std::unordered_map<long, std::vector<long>>();

  long current_ID;
  int u_size = (int)inp_blocksize[OR_X]; // width
  int v_size = (int)inp_blocksize[OR_Z]; // height
  int depth = 1;

  // compute anchors on the zmax plane and write them to the folder and the z max adjacent folder
  for (long pos =0; pos<entries; pos++){
    if (y_max_wall[pos]==y_min_wall[pos]) current_ID = y_max_wall[pos];
    else current_ID = 0;

    if (current_ID!=0){
      if (IDs_present.find(current_ID)==IDs_present.end()){
        IDs_present.insert(current_ID);
        images[current_ID] = new PGMImage(u_size,v_size,depth);
        images[current_ID]->data[pos]=current_ID;
      }
      else{
        images[current_ID]->data[pos]=current_ID;
      }
    }
  }

  for (std::unordered_map <long, PGMImage*>::iterator iter = images.begin(); iter != images.end(); iter++){
    ThinImage(iter->second, iu_centers[iter->first], iv_centers[iter->first]);
  }
  // write to file
  {
    char output_filename_max[4096];
    sprintf(output_filename_max, "%s/output-%04ldz-%04ldy-%04ldx/anchorpoints_computed/%s/%s-Anchors_Comp_YMax-%04ldz-%04ldy-%04ldx.pts", output_dir, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X], prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);

    FILE *maxfp = fopen(output_filename_max, "wb");
    if (!maxfp) { fprintf(stderr, "Failed to open %s\n", output_filename_max); exit(-1); }

    char output_filename_min[4096];
    sprintf(output_filename_min, "%s/output-%04ldz-%04ldy-%04ldx/anchorpoints_computed/%s/%s-Anchors_Comp_YMin-%04ldz-%04ldy-%04ldx.pts", output_dir, block_ind[OR_Z], block_ind[OR_Y]+1, block_ind[OR_X], prefix, prefix, block_ind[OR_Z], block_ind[OR_Y]+1, block_ind[OR_X]);

    FILE *minfp = fopen(output_filename_min, "wb");
    if (!minfp) { fprintf(stderr, "Failed to open %s\n", output_filename_min); exit(-1); }

    // std::cout << "printing to " << output_filename_max << std::endl;
    // std::cout << "printing to " << output_filename_min << std::endl;

    long nsegments = iu_centers.size();
    Writeheader(maxfp, nsegments);
    Writeheader(minfp, nsegments);
    long checksum_max = 0;
    long checksum_min = 0;


    for (std::unordered_map<long, std::vector<long>>::iterator iter = iu_centers.begin(); iter != iu_centers.end(); ++iter) {

      long seg_ID = iter->first;

      long n_anchors =  iu_centers[seg_ID].size();

      long local_index_max[n_anchors];
      long local_index_min[n_anchors];

      if (n_anchors != (long)iv_centers[seg_ID].size()) { fprintf(stderr, "Different number of entries in iu iv segment: %ld\n", seg_ID); exit(-1); }

      if (fwrite(&seg_ID, sizeof(long), 1, maxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_max); exit(-1); }
      if (fwrite(&n_anchors, sizeof(long), 1, maxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_max); exit(-1); }

      if (fwrite(&seg_ID, sizeof(long), 1, minfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_min); exit(-1); }
      if (fwrite(&n_anchors, sizeof(long), 1, minfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_min); exit(-1); }

      for (long pos=0; pos<n_anchors; pos++) {
        long iz = iv_centers[seg_ID][pos];
        long ix = iu_centers[seg_ID][pos];

        long iy = inp_blocksize[OR_Y]-1;
        long up_iv_local_max = IndicesToIndex(ix, iy, iz, inp_sheetsize_block, inp_rowsize_block);
        local_index_max[pos] = up_iv_local_max;
        long up_iv_global_max = IndexLocalToGlobal(up_iv_local_max, block_ind, inp_blocksize, inp_volumesize);
        if (fwrite(&up_iv_global_max, sizeof(long), 1, maxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_max); exit(-1); }
        checksum_max += (up_iv_local_max+up_iv_global_max);


        iy = 0;
        long block_ind_ymax[3] = {block_ind[OR_Z],block_ind[OR_Y]+1,block_ind[OR_X]};
        long up_iv_local_min = IndicesToIndex(ix, iy, iz, inp_sheetsize_block, inp_rowsize_block);
        local_index_min[pos] = up_iv_local_min;
        long up_iv_global_min = IndexLocalToGlobal(up_iv_local_min, block_ind_ymax, inp_blocksize, inp_volumesize);
        if (fwrite(&up_iv_global_min, sizeof(long), 1, minfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_min); exit(-1); }
        checksum_min += (up_iv_local_min+up_iv_global_min);

      }

      for (long pos=0; pos<n_anchors; pos++) {
        if (fwrite(&local_index_max[pos], sizeof(long), 1, maxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_max); exit(-1); }
        if (fwrite(&local_index_min[pos], sizeof(long), 1, minfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_min); exit(-1); }
      }



    }

    if (fwrite(&checksum_max, sizeof(long), 1, maxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_max); exit(-1); }
    if (fwrite(&checksum_min, sizeof(long), 1, minfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_min); exit(-1); }

    fclose(maxfp);
    fclose(minfp);

  }

}

void ProcessXAnchors(const char *prefix, const char* output_dir, long *x_min_wall, long *x_max_wall){

  long entries = (inp_blocksize[OR_Y]*inp_blocksize[OR_Z]);
  std::unordered_set<long> IDs_present = std::unordered_set<long>();
  std::unordered_map <long, PGMImage*> images = std::unordered_map <long, PGMImage*>();
  std::unordered_map<long, std::vector<long>> iu_centers = std::unordered_map<long, std::vector<long>>();
  std::unordered_map<long, std::vector<long>> iv_centers = std::unordered_map<long, std::vector<long>>();

  long current_ID;
  int u_size = (int)inp_blocksize[OR_Y]; // width
  int v_size = (int)inp_blocksize[OR_Z]; // height
  int depth = 1;

  // compute anchors on the zmax plane and write them to the folder and the z max adjacent folder
  for (long pos =0; pos<entries; pos++){
    if (x_max_wall[pos]==x_min_wall[pos]) current_ID = x_max_wall[pos];
    else current_ID = 0;

    if (current_ID!=0){
      if (IDs_present.find(current_ID)==IDs_present.end()){
        IDs_present.insert(current_ID);
        images[current_ID] = new PGMImage(u_size,v_size,depth);
        images[current_ID]->data[pos]=current_ID;
      }
      else{
        images[current_ID]->data[pos]=current_ID;
      }
    }
  }

  for (std::unordered_map <long, PGMImage*>::iterator iter = images.begin(); iter != images.end(); iter++){
    ThinImage(iter->second, iu_centers[iter->first], iv_centers[iter->first]);
  }
  // write to file
  {
    char output_filename_max[4096];
    sprintf(output_filename_max, "%s/output-%04ldz-%04ldy-%04ldx/anchorpoints_computed/%s/%s-Anchors_Comp_XMax-%04ldz-%04ldy-%04ldx.pts", output_dir, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X], prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]);


    FILE *maxfp = fopen(output_filename_max, "wb");
    if (!maxfp) { fprintf(stderr, "Failed to open %s\n", output_filename_max); exit(-1); }

    char output_filename_min[4096];
    sprintf(output_filename_min, "%s/output-%04ldz-%04ldy-%04ldx/anchorpoints_computed/%s/%s-Anchors_Comp_XMin-%04ldz-%04ldy-%04ldx.pts", output_dir, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]+1, prefix, prefix, block_ind[OR_Z], block_ind[OR_Y], block_ind[OR_X]+1);

    FILE *minfp = fopen(output_filename_min, "wb");
    if (!minfp) { fprintf(stderr, "Failed to open %s\n", output_filename_min); exit(-1); }

    // std::cout << "printing to " << output_filename_max << std::endl;
    // std::cout << "printing to " << output_filename_min << std::endl;

    long nsegments = iu_centers.size();
    Writeheader(maxfp, nsegments);
    Writeheader(minfp, nsegments);
    long checksum_max = 0;
    long checksum_min = 0;

    for (std::unordered_map<long, std::vector<long>>::iterator iter = iu_centers.begin(); iter != iu_centers.end(); ++iter) {

      long seg_ID = iter->first;

      long n_anchors =  iu_centers[seg_ID].size();
      if (n_anchors != (long)iv_centers[seg_ID].size()) { fprintf(stderr, "Different number of entries in iu iv segment: %ld\n", seg_ID); exit(-1); }

      if (fwrite(&seg_ID, sizeof(long), 1, maxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_max); exit(-1); }
      if (fwrite(&n_anchors, sizeof(long), 1, maxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_max); exit(-1); }

      if (fwrite(&seg_ID, sizeof(long), 1, minfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_min); exit(-1); }
      if (fwrite(&n_anchors, sizeof(long), 1, minfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_min); exit(-1); }

      long local_index_max[n_anchors];
      long local_index_min[n_anchors];

      for (long pos=0; pos<n_anchors; pos++) {
        long iz = iv_centers[seg_ID][pos];
        long iy = iu_centers[seg_ID][pos];

        long ix = inp_blocksize[OR_X]-1;
        long up_iv_local_max = IndicesToIndex(ix, iy, iz, inp_sheetsize_block, inp_rowsize_block);
        local_index_max[pos] = up_iv_local_max;
        long up_iv_global_max = IndexLocalToGlobal(up_iv_local_max, block_ind, inp_blocksize, inp_volumesize);
        if (fwrite(&up_iv_global_max, sizeof(long), 1, maxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_max); exit(-1); }
        checksum_max += (up_iv_local_max+up_iv_global_max);

        ix = 0;
        long block_ind_xmax[3]= {block_ind[OR_Z],block_ind[OR_Y],block_ind[OR_X]+1};
        long up_iv_local_min = IndicesToIndex(ix, iy, iz, inp_sheetsize_block, inp_rowsize_block);
        local_index_min[pos] = up_iv_local_min;
        long up_iv_global_min = IndexLocalToGlobal(up_iv_local_min, block_ind_xmax, inp_blocksize, inp_volumesize);
        if (fwrite(&up_iv_global_min, sizeof(long), 1, minfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_min); exit(-1); }
        checksum_min += (up_iv_local_min+up_iv_global_min);

      }

      for (long pos=0; pos<n_anchors; pos++) {
        if (fwrite(&local_index_max[pos], sizeof(long), 1, maxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_max); exit(-1); }
        if (fwrite(&local_index_min[pos], sizeof(long), 1, minfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_min); exit(-1); }
      }

    }
    if (fwrite(&checksum_max, sizeof(long), 1, maxfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_max); exit(-1); }
    if (fwrite(&checksum_min, sizeof(long), 1, minfp) != 1) { fprintf(stderr, "Failed to write to %s\n", output_filename_min); exit(-1); }

    fclose(maxfp);
    fclose(minfp);

  }

}

void ThinImage(PGMImage* img, std::vector<long> &iu_centers, std::vector<long> &iv_centers){

  FILE* fp = NULL;

  fp = fopen( "../connectome/ronse_fpta.lut", "rb" );

  unsigned char* lut = NULL;
  lut = (unsigned char*) malloc ( 16777216 * sizeof(unsigned char) );
  if ( !lut ) {
    fprintf(stderr, "Cannot allocate memory for LUT \n"); exit(-1);
    fflush(fp);
    fclose(fp);
  }

  //fread( lut, sizeof(unsigned char), 16777216, fp );
  unsigned long i = 0;
  int read = 0;
  while( !feof(fp) ) {
    unsigned char bits;
    read += fread( &bits, sizeof(unsigned char), 1, fp );
    for ( unsigned int j = 0, k=1; j < 8; j++, k*=2, i++ ) {
      if ( (bits & k) > 0 ) lut[i] = 1; else lut[i] = 0;
    }
  }
  if (read!= 2097152) { fprintf(stderr, "Failed to read LUT \n"); exit(-1); }

  unsigned long nobject = 0;
  unsigned long index = 0;
  for ( int y = 0; y < img->height; y++ ) {
    for ( int x = 0; x < img->width; x++ ) {
      if ( img->data[ index ] != 0 ) {
        nobject++;
        img->data[index] = OBJECT;
      }
      index++;
    }
  }

  img->createBorder(3);

  fpta_thinning( img, lut, NULL );

  img->createBorder(-3);

  unsigned long nskel = 0;
  index = 0;
  for ( int y = 0; y < img->height; y++ ) {
    for ( int x = 0; x < img->width; x++ ) {
      //if ( ((img->data[ index ] & CONTOUR_MASK) != 0) || ((img->data[index] & OBJECT)) != 0 ) {
      if ( (img->data[index] & OBJECT) == OBJECT) {
        nskel++;
        iu_centers.push_back(index%img->width);
        iv_centers.push_back(index/img->width);
      } else {
      }
      index++;
    }
  }

  delete img;

  free( lut );
  fflush(fp);
  fclose(fp);


  // printf( "\n%s :: ", argv[2] );
  // printf( "Number of iterations: %d ", iter );
  // printf( "#object: %d ", nobject );
  // printf( "#skeletal: %d ", nskel );

  // printf( "\n" );
}

// iteration step
unsigned long palagyi_fpta( PGMImage* img, queue<unsigned long> &contour, unsigned char* lut ) {

  unsigned long length = contour.size();

  int w = img->width;

  int env[24] = { -w-1, -w, -w+1, 1, w+1, w, w-1,-1, -2, -w-2, -2*w-2, -2*w-1, -2*w, -2*w+1, -2*w+2, -w+2, 2, w+2, 2*w+2, 2*w+1, 2*w, 2*w-1, 2*w-2, w-2 };
  int n4[4] = { -w, 1, w, -1 };




  queue<unsigned long> deletable;

  for ( unsigned long t = 0; t < length; t++ ) {
    unsigned long p = contour.front();
    contour.pop();
    unsigned long code = 0;
    if ( img->data[p] == CONTOUR_MASK ) {
      //contour.push(p);
      for ( unsigned long i = 0, k = 1; i < 24; i++, k*=2 ) {
        if ( img->data[p+env[i]] != 0 ) {
          code |= k;
        }
      }
      int endpoint = 0;

	  	if ( endpoint == 0 ) {
				if ( lut[code] == 1 ) deletable.push(p); else contour.push(p);
	  	}
    }

  }

  length = deletable.size();

  for ( unsigned long i = 0; i < length; i++ ) {
    unsigned long p = deletable.front();
    deletable.pop();

    img->data[p] = 0;
    for ( int j = 0; j < 4; j++ ) {
      if ( img->data[p+n4[j]] == OBJECT ) {
        img->data[p+n4[j]] = CONTOUR_MASK;
        contour.push( p+n4[j] );
      }
    }

  }

  return length;

}

void fpta_thinning( PGMImage* img, unsigned char *lut, unsigned char *lut2 ) {
  queue<unsigned long> contour;
  int w = img->width, h = img->height;
  unsigned long size = w * h;
  int n4[4] = { -w, 1, w, -1 };

  for ( unsigned long i = 0; i < size; i++ ) {
    if ( img->data[i] == OBJECT ) {
      for ( int j = 0; j < 4; j++ ) {
        if ( img->data[i+n4[j]] == 0 ) {
          contour.push( i );
          img->data[i] = CONTOUR_MASK;
        }
      }
    }
  }

  //printf( "contour size = %d\n", contour.size() );
  int iter = 0;
  unsigned long deleted = 0;
  //while( contour.size() > 0 ) {
  do {
    deleted = palagyi_fpta( img, contour, lut );


    iter++;
  } while( deleted > 0 );

/*
  for ( unsigned long i = 0; i < size; i++ ) {
    if ( (img->data[i] & OBJECT) == OBJECT ) img->data[i] = OBJECT;
  }
*/
  // return iter;

}

void Writeheader(FILE *fp, long &num){

  // write the header parameters to the top of the output file
  int check = 0;
  int size_l = sizeof(long);

  check += fwrite(&(inp_volumesize[OR_Z]), size_l, 1, fp);
  check += fwrite(&(inp_volumesize[OR_Y]), size_l, 1, fp);
  check += fwrite(&(inp_volumesize[OR_X]), size_l, 1, fp);
  check += fwrite(&(inp_blocksize[OR_Z]), size_l, 1, fp);
  check += fwrite(&(inp_blocksize[OR_Y]), size_l, 1, fp);
  check += fwrite(&(inp_blocksize[OR_X]), size_l, 1, fp);
  check += fwrite(&num, size_l, 1, fp);

  if (check != 7) { fprintf(stderr, "Failed to write file in Writeheader\n"); exit(-1); }
}
