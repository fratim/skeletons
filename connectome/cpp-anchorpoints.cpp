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
#include "powersof2.h"
#include <malloc.h>
#include <queue>
#include <string.h>
#include <time.h>

#define CONTOUR_MASK 3
#define OBJECT 1

using namespace std;

class PGMImage {
 public:
   unsigned char* data;
   int width, height, depth;
   PGMImage();
   PGMImage( int width, int height, int depth );
   PGMImage( PGMImage &image );
   ~PGMImage();
   void normalize();
   void binarize(unsigned char threshold = 128);
   int createBorder( int bwidth, unsigned char color = 0 );
   void separateHeader();
};

bool readPGM( char* filename, PGMImage* image ) {
 FILE *fp = NULL;
 char p5[4];
 int width, height, depth;
 char c[10];

 fp = fopen( filename, "rb"  );
 if ( !fp ) {
   printf ( "ERROR: %s can not read\n", filename );
   return false;
 }

 fscanf ( fp, "%s", &p5 );
 fscanf ( fp, "%d", &(image->width) );
 fscanf ( fp, "%d", &(image->height) );
 fscanf ( fp, "%d\n", &(image->depth) );
 //printf ( "Width: %d, Height: %d\n", image->width, image->height );

 if ( image->data != NULL ) {
   //delete[] image->data;
   free ( image->data );
   image->data = NULL;
 }
 image->data = (unsigned char*) malloc ( image->width * image->height * sizeof( unsigned char ));
 unsigned long size =  image->width * image->height;
 if ( image->data != NULL ) {
   for ( unsigned long i = 0; i < size; i++ ) image->data[ i ] = 0;
/*
   for ( unsigned long i = 0; i < size; i++ ) {
       image->data[ i ] = (unsigned char) fgetc( fp );
   }
*/
   fread( (unsigned char* ) image->data, sizeof( unsigned char ), image->width * image->height, fp );
 } else {
     printf ("ERROR: readPGM  -- Not enougth memory\n ");
     return false;
 }

 //printf ( "End of reading\n" );

 fflush(fp);
 fclose(fp);
 return true;
}

bool createPGM(PGMImage* image ) {
 if ( image->data != NULL ) {
   //delete[] image->data;
   free ( image->data );
   image->data = NULL;
 }
 image->data = (unsigned char*) malloc ( image->width * image->height * sizeof( unsigned char ));
 unsigned long size =  image->width * image->height;
 std::cout << "creating dummy pgi" << std::endl <<std::flush;

 if ( image->data != NULL ) {
   image->data[0]=0;
   image->data[1]=1;
   image->data[2]=1;
   image->data[3]=0;
   image->data[4]=0;

   image->data[5]=0;
   image->data[6]=1;
   image->data[7]=1;
   image->data[8]=1;
   image->data[9]=0;

   image->data[10]=0;
   image->data[11]=0;
   image->data[12]=0;
   image->data[13]=0;
   image->data[14]=0;

   image->data[15]=3;
   image->data[16]=0;
   image->data[17]=2;
   image->data[18]=2;
   image->data[19]=0;

   for ( unsigned long i = 0; i < size; i++ ) {
     std::cout << (int)(image->data[i]);
     if ((i+1)%(image->width)==0) std::cout << std::endl;
   }
 } else {
     printf ("ERROR: readPGM  -- Not enougth memory\n ");
     return false;
 }

 printf ( "End of reading\n" );
 return true;
}

bool writePGM( char* filename, PGMImage* image ) {
 FILE *fp = NULL;


 fp = fopen( filename, "wb"  );
 if ( !fp ) {
   printf ( "ERROR: %s can not write\n", filename );
   return false;
 }

 fprintf ( fp, "P5\n%d %d\n%d\n", image->width, image->height, image->depth );

 unsigned long size =  image->width * image->height;
 std::cout << "creating dummy pgi" << std::endl <<std::flush;

 if ( image->data != NULL ) {
   for ( unsigned long i = 0; i < size; i++ ) {
     std::cout << (int)(image->data[i]);
     if ((i+1)%(image->width)==0) std::cout << std::endl;
   }
 }

 if ( image->data != NULL ) {
   fwrite( (unsigned char*) image->data, sizeof( unsigned char ), image->width * image->height, fp );
/*
   unsigned long size = image->width * image->height;

   for ( unsigned long i = 0; i < size; i++ ) {
       fputc( image->data[ i ], fp );
   }
*/

 } else {
     return false;
 }
 fflush(fp);
 fclose(fp);
 return true;
}

PGMImage::PGMImage() {
   this->data = NULL;
   this->width = 0;
   this->height = 0;
   this->depth = 0;
}

PGMImage::PGMImage( int width, int height, int depth ) {
 this->width = width;
 this->height = height;
 this->depth = depth;

 this->data = NULL;
 this->data = ( unsigned char* ) malloc ( width * height * sizeof( unsigned char ) );
 if ( this->data == NULL ) {
   printf ("PGMImage: Allocation error. Not enough memory. \n");
   return;
 }
 int size = width * height;
 for ( int i = 0; i < size; i++ ) {
   this->data[ i ] = 0;
 }


}

PGMImage::PGMImage( PGMImage &image ) {
 this->width = image.width;
 this->height = image.height;
 this->depth = image.depth;

 this->data = NULL;
 this->data = ( unsigned char* ) malloc ( width * height * sizeof( unsigned char ) );
 if ( this->data == NULL ) {
   printf ("PGMImage: Allocation error. Not enough memory. \n");
   return;
 }
 int size = width * height;
 for ( int i = 0; i < size; i++ ) {
   this->data[ i ] = image.data[ i ];
 }

}

PGMImage::~PGMImage() {
   delete[] this->data;
}

void PGMImage::normalize() {
 int size = this->width * this->height;
 int max = 0;
 for ( int i = 0; i < size; i++ ) {
   if ( this->data[ i ] > max ) max = this->data[i];
 }
 if ( max != 0 ) {
   for ( int i = 0; i < size; i++ ) {
     this->data[i] = this->data[i] * 255 / max;
   }
 }
 //this->depth = max;
}

void PGMImage::binarize(unsigned char threshold) {
 int size = this->width * this->height;
 //this->normalize();
 for ( unsigned long i = 0; i < size; i++ ) {
   if ( this->data[i] > threshold ) {
     this->data[i] = 255;
   } else {
     this->data[i] = 0;
   }
 }
}

int PGMImage::createBorder( int bwidth, unsigned char color ) {
 if ( this->data == NULL ) {
   printf( "ERROR: PGMImage::createBorder: data does not exist\n " );
   return -1;
 }
 int new_width, new_height;
 new_width = this->width + 2 * bwidth;
 new_height = this->height + 2 * bwidth;
 unsigned long new_length = new_width * new_height;
 unsigned char *temp;
 temp = ( unsigned char* ) malloc ( new_length * sizeof( unsigned char ) );
 if ( !temp ) {
   printf( "ERROR: PGMImage::createBorder: not enough memory for temp array\n " );
   return -1;
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

/*
 if ( bwidth > 0 ) {
   for ( int y = 0; y < height; y++ ) {
     for ( int x = 0; x < width; x++ ) {
       orig_index = y * width + x;
       // new_index = ( y + bwidth ) * (width + 2*bwidth ) + (x + bwidth);
       int new_x, new_y;
       new_x = x + bwidth;
       new_y = y + bwidth;
       new_index = new_y * new_width + new_x;
       temp[ new_index ] = this->data[orig_index];
     }
   }
 } else {
   for ( int y = 0; y < new_height; y++ ) {
     for ( int x = 0; x < new_width; x++ ) {
       int old_x = x + bwidth;
       int old_y = y + bwidth;
       orig_index = old_y * width + old_x;
       new_index = y * new_width + x;
       temp[ new_index ] = this->data[orig_index];
     }
   }
 }
*/

 free ( data );
 data = temp;
 this->width += 2*bwidth;
 this->height += 2 * bwidth;

}

unsigned long palagyi_fpta( PGMImage* img, queue<unsigned long> &contour, unsigned char* lut );
int fpta_thinning( PGMImage* img, unsigned char *lut, unsigned char *lut2 );
void ThinImages(PGMImage* img);

void ComputeAnchorPoints(const char *prefix, const char* output_dir, long input_blocksize_inp[3], long *z_min_wall, long *z_max_wall, long *y_min_wall, long *y_max_wall, long *x_min_wall, long *x_max_wall){

  long input_blocksize[3];
  input_blocksize[0] = input_blocksize_inp[0];
  input_blocksize[1] = input_blocksize_inp[1];
  input_blocksize[2] = input_blocksize_inp[2];

  long entries = (input_blocksize[OR_Y]*input_blocksize[OR_X]);
  std::unordered_set<long> IDs_present = std::unordered_set<long>();
  std::unordered_map <long, PGMImage*> images = std::unordered_map <long, PGMImage*>();

  long current_ID;
  int width = (int)input_blocksize[OR_X];
  int height = (int)input_blocksize[OR_Y];
  int depth = 1;

  for (long pos =0; pos<entries; pos++){
    if (z_max_wall[pos]==z_min_wall[pos]) current_ID = z_max_wall[pos];
    else current_ID = 0;

    if (current_ID!=0){
      if (IDs_present.find(current_ID)==IDs_present.end()){
        IDs_present.insert(current_ID);
        images[current_ID] = new PGMImage(width,height,depth);
        images[current_ID]->data[pos]=current_ID;
      }
      else{
        images[current_ID]->data[pos]=current_ID;
      }
    }
  }

  for (std::unordered_map <long, PGMImage*>::iterator iter = images.begin(); iter != images.end(); iter++){
    ThinImages(iter->second);
  }

  std::cout << "images found: " << images.size() <<std::endl;

}

void ThinImages(PGMImage* img){

  time_t start_time, end_time, start_thinning, end_thinning;

  start_time = clock();

  FILE* fp = NULL;

	std::cout << "opening LUT" << std::endl << std::flush;
  fp = fopen( "/home/frtim/Documents/Code/skeletons/connectome/ronse_fpta.lut", "rb" );

  unsigned char* lut = NULL;
  lut = (unsigned char*) malloc ( 16777216 * sizeof(unsigned char) );
  if ( !lut ) {
    printf( "ERROR: Can not allocate memory for LUT\n" );
    fflush(fp);
    fclose(fp);
  }

  printf( "width: %d, height = %d\n", img->width, img->height );
  std::cout << "HOSSA A" << std::endl << std::flush;

  //fread( lut, sizeof(unsigned char), 16777216, fp );
  unsigned long i = 0;
  while( !feof(fp) ) {
    unsigned char bits;
    fread( &bits, sizeof(unsigned char), 1, fp );
    for ( unsigned int j = 0, k=1; j < 8; j++, k*=2, i++ ) {
      if ( (bits & k) > 0 ) lut[i] = 1; else lut[i] = 0;
    }
  }

  std::cout << "HOSSA B" << std::endl << std::flush;



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

  start_thinning = clock();

  int iter = fpta_thinning( img, lut, NULL );

  end_thinning = clock();

  img->createBorder(-3);

  unsigned long nskel = 0;
  index = 0;
  for ( int y = 0; y < img->height; y++ ) {
    for ( int x = 0; x < img->width; x++ ) {
      //if ( ((img->data[ index ] & CONTOUR_MASK) != 0) || ((img->data[index] & OBJECT)) != 0 ) {
      if ( (img->data[index] & OBJECT) == OBJECT) {
        nskel++;
        img->data[index] = 255;
      } else {
        img->data[index] = 0;
      }
      index++;
    }
  }

  //
  // writePGM(argv[2], img );

  delete img;


  free( lut );

  fflush(fp);
  fclose(fp);


  end_time = clock();
  // printf( "\n%s :: ", argv[2] );
  printf( "Number of iterations: %d ", iter );
  printf( "#object: %d ", nobject );
  printf( "#skeletal: %d ", nskel );
  printf( "Thinning time: %lf sec ", (double) ( end_thinning - start_thinning ) / CLOCKS_PER_SEC );
  printf( "Elapsed time: %lf sec\n", (double) ( end_time - start_time ) / CLOCKS_PER_SEC );

  printf( "\n" );
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


int fpta_thinning( PGMImage* img, unsigned char *lut, unsigned char *lut2 ) {
  queue<unsigned long> contour;
  int w = img->width, h = img->height;
  unsigned long size = w * h;
  int n4[4] = { -w, 1, w, -1 };

  printf( "Thinning started...\n" );
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
  return iter;

}
