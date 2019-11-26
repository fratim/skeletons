#include <stdio.h>
#include "pgm.h"
// #include "pointenv.h"
#include "powersof2.h"
#include <malloc.h>
#include <vector>
#include <queue>
#include <string.h>
#include <time.h>
#include <iostream>
//#include "fpta.h"

#define CONTOUR_MASK 3
#define OBJECT 1

using namespace std;

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

void help() {
  printf( "Parallel SHRINKING Algorithm (no endpoints preserved)\n" );
  printf( "Usage:  ronse_fp-<ends> input.pgm output.pgm\n" );
  printf( "\n");
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

int main( int argc, char** argv ) {

  time_t start_time, end_time, start_init, end_init, start_thinning, end_thinning;


  start_time = clock();
  if ( argc < 3 )  {
    help();
    return 0;
  }

  FILE* fp = NULL;

	std::cout << "opnening LUTB6" << std::endl << std::flush;
  fp = fopen( "ronse_fpta.lut", "rb" );

  unsigned char* lut = NULL;
  lut = (unsigned char*) malloc ( 16777216 * sizeof(unsigned char) );
  if ( !lut ) {
    printf( "ERROR: Can not allocate memory for LUT\n" );
    fflush(fp);
    fclose(fp);
    return 0;
  }

  //fread( lut, sizeof(unsigned char), 16777216, fp );
  unsigned long i = 0;
  while( !feof(fp) ) {
    unsigned char bits;
    fread( &bits, sizeof(unsigned char), 1, fp );
    for ( unsigned int j = 0, k=1; j < 8; j++, k*=2, i++ ) {
      if ( (bits & k) > 0 ) lut[i] = 1; else lut[i] = 0;
    }
  }

	int width = 5;
	int height = 5;
	int depth = 1;

	PGMImage *img = new PGMImage(width,height,depth);

  createPGM( img );
  printf( "%s\n", argv[2] );
  printf( "width: %d, height = %d\n", img->width, img->height );

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

  end_init = clock();

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


  writePGM(argv[2], img );

  delete img;


  free( lut );

  fflush(fp);
  fclose(fp);


  end_time = clock();
  printf( "\n%s :: ", argv[2] );
  printf( "Number of iterations: %d ", iter );
  printf( "#object: %d ", nobject );
  printf( "#skeletal: %d ", nskel );
  printf( "Init. time: %lf sec ", (double) ( end_init - start_init ) / CLOCKS_PER_SEC );
  printf( "Thinning time: %lf sec ", (double) ( end_thinning - start_thinning ) / CLOCKS_PER_SEC );
  printf( "Elapsed time: %lf sec\n", (double) ( end_time - start_time ) / CLOCKS_PER_SEC );

  printf( "\n" );

  return 0;




}
