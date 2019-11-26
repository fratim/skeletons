#include "pgm.h"
#include <stdio.h>
#include <malloc.h>
#include <iostream>

// #ifndef _PGM_
// #define _PGM_


bool readPGM( char* filename, PGMImage* image );
bool createPGM(PGMImage* image );
bool writePGM( char* filename, PGMImage* image );
PGMImage::PGMImage();
PGMImage::PGMImage( int width, int height, int depth );
PGMImage::PGMImage( PGMImage &image );
PGMImage::~PGMImage();

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

// #endif
