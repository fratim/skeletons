// #ifndef PGM_H
// #define PGM_H

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



bool readPGM( char* filename, PGMImage* image );

bool writePGM( char* filename, PGMImage* image );

bool createPGM(PGMImage* image );



// #endif
