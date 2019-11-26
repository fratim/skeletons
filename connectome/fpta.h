#ifndef __PALAGYI_FPTA_H__
#define __PALAGYI_FPTA_H__

#include "pgm.h"
#include <queue>

using namespace std;

bool generateLUT_FP(  );
bool generateLUT_SI_NE( );
bool generateLUT_SI_SW( );
bool generateLUT_SI_N(  );

bool generateLUT_SF(  );


unsigned long palagyi_fpta( PGMImage* img, queue<unsigned long> &contour, unsigned char* lut );
unsigned long palagyi_SI4( PGMImage* img, queue<unsigned long> &contour, unsigned char* lut );
int fpta_thinning( PGMImage* img, unsigned char *lut );


#endif
