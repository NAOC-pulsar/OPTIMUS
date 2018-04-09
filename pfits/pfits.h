//  Copyright (C) 2015,2016 George Hobbs
// This file is part of the pfits software package
//

/* pfits is free software: you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version. 
 * pfits is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details. 
 * You should have received a copy of the GNU General Public License 
 * along with pfits.  If not, see <http://www.gnu.org/licenses/>. 
*/

/* Definition files for the pfits software */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"

#define FILE_LEN 128
#define MAX_CHANNELS 8192

// Structure containing the PSRFITS header information
typedef struct headerStruct {
  int headerSet; // 1 = set, 0 = unset
  int fileType; // 1 = search mode, 2 = fold mode, 3 = cal
  
  // All files
  int nchan; // Number of channels
  float freq;   // Central frequency for observing band (MHz)
  float bw;     // Observation bandwidth - always positive (MHz)
  int   stt_imjd;   // Integer MJD of observation start
  float stt_smjd; // Number of seconds after imjd for start
  float stt_offs; // Start time offset
  float zeroOff;  // Zero offset
  float chanbw;   // Channel bandwidth - can be positive or negative
  float dm;       // Dispersion measure (cm-3 pc)
  float chanFreq[MAX_CHANNELS]; // Frequency channels (centre of each channel)
  
  char source[FILE_LEN]; // Observation source
  char ra[1024];
  char dec[1024];

  // Search mode files
  //
  int nbits; // Number of bits
  int nsamp; // Number of samples
  int nsub;  // Number of subintegrations
  int nsblk; // Number of samples per subintegration
  float tsamp; // Sample time (s)
  int npol;    // Number of polarisations
  
  // Fold mode files
  int nbin;

} headerStruct;


typedef struct dSetStruct {
  headerStruct *head;
  int headerMemorySet;

  // File information
  int fileSet; // 0 = not set, 1 = set
  char fileName[FILE_LEN]; // File name
  fitsfile *fp;
  int  fileOpen; // 0 = not used, 1 = open, 2 = closed  
} dSetStruct;

typedef struct calibrateStruct {
  float baseline_p0[MAX_CHANNELS];
  float baseline_p1[MAX_CHANNELS];
  float baseline_p2[MAX_CHANNELS];
  float baseline_p3[MAX_CHANNELS];
} calibrateStruct;

// Function definitions
int initialise(dSetStruct **dSet,int debug);
int deallocateMemory(dSetStruct **dSet,int debug);
void errorStop(char *err,dSetStruct *dSet,int debug);
void pfitsLoadHeader(dSetStruct *dSet,int debug);
void pfitsOpenFile(dSetStruct *dSet,int debug);
void setFilename(char *str,dSetStruct *dSet,int debug);

// Loader
void pfits_read1pol_float(float *out,int polNum,dSetStruct *dSet,float t1,float t2,int rangeType,long *nSamples,int *nTimeSamples,int *nFreqSamples,int debugFlag);
void pfits_bytesToFloats(int samplesperbyte,int n,unsigned char *cVals,float *out);
void eightBitsFloat(int eight_bit_number, float *results, int *index);
void fourBitsFloat(int eight_bit_number, float *results, int *index);
void twoBitsFloat(int eight_bit_number, float *results, int *index);
void oneBitFloat(int eight_bit_number, float *results, int *index);
void pfits_read1pol_zeroDM_float(float *out,int polNum,dSetStruct *dSet,float t1,float t2,int rangeType,long *nSamples,int *nTimeSamples,int *nFreqSamples,int debugFlag);

// Calibration
void convertStokes(float p1,float p2,float p3,float p4,float *stokesI,float *stokesQ,float *stokesU,float *stokesV);
void calibrateScalePols(calibrateStruct *cal,float *p1,float *p2,float *p3,float *p4,int n);
