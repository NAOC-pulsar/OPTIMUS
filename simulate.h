#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define MAX_PARAM_FILES 5 // Maximum number of parameter files
#define MAX_RFI 100000    // Maximum number of RFI events
#define MAX_NB_RFI 50     // Maximum number of narrow band RFI signals

typedef struct narowBandRFI {
  float f1; // Low frequency of RFI signal (MHz)
  float f2; // High frequency of RFI signal (MHz)
  float amp; // Amplitude (Jy)
  float t0;  // Start time (seconds since start)
  float duration; // Duration of RFI event (s)
}narrowBandRFI;


// Define the header information
typedef struct header {
  char  format[64]; // Header format
  char  name[128]; // Name of this source
  float f1;         // Frequency of first channel in MHz
  float f2;         // Frequency of last channel in MHz
  int   nchan;      // Number of frequency channels
  float tsamp;      // Sample time in seconds
  float t0;         // Start time in seconds from nominal LST time
  float t1;         // End time in seconds from nominal LST time
  int   imjd;       // Integer MJD
  int   smjd;       // Seconds since imjd
  float stt_offs;   // Offset from imjd + smjd/86400 (sec)
  float raj_rad;    // Right ascension (in radians)
  float decj_rad;   // Declination (in radians)
  int   useAngle;   // 1 = scale with angle to source, 0 = do not scale
  long  initialSeed; // Random number seed
  long  seed;       // Random number seed being used


  
  // System parameters
  float *gain;      // Telescope gain in K/Jy for each frequency channel
  float *tsys;      // System temperature in K for each frequency channel
  int   setGain;    // Has the gain been set? (1=yes, 0=no)
  int   setTsys;    // Has the Tsys been set? (1=yes, 0=no)
  float diameter;   // Telescope diameter in m
  int   surveyType; // 1 = pointed, 2 = drift
  
  // Digitisation
  int   nbits;      // Number of bits
  int   levelSet;   // 0 = no level setting, 1 = set to single value
  // Pulsars
  char predictor[1024]; // Predictor file name
  float p0;
  float dm;
  float width;
  float *flux;      // Flux density in each frequency channel (Jy)
  int   setFlux;
  int   addScatter;
  int   addDMsmear;
  
  // Beam parameters
  float beamRA0;
  float beamDEC0;
  int   setBeam;

  // RFI
  char rfiFile[1024];
  int  nRFI;
  narrowBandRFI *nb_rfi;
  int  n_nb_rfi;

  // Calibration source
  float calAmp;
  float calFreq;
  float calDuty;
} header;

// Function definitions
void simulateWriteHeader(header *head,FILE *fout);
void simulateSetHeaderDefaults(header *head);
void simulateReadHeaderParameters(header *head,FILE *fin);
void simulateDisplayHeaderParameters(header *head);
void simulateReadParamFile(header *head,char *paramFile);
void simulateReleaseMemory(header *head);
