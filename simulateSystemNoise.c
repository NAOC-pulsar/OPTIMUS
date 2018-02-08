// Routine to simulate system noise
// This is assumed to be constant during the observation
// The noise is calculated as a random number chosen from a
// Gaussian distribution.
//
// Currently this is very simple code that simulated a perfectly rectangular pulse
// gcc -O3 -lm -o simulateSystemNoise simulateSystemNoise.c simulate.c T2toolkit.c

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "T2toolkit.h"
#include "simulate.h"

int main(int argc, char *argv[])
{
  header *head;

  double t;
  float val;
  float *scale;
  
  char outName[1024];
  int i;
  FILE *fout,*fout2;
  int debugOut=0;
  double debugOutT1;
  double debugOutT2;
  int debugOutF1;
  int debugOutF2;
  char debugOutName[1024];

  char paramFile[MAX_PARAM_FILES][1024];
  int  useParamFile=0;
  double rad;
  float chanbw;
  int np=2;
  
  head = (header *)malloc(sizeof(header));
  
  simulateSetHeaderDefaults(head);
  strcpy(head->name,"Parkes system");

  // Read input parameters
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-o")==0)
	strcpy(outName,argv[++i]);
      else if (strcmp(argv[i],"-p")==0)
	{
	  strcpy(paramFile[useParamFile],argv[++i]);
	  useParamFile++;
	}
      else if (strcmp(argv[i],"-to")==0) // text output
	{
	  debugOut=1;
	  strcpy(debugOutName,argv[++i]);
	  sscanf(argv[++i],"%lf",&debugOutT1);
	  sscanf(argv[++i],"%lf",&debugOutT2);
	  sscanf(argv[++i],"%d",&debugOutF1);
	  sscanf(argv[++i],"%d",&debugOutF2);
	}

    }

  if (useParamFile>0)
    {
      for (i=0;i<useParamFile;i++)
	simulateReadParamFile(head,paramFile[i]);
    }
  chanbw = fabs(head->f2-head->f1)/head->nchan;
  scale = (float *)malloc(sizeof(float)*head->nchan);
  for (i=0;i<head->nchan;i++)
    {
      // Shall we calculate the scaling factors based on the radiometer equation?
      if (head->setGain == 1 && head->setTsys == 1)
	{
	  rad = head->tsys[i]/head->gain[i]/sqrt(head->tsamp*(chanbw*1e6)*np);
	  scale[i] = rad;
	  printf("Scaling is %g %g %g %g %g\n",scale[i],head->tsys[i],head->gain[i],head->tsamp,chanbw);
	}
      else
	scale[i] = 1;
    }
  
  fout = fopen(outName,"wb");
  if (debugOut==1)
    fout2 = fopen(debugOutName,"w");
  //
  // Write the header information
  simulateWriteHeader(head,fout);
  
  for (t=head->t0;t<head->t1;t+=head->tsamp)
    {
      if (debugOut == 1 && t >= debugOutT1 && t < debugOutT2)
	fprintf(fout2,"%f ",t);
      for (i=0;i<head->nchan;i++)
	{
	  val = scale[i]*TKgaussDev(&(head->seed));
	  fwrite(&val,sizeof(float),1,fout);
	  if (debugOut == 1 && t >= debugOutT1 && t < debugOutT2 &&
	      i >= debugOutF1 && i < debugOutF2)
	    fprintf(fout2,"%f ",val);
	}
      if (debugOut == 1 && t >= debugOutT1 && t < debugOutT2)
	fprintf(fout2,"\n");
    }
  fclose(fout);
  if (debugOut == 1)
    fclose(fout2); 

  simulateReleaseMemory(head);
  free(head);
  free(scale);
}
