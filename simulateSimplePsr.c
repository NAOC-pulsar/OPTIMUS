// gcc -O3 -lm -o simulateSimplePsr simulateSimplePsr.c
// This routine simulates a simple pulsar that is defined by a single periodicity
// (i.e., this routine can not be used when modelling Doppler effects from the Earth's motion or
// from orbital effects).
// This routine also assumes that the dispersion smearing across the band is less than half the pulse period


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "T2toolkit.h"
#include <string.h>
#include "simulate.h"

int main(int argc,char *argv[])
{
  header *head;
  int i,j;
  float amp = 3;
  double *sum;
  float si = -2;
  
  float *profile;
  int npsamp;
  float tval,phase;
  float fval;
  int ival;
  float t;
  float *dt_dm;
  float fref;
  float chanbw;
  
  long pn=0; // Pulse number
  FILE *fout;
  char outName[128] = "pulsar.dat";
  char description[128];
  char format[128];
  int useParamFile=0;
  char paramFile[MAX_PARAM_FILES][1024];
  
  head = (header *)malloc(sizeof(header));
  simulateSetHeaderDefaults(head);

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
    }
  if (useParamFile>0)
    {
      for (i=0;i<useParamFile;i++)
	{
	  printf("Reading parameter file: %s\n",paramFile[i]);
	  simulateReadParamFile(head,paramFile[i]);
	}
    }
  npsamp = (int)(head->p0/head->tsamp);
  dt_dm = (float *)malloc(sizeof(float)*head->nchan);
  profile = (float *)malloc(sizeof(float)*npsamp*head->nchan);
  sum = (double *)malloc(sizeof(double)*head->nchan);
  chanbw = (head->f2-head->f1)/head->nchan;

  fref = 0.5*(head->f2+head->f1);

  for (i=0;i<head->nchan;i++)
    {
      dt_dm[i] = (4.15e-3*head->dm*(pow(fref/1000.0,si)-pow((head->f1+i*chanbw)/1000.0,si)));
      printf("dt_dm[i] = %g %g %g %g %g %g\n",dt_dm[i],head->dm,fref,head->f1,chanbw,si);
    }
  
  // Make a simple profile
  for (i=0;i<npsamp;i++)
    {
      tval = (i)*head->tsamp;
      phase =  tval/head->p0-(int)(tval/head->p0); //+phi0;

      for (j=0;j<head->nchan;j++)
	{
	  profile[i*head->nchan+j] = exp(-pow(phase-0.5+dt_dm[j]/head->p0,2)/2.0/head->width/head->width);
	  if (i==0) 
	    sum[j]=profile[i*head->nchan+j];
	  else
	    sum[j]+=profile[i*head->nchan+j];

	  if (j==0 || j == 50)
	    {
	      printf("profile: %d %g %g ",i,profile[i*head->nchan+j],phase);
	      if (j==50) printf("\n");
	    }
	}
    }
  // Scale the profile amplitude
  for (i=0;i<npsamp;i++)
    {
      for (j=0;j<head->nchan;j++)
	{
	  if (i==0 && head->setFlux==1)
	    printf("Amplitude scaling: %g\n",head->flux[j]/(sum[j]/(double)npsamp));
	  if (head->setFlux==0)
	    profile[i*head->nchan+j]*=amp;
	  else
	    profile[i*head->nchan+j]*=(head->flux[j]/(sum[j]/(double)npsamp));
	}
    }

  fout = fopen(outName,"wb");
  simulateWriteHeader(head,fout);
  
  t=head->t0;
  do {
    // Put down the first pulse
    fwrite(profile,sizeof(float),npsamp*head->nchan,fout);
    t+=head->p0;
    pn++;
  } while (t < head->t1);

  fclose(fout);
  simulateReleaseMemory(head);
  free(profile);
  free(sum);
  free(dt_dm);
  free(head);
}
