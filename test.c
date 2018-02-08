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
  double *sum1;
  float si = -2;
  
  float *profile;
  float *profile1;
  int npsamp;
  int nsamp;
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
  sum = (double *)malloc(sizeof(double)*head->nchan);
  sum1 = (double *)malloc(sizeof(double)*head->nchan);
  chanbw = (head->f2-head->f1)/head->nchan;
  nsamp = (int)(head->t1 - head->t0)/head->tsamp;
  profile = (float *)malloc(sizeof(float)*nsamp*head->nchan);
  profile1 = (float *)malloc(sizeof(float)*npsamp*head->nchan);

  fref = 0.5*(head->f2+head->f1);

  for (i=0;i<head->nchan;i++)
    {
      dt_dm[i] = (4.15e-3*head->dm*(pow(fref/1000.0,si)-pow((head->f1+i*chanbw)/1000.0,si)));
      printf("dt_dm[i]: %g, dm: %g,fref: %g,f1: %g,chanbw: %g,si: %g\n",dt_dm[i],head->dm,fref,head->f1,chanbw,si);
    }
  
  // Make a simple profile
  for (j=0;j<head->nchan;j++)
    {
     sum[j] = 0.;

  //for (i=0;i<npsamp;i++)
      for (i=0;i<nsamp;i++) //from t=0 ~ t=t1
	{
	   tval = (i)*head->tsamp + dt_dm[j];
	   //phase =  tval/head->p0-(int)(tval/head->p0); //+phi0;
	   phase =  fmod(tval/head->p0, 1.);
           if (phase < 0.) phase += 1.;
           phase -= 0.5;

	   profile[i*head->nchan+j] = exp(-pow(phase,2)/2.0/head->width/head->width);
	   sum[j]+=profile[i*head->nchan+j];

    }
  }
  

  // Make a simple profile
  for (i=0;i<npsamp;i++)
    {
      tval = (i)*head->tsamp;
      phase =  tval/head->p0-(int)(tval/head->p0); //+phi0;

      for (j=0;j<head->nchan;j++)
	{
	  profile1[i*head->nchan+j] = exp(-pow(phase-0.5+dt_dm[j]/head->p0,2)/2.0/head->width/head->width);
	  if (i==0) 
	    sum1[j]=profile[i*head->nchan+j];
	  else
	    sum1[j]+=profile[i*head->nchan+j];

	  if (j==0 || j == 50)
	    {
	      printf("profile: %d %g %g ",i,profile[i*head->nchan+j],phase);
	      if (j==50) printf("\n");
	    }
	}
    }

  
  for (j=0;j<head->nchan;j++)
    {
       printf("sum: %g ,nsamp: %d sum1: %g, npsamp: %d\n",sum[j],nsamp,sum1[j],npsamp);
    }

  j=1;
  for (i=0;i<npsamp;i++)
    {
       printf("profile: %g profile1: %g\n",profile[i*head->nchan+j],profile1[i*head->nchan+j]);
    }

	  
  free(profile);
  free(profile1);
  free(sum);
  free(dt_dm);
  free(head);
}
