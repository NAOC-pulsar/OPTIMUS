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
#include <omp.h>

int main(int argc,char *argv[])
{
  header *head;
  int i,j,k;
  float amp = 3;
  double *sum;
  float si = -2;
  
  float width; // smear width 
  float *profile;
  int npsamp;
  int nsamp;
  float tval,phase;
  float fval;
  int ival;
  float t;
  float *dt_dm;
  float fref;
  float chanbw;
  float *chanfref;
  int ngulp = 300000; // the gulp size == 60 s for tsamp = 0.0002s
  int ncap;
  
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
  //
  npsamp = (int)(head->p0/head->tsamp);     //samples in a period
  nsamp = (int)(head->t1 - head->t0)/head->tsamp;  //samples in all file length
  dt_dm = (float *)malloc(sizeof(float)*head->nchan);
  chanfref = (float *)malloc(sizeof(float)*head->nchan);
  sum = (double *)malloc(sizeof(double)*head->nchan);
  chanbw = (head->f2-head->f1)/head->nchan;
  //profile = (float *)malloc(sizeof(float)*npsamp*head->nchan);
  /*profile = (float *)malloc(sizeof(float)*nsamp*head->nchan);*/
  profile = (float *)malloc(sizeof(float)*ngulp*head->nchan);

  //fref = 0.5*(head->f2+head->f1);
  fref = 0.5*(head->f2+head->f2);
  //open file && write header
  fout = fopen(outName,"wb");
  simulateWriteHeader(head,fout);

  // calculate the time smear in each channel 
  for (i=0;i<head->nchan;i++)
  {
      chanfref[i] = (head->f1+i*chanbw);
      //chanfref[i] = (head->f2-i*chanbw);
      dt_dm[i] = (4.15e-3*head->dm*(pow(fref/1000.0,si)-pow(chanfref[i]/1000.0,si)));
      /*printf("dt_dm[i]: %g, dm: %g,fref: %g,f1: %g,chanbw: %g,si: %g\n",dt_dm[i],head->dm,fref,head->f1,chanbw,si);*/
    }
  
  printf("npsamp:%d nsamp:%d chanbw:%f\n",npsamp,nsamp,fabs(chanbw));
  

  // Make a simple profile
  for (k=0;k<=(int)(nsamp/ngulp);k++)
  {
      ncap = nsamp - k*ngulp;
      if (ncap > ngulp) ncap = ngulp;
  
//#pragma omp parallel for default(shared) private(i,j) shared(profile, sum)
      for (j=0;j<head->nchan;j++)
        {
          sum[j] = 0.;

          //width = 2*dt_dm[j]*fabs(chanbw)/(head->f1+(j+0.5)*chanbw) + head->width;chanfref[i]
          width = fabs(2*dt_dm[j]*chanbw/(chanfref[j])) + head->width;
          if (j % 100 == 0 ) printf("j: %d width: %f ,dt_dm: %f ,fref: %f\n",j, width,dt_dm[j],chanfref[j]);

#pragma omp parallel for default(shared) private(i) shared(profile, sum)
          for (i=0;i<ncap;i++) //from t=0 ~ t=t1
            {
               /*printf("*i, j: %d %d\n", i, j);*/
               tval = (i + k*ngulp)*head->tsamp + dt_dm[j];
               //tval = (i + k*ngulp)*head->tsamp;
               phase =  fmod(tval/head->p0, 1.);
               if (phase < 0.) phase += 1.;
               phase -= 0.5;
               
               //profile[(i*head->nchan)+j] = exp(-pow(phase,2)/2.0/head->width/head->width);
               profile[(i*head->nchan)+j] = exp(-pow(phase,2)/2.0/width/width);
               sum[j] += profile[(i*head->nchan)+j];
            }
          /*printf("i: %d ,sum[j]: %g\n",i,sum[j]);*/

#pragma omp parallel for default(shared) private(i) shared(profile, sum)
            for (i=0;i<ncap;i++) //from t=0 ~ t=t1
            {
              if (head->setFlux==0)
                profile[i*head->nchan+j]*=amp;
              else if (sum[j] > 0.)
                profile[i*head->nchan+j]*=(head->flux[j]/(sum[j]/(double)ncap));
            }
        }
        printf("k: %d ,profile[0]: %g\n",k,profile[0]);
    fwrite(profile,sizeof(float),ncap*head->nchan,fout);

  }


  fclose(fout);
  simulateReleaseMemory(head);
  free(profile);
  free(sum);
  free(dt_dm);
  free(head);
}
