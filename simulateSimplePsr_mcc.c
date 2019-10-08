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
#include "ran1.c"
#include "gasdev.c"
#include "time.h"
#include "bessi0.c"

int main(int argc,char *argv[])
{
  header *head;
  int i,j,k;
  float amp = 3;
  double *sum;
  float si = -2;
  float e = 2.718281828 ;
  
  float width; // smear width 
  float *profile;
  float *randarr;
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
  float minfref = 244.140625;      // the profile will be 0 in those channels whose fref < minfref (by mcc)
  /*float minfref = 0.;*/
  int ngulp = 300000; // the gulp size == 60 s for tsamp = 0.0002s
  int ncap;
  int *randidnum;
  float kappa;
  
  long pn=0; // Pulse number
  FILE *fout;
  char outName[128] = "pulsar.dat";
  char description[128];
  char format[128];
  int useParamFile=0;
  char paramFile[MAX_PARAM_FILES][1024];
  int phasenum, nperiods;
  time_t idum = -1 * time(NULL);
  
  
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
  nsamp = (int)ceil((head->t1 - head->t0)/head->tsamp);  //samples in all file length
  nperiods = (int)ceil((head->t1 - head->t0)/head->p0);
  dt_dm = (float *)malloc(sizeof(float)*head->nchan);
  chanfref = (float *)malloc(sizeof(float)*head->nchan);
  sum = (double *)malloc(sizeof(double)*head->nchan);
  chanbw = (head->f2-head->f1)/head->nchan;
  //profile = (float *)malloc(sizeof(float)*npsamp*head->nchan);
  /*profile = (float *)malloc(sizeof(float)*nsamp*head->nchan);*/
  profile = (float *)malloc(sizeof(float)*ngulp*head->nchan);
  randarr = (float *)malloc(sizeof(float)*nperiods);
  randidnum = (int *)malloc(sizeof(int)*ngulp);

  //fref = 0.5*(head->f2+head->f1);
  fref = 0.5*(head->f2+head->f2);

  //open file && write header
  fout = fopen(outName,"wb");

  //change
  //----------------------------------------------------------------------------------------------------
  /*simulateWriteHeader(head,fout);*/

  // calculate the time smear in each channel 
  for (i=0;i<head->nchan;i++)
  {
    chanfref[i] = (head->f1+i*chanbw);
    //chanfref[i] = (head->f2-i*chanbw);
    if (chanfref[i] <= minfref)
      {
        dt_dm[i]=0.;
      }
    else
      {
        dt_dm[i] = (4.15e-3*head->dm*(pow(fref/1000.0,si)-pow(chanfref[i]/1000.0,si)));
      }
     /*printf("dt_dm[i]: %g, dm: %g,fref: %g,f1: %g,chanbw: %g,si: %g\n",dt_dm[i],head->dm,fref,head->f1,chanbw,si);*/
   }
  
  printf("npsamp:%d nsamp:%d chanbw:%f\n",npsamp,nsamp,fabs(chanbw));


  // get random number
  idum = -1 * time(NULL); 

  for (i=0;i<nperiods;i++)
    {
      // float randnum = gasdev(&idum)+1.;
      float randnum = gasdev(&idum)+1.;
      randarr[i] = (randnum+2.12132)/2.12132;
      //2.12132 is 3 sigma
      randarr[i] = randnum > 0 ? randnum : 0;
      printf("%d randarr:%f\n",i,randarr[i]);
    }

  printf("***we passed here!!!\n");
  idum = (int) -1 * time(NULL); //take a new idum

  for (i=0;i<ngulp;i++)
    {
      randidnum[i] = (int) (ran1(&idum) * ngulp);
    }
  idum = (int) -1 * time(NULL); //take a new idum


  // Make a simple profile
  for (k=0;k<=(int)(nsamp/ngulp);k++)
  {
    ncap = nsamp - k*ngulp;
    if (ncap > ngulp) ncap = ngulp;

//#pragma omp parallel for default(shared) private(i,j) shared(profile, sum)
    for (j=0;j<head->nchan;j++)
      {
        sum[j] = 0.;
        //chanfref < minfref the profile =0
        if (chanfref[j] <= minfref)
          {

/*#pragma omp parallel for default(shared) private(i) shared(profile, sum)*/
            for (i=0;i<ncap;i++) //from t=0 ~ t=t1
              {
                if (head->setFlux==0)
                  profile[i*head->nchan+j]*=amp;
                else if (sum[j] == 0.)
                  profile[i*head->nchan+j] = 0.;
              }
          }
        else
          {
            /*width = 2*dt_dm[j]*fabs(chanbw)/(head->f1+(j+0.5)*chanbw) + head->width;chanfref[i]*/
            /*width = fabs(2*dt_dm[j]*chanbw/(chanfref[j])) + head->width;*/
            /*width = head->width;*/
            
            width = (fabs(2*dt_dm[j]*chanbw/(chanfref[j])) + head->width) * 2 * 3.14159265;
            if (j % 100 == 0 ) printf("j: %d width: %f ,dt_dm: %f ,fref: %f\n",j, width,dt_dm[j],chanfref[j]);

/*#pragma omp parallel for default(shared) private(i) shared(profile, sum, randidnum)*/
/*#pragma omp parallel for default(shared) private(i) shared(profile, sum)*/
            for (i=0;i<ncap;i++) //from t=0 ~ t=t1
              {
                /*printf("*i, j: %d %d\n", i, j);*/
                tval = (i + k*ngulp)*head->tsamp + dt_dm[j];
                // tval = (i + k*ngulp)*head->tsamp;
                phase = fmod(tval/head->p0, 1.);
                if (phase < 0.) phase += 1.;
                phase -= 0.5;

                phase *= 2*3.14159265;
                
                phasenum = (int)(tval/head->p0-fmod(tval/head->p0, 1.));
                /*profile[(i*head->nchan)+j] = exp(-pow(phase,2)/2.0/head->width/head->width);*/
                
                /*idx = &randidnum[i];*/

                /*float lognormal = exp(gasdev(&randidnum + i*sizeof(int))) / e;*/
                /*float lognormal = exp(gasdev(&idum)) / e;*/

                /*profile[(i*head->nchan)+j] = exp(-pow(phase,2)/2.0/width/width) * randarr[phasenum] * exp(gasdev(&idum)) / e;*/

                //change to use von Mises profile instead of Gaussian
                
                /*kappa = 1./width/width;*/
                /*if (kappa > 100.){*/
                /*profile[(i*head->nchan)+j] = exp(-pow(phase,2)/2.0/width/width)/width/sqrt(2.* 3.14159265) * randarr[phasenum]* exp(gasdev(&idum)) / e;*/
                /*}*/
                /*else {*/ 
                /*profile[(i*head->nchan)+j] = exp(kappa * cos(phase))/2./3.14159265/bessi0(kappa) * randarr[phasenum] * exp(gasdev(&idum)) / e;*/
                /*}*/

                kappa = 1./width/width;
                if (kappa > 80.){
                profile[(i*head->nchan)+j] = exp(-pow(phase,2)/2.0/width/width)/width/sqrt(2.* 3.14159265)* randarr[phasenum] * exp(gasdev(&idum)) / e ;
                /*profile[(i*head->nchan)+j] = exp(-pow(phase,2)/2.0/width/width);*/
                }
                else {
                profile[(i*head->nchan)+j] = exp(kappa * cos(phase))/2./3.14159265/bessi0(kappa)* randarr[phasenum] * exp(gasdev(&idum)) / e ;
                }

                
          /*if (j == 50){*/
            /*if (i<2*npsamp){*/

                /*printf("profile: %d %g %g %g %g\n ",i,profile[i*head->nchan+j],phase,tval, bessi0(kappa));*/
                      /*}*/
                        /*}*/


                /*profile[(i*head->nchan)+j] = exp(-pow(phase,2)/2.0/width/width) * randarr[phasenum] * exp(gasdev(&randidnum + i*sizeof(int))) / e;*/
                /*profile[(i*head->nchan)+j] = exp(-pow(phase,2)/2.0/width/width) * randarr[phasenum] * lognormal;*/
                /*profile[(i*head->nchan)+j] = exp(-pow(phase,2)/2.0/width/width) * randarr[phasenum] ;*/
                /*profile[(i*head->nchan)+j] = exp(-pow(phase,2)/2.0/width/width) ;*/
                /*profile[(i*head->nchan)+j] = exp(-1.*phase*phase/2.0/width/width) ;*/
                sum[j] += profile[(i*head->nchan)+j];
              }
            /*printf("i: %d ,sum[j]: %g\n",i,sum[j]);*/

  /*Make a simple profile*/
          /*if (j == 50){*/
            /*for (i=0;i<npsamp;i++){*/

                /*[>tval = (i + k*ngulp)*head->tsamp + dt_dm[j];<]*/
                /*tval = (i + k*ngulp)*head->tsamp;*/
                /*phase = fmod(tval/head->p0, 1.);*/
                /*if (phase < 0.) phase += 1.;*/
                /*phase -= 0.5;*/
                
                /*printf("profile: %d %g %g %g\n ",i,profile[i*head->nchan+j],phase,tval);*/
                      /*}*/
                        /*}*/
 

#pragma omp parallel for default(shared) private(i) shared(profile, sum)
             for (i=0;i<ncap;i++) //from t=0 ~ t=t1
               {
                 if (head->setFlux==0)
                   profile[i*head->nchan+j]*=amp;
                 else if (sum[j] > 0.)
                 {
                 /*printf("i,j, (%d, %d),  profile[i*head->nchan+j] %f,  sum[j] %f\n", i, j, profile[i*head->nchan+j], sum[j]);*/
                   //tval = (i + k*ngulp)*head->tsamp + dt_dm[j];
                   profile[i*head->nchan+j]*=(head->flux[j]/(sum[j]/(double)ncap));
                 /*printf("i,j, (%d, %d) head->flux: %f\n", i, j, head->flux[j]);*/
                 }
   }
          }
      }
    /*printf("k: %d ,profile[998]: %g\n; ### sizeof(float) %d",k,profile[1000],*/
            /*sizeof(float));*/

    fwrite(profile,sizeof(float),ncap*head->nchan,fout);

  }



  fclose(fout);
  simulateReleaseMemory(head);
  free(profile);
  free(sum);
  free(dt_dm);
  free(chanfref);
  free(head);
  free(randarr);
  free(randidnum);
}
