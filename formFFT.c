#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <fftw3.h>

// gcc -lm -o formFFT formFFT.c -L../fftw-3.3.4//.libs/ -I../fftw-3.3.4/api/ -lfftw3f

int main(int argc,char *argv[])
{
  char fname[1024];
  char outName[1024];
  FILE *fin;
  long maxVals = 600000;
  float *vals;
  fftwf_complex *outFFT;
  long i;
  long npts=0;
  float dummy;
  fftwf_plan p;
  FILE *fout;
   
  if (!(vals = (float *)malloc(sizeof(float)*maxVals)))
    {
      printf("Unable to allocate memory\n");
      exit(1);
    }
  strcpy(fname,argv[1]);
  strcpy(outName,argv[2]);
  fin = fopen(fname,"r");
  while (!feof(fin))
    {
      if (fscanf(fin,"%f %f",&dummy,&vals[npts])==2)
	{
	  //	  if (dummy > 100 && dummy < 120)
	    npts++;
	}
    }
  npts = 524288;
  fclose(fin);
  printf("Loaded %ld points\n",npts);
  outFFT = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*npts);
  printf("Doing the FFT\n");
  p = fftwf_plan_dft_r2c_1d(npts,vals,outFFT,FFTW_ESTIMATE);  
  fftwf_execute(p);
  fftwf_destroy_plan(p);
  printf("Done the FFT\n");
  fout = fopen(outName,"w");
  for (i=0;i<npts;i++)
    {
      fprintf(fout,"%ld %g\n",i,pow(outFFT[i][0],2)+pow(outFFT[i][1],2));
    }
  free(vals);
  fftwf_free(outFFT);
   

}
