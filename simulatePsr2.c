#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "T2toolkit.h"
#include <string.h>

int main(int argc,char *argv[])
{
  int i,j;
  
  float tsamp = 250e-6;
  float si = -2;
  int nchan=96;
  float f1 = 1518;
  float f2 = 1230;
  float t0 = -1050;
  float t1 = 1050;
  float p0 = 0.3;
  float w = p0/20.0;
  float dm = 50;
  
  float *profile;
  int npsamp = (int)(p0/tsamp);
  float tval,phase;
  float fval;
  int ival;
  float t;
  float dt_dm[nchan];
  float fref;
  float chanbw=(f2-f1)/nchan;
  
  long pn=0; // Pulse number
  FILE *fout;
  char description[128];
  
  
  profile = (float *)malloc(sizeof(float)*npsamp*nchan);
  fref = 0.5*(f2+f1);

  for (i=0;i<nchan;i++)
    {
      dt_dm[i] = (4.15e-3*dm*(pow(fref/1000.0,si)-pow((f1+i*chanbw)/1000.0,si)));
      printf("dt_dm[i] = %g %g %g %g %g %g\n",dt_dm[i],dm,fref,f1,chanbw,si);
    }
  
  // Make a simple profile
  for (i=0;i<npsamp;i++)
    {
      tval = (i)*tsamp;
      phase =  tval/p0-(int)(tval/p0); //+phi0;

      for (j=0;j<nchan;j++)
	{
	  profile[i*nchan+j] = exp(-pow(phase-0.5+dt_dm[j]/p0,2)/2.0/w/w);
	  if (j==0 || j == 50)
	    {
	      printf("profile: %d %g %g ",i,profile[i*nchan+j],phase);
	      if (j==50) printf("\n");
	    }
	}
    }

  fout = fopen("pulsar.dat","wb");
  strcpy(description,"Pulsar");
  fwrite(description,1,128,fout);
  fval = t0; fwrite(&fval,sizeof(float),1,fout);
  fval = t1; fwrite(&fval,sizeof(float),1,fout);
  fval = tsamp; fwrite(&fval,sizeof(float),1,fout);
  fval = f1; fwrite(&fval,sizeof(float),1,fout);
  fval = f2; fwrite(&fval,sizeof(float),1,fout);
  ival = nchan; fwrite(&ival,sizeof(int),1,fout);
  fval = 0.0000; fwrite(&fval,sizeof(float),1,fout); // Right ascension
  fval = 0.0000; fwrite(&fval,sizeof(float),1,fout); // Declination

  t=t0;
  do {
    // Put down the first pulse
    fwrite(profile,sizeof(float),npsamp*nchan,fout);
    t+=p0;
    pn++;
  } while (t < t1);

  fclose(fout);
  free(profile);
}
