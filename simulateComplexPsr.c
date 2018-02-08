// gcc -lm -o simulateComplexPsr simulateComplexPsr.c tempo2pred.c cheby2d.c t1polyco.c 
//
// This routine simulates a pulse train based on a tempo2 predictor file
// This can therefore correctly simulate Doppler shifts caused by the Earth's motion and by orbital effects
//

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "T2toolkit.h"
#include <string.h>
#include "tempo2pred.h"
#include "simulate.h"

int main(int argc,char *argv[])
{
  int i,j;
  header *head;
  float *profile;
  int npsamp;
  float tval,phase;
  float fval;
  float tsamp;
  int ival;
  double t;
  double p0;
  float fref,freq;
  float chanbw;
  int nchan;
  double *sum;
  long pn=0; // Pulse number
  FILE *fout;
  char description[128];
  long double mjd0;
  T2Predictor pred;
  long double phase_ld;
  long double phi_0;
  double amp = 1;
  char outName[128] = "pulsar.dat";
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
  chanbw = (head->f2-head->f1)/head->nchan;
  tsamp = (head->tsamp);
  nchan = (head->nchan);
  
  T2Predictor_Init(&pred);
  T2Predictor_Read(&pred,head->predictor);

  fref = 0.5*(head->f2+head->f1);
  mjd0 = 56000.0;
  
  p0 = 1.0/T2Predictor_GetFrequency(&pred,mjd0,fref);
  phi_0 = T2Predictor_GetPhase(&pred,mjd0,fref);
  npsamp =  (int)(p0/tsamp);
  printf("period from predictor is: %.15g\n",p0);
  
  profile = (float *)malloc(sizeof(float)*npsamp*nchan);
  sum = (double *)malloc(sizeof(double)*head->nchan);

  
  // Make a simple profile
  for (i=0;i<npsamp;i++)
    {
      tval = (i)*tsamp;
      phase =  tval/p0-(int)(tval/p0); //+phi0;

      for (j=0;j<nchan;j++)
	{
	  profile[i*nchan+j] = exp(-pow(phase-0.5,2)/2.0/head->width/head->width);
	  if (i==0) 
	    sum[j]=profile[i*head->nchan+j];
	  else
	    sum[j]+=profile[i*head->nchan+j];
	}
    }
  if (head->addScatter==1)
    {
          float *newPulse;
	  int k;
	  float tscat;
	  int iscat;
	  float fghz;
	  float a,b,c,alpha;
	  float sum1=0,sum2=0;
	  
	  a = -6.46;
	  b = 0.154;
	  c = 1.07;
	  alpha = 4.4;
	  
	  newPulse = malloc(sizeof(float)*(int)(npsamp+1));
	  for (j=0;j<nchan;j++)
	    {
	      fghz = (head->f1+j*chanbw)/1000.0;
	      tscat = pow(10,(a+b*log10(head->dm)+c*log10(head->dm)*log10(head->dm)-alpha*log10(fghz)))/1000.0;
	      iscat = (int)(tscat/head->tsamp);
	      printf("tscat: %g %g %d\n",fghz,tscat,iscat);
	      sum1=0;
	      for (i=0;i<npsamp;i++)
		{
		  newPulse[i] = 0;
		  sum1+=profile[i*nchan+j];
		}
	      if (iscat > 0)
		{
		  sum2=0.0;
		  for (i=0;i<npsamp;i++)
		    {
		      for (k=0;k<200;k++)
			{
			  if ((int)i-(int)k > 0)
			    newPulse[i] += profile[(i-k)*nchan+j]*exp(-k/(float)iscat);
			}
		      sum2+=newPulse[i];
		    }
		  for (i=0;i<npsamp;i++)
		    profile[i*nchan+j] = newPulse[i]*(sum1/sum2);
		}
	    }
	  free(newPulse);
	  
    }

  // Add on dispersion smearing
  {
    float dispSmear;
    float f1samp,fl;
    int   nextra,k,s,c[(int)npsamp+1];
    float *newPulse;
    float si = -2;
    newPulse = malloc(sizeof(float)*(int)(npsamp+1));
    for (j=0;j<nchan;j++)
      {
	for (i=0;i<npsamp;i++)
	  {
	    newPulse[i] = profile[i*nchan+j];
	    c[i] = 1;
	  }
	//	dispSmear = fabs(4.15e-3*dm*(pow((f1+j*chanbw-chanbw/2.0)/1000.0,si)-pow((f1+j*chanbw+chanbw/2.0)/1000.0,si)));
	fl = head->f1+j*chanbw-fabs(chanbw)/2.0;
	f1samp = sqrt(1.0/(1.0/(fl/1000.0)/(fl/1000.0)-head->tsamp/4.13e-3/head->dm))*1000.0;
	nextra = (int)(fabs(chanbw)/fabs(f1samp-fl)+0.5);
	printf("fl f1samp = %ld %g %g %d %g\n",j,fl,f1samp,nextra,chanbw);
	//      }
	fref = head->f1+j*chanbw;
	for (k=0;k<nextra;k++)
	  {
	    dispSmear = (4.15e-3*head->dm*(pow(fref/1000.0,si)-pow((fl+k*fabs(chanbw)/(float)nextra)/1000.0,si)));
	    s = (int)((dispSmear/head->tsamp)+0.5);
	    if (j==nchan-1) printf("dispSmear = %d %g %g %g %d\n",k,dispSmear,fl,fl+k*fabs(chanbw)/(float)nextra,s);
	    
	    for (i=0;i<npsamp;i++)
	      {
		if (((int)i+s) >= 0 && ((int)i+s) < npsamp)
		  {
		    newPulse[i] += profile[((int)i+s)*nchan+j];
		    c[i]++;
		  }
	      }
	  }
	for (i=0;i<npsamp;i++)
	  profile[i*nchan+j] = newPulse[i]/(float)(c[i]);
      }
    free(newPulse);
  }

  
  for (i=0;i<npsamp;i++)
    {
      tval = (i)*tsamp;

      for (j=0;j<nchan;j++)
	{
	  if (j==0 || j == head->nchan-1)
	    {
	      printf("profile: %d %g %g ",i,profile[i*nchan+j],phase);
	      if (j==head->nchan-1) printf("\n");
	    }
	  
	}
    }
  //  exit(1);
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
  printf("Writing data\n");
  for (t=head->t0;t<head->t1;t+=tsamp)
    {
      for (i=0;i<nchan;i++)
	{
	  freq = head->f1+chanbw*i;
	  phase_ld = T2Predictor_GetPhase(&pred,mjd0+(long double)t/86400.0L,freq)-phi_0;
	  phase_ld = (fabsl(phase_ld)-(int)(fabsl(phase_ld)));

	  fwrite(&profile[(int)(phase_ld*npsamp)*nchan+i],sizeof(float),1,fout);
	  //	  if (i==0 && t < head->t0+5)
	  //	    printf("result: %g %d %.5f %g\n",(double)phase_ld,(int)(phase_ld*npsamp),t,profile[(int)(phase_ld*npsamp)*nchan+i]);

	}
    }
  T2Predictor_Destroy(&pred);
  printf("Finishing\n");

  fclose(fout);
  simulateReleaseMemory(head);
  free(profile);
  free(sum);
  free(head);
}
