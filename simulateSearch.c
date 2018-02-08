// 1. Must set lenProc correctly
//
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "T2toolkit.h"

void calculatePulseTrain(float t1,float t2,float dt,float f1,float f2,int nchan);

int main(int argc,char *argv[])
{
  float tsamp; // Sample time in seconds
  float t0,t1; // Start and end time in seconds from nominal LST time
  float timeVal;
  int ncvr=1;
  int rn;
  int nchan=32;
  float f1 = 1300;
  float f2 = 1400;
  int fchan;

  //  tsamp = 1024e-6;
  tsamp = 1e-3;
  t0 = -120;
  t1 = +120;

  // First calculate the pulse train for each pulsar as a function of frequency
  // Either write to disk (in a binary format) or keep in memory
  calculatePulseTrain(t0,t1,tsamp,f1,f2,nchan);
  exit(1);
  for (timeVal = t0;timeVal < t1; timeVal+=tsamp)
    {
      printf("Time: %g\n",timeVal);
      // On moderate time scale:
      // Update:
      // - beam position
      // - angles for each beam to each source


      // For every pulse period for a given pulsar
      // Update scaling factor for the flux density (to deal with nulling etc.)
      // Have a pointer to pulsars that need updating

      
      // Process each receiver
      for (rn=0;rn<ncvr;rn++)
	{
	  // Determine where this receiver is pointing at this time
	  // Use LST
	  
	  // Process each frequency channel
	  for (fchan = 0;fchan < nchan; fchan++)
	    {
	      // Get noise level from the system temperature for
	      // this receiver in this frequency channel

	      // Get sky temperature estimation


	      // Add sources
	      // Determine angle to each source (perhaps only update this on a slower cycle time - use timeVal to update these)
	      // Calculate expected gain
	      

	      // Add all the signals together
	     	    	      
	    }
	  // Digitise the signal
	}
      
      
    }

}

void calculatePulseTrain(float t1,float t2,float dt,float f1,float f2,int nchan)
{
  float p0 = 0.6;
  float w = 0.08;
  float si = -2;
  float timeSpan;
  unsigned long int i,j,k;
  long int ival,kval;
  unsigned long int npulse;
  float s;
  float lenProc = p0*3; 
  double t0,t,deltaT,dm,tval;
  float *signal;
  float amp[nchan],scale,phase,phi0;  
  long seed = TKsetSeed();
  int output = 1; // 1 = ASCII, 2 = binary
  char fname[1024];
  FILE *fout;
  float chanbw = (f2-f1)/nchan;
  float fref = (f2+f1)*0.5;
  float dt_dm[nchan];
  float *pulse;
  float *convolve;
  
  dm = 500;

  // Check processing length
  // Should depend on period, dispersion smearing, scattering, etc.
  lenProc = p0*10; // Should check dispersion time
  pulse = (float *)malloc(sizeof(float)*(int)(lenProc/dt+1)*nchan);
  convolve = (float *)malloc(sizeof(float)*(int)(lenProc/dt+1));
  phi0=0.5;
  
  // Create the pulse shape
  fout = fopen("pulseShape.dat","w");
  for (i=0;i<lenProc/dt;i++)
    {
      tval = (i-lenProc/dt/2.0)*dt;
      phase =  tval/p0; //+phi0;
      for (j=0;j<nchan;j++)
	{
	  if (i*nchan+j > (int)(lenProc/dt+1)*nchan)
	    {
	      printf("THIS should not happen\n");
	      printf("i = %ld\n",i);
	      printf("j = %ld\n",j);
	      printf("nchan = %d\n",nchan);
	      printf("lenProc = %g\n",lenProc);
	      printf("lenProc/dt = %g\n",lenProc/dt);
	      exit(1);
	    }
	  pulse[i*nchan+j] = exp(-pow(phase,2)/2.0/w/w); 
	}

    }
  // Add on scattering tail
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
    
    newPulse = malloc(sizeof(float)*(int)(lenProc/dt+1));
    for (j=0;j<nchan;j++)
      {
	fghz = (f1+j*chanbw)/1000.0;
	tscat = pow(10,(a+b*log10(dm)+c*log10(dm)*log10(dm)-alpha*log10(fghz)))/1000.0;
	iscat = (int)(tscat/dt);
	printf("tscat: %g %g %d\n",fghz,tscat,iscat);
	sum1=0;
	for (i=0;i<lenProc/dt;i++)
	  {
	    newPulse[i] = 0;
	    sum1+=pulse[i*nchan+j];
	  }
	if (iscat > 0)
	  {
	    sum2=0.0;
	    for (i=0;i<lenProc/dt;i++)
	      {
		for (k=0;k<200;k++)
		  {
		    if ((int)i-(int)k > 0)
		      newPulse[i] += pulse[(i-k)*nchan+j]*exp(-k/(float)iscat);
		  }
		sum2+=newPulse[i];
	      }
	    for (i=0;i<lenProc/dt;i++)
	      pulse[i*nchan+j] = newPulse[i]*(sum1/sum2);
	  }
      }
    free(newPulse);
  }

  
  // Add on dispersion smearing
  {
    float dispSmear;
    float f1samp,fl;
    int   nextra,k,s,c[(int)(lenProc/dt)+1];
    float *newPulse;
    newPulse = malloc(sizeof(float)*(int)(lenProc/dt+1));
    for (j=0;j<nchan;j++)
      {
	for (i=0;i<lenProc/dt;i++)
	  {
	    newPulse[i] = pulse[i*nchan+j];
	    c[i] = 1;
	  }
	//	dispSmear = fabs(4.15e-3*dm*(pow((f1+j*chanbw-chanbw/2.0)/1000.0,si)-pow((f1+j*chanbw+chanbw/2.0)/1000.0,si)));
	fl = f1+j*chanbw-fabs(chanbw)/2.0;
	f1samp = sqrt(1.0/(1.0/(fl/1000.0)/(fl/1000.0)-dt/4.13e-3/dm))*1000.0;
	nextra = (int)(fabs(chanbw)/fabs(f1samp-fl)+0.5);
	printf("fl f1samp = %ld %g %g %d %g\n",j,fl,f1samp,nextra,chanbw);

	fref = f1+j*chanbw;
	for (k=0;k<nextra;k++)
	  {
	    dispSmear = (4.15e-3*dm*(pow(fref/1000.0,si)-pow((fl+k*fabs(chanbw)/(float)nextra)/1000.0,si)));
	    s = (int)((dispSmear/dt)+0.5);
	    for (i=0;i<lenProc/dt;i++)
	      {
		if (((int)i+s) >= 0 && ((int)i+s) < lenProc/dt)
		  {
		    newPulse[i] += pulse[((int)i+s)*nchan+j];
		    c[i]++;
		  }
	      }
	  }
	for (i=0;i<lenProc/dt;i++)
	  pulse[i*nchan+j] = newPulse[i]/(float)(c[i]);
      }
    free(newPulse);
  }
  
  // Now print the pulse shape to a file
  for (i=0;i<lenProc/dt;i++)
    {
      fprintf(fout,"%ld ",i);
      for (j=0;j<nchan;j++)
	fprintf(fout,"%g ",pulse[i*nchan+j]);
      fprintf(fout,"\n");      
    }
  fclose(fout);
  
  strcpy(fname,"pulseSeq1.dat");
  
  if (!(signal = (float *)calloc(((t2-t1)/dt*nchan+1),sizeof(float))))
    {
      printf("Unable to allocate memory for calculating the pulse sequence\n");
      exit(1);
    }

  
  timeSpan = (t2-t1)/dt;
  npulse = (int)(timeSpan/p0+0.5);
  printf("Npulse = %ld\n",npulse);

  // Calculate time offsets for the DM
  // Also calculate the flux density
  for (i=0;i<nchan;i++)
    {
      dt_dm[i] = (4.15e-3*dm*(pow(fref/1000.0,si)-pow((f1+i*chanbw)/1000.0,si)));
      amp[i]   = 1.0; //+i*0.1;
      // Can calculate a jitter value here
      printf("dt_dm %ld %g %g.  Amp = %g\n",i,dt_dm[i],dt_dm[i]/dt,amp[i]);
    }
  for (i=0;i<npulse;i++)
    {
      //
      // Can set up dynamic spectrum information here
      // Can also consider drifting subpulses
      //
     
      if (i == 9*npulse/10.)	  printf(".........\n");
      else if (i == 8*npulse/10.) printf("........\n");
      else if (i == 7*npulse/10.) printf(".......\n");
      else if (i == 7*npulse/10.) printf(".......\n");
      else if (i == 6*npulse/10.) printf("......\n");
      else if (i == 5*npulse/10.) printf(".....\n");
      else if (i == 4*npulse/10.) printf("....\n");
      else if (i == 3*npulse/10.) printf("...\n");
      else if (i == 2*npulse/10.) printf("..\n");
      else if (i == npulse/10.)   printf(".\n");
      t0 = t1+i*p0;
      scale = TKranDev(&seed);
      //      printf("i = %ld/%ld %g %g %g %g\n",i,npulse,t0,t0-lenProc,t0+lenProc,dt);
      //      for (t = t0-lenProc;t < t0+lenProc;t += dt)       
      for (k=0;k<lenProc/dt;k++)
      {
	  
	  for (j=0;j<nchan;j++)
	    {
	      //	  if (i==32886)
	      //	    printf("In here %g %g %g\n",t,t0+lenProc,dt);
	      //	  printf("t = %g %g %g %g\n",t,t0,lenProc,dt);
	      //	      if (t >= t1 && t < t2)
		{
		  //		  phase =  (t-t0)/p0+phi0;
		  s = scale*amp[j];
		  //		  ival = ((t-t1-dt_dm[j])/dt)*nchan+j;
		  //		  ival = ((t-t1)/dt)*nchan+j;
		  //		  ival = k;
		  //		  ival -= (int)(dt_dm[j]/dt)*nchan;
		  //	      printf("t = %g ival = %d %g %g %g\n",t,(int)ival,t0,lenProc,dt);
		  ival = ((int)((i*p0+(k*dt-lenProc/2.)-dt_dm[j])/dt))*nchan+j;
		  kval = k*nchan+j;
		  //		  printf("ival = %d\n",(int)ival);
		  if (ival >= 0 && ival < (t2-t1)/dt*nchan)
		    signal[ival]+=s*pulse[kval];
		}
	    }
      }
    }
  printf("Writing output\n");
  fout = fopen(fname,"w");
  
  for (i=0;i<(t2-t1)/dt;i++)
    {
      fprintf(fout,"%ld ",i);
      for (j=0;j<nchan;j++)
	fprintf(fout,"%g ",signal[i*nchan+j]);
      fprintf(fout,"\n");
    }
  fclose(fout);
  free(signal);  
  free(pulse);
  free(convolve);
}




