// I need to set lenProc correctly
//
// Software to simulate a pulse sequence from a pulsar
//
// This software accounts for:
// 1. Dispersion smearing
// 2. A simple model for scattering
//
// The pulse sequence can be defined by:
// 1. A constant period (P), pulse width and amplitude (as a function of frequency)
//
// gcc -lm -o simulatePsr simulatePsr.c T2toolkit.c -O3
//
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "T2toolkit.h"
#include <string.h>

void calculatePulseTrain(float t1,float t2,float dt,float f1,float f2,int nchan,FILE *fout);

int main(int argc,char *argv[])
{
  int i;
  float tsamp; // Sample time in seconds
  float t0,t1; // Start and end time in seconds from nominal LST time
  float timeVal;
  int ncvr=1;
  int rn;
  int nchan=96;
  float f1 = 1518;
  float f2 = 1230;
  int fchan;
  FILE *fout;
  char outName[1024];
  
  // Read input parameters
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-o")==0)
	strcpy(outName,argv[++i]);
    }
  
  //  tsamp = 1024e-6;
  tsamp = 250e-6;
  t0    = -1050;
  t1    = +1050;

  if (!(fout = fopen(outName,"wb")))
    {
      printf("Unable to open output file: %s (-o option)\n",outName);
      exit(1);
    }
  // Write the header information
  {
    char description[128];
    float fval;
    int ival;

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
  }

  
  // First calculate the pulse train for each pulsar as a function of frequency
  // Either write to disk (in a binary format) or keep in memory
  calculatePulseTrain(t0,t1,tsamp,f1,f2,nchan,fout);
  fclose(fout);
}

void calculatePulseTrain(float t1,float t2,float dt,float f1,float f2,int nchan,FILE *fout)
{
  float p0 = 1.1;
  float w = p0/10.0/2.0;
  float si = -2;
  float timeSpan;
  unsigned long int i,j,k;
  long int ival,kval;
  unsigned long int npulse;
  float s;
  float lenProc = p0; //*3; 
  int lenProcPeriod = 1;
  double t0,t,deltaT,dm,tval;
  float *signal;
  float amp[nchan],scale,phase,phi0;  
  long seed = TKsetSeed();
  int output = 1; // 1 = ASCII, 2 = binary
  char fname[1024];
  FILE *fout2;
  float chanbw = (f2-f1)/nchan;
  float fref = (f2+f1)*0.5;
  float dt_dm[nchan];
  float *pulse;
  float *convolve;
  
  dm = 500;

  // Check processing length
  // Should depend on period, dispersion smearing, scattering, etc.
  if (lenProcPeriod == 1)
    lenProc = p0; // Should check dispersion time
  pulse = (float *)malloc(sizeof(float)*(int)(lenProc/dt+1)*nchan);
  convolve = (float *)malloc(sizeof(float)*(int)(lenProc/dt+1));
  phi0=0.5;
  
  // Create the pulse shape
  fout2 = fopen("pulseShape.dat","w");
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
      fprintf(fout2,"%ld ",i);
      for (j=0;j<nchan;j++)
	fprintf(fout2,"%g ",pulse[i*nchan+j]);
      fprintf(fout2,"\n");      
    }
  fclose(fout2);
  
  strcpy(fname,"pulseSeq1.dat");


  // I should load the entire data set into memory
  // Should only write in chunks similar to procLen
  //
  if (!(signal = (float *)calloc(((t2-t1)/dt*nchan+1),sizeof(float))))
    {
      printf("Unable to allocate memory for calculating the pulse sequence\n");
      exit(1);
    }


  // Calculate time offsets for the DM
  // Also calculate the flux density
  for (i=0;i<nchan;i++)
    {
      dt_dm[i] = (4.15e-3*dm*(pow(fref/1000.0,si)-pow((f1+i*chanbw)/1000.0,si)));
      amp[i]   = 0.8; //+i*0.1;
      // Can calculate a jitter value here
      printf("dt_dm %ld %g %g.  Amp = %g\n",i,dt_dm[i],dt_dm[i]/dt,amp[i]);
    }

  timeSpan = (t2-t1)/dt;
  npulse = (int)(timeSpan/p0+0.5);
  printf("Npulse = %ld\n",npulse);

  for (i=0;i<npulse;i++)
    {
      //
      // Can set up dynamic spectrum information here
      // Can also consider drifting subpulses
      //
      if (i == (int)(9*npulse/10.))	  printf(".........\n");
      else if (i == (int)(8*npulse/10.)) printf("........\n");
      else if (i == (int)(7*npulse/10.)) printf(".......\n");
      else if (i == (int)(7*npulse/10.)) printf(".......\n");
      else if (i == (int)(6*npulse/10.)) printf("......\n");
      else if (i == (int)(5*npulse/10.)) printf(".....\n");
      else if (i == (int)(4*npulse/10.)) printf("....\n");
      else if (i == (int)(3*npulse/10.)) printf("...\n");
      else if (i == (int)(2*npulse/10.)) printf("..\n");
      else if (i == (int)(npulse/10.))   printf(".\n");
      else if (i == (int)10.)   printf("X %ld\n",npulse/10);
      else if (i == (int)100.)   printf("XX\n");
      else if (i == (int)1000.)   printf("XXX\n");
      t0 = t1+i*p0;
      scale = 1.0; //TKranDev(&seed);
      //      printf("i = %ld/%ld %g %g %g %g\n",i,npulse,t0,t0-lenProc,t0+lenProc,dt);
      //      for (t = t0-lenProc;t < t0+lenProc;t += dt)       
      for (k=0;k<lenProc/dt;k++)
      {
	kval = k*nchan;
	
	for (j=0;j<nchan;j++)
	  {
	    s = scale*amp[j];
	    ival = ((int)((i*p0+(k*dt-lenProc/2.)-dt_dm[j])/dt))*nchan+j;	 // Does not take long
	    if (ival >= 0 && ival < (t2-t1)/dt*nchan) 
	      {
		//		if (lenProcPeriod==1)
		signal[ival] = s*pulse[kval]; 
		  //		else
		  //		  signal[ival] += s*pulse[kval]; // The += takes a long time. Just = is quick
	      }
	    kval++;
	  }
      }
    }
   
  printf("Writing output\n");
  //  fout2 = fopen(fname,"w");
  
  for (i=0;i<(t2-t1)/dt;i++)
    {
      //      fprintf(fout2,"%ld ",i);
      for (j=0;j<nchan;j++)
	fwrite(&signal[i*nchan+j],sizeof(float),1,fout);
      //	fprintf(fout,"%g ",signal[i*nchan+j]);
      //      fprintf(fout2,"\n");
    }
  //  fclose(fout2);
  free(signal);  
  free(pulse);
  free(convolve);
}
//23 sec
//23 sec
//24 sec
