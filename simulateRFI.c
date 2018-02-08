// Simulate of RFI

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
  float val,valc;
  
  char outName[1024];
  int i,j;
  float freq;
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
  unsigned long *rfiStart;
  unsigned long *rfiEnd;
  float *rfiAmp;
  unsigned long nsamp;
  unsigned long samp;
  float *loadRFIAmp;
  float *loadRFIWidth;
  unsigned long nLoadRFI=0;
  unsigned long rfiNumber;
  FILE *fin;
  int goodtime;
  
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
  rfiStart = (unsigned long *)malloc(sizeof(unsigned long)*head->nRFI);
  rfiEnd = (unsigned long *)malloc(sizeof(unsigned long)*head->nRFI);
  rfiAmp = (float *)malloc(sizeof(float)*head->nRFI);
  loadRFIAmp = (float *)malloc(sizeof(float)*MAX_RFI);
  loadRFIWidth = (float *)malloc(sizeof(float)*MAX_RFI);
  
  printf("Loading RFI events\n");
  fin = fopen(head->rfiFile,"r");
  while (!feof(fin))
    {
      if (fscanf(fin,"%f %f",&loadRFIAmp[nLoadRFI],&loadRFIWidth[nLoadRFI])==2)
	nLoadRFI++;
    }
  fclose(fin);
  printf("Complete loading RFI events: %d\n",(int)nLoadRFI);
  printf("Creating %d RFI signals\n",(int)head->nRFI);
  nsamp = (int)((head->t1-head->t0)/(head->tsamp));
  for (i=0;i<head->nRFI;i++)
    {
      rfiNumber = (int)(TKranDev(&head->seed)*(nLoadRFI-1)+0.5);
      
      rfiStart[i] = TKranDev(&head->seed)*nsamp;
      rfiEnd[i] = rfiStart[i] + loadRFIWidth[rfiNumber]/head->tsamp;
      rfiAmp[i] = loadRFIAmp[rfiNumber];
      if (rfiEnd[i] >= nsamp) rfiEnd[i] = nsamp-1;
    }
  chanbw = fabs(head->f2-head->f1)/head->nchan;
  
  fout = fopen(outName,"wb");
  if (debugOut==1)
    fout2 = fopen(debugOutName,"w");
  //
  // Write the header information
  simulateWriteHeader(head,fout);
  samp = 0;
  for (t=head->t0;t<head->t1;t+=head->tsamp)
    {
      if (debugOut == 1 && t >= debugOutT1 && t < debugOutT2)
	fprintf(fout2,"%f ",t);
      valc = 0;
      // Must speed this up!
      for (j=0;j<head->nRFI;j++)
	{
	  if (samp >= rfiStart[j] && samp <= rfiEnd[j]) {valc=rfiAmp[j]; break;}
	}
      for (i=0;i<head->nchan;i++)
	{
	  freq = head->f1 + i*(head->f2-head->f1)/head->nchan;
	  val = valc;
	  // Check persistent RFI
	  for (j=0;j<head->n_nb_rfi;j++)
	    {
	      goodtime = 1;
	      //	      printf("Checking %g %g %g\n",head->nb_rfi[j].t0,head->nb_rfi[j].duration,t);
	      if (head->nb_rfi[j].t0 >= 0 && (t < head->nb_rfi[j].t0 || t > head->nb_rfi[j].t0 + head->nb_rfi[j].duration))
		  goodtime = 0;
	      if (freq >= head->nb_rfi[j].f1 && freq <= head->nb_rfi[j].f2 && goodtime == 1)
		  val += head->nb_rfi[j].amp;
	    }
	  fwrite(&val,sizeof(float),1,fout);
	  if (debugOut == 1 && t >= debugOutT1 && t < debugOutT2 &&
	      i >= debugOutF1 && i < debugOutF2)
	    fprintf(fout2,"%f ",val);
	}
      if (debugOut == 1 && t >= debugOutT1 && t < debugOutT2)
	fprintf(fout2,"\n");
      samp++;
    }
  fclose(fout);
  if (debugOut == 1)
    fclose(fout2); 

  simulateReleaseMemory(head);
  free(head);
  free(rfiStart);
  free(rfiEnd);
  free(rfiAmp);
  free(loadRFIAmp);
  free(loadRFIWidth);
}

