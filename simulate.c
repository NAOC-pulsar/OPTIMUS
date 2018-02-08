// Routines for the simulation software

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "T2toolkit.h"
#include "simulate.h"

void simulateReadHeaderParameters(header *head,FILE *fin)
{
  fread(&head->format,sizeof(char),64,fin);
  if (strcmp(head->format,"FORMAT 1")==0)
    {
      fread(&head->name,sizeof(char),128,fin);
      fread(&head->t0,sizeof(float),1,fin);
      fread(&head->t1,sizeof(float),1,fin);
      fread(&head->tsamp,sizeof(float),1,fin);
      fread(&head->f1,sizeof(float),1,fin);
      fread(&head->f2,sizeof(float),1,fin);
      fread(&head->nchan,sizeof(int),1,fin);
      fread(&head->raj_rad,sizeof(float),1,fin);
      fread(&head->decj_rad,sizeof(float),1,fin);
      fread(&head->useAngle,sizeof(int),1,fin);
      fread(&head->initialSeed,sizeof(long),1,fin);
      head->seed = head->initialSeed;
    }
  else
    {
      printf("ERROR: Unable to process this file format\n");
      exit(1);
    }
}


void simulateDisplayHeaderParameters(header *head)
{
  printf("Format:       %s\n",head->format);
  printf("Name:         %s\n",head->name);
  printf("t0 (sec):     %f\n",head->t0);
  printf("t1 (sec):     %f\n",head->t1);
  printf("tsamp (sec):  %g\n",head->tsamp);
  printf("f1 (MHz):     %f\n",head->f1);
  printf("f2 (MHz):     %f\n",head->f2);
  printf("nchan:        %d\n",head->nchan);
  printf("RAJ (rad):    %g\n",head->raj_rad);
  printf("DECJ (rad):   %g\n",head->decj_rad);  
  printf("Use angle:    %d\n",head->useAngle);
  printf("Random seed:  %d\n",(int)head->initialSeed);
}

void simulateWriteHeader(header *head,FILE *fout)
{
  float fval;
  long lval;
  int ival;
  char format[64];
  strcpy(format,"FORMAT 1");
  fwrite(format,1,64,fout);
  fwrite(head->name,1,128,fout);
  fval = head->t0; fwrite(&fval,sizeof(float),1,fout);
  fval = head->t1; fwrite(&fval,sizeof(float),1,fout);
  fval = head->tsamp; fwrite(&fval,sizeof(float),1,fout);
  fval = head->f1; fwrite(&fval,sizeof(float),1,fout);
  fval = head->f2; fwrite(&fval,sizeof(float),1,fout);
  ival = head->nchan; fwrite(&ival,sizeof(int),1,fout);
  fval = head->raj_rad; fwrite(&fval,sizeof(float),1,fout); // Right ascension
  fval = head->decj_rad; fwrite(&fval,sizeof(float),1,fout); // Declination
  ival = head->useAngle; fwrite(&ival,sizeof(int),1,fout);
  lval = head->initialSeed; fwrite(&lval,sizeof(long),1,fout);
}

//
// Set some sensible defaults for small output files
//
void simulateSetHeaderDefaults(header *head)
{  
  head->f1       = 1518;
  head->f2       = 1230;
  head->nchan    = 96;
  head->imjd     = 56000;
  head->smjd     = 23456;
  head->stt_offs = 0.1234;
  head->t0       = 0;
  head->t1       = 30;
  head->raj_rad  = 0;
  head->decj_rad = 0;
  head->useAngle = 0;
  head->tsamp    = 250e-6;
  head->seed = head->initialSeed = TKsetSeed();
  strcpy(head->name,"Default");
  strcpy(head->predictor,"t2pred.dat");
  head->setGain = 0;
  head->setTsys = 0;
  head->setFlux = 0;
  head->nbits   = 1;
  head->p0 = 0.3;
  head->dm = 50;
  head->addScatter = 0;
  head->addDMsmear = 0;
  head->width = 0.015;
  head->setBeam = 0;
  head->diameter = 300;
  head->surveyType = 1;
  head->nRFI = 0;
  strcpy(head->rfiFile,"NULL");
  head->nb_rfi = (narrowBandRFI *)malloc(sizeof(narrowBandRFI)*MAX_NB_RFI);
  head->n_nb_rfi = 0;
  head->calAmp = 1; // Fraction of Tsys
  head->calFreq = 1; // Hz
  head->calDuty = 0.5; // Duty cycle
}

void simulateReadParamFile(header *head,char *paramFile)
{
  FILE *fin;
  char line[4096];
  char param[1024];
  float setGain = -1;
  float setTsys = -1;
  float setFlux = -1;
  int i;
  
  if (!(fin = fopen(paramFile,"r")))
    {
      printf("Unable to open file %s\n",paramFile);
      exit(1);
    }
  while (!feof(fin))
    {
      if (fgets(line,4096,fin)!=NULL)
	{
	  line[strlen(line)-1]='\0';
	  sscanf(line,"%s",param);
	  if (strcmp(param,"name:")==0)
	    strcpy(head->name,line+strlen(param)+1);
	  else if (strcmp(param,"predictor:")==0)
	    strcpy(head->predictor,line+strlen(param)+1);
	  else if (strcmp(param,"f1:")==0)
	    sscanf(line+strlen(param)+1,"%f",&(head->f1));
	  else if (strcmp(param,"f2:")==0)
	    sscanf(line+strlen(param)+1,"%f",&(head->f2));
	  else if (strcmp(param,"imjd:")==0)
	    sscanf(line+strlen(param)+1,"%d",&(head->imjd));
	  else if (strcmp(param,"smjd:")==0)
	    sscanf(line+strlen(param)+1,"%d",&(head->smjd));
	  else if (strcmp(param,"stt_offs:")==0)
	    sscanf(line+strlen(param)+1,"%d",&(head->stt_offs));
	  else if (strcmp(param,"t0:")==0)
	    sscanf(line+strlen(param)+1,"%f",&(head->t0));
	  else if (strcmp(param,"t1:")==0)
	    sscanf(line+strlen(param)+1,"%f",&(head->t1));
	  else if (strcmp(param,"tsamp:")==0)
	    sscanf(line+strlen(param)+1,"%f",&(head->tsamp));
	  else if (strcmp(param,"diameter:")==0)
	    sscanf(line+strlen(param)+1,"%f",&(head->diameter));
	  else if (strcmp(param,"surveyType:")==0)
	    sscanf(line+strlen(param)+1,"%d",&(head->surveyType));
	  else if (strcmp(param,"nchan:")==0)
	    sscanf(line+strlen(param)+1,"%d",&(head->nchan));
	  else if (strcmp(param,"raj:")==0)
	    sscanf(line+strlen(param)+1,"%f",&(head->raj_rad));
	  else if (strcmp(param,"decj:")==0)
	    sscanf(line+strlen(param)+1,"%f",&(head->decj_rad));
	  else if (strcmp(param,"p0:")==0)
	    sscanf(line+strlen(param)+1,"%f",&(head->p0));
	  else if (strcmp(param,"dm:")==0)
	    sscanf(line+strlen(param)+1,"%f",&(head->dm));
	  else if (strcmp(param,"calAmp:")==0)
	    sscanf(line+strlen(param)+1,"%f",&(head->calAmp));
	  else if (strcmp(param,"calFreq:")==0)
	    sscanf(line+strlen(param)+1,"%f",&(head->calFreq));
	  else if (strcmp(param,"calDuty:")==0)
	    sscanf(line+strlen(param)+1,"%f",&(head->calDuty));
	  else if (strcmp(param,"scatter:")==0)
	    sscanf(line+strlen(param)+1,"%d",&(head->addScatter));
	  else if (strcmp(param,"dmSmear:")==0)
	    sscanf(line+strlen(param)+1,"%d",&(head->addDMsmear));
	  else if (strcmp(param,"width:")==0)
	    sscanf(line+strlen(param)+1,"%f",&(head->width));	  
	  else if (strcmp(param,"useAngle:")==0)
	    sscanf(line+strlen(param)+1,"%d",&(head->useAngle));
	  else if (strcmp(param,"nbits:")==0)
	    sscanf(line+strlen(param)+1,"%d",&(head->nbits));
	  else if (strcmp(param,"nrfi:")==0)
	    sscanf(line+strlen(param)+1,"%d",&(head->nRFI));
	  else if (strcmp(param,"rfi_file:")==0)
	    strcpy(head->rfiFile,line+strlen(param)+1);	  
	  else if (strcmp(param,"rfi_narrowBand:")==0)
	    {
	      float f1,f2,amp,t0,duration;
	      int nread;
	      nread = sscanf(line+strlen(param)+1,"%f %f %f %f %f",&f1,&f2,&amp,&t0,&duration);
	      head->nb_rfi[head->n_nb_rfi].f1 = f1;
	      head->nb_rfi[head->n_nb_rfi].f2 = f2;
	      head->nb_rfi[head->n_nb_rfi].amp = amp;
	      if (nread == 3)
		{
		  head->nb_rfi[head->n_nb_rfi].t0 = -1;
		  head->nb_rfi[head->n_nb_rfi].duration = -1;
		}
	      else
		{
		  head->nb_rfi[head->n_nb_rfi].t0 = t0;
		  head->nb_rfi[head->n_nb_rfi].duration = duration;
		}
	      
	      (head->n_nb_rfi)++;
	      printf("Loaded %g %g\n",f1,f2);
	    }
	  else if (strcmp(param,"beamRA0:")==0)
	    {sscanf(line+strlen(param)+1,"%f",&(head->beamRA0));head->setBeam=1;}
	  else if (strcmp(param,"beamDEC0:")==0)
	    {sscanf(line+strlen(param)+1,"%f",&(head->beamDEC0)); head->setBeam=1;}
	  else if (strcmp(param,"seed:")==0)
	    {
	      sscanf(line+strlen(param)+1,"%ld",&(head->initialSeed));
	      head->seed = head->initialSeed;
	    }
	  else if (strcmp(param,"gain:")==0)
	    sscanf(line+strlen(param)+1,"%f",&setGain);
	  else if (strcmp(param,"flux:")==0)
	    sscanf(line+strlen(param)+1,"%f",&setFlux);
	  else if (strcmp(param,"tsys:")==0)
	    sscanf(line+strlen(param)+1,"%f",&setTsys);
	}
    }
  
  fclose(fin);

  if (setGain > 0)
    {
      head->gain = (float *)malloc(sizeof(float)*head->nchan);
      for (i=0;i<head->nchan;i++)
	head->gain[i] = setGain;
      head->setGain = 1;
    }
  if (setTsys > 0)
    {
      head->tsys = (float *)malloc(sizeof(float)*head->nchan);
      for (i=0;i<head->nchan;i++)
	head->tsys[i] = setTsys;
      head->setTsys=1;
    }
  if (setFlux > 0)
    {
      head->flux = (float *)malloc(sizeof(float)*head->nchan);
      for (i=0;i<head->nchan;i++)
	head->flux[i] = setFlux;
      head->setFlux = 1;
    }
}

void simulateReleaseMemory(header *head)
{
  if (head->setGain==1)
    free(head->gain);
  if (head->setTsys==1)
    free(head->tsys);
  if (head->setFlux==1)
    free(head->flux);
  free(head->nb_rfi);
}
