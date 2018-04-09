//  Copyright (C) 2015,2016 George Hobbs
// This file is part of the pfits software package
//

/* pfits is free software: you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version. 
 * pfits is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details. 
 * You should have received a copy of the GNU General Public License 
 * along with pfits.  If not, see <http://www.gnu.org/licenses/>. 
*/

//  gcc -lm -o pfits_plot pfits_plot.c pfits_setup.c pfits_loader.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -lcpgplot -lpgplot
//gcc -lm -o pfits_plot pfits_plot.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio 

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>

typedef struct plotStruct {
  int freqTimeRes;  // resolution for frequency-time plots
  float t1,t2;      // Time range
  int rangeType;    // 1 = single subint, 2 = subint ranges, 3 = time range from start of observation
  float minVal;     // Minimum value
  float maxVal;     // Maximum value
  float minx,maxx;  // Plot ranges
  float miny,maxy;
  float ominx,omaxx;
  float ominy,omaxy;

  int windowNum;    // PGPLOT window number
  int xpanelNum;    // PGPLOT panel number (x)
  int ypanelNum;    // PGPLOT panel number (y)
  int colourPlot;   // Colour plot (1 = colour, -1 = grayscale)

  int xPlotType;    // 1 = sample number, 2 = time from start of observation
  long nSamples;  
  int nTimeSamples;      // Number of points
  int nFrequencySamples; // Number of frequency samples
  int polPlot;      // 1 = different pols on different panels. 2 = different pols on same panel (p1,p2,p3,p4)
} plotStruct;

void doPlot(dSetStruct **dSet,int nFiles,plotStruct **plot,int debug);
void initialisePlot(plotStruct *plot);
void loadPlotData(int nFiles,float **plotArr_p0,float **plotArr_p1,float **plotArr_p2,float **plotArr_p3,plotStruct **plot,dSetStruct **dSet,int debug);
void setMinMaxVals(dSetStruct **dSet,plotStruct **plot,int nFiles,int nchan,int nsblk,float **plotArr_p0,float **plotArr_p1,float **plotArr_p2,float **plotArr_p3);
void drawPlot(int nFiles,int fileNum,int pol,dSetStruct **dSet,plotStruct **plot,float **plotArr_p0,
	      float **plotArr_p1,float **plotarr_p2,float **plotArr_p3);
void drawColourMap(dSetStruct *dSet,int pol,plotStruct *plot,float *plotArr_p0,
		   float *plotArr_p1,float *plotArr_p2,float *plotArr_p3);
void drawLinePlot(dSetStruct *dSet,int pol,plotStruct *plot,float *plotArr_p0,
		  float *plotArr_p1,float *plotArr_p2,float *plotArr_p3);


int main(int argc,char *argv[])
{
  dSetStruct **dSet;
  int debug=0;
  int i;
  int nFiles=0,nFiles_t=0;
  plotStruct **plot;
  float sval=-1;
  float t1=-1,t2=-1;
  float s1=-1,s2=-1;
  
  // Count number of files
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	nFiles++;
    }


  plot = (plotStruct **)malloc(sizeof(plotStruct *)*nFiles);
  for (i=0;i<nFiles;i++)
    {
      plot[i] = (plotStruct *)malloc(sizeof(plotStruct));
      initialisePlot(plot[i]);
    }
  
  // Allocate memory for these files
  dSet = (dSetStruct **)malloc(sizeof(dSetStruct *)*nFiles);
  
  // Initialise everything
  for (i=0;i<nFiles;i++) initialise(&(dSet[i]),debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet[nFiles_t++],debug);
      else if (strcmp(argv[i],"-s")==0)
	sscanf(argv[++i],"%f",&sval);
      else if (strcmp(argv[i],"-s1")==0)
	sscanf(argv[++i],"%f",&s1);
      else if (strcmp(argv[i],"-s2")==0)
	sscanf(argv[++i],"%f",&s2);
      else if (strcmp(argv[i],"-t1")==0)
	sscanf(argv[++i],"%f",&t1);
      else if (strcmp(argv[i],"-t2")==0)
	sscanf(argv[++i],"%f",&t2);
    }

  // Open the files and load header information
  for (i=0;i<nFiles;i++)
    {
      pfitsOpenFile(dSet[i],debug);
      pfitsLoadHeader(dSet[i],debug);
    }
  
  // Setup the plots
  if (sval != -1)
    {
      for (i=0;i<nFiles;i++)
	{
	  plot[i]->t2=plot[i]->t1=sval;
	  plot[i]->rangeType=1;
	}
    }
  else if (s1 != -1 && s2 != -1)
    {
      for (i=0;i<nFiles;i++)
	{
	  plot[i]->t1 = s1;
	  plot[i]->t2 = s2;
	  plot[i]->rangeType = 2;
	}
    }
  else if (t1!=-1 && t2==-1) // If a single time is given
    {
      for (i=0;i<nFiles;i++)
	{
	  plot[i]->t2=plot[i]->t1=(int)(t1/(dSet[i]->head->nsblk*dSet[i]->head->tsamp));
	  printf("Set subint = %g\n",plot[i]->t2);
	  plot[i]->rangeType=1;
	}
    }
  else
    {
      for (i=0;i<nFiles;i++)
	{
	  plot[i]->t2=t2;
	  plot[i]->t1=t1;
	  plot[i]->rangeType=3;
	}
    }


  // Do the plot
  doPlot(dSet,nFiles,plot,debug);
 
  // Close the file
  //  pfitsCloseFile(dSet,debug);
  
  // De-allocate the memory
  for (i=0;i<nFiles;i++)
    {
      deallocateMemory(&dSet[i],debug);
      free(plot[i]);
    }
  
  free(dSet);
  free(plot);
}

void doPlot(dSetStruct **dSet,int nFiles,plotStruct **plot,int debug)
{
  float **plotArr_p0;
  float **plotArr_p1;
  float **plotArr_p2;
  float **plotArr_p3;
  
  float tr[6];
  int nsblk;
  int nchan;
  char key;
  float mx,my,mx2,my2;
  int i,j,k,p;
  float dmVal[64];
  float dmVal_time[64];
  int nDM_curve=0;
  int colourPlot=1;
  int xDisplay = 0; // 0 = sample number, 1 = time in seconds from start of observation
  int nWin=0;
  int curWin;
  int nPanelX=0;
  int nPanelY=0;
  char tStr[128];
  int nsamp;
  int nsubLoad;
  
  // Set up defaults and allocate memory
  nchan = dSet[0]->head->nchan;
  nsblk = dSet[0]->head->nsblk;

  nsubLoad = (plot[0]->t2-plot[0]->t1)+1;
  nsamp = nchan*nsblk*nsubLoad;
  
  printf("Source: %s\n",dSet[0]->head->source);
  printf("Npol: %d\n",dSet[0]->head->npol);
  printf("nsblk: %d\n",nsblk);
  printf("nchan: %d\n",nchan);
  printf("nbits: %d\n",dSet[0]->head->nbits);
  plotArr_p0 = (float **)malloc(sizeof(float *)*nFiles);
  plotArr_p1 = (float **)malloc(sizeof(float *)*nFiles);
  plotArr_p2 = (float **)malloc(sizeof(float *)*nFiles);
  plotArr_p3 = (float **)malloc(sizeof(float *)*nFiles);

  for (i=0;i<nFiles;i++)
    {
      plotArr_p0[i] = (float *)malloc(sizeof(float)*nsamp);  
      if (dSet[0]->head->npol==4)
	{
	  plotArr_p1[i] = (float *)malloc(sizeof(float)*nsamp);  
	  plotArr_p2[i] = (float *)malloc(sizeof(float)*nsamp);  
	  plotArr_p3[i] = (float *)malloc(sizeof(float)*nsamp);  
	}
      else if (dSet[0]->head->npol==2)
	{
	  plotArr_p1[i] = (float *)malloc(sizeof(float)*nsamp);  
	  plotArr_p2[i] = (float *)malloc(sizeof(float));  
	  plotArr_p3[i] = (float *)malloc(sizeof(float));  
	}
      else
	{
	  plotArr_p1[i] = (float *)malloc(sizeof(float));  
	  plotArr_p2[i] = (float *)malloc(sizeof(float));  
	  plotArr_p3[i] = (float *)malloc(sizeof(float));  
	}
    }
  // Load the data
  loadPlotData(nFiles,plotArr_p0,plotArr_p1,plotArr_p2,plotArr_p3,plot,dSet,debug);

  // Setup windows and panels
  if (nFiles==1 && dSet[0]->head->npol==1)
    {
      plot[0]->windowNum = 1;
      plot[0]->xpanelNum = 1;
      plot[0]->ypanelNum = 1;
      nWin=1;
      nPanelX=1;
      nPanelY=1;      
    }
  else if (nFiles==1 && dSet[0]->head->npol==4 && plot[0]->polPlot == 1)
    {
      plot[0]->windowNum = 1;
      plot[0]->xpanelNum = 1;
      plot[0]->ypanelNum = 1;
      nWin=1;
      nPanelX=2;
      nPanelY=2;      
    }
  else if (nFiles==1 && dSet[0]->head->npol==2 && plot[0]->polPlot == 1)
    {
      plot[0]->windowNum = 1;
      plot[0]->xpanelNum = 1;
      plot[0]->ypanelNum = 1;
      nWin=1;
      nPanelX=2;
      nPanelY=1;      
    }
  else if (nFiles==1 && dSet[0]->head->npol==4 && plot[0]->polPlot == 2)
    {
      plot[0]->windowNum = 1;
      plot[0]->xpanelNum = 1;
      plot[0]->ypanelNum = 1;
      nWin=1;
      nPanelX=1;
      nPanelY=1;      
    }
  else
    {
      for (i=0;i<nFiles;i++)
	{
	  plot[i]->windowNum = 1;
	  plot[i]->xpanelNum = 1;
	  plot[i]->ypanelNum = i+1;
	}
      nWin=1;
      nPanelX = 1;
      nPanelY = nFiles;
    }


  setMinMaxVals(dSet,plot,nFiles,nchan,nsblk,plotArr_p0,plotArr_p1,plotArr_p2,plotArr_p3);
    
//  do {
      // Go through each window being plotted
      curWin = -1;
      for (k=0;k<nFiles;k++)
	{
	  for (p=0;p<dSet[k]->head->npol;p++)
	    {
	      if (plot[k]->windowNum!=curWin)
		{
		  // Open a window
		  //sprintf(tStr,"%d.eps/CPS",plot[k]->windowNum);
		  sprintf(tStr,"%s.eps/CPS",dSet[0]->head->source);
                  //sprintf(tStr,"?",plot[k]->windowNum);
		  cpgbeg(0,tStr,nPanelX,nPanelY);
		  cpgsch(1.4);  cpgscf(2);  cpgslw(2);
		  cpgask(0);
		  curWin = plot[k]->windowNum;
		}
	      // Find the correct panel
	      cpgpanl(plot[k]->xpanelNum,plot[k]->ypanelNum);
	      if (plot[k]->polPlot==1)
		{
		  if (p==0) drawPlot(nFiles,k,p,dSet,plot,plotArr_p0,plotArr_p1,plotArr_p2,plotArr_p3);
		  if (p==1) {cpgpanl(plot[k]->xpanelNum+1,plot[k]->ypanelNum); drawPlot(nFiles,k,p,dSet,plot,plotArr_p0,plotArr_p1,plotArr_p2,plotArr_p3);}
		  if (p==2) {cpgpanl(plot[k]->xpanelNum,plot[k]->ypanelNum+1);drawPlot(nFiles,k,p,dSet,plot,plotArr_p0,plotArr_p1,plotArr_p2,plotArr_p3);}
		  if (p==3) {cpgpanl(plot[k]->xpanelNum+1,plot[k]->ypanelNum+1); drawPlot(nFiles,k,p,dSet,plot,plotArr_p0,plotArr_p1,plotArr_p2,plotArr_p3);}
		}
	      else if (plot[k]->polPlot==2)
		drawPlot(nFiles,k,p,dSet,plot,plotArr_p0,plotArr_p1,plotArr_p2,plotArr_p3);
	    }
	}
    
    //  cpgcurs(&mx,&my,&key);
     // key = 'q';
//      printf("key = q");
/*
      if (key=='c')
	{
	  for (i=0;i<nFiles;i++)
	    plot[i]->colourPlot*=-1;
	}
      else if (key=='+')
	{
	  // Load the data
	  for (i=0;i<nFiles;i++)
	    {
	      (plot[i]->t1)++;
	      (plot[i]->t2)++;
	    }
	  loadPlotData(nFiles,plotArr_p0,plotArr_p1,plotArr_p2,plotArr_p3,plot,dSet,debug);
	  setMinMaxVals(dSet,plot,nFiles,nchan,nsblk,plotArr_p0,plotArr_p1,plotArr_p2,plotArr_p3);
	}
      else if (key=='-')
	{
	  // Load the data
	  for (i=0;i<nFiles;i++)
	    {
	      (plot[i]->t1)--;
	      (plot[i]->t2)--;
	    }
	  loadPlotData(nFiles,plotArr_p0,plotArr_p1,plotArr_p2,plotArr_p3,plot,dSet,debug);
	  setMinMaxVals(dSet,plot,nFiles,nchan,nsblk,plotArr_p0,plotArr_p1,plotArr_p2,plotArr_p3);
	}
      else if (key=='x')
	{
	  for (i=0;i<nFiles;i++)
	    {	      
	      plot[i]->xPlotType++;
	      if (plot[i]->xPlotType==3)
		plot[i]->xPlotType=1;
	    }
	}
      else if (key=='z')
	{
	  cpgband(2,0,mx,my,&mx2,&my2,&key);
	  for (i=0;i<nFiles;i++)
	    {
	       if (plot[i]->xPlotType==1)
		 {
		   plot[i]->minx = mx;
		   plot[i]->maxx = mx2;
		   plot[i]->miny = my;
		   plot[i]->maxy = my2;
		 }
	       else if (plot[i]->xPlotType==2)
		 {
		   float s1,s2;
		   s1 = dSet[i]->head->nsblk*plot[i]->t1*dSet[i]->head->tsamp;
		   s2 = dSet[i]->head->tsamp;
		   plot[i]->minx = (mx - s1)/s2;
		   plot[i]->maxx = (mx2 - s1)/s2;
		   plot[i]->miny = my;
		   plot[i]->maxy = my2;
		 }
	    }
	}
      else if (key=='u')
	{
	  for (i=0;i<nFiles;i++)
	    {
	      plot[i]->minx = plot[i]->ominx;
	      plot[i]->maxx = plot[i]->omaxx;
	      plot[i]->miny = plot[i]->ominy;
	      plot[i]->maxy = plot[i]->omaxy;
	      
	    }
	}
      else if (key=='A')
	printf("Cursor coordinate is (%g,%g)\n",mx,my);
      else
	printf("Unknown key-stroke %c\n",key);

      else if (key=='c')
	colourPlot*=-1;
      else if (key=='x')
	{
	  xDisplay++;
	  if (xDisplay == 2) xDisplay=0;
	  if (xDisplay == 0)
	    {
	      minx = 0;
	      maxx = nsblk;
	    }
	  else
	    {
	      minx = dSet[0]->head->nsblk*plot->t1*dSet[0]->head->tsamp;
	      maxx = dSet[0]->head->nsblk*(plot->t1+1)*dSet[0]->head->tsamp;
	    }
	  ominx = minx;
	  omaxx = maxx;
	}
      else if (key=='d') // Overlay DM curve
	{
	  //  	  printf("Enter DM ");
	  //	  scanf("%f",&dmVal[nDM_curve++]);

	  printf("Enter DM ");
	  scanf("%f",&dmVal[nDM_curve]);
	  printf("Sample number ");
	  scanf("%f",&dmVal_time[nDM_curve++]);
	  } 
*/
  //  } while (key!='q');
//    } while (key != 'q'); 

  cpgend();
  for (i=0;i<nFiles;i++)
    {
      free(plotArr_p0[i]);
      free(plotArr_p1[i]);
      free(plotArr_p2[i]);
      free(plotArr_p3[i]);
    }
  free(plotArr_p0);
  free(plotArr_p1);
  free(plotArr_p2);
  free(plotArr_p3);
}

void initialisePlot(plotStruct *plot)
{
  plot->freqTimeRes = 2048; // Frequency time plot resolution
  plot->t1 = 1;
  plot->t2 = 1;
  plot->rangeType = 1; // Single subint
  plot->colourPlot = 1;
  plot->xPlotType = 1;
  plot->polPlot = 1;
}



void loadPlotData(int nFiles,float **plotArr_p0,float **plotArr_p1,float **plotArr_p2,float **plotArr_p3,plotStruct **plot,dSetStruct **dSet,int debug)
{
  int i;
  calibrateStruct cal;

  cal.baseline_p0[0] = 145;
  cal.baseline_p1[0] = 145;
  cal.baseline_p2[0] = 185;
  cal.baseline_p3[0] = 205;

  for (i=0;i<nFiles;i++)
    {
      pfits_read1pol_float(plotArr_p0[i],0,dSet[i],plot[i]->t1,plot[i]->t2,plot[i]->rangeType,&(plot[i]->nSamples),&(plot[i]->nTimeSamples),&(plot[i]->nFrequencySamples),debug);    
      
      if (dSet[0]->head->npol==4)
	{
	 pfits_read1pol_float(plotArr_p1[i],1,dSet[i],plot[i]->t1,plot[i]->t2,plot[i]->rangeType,&(plot[i]->nSamples),&(plot[i]->nTimeSamples),&(plot[i]->nFrequencySamples),debug);    
	 pfits_read1pol_float(plotArr_p2[i],2,dSet[i],plot[i]->t1,plot[i]->t2,plot[i]->rangeType,&(plot[i]->nSamples),&(plot[i]->nTimeSamples),&(plot[i]->nFrequencySamples),debug);    
	 pfits_read1pol_float(plotArr_p3[i],3,dSet[i],plot[i]->t1,plot[i]->t2,plot[i]->rangeType,&(plot[i]->nSamples),&(plot[i]->nTimeSamples),&(plot[i]->nFrequencySamples),debug);    

	 // Calibrate/scale
	 //	 calibrateScalePols(&cal,plotArr_p0[i],plotArr_p1[i],plotArr_p2[i],plotArr_p3[i],plot[i]->nSamples);
	}
      else if (dSet[0]->head->npol==2)
	{
	 pfits_read1pol_float(plotArr_p1[i],1,dSet[i],plot[i]->t1,plot[i]->t2,plot[i]->rangeType,&(plot[i]->nSamples),&(plot[i]->nTimeSamples),&(plot[i]->nFrequencySamples),debug);    
	}
    }
} 

void setMinMaxVals(dSetStruct **dSet,plotStruct **plot,int nFiles,int nchan,int nsblk,float **plotArr_p0,float **plotArr_p1,float **plotArr_p2,float **plotArr_p3)
{
  int i,j;

  for (i=0;i<nFiles;i++)
    {
      if (dSet[i]->head->nbits == 1)
	{
	  plot[i]->minVal = -1;
	  plot[i]->maxVal = 1;
	}
      else if (dSet[i]->head->nbits == 2)
	{
	  plot[i]->minVal = 0;
	  plot[i]->maxVal = 4;
	}
      else if (dSet[i]->head->nbits == 4)
	{
	  plot[i]->minVal = 0;
	  plot[i]->maxVal = 16;
	}
      else if (dSet[i]->head->nbits == 8)
	{
	  plot[i]->minVal = 0;
	  plot[i]->maxVal = 255;
	}
      else if (dSet[i]->head->nbits == 16)
	{
	  plot[i]->minVal = -3000; // -32768;
	  plot[i]->maxVal = 3000; //32768;
	}
      else if (dSet[i]->head->nbits == 32)
	{
	  plot[i]->minVal = -1024;
	  plot[i]->maxVal = 1024; // This should be set properly
	}
    
      plot[i]->minx = 0;
      plot[i]->maxx = plot[i]->nTimeSamples;

      if (nchan==1)
	{
	  plot[i]->miny = plot[i]->maxy = plotArr_p0[i][0];
	  for (j=0;j<nsblk;j++)
	    {
	      if (plotArr_p0[i][j] < plot[i]->miny) plot[i]->miny = plotArr_p0[i][j];
	      if (plotArr_p0[i][j] > plot[i]->maxy) plot[i]->maxy = plotArr_p0[i][j];
	    }
	  
	}
      else
	{
	  plot[i]->miny = dSet[i]->head->chanFreq[0];
	  if (nFiles > 1)
	    plot[i]->maxy = dSet[nFiles-1]->head->chanFreq[dSet[nFiles-1]->head->nchan-1];
	  else
	    plot[i]->maxy = dSet[i]->head->chanFreq[dSet[i]->head->nchan-1];
	  //	  plot[i]->miny = 0; plot[i]->maxy=2000;
	}

      plot[i]->ominx = plot[i]->minx;
      plot[i]->omaxx = plot[i]->maxx;
      plot[i]->ominy = plot[i]->miny;
      plot[i]->omaxy = plot[i]->maxy;
    }
}

void drawPlot(int nFiles,int fileNum,int pol,dSetStruct **dSet,plotStruct **plot,float **plotArr_p0,
	      float **plotArr_p1,float **plotArr_p2,float **plotArr_p3)
{
  if (dSet[fileNum]->head->nchan > 1)
    {
      drawColourMap(dSet[fileNum],pol,plot[fileNum],plotArr_p0[fileNum],plotArr_p1[fileNum],plotArr_p2[fileNum],plotArr_p3[fileNum]);
    }
  else
    drawLinePlot(dSet[fileNum],pol,plot[fileNum],plotArr_p0[fileNum],plotArr_p1[fileNum],plotArr_p2[fileNum],plotArr_p3[fileNum]);
}

void drawLinePlot(dSetStruct *dSet,int pol,plotStruct *plot,float *plotArr_p0,
		  float *plotArr_p1,float *plotArr_p2,float *plotArr_p3)
{
  int i;
  float plotX[plot->nTimeSamples];


  for (i=0;i<plot->nTimeSamples;i++)
    {
      plotX[i] = i;
      printf("Res %g\n",plotArr_p0[i]);
    }
  //  plot->maxy = 300;
  //    plot->miny = 6e4;
  cpgsvp(0.1,0.9,0.15,0.9);
  cpgswin(plot->minx,plot->maxx,plot->miny,plot->maxy);
  cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
  if (plot->polPlot==1)
    {
      if (pol==0) cpgline(plot->nTimeSamples,plotX,plotArr_p0);
      else if (pol==1) cpgline(plot->nTimeSamples,plotX,plotArr_p1);
      else if (pol==2) cpgline(plot->nTimeSamples,plotX,plotArr_p2);
      else if (pol==3) cpgline(plot->nTimeSamples,plotX,plotArr_p3);
    }
  else
    {
      cpgline(plot->nTimeSamples,plotX,plotArr_p0);
      cpgsci(2); cpgline(plot->nTimeSamples,plotX,plotArr_p1);
      cpgsci(3); cpgline(plot->nTimeSamples,plotX,plotArr_p2);
      cpgsci(4); cpgline(plot->nTimeSamples,plotX,plotArr_p3);
      cpgsci(1);
    }
}

void drawColourMap(dSetStruct *dSet,int pol,plotStruct *plot,float *plotArr_p0,
		   float *plotArr_p1,float *plotArr_p2,float *plotArr_p3)
{
  float tr[6];
  int nchan = dSet->head->nchan;
  int nsblk = dSet->head->nsblk;
  float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
  float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
  float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
  float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
  char title[128];
  
  if (plot->xPlotType==1)
    {
      tr[0] = 0;  tr[1] = 0;  tr[2] = 1;
    }
  else if (plot->xPlotType==2)
    {
      tr[0] = dSet->head->nsblk*plot->t1*dSet->head->tsamp;
      tr[1] = 0;
      tr[2] = dSet->head->tsamp;
    }
  tr[3] = dSet->head->chanFreq[0];  tr[4] = dSet->head->chanbw;  tr[5] = 0;
  if (plot->xPlotType==1)
    {
      //    cpgenv(plot->minx,plot->maxx,plot->miny,plot->maxy,0,1);
      cpgsvp(0.1,0.9,0.15,0.9);
      cpgswin(plot->minx,plot->maxx,plot->miny,plot->maxy);

      //      cpgswin(plot->minx,plot->maxx,0,1000);
      //        cpgswin(plot->minx,plot->maxx,0,plot->maxy);
      cpgbox("ABCTSN",0,0,"ABCTSN",0,0);
    }
  else if (plot->xPlotType==2)
    {
      //      cpgenv(dSet->head->nsblk*plot->t1*dSet->head->tsamp+dSet->head->tsamp*plot->minx,
      //	     dSet->head->nsblk*plot->t1*dSet->head->tsamp+dSet->head->tsamp*plot->maxx,
      //	     plot->miny,plot->maxy,0,1);
      cpgsvp(0.1,0.9,0.15,0.9);
      cpgswin(dSet->head->nsblk*plot->t1*dSet->head->tsamp+dSet->head->tsamp*plot->minx,dSet->head->nsblk*plot->t1*dSet->head->tsamp+dSet->head->tsamp*plot->maxx,plot->miny,plot->maxy);
      cpgbox("ABCTSN",0,0,"ABCTSN",0,0);


    }
  sprintf(title,"%s (subint %.0f-%.0f)",dSet->fileName,plot->t1,plot->t2);
  
  if (plot->xPlotType==1)
    cpglab("Sample","Frequency (MHz)",title);            
  else if (plot->xPlotType==2)
    cpglab("Time from observation start (s)","Frequency (MHz)",title);            

  if (plot->colourPlot==1)
    {
      int i;
      printf("IN HERE pol = %d %d %d %g %g \n",pol,nchan,plot->nTimeSamples,plot->minVal,plot->maxVal);
      cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
      if (pol==0) cpgimag(plotArr_p0,nchan,plot->nTimeSamples,1,nchan,1,plot->nTimeSamples,plot->minVal,plot->maxVal,tr);
      else if (pol==1) cpgimag(plotArr_p1,nchan,plot->nTimeSamples,1,nchan,1,plot->nTimeSamples,plot->minVal,plot->maxVal,tr);
      else if (pol==2) cpgimag(plotArr_p2,nchan,plot->nTimeSamples,1,nchan,1,plot->nTimeSamples,plot->minVal,plot->maxVal,tr);
      else if (pol==3) cpgimag(plotArr_p3,nchan,plot->nTimeSamples,1,nchan,1,plot->nTimeSamples,plot->minVal,plot->maxVal,tr); 
    }
  else
    {
      if (pol==0) cpggray(plotArr_p0,nchan,plot->nTimeSamples,1,nchan,1,plot->nTimeSamples,plot->minVal,plot->maxVal,tr);
      else if (pol==1) cpggray(plotArr_p1,nchan,plot->nTimeSamples,1,nchan,1,plot->nTimeSamples,plot->minVal,plot->maxVal,tr);
      else if (pol==2) cpggray(plotArr_p2,nchan,plot->nTimeSamples,1,nchan,1,plot->nTimeSamples,plot->minVal,plot->maxVal,tr);
      else if (pol==3) cpggray(plotArr_p3,nchan,plot->nTimeSamples,1,nchan,1,plot->nTimeSamples,plot->minVal,plot->maxVal,tr);
    }
}


/*
      // Plot the DM curves
      {
	float fx[nchan],fy[nchan];
	int j,sampleOff;
	double toff,fref;
	int k,nThick=4;
	fref = dSet[0]->head->chanFreq[0];
	
	for (i=0;i<nDM_curve;i++)
	  {
	    for (k=0;k<nThick;k++)
	      {
		for (j=0;j<nchan;j++)
		  {
		    toff = 4.15e-3*dmVal[i]*(pow(fref/1000.0,-2)-pow(dSet[0]->head->chanFreq[j]/1000.0,-2));
		    sampleOff = toff/dSet[0]->head->tsamp;
		    fx[j] =  dmVal_time[i]-sampleOff+k;
		    fy[j] = dSet[0]->head->chanFreq[j];
		  }
		cpgsci(3);
		cpgline(nchan,fx,fy);
		cpgsci(1);
	      }
	  }
      }

*/
