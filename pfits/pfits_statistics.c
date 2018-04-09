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

//gcc -lm -o pfits_statistcs pfits_statistics.c pfits_statistics.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio 

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include "fitsio.h"

void doPlot(float *sum,float *min,float *max,int totCount,int *histogramTot,float *meanSub,int nchan,dSetStruct *dSet,int nsub);

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  int debug=0;
  int status=0;
  int i,j;
  float *data;
  int subint;
  long nSamples;
  int nTime,nFreq;
  float *sum,*sumSub,*max,*min;
  float *meanSub;
  long totCount;
  int sub0,sub1;
  int pol=0;
  int histogramTot[255];

  for (i=0;i<255;i++)
    histogramTot[i]=0;
 
  sub0 = 0;
  sub1 = 1;
  
  // Initialise everything
  initialise(&dSet,debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-s0")==0)
	sscanf(argv[++i],"%d",&sub0);
      else if (strcmp(argv[i],"-s1")==0)
	sscanf(argv[++i],"%d",&sub1);
      else if (strcmp(argv[i],"-pol")==0)
	sscanf(argv[++i],"%d",&pol);
    }

  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);

  data = (float *)malloc(sizeof(float)*dSet->head->nsblk*dSet->head->nchan);
  sum = (float *)malloc(sizeof(float)*dSet->head->nchan);
  sumSub = (float *)malloc(sizeof(float)*dSet->head->nchan);
  max = (float *)malloc(sizeof(float)*dSet->head->nchan);
  min = (float *)malloc(sizeof(float)*dSet->head->nchan);
  meanSub = (float *)malloc(sizeof(float)*dSet->head->nchan*((sub1-sub0)+1));
  for (i=0;i<dSet->head->nchan;i++)
    sum[i]=sumSub[i]=0.0;
  totCount=0;

  for (subint=sub0;subint<sub1;subint++)
    {
      pfits_read1pol_float(data,pol,dSet,subint,subint,1,&nSamples,&nTime,&nFreq,debug);
      for (i=0;i<dSet->head->nchan;i++)
	sumSub[i]=0.0;
      for (j=0;j<nTime;j++)
	{
	  for (i=0;i<dSet->head->nchan;i++)
	    {
	      //	      printf("Setting %d\n",(int)data[j*nFreq+i]);
	      histogramTot[(int)data[j*nFreq+i]]++;
	      sum[i]+=data[j*nFreq+i];
	      sumSub[i]+=data[j*nFreq+i];
	      if (totCount==0)
		{
		  min[i] = data[j*nFreq+i];
		  max[i] = data[j*nFreq+i];
		}
	      else
		{
		  if (min[i] > data[j*nFreq+i])
		    min[i] = data[j*nFreq+i];
		  if (max[i] < data[j*nFreq+i])
		    max[i] = data[j*nFreq+i];
		}
	    }
		
	  totCount++;
	}
      for (i=0;i<dSet->head->nchan;i++)
	meanSub[i+(subint-sub0)*nFreq] = sumSub[i]/(float)nTime;

    }
  doPlot(sum,min,max,totCount,histogramTot,meanSub,dSet->head->nchan,dSet,(sub1-sub0));
  free(data);
  free(max);
  free(min);
  free(sum);
  free(sumSub);
  free(meanSub);
  // De-allocate the memory
  deallocateMemory(&dSet,debug);
}

void doPlot(float *sum,float *min,float *max,int totCount,int *histogramTot,float *meanSub,int nchan,dSetStruct *dSet,int nsub)
{
  float xVal[nchan];
  float yVal[nchan];
  float yVal2[nchan];
  float yVal3[nchan];
  float histPlotX[255];
  float histPlotY[255];
  int i,j;
  float maxy;
  float miny;
  float maxHist=0;
  int plot=1;
  float mx,my;
  char key;
  
  for (i=0;i<255;i++)
    {
      histPlotX[i] = i;
      histPlotY[i] = histogramTot[i];
      if (maxHist < histPlotY[i]) maxHist = histPlotY[i];
    }
  
  for (i=0;i<nchan;i++)
    {
      xVal[i] = i;
      yVal[i] = sum[i]/(float)totCount;
      yVal2[i] = min[i];
      yVal3[i] = max[i];
      if (i==0)
	{
	  miny = min[i];
	  maxy = max[i];
	}
      else
	{
	  if (miny > min[i]) miny = min[i];
	  if (maxy < max[i]) maxy = max[i];
	}
    }
  cpgbeg(0,"/xs",1,1);
  cpgask(0);
  cpgsch(1.4);
  cpgslw(2);
  cpgscf(2);
  do {
    if (plot==1)
      {
	cpgenv(0,nchan,miny,maxy,0,1);
	cpglab("Frequency channel","Signal","");
	cpgsci(1); cpgbin(nchan,xVal,yVal,0);
	cpgsci(2); cpgbin(nchan,xVal,yVal2,0);
	cpgsci(3); cpgbin(nchan,xVal,yVal3,0);
	cpgsci(1);
      }
    else if (plot==2)
      {
	if (dSet->head->nbits==8)
	  {
	    cpgenv(0,255,0,maxHist+0.1*maxHist,0,1);
	    cpglab("Value","Count","");
	    cpgbin(255,histPlotX,histPlotY,0);
	  }
	else if (dSet->head->nbits==4)
	  {
	    cpgenv(0,20,0,maxHist+0.1*maxHist,0,1);
	    cpglab("Value","Count","");
	    cpgbin(21,histPlotX,histPlotY,0);
	  }
	else if (dSet->head->nbits==2)
	  {
	    cpgenv(0,6,0,maxHist+0.1*maxHist,0,1);
	    cpglab("Value","Count","");
	    cpgbin(7,histPlotX,histPlotY,0);
	  }
	else if (dSet->head->nbits==1)
	  {
	    cpgenv(0,3,0,maxHist+0.1*maxHist,0,1);
	    cpglab("Value","Count","");
	    cpgbin(4,histPlotX,histPlotY,0);
	  }
      }
    else if (plot==3)
      {
	float tr[6];
	float maxV=0;
	float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
	float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
	float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
	float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};

	for (j=0;j<nchan*nsub;j++)
	  {
	    if (maxV < meanSub[j]) maxV=meanSub[j];
	  }
	
	
	tr[0] = 0; tr[1] = 1; tr[2] = 0;
	tr[3] = 0; tr[4] = 0; tr[5] = 1;
	cpgenv(0,nchan,0,nsub,0,1);

	cpglab("Frequency channel","Subintegration number","");
	cpgctab(heat_l,heat_r,heat_g,heat_b,5,1.0,0.5);
	//       	cpgimag(meanSub,nchan,nsub,1,nchan,1,nsub,0,maxV,tr);
	cpggray(meanSub,nchan,nsub,1,nchan,1,nsub,0,maxV,tr);
	cpgaxis("N", 0, nsub, nchan, nsub, dSet->head->chanFreq[0],dSet->head->chanFreq[nchan-1],0.0,0,0.5,0,0.5,-0.5,0.0);
	cpgaxis("N", nchan, 0, nchan, nsub, 0,nsub*dSet->head->nsblk*dSet->head->tsamp,0.0,0,0.5,0,0.5,0.5,0.0);

      }
    cpgcurs(&mx,&my,&key);
    if (key=='1')
      plot=1;
    else if (key=='2')
      plot=2;
    else if (key=='3')
      plot=3;
  } while (key!='q');
  cpgend();

}
