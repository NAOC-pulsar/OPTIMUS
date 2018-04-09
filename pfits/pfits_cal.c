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

//gcc -lm -o pfits_plot pfits_plot.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio 

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"
#include <cpgplot.h>
#include "fitsio.h"
#include "fftw3.h"

void determineOn_Off(dSetStruct *dSet,float *result_p0,float *result_p1,long nSamples);
void measureCalibrationParams(dSetStruct *dSet,long sub0,long sub1,double baseline0,double baseline1,double on0,double on1);
double fortran_mod( double a, double p);


int main(int argc,char *argv[])
{
  dSetStruct *dSet_cal;
  int debug=0;
  int status=0;
  int i,ii;
  long nSamples;
  int nTimeSamples;
  int nFreqSamples;
  float *result_p0,*result_p1;
  long sub0=2; // Ignore the first few subints because of level setting issues
  long sub1;
  double on0,on1;
  double baseline0,baseline1;

  on0 = -115.562;
  on1 = 8.69533;
  baseline0 = 72.336;
  baseline1 = 190.531;
   

  // Initialise everything
  initialise(&dSet_cal,debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-c")==0)
	setFilename(argv[++i],dSet_cal,debug);
    }

  pfitsOpenFile(dSet_cal,debug);
  pfitsLoadHeader(dSet_cal,debug);

  sub1 = dSet_cal->head->nsub;
  result_p0 = (float *)malloc(sizeof(float)*dSet_cal->head->nsblk*dSet_cal->head->nsub);
  result_p1 = (float *)malloc(sizeof(float)*dSet_cal->head->nsblk*dSet_cal->head->nsub);

        pfits_read1pol_zeroDM_float(result_p0,0,dSet_cal,sub0,sub1,2,&nSamples,&nTimeSamples,&nFreqSamples,debug);
        pfits_read1pol_zeroDM_float(result_p1,1,dSet_cal,sub0,sub1,2,&nSamples,&nTimeSamples,&nFreqSamples,debug);
    //  for (i=0;i<nSamples;i++)
  //    printf("result: %ld %g\n",i,result_p0[i]);

  // Should return on0, on1, baseline0 and baseline1 from this function
          determineOn_Off(dSet_cal,result_p0,result_p1,nSamples);

  measureCalibrationParams(dSet_cal,sub0,sub1,baseline0,baseline1,on0,on1);

  
  // De-allocate the memory
  deallocateMemory(&dSet_cal,debug);
  free(result_p0);
  free(result_p1);
}

void measureCalibrationParams(dSetStruct *dSet,long sub0,long sub1,double baseline0,double baseline1,double on0,double on1)
{
  long sub;
  long snum=0;
  int s;
  double check;
  int type=0;
  int status=0;
  int colnum;
  int nchan = dSet->head->nchan;
  int npol = dSet->head->npol;
  int samplesperbyte = 8/dSet->head->nbits;
  int initflag=0;
  unsigned char *cVals;
  unsigned char nval = '0';
  float *floatData;
  float s1;
  double pulseFreq = 11.123;
  double pulsePeriodSamples = (1.0/pulseFreq)/dSet->head->tsamp;
  double pol0_on[nchan];
  double pol0_off[nchan];
  double pol1_on[nchan];
  double pol1_off[nchan];
  double pol2_on[nchan];
  double pol2_off[nchan];
  double pol3_on[nchan];
  double pol3_off[nchan];
  float plotX[nchan];
  float pol0_h[nchan];
  float pol1_h[nchan];
  float pol2_h[nchan];
  float pol3_h[nchan];
  float py0[nchan];
  float py1[nchan];
  float py2[nchan];
  float py3[nchan];
  float minx,maxx;

  int i;
  long nOn=0;
  long nOff=0;
  float miny,maxy;
  float maxy0,maxy1,miny1,miny0,maxy2,miny2;

  int plot=1;
  float mx,my;
  char key;
  int t=1;

  for (i=0;i<nchan;i++)
    {
      pol0_on[i]=0; pol0_off[i]=0;
      pol1_on[i]=0; pol1_off[i]=0;
      pol2_on[i]=0; pol2_off[i]=0;
      pol3_on[i]=0; pol3_off[i]=0;
    }

  floatData = (float *)malloc(sizeof(float)*nchan*npol);
  cVals = (unsigned char *)malloc(sizeof(unsigned char)*nchan*npol/samplesperbyte);

  fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
  if (status){printf("Unable to find the SUBINT table\n"); exit(1);}
  fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum,&status);


  for (sub=sub0;sub<sub1;sub++)
    {
      for (s=0;s<dSet->head->nsblk;s++)
	{
	  s1 = snum-((int)(snum/pulsePeriodSamples+0.5))*pulsePeriodSamples;
	  if (s1 >= on0 && s1 < on1)
	    type=1;
	  else if (s1 >= baseline0 && s1 < baseline1)
	    type=2;
	  else
	    type=3;

	  if (type==1 || type==2)
	    {
	      fits_read_col_byt(dSet->fp,colnum,
				sub+1,
				s*nchan*npol+1,
				nchan*npol/samplesperbyte,
				nval,cVals,&initflag,&status);
	      pfits_bytesToFloats(samplesperbyte,nchan*npol,cVals,floatData);	      
	      for (i=0;i<nchan;i++)
		{
		  if (type==1)
		    {
		      pol0_on[i]+=floatData[i];
		      pol1_on[i]+=floatData[nchan+i];
		      pol2_on[i]+=floatData[2*nchan+i];
		      pol3_on[i]+=floatData[3*nchan+i];
		      //		      printf("Loaded %g %g %g %g\n",pol0_on[i],pol1_on[i],pol2_on[i],pol3_on[i]);
		      nOn++;
		    }
		  else if (type==2)
		    {
		      pol0_off[i]+=floatData[i];
		      pol1_off[i]+=floatData[nchan+i];
		      pol2_off[i]+=floatData[2*nchan+i];
		      pol3_off[i]+=floatData[3*nchan+i];
		      nOff++;
		    }
		}
	    }

	  snum++;
	}
    }
  printf("Channel 100: %g %g\n",pol0_on[100]/(double)nOn,pol0_off[100]/(double)nOff);
  for (i=0;i<nchan;i++)
    {
      plotX[i] = i;
      //      printf("Before %g\n",pol0_on[i]);
      pol0_on[i]/=(double)nOn;
      //      if (i==0)
      //	miny = maxy = pol0_on[i];
      //      printf("Result: %g %d\n",pol0_on[i],nOn);
      pol1_on[i]/=(double)nOn;
      pol2_on[i]/=(double)nOn;
      pol3_on[i]/=(double)nOn;

      pol0_off[i]/=(double)nOff;
      //      printf("offVal = %g\n",pol0_off[i]);
      pol1_off[i]/=(double)nOff;
      pol2_off[i]/=(double)nOff;
      pol3_off[i]/=(double)nOff;

      pol0_h[i] = pol0_on[i]-pol0_off[i];
      pol1_h[i] = pol1_on[i]-pol1_off[i];
      pol2_h[i] = pol2_on[i]-pol2_off[i];
      pol3_h[i] = pol3_on[i]-pol3_off[i];

      /*      if (miny > pol0_off[i]) miny = pol0_off[i];
      if (miny > pol1_off[i]) miny = pol1_off[i];
      if (miny > pol2_off[i]) miny = pol2_off[i];
      if (miny > pol3_off[i]) miny = pol3_off[i];

      if (maxy < pol0_on[i]) maxy = pol0_on[i];
      if (maxy < pol1_on[i]) maxy = pol1_on[i];
      if (maxy < pol2_on[i]) maxy = pol2_on[i];
      if (maxy < pol3_on[i]) maxy = pol3_on[i]; */
      if (i==0)
	miny = maxy = pol0_h[i];

      if (miny > pol0_h[i]) miny = pol0_h[i];
      if (miny > pol1_h[i]) miny = pol1_h[i];
      if (miny > pol2_h[i]) miny = pol2_h[i];
      if (miny > pol3_h[i]) miny = pol3_h[i];

      if (maxy < pol0_h[i]) maxy = pol0_h[i];
      if (maxy < pol1_h[i]) maxy = pol1_h[i];
      if (maxy < pol2_h[i]) maxy = pol2_h[i];
      if (maxy < pol3_h[i]) maxy = pol3_h[i]; 
      
      
    }
    
  printf("Loaded %d on and %d off\n",nOn,nOff);
  minx = 0;
  maxx = nchan;

  do {
    t=1;
    if (plot==1 || plot==2)
      {
 	cpgbeg(0,"/xs",1,1);
	cpgask(0);
	cpgsch(1.4); cpgsfs(2); cpgslw(2);    
      }
    else if (plot==3)
      {
	cpgbeg(0,"/xs",1,3);
	cpgask(0);
	cpgsch(1.4); cpgsfs(2); cpgslw(2);    
      }
    cpgsci(1);

    for (i=0;i<nchan;i++)
      {
	if (plot==1)
	  {
	    py0[i] = pol0_h[i];
	    py1[i] = pol1_h[i];
	    py2[i] = pol2_h[i];
	    py3[i] = pol3_h[i];
	  }
	else if (plot==2)
	  {
	    py0[i] = pol0_h[i]+pol1_h[i]; // This needs checking
	    py1[i] = pol0_h[i]-pol1_h[i];
	    py2[i] = 2*pol2_h[i];
	    py3[i] = 2*pol3_h[i];
	  }
	else if (plot==3)
	  {
	    double stokesI,stokesQ,stokesU,stokesV;
	    double i0,gain,phase;

	    stokesI = pol0_h[i]+pol1_h[i]; // This needs checking
	    stokesQ = pol0_h[i]-pol1_h[i];
	    stokesU = 2*pol2_h[i];
	    stokesV = 2*pol3_h[i];

	    py0[i] = i0 = stokesI;
	    py1[i] = gain = 2*stokesQ/stokesI;
	    py2[i] = phase = atan2(stokesV,stokesU)*180.0/M_PI;
	  }

	if (plot==1 || plot==2)
	  {
	    if (plotX[i] >= minx && plotX[i] <= maxx)
	      {
		if (t==1)
		  {
		    miny = maxy = py0[i];
		    t=2;
		  }
		if (miny > py0[i]) miny = py0[i];
		if (miny > py1[i]) miny = py1[i];
		if (miny > py2[i]) miny = py2[i];
		if (miny > py3[i]) miny = py3[i];
		
		if (maxy < py0[i]) maxy = py0[i];
		if (maxy < py1[i]) maxy = py1[i];
		if (maxy < py2[i]) maxy = py2[i];
		if (maxy < py3[i]) maxy = py3[i];	    
	      }
	  }
	else if (plot==3)
	  {
	    if (plotX[i] >= minx && plotX[i] <= maxx)
	      {
		if (t==1)
		  {
		    maxy0=py0[i];
		    miny0=py0[i];
		    maxy1=py1[i];
		    miny1=py1[i];
		    maxy2=py2[i];
		    miny2=py2[i];
		    t=2;
		  }
		if (py0[i] > maxy0) maxy0=py0[i];
		if (py0[i] < miny0) miny0=py0[i];
		if (py1[i] > maxy1) maxy1=py1[i];
		if (py1[i] < miny1) miny1=py1[i];
		if (py2[i] > maxy2) maxy2=py2[i];
		if (py2[i] < miny2) miny2=py2[i];
	      }
	  }
      }
    printf("maxy = %g %g %g\n",maxy0,maxy1);
    if (plot==1 || plot==2)
      {
	cpgenv(minx,maxx,miny,maxy,0,1);
	cpglab("Frequency channel","","");
	cpgline(nchan,plotX,py0);
	cpgsci(2); cpgline(nchan,plotX,py1);
	cpgsci(3); cpgline(nchan,plotX,py2);
	cpgsci(4); cpgline(nchan,plotX,py3);
      } 
    else if (plot==3)
      {
	cpgenv(minx,maxx,miny0,maxy0,0,1);
	cpgline(nchan,plotX,py0);
	cpgenv(minx,maxx,miny1,maxy1,0,1);
	cpgline(nchan,plotX,py1);
	cpgenv(minx,maxx,miny2,maxy2,0,1);
	cpgline(nchan,plotX,py2);
      }
    cpgcurs(&mx,&my,&key);
    if (key=='1')
      plot=1;
    else if (key=='2')
      plot=2;
    else if (key=='3')
      plot=3;
    else if (key=='z')
      {
	float mx2,my2;

	cpgband(4,0,mx,my,&mx2,&my2,&key);
	minx = mx;
	maxx = mx2;
      }
    else if (key=='u')
      {
	minx = 0;
	maxx = nchan;
      }
    cpgend();
  } while (key!='q');


  free(floatData);
  free(cVals);

}

void determineOn_Off(dSetStruct *dSet,float *result_p0,float *result_p1,long nSamples)
{
  float *plotX;
  float *plotY;
  float *plotOnX;
  float *plotOnY;
  float *plotBaselineX;
  float *plotBaselineY;
  int i;
  float miny,maxy;
  float ominy,omaxy;
  float minx,maxx;
  float mx,my;
  double baseline0=-1,baseline1=-1;
  double on0=-1,on1=-1;
  int setOn=0;
  int setBaseline=0;
  int nOnSamples,nBaselineSamples;
  double pulseFreq = 11.123;
  double pulsePeriodSamples = (1.0/pulseFreq)/dSet->head->tsamp;
  double s1,s2;
  char key;
  double meanOn,meanOff;


  plotX = (float *)malloc(sizeof(float)*nSamples);
  plotY = (float *)malloc(sizeof(float)*nSamples);
  plotOnX = (float *)malloc(sizeof(float)*nSamples);
  plotOnY = (float *)malloc(sizeof(float)*nSamples);
  plotBaselineX = (float *)malloc(sizeof(float)*nSamples);
  plotBaselineY = (float *)malloc(sizeof(float)*nSamples);
  for (i=0;i<nSamples;i++)
    {
      plotX[i] = i;
      plotY[i] = result_p0[i]+result_p1[i];
      if (i==0)
	  miny = maxy = plotY[i];
      else
	{
	  if (miny > plotY[i]) miny = plotY[i];
	  if (maxy < plotY[i]) maxy = plotY[i];
	}
    }
  ominy = miny; 
  omaxy = maxy;
  minx = 0;
  maxx = nSamples;

  cpgbeg(0,"/xs",1,1);
  cpgask(0);
  cpgsch(1.4); cpgslw(2); cpgscf(2);
  do {
    cpgenv(minx,maxx,miny,maxy,0,1);
    cpgline(nSamples,plotX,plotY);

    if (setOn == 1)
      {
	nOnSamples=0;
	meanOn=0;
	for (i=0;i<nSamples;i++)
	  {
	    s1 = i-((int)(i/pulsePeriodSamples+0.5))*pulsePeriodSamples;
	    if (s1 >= on0 && s1 < on1)
	      {
		plotOnX[nOnSamples] = i;
		plotOnY[nOnSamples] = (miny+maxy)/2.0+(maxy-miny)*0.05;
		meanOn+=plotY[i];
		nOnSamples++;
	      }
	  }
	cpgsci(2); cpgpt(nOnSamples,plotOnX,plotOnY,20); cpgsci(1);
	printf("Mean ON value is %g\n",meanOn/(double)nOnSamples);
	printf("On samples are %g %g\n",on0,on1);
      }

    if (setBaseline == 1)
      {
	meanOff=0;
	nBaselineSamples=0;
	for (i=0;i<nSamples;i++)
	  {
	    s1 = i-((int)(i/pulsePeriodSamples+0.5))*pulsePeriodSamples;
	    if (s1 >= baseline0 && s1 < baseline1)
	      {
		plotBaselineX[nBaselineSamples] = i;
		plotBaselineY[nBaselineSamples] = (miny+maxy)/2.0-(maxy-miny)*0.05;
		meanOff+=plotY[i];
		nBaselineSamples++;
	      }
	  }
	cpgsci(4); cpgpt(nBaselineSamples,plotBaselineX,plotBaselineY,20); cpgsci(1);
	printf("Mean OFF value is %g\n",meanOff/(double)nBaselineSamples);
	printf("Baseline samples are %g %g\n",baseline0,baseline1);
      }
    

    cpgcurs(&mx,&my,&key);
    if (key=='z')
      {
	float mx2,my2;
	cpgband(2,0,mx,my,&mx2,&my2,&key);
	minx = mx;
	maxx = mx2;
	miny = my;
	maxy = my2;
      }
    else if (key=='u')
      {
	minx = 0;
	maxx = nSamples;
	miny = ominy;
	maxy = omaxy;
      }
    else if (key=='b')
      {
	float mx2,my2;
	setBaseline=1;
	cpgband(4,0,mx,my,&mx2,&my2,&key);
	baseline0 = mx-((int)(mx/pulsePeriodSamples+0.5))*pulsePeriodSamples;
	baseline1 = mx2-((int)(mx/pulsePeriodSamples+0.5))*pulsePeriodSamples;
	printf("Baseline %g %g\n",baseline0,baseline1);
      }
    else if (key=='o')
      {
	float mx2,my2;
	setOn=1;
	cpgband(4,0,mx,my,&mx2,&my2,&key);
	on0 = mx-((int)(mx/pulsePeriodSamples+0.5))*pulsePeriodSamples;
	on1 = mx2-((int)(mx/pulsePeriodSamples+0.5))*pulsePeriodSamples;
	//        On0 = fortran_mod(on0/pulsePeriodSamples,1);
	//	phaseOn1 = fortran_mod(on1/pulsePeriodSamples,1);
	printf("On %g %g\n",on0,on1);
      }
  } while (key!='q');
  cpgend();
  free(plotX);
  free(plotY);
  free(plotOnX);
  free(plotOnY);
  free(plotBaselineX);
  free(plotBaselineY);
}

double fortran_mod( double a, double p)
{ 
  double ret;
  ret = a - (int)(a/p)*p;
  return ret;
}
