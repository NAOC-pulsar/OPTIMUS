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

// Standard routines for pfits

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pfits.h"
#include "fitsio.h"

/*
 * routine to load 1 polarisation data from the PSRFITS file and return an array of floats
 * note that the PSRFITS format is (nchan,npol,nsblk*nbits/8)
 */
void pfits_read1pol_float(float *out,int polNum,dSetStruct *dSet,float t1,float t2,int rangeType,long *nSamples,int *nTimeSamples,int *nFreqSamples,int debugFlag)
{
  int nchan;
  int nbits;
  int npol;
  int nsblk;
  int samplesperbyte;
  int i,j,i0,i1;
  long ipolpos=0;
  int status=0;
  int colnum;
  unsigned char *cVals;
  unsigned char nval = '0';
  unsigned int *iVals;
  short int *sVals;
  short int n_sval = 0;
  float *fVals;
  float n_fval=0;
  unsigned int n_ival = 0;
  int subint=1;
  int initflag=0;
  long s0,s1;
  long firstSamp,lastSamp;
  float tsamp;
  int scount=0;
  
  nchan = dSet->head->nchan;
  npol  = dSet->head->npol;
  nsblk = dSet->head->nsblk;
  nbits = dSet->head->nbits;
  tsamp = dSet->head->tsamp;
  
  if (rangeType==1)
    {
      s0 = (long)t1;
      s1 = (long)t1;
      firstSamp = 1;
      lastSamp = dSet->head->nsblk;
    }
  else if (rangeType==2)
    {
      s0 = (long)t1;
      s1 = (long)t2;
      firstSamp = 1;
      lastSamp = dSet->head->nsblk;
    }
  else if (rangeType==3)
    {
      s0 = (long)(int)(t1/(tsamp*nsblk));
      s1 = (long)(int)(t2/(tsamp*nsblk));
      firstSamp = (int)((t1-s0*tsamp*nsblk)/tsamp);
      lastSamp = (int)((t2-s1*tsamp*nsblk)/tsamp+0.5);
    }
  if (debugFlag==1) printf("Entering pfits_read1pol_float\n");

  if (nbits <= 8)
    samplesperbyte = 8/nbits;    
  else
    samplesperbyte = 1;
  
  if (dSet->headerMemorySet==1 && dSet->fileOpen==1)
    {
      fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
      if (status){printf("Unable to find the SUBINT table\n"); exit(1);}
      fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum,&status);

      // Allocate enough memory for one group of channels
      if (nbits <= 8)
	cVals = (unsigned char *)malloc(sizeof(unsigned char)*nchan/samplesperbyte); // Note:       reading one polarisation
      else if (nbits==16)
	sVals = (short int *)malloc(sizeof(short int)*nchan);
      else
	fVals = (float *)malloc(sizeof(float)*nchan/samplesperbyte); // Note: only reading one polarisation
	//	iVals = (unsigned int *)malloc(sizeof(unsigned int)*nchan/samplesperbyte); // Note: only reading one polarisation
      for (subint=s0;subint<=s1;subint++)
	{
	  printf("Reading subint %d\n",subint);
	  // Read in 1 subint of data for the specified (polNum) polarisation
	  if (subint==s0)
	    i0 = firstSamp;
	  else
	    i0 = 0;

	  if (subint==s1)
	    i1 = lastSamp;
	  else
	    i1 = nsblk;
	  
	  for (i=i0;i<i1;i++)
	    {
	      if (nbits <= 8)
		{
		  fits_read_col_byt(dSet->fp,colnum,
				    subint+1,
				    1+i*npol*nchan/samplesperbyte+polNum*nchan/samplesperbyte,
				    nchan/samplesperbyte,
				    nval,cVals,&initflag,&status);
		  pfits_bytesToFloats(samplesperbyte,nchan,cVals,out+scount*nchan);
		}
	      else if (nbits == 16)
		{
		  fits_read_col(dSet->fp,TSHORT,colnum,
				subint+1,
				1+i*npol*nchan/samplesperbyte+polNum*nchan/samplesperbyte,
				nchan/samplesperbyte,
				&n_sval,sVals,&initflag,&status);
		  for (j=0;j<nchan;j++)  // Is this really npol??
		    {
		      out[scount*nchan+j]=sVals[j];
		      printf("loaded %d %d %d (%d %d %d %d)\n",scount,j,sVals[j],npol,nchan,samplesperbyte,polNum);
		    }
		}
	      else
		{
		  /*
		    fits_read_col(dSet->fp,TINT,colnum,
		    subint+1,
		    1+i*npol*nchan/samplesperbyte+polNum*nchan/samplesperbyte,
		    nchan/samplesperbyte,
		    &n_ival,iVals,&initflag,&status);
		  */
		  printf("In here\n");
		  fits_read_col(dSet->fp,TFLOAT,colnum,
				subint+1,
				1+i*npol*nchan/samplesperbyte+polNum*nchan/samplesperbyte,
				nchan/samplesperbyte,
				&n_fval,fVals,&initflag,&status);

		  //		  for (j=0;j<npol*nchan;j++)  // Is this really npol??
		  for (j=0;j<nchan;j++)  // Is this really npol??
		    out[scount*nchan+j]=fVals[j];
		    //		    out[scount*nchan+j]=iVals[j];
		  printf("Got here\n");
		}
	      scount++;
	    }
	}
      // Deallocate the memory
      if (nbits <= 8)
	free(cVals);
      else if (nbits==16)
	free(sVals);
      else
	free(fVals);
	//	free(iVals);
    }
  else
    {
      printf("Unable to read file\n");
      exit(1);
    }
  
  if (debugFlag==1) printf("Leaving pfits_read1pol_float\n"); 

  *nSamples = scount*nchan;
  *nFreqSamples = nchan;
  *nTimeSamples = scount;
  printf("Loaded %d\n",*nTimeSamples);
}

/*
 * routine to load 1 polarisation data from the PSRFITS file and return an array of floats summed over frequency
 * note that the PSRFITS format is (nchan,npol,nsblk*nbits/8)
 */
void pfits_read1pol_zeroDM_float(float *out,int polNum,dSetStruct *dSet,float t1,float t2,int rangeType,long *nSamples,int *nTimeSamples,int *nFreqSamples,int debugFlag)
{
  int nchan;
  int nbits;
  int npol;
  int nsblk;
  int samplesperbyte;
  int i,j,i0,i1;
  long ipolpos=0;
  int status=0;
  int colnum;
  unsigned char *cVals;
  unsigned char nval = '0';
  unsigned int *iVals;
  unsigned int n_ival = 0;
  int subint=1;
  int initflag=0;
  long s0,s1;
  long firstSamp,lastSamp;
  float tsamp;
  float *loadVals;
  long scount=0;
  
  nchan = dSet->head->nchan;
  npol  = dSet->head->npol;
  nsblk = dSet->head->nsblk;
  nbits = dSet->head->nbits;
  tsamp = dSet->head->tsamp;
  
  loadVals = (float *)malloc(sizeof(float)*nchan);

  if (rangeType==1)
    {
      s0 = (long)t1;
      s1 = (long)t1;
      firstSamp = 1;
      lastSamp = dSet->head->nsblk;
    }
  else if (rangeType==2)
    {
      s0 = (long)t1;
      s1 = (long)t2;
      firstSamp = 1;
      lastSamp = dSet->head->nsblk;
    }
  else if (rangeType==3)
    {
      s0 = (long)(int)(t1/(tsamp*nsblk));
      s1 = (long)(int)(t2/(tsamp*nsblk));
      firstSamp = (int)((t1-s0*tsamp*nsblk)/tsamp);
      lastSamp = (int)((t2-s1*tsamp*nsblk)/tsamp+0.5);
    }
  if (debugFlag==1) printf("Entering pfits_read1pol_float\n");

  if (nbits <= 8)
    samplesperbyte = 8/nbits;    
  else
    samplesperbyte = 1;
  
  if (dSet->headerMemorySet==1 && dSet->fileOpen==1)
    {
      fits_movnam_hdu(dSet->fp,BINARY_TBL,"SUBINT",1,&status);
      if (status){printf("Unable to find the SUBINT table\n"); exit(1);}
      fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum,&status);

      // Allocate enough memory for one group of channels
      if (nbits <= 8)
	cVals = (unsigned char *)malloc(sizeof(unsigned char)*nchan/samplesperbyte); // Note:       reading one polarisation
      else 
	iVals = (unsigned int *)malloc(sizeof(unsigned int)*nchan/samplesperbyte); // Note: only reading one polarisation
      for (subint=s0;subint<=s1;subint++)
	{
	  printf("Reading subint %d\n",subint);
	  // Read in 1 subint of data for the specified (polNum) polarisation
	  if (subint==s0)
	    i0 = firstSamp;
	  else
	    i0 = 0;

	  if (subint==s1)
	    i1 = lastSamp;
	  else
	    i1 = nsblk;
	  
	  for (i=i0;i<i1;i++)
	    {
	      if (nbits <= 8)
		{
		  fits_read_col_byt(dSet->fp,colnum,
				    subint+1,
				    1+i*npol*nchan/samplesperbyte+polNum*nchan/samplesperbyte,
				    nchan/samplesperbyte,
				    nval,cVals,&initflag,&status);
		  pfits_bytesToFloats(samplesperbyte,nchan,cVals,loadVals);
		  //		  for (j=0;j<nchan;j++)
		  //		    printf("chan: %d %f\n",j,loadVals[j]);
		  //		  exit(1);
		  out[scount]=loadVals[0];
		  for (j=1;j<nchan;j++)
		    out[scount]+=loadVals[j];
		  scount++;
		}
	      else
		{
		  fits_read_col(dSet->fp,TINT,colnum,
				subint+1,
				1+i*npol*nchan/samplesperbyte+polNum*nchan/samplesperbyte,
				nchan/samplesperbyte,
				&n_ival,iVals,&initflag,&status);
		  out[scount] = iVals[0];
		  for (j=1;j<nchan;j++)
		    out[scount]+=iVals[j];
		  scount++;
		}
	    }
	}
      // Deallocate the memory
      if (nbits <= 8)
	free(cVals);
      else
	free(iVals);
    }
  else
    {
      printf("Unable to read file\n");
      exit(1);
    }
  
  if (debugFlag==1) printf("Leaving pfits_read1pol_float\n"); 

  *nSamples = scount;
  *nFreqSamples = 1;
  *nTimeSamples = scount;
  printf("Loaded %d\n",*nTimeSamples);

  free(loadVals);
}


void pfits_bytesToFloats(int samplesperbyte,int n,unsigned char *cVals,float *out)
{
  int i,j;
  int pos=0;
  int index = 0;
  for (i = 0; i < n/samplesperbyte; i++)
    {
      void (*pointerFloat)(int, float *, int*);
      switch (samplesperbyte)
	{
	case 1:
	  pointerFloat = eightBitsFloat;
	  break;

	case 2:
	  pointerFloat = fourBitsFloat;
	  break;

	case 4:
	  pointerFloat = twoBitsFloat;
	  break;

	case 8:
	  pointerFloat = oneBitFloat;
	  break;
	}
      pointerFloat(cVals[i], out, &index);
    }
}


void eightBitsFloat(int eight_bit_number, float *results, int *index)
{
  // 135 subtracted from the original number to make the range from -8 to something
  // 0.5 is then added to remove bias

  //results[(*index)++] = eight_bit_number - 135.5;
  results[(*index)++] = eight_bit_number; // -31.5; // -31.5;
}

void fourBitsFloat(int eight_bit_number, float *results, int *index)
{
  // anding the least significant 4 bits with 0000 1111 will give the (signed) 4-bit number
  // shifting right 4 bits will produce the other (first) 4-bit numbers
  // 0.5 is added to each number to compensate for the bias, making the range -7.5 -> +7.5

  float tempResults[2];
  int i;
  for (i = 0; i < 2; i++) {
    int andedNumber = eight_bit_number & 15;
    tempResults[i] = andedNumber; // - 7.5;
    eight_bit_number = eight_bit_number >> 4;
  }

  for (i = 1; i >= 0; i--)
    results[(*index)++] = tempResults[i];
}

void twoBitsFloat(int eight_bit_number, float *results, int *index)
{
  // anding the least significant 2 bits with 0000 0011 will give the (signed 2-bit number
  // shifting right 2 bits will produce the next (previous) 2-bit number
  // the numbers are adjusted to compensate for bias and to give a rms of ~0.9
  //
  // 1 -> 2.0
  // 0 -> +0.5
  // -1 -> -0.5
  // -2 -> -2.0

  float tempResults[4];
  int i;

  for (i = 0; i < 4; i++) {
    int andedNumber = eight_bit_number & 3;
    switch (andedNumber) {
    case 0:
      tempResults[i] = 0; //-2.5;
      break;
    case 1:
      tempResults[i] = 1; // -0.5;
      break;
    case 2:
      tempResults[i] = 2; //0.5;
      break;
    case 3:
      tempResults[i] = 3; // 2.5;
      break;
    }
    eight_bit_number = eight_bit_number >> 2;
  }

  for (i = 3; i >= 0; i--)
    {
      results[(*index)++] = tempResults[i];
      //      printf("ret: %d %g\n",(*index),tempResults[i]);
    }
}

void oneBitFloat(int eight_bit_number, float *results, int *index)
{
  // anding the least significant bit with 0000 0001 will give the (signed 1-bit number
  // shifting right 1 bit will produce the next (previous) 1-bit number
  //
  // 0.5 is added to each number to compensate for bias
  //

  float tempResults[8];
  int i;

  for (i = 0; i < 8; i++) {
    int andedNumber = eight_bit_number & 1;
    // tempResults[i] = andedNumber ? 0.5 : -0.5;
    tempResults[i] = andedNumber ? 0 : 1;
    eight_bit_number = eight_bit_number >> 1;
  }

  for (i = 7; i >= 0; i--)
    {
      results[(*index)++] = tempResults[i];
      //            printf("This pt: %d %g\n",*index,tempResults[i]);
    }
}


void calibrateScalePols(calibrateStruct *cal,float *p1,float *p2,float *p3,float *p4,int n)
{
  int i;
  float c1,c2,c3,c4;
  for (i=0;i<n;i++)
    {
      p1[i]-=cal->baseline_p0[0];
      p2[i]-=cal->baseline_p1[0];
      p3[i]-=cal->baseline_p2[0];
      p4[i]-=cal->baseline_p3[0];
      convertStokes(p1[i],p2[i],p3[i],p4[i],&c1,&c2,&c3,&c4);
      p1[i] = c1;
      p2[i] = c2;
      p3[i] = c3;
      p4[i] = c4; 
    }
}

void convertStokes(float p1,float p2,float p3,float p4,float *stokesI,float *stokesQ,float *stokesU,float *stokesV)
{
  *stokesI = p1+p2;
  *stokesQ = p1-p2;
  *stokesU = 2*p3;
  *stokesV = 2*p4;
}
