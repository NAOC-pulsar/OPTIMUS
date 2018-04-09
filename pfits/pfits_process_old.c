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

//gcc -lm -o pfits_describe pfits_describe.c pfits_setup.c -I/Users/hob044/hob044/software/cfitsio/include/ -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio 

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "pfits.h"

void processSubint(dSetStruct *dSet,fitsfile *outfptr,int polOut,int polType,int fScrunch,int outFormat,int debug);

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  int debug=0;
  int i;
  char outFile[1024],oname[1024]; // Output filename
  int setOut=0;
  fitsfile *outfptr;
  int status=0;
  int ii,j,k;
  int polOut=-1;
  int polType=-1;
  int fScrunch=0;
  int outFormat = -1; // 0 = int, 1 = 1 bit, 2 = 2 bit, 4 = 4 bit, 8 = 8 bit, -1 = input format
  
  // Initialise everything
  initialise(&dSet,debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-o")==0)
	{
	  strcpy(outFile,argv[++i]);
	  setOut=1;
	}
      else if (strcmp(argv[i],"-p")==0)
	sscanf(argv[++i],"%d",&polOut);
     else if (strcmp(argv[i],"-polT")==0)
	sscanf(argv[++i],"%d",&polType);
     else if (strcmp(argv[i],"-b")==0)
       sscanf(argv[++i],"%d",&outFormat);
     else if (strcmp(argv[i],"-F")==0)
       fScrunch = -1;
    }

  if (setOut==0)
    errorStop("Must provide an output filename using -o\n",dSet,debug);

  // Open the input file
  pfitsOpenFile(dSet,debug);

  // Load header information
  pfitsLoadHeader(dSet,debug);


  // Open the output file
  sprintf(oname,"!%s",outFile);
  if (!fits_create_file(&outfptr,oname,&status))
    {
      // Copy the file
      ii=1;
      while( !fits_movabs_hdu(dSet->fp, ii++, NULL, &status) )
	  fits_copy_hdu(dSet->fp, outfptr, 0, &status);
    }
  status=0;

  // Now process the primary header
  printf("Not updating the primary header\n");
  // Now process the history table
  printf("Not updating the history table\n");
  // Now process the bandpass table
  printf("Not updating the bandpass table\n");


  // Now process the subint table
  processSubint(dSet,outfptr,polOut,polType,fScrunch,outFormat,debug);
  
  fits_close_file(outfptr,&status);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }

  // Close the file
  //  pfitsCloseFile(dSet,debug);
  
  // De-allocate the memory
  free(dSet);
}

void processSubint(dSetStruct *dSet,fitsfile *outfptr,int polOut,int polType,int fScrunch,int outFormat,int debug)
{
  int polIn;
  int status=0;
  int colnum;
  long naxes[3];
  long s,i,j,k;
  float *arr_p0;
  float *arr_p1;
  float *arr_p2;
  float *arr_p3;
  unsigned char *outArrChar;
  unsigned int  *outArrInt;
  int n_p0,n_p1,n_p2,n_p3;
  int n_f0,n_f1,n_f2,n_f3;
  long sp0,sp1,sp2,sp3;
  int nChanOut = dSet->head->nchan;
  
  // Check that we're doing something
  if (polOut!=-1)
    {
      polIn = dSet->head->npol;
      if (polIn!=polOut)
	{	  
	  printf("inPol = %d, polOut = %d\n",polIn,polOut);
	  if (dSet->head->nbits != 8)
	    errorStop("Sorry, polarisation scrunching currently only works for 8-bit data\n",dSet,debug);
	}
    }
  else
    polOut = 4;

  if (fScrunch == -1)
    nChanOut=1;
  
  // Delete the existing DATA column
  fits_get_colnum(outfptr,CASEINSEN,"DATA",&colnum,&status);
  fits_delete_col(outfptr,colnum,&status);

  if (outFormat==-1)
    {

      outFormat = dSet->head->nbits;
      if (outFormat==32) outFormat=0;
    }
  if (outFormat==0)
    printf("Number of bits in output file: 32\n");
  else
    printf("Number of bits in output file: %d\n",outFormat);
  
  // Create a new DATA column
  if (outFormat == 8)
    fits_insert_col(outfptr,colnum,"DATA","1B",&status);
  else if (outFormat == 0)
    fits_insert_col(outfptr,colnum,"DATA","1I",&status);
  else
    errorStop("Sorry cannot process this output format\n",dSet,debug);
    
  // Re-dimension to the correct size
  fits_modify_vector_len(outfptr,colnum,nChanOut*dSet->head->nsblk*polOut,&status);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }
  
  // write tdim
  naxes[0] = nChanOut;
  naxes[1] = polOut;
  naxes[2] = dSet->head->nsblk;	     
  fits_write_tdim(outfptr,colnum,3,naxes,&status);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    } 
  
  // Set array size
  arr_p0 = (float *)malloc(sizeof(float)*dSet->head->nsblk*dSet->head->nchan);
  arr_p1 = (float *)malloc(sizeof(float)*dSet->head->nsblk*dSet->head->nchan);
  arr_p2 = (float *)malloc(sizeof(float)*dSet->head->nsblk*dSet->head->nchan);
  arr_p3 = (float *)malloc(sizeof(float)*dSet->head->nsblk*dSet->head->nchan);  

  if (outFormat == 8)
    outArrChar = (unsigned char *)malloc(sizeof(unsigned char)*dSet->head->nsblk*nChanOut*polOut);
  else if (outFormat == 0)
    outArrInt = (unsigned int *)malloc(sizeof(unsigned int)*dSet->head->nsblk*nChanOut*polOut);

  // Now process each subint at a time
  for (s=1;s<=dSet->head->nsub;s++)
    {
      // Read in the original data
      pfits_read1pol_float(arr_p0,0,dSet,s,s,1,&sp0,&n_p0,&n_f0,debug);
      pfits_read1pol_float(arr_p1,1,dSet,s,s,1,&sp1,&n_p1,&n_f1,debug);
      pfits_read1pol_float(arr_p2,2,dSet,s,s,1,&sp2,&n_p2,&n_f2,debug);
      pfits_read1pol_float(arr_p3,3,dSet,s,s,1,&sp3,&n_p3,&n_f3,debug);    

      for (i=0;i<dSet->head->nsblk*dSet->head->nchan;i++)
	{
	  if (polOut!=-1)
	    {
	      if (polType==-1)
		arr_p0[i] = ((arr_p0[i]+arr_p1[i])/2.0+0.5);
	      else if (polType==0)
		arr_p0[i] = (arr_p0[i]);
	      else if (polType==1)
		arr_p0[i] = (arr_p1[i]);
	      else if (polType==2)
		arr_p0[i] = (arr_p2[i]);
	      else if (polType==3)
		arr_p0[i] = (arr_p3[i]);
	    }	 	  
	}

      // Process the frequency information
      if (fScrunch==-1)
	{
	  int j;
	  float temp0,temp1,temp2,temp3;
	  for (i=0;i<dSet->head->nsblk;i++)
	    {
	      temp0=temp1=temp2=temp3=0.0;
	      for (j=0;j<dSet->head->nchan;j++)
		{
		  temp0 += arr_p0[i*dSet->head->nchan+j];
		  if (polOut==4)
		    {
		      temp1 += arr_p1[i*dSet->head->nchan+j];
		      temp2 += arr_p2[i*dSet->head->nchan+j];
		      temp3 += arr_p3[i*dSet->head->nchan+j];
		    }
		}
	      temp0/=(float)dSet->head->nchan;
	      arr_p0[i] = temp0;
	      if (polOut==4)
		{
		  arr_p1[i] = temp1/(float)dSet->head->nchan;
		  arr_p2[i] = temp2/(float)dSet->head->nchan;
		  arr_p3[i] = temp3/(float)dSet->head->nchan;
		}
	      //	      printf("output: %g\n",temp0);
	    }
	}
      
      // Now write the output data
      for (i=0;i<dSet->head->nsblk;i++)
	{
	  for (j=0;j<nChanOut;j++)
	    {
	      for (k=0;k<polOut;k++)
		{
		  if (outFormat == 8)
		    {
		      if (k==0)
			outArrChar[i*nChanOut*polOut+k*nChanOut+j] = (unsigned char)(arr_p0[i*nChanOut+j]);
		      else if (k==1)
			outArrChar[i*nChanOut*polOut+k*nChanOut+j] = (unsigned char)(arr_p1[i*nChanOut+j]);
		      else if (k==2)
			outArrChar[i*nChanOut*polOut+k*nChanOut+j] = (unsigned char)(arr_p2[i*nChanOut+j]);
		      else if (k==3)
			outArrChar[i*nChanOut*polOut+k*nChanOut+j] = (unsigned char)(arr_p3[i*nChanOut+j]);
		    }
		  else if (outFormat == 0)
		    {
		      if (k==0)
			outArrInt[i*nChanOut*polOut+k*nChanOut+j] = (unsigned int)((arr_p0[i*nChanOut+j]-120)*21);
		      else if (k==1)
			outArrInt[i*nChanOut*polOut+k*nChanOut+j] = (unsigned int)((arr_p1[i*nChanOut+j]-120)*21);
		      else if (k==2)
			outArrInt[i*nChanOut*polOut+k*nChanOut+j] = (unsigned int)((arr_p2[i*nChanOut+j]-120)*21);
		      else if (k==3)
			outArrInt[i*nChanOut*polOut+k*nChanOut+j] = (unsigned int)((arr_p3[i*nChanOut+j]-120)*21);
		    }
		}
	    }
	}
      if (outFormat==8)
	fits_write_col_byt(outfptr,colnum,s,1,dSet->head->nsblk*polOut*nChanOut,outArrChar,&status);
      else if (outFormat == 0)
	fits_write_col(outfptr,TINT,colnum,s,1,dSet->head->nsblk*polOut*nChanOut,outArrInt,&status);
    }

  
  // LOTS OF PROBLEMS:
  // 1. DAT_SCL etc.
  // 2. Should not change the primary NPOL header - should change the history
  // 3. Currently only works for 8bit data
  free(arr_p0);
  free(arr_p1);
  free(arr_p2);
  free(arr_p3);
  if (outFormat == 8)
    free(outArrChar);
  else if (outFormat == 0)  
    free(outArrInt);
  // Now update header information
  fits_movnam_hdu(outfptr,BINARY_TBL,"SUBINT",0,&status);
  fits_update_key(outfptr,TINT,"NPOL",&polOut,NULL,&status);
  fits_update_key(outfptr,TINT,"NCHAN",&nChanOut,NULL,&status);
  if (outFormat==0)
    {
      int nbits=32;
      fits_update_key(outfptr,TINT,"NBITS",&nbits,NULL,&status);
    }
}
