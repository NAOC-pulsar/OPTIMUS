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

void process(dSetStruct *dSet,fitsfile *outfptr,float dm,float fref,int debug);
void processDataSegment(dSetStruct *dSet,unsigned char *inArrChar_1,int n,int nsamp,float *fltArray);
void multMatrix(unsigned char p1,unsigned char p2,unsigned char p3,unsigned char p4,float mat[4][4],
		float *out1,float *out2,float *out3,float *out4);

int main(int argc,char *argv[])
{
  dSetStruct *dSet;
  int debug=0;
  int status=0;
  int i,ii;
  char outFile[1024],oname[1024]; // Output filename
  int setOut=0;
  fitsfile *outfptr;
  float dm=0;
  float fref=-1;
  
  // Initialise everything
  initialise(&dSet,debug);
  
  // Read inputs
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	setFilename(argv[++i],dSet,debug);
      else if (strcmp(argv[i],"-dm")==0)
	sscanf(argv[++i],"%f",&dm);
      else if (strcmp(argv[i],"-fref")==0)
	sscanf(argv[++i],"%f",&fref);
      else if (strcmp(argv[i],"-o")==0)
	{
	  strcpy(outFile,argv[++i]);
	  setOut=1;
	}
    }

  if (setOut==0)
    errorStop("Must provide an output filename using -o\n",dSet,debug);

  pfitsOpenFile(dSet,debug);
  pfitsLoadHeader(dSet,debug);

  // Open the output file
  printf("Making a copy of the input file\n");
  sprintf(oname,"!%s",outFile);
  if (!fits_create_file(&outfptr,oname,&status))
    {
      // Copy the file
      ii=1;
      while( !fits_movabs_hdu(dSet->fp, ii++, NULL, &status) )
	fits_copy_hdu(dSet->fp, outfptr, 0, &status);
    }
  status=0;
  printf("Complete making a copy\n");

  
  // Close the file
  //  pfitsCloseFile(dSet,debug);

  process(dSet,outfptr,dm,fref,debug);

  fits_close_file(outfptr,&status);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }
  
  // De-allocate the memory
  deallocateMemory(&dSet,debug);
}

void process(dSetStruct *dSet,fitsfile *outfptr,float dm,float fref,int debug)
{
  // Dedisperse by calculating each point in turn in the output file, based on the input file
  // Later try ring buffer method and/or FFT method
  int status=0;
  long s,i,j,k;
  long l_samp,l_chan,l_pol;

  unsigned int  *outArrInt;
  unsigned char cSignal;
  int initflag=0;
  
  double toff;
  long   sampleOff;

  unsigned char *outArrChar,*inArrChar_1,*inArrChar_2;
  float *outArrFlt;
  float *fltArray;
  int sub;
  long s1,s2,c,c_in;
  long i_in,j_in,p_in;
  float fChan[dSet->head->nchan];
  float chanbw = dSet->head->chanbw;
  int offChan[dSet->head->nchan];
  int nChanOut;
  int nSubReq;
  int colnum_in,colnum_out;
  long l_sub,r_sub;
  unsigned char nval = 0;
  int f,p;
  unsigned char *posPtr;
  int ptrPos=1,usePtrPos=1;
  long nLoad;
  int curPtr;
  long naxes[3];
  int nbitsOut=32;
  int freqAdd=1;
  int c_out;

  if (freqAdd==0)
    nChanOut = dSet->head->nchan;
  else if (freqAdd==1)
    nChanOut = 1;

  for (i=0;i<dSet->head->nchan;i++)
    fChan[i] = dSet->head->chanFreq[i];
  
  if (fref < 0)
    fref = dSet->head->freq;

  printf("Reference frequency: %g\n",fref);


  fits_update_key(outfptr,TINT,"NBITS",&nbitsOut,NULL,&status);

  
  fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum_in,&status);
  fits_get_colnum(outfptr,CASEINSEN,"DATA",&colnum_out,&status);
  fits_delete_col(outfptr,colnum_out,&status);
  //  fits_insert_col(outfptr,colnum_out,"DATA","1I",&status);
  fits_insert_col(outfptr,colnum_out,"DATA","1E",&status);
  //  fits_modify_vector_len(outfptr_out,colnum,nChanOut*dSet->head->nsblk*polOut,&status);
  fits_modify_vector_len(outfptr,colnum_out,nChanOut*dSet->head->nsblk*dSet->head->npol,&status);

  // write tdim
  naxes[0] = nChanOut;
  naxes[1] = dSet->head->npol; //polOut;
  naxes[2] = dSet->head->nsblk;	     
  fits_write_tdim(outfptr,colnum_out,3,naxes,&status);

  for (i=0;i<dSet->head->nchan;i++)
    {
      fChan[i] = dSet->head->chanFreq[i];
      toff = 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(fChan[i]/1000.0,-2));
      offChan[i] = toff/dSet->head->tsamp;
      printf("offChan[i] = %d\n",offChan[i]);
    }
  // How many samples for the DM smearing across the band?
  toff = 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(fChan[0]/1000.0,-2));
  s1 = toff/dSet->head->tsamp;

  toff = 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(fChan[dSet->head->nchan-1]/1000.0,-2));
  s2 = toff/dSet->head->tsamp;
  printf("s1 and s2 are %ld %ld\n",s1,s2);
  // Calculate how many subints are required
  nSubReq = (int)(fabs((double)(s1-s2)/dSet->head->nsblk)+1);
  printf("nSubReq = %d\n",nSubReq);

  fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum_in,&status);
  fits_get_colnum(outfptr,CASEINSEN,"DATA",&colnum_out,&status);
  
  // Read in the required subints of data

  inArrChar_1 = (unsigned char *)malloc(sizeof(unsigned char)*nSubReq*dSet->head->nsblk*dSet->head->npol*dSet->head->nchan*2);
  inArrChar_2 = (unsigned char *)malloc(sizeof(unsigned char)*nSubReq*dSet->head->nsblk*dSet->head->npol*dSet->head->nchan);
  fltArray = (float *)malloc(sizeof(float)*nSubReq*dSet->head->nsblk*dSet->head->npol*dSet->head->nchan*2);
  posPtr = inArrChar_1;
  
  // Output a single subint at a time
  outArrChar = (unsigned char *)malloc(sizeof(unsigned char)*dSet->head->nsblk*dSet->head->npol*nChanOut);
  outArrFlt = (float *)malloc(sizeof(float)*dSet->head->nsblk*dSet->head->npol*nChanOut);

  l_sub = 0;
  r_sub = 0;
  // Read in the first reqSub*2 amount of data
  for (l_sub=0;l_sub<nSubReq*2;l_sub++)
    {
      fits_read_col_byt(dSet->fp,colnum_in,l_sub+1,
			1,dSet->head->nsblk*dSet->head->npol*dSet->head->nchan
			,nval,inArrChar_1+l_sub*dSet->head->nsblk*dSet->head->npol*dSet->head->nchan,&initflag,&status); 
    }
  processDataSegment(dSet,inArrChar_1,dSet->head->nsblk*dSet->head->npol*dSet->head->nchan*nSubReq*2,nSubReq*2*dSet->head->nsblk,fltArray);

  nLoad = dSet->head->nsblk*dSet->head->npol*dSet->head->nchan*nSubReq*2;
  r_sub = nSubReq*2;
  
  curPtr=0;
  for (l_sub=0;l_sub<dSet->head->nsub;l_sub++)
    {
      printf("Processing subint %ld out of %d\n",l_sub+1,dSet->head->nsub);
      c=0;
      c_out=0;

      if (l_sub!=0)
	{
	  // Read in another subint
	  if (r_sub == dSet->head->nsub-1)
	    printf("Running off the end of the file\n");
	  else
	    {
	      fits_read_col_byt(dSet->fp,colnum_in,r_sub+1,
				1,dSet->head->nsblk*dSet->head->npol*dSet->head->nchan
				,nval,inArrChar_1+curPtr*dSet->head->nsblk*dSet->head->npol*dSet->head->nchan,&initflag,&status);		      
	      processDataSegment(dSet,inArrChar_1+curPtr*dSet->head->nsblk*dSet->head->npol*dSet->head->nchan,dSet->head->nsblk*dSet->head->npol*dSet->head->nchan,dSet->head->nsblk,fltArray+curPtr*dSet->head->nsblk*dSet->head->npol*dSet->head->nchan);
	      r_sub++;
	      curPtr++;
	      if (curPtr==nSubReq*2)
		curPtr=0;
	    }
	}
      
      // Write it out
      for (i=0;i<dSet->head->nsblk;i++)
	{
	  //	  printf("i = %d %d %ld\n",i,dSet->head->nsblk,c);
	  for (p=0;p<dSet->head->npol;p++)
	    {
	      for (j=0;j<dSet->head->nchan;j++)
		{
		  c_in = (c+curPtr*dSet->head->npol*dSet->head->nchan*dSet->head->nsblk)-offChan[j]*dSet->head->npol*dSet->head->nchan;
		  if (c_in >= nLoad) c_in-=nLoad;
		  /*		  if (c_in >= nLoad)
		    {
		      nLoad+=dSet->head->nsblk*dSet->head->npol*dSet->head->nchan;
		      r_sub++;
		      printf("Boo\n");
		      if (ptrPos==1)
			{
			  printf("In here\n");
			  fits_read_col_byt(dSet->fp,colnum_in,r_sub+1,
					    1,dSet->head->nsblk*dSet->head->npol*dSet->head->nchan
					    ,nval,inArrChar_2,&initflag,&status);		      
			}
		      else
			{
			  fits_read_col_byt(dSet->fp,colnum_in,r_sub+1,
					    1,dSet->head->nsblk*dSet->head->npol*dSet->head->nchan
					    ,nval,inArrChar_1,&initflag,&status);		      
			}
			}  */
		  //		  printf("i = %d, c = %ld %ld\n",i,c,c_in);
		  if (c_in < 0) 
		    outArrChar[c_out] = 0;
		  else
		    {
		      //		      if (c_in >= dSet->head->nsblk*dSet->head->npol*dSet->head->nchan)
		      //			{
		      //			  usePtrPos = ptrPos*-1;
		      //			  c_in -= dSet->head->nsblk*dSet->head->npol*dSet->head->nchan;
			  //			  printf("Have: %d %d\n",c_in,dSet->head->nsblk*dSet->head->npol*dSet->head->nchan);
		      //			}
		      //		      else
		      //			usePtrPos = ptrPos;
		      //		      if (c_in > dSet->head->nsblk*dSet->head->npol*dSet->head->nchan)
		      //			{
		      //			  printf("This shouldn't happen\n");
		      //			  exit(1);
		      //			}
		      
		      //		      if (usePtrPos==1)
		      //		outArrChar[c] = fltArray[c_in];
		      if (freqAdd==0)
			outArrFlt[c_out] = fltArray[c_in];
		      else
			{
			  if (j==0) outArrFlt[c_out] = fltArray[c_in]; // THIS IS WRONG BECAUSE OF THE C IN THE C_IN COMMAND
			  else outArrFlt[c_out]+=fltArray[c_in];
			}
			//		      else
			//			outArrChar[c] = inArrChar_2[c_in];
		    }
		  if (freqAdd==0) c_out++;
		  c++;
		}
	      if (freqAdd==1) c_out++;
	    }
	}

      
      //      fits_write_col_byt(outfptr,colnum_out,l_sub+1,1,
      //			 dSet->head->nsblk*dSet->head->npol*dSet->head->nchan,outArrChar,&status);
      printf("Writing floats %d\n",l_sub);
      fits_write_col(outfptr,TFLOAT,colnum_out,l_sub+1,1,dSet->head->nsblk*dSet->head->npol*nChanOut,outArrFlt,&status);
      fits_report_error(stderr,status);
    }
  fits_report_error(stderr,status);
  status=0;
  
  fits_movnam_hdu(outfptr,BINARY_TBL,"SUBINT",0,&status);
  fits_update_key(outfptr,TINT,"NCHAN",&nChanOut,NULL,&status);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }
  free(inArrChar_1);
  free(inArrChar_2);
  free(outArrChar);
  free(fltArray);
  
}

void processDataSegment(dSetStruct *dSet,unsigned char *inArrChar_1,int n,int nsamp,float *fltArray)
{
  int i,j,k;
  long c=0;
  unsigned char *p1;
  unsigned char *p2;
  unsigned char *p3;
  unsigned char *p4;
  float mat[4][4];
  float *out1,*out2,*out3,*out4;

  mat[0][0] = 1;
  mat[1][0] = 1;
  mat[2][0] = 0;
  mat[3][0] = 0;

  mat[0][1] = 1;
  mat[1][1] = -1;
  mat[2][1] = 0;
  mat[3][1] = 0;

  mat[0][2] = 0;
  mat[1][2] = 0;
  mat[2][2] = 2;
  mat[3][2] = 0;

  mat[0][3] = 0;
  mat[1][3] = 0;
  mat[2][3] = 0;
  mat[3][3] = 2; 


  /*  mat[0][0] = 1;
  mat[1][0] = 0;
  mat[2][0] = 0;
  mat[3][0] = 0;

  mat[0][1] = 0;
  mat[1][1] = 1;
  mat[2][1] = 0;
  mat[3][1] = 0;

  mat[0][2] = 0;
  mat[1][2] = 0;
  mat[2][2] = 1;
  mat[3][2] = 0;

  mat[0][3] = 0;
  mat[1][3] = 0;
  mat[2][3] = 0;
  mat[3][3] = 1;*/


  p1 = inArrChar_1;
  p2 = inArrChar_1 + dSet->head->nchan;
  p3 = inArrChar_1 + 2*dSet->head->nchan;
  p4 = inArrChar_1 + 3*dSet->head->nchan;
  
  out1 = fltArray;
  out2 = fltArray + dSet->head->nchan;
  out3 = fltArray + 2*dSet->head->nchan;
  out4 = fltArray + 3*dSet->head->nchan;
  
  for (i=0;i<nsamp;i++)
    {
      for (j=0;j<dSet->head->nchan;j++)
	{           
	  multMatrix(*p1,*p2,*p3,*p4,mat,out1,out2,out3,out4);
	  p1++;
	  p2++;
	  p3++;
	  p4++;
	  out1++;
	  out2++;
	  out3++;
	  out4++;
	}
      p1+=3*dSet->head->nchan;
      p2+=3*dSet->head->nchan;
      p3+=3*dSet->head->nchan;
      p4+=3*dSet->head->nchan;
      out1+=3*dSet->head->nchan;
      out2+=3*dSet->head->nchan;
      out3+=3*dSet->head->nchan;
      out4+=3*dSet->head->nchan;
    }
  
  /*  for (i=0;i<nsamp;i++)
    {

      for (j=0;j<dSet->head->npol;j++)
	{
	  for (k=0;k<dSet->head->nchan;k++)
	    {
	      fltArray[c] = (float)inArrChar_1[c];
	      if (k>100 && k<200) fltArray[c]=0;
	      c++;
	    }
	}
	} */
}

void multMatrix(unsigned char p1,unsigned char p2,unsigned char p3,unsigned char p4,float mat[4][4],
		float *out1,float *out2,float *out3,float *out4)
{
  *out1 = p1*mat[0][0]+p2*mat[1][0]+p3*mat[2][0]+p4*mat[3][0];
  *out2 = p1*mat[0][1]+p2*mat[1][1]+p3*mat[2][1]+p4*mat[3][1];
  *out3 = p1*mat[0][2]+p2*mat[1][2]+p3*mat[2][2]+p4*mat[3][2];
  *out4 = p1*mat[0][3]+p2*mat[1][3]+p3*mat[2][3]+p4*mat[3][3];

}
