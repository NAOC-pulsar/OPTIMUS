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


void dedisperse(dSetStruct *dSet,fitsfile *outfptr,float dm,float fref,int debug);
void dedispRingBuffer(dSetStruct *dSet,fitsfile *outfptr,float dm,float fref,int debug);

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
  sprintf(oname,"!%s",outFile);
  if (!fits_create_file(&outfptr,oname,&status))
    {
      // Copy the file
      ii=1;
      while( !fits_movabs_hdu(dSet->fp, ii++, NULL, &status) )
	fits_copy_hdu(dSet->fp, outfptr, 0, &status);
    }
  status=0;

  
  // Close the file
  //  pfitsCloseFile(dSet,debug);

  dedisperse(dSet,outfptr,dm,fref,debug);

  fits_close_file(outfptr,&status);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }
  
  // De-allocate the memory
  deallocateMemory(&dSet,debug);
}

void dedisperse(dSetStruct *dSet,fitsfile *outfptr,float dm,float fref,int debug)
{
  // Dedisperse by calculating each point in turn in the output file, based on the input file
  // Later try ring buffer method and/or FFT method
  int status=0;
  int colnum_in,colnum_out;
  long s,i,j,k,p,sub;
  long l_sub,l_samp,l_chan,l_pol;

  unsigned char *outArrChar;
  unsigned int  *outArrInt;
  unsigned char cSignal;
  unsigned char nval = 0;
  int initflag=0;
  
  double toff;
  long   sampleOff;

  // Get frequency channels
  float fChan[dSet->head->nchan];
  float chanbw = dSet->head->chanbw;

  int dedispMethod = 3;
  
  // Should read from the file !!!!
  for (i=0;i<dSet->head->nchan;i++)
    fChan[i] = dSet->head->chanFreq[i];
  
  if (fref < 0)
    fref = dSet->head->freq;

  printf("Reference frequency: %g\n",fref);
  
  fits_get_colnum(dSet->fp,CASEINSEN,"DATA",&colnum_in,&status);
  fits_get_colnum(outfptr,CASEINSEN,"DATA",&colnum_out,&status);

  if (dedispMethod == 1)
    { 
      outArrChar = (unsigned char *)malloc(sizeof(unsigned char)*dSet->head->nchan*dSet->head->npol*dSet->head->nsblk);
      printf("Got here\n");
      // Process each output subint in turn
      for (sub=0;sub<dSet->head->nsub;sub++)
	{
	  printf("Processing subint %ld/%d\n",sub+1,dSet->head->nsub);
	  for (s=0;s<dSet->head->nsblk;s++)
	    {
	      for (i=0;i<dSet->head->nchan;i++)
		{
		  // Must find the correct subint, sample, pol and sample to get this data
		  //		  printf("i = %d %d %d\n",(int)i,(int)p,(int)s);
		  //		  double toff;
		  //		  long   sampleOff;
		  //		  long l_sub,l_samp,l_chan,l_pol;
		  toff = 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(fChan[i]/1000.0,-2));
		  sampleOff = toff/dSet->head->tsamp;
		  //		  printf("toff = %g dm = %g sampleOff = %d\n",toff,dm,(int)sampleOff);
		  // Read the relevant byte 		  
		  l_sub = sub;
		  l_samp = s-sampleOff;
		  while (l_samp < 0)
		    {
		      l_sub--;
		      l_samp = dSet->head->nsblk+l_samp;
		    }
		  while (l_samp >= dSet->head->nsblk)
		    {
		      l_sub++;
		      l_samp = l_samp-dSet->head->nsblk;
		    }
		  for (p=0;p<dSet->head->npol;p++)
		    {
		      
		      if (l_sub >= 0 && l_sub < dSet->head->nsub)
			{
			  l_pol = p;		  
			  l_chan = i;
			  
			  fits_read_col_byt(dSet->fp,colnum_in,l_sub+1,
					    1+l_samp*dSet->head->npol*dSet->head->nchan+
					    l_pol*dSet->head->nchan+l_chan
					    ,1,nval,&cSignal,&initflag,&status);		      
			}
		      else
			cSignal=0;
		      
		      outArrChar[s*dSet->head->npol*dSet->head->nchan +
				 p*dSet->head->nchan+i] = cSignal;
		      //		  printf("b\n");
		    }
		}
	    }
	  fits_write_col_byt(outfptr,colnum_out,sub+1,1,
			     dSet->head->nsblk*dSet->head->npol*dSet->head->nchan,
			     outArrChar,&status);
	}
      free(outArrChar);
    }
  else if (dedispMethod == 2)
    {
      double *vals,*ovals;
      double *outFFT;
      long c=0,cc=0;
      fftw_complex *output;
      fftw_plan transform_plan,transform_plan2;

      double real_rotate,ima_rotate,sina,cosina,amp,rot;
      double real,ima;
      double fSpec;
      unsigned char outChar;
      double K = 4149.37759;

      
      vals  = (double *)malloc(sizeof(double)*dSet->head->npol*dSet->head->nsblk*dSet->head->nsub);
      ovals = (double *)malloc(sizeof(double)*dSet->head->npol*dSet->head->nsblk*dSet->head->nsub);
      outFFT = (double *)malloc(sizeof(double)*dSet->head->npol*dSet->head->nsblk*dSet->head->nsub);
      //      output = (fftw_complex*)outFFT;
      output = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*dSet->head->npol*dSet->head->nsblk*dSet->head->nsub);
      for (i=0;i<dSet->head->nchan;i++)
	{
	  printf("Processing channel: %d out of %d\n",(int)i+1,dSet->head->nchan);
	  c=0;
	  for (sub=0;sub<dSet->head->nsub;sub++)
	    {
	      for (s=0;s<dSet->head->nsblk;s++)
		{
		  //		  for (p=0;p<dSet->head->npol;p++)
		  p=0;
		  {
		      fits_read_col_byt(dSet->fp,colnum_in,sub+1,
					1+s*dSet->head->npol*dSet->head->nchan+
					p*dSet->head->nchan+i
					,1,nval,&cSignal,&initflag,&status);		      

		      vals[c++]=cSignal;
		    }
		}
	    }
	  transform_plan = fftw_plan_dft_r2c_1d(c,vals,output,FFTW_ESTIMATE);
	  fftw_execute(transform_plan);
	  // Now inverse FFT


	  for (j=1;j<c/2;j++)
	    {
	      fSpec = j/(c*dSet->head->tsamp);
	      rot = 2.0*M_PI*(K*dm*fSpec)*(1.0/(fChan[i]*fChan[i])-1.0/(fref*fref));

	      real = output[j][0];
	      ima = output[j][1];

	      if (real != 0 || ima != 0)
		{
		  
		  amp=sqrt(real*real+ima*ima);
		  cosina=real/amp;
		  sina=ima/amp;
		  
		  // rotate profile
		  real_rotate=amp*(cosina*cos(rot)-sina*sin(rot));
		  ima_rotate =amp*(sina*cos(rot)+cosina*sin(rot));
		}
	      else
		{
		  real_rotate = 0;
		  ima_rotate = 0;
		}

	      output[j][0] = real_rotate;
	      output[j][1] = ima_rotate;
	    }

	  transform_plan2 = fftw_plan_dft_c2r_1d(c,output,ovals,FFTW_ESTIMATE);
	  fftw_execute(transform_plan2);

	  fftw_destroy_plan(transform_plan);

	  fftw_destroy_plan(transform_plan2);

	  //	  for (i=0;i<c;i++)
	  //	    printf("vals %g %g\n",vals[i],ovals[i]/c);
	    //	    printf("vals %g\n",outFFT[2*i]*outFFT[2*i] + outFFT[2*i+1]*outFFT[2*i+1]);
	  //	  	  printf("Stopping\n");
	  //	  	  	  exit(1);

	  // Writing output
	  printf("Writing output\n");
	  cc=0;
	  for (sub=0;sub<dSet->head->nsub;sub++)
	    {
	      for (s=0;s<dSet->head->nsblk;s++)
		{
		  if (ovals[cc]/c > 255)
		    outChar = 255;
		  else
		    outChar = (unsigned char)(ovals[cc]/c);
		  cc++;
		  
		  for (p=0;p<dSet->head->npol;p++)
		    {

		      //		      printf("Writing to %d %d %d\n",sub+1,1+s*dSet->head->npol*dSet->head->nchan+
		      //			     p*dSet->head->nchan+i,status);
		      fits_write_col_byt(outfptr,colnum_out,sub+1,1+s*dSet->head->npol*dSet->head->nchan+
					 p*dSet->head->nchan+i, 1,
					 &outChar,&status);
		    }
		}
	    }
	  if (status)
	    {
	      fits_report_error(stderr,status);
	      exit(1);
	    }
	}
      free(vals);
    }
  else if (dedispMethod==3)
    dedispRingBuffer(dSet,outfptr,dm,fref,debug);
  
}

void dedispRingBuffer(dSetStruct *dSet,fitsfile *outfptr,float dm,float fref,int debug)
{
  unsigned char *outArrChar,*inArrChar_1,*inArrChar_2;
  int sub,s,i,j;
  long s1,s2,c,c_in;
  long i_in,j_in,p_in;
  double toff;
  float fChan[dSet->head->nchan];
  float chanbw = dSet->head->chanbw;
  int offChan[dSet->head->nchan];
  int nSubReq;
  int colnum_in,colnum_out;
  long l_sub,r_sub;
  int status=0;
  int initflag=0;
  unsigned char nval = 0;
  int f,p;
  unsigned char *posPtr;
  int ptrPos=1,usePtrPos=1;
  long nLoad;
  int curPtr;
  
  for (i=0;i<dSet->head->nchan;i++)
    {
      fChan[i] = dSet->head->chanFreq[i];
      toff = 4.15e-3*dm*(pow(fref/1000.0,-2)-pow(fChan[i]/1000.0,-2));
      offChan[i] = toff/dSet->head->tsamp;
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
  posPtr = inArrChar_1;
  
  // Output a single subint at a time
  outArrChar = (unsigned char *)malloc(sizeof(unsigned char)*dSet->head->nsblk*dSet->head->npol*dSet->head->nchan);

  l_sub = 0;
  r_sub = 0;
  // Read in the first reqSub*2 amount of data
  for (l_sub=0;l_sub<nSubReq*2;l_sub++)
    {
      fits_read_col_byt(dSet->fp,colnum_in,l_sub+1,
			1,dSet->head->nsblk*dSet->head->npol*dSet->head->nchan
			,nval,inArrChar_1+l_sub*dSet->head->nsblk*dSet->head->npol*dSet->head->nchan,&initflag,&status);		      
    }

  nLoad = dSet->head->nsblk*dSet->head->npol*dSet->head->nchan*nSubReq*2;
  r_sub = nSubReq*2;
  
  curPtr=0;
  for (l_sub=0;l_sub<dSet->head->nsub;l_sub++)
    {
      printf("Processing subint %ld out of %d\n",l_sub+1,dSet->head->nsub);
      c=0;

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
		    outArrChar[c] = 0;
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
			outArrChar[c] = inArrChar_1[c_in];
			//		      else
			//			outArrChar[c] = inArrChar_2[c_in];
		    }
		  c++;
		}
	    }
	}

      
      fits_write_col_byt(outfptr,colnum_out,l_sub+1,1,
			 dSet->head->nsblk*dSet->head->npol*dSet->head->nchan,outArrChar,&status);

    }

  free(inArrChar_1);
  free(inArrChar_2);
  free(outArrChar);
}
