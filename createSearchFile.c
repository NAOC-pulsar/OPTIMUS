// Routine to take the output from the simulated data files to produce
// a search-mode file that can be processed using pfits, sigproc or presto
//
// Compile with: gcc -lm -o createSearchFile createSearchFile.c T2toolkit.c -L/Users/hob044/hob044/software/cfitsio/lib -lcfitsio -O3

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "T2toolkit.h"
#include "simulate.h"
#include "fitsio.h"

#define MAX_FILES 10

double psrangle(double centre_long,double centre_lat,double psr_long,double psr_lat);
double beamScaling(double angle,float diameter,float freq);
void digitiseValues(float *chanVals,int nchan,int nbits,unsigned char *chanDigitised,float level);
void createPrimaryHeader(fitsfile *fp,float f1,float f2,header *head);
void deleteUnwantedTables(fitsfile *fp);
void createSubintHeader(fitsfile *fp,float f1,float f2,int nchan,float tsamp,int nbits);
int twoBitDigitise(unsigned char *inArray, unsigned char *outArray, unsigned int n);
double mjd2year(double mjd,int *retDay,int *retMonth,int *retYr);


int main(int argc, char *argv[])
{
  header **head;
  header *outputHeader;
  double t;
  float val;
  
  char outName[1024];
  int i;
  FILE *fout,*fout2;
  int debugOut=0;
  double debugOutT1;
  double debugOutT2;
  int debugOutF1;
  int debugOutF2;
  char debugOutName[1024];

  FILE *fin[MAX_FILES];
  float angle[MAX_FILES];
  float *scaling;
  int nfiles=0;
  int j;
  float readVal;
  float *chanVals;
  int nbits;
  int nsampleByte;
  int nsblk=4096;
  int npol=1;
  unsigned char *chanDigitised;
  fitsfile *fp;
  char fname[1024];
  int status=0;
  int iblk=0;
  int sub=0;
  float *datFreq;
  float *datWts;
  float *datOffs;
  float *datScl;
  float levelSetValue=1;
  double beamRA; 
  double beamDEC;
  float diameter; // Diameter of telescope (m)

  unsigned long angCount=0;
  unsigned long angUpdate;

  char paramFile[MAX_PARAM_FILES][1024];
  char fileNames[MAX_PARAM_FILES][1024];
  int  useParamFile=0;

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
      else if (strcmp(argv[i],"-f")==0) // Input file
	{
	  if (!(fin[nfiles] = fopen(argv[++i],"rb")))
	    {
	      printf("Unable to open file: %s\n",argv[i]);
	      exit(1);
	    }
	  strcpy(fileNames[nfiles],argv[i]);
	  nfiles++;
	}

    }
  outputHeader = (header *)malloc(sizeof(header));
  simulateSetHeaderDefaults(outputHeader);
  if (useParamFile>0)
    {
      for (i=0;i<useParamFile;i++)
	simulateReadParamFile(outputHeader,paramFile[i]);
    }
  nbits = outputHeader->nbits;
  nsampleByte = 8/nbits;
  printf("Nbits = %d level setting = %d\n",nbits,outputHeader->levelSet);
  if (outputHeader->levelSet==1)
    {
      FILE *fin_level;
      char str[1024];
      for (i=0;i<nfiles;i++)
	{
	  sprintf(str,"%s.levelSetting",fileNames[i]);
	 
	  if (!(fin_level = fopen(str,"r")))
	    {
	      // No level setting available
	    }
	  else
	    {
	      printf("Reading level setting from %s\n",str);
	      fscanf(fin_level,"%f",&levelSetValue);
	      fclose(fin_level);
	      break;
	    }
	}
    }
  diameter = outputHeader->diameter;
  
  // Check beam positions
  if (outputHeader->setBeam==1)
    {
      beamRA = outputHeader->beamRA0;
      beamDEC = outputHeader->beamDEC0;
    }
  else
    {
      beamRA = beamDEC = 0;
    }
  
  angUpdate = (unsigned long)(0.1/outputHeader->tsamp);
  printf("Update angle every %ld time samples\n",angUpdate);

  scaling = (float *)malloc(sizeof(float)*MAX_FILES*outputHeader->nchan);
  chanVals = (float *)malloc(sizeof(float)*outputHeader->nchan);
  chanDigitised = (unsigned char *)malloc(sizeof(unsigned char)*outputHeader->nchan/nsampleByte*nsblk);
  datFreq = (float *)malloc(sizeof(float)*outputHeader->nchan);
  datWts = (float *)malloc(sizeof(float)*outputHeader->nchan);
  datOffs = (float *)malloc(sizeof(float)*outputHeader->nchan);
  datScl = (float *)malloc(sizeof(float)*outputHeader->nchan);

  for (i=0;i<outputHeader->nchan;i++)
    {
      datFreq[i] = outputHeader->f1 + (float)i*(outputHeader->f2-outputHeader->f1)/outputHeader->nchan;
      datOffs[i] = 0;
      datWts[i] = 1;
      datScl[i] = 1;
    }
  
  head = (header **)malloc(sizeof(header *)*nfiles);
  for (i=0;i<nfiles;i++)
    {
      head[i] = (header *)malloc(sizeof(header));
      simulateReadHeaderParameters(head[i],fin[i]);
      // Check that the parameters are consistent
      if (head[i]->nchan != outputHeader->nchan)
	{
	  printf("ERROR: File '%s' has different number of channels to the requested output\n",head[i]->name);
	  exit(1);
	}
      if (head[i]->tsamp != outputHeader->tsamp)
	{
	  printf("ERROR: File '%s' has different sampling time to the requested output\n",head[i]->name);
	  exit(1);
	}
      if (head[i]->t1 - head[i]->t0 < outputHeader->t1-outputHeader->t0)
	{
	  printf("ERROR: File '%s' does not have enough data points\n",head[i]->name);
	  exit(1);
	}
      if (head[i]->f1 != outputHeader->f1)
	  printf("WARNING: File '%s' has different f1 to the requested output\n",head[i]->name);
      if (head[i]->f2 != outputHeader->f2)
	  printf("WARNING: File '%s' has different f2 to the requested output\n",head[i]->name);
    }  

  if (strcmp(outName,"default")==0)
    {
      double mjd = outputHeader->imjd + (outputHeader->smjd+outputHeader->stt_offs)/86400.0;
      int day,month,year;
      int hour,min,sec;
      double timeVal;
      int bm=12;
      
      printf("MJD is %g\n",(double)mjd);
      mjd2year(mjd,&day,&month,&year);
      year -= 2000.0;

      timeVal = (mjd - (int)mjd)*86400.0;

      hour = (int)(timeVal/60.0/60.0);
      min  = (int)((timeVal - hour*60*60)/60);
      sec  = (int)(timeVal - hour*60.0*60.0 - min*60.0);
      printf("Hour is %d\n",hour);
      printf("Minute is %d\n",min);
      printf("Filename is sim_%02d%02d%02d_%02d%02d%02d_beam%02d.sf\n",year,month,day,hour,min,sec,bm);
      exit(1);
    }
  
  sprintf(fname,"!%s(psrheader.fits)",outName);
  
  fits_create_file(&fp,fname,&status);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }
  createPrimaryHeader(fp,outputHeader->f1,outputHeader->f2,outputHeader);
  deleteUnwantedTables(fp);
  createSubintHeader(fp,outputHeader->f1,outputHeader->f2,outputHeader->nchan,outputHeader->tsamp,outputHeader->nbits);
  
  if (debugOut==1)
    fout2 = fopen(debugOutName,"w");
  for (t=outputHeader->t0;t<outputHeader->t1;t+=outputHeader->tsamp)
    {
      // Update angles and scaling factors
      // Ensure that we set the angles at the beginning.
      if (angCount==0)
	{
	  for (i=0;i<outputHeader->nchan;i++)
	    {
	      for (j=0;j<nfiles;j++)
		{
		  if (head[j]->useAngle==1 && outputHeader->setBeam == 1)
		    {
		      angle[j] = psrangle(beamRA,beamDEC,head[j]->raj_rad,head[j]->decj_rad); // In degrees
		      scaling[i*nfiles+j] = beamScaling(angle[j],diameter,datFreq[i]); 
		      if (j==1 && i==10)
			printf("Angle is %g %g %g %g\n",t, angle[j],datFreq[i],scaling[i*nfiles+j]);
		    }
		  else {
		    scaling[i*nfiles+j] = 1;
		    // printf("No scalling angle is %g %g i=%d j=%d\n",t, scaling[i*nfiles+j], i, j);
		  }
		}
	    }

	}
      angCount++;
      if (outputHeader->surveyType==2)
	beamRA+=(outputHeader->tsamp/86164.090530833)*2.0*M_PI;
	// GH: changed on 26th April 2017 to use the number of seconds in a sidereal day and not to use cos(dec) here
	//	beamRA+=(outputHeader->tsamp/86400.0*2.0*M_PI*cos(beamDEC));
      if (angCount==angUpdate)
	angCount=0;
      if (debugOut == 1 && t >= debugOutT1 && t < debugOutT2)
	fprintf(fout2,"%f ",t);
      
      for (i=0;i<outputHeader->nchan;i++)
	{
	  chanVals[i]=0;
	  for (j=0;j<nfiles;j++)
	    {
	      fread(&readVal,sizeof(float),1,fin[j]);
	      //	      if (j==0 && i==0)
	      //		printf("scaling: %g\n",scaling[i*nfiles+j]);
	      chanVals[i]+=(scaling[i*nfiles+j]*readVal);
	    }
	} // GH added
      digitiseValues(chanVals,outputHeader->nchan,nbits,chanDigitised+iblk*outputHeader->nchan/nsampleByte,levelSetValue);

      if (debugOut == 1 && t >= debugOutT1 && t < debugOutT2 &&
	  i >= debugOutF1 && i < debugOutF2)
	{
	  for (i=0;i<outputHeader->nchan;i++)
	    {	      
	      fprintf(fout2,"%f ",chanVals[i]);
	    }
	}
      //      fwrite(&chanDigitised,sizeof(unsigned char),nchan/nsampleByte,fout);
      if (debugOut == 1 && t >= debugOutT1 && t < debugOutT2)
	fprintf(fout2,"\n");
      iblk++;
      if (iblk == nsblk)
	{
	  int ival;
	  int colnum;
	  long naxes[4];
	  naxes[0] = outputHeader->nchan;
	  naxes[1] = npol;
	  naxes[2] = nsblk;

	  
	  fits_insert_rows(fp,sub,1,&status);
	  if (status) { fits_report_error(stderr,status); exit(1);}
	  fits_get_colnum(fp,CASEINSEN,(char *)"INDEXVAL",&colnum,&status); fits_report_error(stderr,status);
	  ival = sub+1;  fits_write_col(fp,TINT,colnum,sub+1,1,1,&ival,&status);

	  fits_get_colnum(fp,CASEINSEN,(char *)"DAT_FREQ",&colnum,&status); fits_report_error(stderr,status);
	  fits_modify_vector_len(fp,colnum,outputHeader->nchan,&status); fits_report_error(stderr,status);
	  fits_write_col(fp,TFLOAT,colnum,sub+1,1,outputHeader->nchan,datFreq,&status);
	  
	  fits_get_colnum(fp,CASEINSEN,(char *)"DAT_WTS",&colnum,&status); fits_report_error(stderr,status);
	  fits_modify_vector_len(fp,colnum,outputHeader->nchan,&status); fits_report_error(stderr,status);
	  fits_write_col(fp,TFLOAT,colnum,sub+1,1,outputHeader->nchan,datWts,&status);
	  
	  fits_get_colnum(fp,CASEINSEN,(char *)"DAT_OFFS",&colnum,&status); fits_report_error(stderr,status);
	  fits_modify_vector_len(fp,colnum,outputHeader->nchan,&status); fits_report_error(stderr,status);
	  fits_write_col(fp,TFLOAT,colnum,sub+1,1,outputHeader->nchan,datOffs,&status);
	  
	  fits_get_colnum(fp,CASEINSEN,(char *)"DAT_SCL",&colnum,&status); fits_report_error(stderr,status);
	  fits_modify_vector_len(fp,colnum,outputHeader->nchan,&status); fits_report_error(stderr,status);
	  fits_write_col(fp,TFLOAT,colnum,sub+1,1,outputHeader->nchan,datScl,&status);

	  // Write the data
	  fits_get_colnum(fp,CASEINSEN,(char *)"DATA",&colnum,&status); fits_report_error(stderr,status);
	  fits_modify_vector_len (fp, colnum, (outputHeader->nchan*npol*nsblk), &status); fits_report_error(stderr,status);
	  fits_delete_key(fp, (char *)"TDIM18", &status); // THIS SHOULD NOT BE HARDCODED
	  fits_write_tdim(fp,colnum,3,naxes,&status);fits_report_error(stderr,status);

	  fits_write_col_byt(fp,colnum,sub+1,1,nsblk*npol*outputHeader->nchan/nsampleByte,chanDigitised,&status);
	  if (status) { fits_report_error(stderr,status); exit(1);}

	  iblk=0;
	  sub++;
	}
    }
  fits_close_file(fp,&status);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }

  if (debugOut == 1)
    fclose(fout2);

  for (i=0;i<nfiles;i++)
    {
      fclose(fin[i]);
      simulateReleaseMemory(head[i]);
      free(head[i]);
    }
  simulateReleaseMemory(outputHeader);
  free(head);
  free(outputHeader);
  free(scaling);
  free(chanVals);
  free(chanDigitised);
  free(datFreq);
  free(datWts);
  free(datOffs);
  free(datScl);
}


// Digitise floating values into a small number of bits (currently only 1-bit data has been implemented
//
void digitiseValues(float *chanVals,int nchan,int nbits,unsigned char *chanDigitised,float level)
{
  int nsamplesByte = 8/nbits;
  unsigned char tc;
  double bit_level=0;
  int i,j;
  
  if (nbits==1)
    {
      long n=0;
      for (i=0;i<nchan/nsamplesByte;i++)
        {
          tc=0;
          for (j=0;j<8;j++)
            {
              if (chanVals[n] > bit_level)
                tc = tc | (1 << (7-j));
              n++;
            }
          chanDigitised[i] = tc;
        }
    }
  else if (nbits == 2)
    {
      // Set the levels
      unsigned char inArray[nchan];
      for (i=0;i<nchan;i++)
	{
	  if (chanVals[i] < -80)
	    inArray[i] = 0;
	  else if (chanVals[i] < 0)
	    inArray[i] = 1;
	  else if (chanVals[i] < 80)
	    inArray[i] = 2;
	  else
	    inArray[i] = 3;
			   
	}
      twoBitDigitise(inArray, chanDigitised, nchan);
    }
  else if (nbits == 8)
    {
      for (i=0;i<nchan;i++)
	{
	  chanDigitised[i] = (unsigned char)(chanVals[i]*level+128); 
	  //	  if (i==0) printf("Values are: %g %d\n",chanVals[i],chanDigitised[i]);
	}
    }
  else
    {
      printf("ERROR: Do not know how to output %d bits\n",nbits);
      exit(1);
    }
}

int twoBitDigitise(unsigned char *inArray, unsigned char *outArray, unsigned int n)
{
  /* take a list of chars in inArray which only have the bottom two */
  /* bits set and pack them into outArray with 4 * 2 bits in each   */
  /* char. n is the number of 2 bit chars in inArray.               */
  /* This routine returns 0 if all is OK, otherwise -1              */
  
  unsigned char x = 0;
  unsigned int  i;
  unsigned char *outPtr = outArray;
  unsigned char *inPtr = inArray;
  
  
  for (i = 0; i < n; i++)
    {
      if ((*inPtr & 0x03) != *inPtr)
	{
	  return -1;
	}
      
      x = (x << 2) + *inPtr;
      inPtr++;
      if (((i + 1) % 4) == 0)
	{
	  *outPtr = x;
	  outPtr++;
	  x = 0;
	}
    }
  *outPtr = x;
  return 0;
}


void createPrimaryHeader(fitsfile *fp,float f1,float f2,header *head)
{
  int status=0;
  char str[128];
  float freqc = (f1+f2)/2.0;
  float bw = fabs(f2-f1);
  int stt_imjd=head->imjd;
  float stt_smjd=head->smjd;
  float stt_offs=head->stt_offs;
  int nchan = head->nchan;
  char ra[128] = "17:44:00";
  char dec[128] = "07:47:00";
  float ftemp;
  
  fits_movabs_hdu(fp,1,NULL,&status);
  fits_write_date(fp,&status);
  
  fits_update_key(fp,TFLOAT,"OBSFREQ",&freqc,NULL,&status);
  fits_update_key(fp,TFLOAT,"OBSBW",&bw,NULL,&status);
  fits_update_key(fp,TINT,"OBSNCHAN",&nchan,NULL,&status);
  fits_update_key(fp,TINT,"STT_IMJD",&stt_imjd,NULL,&status);
  fits_update_key(fp,TFLOAT,"STT_SMJD",&stt_smjd,NULL,&status);
  fits_update_key(fp,TFLOAT,"STT_OFFS",&stt_offs,NULL,&status);
  fits_update_key(fp, TSTRING, (char *)"RA", &ra, NULL, &status); 
  fits_update_key(fp, TSTRING, (char *)"DEC", &dec, NULL, &status);
  fits_update_key(fp, TSTRING, (char *)"OBSERVER", (char *)"Zhang Lei", NULL, &status);
  fits_update_key(fp, TSTRING, (char *)"PROJID", (char *)"P888", NULL, &status);
  fits_update_key(fp, TSTRING, (char *)"TELESCOP", (char *)"Parkes", NULL, &status);
  fits_update_key(fp, TSTRING, (char *)"FD_POLN", (char *)"LIN", NULL, &status);
  fits_update_key(fp, TSTRING, (char *)"FRONTEND", (char *)"Frontend", NULL, &status);
  fits_update_key(fp, TSTRING, (char *)"BACKEND", (char *)"Backend", NULL, &status);
  fits_update_key(fp, TSTRING, (char *)"SRC_NAME", (char *)"J1713+0747", NULL, &status);
  fits_update_key(fp, TSTRING, (char *)"OBS_MODE", (char *)"SEARCH", NULL, &status);
  fits_update_key(fp, TSTRING, (char *)"DATE-OBS", (char *)"2016-10-10T03:23:32", NULL, &status);
  fits_update_key(fp, TSTRING, (char *)"TRK_MODE", (char *)"TRACK", NULL, &status);
  ftemp = 0.1; fits_update_key(fp,TFLOAT,"BMAJ",&ftemp,NULL,&status);
  ftemp = 0.1; fits_update_key(fp,TFLOAT,"BMIN",&ftemp,NULL,&status);
  ftemp = 0.0; fits_update_key(fp,TFLOAT,"BPA",&ftemp,NULL,&status);
  ftemp = 1; fits_update_key(fp,TFLOAT,"FD_HAND",&ftemp,NULL,&status);
  ftemp = 0.0; fits_update_key(fp,TFLOAT,"FD_SANG",&ftemp,NULL,&status);
  ftemp = 0.0; fits_update_key(fp,TFLOAT,"FD_XYPH",&ftemp,NULL,&status);
  ftemp = 0.0; fits_update_key(fp,TFLOAT,"CHAN_DM",&ftemp,NULL,&status);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }

}

void deleteUnwantedTables(fitsfile *fp)
{
  int status=0;
  int newHDUtype;
  
  fits_movnam_hdu(fp, BINARY_TBL, (char *)"BANDPASS", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, (char *)"COHDDISP", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, (char *)"PSRPARAM", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, (char *)"POLYCO", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, (char *)"T2PREDICT", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, (char *)"FLUX_CAL", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, (char *)"CAL_POLN", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, (char *)"FEEDPAR", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
}

void createSubintHeader(fitsfile *fp,float f1,float f2,int nchan,float tsamp,int nbits)
{
  int status=0;
  char name[128];
  char str[128];
  char cval[128];
  int npol=1;
  int nsblk=4096;
  float chanbw = (f2-f1)/nchan;
  int itemp;
  
  strcpy(name,"SUBINT");
  fits_movnam_hdu(fp,BINARY_TBL,name,0,&status);fits_report_error(stderr,status);

  if (nbits==1)
    {
      sprintf(cval,"16X");//,nchan*nsblk);
      fits_update_key(fp, TSTRING, (char *)"TFORM18", cval, NULL, &status); // Was 20
      fits_report_error(stderr,status);
    }
  else
    {
      //            sprintf(cval,"%dB",nchan*nsblk);
      sprintf(cval,"2B");
      fits_update_key(fp, TSTRING, (char *)"TFORM18", cval, NULL, &status); // Was 20
      fits_report_error(stderr,status);
    }
  fits_update_key_str(fp,(char *)"INT_TYPE",(char *)"TIME",(char *)"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,(char *)"INT_UNIT",(char *)"SEC",(char *)"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,(char *)"SCALE",(char *)"FluxDen",(char *)"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,(char *)"NBIN",(char *)"1",(char *)"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,(char *)"ZERO_OFF",(char *)"0",(char *)"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,(char *)"SIGNINT",(char *)"0",(char *)"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,(char *)"NSUBOFFS",(char *)"0",(char *)"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,(char *)"POL_TYPE",(char *)"AA+BB",(char *)"",&status);  fits_report_error(stderr,status);


  fits_update_key(fp,TINT,"NPOL",&npol,NULL,&status);
  fits_update_key(fp,TINT,"NBITS",&nbits,NULL,&status);
  fits_update_key(fp,TINT,"NSBLK",&nsblk,NULL,&status);
  fits_update_key(fp,TINT,"NCHAN",&nchan,NULL,&status);
  fits_update_key(fp,TFLOAT,"TBIN",&tsamp,NULL,&status);
  fits_update_key(fp,TFLOAT,"CHAN_BW",&chanbw,NULL,&status);
  itemp = 0; fits_update_key(fp,TINT,"NCHNOFFS",&itemp,NULL,&status);
}

double psrangle(double centre_long,double centre_lat,double psr_long,double psr_lat)
{
  double dlon,dlat,a,c;
  double deg2rad = M_PI/180.0;

  /* Apply the Haversine formula */
  dlon = (psr_long - centre_long);
  dlat = (psr_lat  - centre_lat);
  a = pow(sin(dlat/2.0),2) + cos(centre_lat) *
    cos(psr_lat)*pow(sin(dlon/2.0),2);
  if (a==1)
    c = M_PI/deg2rad;
  else
    c = 2.0 * atan2(sqrt(a),sqrt(1.0-a))/deg2rad;
  return c;
}


// Calculates the reduction in gain caused by the source being offset from the beam pointing direction
// The input angle is in degrees, the diameter is the telescope diameter (m) and freq is the observing frequency (MHz).  The output value is the scaling factor
//
double beamScaling(double angle,float diameter,float freq)
{
  double radangle = angle*M_PI/180.0;
  double tt;
  double lambda = 3.0e8/(freq*1.0e6);
  double ang;
  ang = 1.22*lambda/diameter;
  
  tt = radangle*M_PI/ang;
  if (tt == 0)
    return 1;
  else
    return pow(sin(tt)/tt,2);
}

double mjd2year(double mjd,int *retDay,int *retMonth,int *retYr)
{
  double jd,fjd,day;
  int ijd,b,c,d,e,g,month,year;
  //  int retYr,retDay,stat;
  int stat;
  
  jd = mjd + 2400000.5;
  ijd = (int)(jd+0.5);
  fjd = (jd+0.5)-ijd;
  if (ijd > 2299160)
    {
      int a;
      a = (int)((ijd-1867216.25)/36524.25);
      b = ijd + 1 + a - (int)(a/4.0);
    }
  else
    b = ijd;
  
  c = b + 1524;
  d = (int)((c - 122.1)/365.25);
  e = (int)(365.25*d);
  g = (int)((c-e)/30.6001);
  day = c-e+fjd-(int)(30.6001*g);
  if (g<13.5)
    month = g-1;
  else
    month = g-13;
  if (month>2.5)
    year = d-4716;
  else
    year = d-4715;
  //  slaCalyd(year, month, (int)day, retYr, retDay, &stat);
  *retDay = (int)day;
  *retYr = year;
  *retMonth = month;
  
  //  return  retYr+(retDay+(day-(int)day))/365.25;
}
