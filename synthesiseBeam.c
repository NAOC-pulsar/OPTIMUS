#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

int main(int argc,char *argv[])
{
  FILE *fin;
  FILE *fout;
  int i;
  float *vals1;
  float *vals2;
  float *vals3;
  float *timeVals;
  
  long n=0;

  long maxVals = 600000;
  float synthBeam;
  double gain1,gain2,gain3;
  
  if (!(vals1 = (float *)malloc(sizeof(float)*maxVals)))
    {
      printf("Unable to allocate memory\n");
      exit(1);
    }
  if (!(vals2 = (float *)malloc(sizeof(float)*maxVals)))
    {
      printf("Unable to allocate memory\n");
      exit(1);
    }
  if (!(vals3 = (float *)malloc(sizeof(float)*maxVals)))
    {
      printf("Unable to allocate memory\n");
      exit(1);
    }
  if (!(timeVals = (float *)malloc(sizeof(float)*maxVals)))
    {
      printf("Unable to allocate memory\n");
      exit(1);
    }

  // Read beam1
  fin = fopen(argv[1],"r");
  n=0;
  while (!feof(fin))
    {
      if (fscanf(fin,"%f %f",&timeVals[n],&vals1[n])==2)
	n++;
    }
  fclose(fin);

  // Read beam2
  fin = fopen(argv[2],"r");
  n=0;
  while (!feof(fin))
    {
      if (fscanf(fin,"%f %f",&timeVals[n],&vals2[n])==2)
	n++;
    }
  fclose(fin);

  // Read beam3
  fin = fopen(argv[3],"r");
  n=0;
  while (!feof(fin))
    {
      if (fscanf(fin,"%f %f",&timeVals[n],&vals3[n])==2)
	n++;
    }
  fclose(fin);

  // Form beam and write out
  fout = fopen("synthBeam.dat","w");
  for (i=0;i<n;i++)
    {
      gain1 = pow(sin((timeVals[i]-55)/4.)/((timeVals[i]-55)/4.),2);
      gain2 = pow(sin((timeVals[i]-83)/4.)/((timeVals[i]-83)/4.),2);
      gain3 = pow(sin((timeVals[i]-110)/4.)/((timeVals[i]-110)/4.),2);
      synthBeam = gain1*vals1[i]+gain2*vals2[i]+gain3*vals3[i];
      fprintf(fout,"%g %g\n",timeVals[i],synthBeam);
    }
  fclose(fout);
 
  
  free(timeVals);
  free(vals1);
  free(vals2);
  free(vals3);

}
