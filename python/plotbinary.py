import numpy as np 
import fitsio
import os
import datetime
import time
import sys
from array import array
from pylab import *
import astropy.io.fits as pyfits
from scipy import interpolate




if (len(sys.argv)<1):
  print 'too few inputs!'
  sys.exit()
else:
  print 'input seems OK'



print 'record start time:'
starttime=datetime.datetime.now()
print starttime


rowdatafile=sys.argv[1]



nsblk=4096
nchan=4096

simdata=np.zeros((64,nsblk,1,nchan,1))

rowdata=np.fromfile(rowdatafile,dtype=np.float32,count=-1)

print "rowdata.shape:",rowdata.shape
print "simdata.shape",simdata.shape
print "simdata.shape",simdata[0,0,0,:,0].shape

temp=np.zeros(64*nsblk)

print "start reshape file",datetime.datetime.now()
simdata = np.fromfile(rowdatafile,dtype=np.float32,count=-1).reshape((64,nsblk,1,nchan,1),order='C')
print "end reshape file",datetime.datetime.now()

for i in range(64):
    temp[i*4096:(i+1)*4096]=simdata[i,:,0,1000,0]
#temp = simdata[1,:,0,1000,0]

print np.max(temp)
plot(temp)
show()


#temp1=np.zeros(64*nsblk)
#print "start reshape file",datetime.datetime.now()
#simdata=np.zeros((64,nsblk,1,nchan,1))
#
#rowdata=np.fromfile(rowdatafile,dtype=np.float32,count=-1)
#
#print "rowdata.shape:",rowdata.shape
#print "simdata.shape",simdata.shape
#print "simdata.shape",simdata[0,0,0,:,0].shape
#for i in range(64):
#    for j in range(nsblk):
#        simdata[i,j,0,:,0]=rowdata[(i*4096+j)*nchan:(i*4096+j+1)*nchan]
#
#print "end reshape file",datetime.datetime.now()
#
#
#for i in range(64*nsblk):
#    temp1[i]=rowdata[i*nchan+1000]-temp[i]
##
##
#plot(temp1)
###plot(data.sum(axis=0))
###plot(data.sum(axis=1))
#show()


