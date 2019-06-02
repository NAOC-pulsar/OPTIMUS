import os
import sys
import time
import math
import pywt
import fitsio
import datetime
import numpy as np 
from pylab import *
from array import array
from decimal import Decimal
import astropy.io.fits as pyfits

import matplotlib as mpl
import matplotlib.pyplot as plt
#from scipy import interpolate
#import pandas as pd

#---------------
#add drifting function
#beam weight

def beamWeight(filename, scale=1, percent=0, flag=1):
    #settings
    C = 299794580.0
    D = 300.0
    #Pi = 3.141592654
    Pi = math.pi
    Omega = 2.*Pi/3600./24.
    Filename = filename
    Scale = scale
    Percent = percent
    Flag = flag
    #read file 
    hdulist = fitsio.FITS(Filename, 'r')
    hdu1 = hdulist[1]
    hdu0 = hdulist[0]
    header1 = hdu1.read_header()
    header0 = hdu0.read_header()
    data1=hdulist[1][:]
    #Time_subint = header0['TBIN']*header1['NSBLK']
    Time_subint = header1['TBIN']
    Freq_chann = header1['CHAN_BW']
    Nsblk=header1['NSBLK']
    #Freq_start = header0['OBSFREQ']-(header0['OBSNCHAN']/2.)*Freq_chann
    print header0['OBSFREQ'], header0['OBSNCHAN'], Freq_chann, Time_subint, Nsblk
    Freq_start = header0['OBSFREQ']
    Nsub = np.array(data1['DAT_WTS']).shape[0]
    Nchan = np.array(data1['DAT_WTS']).shape[1]
    Tau = round(Time_subint * Nsblk * Nsub * Percent)	
    Weight = np.zeros((Nsblk*Nsub,Nchan))
    
    for isub in range(0,Nsblk*Nsub):
        Time = Time_subint * isub
        for ichan in range(0,Nchan):
            Freq = Freq_start + ichan*Freq_chann
            if (Flag == 1):
                Weight[isub,ichan] = 1
            else :
                Eta = -2.*Pi*Omega*(Time-Tau)/(1.22*(C/Freq/1000000.)/D)
                #print Freq, Eta, Time, Tau
                if ((Time-Tau) == 0):
                    Weight[isub,ichan] = Scale*1.
                else:
                    Weight[isub,ichan] = Scale*(math.sin(Eta)/Eta)**2.
#                    Weight[isub,ichan] = Scale*math.cos(Dec*2.*Pi/360.)*(math.sin(Eta)/Eta)**2.
    hdulist.close()
    #Weight = Weight.reshape((Nsub,Nsblk,1,Nchan,1),order='C')
    #print 'DAT_WTS[',int(Nsub*Percent),',',Nchan,']=', data1['DAT_WTS'][int(Nsub*Percent),Nchan-1], '->', Dat_wts[int(Nsub*Percent),Nchan-1]
    #print 'DAT_WTS[',Nsub,',',Nchan,']=', data1['DAT_WTS'][Nsub-1,Nchan-1], '->', Dat_wts[Nsub-1,Nchan-1]
    print Weight.shape, np.min(Weight), np.max(Weight)
    return Weight


if __name__ == '__main__':

    Flag = 0
    Scale = 1
    Percent = 0.6
    infile = 'FP20180104_0-1GHz_Dec-0.5_drifting_0779.fits'
    Weight = beamWeight(infile, scale=1, percent=0, flag=1)
    Weight = Weight.reshape(-1)
    f=open('noBeamWeight.dat','wb')
    f.write(Weight)
    f.close()


