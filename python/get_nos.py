from __future__ import division
import numpy as np
from pylab import *
import pywt
from scipy import interpolate
import pandas as pd
import copy as cp

output_nos=[]
output_bandpass=[]
def smooth(sig,threshold = 3, level=6, wavelet='db8'):
        sigma = sig.std()
        dwtmatr = pywt.wavedec(data=sig, wavelet=wavelet, level=level)
        denoised = dwtmatr[:]
        denoised[1:] = [pywt.threshold(i, value=threshold*sigma, mode='soft') for i in dwtmatr[1:]]
        smoothed_sig = pywt.waverec(denoised, wavelet, mode='sp1')[:sig.size]
        noises = sig - smoothed_sig
        return smoothed_sig, noises

bandpass = np.loadtxt("FP20180126_0-1GHz_Dec+46.6_drifting_0804.txt")
print bandpass
idxbad_chan=[]
inf=float('inf')
print "bandpass.shape",bandpass.shape
print "length(bandpass)",len(bandpass)
idxgood_num_=[]
idxbad_chan_=[]
copy=cp.deepcopy(bandpass)
print copy
sig,nos = smooth(copy)
print sig.shape
print nos.shape
output_nos.append(nos)
idxarr = np.arange(copy.size)
tck = interpolate.splrep(idxarr,sig)
xnew=np.arange(1,4097,1)
y = interpolate.splev(xnew, tck)
y2 = abs((y-copy)/y)
y2[abs(y2)==inf]=0
idxgood = idxarr[y2<0.3]
idxbad = idxarr[y2>=0.3]
print idxgood,idxbad
a=idxgood.size
idxgood_num_.append(a)
copy[idxbad]=1
idxbad_chan_.append(copy)

print copy
print bandpass
#plot(bandpass)
#plot(idxgood,bandpass[idxgood],'.')
plot(sig,'.')
plot(nos,'.')
show()
