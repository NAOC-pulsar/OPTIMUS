#!/usr/bin/python
import numpy as np
import fitsio
import sys, os

"""
usage: %s filename freqlow freqhigh
""" % sys.argv[0]

filename = sys.argv[1]
freqlow = float(sys.argv[2])
freqhigh = float(sys.argv[3])
outname = sys.argv[4]

#read fits file
data0, hdr0 = fitsio.read(filename, ext=0, header=True)

#print hdr0['OBSNCHAN'], hdr0['OBSFREQ'], hdr0['OBSBW']
print 'infile %s\nfrequency range %s, %s \noutfile %s' % (filename, freqlow, freqhigh, outname)

data, hdr = fitsio.read(filename, ext=1, header=True)

freq = data['DAT_FREQ'][0]
idx = np.nonzero(np.logical_and((freq > freqlow), (freq < freqhigh)))[0]

idxlow = min(idx)
idxhigh = max(idx)
newfreq = freq[idx]
#print idxlow, idxhigh
#print newfreq
nchan = newfreq.size
print 'selected ', nchan, ' number of frequency channels.'
avefreq = newfreq.mean()
bandwidth = newfreq.max() - newfreq.min() + hdr['CHAN_BW']

hdr['NCHAN'] = nchan
nrow = len(data)
dtype = eval(str(data.dtype))
newdtype = []

for dt in dtype:
    #print dt
    ndt = list(dt)
    key = ndt[0]
    if key == 'DAT_FREQ':
        ndt[2] = (nchan,) 
    elif key == 'DAT_WTS':
        ndt[2] = (nchan,) 
    elif key == 'DAT_OFFS':
        ndt[2] = (nchan,) 
    elif key == 'DAT_SCL':
        ndt[2] = (nchan,) 
    elif key == 'DATA':
        ndt[2] = (4096, 1, nchan, 1) 
    newdtype.append(tuple(ndt))

#print newdtype
newdata = np.zeros(nrow, dtype=newdtype)
for ndt in newdtype:
    key = ndt[0]
    if key == 'DAT_FREQ':
        newdata[key] = np.squeeze(data[key][:,idx])
    elif key == 'DAT_WTS': 
        newdata[key] = np.squeeze(data[key][:,idx])
    elif key == 'DAT_OFFS': 
        newdata[key] = np.squeeze(data[key][:,idx])
    elif key == 'DAT_SCL': 
        newdata[key] = np.squeeze(data[key][:,idx])
    elif key == 'DATA':
        print newdata[key].shape, data[key][:, :, :, idx, :].shape
        newdata[key] = data[key][:, :, :, idx, :]
    else:
        newdata[key] = data[key]

'''write to files'''
hdr0['OBSNCHAN'] = nchan
hdr0['OBSFREQ'] = avefreq
hdr0['OBSBW'] = bandwidth
hdr['POL_TYPE'] = 'AA+BB'
fitsio.write(outname, data0, header=hdr0, clobber=True)
fitsio.write(outname, newdata, header=hdr, extname='SUBINT')


