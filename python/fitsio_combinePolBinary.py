import numpy as np 
import fitsio
import os
import datetime
import time
import sys
from array import array
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import *
import astropy.io.fits as pyfits
from decimal import Decimal
import pywt
#from scipy import interpolate
#import pandas as pd

# get bandpass from fits file
def get_bandpass(filename):
    secperday = 3600 * 24
    infile=sys.argv[1]
#    print filename
    try:
        hdulist = pyfits.open(filename)
        hdu0 = hdulist[0]
        hdu1 = hdulist[1]
        data1 = hdu1.data['data']
        header1 = hdu1.header
        fchannel = hdulist['SUBINT'].data[0]['DAT_FREQ']
        fch1 = fchannel[0]
        obsfreq = hdu0.header['OBSFREQ']
        obsnchan = hdu0.header['OBSNCHAN']
        obsbw = hdu0.header['OBSBW']
        fmin = obsfreq - obsbw/2.
        fmax = obsfreq + obsbw/2.
        nf = obsnchan
        df = hdu1.header['CHAN_BW']
        tsamp = hdu1.header['TBIN']
        nsubint = hdu1.header['NAXIS2']
        samppersubint = int(hdu1.header['NSBLK'])
        nsamp = nsubint * samppersubint
        sourename = hdu0.header['SRC_NAME']
        ra = hdu0.header['RA']
        dec = hdu0.header['DEC']
        subintoffset = hdu1.header['NSUBOFFS']
        tstart = "%.13f" % (Decimal(hdu0.header['STT_IMJD']) + Decimal(hdu0.header['STT_SMJD'] + tsamp * samppersubint * subintoffset )/secperday )
        nbits = hdu0.header['BITPIX']
        header = hdu0.header + hdu1.header
        dtype = ''
        a,b,c,d,e = data1.shape
        if c > 1:
            data = data1[:,:,0,:,:].squeeze().reshape((-1,d))
        else:
            data = data1.squeeze().reshape((-1,d))
        l, m = data.shape
        data = data.reshape(l/64, 64, d).sum(axis=1)
        data1 = np.sum(data,axis=0)
    except:
        print ('Error')
    else:
        return data1


#smooth the bandpass
#need change#
def smooth_bandpass(bandpass):
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
#    print "bandpass.shape",bandpass.shape 
    idxbad_chan=[]
    inf=float('inf')
#    print "bandpass.shape",bandpass.shape
#    print "length(bandpass)",len(bandpass)
    idxgood_num_=[]
    idxbad_chan_=[]
    sig,nos = smooth(bandpass)
    return sig




######################################################################
# 2018/04/04 #  mcc       output: *.fits (pol averaged obs data + binary data)
######################################################################

#get input files
if (len(sys.argv)<3):
  print 'too few inputs!'
  print 'example:'
  print 'usage: python fits_combinePolBinary.py fitsfile binaryfile combinefile'
  sys.exit()
else:
  print 'input seems OK'



print 'record start time:'
starttime=datetime.datetime.now()
print starttime


infile=sys.argv[1]
rowdatafile=sys.argv[2]
outfile=sys.argv[3]

#u19700101=62135683200.0

#==============================================================
#get data in fits file
fits=fitsio.FITS(infile)

hdu0 = fits[0]
header0 = hdu0.read_header()

hdu1 = fits[1]
header1 = hdu1.read_header()

hdrver=header0['HDRVER']
date=header0['DATE']
ant_x=header0['ANT_X']
ant_y=header0['ANT_Y']
ant_z=header0['ANT_Z']
obsfreq=header0['OBSFREQ']
obsbw=header0['OBSBW']
nchan=header0['OBSNCHAN']
ra=header0['RA']
dec=header0['DEC']
bmaj=header0['BMAJ']
bmin=header0['BMIN']
date_obs=header0['DATE-OBS']
stt_imjd=header0['STT_IMJD']
stt_smjd=header0['STT_SMJD']
stt_offs=header0['STT_OFFS']
stt_lst=header0['STT_LST']
nsuboffs=header1['NSUBOFFS']
nchnoffs=header1['NCHNOFFS']
nsblk=header1['NSBLK']
npol=header1['NPOL']
tbin=header1['TBIN']
chan_bw=header1['CHAN_BW']
print 'NSBLK: ',nsblk,'sample time(s): ',tbin,'channel width(MHz): ',chan_bw

nline=header1['NAXIS2']
pol_type=header1['POL_TYPE']



rmcomm='rm -f '+outfile
os.system(rmcomm)
print rmcomm,outfile
fitsout=fitsio.FITS(outfile,'rw')

dataout = np.zeros(1,dtype=[('TSUBINT','float64'),('OFFS_SUB','float64'),('LST_SUB','float64'),('RA_SUB','float64'),('DEC_SUB','float64'),('GLON_SUB','float64'),('GLAT_SUB','float64'),('FD_ANG','float32'),('POS_ANG','float32'),('PAR_ANG','float32'),('TEL_AZ','float32'),('TEL_ZEN','float32'),('DAT_FREQ','float32',(nchan)),('DAT_WTS','float32',(nchan)),('DAT_OFFS','float32',(nchan)),('DAT_SCL','float32',(nchan)),('DATA','uint8',(nsblk,1,nchan,1))])


#==============================================================
#get bandpass and smooth it
bandpass=get_bandpass(infile)
smoothBandpass=smooth_bandpass(bandpass)
#smoothBandpass=smoothBandpass/max(smoothBandpass)
smoothBandpass=smoothBandpass/(52.4288/0.0002)
print "normalized ",'max(smoothBandpass):',np.max(smoothBandpass),'min(smoothBandpass):',np.min(smoothBandpass)
#plot(bandpass,',')
#plot(smoothBandpass,',')
#show()



#==============================================================
# read binanry data and use smoothed bandpass changing the profile in each channel
#binary recorded from float32 to uint8
simdata=np.zeros((64,nsblk,1,nchan,1))
rowdata=np.fromfile(rowdatafile,dtype=np.float32,count=-1)

#print "binarydata.shape:",rowdata.shape
#print "simdata.shape",simdata.shape

for i in range(64):
    for j in range(nsblk):
        simdata[i,j,0,:,0]=rowdata[(i*4096+j)*nchan:(i*4096+j+1)*nchan]

print 'simdata.dtype',simdata.dtype,'simadata.max',np.max(simdata),'simdata.min',np.min(simdata)

for i in range(nchan):
    simdata[:,:,0,i,0]=simdata[:,:,0,i,0]*smoothBandpass[i]

print 'simdata.dtype',simdata.dtype,'simadata.max',np.max(simdata),'simdata.min',np.min(simdata)
simdata[simdata > 127.] = 127.
simdata[simdata < 0.] = 0.

simdata = simdata.astype(uint8)
print 'simdata.dtype',simdata.dtype,'simadata.max',np.max(simdata),'simdata.min',np.min(simdata)

#simdata = np.fromfile(rowdatafile,dtype=np.float32,count=-1).reshape((64,nsblk,1,nchan,1))




#=======================================================================
rowindex=0
dataout['TSUBINT'][0]=fits[1].read(rows=[rowindex], columns=['TSUBINT'])[0][0]
dataout['OFFS_SUB'][0]=fits[1].read(rows=[rowindex], columns=['OFFS_SUB'])[0][0]
dataout['LST_SUB'][0]=fits[1].read(rows=[rowindex], columns=['LST_SUB'])[0][0]
dataout['RA_SUB'][0]=fits[1].read(rows=[rowindex], columns=['RA_SUB'])[0][0]
dataout['DEC_SUB'][0]=fits[1].read(rows=[rowindex], columns=['DEC_SUB'])[0][0]
dataout['GLON_SUB'][0]=fits[1].read(rows=[rowindex], columns=['GLON_SUB'])[0][0]
dataout['GLAT_SUB'][0]=fits[1].read(rows=[rowindex], columns=['GLAT_SUB'])[0][0]
dataout['FD_ANG'][0]=fits[1].read(rows=[rowindex], columns=['FD_ANG'])[0][0]
dataout['POS_ANG'][0]=fits[1].read(rows=[rowindex], columns=['POS_ANG'])[0][0]
dataout['PAR_ANG'][0]=fits[1].read(rows=[rowindex], columns=['PAR_ANG'])[0][0]
dataout['TEL_AZ'][0]=fits[1].read(rows=[rowindex], columns=['TEL_AZ'])[0][0]
dataout['TEL_ZEN'][0]=fits[1].read(rows=[rowindex], columns=['TEL_ZEN'])[0][0]
dataout['DAT_FREQ'][0]=fits[1].read(rows=[rowindex], columns=['DAT_FREQ'])[0][0][0:nchan]
dataout['DAT_WTS'][0]=fits[1].read(rows=[rowindex], columns=['DAT_WTS'])[0][0][0:nchan]
dataout['DAT_OFFS'][0][0:nchan]=fits[1].read(rows=[rowindex], columns=['DAT_OFFS'])[0][0][0:nchan]
dataout['DAT_SCL'][0][0:nchan]=fits[1].read(rows=[rowindex], columns=['DAT_SCL'])[0][0][0:nchan]

data=fits[1].read(rows=[rowindex], columns=['DATA'])
print 'data.dtype',data[0][0].dtype,'data.max',np.max(data[0][0]),'data.min',np.min(data[0][0])
#==============================================================
#for subint 1 : add the binary data and real obs data
for subindex in range(nsblk):
    #temp = (data[0][0][subindex,0,:,0]+data[0][0][subindex,1,:,0])/4
    temp = data[0][0][subindex,0,:,0]/2
    temp[temp> 127] = 127
    temp[temp< 0 ] = 0
    dataout['DATA'][0][subindex,0,:,0] = temp +simdata[rowindex,subindex,0,:,0]
fitsout.write(dataout)
#=======================================================================

fitsout[0].write_key('HDRVER',hdrver,comment="")
fitsout[0].write_key('FITSTYPE','PSRFITS',comment="FITS definition ")
fitsout[0].write_key('DATE',date,comment="")
fitsout[0].write_key('OBSERVER','FAST_TEAM',comment="Observer name")
fitsout[0].write_key('PROJID','Drift',comment="Project name")
fitsout[0].write_key('TELESCOP','FAST',comment="Telescope name")
fitsout[0].write_key('ANT_X',ant_x,comment="")
fitsout[0].write_key('ANT_Y',ant_y,comment="")
fitsout[0].write_key('ANT_Z',ant_z,comment="")
fitsout[0].write_key('FRONTEND','WIDEBAND',comment="Frontend ID")
fitsout[0].write_key('NRCVR',1,comment="")
fitsout[0].write_key('FD_POLN','LIN',comment="LIN or CIRC")
fitsout[0].write_key('FD_HAND',1,comment="")
fitsout[0].write_key('FD_SANG',0.,comment="")
fitsout[0].write_key('FD_XYPH',0.,comment="")
fitsout[0].write_key('BACKEND','ROACH',comment="Backend ID")
fitsout[0].write_key('BECONFIG','N/A',comment="")
fitsout[0].write_key('BE_PHASE',1,comment="")
fitsout[0].write_key('BE_DCC',0,comment="")
fitsout[0].write_key('BE_DELAY',0.,comment="")
fitsout[0].write_key('TCYCLE',0.,comment="")
fitsout[0].write_key('OBS_MODE','SEARCH',comment="(PSR, CAL, SEARCH)")
fitsout[0].write_key('DATE-OBS',date_obs,comment="Date of observation")
fitsout[0].write_key('OBSFREQ',obsfreq,comment="[MHz] Bandfrequency")
fitsout[0].write_key('OBSBW',obsbw,comment="[MHz] Bandwidth")
fitsout[0].write_key('OBSNCHAN',nchan,comment="Number of channels")
fitsout[0].write_key('CHAN_DM',0.,comment="")
fitsout[0].write_key('SRC_NAME','Drift',comment="Source or scan ID")
fitsout[0].write_key('COORD_MD','J2000',comment="")
fitsout[0].write_key('EQUINOX',2000.,comment="")

fitsout[0].write_key('RA',ra,comment="")
fitsout[0].write_key('DEC',dec,comment="")
fitsout[0].write_key('BMAJ',bmaj,comment="[deg] Beam major axis length")
fitsout[0].write_key('BMIN',bmin,comment="[deg] Beam minor axis length")
fitsout[0].write_key('BPA',0.,comment="[deg] Beam position angle")
fitsout[0].write_key('STT_CRD1','00:00:00.00',comment="")
fitsout[0].write_key('STT_CRD2','00:00:00.00',comment="")
fitsout[0].write_key('TRK_MODE','TRACK',comment="")
fitsout[0].write_key('STP_CRD1','00:00:00.00',comment="")
fitsout[0].write_key('STP_CRD2','00:00:00.00',comment="")
fitsout[0].write_key('SCANLEN',0.,comment="")
fitsout[0].write_key('FD_MODE','FA',comment="")
fitsout[0].write_key('FA_REQ',0.,comment="")
fitsout[0].write_key('CAL_MODE','OFF',comment="")
fitsout[0].write_key('CAL_FREQ',0.,comment="")
fitsout[0].write_key('CAL_DCYC',0.,comment="")
fitsout[0].write_key('CAL_PHS',0.,comment="")
fitsout[0].write_key('STT_IMJD',stt_imjd,comment="Start MJD (UTC days) (J - long integer)")
fitsout[0].write_key('STT_SMJD',stt_smjd,comment="[s] Start time (sec past UTC 00h) (J)")
fitsout[0].write_key('STT_OFFS',stt_offs,comment="[s] Start time offset (D)")
fitsout[0].write_key('STT_LST',stt_lst,comment="[s] Start LST (D)")

fitsout[1].write_key('INT_TYPE','TIME',comment="Time axis (TIME, BINPHSPERI, BINLNGASC, etc)")
fitsout[1].write_key('INT_UNIT','SEC',comment="Unit of time axis (SEC, PHS (0-1),DEG)")
fitsout[1].write_key('SCALE','FluxDen',comment="")
fitsout[1].write_key('NPOL',1,comment="Nr of polarisations")
fitsout[1].write_key('POL_TYPE','AABB',comment="Polarisation identifier")
fitsout[1].write_key('TBIN',tbin,comment="[s] Time per bin or sample")
fitsout[1].write_key('NBIN',1,comment="")
fitsout[1].write_key('NBIN_PRD',0,comment="Nr of bins/pulse period (for gated data)")
fitsout[1].write_key('PHS_OFFS',0.0,comment="Phase offset of bin 0 for gated data")
fitsout[1].write_key('NBITS',8,comment="Nr of bits/datum ")
fitsout[1].write_key('NSUBOFFS',nsuboffs,comment="Subint offset ")
fitsout[1].write_key('NCHNOFFS',nchnoffs,comment="Channel/sub-band offset for split files")
fitsout[1].write_key('NCHAN',nchan,comment="Number of channels")
fitsout[1].write_key('CHAN_BW',chan_bw,comment="[MHz] Channel/sub-band width")
fitsout[1].write_key('NSBLK',nsblk,comment="Samples/row ")
fitsout[1].write_key('EXTNAME','SUBINT  ',comment="name of this binary table extension")
fitsout[1].write_key('EXTVER',1,comment="")


#==============================================================
#for subint 2-64 : add the binary data and real obs data
for rowindex in range(1,nline):
    dataout['TSUBINT'][0]=fits[1].read(rows=[rowindex], columns=['TSUBINT'])[0][0]
    dataout['OFFS_SUB'][0]=fits[1].read(rows=[rowindex], columns=['OFFS_SUB'])[0][0]
    dataout['LST_SUB'][0]=fits[1].read(rows=[rowindex], columns=['LST_SUB'])[0][0]
    dataout['RA_SUB'][0]=fits[1].read(rows=[rowindex], columns=['RA_SUB'])[0][0]
    dataout['DEC_SUB'][0]=fits[1].read(rows=[rowindex], columns=['DEC_SUB'])[0][0]
    dataout['GLON_SUB'][0]=fits[1].read(rows=[rowindex], columns=['GLON_SUB'])[0][0]
    dataout['GLAT_SUB'][0]=fits[1].read(rows=[rowindex], columns=['GLAT_SUB'])[0][0]
    dataout['FD_ANG'][0]=fits[1].read(rows=[rowindex], columns=['FD_ANG'])[0][0]
    dataout['POS_ANG'][0]=fits[1].read(rows=[rowindex], columns=['POS_ANG'])[0][0]
    dataout['PAR_ANG'][0]=fits[1].read(rows=[rowindex], columns=['PAR_ANG'])[0][0]
    dataout['TEL_AZ'][0]=fits[1].read(rows=[rowindex], columns=['TEL_AZ'])[0][0]
    dataout['TEL_ZEN'][0]=fits[1].read(rows=[rowindex], columns=['TEL_ZEN'])[0][0]
    dataout['DAT_FREQ'][0]=fits[1].read(rows=[rowindex], columns=['DAT_FREQ'])[0][0][0:nchan]
    dataout['DAT_WTS'][0]=fits[1].read(rows=[rowindex], columns=['DAT_WTS'])[0][0][0:nchan]
    dataout['DAT_OFFS'][0][0:nchan]=fits[1].read(rows=[rowindex], columns=['DAT_OFFS'])[0][0][0:nchan]
    dataout['DAT_SCL'][0][0:nchan]=fits[1].read(rows=[rowindex], columns=['DAT_SCL'])[0][0][0:nchan]

    data=fits[1].read(rows=[rowindex], columns=['DATA'])
    for subindex in range(nsblk):
        #temp = (data[0][0][subindex,0,:,0]+data[0][0][subindex,1,:,0])/4 
        temp = data[0][0][subindex,0,:,0]/2 
        temp[temp> 127] = 127
        temp[temp< 0 ] = 0
        dataout['DATA'][0][subindex,0,:,0] = temp+simdata[rowindex,subindex,0,:,0]
    fitsout[-1].append(dataout)

#fitsout.write(dataout)
#==============================================================


fitsout.close()
print '--------------------------------------------'
print '             Finished!                      '


endtime=datetime.datetime.now()
print 'START:',starttime
print 'END:',endtime
duration=endtime-starttime
print 'DURATION:',duration.seconds,' sec'
