import numpy as np
import fitsio
import sys
import matplotlib.pyplot as plt
import time
from pylab import *
import os
import pandas as pd
import progressbar

if len(sys.argv) < 2:
    print '\033[1;31;47mNo action specified.\033[0m'
    sys.exit()
elif len(sys.argv) == 2:
    if sys.argv[1].startswith('-'):
        option = sys.argv[1][1:]
        # fetch sys.argv[1] but without the first two characters
        if option == 'version': 
            print 'Version 18/05/22 by mcc'
            sys.exit()
        elif option == 'help':  
            print \
    '''
usage:
    python fastFitsCut.py [-f1 startchannel] [-f2 endchannel] [-t1 startsubint] [-t2 endsubint] FAST.fits
    example:
    python fastFitsCut.py -f1 1000 -f2 3000 -t1 0 -t2 31 FAST.fits
    python fastFitsCut.py -f1 1000 -f2 3000 -t1 0 FAST.fits
    (if someone lost, the input will be set to a default num.)
    output files:
    FAST_cut_AA.fits FAST_cut_BB.fits FAST_cut_total.fits
Options include:
    -version : Prints the version number
    -help    : Display this help
    '''
    else :
        print \
    '''\033[1;31;47m
    Wrong input.
    You must given a ferq or subint number.
    \033[0m'''
    sys.exit()
else :
    # read file
    filename=sys.argv[-1]
    fits=fitsio.FITS(filename)
    h0 = fits[0].read_header()
    h1 = fits[1].read_header()
    freq=h0['OBSFREQ']
    nchan_origin=h0['OBSNCHAN']
    widthfreq=h0['OBSBW']
    
    tsample=h1['TBIN']
    nsblk=h1['NSBLK']
    npol=h1['NPOL']
    chan_bw=h1['CHAN_BW']

    #default setting
    startfreq = 0
    endfreq = int(nchan_origin)-1
    startn = 0
    endn = int(h1['NAXIS2']) - 1 # 0 ~ max-1 <==> 1 ~ max

    #set channel num and subint
    for i in range(len(sys.argv)):
        if sys.argv[i].startswith('-'):
            option = sys.argv[i][1:]
            if option == 'f1':
                startfreq=int(sys.argv[i+1])
            if option == 'f2':
                endfreq=int(sys.argv[i+1])
            if option == 't1':
                startn = int(sys.argv[i+1])
            if option == 't2':
                endn = int(sys.argv[i+1])


outname1=filename[0:-5]+'_cut_AA.fits'
outname2=filename[0:-5]+'_cut_BB.fits'
outname3=filename[0:-5]+'_cut_total.fits'

print 'input:\n    startchan: {startfreq}\n    endchan: {endfreq}\n    startsubint: {startn}\n    endsubint: {endn}\n'.format(startfreq = startfreq, endfreq = endfreq, startn = startn, endn = endn)
print 'out put files:\n    %s\n    %s\n    %s' %(outname1,outname2,outname3)

chnum=endfreq-startfreq+1
linenum=endn-startn+1
nrow=linenum

nchan=chnum


specs=np.zeros((nrow*nsblk,nchan))
specs_av=np.zeros((nrow,nchan))


command1='rm -f '+outname1
command2='rm -f '+outname2
command3='rm -f '+outname3
os.system(command1)
os.system(command2)
os.system(command3)



data3 = np.zeros(1,dtype=[('TSUBINT','float64'),('OFFS_SUB','float64'),('LST_SUB','float64'),('RA_SUB','float64'),('DEC_SUB','float64'),('GLON_SUB','float64'),('GLAT_SUB','float64'),('FD_ANG','float32'),('POS_ANG','float32'),('PAR_ANG','float32'),('TEL_AZ','float32'),('TEL_ZEN','float32'),('DAT_FREQ','float32',(nchan)),('DAT_WTS','float32',(nchan)),('DAT_OFFS','float32',(2*nchan)),('DAT_SCL','float32',(2*nchan)),('DATA','uint8',(nsblk,1,nchan,1))])

data4 = np.zeros(1,dtype=[('TSUBINT','float64'),('OFFS_SUB','float64'),('LST_SUB','float64'),('RA_SUB','float64'),('DEC_SUB','float64'),('GLON_SUB','float64'),('GLAT_SUB','float64'),('FD_ANG','float32'),('POS_ANG','float32'),('PAR_ANG','float32'),('TEL_AZ','float32'),('TEL_ZEN','float32'),('DAT_FREQ','float32',(nchan)),('DAT_WTS','float32',(nchan)),('DAT_OFFS','float32',(2*nchan)),('DAT_SCL','float32',(2*nchan)),('DATA','uint8',(nsblk,1,nchan,1))])

data5 = np.zeros(1,dtype=[('TSUBINT','float64'),('OFFS_SUB','float64'),('LST_SUB','float64'),('RA_SUB','float64'),('DEC_SUB','float64'),('GLON_SUB','float64'),('GLAT_SUB','float64'),('FD_ANG','float32'),('POS_ANG','float32'),('PAR_ANG','float32'),('TEL_AZ','float32'),('TEL_ZEN','float32'),('DAT_FREQ','float32',(nchan)),('DAT_WTS','float32',(nchan)),('DAT_OFFS','float32',(2*nchan)),('DAT_SCL','float32',(2*nchan)),('DATA','uint8',(nsblk,1,nchan,1))])

fitsout1=fitsio.FITS(outname1,'rw')
fitsout2=fitsio.FITS(outname2,'rw')
fitsout3=fitsio.FITS(outname3,'rw')



ch=np.array(range(nchan))
nu=(ch*1.0/(nchan-1)-0.5)*widthfreq+freq

bar = progressbar.ProgressBar()

for index in bar(range(nrow)):
    rowindex=index+startn

    bar.update(index+1)
    #print "processing subint: %s" %(index+startn)

    data = fits[1].read(rows=[rowindex], columns=['DATA'])

    tempspec1=data[0][0][:,0,:,0]
    tempspec2=data[0][0][:,1,:,0]
    #print tempspec1.shape, tempspec2.shape
    data3['TSUBINT'][0]=fits[1].read(rows=[rowindex], columns=['TSUBINT'])[0][0]
    data3['OFFS_SUB'][0]=fits[1].read(rows=[rowindex], columns=['OFFS_SUB'])[0][0]
    data3['LST_SUB'][0]=fits[1].read(rows=[rowindex], columns=['LST_SUB'])[0][0]
    data3['RA_SUB'][0]=fits[1].read(rows=[rowindex], columns=['RA_SUB'])[0][0]
    data3['DEC_SUB'][0]=fits[1].read(rows=[rowindex], columns=['DEC_SUB'])[0][0]
    data3['GLON_SUB'][0]=fits[1].read(rows=[rowindex], columns=['GLON_SUB'])[0][0]
    data3['GLAT_SUB'][0]=fits[1].read(rows=[rowindex], columns=['GLAT_SUB'])[0][0]
    data3['FD_ANG'][0]=fits[1].read(rows=[rowindex], columns=['FD_ANG'])[0][0]
    data3['POS_ANG'][0]=fits[1].read(rows=[rowindex], columns=['POS_ANG'])[0][0]
    data3['PAR_ANG'][0]=fits[1].read(rows=[rowindex], columns=['PAR_ANG'])[0][0]
    data3['TEL_AZ'][0]=fits[1].read(rows=[rowindex], columns=['TEL_AZ'])[0][0]
    data3['TEL_ZEN'][0]=fits[1].read(rows=[rowindex], columns=['TEL_ZEN'])[0][0]
    data3['DAT_FREQ'][0]=fits[1].read(rows=[rowindex], columns=['DAT_FREQ'])[0][0][startfreq:endfreq+1]
    data3['DAT_WTS'][0]=fits[1].read(rows=[rowindex], columns=['DAT_WTS'])[0][0][startfreq:endfreq+1]
    
    data3['DAT_OFFS'][0][0:nchan]=fits[1].read(rows=[rowindex], columns=['DAT_OFFS'])[0][0][startfreq:endfreq+1]
    data3['DAT_OFFS'][0][nchan:2*nchan]=fits[1].read(rows=[rowindex], columns=['DAT_OFFS'])[0][0][startfreq+nchan_origin:endfreq+1+nchan_origin]
    data3['DAT_SCL'][0][0:nchan]=fits[1].read(rows=[rowindex], columns=['DAT_SCL'])[0][0][startfreq:endfreq+1]
    data3['DAT_SCL'][0][nchan:2*nchan]=fits[1].read(rows=[rowindex], columns=['DAT_SCL'])[0][0][startfreq+nchan_origin:endfreq+1+nchan_origin]
    #data reshape
    data3['DATA'][0][:,0,:,0]=tempspec1[:, startfreq:endfreq+1]

    data4['TSUBINT'][0]=fits[1].read(rows=[rowindex], columns=['TSUBINT'])[0][0]
    data4['OFFS_SUB'][0]=fits[1].read(rows=[rowindex], columns=['OFFS_SUB'])[0][0]
    data4['LST_SUB'][0]=fits[1].read(rows=[rowindex], columns=['LST_SUB'])[0][0]
    data4['RA_SUB'][0]=fits[1].read(rows=[rowindex], columns=['RA_SUB'])[0][0]
    data4['DEC_SUB'][0]=fits[1].read(rows=[rowindex], columns=['DEC_SUB'])[0][0]
    data4['GLON_SUB'][0]=fits[1].read(rows=[rowindex], columns=['GLON_SUB'])[0][0]
    data4['GLAT_SUB'][0]=fits[1].read(rows=[rowindex], columns=['GLAT_SUB'])[0][0]
    data4['FD_ANG'][0]=fits[1].read(rows=[rowindex], columns=['FD_ANG'])[0][0]
    data4['POS_ANG'][0]=fits[1].read(rows=[rowindex], columns=['POS_ANG'])[0][0]
    data4['PAR_ANG'][0]=fits[1].read(rows=[rowindex], columns=['PAR_ANG'])[0][0]
    data4['TEL_AZ'][0]=fits[1].read(rows=[rowindex], columns=['TEL_AZ'])[0][0]
    data4['TEL_ZEN'][0]=fits[1].read(rows=[rowindex], columns=['TEL_ZEN'])[0][0]
    data4['DAT_FREQ'][0]=fits[1].read(rows=[rowindex], columns=['DAT_FREQ'])[0][0][startfreq:endfreq+1]
    data4['DAT_WTS'][0]=fits[1].read(rows=[rowindex], columns=['DAT_WTS'])[0][0][startfreq:endfreq+1]
    
    data4['DAT_OFFS'][0][0:nchan]=fits[1].read(rows=[rowindex], columns=['DAT_OFFS'])[0][0][startfreq:endfreq+1]
    data4['DAT_OFFS'][0][nchan:2*nchan]=fits[1].read(rows=[rowindex], columns=['DAT_OFFS'])[0][0][startfreq+nchan_origin:endfreq+1+nchan_origin]
    data4['DAT_SCL'][0][0:nchan]=fits[1].read(rows=[rowindex], columns=['DAT_SCL'])[0][0][startfreq:endfreq+1]
    data4['DAT_SCL'][0][nchan:2*nchan]=fits[1].read(rows=[rowindex], columns=['DAT_SCL'])[0][0][startfreq+nchan_origin:endfreq+1+nchan_origin]
    #data reshape
    #    data4['DATA'][0][subindex,0,:,0]=tempspec2[startfreq:endfreq+1]
    data4['DATA'][0][:,0,:,0]=tempspec2[:, startfreq:endfreq+1]


    data5['TSUBINT'][0]=np.array(fits[1].read(rows=[rowindex], columns=['TSUBINT']))[0][0]
    data5['OFFS_SUB'][0]=np.array(fits[1].read(rows=[rowindex], columns=['OFFS_SUB']))[0][0]
    data5['LST_SUB'][0]=np.array(fits[1].read(rows=[rowindex], columns=['LST_SUB']))[0][0]
    data5['RA_SUB'][0]=np.array(fits[1].read(rows=[rowindex], columns=['RA_SUB']))[0][0]
    data5['DEC_SUB'][0]=np.array(fits[1].read(rows=[rowindex], columns=['DEC_SUB']))[0][0]
    data5['GLON_SUB'][0]=np.array(fits[1].read(rows=[rowindex], columns=['GLON_SUB']))[0][0]
    data5['GLAT_SUB'][0]=np.array(fits[1].read(rows=[rowindex], columns=['GLAT_SUB']))[0][0]
    data5['FD_ANG'][0]=np.array(fits[1].read(rows=[rowindex], columns=['FD_ANG']))[0][0]
    data5['POS_ANG'][0]=np.array(fits[1].read(rows=[rowindex], columns=['POS_ANG']))[0][0]
    data5['PAR_ANG'][0]=np.array(fits[1].read(rows=[rowindex], columns=['PAR_ANG']))[0][0]
    data5['TEL_AZ'][0]=np.array(fits[1].read(rows=[rowindex], columns=['TEL_AZ']))[0][0]
    data5['TEL_ZEN'][0]=np.array(fits[1].read(rows=[rowindex], columns=['TEL_ZEN']))[0][0]
    data5['DAT_FREQ'][0]=np.array(fits[1].read(rows=[rowindex], columns=['DAT_FREQ']))[0][0][startfreq:endfreq+1]
    data5['DAT_WTS'][0]=np.array(fits[1].read(rows=[rowindex], columns=['DAT_WTS']))[0][0][startfreq:endfreq+1]
    
    data5['DAT_OFFS'][0][0:nchan]=np.array(fits[1].read(rows=[rowindex], columns=['DAT_OFFS']))[0][0][startfreq:endfreq+1]
    data5['DAT_OFFS'][0][nchan:2*nchan]=np.array(fits[1].read(rows=[rowindex], columns=['DAT_OFFS']))[0][0][startfreq+nchan_origin:endfreq+1+nchan_origin]
    data5['DAT_SCL'][0][0:nchan]=np.array(fits[1].read(rows=[rowindex], columns=['DAT_SCL']))[0][0][startfreq:endfreq+1]
    data5['DAT_SCL'][0][nchan:2*nchan]=np.array(fits[1].read(rows=[rowindex], columns=['DAT_SCL']))[0][0][startfreq+nchan_origin:endfreq+1+nchan_origin]
    #data reshape
    #    data5['DATA'][0][subindex,0,:,0]=(tempspec1[startfreq:endfreq+1]+tempspec2[startfreq:endfreq+1])/2
    data5['DATA'][0][:,0,:,0] = (tempspec1[:, startfreq:endfreq+1] + tempspec2[:, startfreq:endfreq+1])/2


    if(index==0):
        outfreq=((startfreq+endfreq)*1.0/2+1.0)/((nchan-1.0)*1.0/2+1.0)*freq
        hdr=fitsio.FITSHDR()
        fitsout1.write(data3,header=hdr)
        fitsout1[0].write_key('HDRVER', h0['HDRVER'], comment="")
        fitsout1[0].write_key('FITSTYPE', h0['FITSTYPE'], comment="")
        fitsout1[0].write_key('DATE', h0['DATE'], comment="")
        fitsout1[0].write_key('OBSERVER', h0['OBSERVER'], comment="")
        fitsout1[0].write_key('PROJID', h0['PROJID'], comment="")
        fitsout1[0].write_key('TELESCOP', h0['TELESCOP'], comment="")
        fitsout1[0].write_key('ANT_X', h0['ANT_X'], comment="")
        fitsout1[0].write_key('ANT_Y', h0['ANT_Y'], comment="")
        fitsout1[0].write_key('ANT_Z', h0['ANT_Z'], comment="")
        fitsout1[0].write_key('FRONTEND', h0['FRONTEND'], comment="")
        fitsout1[0].write_key('NRCVR', 1, comment="")
        fitsout1[0].write_key('FD_POLN', h0['FD_POLN'], comment="")
        fitsout1[0].write_key('FD_HAND', h0['FD_HAND'], comment="")
        fitsout1[0].write_key('FD_SANG', h0['FD_SANG'], comment="")
        fitsout1[0].write_key('FD_XYPH', h0['FD_XYPH'], comment="")
        fitsout1[0].write_key('BACKEND', h0['BACKEND'], comment="")
        fitsout1[0].write_key('BECONFIG', h0['BECONFIG'], comment="")
        fitsout1[0].write_key('BE_PHASE', h0['BE_PHASE'], comment="")
        fitsout1[0].write_key('BE_DCC', h0['BE_DCC'], comment="")
        fitsout1[0].write_key('BE_DELAY', h0['BE_DELAY'], comment="")
        fitsout1[0].write_key('TCYCLE', h0['TCYCLE'], comment="")
        fitsout1[0].write_key('OBS_MODE', h0['OBS_MODE'], comment="")
        fitsout1[0].write_key('DATE-OBS', h0['DATE-OBS'], comment="")
        
        fitsout1[0].write_key('OBSFREQ', outfreq, comment="")
        fitsout1[0].write_key('OBSBW', chnum*1.0, comment="")
        fitsout1[0].write_key('OBSNCHAN', chnum, comment="")
        fitsout1[0].write_key('CHAN_DM', h0['CHAN_DM'], comment="")
        fitsout1[0].write_key('SRC_NAME', h0['SRC_NAME'], comment="")
        fitsout1[0].write_key('COORD_MD', h0['COORD_MD'], comment="")
        fitsout1[0].write_key('EQUINOX', h0['EQUINOX'], comment="")
        fitsout1[0].write_key('RA', h0['RA'], comment="")
        fitsout1[0].write_key('DEC', h0['DEC'], comment="")
        fitsout1[0].write_key('BMAJ', h0['BMAJ'], comment="")
        fitsout1[0].write_key('BMIN', h0['BMIN'], comment="")
        fitsout1[0].write_key('BPA', h0['BPA'], comment="")
        fitsout1[0].write_key('STT_CRD1', h0['STT_CRD1'], comment="")
        fitsout1[0].write_key('STT_CRD2', h0['STT_CRD2'], comment="")
        fitsout1[0].write_key('TRK_MODE', h0['TRK_MODE'], comment="")
        fitsout1[0].write_key('STP_CRD1', h0['STP_CRD1'], comment="")
        fitsout1[0].write_key('STP_CRD2', h0['STP_CRD2'], comment="")
        fitsout1[0].write_key('SCANLEN', h0['SCANLEN'], comment="")
        fitsout1[0].write_key('FD_MODE', h0['FD_MODE'], comment="")
        fitsout1[0].write_key('FA_REQ', h0['FA_REQ'], comment="")
        fitsout1[0].write_key('CAL_MODE', h0['CAL_MODE'], comment="")
        fitsout1[0].write_key('CAL_FREQ', h0['CAL_FREQ'], comment="")
        fitsout1[0].write_key('CAL_DCYC', h0['CAL_DCYC'], comment="")
        fitsout1[0].write_key('CAL_PHS', h0['CAL_PHS'], comment="")
        fitsout1[0].write_key('STT_IMJD', h0['STT_IMJD'], comment="")
        fitsout1[0].write_key('STT_SMJD', h0['STT_SMJD'], comment="")
        fitsout1[0].write_key('STT_OFFS', h0['STT_OFFS'], comment="")
        fitsout1[0].write_key('STT_LST', h0['STT_LST'], comment="")
        fitsout1[1].write_key('INT_TYPE','TIME', comment="")
        fitsout1[1].write_key('INT_UNIT','SEC', comment="")
        fitsout1[1].write_key('SCALE','FluxDec', comment="")
        fitsout1[1].write_key('NPOL', 1, comment="")
        fitsout1[1].write_key('POL_TYPE','AA', comment="")
        fitsout1[1].write_key('TBIN', tsample, comment="")
        fitsout1[1].write_key('NBIN',1, comment="")
        fitsout1[1].write_key('NBIN_PRD',0, comment="")
        fitsout1[1].write_key('PHS_OFFS',0.0, comment="")
        fitsout1[1].write_key('NBITS',8, comment="")
        fitsout1[1].write_key('NSUBOFFS',0, comment="")
        fitsout1[1].write_key('NCHAN',chnum, comment="")
        fitsout1[1].write_key('CHAN_BW',chan_bw, comment="")
        fitsout1[1].write_key('NCHNOFFS',0, comment="")
        fitsout1[1].write_key('NSBLK', nsblk, comment="")
        fitsout1[1].write_key('EXTNAME','SUBINT  ',comment="")

        fitsout2.write(data4,header=h0)
        fitsout2[0].write_key('HDRVER', h0['HDRVER'], comment="")
        fitsout2[0].write_key('FITSTYPE', h0['FITSTYPE'], comment="")
        fitsout2[0].write_key('DATE', h0['DATE'], comment="")
        fitsout2[0].write_key('OBSERVER', h0['OBSERVER'], comment="")
        fitsout2[0].write_key('PROJID', h0['PROJID'], comment="")
        fitsout2[0].write_key('TELESCOP', h0['TELESCOP'], comment="")
        fitsout2[0].write_key('ANT_X', h0['ANT_X'], comment="")
        fitsout2[0].write_key('ANT_Y', h0['ANT_Y'], comment="")
        fitsout2[0].write_key('ANT_Z', h0['ANT_Z'], comment="")
        fitsout2[0].write_key('FRONTEND', h0['FRONTEND'], comment="")
        fitsout2[0].write_key('NRCVR', 1, comment="")
        fitsout2[0].write_key('FD_POLN', h0['FD_POLN'], comment="")
        fitsout2[0].write_key('FD_HAND', h0['FD_HAND'], comment="")
        fitsout2[0].write_key('FD_SANG', h0['FD_SANG'], comment="")
        fitsout2[0].write_key('FD_XYPH', h0['FD_XYPH'], comment="")
        fitsout2[0].write_key('BACKEND', h0['BACKEND'], comment="")
        fitsout2[0].write_key('BECONFIG', h0['BECONFIG'], comment="")
        fitsout2[0].write_key('BE_PHASE', h0['BE_PHASE'], comment="")
        fitsout2[0].write_key('BE_DCC', h0['BE_DCC'], comment="")
        fitsout2[0].write_key('BE_DELAY', h0['BE_DELAY'], comment="")
        fitsout2[0].write_key('TCYCLE', h0['TCYCLE'], comment="")
        fitsout2[0].write_key('OBS_MODE', h0['OBS_MODE'], comment="")
        fitsout2[0].write_key('DATE-OBS', h0['DATE-OBS'], comment="")

        fitsout2[0].write_key('OBSFREQ', outfreq, comment="")
        fitsout2[0].write_key('OBSBW', chnum*1.0, comment="")
        fitsout2[0].write_key('OBSNCHAN', chnum, comment="")
        fitsout2[0].write_key('CHAN_DM', h0['CHAN_DM'], comment="")
        fitsout2[0].write_key('SRC_NAME', h0['SRC_NAME'], comment="")
        fitsout2[0].write_key('COORD_MD', h0['COORD_MD'], comment="")
        fitsout2[0].write_key('EQUINOX', h0['EQUINOX'], comment="")
        fitsout2[0].write_key('RA', h0['RA'], comment="")
        fitsout2[0].write_key('DEC', h0['DEC'], comment="")
        fitsout2[0].write_key('BMAJ', h0['BMAJ'], comment="")
        fitsout2[0].write_key('BMIN', h0['BMIN'], comment="")
        fitsout2[0].write_key('BPA', h0['BPA'], comment="")
        fitsout2[0].write_key('STT_CRD1', h0['STT_CRD1'], comment="")
        fitsout2[0].write_key('STT_CRD2', h0['STT_CRD2'], comment="")
        fitsout2[0].write_key('TRK_MODE', h0['TRK_MODE'], comment="")
        fitsout2[0].write_key('STP_CRD1', h0['STP_CRD1'], comment="")
        fitsout2[0].write_key('STP_CRD2', h0['STP_CRD2'], comment="")
        fitsout2[0].write_key('SCANLEN', h0['SCANLEN'], comment="")
        fitsout2[0].write_key('FD_MODE', h0['FD_MODE'], comment="")
        fitsout2[0].write_key('FA_REQ', h0['FA_REQ'], comment="")
        fitsout2[0].write_key('CAL_MODE', h0['CAL_MODE'], comment="")
        fitsout2[0].write_key('CAL_FREQ', h0['CAL_FREQ'], comment="")
        fitsout2[0].write_key('CAL_DCYC', h0['CAL_DCYC'], comment="")
        fitsout2[0].write_key('CAL_PHS', h0['CAL_PHS'], comment="")
        fitsout2[0].write_key('STT_IMJD', h0['STT_IMJD'], comment="")
        fitsout2[0].write_key('STT_SMJD', h0['STT_SMJD'], comment="")
        fitsout2[0].write_key('STT_OFFS', h0['STT_OFFS'], comment="")
        fitsout2[0].write_key('STT_LST', h0['STT_LST'], comment="")
        fitsout2[1].write_key('INT_TYPE','TIME', comment="")
        fitsout2[1].write_key('INT_UNIT','SEC', comment="")
        fitsout2[1].write_key('SCALE','FluxDec', comment="")
        fitsout2[1].write_key('NPOL', 1, comment="")
        fitsout2[1].write_key('POL_TYPE','AA', comment="")
        fitsout2[1].write_key('TBIN', tsample, comment="")
        fitsout2[1].write_key('NBIN',1, comment="")
        fitsout2[1].write_key('NBIN_PRD',0, comment="")
        fitsout2[1].write_key('PHS_OFFS',0.0, comment="")
        fitsout2[1].write_key('NBITS',8, comment="")
        fitsout2[1].write_key('NSUBOFFS',0, comment="")
        fitsout2[1].write_key('NCHAN',chnum, comment="")
        fitsout2[1].write_key('CHAN_BW',chan_bw, comment="")
        fitsout2[1].write_key('NCHNOFFS',0, comment="")
        fitsout2[1].write_key('NSBLK', nsblk, comment="")
        fitsout2[1].write_key('EXTNAME','SUBINT  ',comment="")



        fitsout3.write(data5,header=hdr)
        fitsout3[0].write_key('HDRVER', h0['HDRVER'], comment="")
        fitsout3[0].write_key('FITSTYPE', h0['FITSTYPE'], comment="")
        fitsout3[0].write_key('DATE', h0['DATE'], comment="")
        fitsout3[0].write_key('OBSERVER', h0['OBSERVER'], comment="")
        fitsout3[0].write_key('PROJID', h0['PROJID'], comment="")
        fitsout3[0].write_key('TELESCOP', h0['TELESCOP'], comment="")
        fitsout3[0].write_key('ANT_X', h0['ANT_X'], comment="")
        fitsout3[0].write_key('ANT_Y', h0['ANT_Y'], comment="")
        fitsout3[0].write_key('ANT_Z', h0['ANT_Z'], comment="")
        fitsout3[0].write_key('FRONTEND', h0['FRONTEND'], comment="")
        fitsout3[0].write_key('NRCVR', 1, comment="")
        fitsout3[0].write_key('FD_POLN', h0['FD_POLN'], comment="")
        fitsout3[0].write_key('FD_HAND', h0['FD_HAND'], comment="")
        fitsout3[0].write_key('FD_SANG', h0['FD_SANG'], comment="")
        fitsout3[0].write_key('FD_XYPH', h0['FD_XYPH'], comment="")
        fitsout3[0].write_key('BACKEND', h0['BACKEND'], comment="")
        fitsout3[0].write_key('BECONFIG', h0['BECONFIG'], comment="")
        fitsout3[0].write_key('BE_PHASE', h0['BE_PHASE'], comment="")
        fitsout3[0].write_key('BE_DCC', h0['BE_DCC'], comment="")
        fitsout3[0].write_key('BE_DELAY', h0['BE_DELAY'], comment="")
        fitsout3[0].write_key('TCYCLE', h0['TCYCLE'], comment="")
        fitsout3[0].write_key('OBS_MODE', h0['OBS_MODE'], comment="")
        fitsout3[0].write_key('DATE-OBS', h0['DATE-OBS'], comment="")

        fitsout3[0].write_key('OBSFREQ', outfreq, comment="")
        fitsout3[0].write_key('OBSBW', chnum*1.0, comment="")
        fitsout3[0].write_key('OBSNCHAN', chnum, comment="")
        fitsout3[0].write_key('CHAN_DM', h0['CHAN_DM'], comment="")
        fitsout3[0].write_key('SRC_NAME', h0['SRC_NAME'], comment="")
        fitsout3[0].write_key('COORD_MD', h0['COORD_MD'], comment="")
        fitsout3[0].write_key('EQUINOX', h0['EQUINOX'], comment="")
        fitsout3[0].write_key('RA', h0['RA'], comment="")
        fitsout3[0].write_key('DEC', h0['DEC'], comment="")
        fitsout3[0].write_key('BMAJ', h0['BMAJ'], comment="")
        fitsout3[0].write_key('BMIN', h0['BMIN'], comment="")
        fitsout3[0].write_key('BPA', h0['BPA'], comment="")
        fitsout3[0].write_key('STT_CRD1', h0['STT_CRD1'], comment="")
        fitsout3[0].write_key('STT_CRD2', h0['STT_CRD2'], comment="")
        fitsout3[0].write_key('TRK_MODE', h0['TRK_MODE'], comment="")
        fitsout3[0].write_key('STP_CRD1', h0['STP_CRD1'], comment="")
        fitsout3[0].write_key('STP_CRD2', h0['STP_CRD2'], comment="")
        fitsout3[0].write_key('SCANLEN', h0['SCANLEN'], comment="")
        fitsout3[0].write_key('FD_MODE', h0['FD_MODE'], comment="")
        fitsout3[0].write_key('FA_REQ', h0['FA_REQ'], comment="")
        fitsout3[0].write_key('CAL_MODE', h0['CAL_MODE'], comment="")
        fitsout3[0].write_key('CAL_FREQ', h0['CAL_FREQ'], comment="")
        fitsout3[0].write_key('CAL_DCYC', h0['CAL_DCYC'], comment="")
        fitsout3[0].write_key('CAL_PHS', h0['CAL_PHS'], comment="")
        fitsout3[0].write_key('STT_IMJD', h0['STT_IMJD'], comment="")
        fitsout3[0].write_key('STT_SMJD', h0['STT_SMJD'], comment="")
        fitsout3[0].write_key('STT_OFFS', h0['STT_OFFS'], comment="")
        fitsout3[0].write_key('STT_LST', h0['STT_LST'], comment="")
        fitsout3[1].write_key('INT_TYPE','TIME', comment="")
        fitsout3[1].write_key('INT_UNIT','SEC', comment="")
        fitsout3[1].write_key('SCALE','FluxDec', comment="")
        fitsout3[1].write_key('NPOL', 1, comment="")
        fitsout3[1].write_key('POL_TYPE','AA', comment="")
        fitsout3[1].write_key('TBIN', tsample, comment="")
        fitsout3[1].write_key('NBIN',1, comment="")
        fitsout3[1].write_key('NBIN_PRD',0, comment="")
        fitsout3[1].write_key('PHS_OFFS',0.0, comment="")
        fitsout3[1].write_key('NBITS',8, comment="")
        fitsout3[1].write_key('NSUBOFFS',0, comment="")
        fitsout3[1].write_key('NCHAN',chnum, comment="")
        fitsout3[1].write_key('CHAN_BW',chan_bw, comment="")
        fitsout3[1].write_key('NCHNOFFS',0, comment="")
        fitsout3[1].write_key('NSBLK', nsblk, comment="")
        fitsout3[1].write_key('EXTNAME','SUBINT  ',comment="")
    else:
        fitsout1[-1].append(data3)
        fitsout2[-1].append(data4)
        fitsout3[-1].append(data5)

fits.close()
fitsout1.close()
fitsout2.close()
fitsout3.close()
print "------------------Finish------------------"
