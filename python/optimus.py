#!/usr/bin/python
# -*- coding: UTF-8 -*-
import numpy as np
import math
import os, sys, glob, re
from commands import getoutput
import time
import datetime

#return string list
class telescope():
    "aaa"
    def __init__(self, name = 'FAST system', f1 = 0.125, f2 = 1024.125, nchan = 4096, t0 = 0, t1 = 52.4288, gain = 1, tsys = 30, raj = 0, decj = 0, useAngel = 0, tsamp = 0.0002, digitiser = 8):
        self.name = name
        self.f1 = f1
        self.f2 = f2 
        self.nchan = nchan
        self.t0 = t0
        self.t1 = t1
        self.gain = gain 
        self.tsys = tsys
        self.raj =  raj
        self.decj =  decj
        self.useAngel = useAngel
        self.tsamp = tsamp
        self.digitiser = digitiser

    def __str__(self):
        params=['name','f1', 'f2', 'nchan', 't0','t1', 'gain', 'tsys', 'raj', 'decj', 'useAngel', 'tsamp', 'digitiser']
        return 'name: %s\nf1: %.3f \nf2: %.3f\nnchan: %.0f\nt0: %.5f\nt1: %.5f\ngain: %.2f\ntsys: %.2f\nraj: %.9f\ndecj: %.9f\nuseAngle: %.0f\ntsamp: %.6f\ndigitiser: %.0f\n' %(self.name, self.f1, self.f2, self.nchan, self.t0, self.t1, self.gain, self.tsys, self.raj, self.decj, self.useAngel, self.tsamp, self.digitiser)
    __repr__ = __str__
    
    def __doc__(self):
        'usage: import optimus'



#set pulsar params
#return string list
class pulsar():
    #'name: J1950+30\n', 'p0: 0.2\n', 'dm: 50.0\n', 'raj: 4.510914803\n', 'decj: 0.13602659 \n', 'width: 0.01\n', 'flux: 0.000008\n', 'useAngle: 0\n'
    def __init__(self, name = 'J1950+30', period = '0.2', dm = '50', raj = 4.510914803, decj = 0.13602659, width = 0.01, flux = '0.001', useAngel = 0):
        self.name = name
        self.dm = dm
        self.period = period
        self.raj = raj
        self.decj = decj
        self.width = width
        self.flux = flux
        self.useAngel = useAngel
    def __str__(self):
        #PSRparams=['name: J1950+30\n', 'p0: 0.2\n', 'dm: 50.0\n', 'raj: 4.510914803\n', 'decj: 0.13602659 \n', 'width: 0.01\n', 'flux: 0.000008\n', 'useAngle: 0\n']
        #return 'name: %s\np0: %.9f \ndm: %.2f\nraj: %.9f\ndecj: %.9f\nflux: %.9f\nwidth: %.4f\nuseAngle: %.0f' %(self.name, self.period, self.dm, self.raj, self.decj, self.flux, self.width,self.useAngel)
        return 'name: %s\np0: %s \ndm: %s\nraj: %.9f\ndecj: %.9f\nflux: %s\nwidth: %.4f\nuseAngle: %.0f\n' %(self.name, self.period, self.dm, self.raj, self.decj, self.flux, self.width,self.useAngel)
    __repr__ = __str__


#class defaultparams():
#    def __init__(self, name = ""):
#        self.name = name
#        self.data_dic = {}
#        self.index = -1
    


#write paramfiles
class writeParamFile() :
    def __init__(self,filename,paramStruct):
        # write params
        f=open(filename,'w')
        f.writelines(str(paramStruct))
        f.close()

#calculate the width of the pulsar.
#input period (s) meanFref (MHz) randomNum (-1~1)
def pulseWidth(period,meanFref,randomNum):
    meanFref = meanFref/1000.
    B=1.*1E12
    R_LightCylinder = 4.77*10000.*period
    periodDot = ((B/(3.2*1E19))**2)/period
    R_EmisionBeam = 400.*pow(meanFref,-0.26)*pow((periodDot*(1E15)),0.07)*pow(period,0.3)
    rho = 86.*pow(R_EmisionBeam/R_LightCylinder,0.5)
    rho = rho/360.
    w = 2*rho*pow((1.-randomNum**2),0.5)
    return w

#if __name__ == "__main__":
#    print '''
#usage: import optimus
#
#    '''
#print optimus.telescope()
#print optimus.pulsar()
#print optimus.writeParamFile(outputFilename, optimus.telescope())
#print optimus.writeParamFile(outputFilename, optimus.pulsar(dm = dm, flux = fu, period = p, width = width))
#    '
#set telescope params
