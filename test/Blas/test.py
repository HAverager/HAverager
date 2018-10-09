#!/usr/bin/env python

from numpy import *


def gen(ndata=1000,nsys=1000,nset=2):
    ''' Generate data ''' 

    stat_s = random.normal(0.,1.,(ndata,nset))
    syst_b = random.normal(0.,1.,nsys)
    syst_g = random.normal(0.,1.,(nsys,ndata,nset))

    cent = 100.+stat_s+sum(syst_b[:,newaxis,newaxis]*syst_g,axis=0)
    stat = ones( (ndata,nset) )
    return cent,stat,syst_g
    
# Set path of the averager
import sys
sys.path.append('../bin')
import averager

generate = 0

nm = 2000
ns = 2000
if generate>0:
    ce,er,sy = gen(nm,ns,2)


    ce.tofile("ce.dat")
    er.tofile("er.dat")
    sy.tofile("sy.dat")

    print sy.shape
    exit(0)

ce = fromfile("ce.dat").reshape(nm,2)
er = fromfile("er.dat").reshape(nm,2)
sy = fromfile("sy.dat").reshape(ns,nm,2)

#initialization (optional information)
averager.avin.initvariables()
averager.avin.setoutputfolder('./o')
averager.avin.initeration = 3
#perform averaging
dataAv,statAv,systAv = averager.average(ce,er,sy)

print dataAv
