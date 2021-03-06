#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SSVEP project

import numpy as np
import pylab as py
import scipy.signal as ss
from ssvepfun import *
import os
import pywt # do dekompozycji falkowej
import sys

FOLDER='data63/'
try:
	NAZWA= FOLDER+sys.argv[1]
except Exception:
	NAZWA=FOLDER+'h6030_s06_mg_o1'

o1_bez=np.load(NAZWA+'_bez.npy')
o1_sty=np.load(NAZWA+'_stym.npy')

fs=512


def usrednione_widmo(arr):
	win=np.ones(len(arr[0]))
	l=0
	Flag=1
	for trial in arr:
		(P,fv)=periodogram(trial,win,fs)
		l+=1
		if Flag==1:
			wek=np.zeros(len(P))
			Flag=0
		wek+=P
	return wek/l,fv

def usrednione_sygnaly():
	sred_r,sred_s=np.average(o1_bez,axis=0),np.average(o1_sty,axis=0)	
	Ps,fv=periodogram(sred_s,np.ones(len(sred_s)),fs)
	Pr,fv=periodogram(sred_r,np.ones(len(sred_r)),fs)
	ryspary(fv,Pr,fv,Ps)


def histogramy():
	py.subplot(211)
	py.hist(wr,bins=40)
	py.xlim([0,2.5e6])
	py.subplot(212)
	py.hist(ws,bins=40)
	py.xlim([0,2.5e6])
	py.show()

def spektrogramy():
	sred_r,sred_s=np.average(o1_bez,axis=0),np.average(o1_sty,axis=0)
	py.subplot(121)
	py.specgram(sred_r, Fs=fs, scale_by_freq=True)
	py.title('ref')
	py.subplot(122)
	py.specgram(sred_s, Fs=fs, scale_by_freq=True)
	py.title('stym')
	py.show()

def wspkor(ar1,ar2):
	'''oblicza wspolczynnik korelacji dla dwu prob takiej samej wielkosci'''
	return np.corrcoef(ar1,ar2,ddof=1)[0][1]

rm,ws,wr=oblicz_roznice_mocy(o1_bez, o1_sty,wiecej=8)

py.plot(np.arange(len(rm)),rm,'go')
py.xlabel('nr trialu')
py.ylabel('roznica mocy')
py.savefig(os.getcwd()+'/result/'+NAZWA[-15:]+'moctrial.png')
#py.show()

#uw_ref,fv1=usrednione_widmo(o1_bez)
#uw_sty,fv2=usrednione_widmo(o1_sty)

#ryspary(fv1,uw_ref,fv2,uw_sty,lim_y=np.max(uw_sty)+1e4)

#print 'wspolczynnik korelacji', wspkor(rm,np.arange(len(rm)))
#spektrogramy()

def fazy():
	pha,amp=rob_liste_fazy(o1_sty)
	py.figure()
	py.polar(pha,amp,'o')
	py.savefig(os.getcwd()+'/result/'+NAZWA[-15:]+'faz.png')
#fazy()

