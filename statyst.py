#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SSVEP project - statystyka

import numpy as np
import pylab as py
import scipy.signal as ss
import scipy.stats as st
from ssvepfun import *
import sys

try:
	NAZWA='data/'+sys.argv[1]
except Exception:
	NAZWA='data/h5_s05_ml_o1'


o1_bez=np.load(NAZWA+'_bez.npy')
o1_sty=np.load(NAZWA+'_stym.npy')

fs=512

def randsample(x):
	'''zwraca wektor z losowo wybranymi elementami wektora x. Losowanie odbywa się z powtórzeniami''' 
	n=len(x)
	ind = np.random.randint(n, size = n)
	y = x[ind]
	return y

rm,ws,wr=oblicz_roznice_mocy(o1_bez, o1_sty,wiecej=8)
print 'wczytano'

N_rep=1e4
r_rf = np.zeros((N_rep,1))
for i in np.arange(N_rep):
	x_rf = randsample(wr)
	r_rf[i] = np.mean(x_rf)

print 'zakonczono losowanie'

ci_rf = st.scoreatpercentile(r_rf, 95)
print 'Z bootstrapu przedzial gorny od: ', ci_rf

def zhistem():
	py.hist(ws)
	py.hist(wr)
	py.plot(ci_rf,100,'ro')
	py.show()

zhistem()

print 'wieksze\nz ref', np.sum(wr>ci_rf)
print 'z sty %d -procetowo %d' %(np.sum(ws>ci_rf),np.sum(ws>ci_rf)*1e2/len(ws))
print 'mniejsze\nz ref', np.sum(wr<ci_rf)
print 'z sty %d -procetowo %d' %(np.sum(ws<ci_rf),np.sum(ws<ci_rf)*1e2/len(ws))

print np.sum(wr>ci_rf)*1./len(wr)
