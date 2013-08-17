#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SSVEP project

import numpy as np
import pylab as py
import scipy.signal as ss
from ssvepfun import *

#np.set_printoptions(threshold='nan')

def filtracja(kan,fs, czest):
	n=1
	[b,a]=ss.butter(n, [48./(fs/2.), 52./(fs/2.)],btype='bandstop')
	y1 = ss.filtfilt(b,a,kan)
	#[b,a]=ss.butter(n, [148./(fs/2.), 152./(fs/2.)],btype='bandstop')
	#y2 = ss.filtfilt(b,a,y1)

	[bx,ax]=ss.butter(n, czest/(fs/2.),btype='highpass')
	y1a = ss.filtfilt(bx,ax,y1) 
	return y1a

dane=np.fromfile('data63/h6030_s01_mn.raw','float32')
chanNr= 24#liczba kanalow
fs=512.

ch_o1=18
ear1=20
ch_diod=22

o1=dane[ch_o1-1::chanNr]
diod=dane[ch_diod-1::chanNr]

b1,b2=3550, 4e3
first20=filtracja(o1[b1*fs:fs*b2],fs,1)
d20=diod[b1*fs:fs*b2]
wk=np.where(d20>8e3)[0]
#wk2=np.where(diod>6e3)[0]
print 'dlugosc',len(diod)*1./(fs*60)
print 'mx',np.max(d20)
print 'len f20', len(first20), np.max(first20)

def rysownik():
	py.subplot(211)
	py.title('o1')
	py.plot(b1+np.arange(len(first20))/fs,first20)
	py.ylim([-2e3,2e3])
	py.subplot(212)
	py.title('dioda')
	py.plot(b1+np.arange(len(d20))/fs,d20, '.')
	py.ylim([0,200e3])
	#py.figure()
	#py.plot(d20[wk],'ro')
	py.show()

rysownik()

def idxStym(arr,odl=500):
	'''zwraca indeksy rozpoczecia stymulacji
	arr - '''
	bc=np.hstack((arr[0],arr[:-1]))
	idx=np.where(arr-bc>odl)[0]
	return arr[idx]

def idxBreak(arr,odl=500):
	bc=np.hstack((arr[1:],arr[-1]))
	idx=np.where(bc-arr>odl)[0]
	return arr[idx]


def rys2():
	tryin1=idxStym(wk)
	tryin2=idxBreak(wk)
	py.plot(np.arange(len(d20))/fs,d20)
	py.plot(tryin1/fs,np.ones(len(tryin1))*10e3,'ro')
	py.plot(tryin2/fs,np.ones(len(tryin2))*10e3,'gx')
	py.show()

rys2()
def tnij(data,idxwyc,length=12*fs):
	ls=[]
	for e,i in enumerate(idxwyc):
		ls.append(data[i:i+length])
	return np.array(ls)

#ad=idxBreak(wk2)
#add=idxStym(wk2)
#arrTrial_diod=tnij(diod,ad)

#py.plot(arrTrial_diod[14])
#py.show()

def periodogram2(sygnal,Fs,okno):
 s=sygnal*okno
 N_fft=len(sygnal)
 S=np.fft.fft(s,N_fft)
 P=S*S.conj()
 P=P.real/np.sum(okno**2)
 F=np.fft.fftfreq(N_fft,1./Fs)
 return (np.fft.fftshift(P),np.fft.fftshift(F))

#NAZWA='data/h5_s01_mn_o1'
#o1_bez=np.load(NAZWA+'_stym.npy')

#F1,fq1=periodogram(o1_bez[0],fs,np.hamming(len(o1_bez[0])))
#F2,fq2=periodogram2(o1_bez[0],fs,np.hanning(len(o1_bez[0])))
#ryspary(fq1,F1,fq2,F2)
