#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SSVEP project - funkcje

import numpy as np
import pylab as py
import scipy.signal as ss


def filtracja(kan,fs, czest):
	n=1
	[bx,ax]=ss.butter(n, czest/(fs/2.),btype='highpass')
	y1a = ss.filtfilt(bx,ax,kan)
	
	[b,a]=ss.butter(n, [48./(fs/2.), 52./(fs/2.)],btype='bandstop')
	y1 = ss.filtfilt(b,a,y1a)
 
	return y1

def idxStym(arr, odc=6e3,odl=500):
	'''zwraca indeksy rozpoczecia stymulacji
	arr - sygnal z FOTODIODY, odc - granica odcieca, odl - odleglosc miedzy stymulacjiami'''
	wk=np.where(arr>odc)[0]
	bc=np.hstack((wk[0],wk[:-1]))
	idx=np.where(wk-bc>odl)[0]+3 # +3 bo wczesniej jakis efekt brzegowy na diodzie
	return wk[idx]
	
def idxBreak(arr,odc=6e3,odl=500):
	'''zwraca indeksy rozpoczecia przerwy
	arr - sygnal z FOTODIODY, odc - granica odcieca, odl - odleglosc miedzy stymulacjiami'''
	wk=np.where(arr>odc)[0]
	bc=np.hstack((wk[1:],wk[-1]))
	idx=np.where(bc-wk>odl)[0]
	return wk[idx]

def tnij(data,idxwyc,length):
	'''wycina z macierzy data fragmenty zaczynajace sie od momentu idxwyc 
	i trawajace length w probkach'''
	ls=[]
	for i in idxwyc:
		ls.append(data[i:i+length])
	return np.array(ls)

def periodogram(s,w,fs,k=1):
	N=len(s)
	w*=1./np.sqrt(np.dot(w,w))
	F=np.fft.fft(s*w,int(k*N))
	ft=np.fft.fftshift(abs(F))
	fv=np.fft.fftshift(np.fft.fftfreq(len(ft),1./fs))
	return (ft**2,fv)


def scalogram(data):
	bottom = 0

	vmin = min(map(lambda x: min(abs(x)), data))
	vmax = max(map(lambda x: max(abs(x)), data))

	py.gca().set_autoscale_on(False)

	for row in range(0, len(data)):
		scale = 2.0 ** (row - len(data))

		py.imshow(
			np.array([abs(data[row])]),
			interpolation = 'nearest',
			vmin = vmin,
			vmax = vmax,
			extent = [0, 1, bottom, bottom + scale])

		bottom += scale

def przygotuj_sygn(danne, kanal, Fs=512, chanNr=24, ear1=20):
	'''dannne - caly plik z danymi, kanal - nr kanalu ktory chcemy wyciagnac
	Fs - wiadomo, chanNr- calkowita liczba kanalow, ear1 - nr montazu 
	pierwszej elektrody z ucha'''
	kan=danne[kanal-1::chanNr]
	ucho1=danne[ear1-1::chanNr]
	ucho2=danne[ear1::chanNr]
	
	kan=kan-(ucho1+ucho2)/2
	kan=filtracja(kan,Fs,0.5)
	kan-=np.mean(kan)
	return kan

def ind_cz(fq,step=15):
	"""Zwraca indeksy cz. podstawowej (step) i jej harmonicznych"""
	indlis=[]
	for j in range(0,20,step)[1:]:
		idx=np.where(fq>=(j-0.1))[0][0]
		indlis.append(idx)
	return indlis

def rob_liste_mocy(arr,fs=512):
	Flg=1
	win=np.ones(len(arr[0]))
	wek=np.zeros(len(arr))
	for e,tr in enumerate(arr):
		(P,fv)=periodogram(tr,win,fs)
		if Flg==1:
			ind=ind_cz(fv)
			Flg=0
		wek[e]=P[ind]
	return wek

def rob_liste_fazy(arr,fs=512):
	"""zwraca wartosci fazy dla odpowiedniego trialu"""
	Flg=1
	phase=np.zeros(len(arr))
	r=np.zeros(len(arr))
	for e,tr in enumerate(arr):
		F=np.fft.fftshift(np.fft.fft(tr))
		fv=np.fft.fftshift(np.fft.fftfreq(len(F),1./fs))
		if Flg==1:
			ind=ind_cz(fv)
			Flg=0
		phase[e]=np.angle(F[ind])
		r[e]=abs(F[ind])
	return phase,r

def oblicz_roznice_mocy(arr_bez,arr_sty,wiecej=0):
	rpc=np.zeros(len(arr_bez))
	wek_r,wek_s=rob_liste_mocy(arr_bez),rob_liste_mocy(arr_sty)
	if np.any(wek_r>2e7):
		wek_r[np.where(wek_r>2e7)[0]]=0
		wek_s[np.where(wek_r>2e7)[0]]=0
	
	rpc=wek_s-wek_r


	if wiecej==0: return rpc
	else: return rpc,wek_s, wek_r

def ryspary(osx,a1,osx2,a2,pokaz=True,lim_y=None):
	py.subplot(211)
	if osx != None: py.plot(osx,a1)
	else: py.plot(a1)
	if lim_y != None: py.ylim([0,lim_y])
	py.subplot(212)
	if osx2 != None: py.plot(osx2,a2)
	else: py.plot(a2)
	if lim_y != None: py.ylim([0,lim_y])
	
	if pokaz: py.show()

def ryspary_brzeg(osx,a1,osx2,a2,pokaz=True,lim_y=None):
	py.subplot(211)
	if osx != None:
		py.plot(osx[len(osx)/2:],a1[len(a1)/2:])
		py.xlim([0,100])
		py.xticks(np.arange(0,100,5))
	else:
		py.plot(a1[:len(a1)/2])
	if lim_y != None: py.ylim([0,lim_y])
	py.subplot(212)
	if osx2 != None:
		py.plot(osx2[len(osx2)/2:],a2[len(a2)/2:])
		py.xlim([0,100])
		py.xticks(np.arange(0,100,5))
	else:
		py.plot(a2[:len(a2)/2])
	if lim_y != None: py.ylim([0,lim_y])
	
	if pokaz: py.show()
