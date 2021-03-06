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

def ind_cz(fq,P=None, step=15, rfval=None):
	"""troche sie funkcja skomplikowala ale musialem jakos obejsc
	problem roznych indeksow wektora czestosci z ref i sty i faktu
	ze pik stymulacji nie zawsze jest w 15 Hz...
	
	Zwraca indeksy cz. podstawowej (step) i jej harmonicznych"""
	if P==None:
		indlis=[]
		for j in range(0,20,step)[1:]:
			idx=np.where(fq>(j-0.01))[0][0]
			indlis.append(idx)
		return indlis
	else:
		if rfval==None:
			czp,czk=step-0.2,step+0.2
			ind_czp, ind_czk = np.where(fq>=czp)[0][0],np.where(fq>=czk)[0][0]
			ind_of_sty=np.argmax(P[ind_czp:ind_czk])
			#print fq[ind_czp+ind_of_sty],
			return ind_czp+ind_of_sty
		else:
			idx=np.where(fq>=rfval)[0][0]
			return idx

def rob_liste_mocy(arr,fs=512,czest=15):
	'''dla podanego wektora z trialami wyrzuca wektor
	wartosci mocy obliczonej periodogramem dla czestotliwosci'''
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

def lista_mocy_stym(arr,fs=512,czest=15):
	win=np.ones(len(arr[0])) #okno
	wek=np.zeros(len(arr))
	ind_of_mx=wek.copy()
	for e,trial in enumerate(arr):
		(P,fv)=periodogram(trial,win,fs)
		idx=ind_cz(fv,P,step=czest)
		ind_of_mx[e]=fv[idx]
		#print fv[idx],
		wek[e]=P[idx]
	return wek, ind_of_mx

def lista_mocy_ref(arr,ind_mx,fs=512,czest=15):
	win=np.ones(len(arr[0])) #okno
	wek=np.zeros(len(arr))
	for e,trial in enumerate(arr):
		(P,fv)=periodogram(trial,win,fs)
		idx=ind_cz(fv,P,step=czest,rfval=ind_mx[e])
		wek[e]=P[idx]
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
	'''po zadaniu arr_bez - wektora z trialami bez stymulacji i
	arr_sty wektora z symulacja wyrzuca wartosci roznicy dla czestosci
	ustawionej w funkcji rob_liste_mocy 
	!!!!!
	lub jesli uzyto lista_mocy_... to 
	dla wiecej=True podaje tez wektory z wartoscia mocy stym i ref'''
	rpc=np.zeros(len(arr_bez))
	wek_s, ind_mx = lista_mocy_stym(arr_sty)
	wek_r = lista_mocy_ref(arr_bez,ind_mx)
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
