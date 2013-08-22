#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SSVEP project

import numpy as np
import pylab as py
import scipy.signal as ss
from ssvepfun import *
import os

#### wlasciwosci
NAZWA_PLIKU='h6030_s05_dk'
FOLDER='data63'
ELEKTRODA='o1' # <- nazwa elektrody
ch_o=18# 18 - O1, 19 - O2, 15 - Pz  <---nr elektrody z pliku XML
SKALA=3e6
OBRAZKI=0
SAVE=0 # zapis pocietych danych do pliku

####

dane=np.fromfile(FOLDER+'/'+NAZWA_PLIKU+'.raw','float32')
chanNr= 24#liczba kanalow
fs=512.
ear1=20
ch_diod=22
diod=dane[ch_diod-1::chanNr]

def przygotuj_sygnal(danne, kanal, Fs=512, chanNr=24, ear1=20):
	'''dannne - caly plik z danymi, kanal - nr kanalu ktory chcemy wyciagnac
	Fs - wiadomo, chanNr- calkowita liczba kanalow, ear1 - nr montazu 
	pierwszej elektrody z ucha'''
	kan=danne[kanal-1::chanNr]
	ucho1=danne[ear1-1::chanNr]
	ucho2=danne[ear1::chanNr]
	
	kan=kan-(ucho1+ucho2)/2
	kan=filtracja(kan,Fs,3)
	kan-=np.mean(kan)
	return kan

# w tym momencie dane z elektrody odfiltrowane, odjete uszy, i odjeta srednia
elek=przygotuj_sygnal(dane,ch_o)
przerwy=idxBreak(diod)
stym=idxStym(diod)

elek_bez=tnij(elek,przerwy+0.15*fs,29.85*fs)
elek_sty=tnij(elek,stym,60.15*fs)

def zapisz():
	np.save(FOLDER+'/'+NAZWA_PLIKU+'_'+ELEKTRODA+'_bez',elek_bez)
	np.save(FOLDER+'/'+NAZWA_PLIKU+'_'+ELEKTRODA+'_stym',elek_sty)
	print 'zapisano'

if SAVE==True: zapisz()


def zobacz_pojedynczy_trial( PZ=240 ):
	'PZ - numer realizacji, okno prostokatne'

	win=np.ones(elek_bez[PZ].size)
	wins=np.ones(elek_sty[PZ].size)
	(P,fv)=periodogram(elek_bez[PZ],win,fs)
	(Ps,fvs)=periodogram(elek_sty[PZ],wins,fs)

	ryspary_brzeg(fv,P,fvs,Ps)

def zapis_obr():
	''''zapisuje obrazki z periodogramu oknem prostakatnym sygnalu
	z elek_sty i elek_bez w wpliku zadanym w NEWPATH'''
	win=np.ones(elek_bez[0].size)
	wins=np.ones(elek_sty[0].size)
	for PZ in xrange(len(elek_bez)-1):
		(P,fv)=periodogram(elek_bez[PZ],win,fs)
		(Ps,fvs)=periodogram(elek_sty[PZ],wins,fs)
		py.figure()
		ryspary_brzeg(fv,P,fvs,Ps,pokaz=False,lim_y=SKALA)
		
		newpath = os.getcwd()+FOLDER+'/obr/'+NAZWA_PLIKU+ELEKTRODA
		if not os.path.exists(newpath): os.makedirs(newpath)
		py.savefig(os.getcwd()+FOLDER+'/obr/'+NAZWA_PLIKU+ELEKTRODA+'/'+str(PZ)+'.png')

		print PZ,

if OBRAZKI==True: zapis_obr()
