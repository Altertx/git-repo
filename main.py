#!/usr/bin/python
# -*- coding: UTF-8 -*-

import numpy as np
from geo_v1 import *

import sys
sys.path.append('C:\\Users\\BRONX\\SEM_4\\INF 2\\Projekt_geo')

geo = Transformacje(model = "wgs84")

plik = "wsp_inp.txt"
# odczyt z pliku: https://docs.scipy.org/doc/numpy-1.15.1/reference/generated/numpy.genfromtxt.html
tablica = np.genfromtxt(plik, delimiter=',', skip_header = 4)

PL = []
for krotka in tablica:
    PL.append(geo.xyz2plh(krotka[0], krotka[1], krotka[2]))
PL_arr = np.array(PL)

XY_GK_2000 = []
for krotka in PL:
    XY_GK_2000.append(geo.pl2xy(krotka[0], krotka[1], geo.L0_2000))
# nie tworze arraya XY_GK_2000 bo nie bede wyswietlal tych wartosci

XY_2000 = []
for krotka in XY_GK_2000:
    XY_2000.append(geo.u2000(krotka[0], krotka[1]))
XY_2000_arr = np.array(XY_2000)

XY_GK_1992 = []
for krotka in PL:
    XY_GK_1992.append(geo.pl2xy(krotka[0], krotka[1], geo.L0_1992))
# nie tworze arraya XY_GK_1992 bo nie bede wyswietlal tych wartosci

XY_1992 = []
for krotka in XY_GK_1992:
    XY_1992.append(geo.u1992(krotka[0], krotka[1]))
XY_1992_arr = np.array(XY_1992)



wszystko = np.column_stack([PL_arr,XY_2000_arr,XY_1992_arr])
# zapis: https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.savetxt.html
np.savetxt("wsp_out.txt", wszystko, delimiter=',', fmt = ['%11.8f', '%12.8f', '%8.3f' , '%12.3f' , '%12.3f' , '%11.3f' , '%11.3f'], header = 'konwersja współrzednych geodezyjnych \n\\\\ Albert Kalinowski')

    
# np.savetxt("wsp_out.txt" , wyniki , delimiter = ',',fmt = ['%10.2f', '%10.2f', '%10.3f'])
