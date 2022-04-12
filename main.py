import numpy as np
from geo_v1 import *

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


dX = []
i = 1
while i < len(tablica):
    dX.append(tablica[i]-tablica[0])
    i+=1
dX_arr = np.array(dX)
dX_arr = np.transpose(dX_arr)

dx = []
a = np.hsplit(dX_arr,11)
for i in a:
    dx.append(geo.xyz2neu(PL_arr[0,0], PL_arr[0,1], i))
    
dx_arr = np.hstack((dx[0],dx[1],dx[2],dx[3],dx[4],dx[5],dx[6],dx[7],dx[8],dx[9],dx[10]))
dx_arr = np.transpose(dx_arr) # NEU

s_Az_z = []
for krotka in dx_arr:
    s_Az_z.append(geo.s_azimuth_elevation(krotka[0], krotka[1], krotka[2]))
s_Az_z_arr = np.array(s_Az_z)


zeros = np.array([[0,0,0]])
dx_arr = np.concatenate((zeros,dx_arr), axis = 0)
s_Az_z_arr = np.concatenate((zeros,s_Az_z_arr), axis = 0)


results = np.column_stack([PL_arr,XY_2000_arr,XY_1992_arr,dx_arr,s_Az_z_arr])
# zapis: https://docs.scipy.org/doc/numpy-1.15.0/reference/generated/numpy.savetxt.html
np.savetxt("wsp_out.txt", results , delimiter=',', fmt = ['%11.8f', '%12.8f', '%8.3f' , '%12.3f' , '%12.3f' , '%11.3f' , '%11.3f' , '%7.3f' , '%7.3f' , '%6.3f' , '%6.3f' , '%13.8f' , '%12.8f'],
header = 'konwersja współrzednych geodezyjnych \nAlbert Kalinowski\n\nfi[dec_d]   la[dec_d]     h[m]'
'     X2000[m]    Y2000[m]     X1992[m]   Y1992[m]      N[m]    E[m]   U[m]   s[m]    Az[dec_d]    z[dec_d]')


