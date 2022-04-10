# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 13:10:35 2022

@author: Ola
"""

import numpy as np
from proj import *

plik = "wsp_inp.txt"
tablica = np.genfromtxt(plik, delimiter=',', skip_header = 4)

T = Transformacje(model = "grs80")

blh = []
u92 = []
for i in range(len(tablica)):
    b, l, h = T.xyz2blh_hirvonen(tablica[i][0], tablica[i][1], tablica[i][2])
    x92, y92 = T.u1992(b, l)
    
    blh.append([b, l, h])
    u92.append([x92, y92])
    

wsp_out = np.hstack([blh, u92])    
np.savetxt("wsp_out.txt", wsp_out, delimiter=',  ', fmt = ['%10.2f', '%10.2f', '%10.3f', '%10.3f', '%10.3f'], header = 'Konwersja współrzędnych geodezyjnych \\ Aleksandra Skolimowska \n blh      układ 1992')
  