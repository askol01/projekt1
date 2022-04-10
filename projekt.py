# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 13:10:35 2022

@author: Ola
"""

import numpy as np
from proj import *

plik = "wsp_inp.txt"
tablica = np.genfromtxt(plik, delimiter=',', skip_header = 4)

T_grs80 = Transformacje(model = "grs80")


while True:
    print('\nWybierz opcję')
    print('\nA - Transformacja współrzędnych ortokartezjańskich X, Y, Z na współrzędne geodezyjne fi, lambda, h')
    print('B - Transformacja współrzędnych geodezyjnych fi, lambda, h na współrzędne ortokartezjańskie X, Y, Z')
    print('C - Transformacja współrzednych geodezyjnych do układu 1992')
    print('D - Transformacja współrzednych geodezyjnych do układu 2000')
    print('E - Transformacja współrzędnych ortokartezjańskich X, Y, Z na współrzędne topocentryczne N, E, U')
    print('F - Obliczenie odległosci 2D i 3D')
    print('G - Obliczenie kąta azymutu i kąta elewacji')
    print('H - Zakończ')
    
    opcja = input('Wybierz opcję = ')
    
    if opcja == 'H':
        break
    
    elif opcja == 'A':
        print('\nTransformacja współrzędnych ortokartezjańskich X, Y, Z na współrzędne geodezyjne fi, lambda, h')
        print(' fi  |  lambda   |   h')
        blh = []
        for i in range(len(tablica)):
            b, l, h = T_grs80.xyz2blh_hirvonen(tablica[i][0], tablica[i][1], tablica[i][2])
            blh.append([b, l, h])
        blh = np.array(blh)
        print(blh)
        
    elif opcja == 'B':
        print('\nTransformacja współrzędnych geodezyjnych fi, lambda, h na współrzędne ortokartezjańskie X, Y, Z')
        print(' X  |  Y   |   Z')
        XYZ = []
        for i in range(len(tablica)):
            b, l, h = T_grs80.xyz2blh_hirvonen(tablica[i][0], tablica[i][1], tablica[i][2])
            X, Y, Z = T_grs80.blh2xyz_odwrotny(b, l, h)
            XYZ.append([X, Y, Z])
        XYZ = np.array(XYZ)
        print(XYZ)
        
    elif opcja == 'C':
        print('\nTransformacja współrzednych geodezyjnych do układu 1992')
        print(' x92  |  y92  ')
        u92 = []
        for i in range(len(tablica)):
            b, l, h = T_grs80.xyz2blh_hirvonen(tablica[i][0], tablica[i][1], tablica[i][2])
            x92, y92 = T_grs80.u1992(b, l)
            u92.append([x92, y92])
        u92 = np.array(u92)
        print(u92)
        
    elif opcja == 'D':
        print('\nTransformacja współrzednych geodezyjnych do układu 2000')
        print(' x00  |  y00  ')
        u00 = []
        for i in range(len(tablica)):
            b, l, h = T_grs80.xyz2blh_hirvonen(tablica[i][0], tablica[i][1], tablica[i][2])
            x00, y00 = T_grs80.u2000(b, l)
            u00.append([x00, y00])
        u00 = np.array(u00)
        print(u00)
    
    elif opcja == 'E':
        print('\nTransformacja współrzędnych ortokartezjańskich X, Y, Z na współrzędne topocentryczne N, E, U')
        print(' N  |  E  |  U')
        NEU = []
        X0 = 0
        Y0 = 0
        Z0 = 0
        for i in range(len(tablica)):
            n,e,u = T_grs80.NEU(X0, Y0, Z0, tablica[i][0], tablica[i][1],tablica[i][2])
            NEU.append([n,e,u])
        NEU = np.array(NEU)
        print(NEU)
    
    elif opcja == 'F':
        print('\nObliczenie odległosci 2D i 3D')
        print(' 2D  |  3D  ')
        X0 = 0
        Y0 = 0
        Z0 = 0
        d2_d3 = []
        for i in range(len(tablica)):
            d2, d3 = T_grs80.odl_2d_3d(X0, Y0, Z0, tablica[i][0], tablica[i][1],tablica[i][2])
            d2_d3.append([d2, d3])
        d2_d3 = np.array(d2_d3)
        print(d2_d3)
    
    elif opcja == 'G':
        print('\nObliczenie kąta azymutu i kąta elewacji')
        print(' Azymut  |  Elewacja  ')
        X0 = 0
        Y0 = 0
        Z0 = 0
        az_el = []
        for i in range(len(tablica)):
            az, el = T_grs80.azym_elew(X0, Y0, Z0, tablica[i][0], tablica[i][1],tablica[i][2])
            az_el.append(az, el)
        az_el = np.array(az_el)
        print(az_el)
    else:
        print('\nPodano złą opcję, wybierz jeszcze raz.')

    

#wsp_out = np.hstack([blh, u92])    
#np.savetxt("wsp_out.txt", wsp_out, delimiter=',  ', fmt = ['%10.2f', '%10.2f', '%10.3f', '%10.3f', '%10.3f'], header = 'Konwersja współrzędnych geodezyjnych \\ Aleksandra Skolimowska \n blh      układ 1992')
  