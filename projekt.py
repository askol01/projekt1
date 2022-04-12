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
wsp_out = f'Konwersja współrzędnych geodezyjnych \\ Aleksandra Skolimowska\n'

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
        wsp_out += f'       fi          |          lam            |          h \n'
        blh = []
        for i in range(len(tablica)):
            b, l, h = T_grs80.xyz2blh_hirvonen(tablica[i][0], tablica[i][1], tablica[i][2])
            blh.append([b, l, h])
        blh = np.array(blh)
        for i in range(len(blh)):
            wsp_out += f'{blh[i][0]}  |  {blh[i][1]} |   {blh[i][2]} \n'
        print(blh)
        
    elif opcja == 'B':
        print('\nTransformacja współrzędnych geodezyjnych fi, lambda, h na współrzędne ortokartezjańskie X, Y, Z')
        print(' X  |  Y   |   Z')
        wsp_out += f'        X          |         Y          |         Z \n'
        XYZ = []
        for i in range(len(tablica)):
            b, l, h = T_grs80.xyz2blh_hirvonen(tablica[i][0], tablica[i][1], tablica[i][2])
            X, Y, Z = T_grs80.blh2xyz_odwrotny(b, l, h)
            XYZ.append([X, Y, Z])
        XYZ = np.array(XYZ)
        for i in range(len(XYZ)):
            wsp_out += f'{XYZ[i][0]}  |  {XYZ[i][1]} |   {XYZ[i][2]} \n'
        print(XYZ)
        
    elif opcja == 'C':
        print('\nTransformacja współrzednych geodezyjnych do układu 1992')
        print(' x92  |  y92  ')
        wsp_out += f'     x92      |       y92 \n'
        u92 = []
        for i in range(len(tablica)):
            b, l, h = T_grs80.xyz2blh_hirvonen(tablica[i][0], tablica[i][1], tablica[i][2])
            x92, y92 = T_grs80.u1992(b, l)
            u92.append([x92, y92])
        u92 = np.array(u92)
        for i in range(len(u92)):
            wsp_out += f'{u92[i][0]}  |  {u92[i][1]} \n'
        print(u92)
        
    elif opcja == 'D':
        print('\nTransformacja współrzednych geodezyjnych do układu 2000')
        print(' x00  |  y00  ')
        wsp_out += f'     x00      |       y00  \n'
        u00 = []
        for i in range(len(tablica)):
            b, l, h = T_grs80.xyz2blh_hirvonen(tablica[i][0], tablica[i][1], tablica[i][2])
            x00, y00 = T_grs80.u2000(b, l)
            u00.append([x00, y00])
        u00 = np.array(u00)
        for i in range(len(u00)):
            wsp_out += f'{u00[i][0]}  |  {u00[i][1]} \n'
        print(u00)
    
    elif opcja == 'E':
        print('\nTransformacja współrzędnych ortokartezjańskich X, Y, Z na współrzędne topocentryczne N, E, U')
        print(' N  |  E  |  U')
        wsp_out += f'        N           |           E           |          U  \n'
        NEU = []
        X0 = 0
        Y0 = 0
        Z0 = 0
        for i in range(len(tablica)):
            n,e,u = T_grs80.NEU(X0, Y0, Z0, tablica[i][0], tablica[i][1],tablica[i][2])
            NEU.append([n,e,u])
        NEU = np.array(NEU)
        for i in range(len(NEU)):
            wsp_out += f'{NEU[i][0]}  |  {NEU[i][1]} |   {NEU[i][2]} \n'
        
        print(NEU)
    
    elif opcja == 'F':
        print('\nObliczenie odległosci 2D i 3D')
        print(' 2D  |  3D  ')
        wsp_out += f'       2D         |         3D  \n'
        X0 = 0
        Y0 = 0
        Z0 = 0
        d2_d3 = []
        for i in range(len(tablica)):
            d2, d3 = T_grs80.odl_2d_3d(X0, Y0, Z0, tablica[i][0], tablica[i][1],tablica[i][2])
            d2_d3.append([d2, d3])
        d2_d3 = np.array(d2_d3)
        for i in range(len(d2_d3)):
            wsp_out += f'{d2_d3[i][0]}  |  {d2_d3[i][1]} \n'
        print(d2_d3)
    
    elif opcja == 'G':
        print('\nObliczenie kąta azymutu i kąta elewacji')
        print(' Azymut  |  Elewacja  ')
        wsp_out += f'       Azymut        |        Elewacja  \n'
        X0 = 0
        Y0 = 0
        Z0 = 0
        az_el = []
        for i in range(len(tablica)):
            az, el = T_grs80.azym_elew(X0, Y0, Z0, tablica[i][0], tablica[i][1],tablica[i][2])
            az_el.append([az, el])
        az_el = np.array(az_el)
        for i in range(len(az_el)):
            wsp_out += f'{az_el[i][0]}  |  {az_el[i][1]} \n'
        print(az_el)
    else:
        print('\nPodano złą opcję, wybierz jeszcze raz.')

    
zapis_plik = open('wsp_out.txt', 'a')
zapis_plik.write(wsp_out)
zapis_plik.close()
