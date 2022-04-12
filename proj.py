# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 13:30:56 2022

@author: Ola
"""
import math
import numpy as np

class Transformacje:
    
    def __init__(self, model: str = "wgs84"):
        '''
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2

        '''

        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b) / self.a
        self.ecc2 = (2 * self.flattening - self.flattening ** 2)

    def xyz2blh_hirvonen(self, X, Y, Z):
        '''
        Algorytm Hirvonena - algorytm służący do transformacji współrzędnych 
        ortokartezjańskich X, Y, Z na współrzędne geodezyjne fi (szerokosć), 
        lambda (długosć), h (wysokosć). Jest to proces iteracyjny. W wyniku 
        3-4-krotnego powtarzania procedury można przeliczyć współrzędne na 
        poziomie dokładności 1 cm.

        Parameters
        ----------
        X : FLOAT
            - wartosc X w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        Y : FLOAT
            - wartosc Y w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        Z : FLOAT
            - wartosc Z w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]

        Returns
        -------
        fi : FLOAT
            - wartosć szerokosci geodezyjnej, liczba zmiennoprzecinkowa [stopnie]
        lam : FLOAT
            - wartosć dlugosci geodezyjnej, liczba zmiennoprzecinkowa [stopnie]
        h : FLOAT
            - wartosć wysokosci elipsoidalnej, liczba zmiennoprzecinkowa [m]

        '''
        r = math.sqrt(X**2 + Y**2)
        fi_n = math.atan(Z/(r*(1-self.ecc2)))
        eps = 0.000001/3600 *math.pi/180 # radiany
        fi = fi_n*2
        while np.abs(fi_n - fi) > eps:
            fi = fi_n
            N = self.a/np.sqrt(1-self.ecc2*np.sin(fi_n)**2)
            h = r/np.cos(fi_n) -N
            fi_n = math.atan(Z/(r*(1-self.ecc2*(N/(N + h)))))
        lam = math.atan(Y/X)
        h = r/np.cos(fi_n) -N
        fi = math.degrees(fi)
        lam = math.degrees(lam)
        return fi, lam, h
    
    def blh2xyz_odwrotny(self, fi, lam, h):
        '''
        ALgorytm odwrotny do algorytmu Hirvonena. Funkcja wykonuje transformację 
        współrzędnych krzywoliniowych do układu współrzędnych ortokartezjańskich.

        Parameters
        ----------
        fi : FLOAT
            - wartosć szerokosci geodezyjnej, liczba zmiennoprzecinkowa [stopnie]
        lam : FLOAT
            - wartosć dlugosci geodezyjnej, liczba zmiennoprzecinkowa [stopnie]
        h : FLOAT
            - wartosć wysokosci elipsoidalnej, liczba zmiennoprzecinkowa [m]

        Returns
        -------
        X : FLOAT
            - wartosc X w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        Y : FLOAT
            - wartosc Y w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        Z : FLOAT
            - wartosc Z w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]

        '''
        fi = np.radians(fi)
        lam = np.radians(lam)
        N = self.a/math.sqrt(1-self.ecc2*math.sin(fi)**2)
        X = (N + h) * math.cos(fi) * math.cos(lam)
        Y = (N + h) * math.cos(fi) * math.sin(lam)
        Z = (N*(1-self.ecc2) + h) * math.sin(fi)
        return X, Y, Z
    
    def u1992(self, fi, lam):
        '''
        Funkcja wykonująca transformacje współrzednych krzywoliniowych do układu 
        płaskiego 1992. Układ współrzędnych 1992 jest to układ współrzędnych 
        płaskich prostokątnych oparty na odwzorowaniu Gaussa-Krügera dla 
        elipsoidy GRS80 w jednej dziesięciostopniowej strefie.

        Parameters
        ----------
        fi : FLOAT
            - wartosć szerokosci geodezyjnej, liczba zmiennoprzecinkowa [stopnie]
        lam : FLOAT
            - wartosć dlugosci geodezyjnej, liczba zmiennoprzecinkowa [stopnie]

        Returns
        -------
        x92 : FLOAT
            -  wartosc X w układzie 1992, liczba zmiennoprzecinkowa [m]
        y92 : FLOAT
            -  wartosc Y w układzie 1992, liczba zmiennoprzecinkowa [m]

        '''
        fi = np.radians(fi)
        lam = np.radians(lam)
        m_0 = 0.9993
        N = self.a/(np.sqrt(1-self.ecc2 * np.sin(fi)**2))
        t = np.tan(fi)
        e_2 = self.ecc2/(1-self.ecc2)
        n2 = e_2 * (np.cos(fi))**2
        
        lam_0 = math.radians(19)
        l = lam - lam_0
        
        A0 = 1 - (self.ecc2/4) - ((3*(self.ecc2**2))/64) - ((5*(self.ecc2**3))/256)   
        A2 = (3/8) * (self.ecc2 + ((self.ecc2**2)/4) + ((15 * (self.ecc2**3))/128))
        A4 = (15/256) * (self.ecc2**2 + ((3*(self.ecc2**3))/4))
        A6 = (35 * (self.ecc2**3))/3072 
        
    
        sig = self.a * ((A0*fi) - (A2*np.sin(2*fi)) + (A4*np.sin(4*fi)) - (A6*np.sin(6*fi)))
        
        x = sig + ((l**2)/2) * N *np.sin(fi) * np.cos(fi) * (1 + ((l**2)/12) * ((math.cos(fi))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((l**4)/360) * ((math.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*math.cos(fi)) * (1 + ((((l**2)/6) * (math.cos(fi))**2) * (1-t**2+n2)) +  (((l**4)/(120)) * (math.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
        
        x92 = round(m_0*x - 5300000, 3)
        y92 = round(m_0*y + 500000, 3)
        return x92, y92 
    
    def u2000(self, fi, lam):
        '''
        Funkcja wykonująca transformacje współrzednych krzywoliniowych do układu 
        płaskiego 2000. Układ współrzędnych 2000 jest to układ współrzędnych 
        płaskich prostokątnych oparty na odwzorowaniu Gaussa-Krügera dla elipsoidy 
        GRS 80 w czterech trzystopniowych strefach o południkach osiowych 15°E, 
        18°E, 21°E i 24°E, oznaczonych kolejno numerami – 5, 6, 7 i 8. Skala 
        długości odwzorowania na południkach osiowych wynosi m0 = 0,999923.

        Parameters
        ----------
        fi : FLOAT
            - wartosć szerokosci geodezyjnej, liczba zmiennoprzecinkowa [stopnie]
        lam : FLOAT
            - wartosć dlugosci geodezyjnej, liczba zmiennoprzecinkowa [stopnie]

        Returns
        -------
        x00 : FLOAT
            -  wartosc X w układzie 2000, liczba zmiennoprzecinkowa [m]
        y00 : FLOAT
            -  wartosc Y w układzie 2000, liczba zmiennoprzecinkowa [m]

        '''
        fi = np.radians(fi)
        lam = np.radians(lam)
        m = 0.999923
        N = self.a/math.sqrt(1-self.ecc2*math.sin(fi)**2)
        t = np.tan(fi)
        e_2 = self.ecc2/(1-self.ecc2)
        n2 = e_2 * (np.cos(fi))**2
        
        lam = math.degrees(lam)
        if lam>13.5 and lam <16.5:
            s = 5
            lam0 = 15
        elif lam>16.5 and lam <19.5:
            s = 6
            lam0 = 18
        elif lam>19.5 and lam <22.5:
            s = 7
            lam0 = 21
        elif lam>22.5 and lam <25.5:
            s = 8
            lam0 = 24
            
        lam = math.radians(lam)
        lam0 = math.radians(lam0)
        l = lam - lam0
    
        A0 = 1 - (self.ecc2/4) - ((3*(self.ecc2**2))/64) - ((5*(self.ecc2**3))/256)   
        A2 = (3/8) * (self.ecc2 + ((self.ecc2**2)/4) + ((15 * (self.ecc2**3))/128))
        A4 = (15/256) * (self.ecc2**2 + ((3*(self.ecc2**3))/4))
        A6 = (35 * (self.ecc2**3))/3072 
        
    
        sig = self.a * ((A0*fi) - (A2*np.sin(2*fi)) + (A4*np.sin(4*fi)) - (A6*np.sin(6*fi)))
        
        x = sig + ((l**2)/2) * N *np.sin(fi) * np.cos(fi) * (1 + ((l**2)/12) * ((math.cos(fi))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((l**4)/360) * ((math.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*math.cos(fi)) * (1 + ((((l**2)/6) * (math.cos(fi))**2) * (1-t**2+n2)) +  (((l**4)/(120)) * (math.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
        
        x00 = round(m * x, 3) 
        y00 = round(m * y + (s*1000000) + 500000, 3)
         
        return x00, y00
    
    
    
    def NEU(self, X0, Y0, Z0, X, Y, Z):
        '''
        Funkcja wykonuje transformację współrzędnych ortokartezjańskich do układu 
        współrzędnych topocentrycznych. Układ topocentryczny jest to układ współrzędnych 
        ze środkiem znajdującym się w miejscu obserwacji. 

        Parameters
        ----------
        X0 : FLOAT
            - wartosc X punktu początkowego w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        Y0 : FLOAT
            - wartosc Y punktu początkowego w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        Z0 : FLOAT
            - wartosc Z punktu początkowego w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        X : FLOAT
            - wartosc X w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        Y : FLOAT
            - wartosc Y w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        Z : FLOAT
            - wartosc Z w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]

        Returns
        -------
        N : FLOAT
            DESCRIPTION.
        E : FLOAT
            DESCRIPTION.
        U : FLOAT
            DESCRIPTION.

        '''
        dX = X - X0
        dY = Y - Y0
        dZ = Z- Z0
            
        fi, lam, h = self.xyz2blh_hirvonen(X, Y, Z)
            
        R = np.array([[-np.sin(fi)*np.cos(lam), -np.sin(fi)*np.sin(lam), np.cos(fi)],
                      [-np.sin(lam), np.cos(lam), 0],
                      [np.cos(fi)*np.cos(lam), np.cos(fi)*np.sin(lam), np.sin(fi)]])
        
        d = np.array([[dX],
                     [dY],
                     [dZ]])
        NEU = R*d
        
        n = NEU[0,0]
        e = NEU[1,0]
        u = NEU[2,0]
        return n, e, u
            
            
    def odl_2d_3d(self, X0, Y0, Z0, X, Y, Z):
        '''
        Funkcja oblicza odległosć 2D i 3D na podstawie współrzędnych w układzie 
        ortokartezjańskim.

        Parameters
        ----------
        X0 : FLOAT
            - wartosc X punktu początkowego w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        Y0 : FLOAT
            - wartosc Y punktu początkowego w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        Z0 : FLOAT
            - wartosc Z punktu początkowego w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        X : FLOAT
            - wartosc X w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        Y : FLOAT
            - wartosc Y w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        Z : FLOAT
            - wartosc Z w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]

        Returns
        -------
        d2 : FLOAT
            - odległosć w przestrzeni dwuwymiarowej, liczba zmiennoprzecinkowa [m]
        d3 : FLOAT
            - odległosć w przestrzeni trójwymiarowej, liczba zmiennoprzecinkowa [m]

        '''
        dX = X - X0
        dY = Y - Y0
        dZ = Z - Z0
        
        d2 = math.sqrt(dX**2 + dY**2)
        d3 = math.sqrt(dX**2 + dY**2 + dZ**2)
        return d2, d3
    

    
    def azym_elew(self, X0, Y0, Z0, X, Y, Z):
        '''
        Funkcja oblicza kąt azymutu i kąt elewacji na podstawie współrzędnych w 
        układzie ortokartezjańskim.

        Parameters
        ----------
        X0 : FLOAT
            - wartosc X punktu początkowego w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        Y0 : FLOAT
            - wartosc Y punktu początkowego w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        Z0 : FLOAT
            - wartosc Z punktu początkowego w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        X : FLOAT
            - wartosc X w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        Y : FLOAT
            - wartosc Y w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]
        Z : FLOAT
            - wartosc Z w układzie ortokartezjańskim, liczba zmiennoprzecinkowa [m]

        Returns
        -------
        azymut : FLOAT
            - wartosć kąta azymutu, liczba zmiennoprzecinkowa [stopnie]
        elewacja : FLOAT
            - wartosć kąta elewacji, liczba zmiennoprzecinkowa [stopnie]

        '''
 
        N,E,U = self.NEU(X0, Y0, Z0, X, Y, Z)
        
        hz = np.sqrt(E**2 + N**2)
        el = np.sqrt(E**2 + N**2 + U**2)
        azymut = math.atan2(E, N)
        if azymut < 0:
            azymut = azymut + 2*np.pi
        
        elewacja = math.atan2(U, hz)
        
        azymut =  np.degrees(azymut)
        elewacja = np.degrees(elewacja)
        return azymut, elewacja
    
    

