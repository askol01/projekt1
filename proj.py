# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 13:30:56 2022

@author: Ola
"""
import math
import numpy as np

class Transformacje:
    
    def __init__(self, model: str = "wgs84"):

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
        N = self.a/math.sqrt(1-self.ecc2*math.sin(fi)**2)
        X = (N + h) * math.cos(fi) * math.cos(lam)
        Y = (N + h) * math.cos(fi) * math.sin(lam)
        Z = (N*(1-self.ecc2) + h) * math.sin(fi)
        return(X, Y, Z)
    
    def u1992(self, fi, lam):
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
         
        return(x00, y00)
    
    
    def fi_lam2neu(self, fi, lam):
        n = np.array([-np.sin(fi)*np.cos(lam), -np.sin(fi)*np.sin(lam), np.cos(fi)])#dY
        e = np.array([-np.sin(lam), np.cos(lam), 0])#dX
        u = np.array([np.cos(fi)*np.cos(lam), np.cos(fi)*np.sin(lam), np.sin(fi)])#dZ
        return n, e, u
    
    def NEU(self, X0, Y0, Z0, X, Y, Z, typ: str = "XYZ"):
        if typ == "XYZ":
            dX = X - X0
            dY = Y - Y0
            dZ = Z- Z0
            
            fi, lam, h = self.xyz2blh_hirvonen(X, Y, Z)
            
            n, e, u = self.fi_lam2neu(fi, lam)
            
            N = n * dY
            E = e * dX
            U = u * dZ
        return N, E, U
            
            
    def odl_2d(self, X1, Y1, X2, Y2):
        dX = X2 - X1
        dY = Y2 - Y1
        s = math.sqrt(dX**2 + dY**2)
        return s
    
    def odl_3d(self, X1, Y1, Z1, X2, Y2, Z2):
        dX = X2 - X1
        dY = Y2 - Y1
        dZ = Z2 - Z1
        s = math.sqrt(dX**2 + dY**2 + dZ**2)
        return s
    
    def azymut(self, X1, Y1, X2, Y2):
        dX = X2 - X1
        dY = Y2 - Y1
        
        A = np.arctan(dY/dX)
        Az = np.degrees(A)
        
        if dX < 0 and dY < 0:
            Azymut = 180 + Az
        elif dX > 0 and dY > 0:
            Azymut = Az
        elif dX < 0 and dY > 0:
            Azymut = 180 - Az
        elif dX > 0 and dY < 0:
            Azymut = 360 - Az
        elif dX == 0 and dY > 0:
            Azymut = 90
        elif dX < 0 and dY == 0:
            Azymut = 180
        elif dX == 0 and dY < 0:
            Azymut = 270
        elif dX > 0 and dY == 0:
            Azymut = 0
         
        d = math.sqrt((dX**2) + (dY**2))
        
        return (Azymut, d)
    
    

#np.genfromtxt()