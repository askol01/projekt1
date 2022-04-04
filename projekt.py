# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 13:10:35 2022

@author: Ola
"""

import numpy as np
from proj import *

plik = "wsp_inp.txt"
tablica = np.genfromtxt(plik, delimiter=',', skip_header = 4)

T = Transformacje()

