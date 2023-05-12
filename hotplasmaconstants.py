#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 11:28:27 2022

@author: Bodhi Biswas
"""

#contants
c=2.998e8
e=1.60217662e-19
eps = 8.8541878176e-12 
me = 9.10938356e-31
u = 1.660539040e-27  # unified mass unit
#mi = 1.00794*u
mi = 2.01410178*u    # Deuterium
#mim = 39.948*u      # Ar (39.948)
#mim = 1.00794*u           # H2
mim = 12.011*u           # C12
Ze = -1.
Zi = 1.
Zim = 6.
mu0 = 1.256637e-6
impfrac=0.00
lambda_res = 50 #used for lin-liu CD model
theta_res = 51 #odd integer. used for lin-liu CD model
