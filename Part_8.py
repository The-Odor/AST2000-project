# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 17:42:14 2017

@author: theod
"""
v = 0.902 #c
Lplan = 500
Lship = 1200

#%% 3A

tA = 0
tB = 8.5 #milliseconds

tC = (tA + Lplan)/(1-v)
print tC #correct!

tD = (tB + Lplan)/(1+v)
print tD