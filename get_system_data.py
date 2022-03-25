# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 20:34:55 2017

@author: theod
"""
from ast2000solarsystem_27 import AST2000SolarSystem

seed = 20229

instance = AST2000SolarSystem(seed)

function = raw_input("What do you want from SS today? ")


if "," in function:
    functions = function.split(",")
    for function in functions:

        operation = "instance." + function
        print        
        try:
            print eval(operation)
        except:
            print eval(operation + "()")

            
else:
        operation = "instance." + function
        print
        try:
            print eval(operation)
        except:
            print eval(operation + "()")

