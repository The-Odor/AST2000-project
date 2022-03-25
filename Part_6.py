# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 21:31:10 2017

@author: theod
"""
from GOOD_BIG_SPECIAL_PROJECT import Planet_Class, print_time
from time import time

time_start = time()

scan = Planet_Class()
#    scan.show_spectrum()
#    print "modelling"
#    scan.model_spectrum()
#    print "modelled"

    #
    #    G = scan.G
    #    M = scan.syst.mass[3] * 1.989e30
    #
    #    satt_plan = array([-1502303.28257189, -18720571.26153056])
    #    awe  = array([3714.7509967 ,  -307.302641468])
    #    awer = norm(awe)
    #    x = vec_abs(satt_plan)
    #    v = sqrt((G*M)/x)
    #    print v, awer
    #
    #    print norm(array([-1512351.80387289, -18697059.86392402]))
    #    print norm(array([-8.36388178e+06,  -1.03401866e+08]))


#    scan.xyz_to_angular()
scan.land()
#print scan.parachute
#    scan.simulate_orbit()



print_time(time_start)
