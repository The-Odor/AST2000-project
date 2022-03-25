# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 15:03:38 2017

@author: Bruker
"""
from time import time
from GOOD_BIG_SPECIAL_PROJECT import Rocket_Class, print_time
time1 = time()
#from ast2000solarsystem_27_v5 import AST2000SolarSystem
#syst = AST2000SolarSystem(20229)
#syst.check_planet_positions
#raise Exception

L, T, N, m, k = 1E-6, 1E5, int(1E5), 2*1.67E-27, 1.381E-23
delt, n = 1E-9, int(1E3)
dt, G   = delt/float(n), 6.67408E-11
SAT     = 1100

rocket = Rocket_Class(L, T, N, m, k, SAT, delt, n, dt, G, 20229)
rocket.engine_init_posvel()
rocket.engine_posvel(short = True)
rocket.engine_find_momentum(hole = True, hole_size = L/2.)
rocket.find_acc() 



print "\nparticle escapees", rocket.partic.sum()

print "momentum escapees", rocket.moment.sum()

V = rocket.Vel_part_mean_func(T, m, k)
print "Change in speed", rocket.d_speed
print "average velocit", V
print "expected escape", ((V*delt)/(L))*N/2.
print "\n\n\n"

#rocket.TakeOff(graph = True, boxes = 6.07e12, fuel_ini = 3741)
rocket.TakeOff(graph = True, boxes = 1e13, fuel_ini = 7205)
print rocket.fuel_post_launch

print_time(time1)

raise Exception
rocket.TakeOff_test(SAT_override = True, SAT_ = 2100)
print "Running rocket.test_launch(), wish me lack"
rocket.test_launch()




