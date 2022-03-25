# -*- coding: utf-8 -*-
"""
Created on Thu Oct 05 22:19:53 2017

@author: theod
"""
from time import time
from GOOD_BIG_SPECIAL_PROJECT import Sattelite_Class, Rocket_Class, print_time
time1 = time()


satt = Sattelite_Class()
L, T, N, m, k = 1E-6, 1E5, int(1E5), 2*1.67E-27, 1.381E-23
delt, n = 1E-9, int(1E3)
dt, G   = delt/float(n), 6.67408E-11
SAT     = 1100
seed    = 20229
rock = Rocket_Class(L, T, N, m, k, SAT, delt, n, dt, G, seed)

rock.engine_init_posvel()
rock.engine_posvel(short = True)
rock.engine_find_momentum(hole = True, hole_size = L/2.)
rock.find_acc()
rock.TakeOff(graph = False)
print "now running satt class"
satt.find_solar()
print satt.solar_panel_size
satt.launch = rock.launch
satt.start  = rock.launch_point
satt.fuel   = rock.fuel_post_launch
satt.SAT, satt.k, satt.T, satt.m = rock.SAT, rock.k, rock.T, rock.m
print "now running find_sat_launch"
satt.satt_launch(start = 3.27243400, time = 3.4, dt = 1E-4, graph = True)
satt.find_satt_launch(final_graph = False, tries_ = 100, inclines = 7, dt_ = 1E-4)
print satt.iniboost
raise Exception
satt.satt_launch(start = 0, time = 0.5, dt = 1E-4, graph = True, iniboost = "own")

timespan = int(time()-time1)
time_hours = timespan/(60**2)
time_minutes = (timespan/60)%60
time_seconds = (timespan%60)
print """\n\n\nTime taken: {0} hours
        {1} minutes
        {2} seconds""".format(time_hours, time_minutes, time_seconds)

