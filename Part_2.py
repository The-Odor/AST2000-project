# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 16:01:58 2017

@author: Bruker
"""
from time import time
from GOOD_BIG_SPECIAL_PROJECT import System_Class
time1 = time()


plan = System_Class()
plan.move_system(graph = False, make_new = True, save_file = True, time = 20, dt = 1e-5, steps_saved = "own") #Theodor, they saved 100 000 btw
raise Exception("you generated the system, good job! Now you can run the GOOD_BIG_SPECIAL_PROJECT.py*")
plan.make_xml(make_new = False)
plan.check_energy()
plan.read_stars()
plan.plot_rotavel()


time2 = time()

timespan = int(time2-time1)
time_hours = timespan/(60**2)
time_minutes = (timespan/60)%60
time_seconds = (timespan%60)
print """\n\n\nTime taken: {0} hours
            {1} minutes
            {2} seconds""".format(time_hours, time_minutes, time_seconds)






