# -*- coding: utf-8 -*-
"""
Created on Tue Nov 07 00:55:08 2017

@author: theod
"""
from time import time
from GOOD_BIG_SPECIAL_PROJECT import Travel_Class, print_time, vec_abs, rotate_angle_2d, Sattelite_Class, Rocket_Class
from numpy import array, sqrt, pi, dot
tau = 2*pi

time1 = time()


trav = Travel_Class()
trav.set_stuff(3.27142313, debug = False)
print
position = array([-0.224963283026 , 0.878248735469])
velocity = array([-6.02607206493 , 9.64662891945])
#trav.Space_align(array([-0.23622308617 , 0.899294214964]), array([-5.62262766726 , 10.525217190]), time = 3.347)

trav.new_satt_launch(3.349, position, velocity)
print trav.min_dist

#r = 0.00296162375604
G = trav.G
M = trav.syst.mass[3]

satt_plan = array([-0.247213219019 , 0.922308864884]) - array([-0.246899030046 , 0.921693446724])
x = vec_abs(satt_plan)
v = sqrt((G*M)/x)
#print v, "speed"
print "distance away from planet:", x
satt_plan = satt_plan/vec_abs(satt_plan)

boost = v * rotate_angle_2d(satt_plan, tau/4 )#+ tau/8) * 0.5
plan_vel = array([-4.12497774614 , -1.10262753071])
sat_vel  = array([-5.24957617374 , 11.5518114734])

actu_boost = boost - sat_vel + plan_vel
print actu_boost, "actual boost"





satt_plan = array([ -0.935672284616 , 0.0826882008398  ]) - array([-0.935856435569 , 0.0824711688509 ])
x = vec_abs(satt_plan)
v = sqrt((G*M)/x)
print v, "speed"
print "distance away from planet:", x
satt_plan = satt_plan/vec_abs(satt_plan)

boost = v * rotate_angle_2d(satt_plan, tau/4. )#+ tau/8.) * 0.5
print dot(boost, satt_plan)
print vec_abs(boost)/v, "aaaaaaaaaaaaa"
plan_vel = array([-0.283614826025 , -4.32845512003])
sat_vel  = array([0.493210702456 , -5.79653309944])

actu_boost = boost - sat_vel + plan_vel
print actu_boost, "actual boost"





trav.new_satt_launch(3.35, position, velocity, False)
satell = array([-0.936589097109 , 0.0774060137645])
planet = array([-0.93612981763 , 0.0781418275745])
print vec_abs(planet-satell)/x

#from numpy import linspace
#for i in linspace(3.5, 3.62, 13):
#    trav.new_satt_launch(i, position, velocity, False)


print_time(time1)




#launch
#boost 3.27243400 2.84814952 3.94035904
#
#boost 3.284612913 -4.42649216  2.32511658
#boost 3.345507478 0.43570066 -0.90009163
#boost 3.34 0.90148639  0.43280745
#
#boost 3.3494215 1.80290047 -10.30902573

#MOTHAFUCKN PIXEL PERFECT,                 but apparently too far away :(


#launch
#boost 3.27243400 2.84814952 3.94035904
#
#boost 3.284612913 -4.42649216  2.32511658
#boost 3.345507478 0.43570066 -0.90009163
#boost 3.34 0.90148639  0.43280745
#boost 3.345 0.39297672  0.91954842
#
#boost 3.349 1.98120954 -11.8138806
#
# WE FUCKING NAILED IT EXCEPT NOT BECAUSE I MISSED IT BY LIKE  NOTHING GRRR



#launch
#boost 3.27243400 2.84814952 3.94035904
#
#boost 3.284612913 -4.42649216  2.32511658
#boost 3.345507478 0.43570066 -0.90009163
#boost 3.34 0.90148639  0.43280745
#boost 3.345 0.39297672  0.91954842
#boost 3.347 0.11638053  0.9932047
#
#boost 3.349 0.42281539 -13.01271977

#IT'S A WORKING, IT'S A WOOOOORKIIIINNGGGGG

