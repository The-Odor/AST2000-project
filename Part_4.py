# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 14:37:47 2017

@author: theod
"""
from time import time
from numpy import pi; tau = 2*pi
from GOOD_BIG_SPECIAL_PROJECT import Orientation_Class, print_time
time1 = time()
p1 = 640
p2 = 480

orient = Orientation_Class()
a = orient.get_projection(alpha = 70, pro_amo = 360, pixels = (p2, p1)) #(360L, 640L, 480L, 3L)

for i in range(len(a)):
    break
    i += 1
    print "now saving image", i
    orient.save_img(a[i-1], "pictures\img%.3i.png"%i)


samples = ["0000", "0200", "0435", "0911", "1400", "1900"]
for samp in samples:
    break
    orient.get_direction("Skynet samples\sample{}.png".format(samp), show_data = False)
phile = "find_orient.png"
orient.get_direction(phile, show_data = True)

    


phi = (339.344868 * tau/360., 216.765374 * tau/360.)
#delam = (-0.018491838551, -0.018272021659)
#orient.get_velocity(phi, delam)
delam = ( 0.073343844630,-0.080813638011)
delam = (-0.080813638011, 0.073343844630)
delam = ( 0.027060240534, 0.049484564636)
orient.get_velocity(phi, delam)

from numpy import load
infile = open("pos.npy", "rb")
pos = load(infile)
infile.close()
t = 0.000037

y = orient.get_position_numeric(pos, t)
print y



print_time(time1)
