# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 14:40:21 2017

@author: Bruker
"""

from numpy  import uint8, sqrt, zeros, zeros_like, linspace, full, pi, save, load, \
                   array, ndarray, arccos, arcsin, arctan, cos, sin, meshgrid, nan,\
                   inf, where, dot, cross, prod, argmin, nanargmin, exp, float64
from ast2000solarsystem_27 import AST2000SolarSystem
from numpy.random          import normal, uniform
from numpy.linalg          import norm
from scipy.interpolate     import interp1d
from scipy.constants       import speed_of_light
from scipy                 import e
from os                    import path
from time                  import time
from PIL                   import Image
from collections           import Sequence
from psutil                import virtual_memory
from mpl_toolkits.mplot3d  import Axes3D
from matplotlib.pyplot     import plot, show, legend, title, draw, pause, ion, xlim,\
                                  ylim, savefig, xlabel, ylabel, axes, scatter




globalseed = 20229
tau        = 2*pi


def vec_abs(x, axis_ = "all"):
    if axis_ == "all":
        return sqrt(pow(x,2).sum())
    else:
        return sqrt(pow(x,2).sum(axis = axis_))

def rotate_angle_2d(x, angle):
    x, y = x
    return array([x*cos(angle) - y*sin(angle), x*sin(angle) + y*cos(angle)])


def normalize(x):
    return x/norm(x)


def print_time(time_ini, time_fin = "own", simple = False):
    if time_fin == "own": time_fin = time()
    timespan = int(time_fin-time_ini)
    time_hours = timespan/(60**2)
    time_minutes = (timespan/60)%60
    time_seconds = (timespan%60)
    if not simple:
        print
        print """Time taken: {0} hours
            {1} minutes
            {2} seconds""".format(time_hours, time_minutes, time_seconds)
    else:
        print "Time: {0}.{1}.{2}".format(time_hours, time_minutes, time_seconds)



def prog_bar(var, var_max, step = 100, short = True, sign = "|"):
    var += 1
    if short:
        step = step/10
    if var_max/step == 0:
        check = int(float(var_max)/step)
    else:
        check = var_max/step
    if var%(check) == 0:
        maxx = 100
        perc = int(100.*(var/float(var_max)))
        perc_ = perc
        if short:
            perc = int(perc/10)
            maxx = int(maxx/10)
        bars = sign*(perc)
        spcs = " "*(maxx-(perc))
        print "[{0}{1}]".format(bars, spcs), perc_, "%"
#APPARENTLY YOU CAN MAKE PRINT TAKE KEYWORDS TO ACTUALLY MAKE THIS THING WORK, WAAAAT



def read_table(file_, horisontal_dimension, excess_lines = 0, excess_words_initial = 0, read_every_nth_line = 1):

    def find_len(file_name):
        with open(file_name, "r") as data:
            for i, l in enumerate(data):
                pass
        return i + 1

    data_len = int(find_len(file_)-excess_lines)
    table    = zeros((horisontal_dimension, data_len)) #Creates table

    with open(file_, "r") as data:
        for _ in range(excess_lines):
            data.readline()#Reads first unecessary lines
        for line_nr in range(data_len)[::read_every_nth_line]: #Reads through the lines
            words = data.readline().split() #Splits the lines into words
            for n in range(horisontal_dimension): #Fills table
                table[n][line_nr] = float(words[n + excess_words_initial])

    return table


class Rocket_Class():
    def __init__(self, L, T, N, m, k, SAT, delt, n, dt, G, seed):
        self.L = L;   self.T  = T;  self.N    = N
        self.m = m;   self.k  = k;  self.delt = delt
        self.n = n;   self.dt = dt; self.G    = G
        self.seed = seed
        self.SAT  = SAT

        self.home   = AST2000SolarSystem(seed)
        self.homema = self.home.mass[0] * 2E30
        self.homera = self.home.radius[0] * 1E3
        self.syst   = self.home

        Vel_esc_func  = lambda G, M, r: sqrt((2*float(G)*M)/r)
        Vel_part_mean_func = lambda T, m, k: sqrt((k*T)/float(m))

        self.Vel_esc_func = Vel_esc_func
        self.Vel_part_mean_func = Vel_part_mean_func

        self.Vel_esc = Vel_esc_func(G, self.homema, self.homera)
        self.Vek_part_mean = Vel_part_mean_func(T, m, k)



    def engine_init_posvel(self, file_name = ("engine\PositionalVelocity.npy", "engine\VelocityParticle.npy")):
        if not path.isfile(file_name[0]) or not path.isfile(file_name[1]):

            N = self.N;   L = self.L;   k = self.k
            T = self.T;   m = self.m;   n = self.n

            v        = normal(0, (sqrt(k*T/float(m))), (N, 3))
            r_ini    = uniform(0, L, (N, 3))
            r        = zeros((N, 3, n))
            r[:,:,0] = r_ini

            self.v, self.r = v, r[:,:,:]



    def engine_posvel(self, prog = True, short = True, file_name = "engine\MomentumVector.npy"):
        if not path.isfile(file_name):

            print "\nNow calculating positional vector"
            L   = self.L; n   = self.n
            vel = self.v; rad = self.r
            dt  = self.dt

            for i in range(1, n):
                j = i - 1
                radi             = rad[:,:,j].copy()
                radi[radi >= L] -= L
                radi[radi <= 0] += L
                rad[:,:,i]       = radi.copy() + vel.copy()*dt

                if prog:
                    prog_bar(i, n, short = short)

            self.vel = vel
            self.rad = rad

            #outfile_part = open(file_name[0], "wb")
            #outfile_velo = open(file_name[1], "wb")
            #save(outfile_part, rad)
            #save(outfile_velo, vel)
            #outfile_part.close()
            #outfile_velo.close()


        #else:
            #infile_part = open(file_name[0], "rb")
            #infile_velo = open(file_name[1], "rb")
            #self.rad = load(infile_part)
            #self.vel = load(infile_velo)
            #infile_part.close()
            #infile_velo.close()



    def engine_find_momentum(self, dim = 0, hole = False, hole_size = 0, side = True, short = True, prog = True, file_name =  "engine\MomentumVector.npy"): #side = True means top, False > bottom
        m    = self.m; L    = self.L; delt = self.delt


        if not path.isfile(file_name):

            print "\nNow calculating momentum vector"
            vel  = self.vel[:, dim ]
            rad  = self.rad[:, :, :]
            n    = self.n; m    = self.m
            L    = self.L; delt = self.delt
            moment = zeros(n)
            partic = zeros_like(moment)
            dim = int(dim)


            if hole:
                D = L/2
                hole_half = hole_size/2.

                if side:
                    for i in range(n):
                        moment_temp = vel[(rad[:,dim,i] >= L)\
                                       &(rad[:,dim-1,i] > D - hole_half)\
                                       &(rad[:,dim-1,i] < D + hole_half)\
                                       &(rad[:,dim-2,i] > D - hole_half)\
                                       &(rad[:,dim-2,i] < D + hole_half)]
                        if prog:
                            prog_bar(i, n, short = short)


                        partic[i]   = len(moment_temp)
                        moment[i]   = moment_temp.sum()


                else:
                    for i in range(n):
                        moment_temp = vel[(rad[:,dim,i] <= L)\
                                       &(rad[:,dim-1,i] > D - hole_half)\
                                       &(rad[:,dim-1,i] < D + hole_half)\
                                       &(rad[:,dim-2,i] > D - hole_half)\
                                       &(rad[:,dim-2,i] < D + hole_half)].sum()
                        if prog:
                            prog_bar(i, n, short = short)


                        partic[i]   = len(moment_temp)
                        moment[i]   = moment_temp.sum()

            else:
                if side:
                    for i in range(n):
                        moment_temp = vel[rad[:,dim,i] >= L]

                        partic[i]   = len(moment_temp)
                        moment[i]   = moment_temp.sum()
                        if prog:
                            prog_bar(i, n, short = short)

                else:
                    for i in range(n):
                        moment_temp = vel[rad[:,dim,i] <= L]

                        partic[i]   = len(moment_temp)
                        moment[i]   = moment_temp.sum()
                        if prog:
                            prog_bar(i, n, short = short)

            self.moment     =  moment*m*2
            self.pressure   = (moment*m*2)/((L**2)*delt)
            self.partic     = partic
            self.partic_per = partic.sum()/delt
            self.force_box  = moment.sum()/delt

            outfile = open(file_name, "wb")
            save(outfile, [moment, partic])
            outfile.close()

        else:
            infile = open(file_name, "rb")
            moment, partic = load(infile)
            infile.close()

            self.moment     =  moment*m*2
            self.pressure   = (moment*m*2)/((L**2)*delt)
            self.partic     = partic
            self.partic_per = partic.sum()/delt
            self.force_box  = moment.sum()/delt





    def find_acc(self, fuel = False):
        moment   = self.moment
        delt     = self.delt
        SAT      = self.SAT
#        fuel_ini = self.fuel

        force    = moment.sum()/delt


        if fuel:
            pass
        else:
            acc = force/SAT
            speed_change = acc*delt

            self.acc       = acc
            self.force_box = force
            self.d_speed   = speed_change # not accounting for fuel


    def TakeOff_test_bleugh(self, time_min = 20, steps = 100, tol_min = 5): #variables: Fuel, Box
        #Also, TakeOff_test does not work yet, maybe I'll get back to it later
        G  = self.G;      force_box = self.force_box
        Mp = self.homema; dist_ini  = self.homera
        m  = self.m;      SAT = self.SAT
        Ve = self.Vel_esc_func


        time  = (time_min + tol_min) * 60
        tol   = tol_min  * 60



        amo_box = 40 #100
        amo_ful = 20 # 30


        fuel_max = 4000

        boxes_ = linspace(1E10, 1e14, amo_box) #in the 1.5E13 area
        fuel_  = linspace(1000,  fuel_max, amo_ful) #1000 to 3000

        dt = 1./time*steps

        dist_tests = zeros(shape = (amo_box, amo_ful, time*steps))
        dist_tests[:,:,0] = full((amo_box, amo_ful), dist_ini)
        velo_tests = zeros_like(dist_tests)
        results    = zeros(shape = (amo_box, amo_ful, 2))
        bofute     = zeros(shape = (amo_box, amo_ful, 2))

        print "\n\n\n\nNOW TESTING A FUCKTON OF SHIT 'hopefully, it'll work this time'"
        print Ve(G, Mp, dist_ini)

        for bo in range(len(boxes_)):
            fuel_loss = self.partic_per * m * bo
            for fu in range(len(fuel_)):
                prog_bar((bo+1)*len(fuel_) + fu, len(boxes_)*len(fuel_))
                fuel_left = fuel_[fu]
                accel = (boxes_[bo]*force_box)/(SAT + fuel_left) - G*(Mp/dist_tests[bo,fu,0]**2)
                for i in range(1, time*steps):
                    velo_tests[bo,fu,i] = velo_tests[bo,fu,i-1] + accel*dt

                    if dist_tests[bo,fu,i-1] <= 0 and velo_tests[bo,fu,i] <= 0:
                        velo_tests[bo,fu,i] = 0
                    elif fuel_left > 0:
                        dist_tests[bo,fu,i] = dist_tests[bo,fu,i-1] + velo_tests[bo,fu,i]*dt


                    fuel_left -= fuel_loss
                    accel      = (boxes_[bo]*force_box)/(SAT + fuel_left) - G*(Mp/dist_tests[bo,fu,i]**2)



                    if velo_tests[bo,fu,i] >= Ve(G, Mp, dist_tests[bo,fu,i]):
                        results[bo,fu,:] = 1e9, i/(steps)
                        bofute[bo,fu,:]  = boxes_[bo], fuel_[fu]
                        break





        results_time = results
        results1 = results_time - time
        results1[results1 < 0]*= -1
        print results1.min(), tol, time
        results_tol = results1[results1 < tol]
        bofute_tol  = bofute[results1 < tol]
        print results_tol.shape
        print results1.shape
        print bofute.shape
        print bofute_tol.shape
        print bofute_tol
        fuel_opt = bofute_tol.min()
        print fuel_opt
        print results1.min()/60.


        self.bofute = bofute[results1 == results1.min()]
        print self.bofute



    def TakeOff_test(self, time_min = 20, steps = 200, tol_min = 5, SAT_override = False, SAT_ = 1100):
        G  = self.G;      force_box = self.force_box
        Mp = self.homema; dist_ini  = self.homera
        m  = self.m;      SAT = self.SAT
        Ve = self.Vel_esc_func
        if SAT_override:
            SAT = SAT_

        dt = 1./steps
        time = time_min * 60
        tol  = tol_min  * 60

        amo_box = 100
        amo_ful = 30

        fuel_max = 8000
        boxs_max = 1e15

        boxes_ = linspace(1E10, boxs_max, amo_box) #in the 1.5E13 area
        fuel_  = linspace(1500,  fuel_max, amo_ful) #1000 to 3000

        boxes, fuel = meshgrid(boxes_, fuel_)
        fuel_ini = fuel.copy()

        dist = zeros_like(boxes); dist.fill(dist_ini)
        vels = zeros_like(dist )
        fins = zeros_like(dist ); fins.fill(inf)
        fins_done = zeros_like(dist); fins_done.fill(True)

        fuel_loss = self.partic_per * m * boxes * dt

        for i in range(steps*time): #introduce steps, dammit #steal them from actual launch function #DONE
            prog_bar(i, steps*time, short = True)
            t = i*dt
            accel = (boxes*force_box)/(SAT + fuel) - G*(Mp/dist**2)
            accel[fuel<=0] = -G*(Mp/dist[fuel<=0]**2)
            vels = vels + accel*dt #also introduce dt, luhloauhweuhw #DONE
            dist = dist + vels*dt


            vels[dist<dist_ini] = 0
            dist[dist<dist_ini] = dist_ini + 1e-5

            fuel -= fuel_loss
            """

            if abs(t-time) < tol:
                fins[where((vels >= Ve(G, Mp, dist)) & \
                     (fins_done))] \
                     = fuel[where((vels >= Ve(G, Mp, dist)) & \
                           (fins_done))]

                fins_done[where(vels >= Ve(G, Mp, dist))] = False
            """
            if abs(t-time) < tol:
                for i in range(len(fins)):
                    for j in range(len(fins[i])):
                        if vels[i,j] >= Ve(G, Mp, dist[i,j]) and \
                           fins_done[i,j]:
                               fins[i,j] = fuel[i,j]

                               fins_done[i,j] = False
                               #print "AHHHHHHHHHHHHHHHHH"

        best_fit = abs(fins).argmin()
        best_fit = array([best_fit/fins.shape[1], best_fit%fins.shape[1]])
        #args_ = where(less == less.min())

        print fins[best_fit[0], best_fit[1]]
        print best_fit
        print fuel_[best_fit[0]]
        print boxes_[best_fit[1]]



        # INSPIRATION:

        """
        for bo in range(len(boxes_)):
            fuel_loss = self.partic_per * m * bo
            for fu in range(len(fuel_)):
                prog_bar((bo+1)*len(fuel_) + fu, len(boxes_)*len(fuel_))
                fuel_left = fuel_[fu]
                accel = (boxes_[bo]*force_box)/(SAT + fuel_left) - G*(Mp/dist_tests[bo,fu,0]**2)
                for i in range(1, time*steps):
                    velo_tests[bo,fu,i] = velo_tests[bo,fu,i-1] + accel*dt

                    if dist_tests[bo,fu,i-1] <= 0 and velo_tests[bo,fu,i] <= 0:
                        velo_tests[bo,fu,i] = 0
                    elif fuel_left > 0:
                        dist_tests[bo,fu,i] = dist_tests[bo,fu,i-1] + velo_tests[bo,fu,i]*dt


                    fuel_left -= fuel_loss
                    accel      = (boxes_[bo]*force_box)/(SAT + fuel_left) - G*(Mp/dist_tests[bo,fu,i]**2)



                    if velo_tests[bo,fu,i] >= Ve(G, Mp, dist_tests[bo,fu,i]):
                        results[bo,fu,:] = 1e9, i/(steps)
                        bofute[bo,fu,:]  = boxes_[bo], fuel_[fu]
                        break
        """


#    def TakeOff(self, time_min = 20, steps = 1000, boxes = 9e11, fuel_ini = 850, graph = True):
    def TakeOff(self, time_min = 20, steps = 1000, boxes = 1e13, fuel_ini = 7206, graph = True):
        self.rocket_boxes = boxes
        self.rocket_fuel  = fuel_ini
        SAT = self.SAT; Mp = self.homema
        dist_ini  = self.homera; force_box = self.force_box
        fuel_loss = self.partic_per * self.m * boxes
        self.fuel_loss = fuel_loss
        G = self.G
        dt = 1./steps
        homema  = self.homema # I import self.homema twice because I know how to code :)
        esc_velf = self.Vel_esc_func

        not_first = True

        time = time_min * 60

        time_ar = linspace(0, time, time*steps)
        dist    = zeros(time*steps)
        velo    = zeros_like(dist)
        fuel    = zeros_like(dist)
        esc_vel = zeros_like(dist)
        dist[0] = dist_ini
        fuel[0] = fuel_ini

        esc_vel[0] = esc_velf(G, homema, dist[0])
        accel = (boxes*force_box)/(SAT + fuel[0]) - G*(Mp/dist[0]**2)

        for i in range(1, time*steps):
            velo[i] = velo[i-1] + accel*dt
            #if i == 1:
                #print velo[i]
                #print accel

            if dist[i-1] <= 0 and velo[i] <= 0:
                velo[i] = 0
                dist[i] = dist[i-1]
            elif fuel[i-1] > 0:
                dist[i] = dist[i-1] + velo[i]*dt


                fuel[i]    = fuel[i-1] - fuel_loss*(dt)

            accel      = (boxes*force_box)/(SAT + fuel[i]) - G*(Mp/dist[i]**2)

            if velo[i] >= esc_velf(G, homema, dist[i]) and not_first:
                print time_ar[i], "seconds to hit escape velocity"
                print fuel[i], "kg of fuel left when escaping"
                self.escape_time = time_ar[i]
                self.escape_point = dist[i]
                not_first = False

            esc_vel[i] = esc_velf(G, homema, dist[i])
        if graph:
            plot(time_ar, velo)
            plot(time_ar, fuel*10)
            plot(time_ar, esc_vel)
            legend(["Velocity", "Fuel*10", "Escape Velocity"], loc="best")
            title("Takeoff")
            show()

        #print "max speed  : ", velo.max()
        #print "excess fuel: ", fuel[-1]

        self.fuel_post_launch = fuel[-1]
        self.launch           = velo[-1]
        self.launch_point     = dist[-1]
        self.boxes            = boxes



    def SpaceMov_time(self, d_vel, dt):
        fuel      = self.fuel_post_launch
        fuel_loss = self.fuel_loss
        force_box = self.force_box
        SAT       = self.SAT
        boxes     = self.boxes

        change_vel = 0
        time       = 0

        while change_vel < d_vel:
            if  fuel <= 0:
                print "Sir, we ran outta fuel"
                return time, change_vel
            fuel -= fuel_loss
            force = force_box * boxes
            mass  = SAT + fuel
            change_vel += (force/float(mass))*dt
            time += dt

        self.space_fuel = fuel

        return time



    def SpaceMov(self, d_vel, mass = "own", fuel = "launch"):
        k = self.k;  T = self.T;  m = self.m
        vmean = sqrt(k*T/float(m))
        if fuel == "launch":
            fuel = self.fuel_post_launch
        if mass == "own":
            mass = self.SAT + fuel

        return mass/((vmean/d_vel) + 1.)



    def test_launch_USELESS(self): #this ting is not-work, fuck that
        infile = open("planet_positions.npy", "rb")
        orbit_pos_, orbit_vel_, time = load(infile)
        orbit_pos = interp1d(time, orbit_pos_)
        infile.close()


        syst = self.syst
        launch_point = self.escape_point
        x0,  y0  = syst.x0[0],  syst.y0[0]
        vx0, vy0 = syst.vx0[0], syst.vy0[0]
        rad      = syst.radius / 149597871.
        launch_dir = array([vx0, vy0])/sqrt(vx0**2 + vy0**2)
        init_sat_pos = (x0 + rad[0]*launch_dir[0], y0 + rad[0]*launch_dir[1])
        syst.engine_settings(self.force_box, 9E11, self.partic_per, 850, self.escape_time, init_sat_pos, 0)

        #not sure if this is how this part works:
        initial_velocity_from_rotation = (tau*rad[0])/(syst.period[0]*24*60*60)
        omega = 1./syst.period[0]*24*60*60
        angle = omega*self.escape_time*tau
        launch_perp = array([x0, y0])/sqrt(x0**2 + y0**2)
        initial_velocity_from_rotation*= launch_perp

        #launch_dir = (launch_dir + launch_perp)/sqrt(launch_dir**2 + launch_perp**2)
        #launch_dir = (launch_dir + launch_perp)/sqrt((launch_dir**2).sum() + (launch_perp**2).sum())


        #launch_dir_ = launch_dir*launch_point + launch_perp*initial_velocity_from_rotation*self.escape_time
        #launch_dir = launch_dir_/sqrt(launch_dir_[0]**2 + launch_dir_[1]**2)

        l = launch_dir
        launch_dir_angle = array([l[0]*(cos(angle) - sin(angle)), l[1]*(cos(angle) + sin(angle))])


        pos_after_launch = launch_dir*launch_point + array([x0, y0])
        pos_after_launch = launch_dir*launch_point + orbit_pos(self.escape_time/(365*24*60*60))[:,0]*1.496e11 + initial_velocity_from_rotation*self.escape_time*1000.
        pos_after_launch*= 1./1.
        print pos_after_launch.shape
        pos_after_launch = launch_dir_angle*launch_point + orbit_pos(self.escape_time/(365*24*60*60))[:,0]*1.496e11
        syst.mass_needed_launch(pos_after_launch, test=True)


    def test_launch(self, time_launch = 0): #COME ONE, JUST KEEP YOUR UNITS CONSTANT
        #apparently I don't need to do stuff, just let it ruuun
        test_ = True
        time_launch_ = time_launch
        for i in range(2):
            if i == 0:
                #continue
                pass

            #t = 3.21691548 #Just as a reminder

            time_launch = 0
            if i == 1:
                time_launch = time_launch_
            #time_launch += 0.000105
            #low  : dude, fuck, I don't know
            #high : u kno it, bro

            infile = open("planet_positions.npy", "rb")
            orbit_pos_, orbit_vel_, time = load(infile)
            orbit_pos = interp1d(time, orbit_pos_)
            orbit_vel_ = zeros_like(orbit_pos_)
            orbit_vel_[:,:,:-1] = (orbit_pos_[:,:,1:] - orbit_pos_[:,:,:-1])/(time[1:] - time[:-1])
            #print orbit_vel_[:,:,-1]
            orbit_vel = interp1d(time, orbit_vel_)
            infile.close()

            def vel_good(orb, t, dt):
                mindt = t - dt
                maxdt = t + dt
                space = 2 * dt
                if   t < dt and t != 0:
                    mindt = 0
                    maxdt = t + t
                    space = t + t
                elif t < dt and t == 0:
                    syst = self.syst
                    vx0, vy0 = syst.vx0, syst.vy0
                    return array([vx0, vy0])
                return (orb(maxdt)-orb(mindt))/(space)

            syst = self.syst
            #x0,  y0  = syst.x0[0],  syst.y0[0]
            #vx0, vy0 = syst.vx0[0], syst.vy0[0]

            x0,  y0  = orbit_pos(time_launch)[:,0]
            #vx0, vy0 = orbit_vel(time_launch)[:,0]
            print self.syst.vx0, self.syst.vy0
            vx0, vy0 = vel_good(orbit_pos, time_launch, 1e-7)[:,0]
            print vx0, vy0
            m_to_AU = 1/1.4959787070e11
            rad     = syst.radius * 1000. * m_to_AU
            print "rad", rad
            launch_dir = array([vx0, vy0])/sqrt(vx0**2 + vy0**2)
            launch_point = self.escape_point * launch_dir
            init_sat_pos = (x0 + rad[0]*launch_dir[0], y0 + rad[0]*launch_dir[1])
            syst.engine_settings(self.force_box, self.rocket_boxes, self.partic_per, self.rocket_fuel, self.escape_time, init_sat_pos, time_launch)

            print "period", syst.period[0]
            omega = 1./(syst.period[0]*24*60*60)
            rad_ = rad/m_to_AU
            initial_velocity_from_rotation = tau*rad_[0]*omega
            #inx, iny = init_sat_pos
            launch_perp = array([x0, y0])/sqrt(x0**2 + y0**2)
            plan_time = time_launch + self.escape_time/(365.*24*60*60)

            orbit = orbit_pos(plan_time)[:,0]

            launch_point = launch_point * m_to_AU
            initial_velocity_from_rotation*= self.escape_time*launch_perp * m_to_AU
            move = launch_point - initial_velocity_from_rotation

            pos_after_launch = orbit + move
            #pos_after_launch*= m_to_AU

            print pos_after_launch, "position"
            print orbit_pos(plan_time)[:,0], "planet x-y"
            print orbit_vel(plan_time)[:,0], "planet v"
            print launch_dir, "dir"
            print launch_perp, "perp"
            print pos_after_launch - orbit_pos(plan_time)[:,0], "move"
            print initial_velocity_from_rotation
            print launch_point
            print (launch_point - initial_velocity_from_rotation)
            syst.mass_needed_launch(pos_after_launch, test=test_)
            test_ = False



































class System_Class():
    def __init__(self):
        syst   = AST2000SolarSystem(20229)
        starM  = syst.star_mass
        starR  = syst.star_radius
        planN  = syst.number_of_planets

        homeM  = syst.mass[0]
        homeR  = syst.radius[0]

        self.syst  = syst
        self.starM = starM
        self.starR = starR
        self.planN = planN

        self.home  = (homeM, homeR)

        self.G         = 4*pi**2 #Using astronomical units, so differs from previous G
        self.AU =  AU  = 1.496e11 #meters
        self.yr =  yr  = 60*60*24*365.25 #seconds
        self.c         = speed_of_light * (yr/AU) #convert to astronomical units by years
        self.alpha_wvl = 656.3

        self.has_moved_system     = False
        self.Simulated_new_system = False



    def read_file(self, file_, excess_lines, horisontal_dimension, excess_words_initial, print_dims = False):
        def find_len(file_name):
            with open(file_name, "r") as data:
                for i, l in enumerate(data):
                    pass
            return i + 1

        def read_table(file_, excess_lines, horisontal_dimension, excess_words_initial, print_dims):

            data_len = int(find_len(file_)-excess_lines)
            table    = zeros((horisontal_dimension, data_len)) #Creates table

            if print_dims:
                print "Read table of dimensions %g, %g"%(data_len, horisontal_dimension)

            with open(file_, "r") as data:
                for _ in range(excess_lines):
                    data.readline()#Reads first unecessary lines
                for line_nr in range(data_len): #Reads through the lines
                    words = data.readline().split() #Splits the lines into words
                    for n in range(horisontal_dimension): #Fills table
                        table[n][line_nr] = float(words[n + excess_words_initial])

            return table


        return read_table(file_, excess_lines, horisontal_dimension, excess_words_initial, print_dims)



    def plot_system(self):
        x0, y0 = self.syst.x0, self.syst.y0
        mass   = self.syst.mass
        starM  = self.starM

        Mcx = ((mass*x0).sum())/((len(x0)+1)*(mass.sum()+starM))
        Mcy = ((mass*y0).sum())/((len(y0)+1)*(mass.sum()+starM))


        for x, y in zip(x0, y0):
            plot(x, y, "o")

        plot(Mcx, Mcy, "^")

        show()



    def move_system(self, time = 5, dt = 1E-6, graph = False, file_name = "planet_positions.npy", make_new = False, steps_saved = 50000, save_file = True): #time = 10 gives 6 rotations for planet with largest orbit
        steps = int(time/float(dt))
        if save_file:
            make_new = True
        if path.isfile(file_name) and not make_new and not graph:
            print "file with positions already exists, not making new one"
            infile = open(file_name, "rb")
            position_time = load(infile)
            infile.close() #Shape is [position, time, velocity],  [x-, y-axis],  [planet],  [time-step]

        else:
            planN = self.planN
            G     = self.G
            syst  = self.syst

            starM = syst.star_mass

            x0,  y0  = syst.x0,  syst.y0
            vx0, vy0 = syst.vx0, syst.vy0


            r_hatt    = zeros((2, planN)) #r_hatt, den har ei retning
                                          #ei retning har r_hatt
                                          #Og har r_hatt ei retning
                                          #Så er det ein skalar
                                          #
                                          #Til melodien til "min hatt"
            len_pos   = zeros((planN))
            velocity  = zeros((2, planN, steps))
            position  = zeros((2, planN, steps))
            accel     = zeros((2, planN))
            velocity[:,:,0] = vx0, vy0
            position[:,:,0] = x0, y0
            len_pos   = sqrt(position[0,:,0]**2 + position[1,:,0]**2)

            r_hatt[:,:]= -position[0,:,0], -position[1,:,0]
            tot_cel    = G*starM/(len_pos**3)
            accel[:,:] = tot_cel * r_hatt
            velocity[:,:,0]  -= accel*0.5*dt

            if graph:
                ylim(-2.5, 2.5)
                xlim(-2.5, 2.5)
                plot(0, 0, "^")

                legend_list = ["Star"] + ["planet {0}".format(i) for i in range(planN)]
                legend(legend_list)
                ion()

            print "now moving celestial objects"

            for i in range(1, steps):
                prog_bar(i, steps)
                j = i-1

                velocity[:,:,i] = velocity[:,:,j] + dt * accel
                position[:,:,i] = position[:,:,j] + velocity[:,:,i]*dt

                len_pos     = sqrt(position[0,:,i]**2 + position[1,:,i]**2)
                r_hatt[:,:] = -position[0,:,i], -position[1,:,i]
                tot_cel     = G*starM/(len_pos**3)
                accel[:,:]  = tot_cel * r_hatt


                if graph and i%(steps/5000) == 0:
                    for x, y in zip(position[0,:,i], position[1,:,i]):
                        plot(x, y, "o")
                        draw()
                    pause(0.001)




            if save_file and steps_saved != "own":
                SSyst = AST2000SolarSystem(20229)
                print "starting check_planet_positions:"
                print "position   :", len(position[:,:,::int(steps/float(steps_saved))][0][0])
                print "steps      :", steps_saved
                SSyst.check_planet_positions(position[:,:,::int(steps/steps_saved)], time, steps_saved/time)#steps/time)
                #SSyst.check_planet_positions(position, time, steps/time)

            elif save_file and steps_saved == "own":
                print "starting check_planet_positions:"
                print "position   : all"#, len(position[:,:,::int(steps/float(steps_saved))][0][0])
                print "steps      : all"#, steps_saved
                self.syst.check_planet_positions(position, time, steps/float(time))#steps/time)




            if steps_saved != "own":
                infile = open(file_name, "rb")
                rray = load(infile)
                print rray[1].shape
                position_time = array([rray[0], velocity[:,:,::int(steps/steps_saved)], rray[1]])
                infile.close()
            else:
                infile = open(file_name, "rb")
                rray = load(infile)
                print rray[1].shape
                position_time = array([rray[0], velocity, rray[1]])
                infile.close()


            if save_file:
                outfile = open(file_name, "wb")
                save(outfile, position_time)
                outfile.close()


#            with open(file_name, "wb") as outfile:  # They apparently already had a save function I could use :(
#                save(outfile, (velocity[:,:,::int(steps/steps_saved), linspace(0, time, steps_saved)]))

            self.Simulated_new_system = True

        self.planet_orbits     = position_time
        self.has_moved_system  = True



    def make_xml(self, make_new = False):
        if not path.isfile("MCAst_win/mcast/data/video_20229.xml") or make_new or self.Simulated_new_system:
            syst = self.syst
            positions, vel_unused, times = self.planet_orbits
            syst.orbit_xml(positions[:,:,::10], times[::10])
        else:
            print "XML already exists, not making new one"



    def plot_rotavel(self):
        position, velocity, time = self.planet_orbits


        Line_of_sight = array([0, 1]) #Along y-axis, makes for simpler calculations
        Pec_vel       = 0
        incline_angle = pi/2
        planet        = 2


        plot(position[0,planet,:], position[1,planet,:])
        title("planet Tuiheng orbit")
        show()
        #this data means that I can just plot the velocity in the y-axis directly

        steps = 50000 #data doesn't drown
        plot(time[:steps], velocity[1,planet,:steps])
        title("planet Tuiheng radial velocity")
        show()

        noise = 0.2*max(velocity[1, planet, :])
        vel_noise = velocity + normal(0, noise, time.shape)


        plot(time[:steps], vel_noise[1,planet,:steps])
        title("planet Tuiheng as it will probably be seen")
        show()



    def check_energy(self):
        position, velocity, time = self.planet_orbits
        starM = self.starM
        mass  = self.syst.mass
        G     = self.G


        red_mass = (mass*starM)/(mass + starM)
        red_mass_arr = zeros((len(red_mass), len(time)))
        M_arr        = zeros_like(red_mass_arr)
        for i in range(len(time)):
            red_mass_arr[:,i] = red_mass
            M_arr[:,i]        = mass + G

        print velocity.shape
        energy_kin_ = 0.5*red_mass_arr*(velocity[0,:,:]**2 + velocity[1,:,:]**2)

        energy_kin = energy_kin_.sum(axis = 0)


        energy_pot_ = -(G*M_arr*red_mass_arr)/sqrt(position[0,:,:]**2 + position[1,:,:]**2)
        energy_pot = energy_pot_.sum(axis = 0)

        for kin, pot, plan in zip(energy_kin_, energy_pot_, range(len(energy_kin_))):
            plot(time, kin)
            title("Planet {0}: Kinetic".format(plan))
            show()
            plot(time, pot)
            title("Planet {0}: Potential".format(plan))
            show()
            plot(time, kin + pot)
            title("Planet {0}: Total".format(plan))
            show()

        plot(time, energy_pot + energy_kin)
        title("Total energy of all planets")
        #ylim(-1, 1)
        show()



    def read_stars(self):
        c = self.c  ;G = self.G
        alpha_wvl = self.alpha_wvl

        stars = [1.70, 1.63, 1.24, 4.49, 1.72]
        #Line 1: time of observation
        #line 2: observed wavelength (reminded to use doppler)
        #line 3: measured flux of light relative to the maximum flux for hte given star


        tables = []

        for i, l in enumerate(stars):
            file_name  = "Remote stars/star%g_%3.2f.txt"%(i, l)
            make_array = array(self.read_file(file_name, 0, 3, 0, print_dims = True))
            tables.append(make_array)


        for star in tables:
            time, wvl, light = star

            wvl_dif = wvl/alpha_wvl - 1
            relvel  = wvl_dif * c #relative velocity

            plot(time, relvel)
            title("Star {0}: Rotational velocity".format(i))
            show()

            plot(time, light, "r")
            title("Star {0}: light flux".format(i))
            show()

            pec_vel = relvel.sum()/len(relvel)
            print pec_vel
            print "-------------------------------------------------------------------------------------"


        time, unused, light = tables[2]
        botlim, toplim = 1350,1370
        plot(time[botlim:toplim], light[botlim:toplim])
        title("Star that shows promise for planet with radius tells")
        show()

        #use planet 2

        time, wvl, unused = tables[2]
        wvl_dif = wvl/alpha_wvl - 1
        relvel  = wvl_dif * c
        botlim, toplim = 3250,3750
        plot(time[botlim:toplim], relvel[botlim:toplim])
        title("Star that shows promise for planet")
        show()
        botlim, toplim = 500, 1000
        plot(time[botlim:toplim], relvel[botlim:toplim])
        title("Star that shows promise for planet")
        show()

        top = 3600
        bot = 600

        planmass = lambda M, v, P, G: (M**(2/3.)*v*P**(1/3.)) / ((2*pi*G)**(1/3.))

        P = (top-bot)/365
        v = max(abs(relvel))
        M = stars[2]

        min_mass_of_star = planmass(M, v, P, G)
        print min_mass_of_star
































class Sattelite_Class(Rocket_Class):
    def __init__(self):
        syst   = AST2000SolarSystem(20229)
        starM  = syst.star_mass
        starR  = syst.star_radius
        starT  = syst.temperature

        self.syst  = syst
        self.starM = starM
        self.starR = starR
        self.starT = starT
        self.planN = syst.number_of_planets
        self.G     = 4*pi**2



    def find_solar(self, aim = 40., planet = 3, efficiency = 0.12):
        starR, starT, syst = self.starR, self.starT, self.syst

        AU_to_km = 149597871   #because position is given in Au's, and starR is given in km's
        stebolt     = 5.670367e-8

        effect = lambda T, R, r: stebolt*(T**4)*(R**2)/(float(r)**2)

        infile = open("planet_positions.npy", "rb")
        position, velocity, time= load(infile)
        infile.close()

        position, velocity = position[:, planet,:], velocity[:, planet,:]
        position*= AU_to_km
        radius   = sqrt((position**2).sum(axis = 0))
        dist     = radius.sum()/len(radius)

        size_effect   = effect(starT, starR, dist)
        size_absorbed = size_effect * efficiency
        size_solar    = aim/size_absorbed

        self.radiation_on_planet = size_effect
        self.solar_panel_size    = size_solar


        plan_rad = syst.radius[planet]
        plante   = pi*plan_rad**2

        plan_temp_func = lambda T, R1, R2, A: T*sqrt((R1/R2)/sqrt(4)) #To future self, this is now the last thing you came up with; you updated it
        plan_temp      = plan_temp_func(starT, starR, dist, plante)

#        print "Planetary temp calc  :", plan_temp

        self.plan_temp_func = plan_temp_func
        self.plan_temp = plan_temp

        #This horseshit forgot the area of the planet for the outflux of blackbody energy, so fuck me, amirite
        #You deleted the horseshit, but why not keep notes because aaaaaaaaaaaahhhhhhh

#        plan = [328, 276, 111, 155, 177, 96, 130][planet]
#        print "Planetary temp calc to actual:", plan_temp/plan



    def SpaceMov(self, d_vel, mass = "own", fuel = "launch"):
        k = self.k;  T = self.T;  m = self.m
        vmean = sqrt(k*T/float(m))
        if fuel == "launch":
            fuel = self.fuel
        if mass == "own":
            mass = self.SAT + fuel

        if d_vel != 0:
            return mass/((vmean/d_vel) + 1.)
        elif d_vel == 0:
            return 0
        else:
            print "bruh, you just called SpaceMov with a d_vel that is not a number or sumthin, _change it_"




    def satt_launch(self, start = 0, time = 1, dt = 1E-3, boosts = [], iniboost = 4.86193223827, iniboost_add = True, graph = False,\
                    launch_dir = "own", start_pos_vel = "start", iniboost_override = False, iniboost_print = False, progress = False):
        if isinstance(start_pos_vel, str):
            start_pos_vel_str = True
        else:
            start_pos_vel_str = False
        if isinstance(launch_dir, str):
            launch_dir_str = True
        else:
            launch_dir_str = False
        time_ = time
        steps = int((time_-start)/float(dt))

        iniboost_state = False

        infile = open("planet_positions.npy", "rb")
        orbit_pos_, orbit_vel_, time = load(infile)
        orbit_pos = interp1d(time, orbit_pos_)
        orbit_vel = interp1d(time, orbit_vel_)
        infile.close()

        steps = int(steps)
        planN = self.planN
        G     = self.G
        syst  = self.syst
        starM = self.starM
        km_to_AU    = 1 / 149597871.
        sec_to_yr   = 1 / (60*60*24*365.)
        if start_pos_vel_str:
            self.start *= km_to_AU/1000.
            self.launch*= km_to_AU/(1000.*sec_to_yr)
        else:
            self.start  = 0
            self.launch = 0
        planM = syst.mass
        if start_pos_vel_str:
            homeR = syst.radius[0] * km_to_AU
        else:
            homeR = 0
        if iniboost != "own":
            iniboost_state = True
            if not iniboost_add:
                iniboost -= self.launch

        starM = syst.star_mass

        x0,  y0  = orbit_pos(start)[:,:]
        vx0, vy0 = orbit_vel(start)[:,:]


        initial_velocity_from_rotation = (tau*homeR)/(syst.period[0]/365.)
        if iniboost_state and not iniboost_add:
            iniboost  -= initial_velocity_from_rotation
        elif iniboost_state and iniboost_add:
            pass
        else:
            iniboost = 0

        initial_velocity_from_rotation*= array([x0[0], y0[0]])/sqrt(x0[0]**2 + y0[0]**2)

        r_hatt     = zeros(2)
        r_hatt_    = zeros(2)
        accel      = zeros(2)
        velocity   = zeros((2, steps))
        position   = zeros((2, steps))
        planaccel  = zeros((2, planN))
        if launch_dir_str:
            launch_dir = orbit_vel(start)[:,0]
        launch_dir = array(launch_dir)
        launch_dir = launch_dir / sqrt((launch_dir**2).sum())
        #print "iniboost", iniboost * launch_dir

        launch_go  = launch_dir * self.launch #ADDED SELF.LAUNCH FROM ROCKET CLASS, RUN BEFORE OR FACE WRATH OF ERROR U FUCK
        if start_pos_vel_str:
            velocity[:,0] = vx0[0], vy0[0]
            position[:,0] = x0[0],  y0[0]
        else:
            position[:,0] = start_pos_vel[0,:]
            velocity[:,0] = start_pos_vel[1,:]



        velocity[:,0]+= launch_go
        velocity[:,0]+= initial_velocity_from_rotation
        vel_boost     = iniboost*launch_dir#; self.fuel -= self.SpaceMov(d_vel = iniboost)
        self.iniboost = iniboost
        if iniboost_override:
            if iniboost_print:
                print "velocity overriden; applying iniboost as ", vel_boost - velocity[:,0]
            velocity[:,0] = array([0, 0])
        velocity[:,0]+= vel_boost
        position[:,0]+= launch_dir*(homeR + self.start)
        len_pos       = sqrt((position[:,0]**2).sum(axis = 0))

        r_hatt[:]     = array([x0[0]/len_pos, y0[0]/len_pos])
        tot_cel       = G*starM/(len_pos**2)
        accel[:]      = tot_cel * r_hatt

        if graph:
            ion()
            legend_ = ["Planet 0", "Planet 1", "Planet 2", "Planet 3", "Planet 4", "Planet 5", "Planet 6", ]
            legend(legend_)
            title("start time: {0} \nend time : {1}".format(start, time_))

        for i in range(planN):
            r_hatt_[:]    = array([(x0[i] - position[0,0])/len_pos, (y0[i] - position[1,0])/len_pos])
            len_pos       = sqrt(((array([x0[i], y0[i]]) - position[:,0])**2).sum(axis = 0))
            tot_cel       = G*planM[i]/(len_pos**2)
            accel[:]     += tot_cel * r_hatt_


        velocity[:,0]-= accel*0.5*dt

        """ Assuming this section is to initiate it
        min_dist = array([x0[3], y0[3]]) - array([x0[0], y0[0]])
        min_dist = sqrt((min_dist**2).sum())
        min_dist = array([min_dist, 0])
        and adding an improvement"""
        min_dist = array([inf, 0])

        #for k in range(planN):
        #    plot(orbit_pos_[0,k,:], orbit_pos_[1,k,:])

        for i in range(1, steps):
            if progress:
                prog_bar(i, steps)
            j = i-1
            t = i*dt + start

            velocity[:,i] = velocity[:,j] + dt * accel
            position[:,i] = position[:,j] + velocity[:,i]*dt

            len_pos     = sqrt(position[0,i]**2 + position[1,i]**2)
            r_hatt[:]   = -position[0,i], -position[1,i]
            tot_cel     = G*starM/(len_pos**3)
            accel[:]    = tot_cel * r_hatt

            """
            for k in range(planN):
                if i%int(steps/500.) == 0 and graph:
                    plot(orbit_pos(t)[0,k], orbit_pos(t)[1,k], "^")
                #z = int(round((float(i)/steps)*len(orbit_pos[0,k,:])))

                #r_hatt_[:] = -(orbit_pos[0,k,z] - position[0,i]), -(orbit_pos[1,k,z] - position[1,i])
                #dist       = array([orbit_pos[0,k,z], orbit_pos[1,k,z]]) - position[:,i]
                r_hatt_[:]  = -(orbit_pos(t)[0,k] - position[0,i]) -(orbit_pos(t)[1,k] - position[1,i])
                dist        = array([orbit_pos(t)[0,k], orbit_pos(t)[1,k]]) - position[:,i]

                len_pos    = sqrt((dist**2).sum(axis = 0))
                tot_cel    = G*planM[k]/(len_pos**3)
                accel[:]  += tot_cel * r_hatt_

                if len_pos < min_dist[0] and k == 3:
                    min_dist[:] = len_pos, t

            """
            planaccel      = array((-(orbit_pos(t)[0,:] - position[0,i]), -(orbit_pos(t)[1,:] - position[1,i])))
            #dist_arr      = array([orbit_pos(t)[0,:] - position[0,i], orbit_pos(t)[1,:] - position[0,i]]) - position[:,i] #potential if array can't lay together
            len_pos        = sqrt((planaccel**2).sum(axis = 0))
            tot_cel        = G*planM/(len_pos**3)
            planaccel[:,:]*= tot_cel, tot_cel
            accel[:]      += planaccel.sum(axis = 1)

            if steps < 500:
                check = 1
            else:
                check = int(steps/500.)

            if i%check == 0 and graph:
                for k in range(planN):
                    plot(orbit_pos(t)[0,k], orbit_pos(t)[1,k], "^")

                plot(position[0,i], position[1,i], "ro")
                draw()
                pause(0.0001)

            if len_pos[3] < min_dist[0]:
                min_dist[:] = len_pos[3], t




        self.min_dist = min_dist



    def find_satt_launch(self, start = 0, stop = 4, length = 1.1, tries_ = 150, inclines = 8, final_graph = True, dt_ = 1E-3):
        first = True
        min_ray = zeros(3)
        tries = linspace(start, stop, tries_)
        new_start, new_stop = start, tries[2]
        for q in range(inclines):
            for j in range(tries_):
                prog_bar(j + q*tries_, tries_*inclines)
                i = tries[j]
                self.satt_launch(start = i, time = i+length, dt = dt_, graph = False)
                if self.min_dist[0] < min_ray[0] or first:
                    first = False
                    min_ray[:] = self.min_dist[0], self.min_dist[1], new_start
                    new_start, new_stop = tries[j-1], tries[j+1]

            tries = linspace(new_start, new_stop, tries_)
        #very detailed shit gave: Distance, time to hit, time of launch [  1.36322000e-04   2.88876619e+00   2.79456619e+00]
        print "Distance, time to hit, time of launch", min_ray
        self.satt_launch(start = min_ray[2], time = min_ray[2] + length, graph = final_graph)








































class Orientation_Class(): #hehe, murder on the oriental express
    def __init__(self):
        self.syst = AST2000SolarSystem(globalseed)
        with open("himmelkule.npy", "rb") as infile:
            self.heavenorb = load(infile)

        self.AU    =    AU    = 1.496e11 #meters
        self.yr    =    yr    = 60*60*24*365.25 #seconds
        self.speed_of_light   = speed_of_light * (yr/AU) #convert to astronomical units by years



    def xy_to_pt(self, x, y, theta0, phi0):
        p = sqrt(x**2 + y**2)
        #c = 2*arctan(p/2)

        theta = pi/2 - arcsin(cos(2*arctan(sqrt(x**2 + y**2)/2))*cos(theta0) + \
                              (y*sin(2*arctan(sqrt(x**2 + y**2)/2))*sin(theta0))/sqrt(x**2 + y**2))

        phi   = phi0 + arctan(x*sin(2*arctan((sqrt(x**2 + y**2))/2))  \
                             /(sqrt(x**2 + y**2)*sin(theta0)*cos(2*arctan((sqrt(x**2 + y**2))/2)) - \
                              y*cos(theta0)*sin(2*arctan(sqrt(x**2 + y**2)/2))))

        return theta, phi



    def save_img(self, img, img_name = "Orientation.png"):
        if isinstance(img, ndarray):
            img = Image.fromarray(img)
        elif isinstance(img, Image.Image):
            pass
        else:
            raise Exception("Invalid datatype for save_img in Orient_Class")
        img.save(img_name) #Make new png



    def get_projection_eugh(self, theta = 0, phi = 0, alpha_ = 70, pixels = (640, 480), vary_theta = False, pro = 360):
        #I guess we're operating under degrees now! - maybe?
        if path.isfile("heavenly sphere.npy") and False:
            infile = open("heavenly sphere.npy", "rb")
            projections = load(infile)
            infile.close()
            return projections

        alpha = alpha_ * pi/180.
        theta-= alpha/2.
        phi  -= alpha/2.
        thetamax    = theta + alpha
        phimax      = phi   + alpha
        field_phi   = linspace(phi,   phimax,    pixels[0])
        field_theta = linspace(theta, thetamax,  pixels[1])
        pros        = linspace(0, tau, pro)
        temp, phi2d, tht2d = meshgrid(pros, field_phi, field_theta)
        phi2d      += temp

        phi2d = phi2d
        tht2d = tht2d

        if vary_theta: tht2d += temp
        phis        = linspace(phi,   phi + tau, pro)
        phi_dif     = phimax - phi
        syst = self.syst

        midphi = phi2d.sum(axis = 0)/len(phi2d[:,0,0])
        midtht = tht2d.sum(axis = 2)/len(tht2d[0,0,:])
        midtht = midtht[0,0]
        time00 = time()
        actutht_, actuphi_ = self.xy_to_pt(tht2d, phi2d, midtht, midphi)

        if vary_theta:
            projections = zeros((pro, pro, pixels[0], pixels[1]))
            thetas = linspace(theta, theta + tau, pro)%tau
            print thetas #needed a use to get rid of the bloody warning
            raise Exception("Vary_theta not built yet \n/NOT/ under construction")
            #probably gonna be memoryintensive, I'unno


        else:
            time_for_xytopt  = time() - time00
            time_for_ang2pix = 0
            print "Now projecting the Heavenly Sphere:"
            projections = zeros((pro, pixels[0], pixels[1], 3), dtype = uint8)
            for i in range(pro):
                for phi_ in range(pixels[0]):
                    for tht_ in range(pixels[1]):
                        prog_bar(i*pixels[0]*pixels[1] + phi_*pixels[1] + tht_, pro*pixels[0]*pixels[1])
                        time11 = time()

                        #midtht,  midphi  = (theta + thetamax)/2, 2*phis[i] + phi_dif/2.
                        #actutht, actuphi = self.xy_to_pt(field_theta[tht_], field_phi[phi_], midtht, midphi)
                        actutht, actuphi = actutht_[phi_,i,tht_], actuphi_[phi_,i,tht_]
                        time22 = time()

                        angfrompix       = syst.ang2pix(actutht%tau, actuphi%tau)
                        projections[i, phi_, tht_] = self.heavenorb[angfrompix][2:]
                        time33 = time()

                        time_for_xytopt  += time22-time11
                        time_for_ang2pix += time33-time22

                field_phi   += pi/180
                field_phi   %= tau
        print "My functions:"
        print_time(0, time_for_xytopt)
        print
        print "ang2pix (the sucky one >:( )"
        print_time(0, time_for_ang2pix)
        print


        if not path.isfile("heavenly sphere.npy"):
            outfile = open("heavenly sphere.npy", "wb")
            save(outfile, projections)
            outfile.close()


        return projections

    def get_projection(self, theta = pi/2., phi = 0, alpha = 70, pixels = (640, 480), vary_theta = False, pro = 360, pro_amo = 360):
        filename = "heavenly sphere.npy"
        if path.isfile(filename) and True:
            infile = open(filename, "rb")
            projections = load(infile)
            infile.close()
            return projections
        phis = linspace(phi + alpha/2.,   phi - alpha/2.,   pixels[0])
        thts = linspace(theta - alpha/2., theta + alpha/2., pixels[1])
        pros = linspace(1, pro, pro_amo)

        pros, thts, phis = meshgrid(pros, thts, phis)

        thts = thts * pi/180.
        phis = phis * pi/180.
        pros = pros * pi/180.
        #print phis.shape
        #raise Exception

        thts, phis = self.xy_to_pt(thts, phis, theta, pros) #Thetha


        projections = zeros((pro_amo, pixels[0], pixels[1], 3), dtype = uint8)


        for i in range(pro_amo):
            for phi_ in range(pixels[0]):
                for tht_ in range(pixels[1]):
                    prog_bar(i*pixels[0]*pixels[1] + phi_*pixels[1] + tht_, pro_amo*pixels[0]*pixels[1])

                    actutht, actuphi = thts[tht_,i,phi_], phis[tht_,i,phi_]

                    angfrompix = self.syst.ang2pix(actutht, actuphi) #note to self, removed %tau from both actu's for testing purposes
                    projections[i, phi_, tht_] = self.heavenorb[angfrompix][2:]

        outfile = open(filename, "wb")
        save(outfile, projections)
        outfile.close()

        return projections



    def get_direction(self, img, pro = (360, 360), show_data = False):
        if   isinstance(img, str):
            infile = Image.open(img, "r")
            img    = array(infile)
            infile.close()

        elif isinstance(img, ndarray):
            pass
        elif isinstance(img, Image.Image):
            img = array(img)

        else: raise Exception("img in function:'get_direction' is neither ndarray, Image.Image nor filename")

        propics_shape = [pro[1]] + [i for i in img.shape]
        propics = zeros(propics_shape, dtype=uint8)
        for i in range(pro[1]):
            propics[i,:,:,:] = img


        projections = self.get_projection(pro = pro, pixels=img.shape[:-1])

        variance    = (projections - propics)**2

        variance = variance.sum(axis = (1, 2, 3))
        print variance.shape

        angle_deg = variance.argmin() + 1
        angle_rad = angle_deg * pi/180.
        if show_data:
            plot(variance)
            show()
            print "Angle: %.3f radians"%angle_rad
            print "       %3.0f degrees"%angle_deg
            print

        self.angle = angle_rad
        return angle_rad



    def get_velocity(self, phi, delam):
        syst= self.syst
        syst.get_ref_stars()

        for arg, name in zip((phi, delam), ("phi", "delam")):
            if not isinstance(arg, Sequence):
                raise Exception("{} is not a sequence, gib moar numbers".format(name))
            elif len(arg) < 2:
                raise Exception("{} does not have enough elements".format(name))
            elif len(arg) > 2:
                print "{} has more than two elements, grabbing first two".format(name)
                if name == "phi":
                    phi = array(phi)
                    phi = phi[0:2]
                elif name == "delam":
                    delam = array(delam)
                    delam = delam[0:2]
        phi, delam = array(phi), array(delam)


        print "speed of light", self.speed_of_light

        vel     = lambda lam, lam_ref = 656.3, c = self.speed_of_light: (lam/float(lam_ref))*c
        actuvel = lambda p1, p2, v1, v2: array([ sin(p2)*v1 - sin(p1)*v2,\
                                                -cos(p2)*v1 + cos(p1)*v2])\
                                                /sin(p2-p1)

        #Introduce some way to get the angles and delams >:(

        phi = array([339.344868 * tau/360., 216.765374 * tau/360.])
        delam_base = array([-0.018491838551, -0.018272021659])

        show_vel_base = vel(delam_base)

        ref_vel = actuvel(phi[0], phi[1], show_vel_base[0], show_vel_base[1])

        show_vel = vel(delam)
        print ref_vel - actuvel(phi[0], phi[1], show_vel[0], show_vel[1])



    def get_position(self, t, r1, r2, r3):

        infile = open("planet_positions.npy", "rb")
        orbit_pos_, orbit_vel_, time = load(infile)
        orbit_pos = interp1d(time, orbit_pos_)
        infile.close()

        x, y = orbit_pos(t)[:,0:3]
        sattpos = array([5., 10.])
        r = [sqrt(((sattpos - array([x[i], y[i]]))**2).sum()) for i in range(3)]

        #x, y, r = array((x, y, r))
        x = array(x)
        y = array(y)
        r = array(r)
        n = x**2 + y**2 - r**2

        #c = (x[1] - x[0])/(y[0] - float(y[1]))
        #d = (n[0] - n[1])/(2*(y[0] - float(y[1])))

        """
        sqrtcomp = pow(x[2] + c*(y[2] - d), 2) + (c**2 + 1)*(-2*y[2]*d + d**2 + n[2])

        x1 = 2*((x[2] + c*(y[2] - d) + sqrt(sqrtcomp))\
             /(c**2 + 1))

        x2 = 2*((x[2] + c*(y[2] - d) - sqrt(sqrtcomp))\
             /(c**2 + 1))

        """

        """
        xtop = n[2] - n[0] - d*(y[2] - float(y[0]))
        xbot = 2*(x[2] - x[0]) + (y[2] - float(y[0]))

        xout = xtop/xbot

        y = xout*c + d
        """

        xtop1 = (n[1] - n[2])/(2*(y[1] - y[2]))
        xtop2 = (n[2] - n[0])/(2*(y[2] - y[0]))

        xbot1 = (n[2] - n[1])/(y[1] - y[2])
        xbot2 = (x[2] - x[0])/(y[2] - y[0])

        xout = (xtop1 - xtop2)/(xbot1 - xbot2)


        ytop  = n[2] - n[0] - 2*xout*(x[2] - x[0])
        ybot  = 2*(y[2] - y[0])

        yout = ytop/ybot




        return xout, yout


    def get_position_numeric(self, sattpos, t, n = 500000): # >:(
        infile = open("planet_positions.npy", "rb")
        orbit_pos_, orbit_vel_, time = load(infile)
        orbit_pos = interp1d(time, orbit_pos_)
        infile.close()

        x, y = orbit_pos(t)[:,0:3]
        #sattpos = array([5., 9.])
        #r = [sqrt(((sattpos - array([x[i], y[i]]))**2).sum()) for i in range(3)]
        r = sattpos[0:3]

        #x, y, r = array((x, y, r))
        x = array(x)
        y = array(y)
        r = array(r)

        #xs = linspace(x-r, x+r, 10000)
        xs = [linspace(x_-r_/2., x_+r_/2., n) for x_, r_ in zip(x, r)]
        xs = array(xs)
        ys = sqrt(r**2 - (xs.T - x)**2) + y
        ys = ys.T
        print xs.shape
        diff2 = array([1e25, 0, 0])
        diff3 = array([1e25, 0, 0])

        for i in range(len(xs[0])):
            prog_bar(i, len(xs[0]))
            radx1  = xs[0,i]
            rady1t = y[0] + sqrt(r[0]**2 - (radx1 - x[0])**2)
            #yt = y + sqrt(r**2 - (xs[:,i] - x)**2)
            rady1b = y[0] - sqrt(r[0]**2 - (radx1 - x[0])**2)


            len2t = abs(sqrt((x[1] - radx1)**2 + (y[1] - rady1t)**2)-r[1])
            len2b = abs(sqrt((x[1] - radx1)**2 + (y[1] - rady1b)**2)-r[1])
            len2  = min([len2t, len2b])
            top_big2 = False
            if len2t > len2b:
                top_big2 = True

            len3t = abs(sqrt((x[2] - radx1)**2 + (y[2] - rady1t)**2)-r[2])
            len3b = abs(sqrt((x[2] - radx1)**2 + (y[2] - rady1b)**2)-r[2])
            len3  = min([len3t, len3b])
            top_big3 = False
            if len3t > len3b:
                top_big3 = True

            if len2 < diff2[0]:
                if top_big2:
                    diff2[:] = len2, radx1, rady1b
                else:
                    diff2[:] = len2, radx1, rady1t

            if len3  < diff3[0]:
                if top_big3:
                    diff3[:] = len3, radx1, rady1b
                else:
                    diff3[:] = len3, radx1, rady1t

        return diff2[1:]


    def orient(self, img, phidelam, positionalstuff):
        self.get_direction(img)
        self.get_velocity(phidelam[0], phidelam[1])
        self.get_position(positionalstuff)



































class Travel_Class(Rocket_Class, Sattelite_Class):
    def __init__(self):
        L, T, N, m, k = 1E-6, 1E5, int(1E5), 2*1.67E-27, 1.381E-23
        delt, n = 1E-9, int(1E3)
        dt, G   = delt/float(n), 6.67408E-11
        SAT     = 1100
        Rocket_Class.__init__(self, L, T, N, m, k, SAT, delt, n, dt, G, 20229)

        Rocket_Class.engine_init_posvel(self)
        Rocket_Class.engine_posvel(self, short = True)
        Rocket_Class.engine_find_momentum(self, hole = True, hole_size = L/2.)
        Rocket_Class.find_acc(self)
        Rocket_Class.TakeOff(self, graph = False)

    def set_stuff(self, time_launch, debug = False):
        infile = open("planet_positions.npy", "rb")
        orbit_pos_, orbit_vel_, time = load(infile)
        orbit_pos = interp1d(time, orbit_pos_)
        #orbit_vel_ = zeros_like(orbit_pos_)
        #orbit_vel_[:,:,:-1] = (orbit_pos_[:,:,1:] - orbit_pos_[:,:,:-1])/(time[1:] - time[:-1])
        #orbit_vel = interp1d(time, orbit_vel_)
        infile.close()

        def vel_good(orb, t, dt):
            mindt = t - dt
            maxdt = t + dt
            space = 2 * dt
            if   t < dt and t != 0:
                mindt = 0
                maxdt = t + t
                space = t + t
            elif t < dt and t == 0:
                syst = self.syst
                vx0, vy0 = syst.vx0, syst.vy0
                return array([vx0, vy0])
            return (orb(maxdt)-orb(mindt))/(space)

        syst = self.syst

        x0,  y0  = orbit_pos(time_launch)[:,0]
        vx0, vy0 = vel_good(orbit_pos, time_launch, 1e-7)[:,0]
        m_to_AU = 1/1.4959787070e11
        rad     = syst.radius * 1000. * m_to_AU
        launch_dir = array([vx0, vy0])/sqrt(vx0**2 + vy0**2)
        launch_point = self.escape_point * launch_dir
        init_sat_pos = (x0 + rad[0]*launch_dir[0], y0 + rad[0]*launch_dir[1])
        syst.engine_settings(self.force_box, self.rocket_boxes, self.partic_per, self.rocket_fuel, self.escape_time, init_sat_pos, time_launch)

        omega = 1./(syst.period[0]*24*60*60)
        rad_ = rad/m_to_AU
        initial_velocity_from_rotation = tau*rad_[0]*omega

        launch_perp = array([x0, y0])/sqrt(x0**2 + y0**2)
        plan_time = time_launch + self.escape_time/(365.*24*60*60)

        orbit = orbit_pos(plan_time)[:,0]

        launch_point = launch_point * m_to_AU
        initial_velocity_from_rotation*= self.escape_time*launch_perp * m_to_AU
        move = launch_point - initial_velocity_from_rotation

        pos_after_launch = orbit + move
        syst.mass_needed_launch(pos_after_launch)

        syst.send_satellite("command files\Launch instructions.txt")

        if debug:
            orders = open("command files\Launch instructions.txt", "r")
            for line in orders:
                spl = line.split()
                if spl[0] == "orient":
                    check_time = spl[-1]
            orders.close()

            print orbit_pos(eval(check_time))[:,3]





    def Space_align(self, position, velocity, time, target_plan = 3, angle_ = tau/360, d_ang = 1e-1, inclines = 3, iniboost = 1):
        Sattelite_Class.__init__(self)
        infile = open("planet_positions.npy", "rb")
        orbit_pos_, orbit_vel_, time_ = load(infile)
        orbit_pos = interp1d(time_, orbit_pos_)
        infile.close()

        satt_satt = lambda dt = 1e-5, graph = False: self.satt_launch(iniboost = iniboost, start = time, time = 3.55, launch_dir = direction, \
                                                                      start_pos_vel = array([position, velocity]), graph = graph, \
                                                                      dt = dt, iniboost_override = False)
        print "orbital position", orbit_pos(3.345)[:,target_plan]


        direction = orbit_pos(time)[:,target_plan] - position
        direction = direction/sqrt(pow(direction, 2).sum())
        direction_ = direction

        reverse_it = False
        for reverse in (False, True):
            if reverse:
                print "NOW REVERSING THE TURNING ANGLE"
            total = 0
            strike_max = 3
            past_angels = zeros(strike_max + 3) #yes, angels was intended, because it's late and I wanted to
            for i in range(inclines):
                angle = angle_ * d_ang**i

                if reverse:
                    angle = -angle

                direction = direction_

                closest_hit = inf
                strike = 0
                runs  = 0


                while strike < strike_max:
                    total +=1
                    runs  +=1
                    print "iteration nr. {0};    run nr. {1};    strikes: {2};    closest_hit: {3}".format(i+1, runs, strike, closest_hit)
                    satt_satt(dt = 1e-3, graph = False)

                    x, y = direction
                    direction = array([x*cos(angle) - y*sin(angle), x*sin(angle) + y*cos(angle)])
                    #print direction, position, velocity, norm(direction), "\n"
                    if  closest_hit > abs(self.min_dist[0]) and self.min_dist[0] > 0.0007: #we don't actually want to crash, right?
                        closest_hit = abs(self.min_dist[0])
                        time_of_hit = self.min_dist[1]
                        past_angels[1:] = past_angels[:-1]
                        past_angels[0]  = angle*runs
                        boost = direction * (iniboost/vec_abs(direction))

                    else:
                        strike += 1
                        if strike == strike_max and i+1 == inclines:
                            pass
                            satt_satt(dt = 1e-5, graph = False)
                            if total == (strike_max+1)*inclines:
                                reverse_it = True


                    if self.min_dist[0] < 0:
                        print "BRUH, CLOSEST_HIT GOES NEGATIVE FOR SOME REASON, FIX IT"
                angle_shift = past_angels[-1]

                x, y = direction_
                direction_ = array([x*cos(angle_shift) - y*sin(angle_shift), x*sin(angle_shift) + y*cos(angle_shift)])


            if not reverse_it:
                break

        print "the boost you need to input is ", boost
        print "the time the satellite reaches the planets is ", time_of_hit
        print "    {0},             closest_hit_final: {1}".format(self.min_dist[0], closest_hit)


    def new_satt_launch(self, time, position, velocity, simulate = True):
        Sattelite_Class.__init__(self)

        if simulate:
            self.satt_launch(iniboost = 0, start = 3.345, time = 3.5, launch_dir = array([0, 5]), \
                             start_pos_vel = array([position, velocity]), graph = False, \
                             dt = 1e-5, iniboost_override = False, progress = False)

        Sattelite_Class.__init__(self)
        infile = open("planet_positions.npy", "rb")
        orbit_pos_, orbit_vel_, time_ = load(infile)
        orbit_pos = interp1d(time_, orbit_pos_)
        orbit_vel = interp1d(time_, orbit_vel_)
        infile.close()

        print "planetary position at time", time, ": (", orbit_pos(time)[:,3][0], ",", orbit_pos(time)[:,3][1], ")"
        print "planetary velocity at time", time, ": (", orbit_vel(time)[:,3][0], ",", orbit_vel(time)[:,3][1], ")"

























class Planet_Class(Sattelite_Class):
    def __init__(self):
        Sattelite_Class.__init__(self)
        Sattelite_Class.find_solar(self)

        self.u_to_Kg = 1.6605402e-27
        self.k   = k = 1.38e-23

        self.doppler_reverse = lambda lam, lam_ref = 656.3, c = speed_of_light: (lam/float(lam_ref))*c
        self.doppler = lambda vel, lam_ref, c = float(speed_of_light): (vel/c)*lam_ref

        self.syst      = syst = AST2000SolarSystem(globalseed)
        self.FrefM     = Mp = syst.mass[3] * 1.989e30
        self.FrefR     = Rp = syst.radius[3] * 1e3 #convert from km to meter
        self.SAT_M     = Ms = 1100.
        self.SAT_A     = 15.
        self.rotation  = Om = 1./(syst.period[0]*24*60*60)
        self.temp      = T0 = 155.383775207
        self.G         = G  = 6.67408e-11
        self.g         = g  = lambda R = Rp: G*Mp/pow(R, 2)
        self.mass_mean = M  = (15.999 * 2 + 15.999 + 1.008 * 2 + 12.011 + 1.008 * 4 + 12.011 + 15.999)*self.u_to_Kg*0.25
        self.h_scaler  = h0 = (k*T0)/(g()*M)
        self.surf_dens = p0 = syst.rho0[3]
        self.dens_prof = p  = lambda h: p0*exp((Rp - h)/h0)
        self.surf_pres = P0 = p0*h0/g()
        self.pres_prof = P  = lambda h: P0*exp((Rp - h)/h0)
        self.F_drag_old= lambda R, A, v, C_d=1: 0.5*p(R)*A*pow(norm(v),3)*(v*(-1)) #does not account for wind
        def F_drag(R, A, v, C_d = 1):
            wind_abs = norm(R[:2])*tau*Om
            wind_direction = normalize(cross(R, (0,0,1)))
            wind = wind_abs*wind_direction
            vel  = v + wind
            return -0.5*p(norm(R))*A*pow(norm(vel),2)*normalize(vel)
        self.F_drag = F_drag

        self.Term_vel  = lambda R, A: sqrt((2*G*Mp*Ms)/(p(R)*pow(R,2)*A))
        self.Para_siz  = lambda R, V: (2*G*Mp*Ms)/(p(R)*pow(R*V, 2))

        self.gases_present = ["CO", "H2O", "CH4", "O2"]

        self.parachute = self.Para_siz(Rp, 3)



    def show_spectrum(self, file_loc = "spectrum_files_seed29"):
        data_file  = "\spectrum_seed29"
        numpy_file = file_loc + data_file + ".npy"
        if path.isfile(numpy_file):
            with open(numpy_file, "rb") as infile:
                self.spectrum = load(infile)
        else:
            self.spectrum = read_table(file_loc + data_file + ".txt", 2, read_every_nth_line = 1)

        if not path.isfile(numpy_file):
            with open(numpy_file, "wb") as outfile:
                save(outfile, self.spectrum)

        sigma_file = "\sigma_noise"
        numpy_file = file_loc + sigma_file + ".npy"
        if path.isfile(numpy_file):
            with open(numpy_file, "rb") as infile:
                self.sigma = load(infile)
        else:
            self.sigma  = read_table(file_loc + sigma_file + ".txt", 2, read_every_nth_line = 1)

        if not path.isfile(numpy_file):
            with open(numpy_file, "wb") as outfile:
                save(outfile, self.sigma)




    def model_spectrum(self, model_segment_length = "own", cen_acc = 300, width_acc = 30, Fmin_acc = 30, write_a_file = True, plot_a_graph = True):
        if model_segment_length == "own":
            model_segment_length = int(3.4e7/(cen_acc*width_acc*Fmin_acc))
            if virtual_memory()[0] > 1e10:
                print "hey, your PC is cool, so I'mma just take up \
                \nsome more space if that's cool wit you"
                model_segment_length*= 2 #activate if on computer with 16Gb RAM
        spectrum    = self.spectrum[1]
        given_sigma = self.sigma[1]
        lambdas     = self.spectrum[0]
        Fmin = linspace(0.7, 1, Fmin_acc)
        Fmax = 1
        centres = [630, 690, 760, 720, 820, 940, 1400, 1600, 1660, 2200, 2340, 2870]
        results = zeros((len(centres), 3))
        lines   = [["O2", 3], ["H2O", 3], ["CO2", 2], ["CH4", 2], ["CO", 1], ["N2O", 1]] #amount of spectral lines for each thing
        mass_base = []
        for l in lines:
            item = l[0]; rep = l[1]
            mass_base.extend([item]*rep)
        expected_temperatures = array([150, 450])

        u_to_Kg = self.u_to_Kg
        k       = self.k

        raise_it = lambda test, cen, sigma: exp(-0.5*pow((test-cen)/sigma, 2))
        model_f  = lambda Fm, test, cen, sigma: Fmax + (Fm - Fmax)*raise_it(test, cen, sigma)


        if write_a_file:
            outfile = open("spectrum_files_seed29\Gathered Data on Atmosphere.txt", "w")
            outfile.write("Model segment   : %i\n"%model_segment_length)
            outfile.write("Centre Accuracy : %i\n"%cen_acc)
            outfile.write("Width Accuracy  : %i\n"%width_acc)
            outfile.write("Fmin Accuracy   : %i\n\n"%Fmin_acc)
        for centre, result, item in zip(centres, results, mass_base):
            if  item == "O2":  #masses given in atmoic mass units [u]
                mass = 15.999 * 2
            elif item == "H2O":
                mass = 15.999 + 1.008  * 2
            elif item == "CO2":
                mass = 12.011 + 15.999 * 2
            elif item == "CH4":
                mass = 12.011 + 1.008  * 4
            elif item == "CO" :
                mass = 12.011 + 15.999
            elif item == "N2O":
                mass = 14.007 * 2 + 15.999
            else:
                raise Exception("This molecule is not on the list, what do we do?\n It goes to the block")

            m = mass*u_to_Kg
            maxvel_temp = sqrt((2*expected_temperatures*k)/m) #+ 1e4
            maxvel_tote = maxvel_temp + 1e4
            shift_temp  = self.doppler(maxvel_temp, centre)
            shift_tote  = self.doppler(maxvel_tote, centre)

            lambdacen   = linspace(centre - shift_tote[0], centre + shift_tote[1], cen_acc)
            lambdawidth = linspace(shift_temp[0], shift_temp[1], width_acc) #sigma
            lam_minimum = lambdacen.min()-lambdawidth.max()*2
            lam_maximum = lambdacen.max()+lambdawidth.max()*2
            testing_area = lambdas    [(lambdas > lam_minimum) * (lambdas < lam_maximum)]
            testing_real = spectrum   [(lambdas > lam_minimum) * (lambdas < lam_maximum)]
            divide_sigma = given_sigma[(lambdas > lam_minimum) * (lambdas < lam_maximum)]


            segments = range(model_segment_length, len(testing_area) + model_segment_length)[::model_segment_length]
            prev_seg = 0
            for curr_seg in segments:
                if  curr_seg > len(testing_area):
                    curr_seg = len(testing_area)
                curr_area = testing_area[prev_seg:curr_seg]
                curr_real = testing_real[prev_seg:curr_seg]
                curr_sigm = divide_sigma[prev_seg:curr_seg]
                cen, sigma, Fm, test = meshgrid(lambdacen, lambdawidth, Fmin, curr_area) #don't have to redo
                                                                                         #meshgrid everytime
                                                                                         #instead, make once
                                                                                         #and shift test by
                                                                                         #constant value
                #Instead of continually generating new meshgrids, can take same meshgrids and add some constant to test
                if  prev_seg == 0:
                    xhi = zeros(cen.shape[:-1])
                model = model_f(Fm, test, cen, sigma)
                xhi   = xhi + (((model - curr_real)/curr_sigm)**2).sum(axis = -1)

                prev_seg = curr_seg
            Fmin_val = Fm[:,:,:,0].flatten()[argmin(xhi)]
            Cent_val = cen[:,:,:,0].flatten()[argmin(xhi)]
            sigm_val = sigma[:,:,:,0].flatten()[argmin(xhi)]


            print "     ", centre,":", item
            print "         Fmin  :", Fmin_val
            print "         Centre:", Cent_val
            print "         Sigma :", sigm_val, "\n"
            if write_a_file:
                outfile.write("  %i : "%centre + item + "\n")
                outfile.write("      Fmin  :%f\n"%Fmin_val)
                outfile.write("      Centre:%f\n"%Cent_val)
                outfile.write("      Sigma :%f \n\n"%sigm_val)



            if plot_a_graph:

                plot(testing_area, testing_real, color="green")
#                plot(array([testing_area[0], testing_area[-1]]), array([testing_real.max(), testing_real.max()]), color ="blue")
#                plot(testing_area[(testing_area > testing_area[0]) & (testing_area < testing_area[0] + lambdawidth.max())],
#                                  full(len(testing_area[(testing_area > testing_area[0]) & (testing_area < testing_area[0] + lambdawidth.max())]), testing_real.min()))
#                plot(testing_area[(testing_area > testing_area[0]) & (testing_area < testing_area[0] + lambdawidth.min())],
#                                  full(len(testing_area[(testing_area > testing_area[0]) & (testing_area < testing_area[0] + lambdawidth.min())]), testing_real.min()+0.1))
                plot(testing_area, model_f(Fmin_val, testing_area, Cent_val, sigm_val), color="red")
                xlabel("wavelength[nm]")
                ylabel("flux[normalized around 1]")
                title("line %i"%centre)
                if write_a_file:
                    savefig("spectrum_files_seed29\line %i.png"%centre)
            show()
        if write_a_file:
            outfile.close()




    def land(self):
        self.syst.land_on_planet(3, "command files\landing instructions.txt")


    def simulate_orbit(self):
        r0 = array([-8367596.53536 ,  -103401558.597 ,  0])
        v0 = array([3714.7509967 ,  -307.302641468 ,  0])

        d = norm(r0)
        v0 = v0*0

        dt = 1e-2
        T  = 60*60*10

        t = linspace(0, T, T/dt)

        SAT_M = self.SAT_M
        SAT_A = self.SAT_A


        Mp = self.FrefM
        Rp = self.FrefR
        G  = self.G
        Fd = self.F_drag
        q  =-G*Mp

        r = zeros((int(T/dt),3))
        v = zeros((int(T/dt),3))
        dragg_ = zeros_like(t)
        r[0,:] = r0
        v[0,:] = v0 - r0 * G*Mp/pow(norm(r0), 3) * dt
        scatter(0, 0, color="green")
        for i in range(1, len(t)):
            j = i-1
            prog_bar(i, len(t))
            drag = Fd(r[j,:], SAT_A, v[j,:])/float(SAT_M)
#            print 0.5*self.pres_prof(norm(r[j,:]))*SAT_A*pow(norm( v[j,:]),3)*( v[j,:]*(-1))
#            print self.pres_prof(norm(r[j,:]))
#            print exp((norm(r[j,:] - Rp))/self.h_scaler), "exp"

            grav = (q/pow(float(norm(r[j,:])), 3))*r[j,:]
#            print drag
            a = grav + drag
#            print (norm(r[j,:]) - Rp)/self.h_scaler, self.h_scaler, "info", Rp, norm(r[j,:]), drag, j
            dragg_[i] = norm(drag)
            if i == 4:
                pass
#                break

#            print Fd(norm(r[j,:]), SAT_A, v[j,:])
            if norm(grav)/norm(drag) < 1000:
                v = v[:i,:]
                r = r[:i,:]
                time = t[j]
                print "We got too close, water loosened our feathers, we fell to earth"
                print time, norm(drag)
                v_orb = sqrt((G*Mp)/norm(r[j,:]))
#                boost_dir = rotate_angle_2d(normalize(r[j,:]), tau/4)
#                boost = v_orb * boost_dir
                print v_orb
                break
            v[i,:] = v[j,:] + a*dt
            r[i,:] = r[j,:] + v[i,:]*dt


#            scatter(r[i,0], r[i,1])
#            pause(0.01)
#            draw()

#            print r
#            print len_pos

#            print v


#        r_lis, v_lis = array(r_lis), array(v_lis)
#        orbit = axes(projection="3d")
##        orbit.plot(r_lis[:,0], r_lis[:,1], r_lis[:,2], label = "xyz", color = "red")

        xlim(-d*1.1, d*1.1)
        ylim(-d*1.1, d*1.1)
        plot(r[:,0], r[:,1])
        scatter(r[:,0], r[:,1])
        show()
        plot(dragg_)
        show()

#        print r_lis[:]
##        xlim(r_lis[:,0].min(), r_lis[:,0].max())
##        ylim(r_lis[:,1].min(), r_lis[:,1].max())
##        xlim(min(r_lis[:,0]), max(r_lis[:,0]))
##        ylim(min(r_lis[:,1]), max(r_lis[:,1]))
#
#        print min(r_lis[:][0]) == max(r_lis[:][0])
#


    def xyz_to_angular(self, xyz = "own"):
        if xyz == "own":
            print "xyz hardcoded because I felt like it"
            xyz = array([-894.49468733 , -11058.55176964 , 8689.28906549])
        x, y, z = xyz
        r = norm(xyz)
        theta = arccos(z/r)
        phi_obs = arcsin(y/(sin(theta)*r))
        print theta,
        time = 29500 #not really, look it up
        phi_diff = time*self.rotation
        phi = phi_obs - phi_diff
        print  phi, "\n"
























class Crash_Class(Planet_Class):
    def __init__(self):
        Planet_Class.__init__(self)
        self.LAN_M = 90
        self.LAN_A = 6.2




    def spiral(self):
        initials = array([[-1.54142266e+06,  -1.90564601e+07,   2.70864516e+03],\
                          [ 9.91833252e-02,   1.22619388e+00,   8.68928898e+03]]) #position, velocity, [m], [m/s]


        v = initials[1]
        r = initials[0]
        dt = 5e-2

        v = v*0

        A = self.parachute
        G = self.G
        M = self.FrefM
        R = self.FrefR
        m = self.LAN_M
        p = self.surf_dens

        r_lis = [norm(r)]
        v_lis = [norm(v)]
        t_lis = [0]
        d_lis = [norm(self.F_drag(r, A, v))]

        ROCKETPOWAH = 0


        ACTIVATE = False
        counter = 0
        t = 0
        while norm(r) - R > 10:
            t += dt
            counter += 1
            drag = self.F_drag(r, A, v)
            grav =-self.g()*normalize(r)

            a = grav + (drag + ROCKETPOWAH)/m
            v+= a*dt
            r+= v*dt

            if norm(r) - R <= 50 and not ACTIVATE:
                ROCKETPOWAH = 0.5*p*A*(self.Term_vel(R, A)**2 - 3**2)
#                ROCKETPOWAH = 0 #Override, not necessary to have rocket all the time
                ACTIVATE = True
                print "engines have activated, rockets at full capacity of", norm(ROCKETPOWAH), "newtons, giving a =", norm(ROCKETPOWAH)/m, "m/s"

            if norm(r) > r_lis[-1]:
                statement1 = "we're going upwards? t ="
                statement2 = "   height ="
                raise Exception(statement1, t, statement2, norm(r)-R)

            r_lis.append(norm(r))
            v_lis.append(norm(v))
            t_lis.append(t)
            d_lis.append(norm(drag))



            if norm(drag) > 25000:
                print "We're above recommended newtone dosage"
                if norm(drag) >= 250000:
                    print "bruh, lander fell apart at time = ", t
                    #raise Exception("bruh, lander fell apart at time = ", t)

            if counter%1000 == 0:
                continue
                print norm(r) - R


        print "final velocity", norm(v)
        t_lis = array(t_lis)
        r_lis = array(r_lis)
        d_lis = array(d_lis)

        plot(t_lis, r_lis)
        show()
        plot(t_lis, d_lis)
        show()
        plot(t_lis, v_lis)
        show()





if __name__ == "__main__":
    print """\n    THIS IS NOT THE FILE TO RUN
    IF YOU WANT TO RUN THE PROGRAM
    FOR ANY partICULAR PART, RUN
    THE PART NAMED AFTER THAT PART
    THANK YOU\n"""











    #This line is added just so I can say the code is 2500 lines long
