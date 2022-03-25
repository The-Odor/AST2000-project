# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 23:06:46 2017

@author: theod
"""
# mostly just calculations
from numpy import pi
from numpy import sqrt
T = 4097.8720398
R = 401034.446737 * 1000
s = 5.6704e-8
lumi = 4*pi*R**2*s*T**4

lumi_rel = lumi/3.828e26

sun_lumi = 3.828e26

sun_mass = 1.989e30
sun_temp = 5778
sun_radi = 695700 * 1000
sol_mass = 0.45144
sol_temp = T/sun_temp



lumi_from_mass = sol_mass**4
lumi_from_temp = sol_temp**2

print lumi_rel
print lumi_from_mass
print lumi_from_temp
k = 1.38064852e-23
G = 6.67408e-11
Jeans_len = lambda M, T, mu, k_ = k, G_ = G: (G_*mu*M)/(5*k_*T)

mu = 1.6737236e-27*0.75 + 6.6464764e-27*0.25

Jeans = Jeans_len(sun_mass, 10, mu)
lumi_Jeans = 4*pi*Jeans**2*s*10**2
print lumi_Jeans/sun_lumi



#
#def Formal_jeans():
#    T = 10
#    k = 1.38064852e-23
#    G = 6.67408e-11
#    mu = 1.6737236e-27*0.75 + 6.6464764e-27*0.25
#    rho = (3*sun_mass)/(3*pi*)
#    formal_jeans = (15*k*T)/(4*pi*G*mu*rho)

sol_mass = 0.451449687213
rho_sol = (3*sol_mass*sun_mass)/(4*pi*R**3)
mu_prot = 1.67e-27
T_core = T + ((2.*mu_prot)/(3*k))*G*pi*rho_sol*R**2

print T_core

m_prot = mu_prot
m_heli = 6.6464764e-27
m_CNO  = 1.9944e-26 + 2.32586e-26 + 2.6566e-26

m_sum = m_prot + m_heli + m_CNO
X_H = m_prot/m_sum
X_He = m_heli/m_sum
X_CNO = m_CNO/m_sum
print
print X_H
print X_He
print X_CNO
print 

T6 = 9.; T6 = T_core / 1e6
e0pp = 1.08e-12
e0CNO = 8.24e-31
epp = e0pp*(X_H**2)*rho_sol*(T6**4)
eCNO = e0CNO*X_H*X_CNO*rho_sol*(T6**20)
print epp, eCNO

M_core = rho_sol*(4./3)*pi*pow(0.2*R, 3)

RADILUMI = M_core*(epp + eCNO)

print RADILUMI

radiation_temperature = sqrt(sqrt((RADILUMI)/(s*4*pi*R**2)))
print radiation_temperature, R

sol_lifetime = 1/pow(sol_mass, 3)
print sol_lifetime

def white_dwarf():
    import numpy as np
    
    pi = np.pi
    h  = 6.626e-34#planck constant
    ZA = 0.5
#    G  = G
    mH = mu_prot
    me = 9.109e-31#electron mass
    Mc = 1.44*sun_mass#chandrasekhar
    M  = Mc*(sol_mass/8.)

    factors = np.zeros(4)
    factors[0] = pow(3/(2*pi), 4/3.)
    factors[1] = pow(h,2)/(20*me*G)
    factors[2] = pow(ZA/mH, 5/3.)
    factors[3] = pow(M, -1/3.)
    
    RWD = factors.prod()
    
    print "litre", (M/((4/3.)*pi*pow(RWD,3)))*1000
    
    print "acceleration", G*M/pow(RWD, 2)
    
    
    
    return RWD
    
print white_dwarf()