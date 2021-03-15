"""Object name: Resistance
   Function name: serial_sum(R,nori,nend), performs serial sum of a resistance object list from nori to nend
   Function name: parallel_sum(R,nori,nend), performs parallel sum of a resistance object list from nori to nend
   """
### definition of thermal resistance ###
from sympy.interactive import printing
printing.init_printing(use_latex='mathjax')


from IPython.display import display,Image, Latex
import numpy as np
import math
import scipy.constants as sc

import sympy as sym
#from sympy import *

class Resistance(object):
    """ Defines thermal resistances for conduction, convection and radiation heat transfer. 
        First define the object attached with class with the name used in the thermal circuit
        and the units, which can only be 'W', 'W/m' or 'W/m^2'
        Second use self.conduction, self.convection or self.radiation to calculate your 
        resistance. Each mode requires different arguments:
        
        from Libraries import HT_thermal_resistance as res
        R = []
        R.append(res.Resistance("$label$", "units")) where units = 'W', 'W/m' or 'W/m^2'
        
        then
        
        For conduction, there are 3 options:
        
        - R.cond_plane(k, L, A = 1.0) for planar conduction: k is the thermal conductivity,
                L is the thickness of the wall, and A is the optional surface area (=1 by default)
        - R.cond_cylinder(k , ra, rb, L = 1.0, angle = 2.*math.pi) for conduction in a 
                cylindrical shell between the radii ra (internal) and rb (external). L is the length
                of the shell (optional, default = 1) and angle is angular dimension of shell, also 
                optional and set to a full revolution by default (2 pi)
        - R.cond_sphere(k, ra, rb, scale = 1.0) for conductuion within a spherical shell bounded by radii ra and rb
            ra < rb. The optional parameter scale allows to calculate the thermal resistance for a fraction
            of a spherical shell. For instance a cornea is about 1/3 of spherical shell, so scale = 1./3.
        
        Convection:
        - R.convection(h, A = 1.0), where h is the convection coefficient (W/m^2K) and A is 
        the surface area (optional, default is unit surface aera 1 m^2)
        
        Radiation:
        - R.radiation(eps, T_s, T_sur, A = 1.0), where eps is the permissivity of the material, T_s
        the surface temperature, T_sur the far away surface temperature, A the surface area (optional, 
        by default A is the unit surface area 1 m^2).
        
        Contact:
        
        - R.contact(R,A,R_name= "R_{t}",A_name = "A",T_a_name = "T_a",Tb_name = "T_b"), where R is the contact resistance, typically obtained from a table
        A is the surface area
        The minimum number of arguments are:
        R.contact(R,A)
        
        R.display_equation(index) displays the heat flux/rate equations for a given resistance. index is the number of 
        your resistance (you specify)
        
        Outputs:
        - R[i].R the resistance of element i, R[i].h the convection or radiation coefficient.
        
        Functions include
        R_tot = res.serial_sum(R,first_resistance,last_resistance) sums serial resistance
        R_tot = res.parallel_sum(R,first_resistance,last_resistance) sums parallel resistance
        
        
        
        """
    def __init__(self,name,units):
        self.name = name
        self.units = units
    def cond_plane(self, k, L, A = 1.0):
        self.mode = "conduction"
        self.geometry = "planar"
        self.k = k
        if k <= 0.:
            print("problem with the definition of thermal conductivity")
        self.L = L
        self.A = A
        self.R = self.L / (self.k * self.A)
    def cond_cylinder(self, k , ra, rb, L = 1.0, angle = 2.*math.pi):
        self.mode = "conduction"
        self.geometry = "cylindrical"
        self.k = k
        if k <= 0.:
            print("problem with the definition of thermal conductivity")
        self.ra = ra
        self.rb = rb
        if ra*rb <= 0.:
            print("problem with the definition of radii")
        self.L = L
        self.angle = angle
        self.R = np.log(rb/ra)/(angle*L*k)
    def cond_sphere(self, k, ra, rb, scale = 1.0):
        self.mode = "conduction"
        self.geometry = "spherical"
        self.k = k
        if k <= 0.:
            print("problem with the definition of thermal conductivity")
        self.ra = ra
        self.rb = rb
        if ra*rb <= 0.:
            print("problem with the definition of radii")
        self.R = (1./r_a-1./r_b)/(scale*4.*math.pi*k)   
    def convection(self, h, A = 1.0):
        self.mode = 'convection'
        self.geometry = "whatever"
        self.R = 1./(h*A)
        self.A = A
        self.h = h
    def radiation(self,eps,T_s,T_sur, A = 1.0):
        self.R = 1./(eps*sc.sigma*(T_s+T_sur)*(T_s**2+T_sur**2)*A)
        self.mode = 'radiation'
        self.geometry = "whatever"
        self.A = A
        self.h = eps*sc.sigma*(T_s+T_sur)*(T_s**2+T_sur**2)
    def contact(self, R, A=1.0):
        self.R = R/A
        self.geometry = 'whatever'
        self.mode = 'contact'
        
    
        
### summation of thermal resistance (R is a vector) ###
def serial_sum(R,nori,nend):
    sum = 0.
    for i in range(nori,nend+1):
        sum += R[i].R
    return sum

def parallel_sum(R,nori,nend):
    sum = 0.
    for i in range(nori,nend+1):
        sum += 1./R[i].R
    return 1./sum


