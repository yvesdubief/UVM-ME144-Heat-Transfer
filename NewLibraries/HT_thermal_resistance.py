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
        
        R.conduction(geo, k, thickness = 0.0, A = 1.0, L_pipe = 1.0, r_a = 0., r_b = 0.,k_name = "k",\
                   thickness_name = "L", L_pipe_name = "L", r_a_name = "r_a",r_b_name = "r_b",A_name = "A",\
                   T_a_name = "T_a",T_b_name = "T_b"),
        where geo can only be 'plane','cylindrical' or 'spherical'
        The minimum number of arguments are:
        R.conduction("plane", k, thickness = a) for heat flux (where a>0)
        R.conduction("plane", k, thickness = a, A = lengthorarea) for heat rate by unit length (m) or for heat rate (m^2)
        R.conduction("cylindrical",k, r_a = a, r_b = b) for heat rate per unit length of the pipe
        R.conduction("cylindrical",k, r_a = a, r_b = b, L_pipe = L) for heat rate 
        R.conduction("spherical",k, r_a = a, r_b = b) for heat rate
        
        thickness is the thickness of the material for plane conduction.
        r_a is the inner radius of the cylinder/sphere, r_b is the outer radius of the cylinder/sphere.
        A is the surface area of the system for plane conduction
        L_pipe is the pipe length for cylindrical conduction.
        All arguments ending with _name are used to write heat flux/rate equations(they are strings preferably 
        LaTeX formatted, without $$)
        
        R.convection(h,A,h_name = "h",A_name = "A",T_a_name = "T_a",Tb_name = "T_b"), where h is the convection coefficient (W/m^2K) and A is 
        the surface area. All arguments ending with _name are used to write the flux equations(they are strings 
        preferably LaTeX formatted)
        The minimum number of arguments are:
        R.convection(h,A)
        
        R.radiation(eps,T_s,T_sur,A,h_name = "h_r",A_name = "A",Ts_name = "T_s",Tsur_name = "T_{sur}"), where eps is the permissivity of the material, T_s
        the surface temperature, T_sur the far away surface temperature, A the surface area.
        The minimum number of arguments are:
        R.radiation(eps,T_s,T_sur,A)
        
        R.contact(R,A,R_name= "R_{t}",A_name = "A",T_a_name = "T_a",Tb_name = "T_b"), where R is the contact resistance, typically obtained from a table
        A is the surface area
        The minimum number of arguments are:
        R.contact(R,A)
        
        R.display_equation(index) displays the heat flux/rate equations for a given resistance. index is the number of 
        your resistance (you specify)
        
        Outputs include R[i].R the resistance of element i, R[i].h the convection or radiation coefficient.
        
        Functions include
        R_tot = res.serial_sum(R,first_resistance,last_resistance) sums serial resistance
        R_tot = res.parallel_sum(R,first_resistance,last_resistance) sums parallel resistance
        
        
        
        """
    def __init__(self,name,units):
        self.name = name
        self.units = units
    def conduction(self,geo, k, thickness = 0.0, A = 1.0, L_pipe = 1.0, r_a = 0., r_b = 0.,k_name = "k",\
                   thickness_name = "L", L_pipe_name = "L", r_a_name = "r_a",r_b_name = "r_b",A_name = "A",\
                   T_a_name = "T_a",T_b_name = "T_b"):
        self.geometry = geo
        self.k = k
        self.mode = 'conduction'
        self.k_name = k_name
        self.thickness = thickness
        self.L_pipe = L_pipe
        self.r_a = r_a
        self.r_b = r_b
        self.r_a_name = r_a_name
        self.r_b_name = r_b_name
        self.surface_name = A_name
        if geo == 'plane':
            if thickness == 0.:
                print("Warning you need to input thickness = a (where a>0) for plane conduction")
            self.surface_scale = A
            self.surface_name = A_name
        elif geo == 'cylindrical':
            self.surface_scale = L_pipe
            self.surface_name = L_pipe_name
        if geo != 'plane':
            if r_a == 0. or r_b == 0.:
                print("Warning you need to input r_a = a, r_b = b (where a,b>0) for cylindrical or spherical conduction")
        self.T_a_name = T_a_name
        self.T_b_name = T_b_name
        self.thickness_name = thickness_name
        self.L_pipe_name = L_pipe_name
        if self.geometry == 'plane':
            self.R = thickness/(k*A)
        elif self.geometry == 'cylindrical':
            print(self.r_a,self.r_b,self.L_pipe,self.k)
            if r_b == 0.:
                print("Warning rb must be specified for cylindrical geometries")
            self.R = np.log(r_b/r_a)/(2.*math.pi*L_pipe*k)
        elif self.geometry == 'spherical':
            if rb == 0.:
                print("Warning rb must be specified for spherical geometries")
            self.R = (1./r_a-1./r_b)/(4.*math.pi*k)
        else :
            print("geometry is not plane, cylindrical or spherical, cannot compute")
    def convection(self,h,A,h_name = "h",A_name = "A",T_a_name = "T_a",T_b_name = "T_b"):
        self.mode = 'convection'
        self.R = 1./(h*A)
        self.surface_scale = A
        self.h = h
        self.h_name = h_name
        self.surface_name = A_name
        self.T_a_name = T_a_name
        self.T_b_name = T_b_name
    def radiation(self,eps,T_s,T_sur,A,h_name = "h_r",A_name = "A",Ts_name = "T_s",Tsur_name = "T_{sur}"):
        self.R = 1./(eps*sc.sigma*(T_s+T_sur)*(T_s**2+T_sur**2)*A)
        self.mode = 'radiation'
        self.surface_scale = A
        self.h = eps*sc.sigma*(T_s+T_sur)*(T_s**2+T_sur**2)
        self.surface_name = A_name
        self.h_name = h_name
        self.Ts_name = Ts_name
        self.Tsur_name = Tsur_name
    def contact(self,R,A,R_name= "R_{t}",A_name = "A",T_a_name = "T_a",T_b_name = "T_b"):
        self.R = R/A
        self.mode = 'contact'
        self.R_name = R_name
        self.surface_scale = A
        self.surface_name = A_name
        self.T_a_name = T_a_name
        self.T_b_name = T_b_name
        
    def display_equation(self,index):

        Tasym = sym.symbols(self.T_a_name)
        Tbsym = sym.symbols(self.T_b_name)
        if self.units == 'W':
            Asym = sym.symbols(self.surface_name)
            namesym = "q_"+str(index)
        elif self.units == 'W/m':
            Asym = sym.symbols(self.surface_name)
            namesym = "q'_"+str(index)
        elif self.units == 'W/m^2':
            namesym = "q''_"+str(index)
        else:
            print('units are not properly defined')
        qsym = sym.symbols(namesym)
        Rsym = sym.symbols(self.name[1:-1])
        eq = sym.Eq(qsym,(1/Rsym)*(Tasym-Tbsym))
        if self.mode == 'conduction':
            thicksym = sym.symbols(self.thickness_name)
            rasym = sym.symbols(self.r_a_name)
            rbsym = sym.symbols(self.r_b_name)
            cstsym = sym.symbols(self.k_name)
            if self.geometry == 'plane':
                if self.units != 'W/m^2':
                    eq1 = sym.Eq(qsym,cstsym*Asym/thicksym*(Tasym-Tbsym))
                else:
                    eq1 = sym.Eq(qsym,cstsym/thicksym*(Tasym-Tbsym))
            elif self.geometry == 'cylindrical':
                if self.units == 'W':
                    eq1 = sym.Eq(qsym,2*sym.pi*cstsym/sym.log(rbsym/rasym)*Asym*(Tasym-Tbsym))
                else:
                    eq1 = sym.Eq(qsym,2*sym.pi*cstsym/sym.log(rbsym/rasym)*(Tasym-Tbsym))
            elif self.geometry == 'spherical':
                eq1 = sym.Eq(qsym,4*sym.pi*cstsym/(1/rasym-1/rbsym)*(Tasym-Tbsym))
                
        elif self.mode == 'convection':
            cstsym = sym.symbols(self.h_name)
            if self.units == 'W/m^2':
                eq1 = sym.Eq(qsym,cstsym*(Tasym-Tbsym))
            else:
                eq1 = sym.Eq(qsym,cstsym*Asym*(Tasym-Tbsym))
        elif self.mode == 'radiation':
            cstsym = sym.symbols(self.h_name)
            if self.units == 'W/m^2':
                eq1 = sym.Eq(qsym,cstsym*(Tasym-Tbsym))
            else:
                eq1 = sym.Eq(qsym,cstsym*Asym*(Tasym-Tbsym))
        elif self.mode == 'contact':
            cstsym = sym.symbols(self.R_name)
            if self.units == 'W/m^2':
                eq1 = sym.Eq(qsym,cstsym*(Tasym-Tbsym))
            else:
                eq1 = sym.Eq(qsym,(Asym/cstsym)*(Tasym-Tbsym))
        
        return display(eq,eq1)
        
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

