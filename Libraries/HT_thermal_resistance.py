### definition of thermal resistance ###
from __future__ import division
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
        .conduction(self,geo,k,r_a,r_b,A,r_a_name,r_b_name,A_name,Ta_name,Tb_name), where geo can only be 'plane',
        'cylindrical' or 'spherical', r_a is the first length and the only that matters in case
        of planar conduction, r_b is the second length. r_b must be an input even for planar
        (could be 0.). A is the surface area of the system for plane conduction, or the length of the cylinder
        for cylindrical conduction. Set to 0 if spherical. All arguments ending with _name are
        used to write the flux equations(they are strings preferably LaTeX formatted)
        """
    def __init__(self,name,units):
        self.name = name
        self.units = units
    def conduction(self,geo,k,ra,rb,A,k_name,ra_name,rb_name,A_name,Ta_name,Tb_name):
        self.geometry = geo
        self.mode = 'conduction'
        self.k_name = k_name
        self.ra_name = ra_name
        self.rb_name = rb_name
        self.surface_name = A_name
        self.surface_scale = A
        self.Ta_name = Ta_name
        self.Tb_name = Tb_name
        if self.geometry == 'plane':
            self.R = ra/(k*A)
        elif self.geometry == 'cylindrical':
            self.R = np.log(rb/ra)/(2.*math.pi*self.surface_scale*k)
        elif self.geometry == 'spherical':
            self.R = (1./ra-1./rb)/(4.*math.pi*k)
        else :
            print("geometry is not plane, cylindrical or spherical, cannot compute")
    def convection(self,h,A,h_name,A_name,Ta_name,Tb_name):
        self.mode = 'convection'
        self.R = 1./(h*A)
        self.surface_scale = A
        self.h = h
        self.h_name = h_name
        self.surface_name = A_name
        self.Ta_name = Ta_name
        self.Tb_name = Tb_name
    def radiation(self,eps,T_s,T_sur,A,h_name,A_name,Ta_name,Tb_name):
        self.R = 1./(eps*sc.sigma*(T_s+T_sur)*(T_s**2+T_sur**2)*A)
        self.mode = 'radiation'
        self.surface_scale = A
        self.h = eps*sc.sigma*(T_s+T_sur)*(T_s**2+T_sur**2)
        self.surface_name = A_name
        self.h_name = h_name
        self.Ta_name = Ta_name
        self.Tb_name = Tb_name
    def contact(self,R,A,R_name,A_name,Ta_name,Tb_name):
        self.R = R/A
        self.mode = 'contact'
        self.R_name = R_name
        self.surface_scale = A
        self.surface_name = A_name
        self.Ta_name = Ta_name
        self.Tb_name = Tb_name
        
    def display_equation(self,index):

        Tasym = sym.symbols(self.Ta_name)
        Tbsym = sym.symbols(self.Tb_name)
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
            rasym = sym.symbols(self.ra_name)
            rbsym = sym.symbols(self.rb_name)
            cstsym = sym.symbols(self.k_name)
            if self.geometry == 'plane':
                if self.units != 'W/m^2':
                    eq1 = sym.Eq(qsym,cstsym*Asym/rasym*(Tasym-Tbsym))
                else:
                    eq1 = sym.Eq(qsym,cstsym/rasym*(Tasym-Tbsym))
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

