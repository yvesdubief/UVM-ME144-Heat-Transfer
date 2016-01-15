### definition of thermal resistance ###

import numpy as np
import math
import scipy.constants as sc

class Resistance(object):
    """ Defines thermal resistances for conduction, convection and radiation heat transfer. 
        First define the mode when defining your object
        Second use self.conduction, self.convection or self.radiation to calculate your resistance """
    def __init__(self,name):
        self.name = name
    def conduction(self,geo,k,r_a,r_b,A):
        self.geometry = geo
        if self.geometry == 'plane':
            self.R = r_a/(k*A)
        elif self.geometry == 'cylindrical':
            self.R = np.log(r_b/r_a)/(2.*math.pi*L*k)
        elif self.geometry == 'spherical':
            self.R = (1./r_a-1./r_b)/(4.*math.pi*k)
        else :
            print("geometry is not plane, cylindrical or spherical, cannot compute")
    def convection(self,h,A):
        self.R = 1./(h*A)
    def radiation(eps,T_s,T_infty,A):
        self.R = 1./(eps*sc.sigma*(T_s+T_infty)*(T_s**2+T_infty**2)*A)
        
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