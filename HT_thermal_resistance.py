### definition of thermal resistance ###

import numpy as np
import math
import scipy.constants as sc

def R_th_cond_plane(k,L,A):
    return L/(k*A)

def R_th_cond_cylindrical(k,r_a,r_b,L):
    return np.log(r_b/r_a)/(2.*math.pi*L*k)


def R_th_cond_spherical(k,r_a,r_b,L):
    return (1./r_a-1./r_b)/(4.*math.pi*k)

def R_th_convection(h,A):
    return 1./(h*A)

def R_th_radiation(eps,T_s,T_infty,A):
    return 1./(eps*sc.sigma*(T_s+T_infty)*(T_s**2+T_infty**2)*A)

### summation of thermal resistance (R is a vector) ###
def serial_sum(R):
    return np.sum(R)

def parallel_sum(R):
    return 1./(np.sum(1./R))