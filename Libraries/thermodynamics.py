import numpy as np
import scipy
import scipy.optimize
from scipy.constants.constants import C2K
from scipy.constants.constants import K2C
from scipy.constants.constants import F2K
from scipy.constants.constants import K2F
import scipy.constants as sc

def interpolate_table(target,index,xquantity,yquantity):
    return yquantity[index] + \
                (yquantity[index+1]-yquantity[index])* \
                (target-xquantity[index])/(xquantity[index+1]-xquantity[index])

def water_look_up_1atm(T_o):
    T,p,rho,Cp,mu,k = \
    np.genfromtxt('tables/water1atm.csv', delimiter=',', skip_header = 1, unpack=True, dtype=None)
    Ntab = len(T)
    Cp *= 1e3
    nu = mu/rho 
    alpha = k/(rho*Cp)
    Pr = nu/alpha
    dT = T[1] - T[0]
    i = int((T_o-T[0])/dT)
    if (i == Ntab - 1):
        i == Ntab - 2
    rho_o = interpolate_table(T_o,i,T,rho)
    Cp_o = interpolate_table(T_o,i,T,Cp)
    mu_o = interpolate_table(T_o,i,T,mu)
    k_o = interpolate_table(T_o,i,T,k)
    nu_o = interpolate_table(T_o,i,T,nu)
    alpha_o = interpolate_table(T_o,i,T,alpha)
    Pr_o = interpolate_table(T_o,i,T,Pr)
    # compute beta from -rho(d rho/dT)
    beta = -(1./rho)*np.gradient(rho)/dT
    beta_o = interpolate_table(T_o,i,T,beta)
    return rho_o,Cp_o,mu_o,k_o,nu_o,alpha_o,Pr_o,beta_o

def air_look_up_1atm(T_o):
    T,rho,Cp,k,nu,beta,Pr = \
    np.genfromtxt('tables/air1atm.csv', delimiter=',', skip_header = 1, unpack=True, dtype=None)
    Ntab = len(T)
    T = C2K(T)
    Cp *= 1e3
    nu *= 1e-6
    mu = rho*nu
    alpha = k/(rho*Cp)
    Pr = nu/alpha
    i = 0
    while (T[i] < T_o) and (i<Ntab):
        i += 1
    i -=1
    if (i == Ntab - 1):
        i == Ntab - 2
    rho_o = interpolate_table(T_o,i,T,rho)
    Cp_o = interpolate_table(T_o,i,T,Cp)
    k_o = interpolate_table(T_o,i,T,k)
    nu_o = interpolate_table(T_o,i,T,nu)
    mu_o = interpolate_table(T_o,i,T,mu)
    alpha_o = interpolate_table(T_o,i,T,alpha)
    Pr_o = interpolate_table(T_o,i,T,Pr)
    beta = 1./T_o
    return rho_o,Cp_o,mu_o,k_o,nu_o,alpha_o,Pr_o,beta