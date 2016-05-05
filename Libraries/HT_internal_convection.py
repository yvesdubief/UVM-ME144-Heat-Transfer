""" 
Object name: PipeFlow
"""
import numpy as np
import scipy
import scipy.optimize

class PipeFlow(object):
    """ Definition of boundary layer thickness, friction coefficient, Nusselt number (both local and average)
        as a function of the regime.
        import HT_external_convection.py as extconv
        
        bl =extconv.FlatPlate(regime,thermal_bc,U_infty,nu,alpha_f,L,xi,Re_xc)
        where 
        regime = 'laminar' or 'turbulent' or 'mixed', 
        thermal_bc = 'isothermal', 'heat flux', 'unheated starting length',
        U_infty is the free stream velocity,
        nu the fluid viscosity,
        alpha the fluid thermal diffusivity,
        L length of the plate
        xi unheated started length
        Re_xc critical Reynolds number for transition laminar to turbulence
        
        bl.local(x) calculates the local Re, Cf, Nu at x based on thermal_bc
        
        bl.average(x) calculates the average Cf, Nu over a length x from the leading edge
        
        """
    def __init__(self, Re=0.0,eps = 0.0, Um = 0.0, D = 0.0, mdot = 0.0, nu = 0.0 ):
        if Re == 0.0:
            if Um != 0.0 and nu != 0.0:
                Re = Um*D/nu
            else:
                print("Warning if Re == 0, Um, D and nu must be specified")
        self.Re = Re
        self.eps = eps
        self.Um = Um
        self.D = D
        if D == 0.:
            print("Warning D must be different ")
        self.mdot = mdot
        self.nu = nu
        
        if eps == 0.0:
            print("Pipe wall is assumed to be hydrodynamically smooth")
            
    def pressure_drop_pipe(self,f = 0.0,L = 0.0,rho = 0.0, Um = 0.0):
        self.pressure_drop = f*(L/self.D)*(rho*Um**2)/2.


    def f_laminar(self, Re = 0.0):
        if Re == 0. and self.Re !=0:
            Re = self.Re
        else:
            print("Warning Reynolds number is not defined")
        self.f = 64./Re

    def f_turbulent(self,Re = 0.0, eps = 0.0):
        if Re == 0. and self.Re !=0:
            Re = self.Re
        else:
            print("Warning Reynolds number is not defined")
        e = eps
     
        f_0 = (0.790*np.log(Re)- 1.64)**(-2.)
        if (e > 0.):
            f_1 = 1./(-2.0*np.log10(e/3.71))**2
        else:
            f_1 = f_0
        f_guess = np.max(f_0,f_1)
        #f_guess = 0.04
        def f_tmp(x):
            y = (-2*np.log10((2.51/(Re*np.sqrt(x))) + (e/(3.71))) - 1.0/np.sqrt(x))
            return y
        y = scipy.optimize.fsolve(f_tmp, f_guess)
        self.f = y
        
    def h_laminar_isothermal(self,k):
        self.h = 3.66*k/self.D

    def h_laminar_isoflux(self,k):
        self.h = 4.36*k/D

    def Dittus_Boelter(self,mode,Pr,Re = 0.0):
        if Re == 0. and self.Re !=0:
            Re = self.Re
        else:
            print("Warning Reynolds number is not defined")
        if (mode == 'heating'):
            n = 0.4
        elif (mode == 'cooling'):
            n = 0.3
        self.Nu = 0.023*Re**(4./5.)*Pr**n

    def Nu_turbulent_Sieder_Tate(self,Pr,mu,mu_s, Re = 0.0):
        if Re == 0. and self.Re !=0:
            Re = self.Re
        else:
            print("Warning Reynolds number is not defined")
        self.Nu = 0.027*Re**(4./5.)*Pr*(1./3.)*(mu/mu_s)**0.14

    def Nu_turbulent_Gnielinski(self,Pr,f = 0.0,Re = 0.0):
        if Re == 0. and self.Re !=0:
            Re = self.Re
        else:
            print("Warning Reynolds number is not defined")
        self.Nu = (f/8.)*(Re-1000.)*Pr/(1+12.7*(f/8.)**0.5*(Pr**(2./3.)-1.))

    def Nu_turbulent_Skupinski(self,Pr, Re = 0.0):
        if Re == 0. and self.Re !=0:
            Re = self.Re
        else:
            print("Warning Reynolds number is not defined")
        self.Nu = 4.82+0.0185*(Re*Pr)**0.827

    def Nu_turbulent_Seban(self,Pr, Re = 0.0):
        if Re == 0. and self.Re !=0:
            Re = self.Re
        else:
            print("Warning Reynolds number is not defined")
        self.Nu = 5.0+0.025*(Re*Pr)**0.8

def log_mean_temperature(T_s,T_o,T_i):
    if (T_s < min(T_o,T_i)):
        DT_o = T_o-T_s
        DT_i = T_i-T_s
    elif (T_s > max(T_o,T_i)):
        DT_o = T_s-T_o
        DT_i = T_s-T_i
    return (DT_o-DT_i)/np.log(DT_o/DT_i)


def T_mx_Ts_constant(T_s,T_mi,P,mdot,Cp,hbar,x):
    return T_s-(T_s-T_mi)*np.exp(-P*x*hbar/(mdot*Cp))

def T_mo_T_infty(T_infty,T_mi,P,L,mdot,Cp,R_tot):
    return T_infty-(Tinfty-T_mi)*np.exp(-1/(mdot*Cp*Rtot))
