from __future__ import division
from sympy.interactive import printing
printing.init_printing(use_latex='mathjax')


from IPython.display import display,Image, Latex
import numpy as np
import math
import scipy.constants as sc

import sympy as sym
#from sympy import *

class FlatPlate(object):
    """ Definition of boundary layer thickness, friction coefficient, Nusselt number (both local and average)
        as a function of the regime.
        import HT_external_convection_flat_plate.py as flatplate
        
        bl =flatplate.FlatPlate(regime,thermal_bc,U_infty,nu,alpha_f,L,xi,Re_xc)
        where 
        regime = 'laminar' or 'turbulent' or 'mixed', 
        thermal_bc = 'isothermal', 'heat flux', 'unheated starting length',
        U_infty is the free stream velocity,
        nu the fluid viscosity,
        alpha the fluid thermal diffusivity,
        L length of the plate
        xi unheated started length
        Re_xc critical Reynolds number for transition laminar to turbulence
        
        """
    def __init__(self,regime,thermal_bc,U_infty,nu,alpha,L,xi,Re_xc):
        self.regime = regime
        self.thermal_bc = thermal_bc
        self.U_infty = U_infty
        self.nu = nu
        self.alpha = alpha
        self.Pr = self.nu/self.alpha
        self.L = L
        self.xi = xi
        self.Re_xc = Re_xc
        self.Re_L = self.L*self.U_infty/self.nu
        self.x_c = self.Re_xc*self.nu/self.U_infty
        if self.regime != "laminar" and self.regime and "turbulent" and self.regime != "mixed":
            print("regime is not properly defined")
        if self.thermal_bc != "isothermal" and self.thermal_bc != "heat flux" and self.thermal_bc != "unheated starting length":
            print("thermal boundary condition is not properly defined")
        if self.Re_L > self.Re_xc and self.regime == "laminar":
            print("The end plate Reynolds number is larger than the critical Reynolds number")
    def local(self,x):
        self.x = x
        self.Re_x = self.U_infty*self.x/self.nu
        if x == 0.:
            self.delta_x = 0.
            self.delta_Tx = 0.
            self.C_fx = 0.
            self.Nu_x = 0.
        else:
            if self.regime == "laminar":
                self.delta_x = 5.0*self.x/np.sqrt(self.Re_x)
                self.C_fx = 0.664*np.power(self.Re_x,-1./2.)
                if self.thermal_bc == "isothermal":
                    self.Nu_x = 0.332*np.power(self.Re_x,1./2.)*np.power(self.Pr,1./3.)
                elif self.thermal_bc == "heat flux":
                    self.Nu_x = 0.453*np.power(self.Re_x,1./2.)*np.power(self.Pr,1./3.)
                elif self.thermal_bc == "unheated starting length":
                    Re_xi = self.xi*self.U_infty/self.nu 
                    self.Nu_x = 0.332*np.power(self.Re_xi,1./2.)*np.power(self.Pr,1./3.)/\
                            np.power(1.-np.power(self.xi/self.x,3./4.),1./3.)
            elif self.regime == "turbulent":
                self.delta_x = 0.37*self.x*np.power(self.Re_x,-1./5.)
                self.C_fx = 0.0592*np.power(self.Re_x,-1./5.)
                if self.thermal_bc == "isothermal":
                    self.Nu_x = 0.0296*np.power(self.Re_x,4./5.)*np.power(self.Pr,1./3.)
                elif self.thermal_bc == "heat flux":
                    self.Nu_x = 0.0296*np.power(self.Re_x,4./5.)*np.power(self.Pr,1./3.)
                elif self.thermal_bc == "unheated starting length":
                    Re_xi = self.xi*self.U_infty/self.nu
                    self.Nu_x = 0.0296*np.power(self.Re_xi,4./5.)*np.power(self.Pr,1./3.)/\
                            np.power(1.-np.power(self.xi/self.x,9./10.),1./9.)
            elif self.regime == "mixed":
                if self.x < self.x_c:
                    self.delta_x = 5.0*self.x/np.sqrt(self.Re_x)
                    self.C_fx = 0.664*np.power(self.Re_x,-1./2.)
                    if self.thermal_bc == "isothermal":
                        self.Nu_x = 0.332*np.power(self.Re_x,1./2.)*np.power(self.Pr,1./3.)
                    elif self.thermal_bc == "heat flux":
                        self.Nu_x = 0.453*np.power(self.Re_x,1./2.)*np.power(self.Pr,1./3.)
                    elif self.thermal_bc == "unheated starting length":
                        Re_xi = self.xi*self.U_infty/self.nu 
                        self.Nu_x = 0.332*np.power(self.Re_xi,1./2.)*np.power(self.Pr,1./3.)/\
                            np.power(1.-np.power(self.xi/self.x,3./4.),1./3.)
                else:
                    self.delta_x = 0.37*self.x*np.power(self.Re_x,-1./5.)
                    self.C_fx = 0.0592*np.power(self.Re_x,-1./5.)
                    if self.thermal_bc == "isothermal":
                        self.Nu_x = 0.0296*np.power(self.Re_x,4./5.)*np.power(self.Pr,1./3.)
                    elif self.thermal_bc == "heat flux":
                        self.Nu_x = 0.0296*np.power(self.Re_x,4./5.)*np.power(self.Pr,1./3.)
                    elif self.thermal_bc == "unheated starting length":
                        Re_xi = self.xi*self.U_infty/self.nu
                        self.Nu_x = 0.0296*np.power(self.Re_xi,4./5.)*np.power(self.Pr,1./3.)/\
                            np.power(1.-np.power(self.xi/self.x,9./10.),1./9.)
            self.delta_Tx = self.delta_x*np.power(self.Pr,-1./3.)
    def average(self,x):
        self.x = x
        self.Re_x = self.U_infty*self.x/self.nu
        if x == 0.:
            print("The length cannot be zero")
        if self.regime == "laminar":
            self.C_fave = 1.328*np.power(self.Re_x,-1./2.)
            if self.thermal_bc == "isothermal" or self.thermal_bc == "heat flux":
                self.Nu_ave = 0.664*np.power(self.Re_x,1./2.)*np.power(self.Pr,1./3.)
            elif self.thermal_bc == "unheated starting length":
                p = 2.
                self.Nu_ave = 0.664*np.power(self.Re_xi,1./2.)*np.power(self.Pr,1./3.)*\
                              x/(x-xi)*np.power(1.-np.power(xi/x,(p+1.)/(p+2.)),p/(p+1.))
        elif self.regime == "turbulent":
            self.C_fave = 0.074*np.power(self.Re_x,-1./5.)
            if self.thermal_bc == "isothermal" or self.thermal_bc == "heat flux":
                self.Nu_ave = 0.037*np.power(self.Re_x,4./5.)*np.power(self.Pr,1./3.)
            elif self.thermal_bc == "unheated starting length":
                p = 8.
                self.Nu_ave = 0.664*np.power(self.Re_xi,1./2.)*np.power(self.Pr,1./3.)*\
                              x/(x-xi)*np.power(1.-np.power(xi/x,(p+1.)/(p+2.)),p/(p+1.))
        elif self.regime == "mixed":
            A = 0.037*np.power(self.Re_xc,4./5.)-0.664*np.power(self.Re_xc,1./2.)
            
            self.C_fave = 0.074*np.power(self.Re_x,-1./5.) - 2.*A/self.Re_x
            self.Nu_ave = (0.037*np.power(self.Re_x,4./5.) - A)*np.power(self.Pr,1./3.)
                
        
            