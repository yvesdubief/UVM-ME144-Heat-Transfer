"""Object: ExtSurfaces"""
from sympy.interactive import printing
printing.init_printing(use_latex='mathjax')


from IPython.display import display,Image, Latex
import numpy as np
import math
import scipy.constants as sc

import sympy as sym
#from sympy import *

class ExtSurfaces(object):
    """ Defines temperature distribution, heat rate for constant cross sectional area fins.
        from Libraries import HT_conduction_extended_surfaces as condext
        
        fin = condext.ExtSurfaces(T_b,T_infty,T_L,k,h,P,Ac,L)
            calculates fin.m, fin.M which are constants used in flux calculation. Also provides
            fin.theta_b,.theta_L,.T_b,.T_infty,.T_L,.h,.k,.h,.P,.Ac,.L,.Af(fin exposed surface area)
        fin.heat_rate(bc) calculate the heat rate for bc="convection", "adiabatic", "isothermal", "infinite"
            The ouptuts are fin.q_f, fin.effectiveness, fin.resistance, fin.efficiency
        fin.temperature(bc,x) calculates the temperature as a function of bc and the location x
            The output is fin.theta_over_theta_b
        fin.equations(T_b_name,T_infty_name,T_L_name,k_name,h_name,P_name,Ac_name,L_name) writes all the equations for you
            you need to run fin.heat_rate first.
    """
    def __init__(self,T_b,T_infty,T_L,k,h,P,Ac,L):
        self.T_b = T_b
        self.T_infty = T_infty
        self.T_L = T_L
        theta_b = T_b-T_infty
        theta_L = T_L-T_infty
        self.theta_b = T_b-T_infty
        self.theta_L = T_L-T_infty
        self.k = k
        self.h = h
        self.P = P
        self.Ac = Ac
        self.L = L
        self.Af = self.P*self.L
        m = np.sqrt(self.h*self.P/(self.k*self.Ac))
        self.m = m
        M = np.sqrt(self.h*self.P*self.k*self.Ac)*self.theta_b
        self.M = M
    def heat_rate(self,bc):
        self.bc = bc
        it_works = True
        if self.bc == "convection":
            self.q_f = self.M*(np.sinh(self.m*self.L) + (self.h/(self.m*self.k))*np.cosh(self.m*self.L))/\
                    (np.cosh(self.m*self.L) + (self.h/(self.m*self.k))*np.sinh(self.m*self.L))
        elif self.bc == "adiabatic":
            self.q_f = self.M*np.tanh(self.m*self.L)
        elif self.bc == "isothermal":
            self.q_f = self.M*np.cosh(self.m*self.L - self.theta_L/self.theta_b)/np.sinh(self.m*self.L)
        elif self.bc == 'infinite':
            self.q_f = self.M
        else:
            print("boundary condition is not properly defined")
            it_works = False
        if it_works:
            self.effectiveness = self.q_f/(self.h*self.Ac*self.theta_b)
            self.Resistance = self.theta_b/self.q_f
            self.efficiency = self.q_f/(self.h*self.Af*self.theta_b)
        
            
    def temperature(self,bc,x):
        self.bc = bc
        if self.bc == "convection":
            self.theta_over_theta_b = (np.cosh(self.m*(self.L-x)) + (self.h/(self.m*self.k))*np.sinh(self.m*(self.L-x)))/\
                    (np.cosh(self.m*self.L) + (self.h/(self.m*self.k))*np.sinh(self.m*self.L))
        elif self.bc == "adiabatic":
            self.theta_over_theta_b = np.cosh(self.m*(self.L-x))/np.cosh(self.m*self.L)
        elif self.bc == "isothermal":
            self.theta_over_theta_b = ((self.theta_L/self.theta_b)*np.sinh(self.m*self.L)+np.sinh(self.m*self.L - x))\
                                        /np.sinh(self.m*self.L)
        elif self.bc == 'infinite':
            self.theta_over_theta_b = np.exp(-self.m*x)
        else:
            print("boundary condition is not properly defined")
        self.T_x = self.T_infty + self.theta_over_theta_b*self.theta_b
        
    def equations(self,T_b_name,T_infty_name,T_L_name,k_name,h_name,P_name,Ac_name,L_name):
        T_x_sym = sym.symbols(r"T(x)")
        T_b_sym = sym.symbols(T_b_name)
        T_infty_sym = sym.symbols(T_infty_name)
        T_L_sym = sym.symbols(T_L_name)
        theta_sym = sym.symbols(r"\theta(x)")
        eq_theta = sym.Eq(theta_sym,T_x_sym-T_infty_sym)
        theta_b_sym = sym.symbols(r"\theta_b")
        eq_theta_b = sym.Eq(theta_b_sym,T_b_sym-T_infty_sym)
        theta_L_sym = sym.symbols(r"{\theta}_L")
        eq_theta_L = sym.Eq(theta_L_sym,T_L_sym-T_infty_sym)
        k_sym = sym.symbols(k_name)
        h_sym = sym.symbols(h_name)
        P_sym = sym.symbols(P_name)
        Ac_sym = sym.symbols(Ac_name)
        L_sym = sym.symbols(L_name)
        m_sym = sym.symbols(r"m")
        M_sym = sym.symbols(r"M")
        x_sym = sym.symbols(r"x")
        eq_m2 = sym.Eq(m_sym**2,h_sym*P_sym/(k_sym*Ac_sym))
        eq_M = sym.Eq(M_sym,sym.sqrt(h_sym*P_sym*k_sym*Ac_sym)*theta_b_sym)
        q_f_sym = sym.symbols(r"q_f")
        if self.bc == 'convection':
            eq_q = sym.Eq(q_f_sym,M_sym*(sym.sinh(m_sym*L_sym) + (h_sym/(m_sym*k_sym))*sym.cosh(m_sym*L_sym))/\
                    (sym.cosh(m_sym*L_sym) + (h_sym/(m_sym*k_sym))*sym.sinh(m_sym*L_sym)))
            eq_temp = sym.Eq(theta_sym/theta_b_sym,(sym.cosh(m_sym*(L_sym-x_sym)) + (h_sym/(m_sym*k_sym))*sym.sinh(m_sym*(L_sym-x_sym)))/\
                    (sym.cosh(m_sym*L_sym) + (h_sym/(m_sym*k_sym))*sym.sinh(m_sym*L_sym)))
        elif self.bc == "adiabatic":
            eq_q = sym.Eq(q_f_sym,M_sym*sym.tanh(m_sym*L_sym))
            eq_temp = sym.Eq(theta_sym/theta_b_sym,sym.cosh(m_sym*(L_sym-x_sym))/sym.cosh(m_sym*L_sym))
        elif self.bc == "isothermal":
            eq_q = sym.Eq(q_f_sym,M_sym*sym.cosh(m_sym*L_sym - theta_L_sym/theta_b_sym)/sym.sinh(m_sym*L_sym))
            eq_temp = sym.Eq(theta_sym/theta_b_sym,((theta_L_sym/theta_b_sym)*sym.sinh(m_sym*L_sym)+sym.sinh(m_sym*L_sym - x_sym))\
                                        /sym.sinh(m_sym*L_sym))
        elif self.bc == 'infinite':
            eq_q = sym.Eq(q_f_sym,M_sym)
            eq_temp = sym.Eq(theta_sym/theta_b_sym,sym.exp(-m_sym*x_sym))
        else:
            print("boundary condition is not properly defined")
        return display(eq_theta,eq_theta_b,eq_theta_L,eq_m2,eq_M,eq_q,eq_temp)
            
            
        