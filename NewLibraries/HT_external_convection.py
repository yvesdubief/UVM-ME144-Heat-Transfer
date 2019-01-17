""" Object name 1: FlatPlate
    Object name 2: CircularCylinder
    Object name 3: NoncircularCylinder
    Object name 4: BankofTubes
"""

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
        import HT_external_convection.py as extconv
        
        bl =extconv.FlatPlate(regime,thermal_bc,U_infty,nu,alpha,L,xi=0.0,Re_xc=5e5)
            where regime = 'laminar' or 'turbulent' or 'mixed', 
            thermal_bc = 'isothermal', 'heat flux', 'unheated starting length',
            U_infty is the free stream velocity,
            nu the fluid viscosity,
            alpha the fluid thermal diffusivity,
            L length of the plate
            xi unheated started length (only applies of using unheated starting length)
            Re_xc critical Reynolds number for transition laminar to turbulence
            
        output: bl.Re_L Reynolds at the trailing edge of the plate (x=L)
        
        bl.local(x) calculates the local Re (bl.Re_x), Cf (bl.Cf_x), Nu (bl.Nu_x) and velocity
                    thermal boundary layer thicknesses (bl.delta_x and bl.delta_Tx) at x based on thermal_bc
        
        bl.average(x) calculates the average Cf (bl.C_fave), Nu (bl.Nu_ave) over a length x from the leading edge
        
        """
    def __init__(self,regime,thermal_bc,U_infty,nu,alpha,L,xi=0.0,Re_xc=5e5):
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
            print("Warning: regime is not properly defined")
        if self.thermal_bc != "isothermal" and self.thermal_bc != "heat flux" and self.thermal_bc != "unheated starting length":
            print("Warning: thermal boundary condition is not properly defined")
        if self.Re_L > self.Re_xc and self.regime == "laminar":
            print("Warning: The end plate Reynolds number is larger than the critical Reynolds number, consider 'mixed' regime instead")
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
                    self.Re_xi = self.xi*self.U_infty/self.nu 
                    self.Nu_x = 0.332*np.power(self.Re_x,1./2.)*np.power(self.Pr,1./3.)/\
                            np.power(1.-np.power(self.xi/self.x,3./4.),1./3.)
            elif self.regime == "turbulent":
                self.delta_x = 0.37*self.x*np.power(self.Re_x,-1./5.)
                self.C_fx = 0.0592*np.power(self.Re_x,-1./5.)
                if self.thermal_bc == "isothermal":
                    self.Nu_x = 0.0296*np.power(self.Re_x,4./5.)*np.power(self.Pr,1./3.)
                elif self.thermal_bc == "heat flux":
                    self.Nu_x = 0.0296*np.power(self.Re_x,4./5.)*np.power(self.Pr,1./3.)
                elif self.thermal_bc == "unheated starting length":
                    self.Re_xi = self.xi*self.U_infty/self.nu
                    self.Nu_x = 0.0296*np.power(self.Re_x,4./5.)*np.power(self.Pr,1./3.)/\
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
                        self.Re_xi = self.xi*self.U_infty/self.nu 
                        self.Nu_x = 0.332*np.power(self.Re_x,1./2.)*np.power(self.Pr,1./3.)/\
                            np.power(1.-np.power(self.xi/self.x,3./4.),1./3.)
                else:
                    self.delta_x = 0.37*self.x*np.power(self.Re_x,-1./5.)
                    self.C_fx = 0.0592*np.power(self.Re_x,-1./5.)
                    if self.thermal_bc == "isothermal":
                        self.Nu_x = 0.0296*np.power(self.Re_x,4./5.)*np.power(self.Pr,1./3.)
                    elif self.thermal_bc == "heat flux":
                        self.Nu_x = 0.0296*np.power(self.Re_x,4./5.)*np.power(self.Pr,1./3.)
                    elif self.thermal_bc == "unheated starting length":
                        self.Re_xi = self.xi*self.U_infty/self.nu
                        self.Nu_x = 0.0296*np.power(self.Re_x,4./5.)*np.power(self.Pr,1./3.)/\
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
                self.Re_xi = self.xi*self.U_infty/self.nu
                self.Nu_ave = 0.664*np.power(self.Re_x,1./2.)*np.power(self.Pr,1./3.)*\
                              x/(x-self.xi)*np.power(1.-np.power(self.xi/x,(p+1.)/(p+2.)),p/(p+1.))
        elif self.regime == "turbulent":
            self.C_fave = 0.074*np.power(self.Re_x,-1./5.)
            if self.thermal_bc == "isothermal" or self.thermal_bc == "heat flux":
                self.Nu_ave = 0.037*np.power(self.Re_x,4./5.)*np.power(self.Pr,1./3.)
            elif self.thermal_bc == "unheated starting length":
                p = 8.
                self.Re_xi = self.xi*self.U_infty/self.nu
                self.Nu_ave = 0.664*np.power(self.Re_x,1./2.)*np.power(self.Pr,1./3.)*\
                              x/(x-self.xi)*np.power(1.-np.power(self.xi/x,(p+1.)/(p+2.)),p/(p+1.))
        elif self.regime == "mixed":
            A = 0.037*np.power(self.Re_xc,4./5.)-0.664*np.power(self.Re_xc,1./2.)
            
            self.C_fave = 0.074*np.power(self.Re_x,-1./5.) - 2.*A/self.Re_x
            self.Nu_ave = (0.037*np.power(self.Re_x,4./5.) - A)*np.power(self.Pr,1./3.)
                
class CircularCylinder(object):
    """ Nusselt correlations for cylinders
    import HT_external_convection.py as extconv
        
    bluff_body =extconv.CircularCylinder(correlation,Re,Pr,Pr_s = 0.0) 
        where Re, Pr, and Pr_s are the Reynolds number, Prandtl number of the flow and surface Prandtl numbers, respectively. If using Hilpert of Churchill Bernstein correlations,
        Re and Pr must be defined at film temperature, Pr_s can be set to anything since it is not used. 
        If using Zukauskas, Re and Pr are defined at temperature at infinity.
        correlation may be 'Hilpert', 'Churchill-Bernstein', 'Zukauskas'
        Example:
        bluff_body = extconv.CircularCylinder('Hilpert',Re,Pr)
        bluff_body = extconv.CircularCylinder('Churchill-Bernstein',Re,Pr)
        bluff_body = extconv.CircularCylinder('Zukauskas',Re,Pr,Pr_s = xx)
        
        Output: bluff_body.Nu average Nusselt number also bluff_body.correlation, bluff_body.Re, bluff_body.Pr, bluff_body.Pr_s
    
    bluff_body.correlation('Name of the correlation')
    Name of the correlation may be 'Hilpert', 'Churchill-Bernstein', 'Zukauskas'
    
    
    """
    def __init__(self,correlation,Re,Pr,Pr_s = 0.0):
        self.correlation = correlation
        self.Re = Re
        self.Pr = Pr
        self.Pr_s = Pr_s
        if correlation == "Zukauskas" and Pr_s == 0.0:
            print("Warning: Zukauskas correlation requires Pr_s")
        if self.correlation == "Hilpert":
            if self.Re < 0.4:
                print("Warning, Reynolds number too low for Hilpert Correlation")
                self.Nu = 0.
            elif  self.Re < 4.:
                C = 0.989
                m = 0.33
            elif  self.Re < 40:
                C = 0.911
                m = 0.385
            elif self.Re < 4000:
                C = 0.683
                m = 0.466
            elif self.Re < 40000.:
                C = 0.193
                m = 0.618
            elif self.Re <= 400000.:
                C = 0.027
                m = 0.805
            else :
                print("Warning Reynolds number is too high for the Hilpert Correlation")
                self.Nu = 0.
            if self.Re >= 0.4 and self.Re <= 400000.:
                self.Nu = C * self.Re**m * self.Pr**(1./3.)
        elif self.correlation == "Churchill-Bernstein":
            if (self.Re*self.Pr < 0.2):
                print("Warning: Product RePr lower than acceptable limit for Churchill Bernstein Correlation")
                self.Nu = 0.
    
            else:
                self.Nu = 0.3+(0.62*self.Re**(0.5)*self.Pr**(1./3.)) \
                  /(1.+(0.4/self.Pr)**(2./3.))**(1./4.) \
                *(1.+(self.Re/282000.)**(5./8.))**(4./5.)
        elif self.correlation == "Zukauskas":
            if (self.Pr <= 10):
                n = 0.37
            else:
                n = 0.36
            if (self.Re < 1.) and (self.Re > 1.e6):
                print("Warning Reynolds number out of bounds for the Zukauskas Correlation")
                self.Nu = 0.
            else:
                if (self.Re < 40.):
                    C = 0.75
                    m = 0.4
                elif (self.Re < 1000.):
                    C = 0.51
                    m = 0.5
                elif (self.Re < 2.e5):
                    C = 0.26
                    m = 0.6
                else:
                    C = 0.076
                    m = 0.7
                self.Nu = C*self.Re**m*self.Pr**n*(self.Pr/self.Pr_s)**(1./4.)

class NonCircularCylinder(object):
    """ Nusselt correlations for  cylinders with non circular cross-sections.
    import HT_external_convection.py as extconv
        
    bluff_body =extconv.NonCircularCylinder(geometry,Re,Pr) where 
    geometry = "angled square" square with stagnation point on one of its edges
               "square" square with stagnation point at the center of one of its faces
               "angled hexagon" hexagon with stagnation point on one of its edges
               "hexagon" hexagon with stagnation point at the center of one of its faces
               "thin plate" thin plate perpendicular to the flow
    Re: Reynolds number at film temperature
    Pr: Prandtl number at film temperature
    
    Output: bluff_body.Nu, bluff_body.Nu_front, bluff_body.Nu_back, the last two are for thin plate only
    also bluff_body.geometry, bluff_body.Re, bluff_body.Pr
    
    """
    def __init__(self,geometry,Re,Pr):
        self.geometry = geometry
        self.Re = Re
        self.Pr = Pr
        if self.geometry == "angled square":
            self.Nu_front = np.inf
            self.Nu_back = np.inf
            if self.Re < 6000:
                print("Warning, Reynolds number too low for Hilpert Correlation")
                self.Nu = np.inf
            elif  self.Re <= 60000.:
                C = 0.304
                m = 0.59
                self.Nu = C * self.Re**m * self.Pr**(1./3.)
            else :
                print("Warning Reynolds number is too high for the Hilpert Correlation")
                self.Nu = np.inf
        elif self.geometry == "square":
            self.Nu_front = np.inf
            self.Nu_back = np.inf
            if self.Re < 5000:
                print("Warning, Reynolds number too low for Hilpert Correlation")
                self.Nu = np.inf
            elif  self.Re <= 60000.:
                C = 0.158
                m = 0.66
                self.Nu = C * self.Re**m * self.Pr**(1./3.)
            else :
                print("Warning Reynolds number is too high for the Hilpert Correlation")
                self.Nu = np.inf
        elif self.geometry == "angled hexagon":
            self.Nu_front = np.inf
            self.Nu_back = np.inf
            if self.Re < 4500:
                print("Warning, Reynolds number too low for Hilpert Correlation")
                self.Nu = np.inf
            elif  self.Re <= 90700.:
                C = 0.150
                m = 0.638
                self.Nu = C * self.Re**m * self.Pr**(1./3.)
            else :
                print("Warning Reynolds number is too high for the Hilpert Correlation")
                self.Nu = np.inf
        elif self.geometry == "hexagon":
            self.Nu_front = np.inf
            self.Nu_back = np.inf
            if self.Re < 5200:
                print("Warning, Reynolds number too low for Hilpert Correlation")
                self.Nu = np.inf
            elif  self.Re <= 20400.:
                C = 0.164
                m = 0.638
                self.Nu = C * self.Re**m * self.Pr**(1./3.)
            elif  self.Re <= 105000.:
                C = 0.039
                m = 0.78
                self.Nu = C * self.Re**m * self.Pr**(1./3.)
            else :
                print("Warning Reynolds number is too high for the Hilpert Correlation")
                self.Nu = np.inf
        elif self.geometry == "thin plate":
            self.Nu = np.inf
            if self.Re < 10000:
                print("Warning, Reynolds number too low for Hilpert Correlation")
                self.Nu_front = np.inf
            elif  self.Re <= 50000.:
                C = 0.667
                m = 0.5
                self.Nu_front = C * self.Re**m * self.Pr**(1./3.)
            else :
                print("Warning Reynolds number is too high for the Hilpert Correlation for Nu_front")
                self.Nu_back = np.inf
            if self.Re < 7000:
                print("Warning, Reynolds number too low for Hilpert Correlation")
                self.Nu_back = np.inf
            elif  self.Re <= 80000.:
                C = 0.191
                m = 0.667
                self.Nu_back = C * self.Re**m * self.Pr**(1./3.)
            else :
                print("Warning Reynolds number is too high for the Hilpert Correlation for Nu_front")
                self.Nu_back = np.inf
                
class BankofTubes(object):
    """ Nusselt correlations for flow across banks of tubes
    import HT_external_convection.py as extconv
        
    bank =extconv.BankofTubes(arrangement,V_i,D,nu,Pr,Pr_s,S_L,S_T,N_L,N_T) where 
    arrangement = "aligned" tubes are aligned in row and column
                  "staggered" tubes are staggered from one row to the next
    V_i: Inlet velocity
    Pr: Prandtl number at arithmetic mean temperature
    Pr_s: Prandtl number at surface temperature
    S_L: tube center to tube center separation  between two consecutive rows (perpendicular to the flow)
    S_T: tube center to tube center separation  between two consecutive rows (aligned with the flow)
    N_L: number of rows perpendicular to flow
    N_T: number of rows aligned with flow if unknown giev your best guess. 
    
    Output: bank.Nu: average Nusselt number
            bank.arrangement,.Re,.Pr,.Pr_s,.S_L,.S_T,.N_L,.N_T,.N=self.N_L*self.N_T
    
    Functions:        
    bank.heat_rate(hbar,D,T_s,T_i,T_o) returns the heart per tube length based on the average convection coefficient
    bank.Vmax
    hbar = bank.Nu*k/D, the tube diameter and the log mean temperature obtained from the function
    Delta_T_lm(T_s,T_i,T_o) where T_s is the surface temperature, T_i is the inlet temperature, T_o the outlet
    The outlet temperature is calculated from the function
    
    T_o = extconv.temperature_outlet_tube_banks(T_s,T_i,D,N,N_T,hbar,rho_i,V_i,S_T,Cp_i)
    
    Note that the density and specific heat are estimated at T_i
    
    Other functions include:
    N_L = extconv.N_L_for_given_To(T_s,T_i,T_o,D,hbar,rho_i,V_i,S_T,Cp_i) which gives the Number of rows for a given T_o
    Note that the density and specific heat are estimated at T_i
    
    and
    
    T_lm = extconv.Delta_T_lm(T_s,T_i,T_o) which calculates the log mean
    
    self.Vmax
    """
    def __init__(self,arrangement,V_i,D,nu,Pr,Pr_s,S_L,S_T,N_L,N_T):
        self.arrangement = arrangement
        self.Pr = Pr
        self.Pr_s = Pr_s
        self.S_L = S_L
        self.S_T = S_T
        self.N_L = N_L
        self.N_T = N_T
        self.N = N_L*N_T
        self.D = D
        if self.arrangement == 'aligned':
            self.Vmax = self.S_T*V_i/(self.S_T-D)
        elif self.arrangement == 'staggered':
            self.S_D = np.sqrt(self.S_L**2+(self.S_T/2.)**2)
            self.Vmax = self.S_T*V_i/(2.*(self.S_D-D))
        Re = self.Vmax*self.D/nu
        self.Re = Re
        self.Nu = np.inf
        Corr_aligned = np.array([0.70,0.80,0.86,0.90,0.92,0.94,0.95,0.96,0.96,0.97,0.97,0.97,0.98,0.99,0.99,0.99,0.99,0.99,0.99])
        Corr_staggered = np.array([0.64,0.76,0.84,0.89,0.92,0.94,0.95,0.96,0.96,0.97,0.97,0.97,0.98,0.99,0.99,0.99,0.99,0.99,0.99])
        if (N_L < 20):
            if arrangement == 'aligned':
                Corr = Corr_aligned[N_L-1]
            elif arrangement == 'staggered':
                Corr = Corr_staggered[N_L-1]
        else:
            Corr = 1.
        if (Re < 10.):
            print('Warning: Re is out of bounds')
        if (Re >= 10.) and (Re <= 100.):
            if arrangement == 'aligned':
                C = 0.8
                m = 0.4
            elif arrangement == 'staggered':
                C = 0.9
                m = 0.4
            self.Nu = Corr*C*Re**m*Pr**(0.36)*(Pr/Pr_s)**(1./4.)
        elif (Re > 100.) and (Re <= 1000.):
            C = 0.51
            m = 0.
            self.Nu = Corr*C*Re**m*Pr**(0.36)*(Pr/Pr_s)**(1./4.)
        elif (Re > 1000.) and (Re <= 2.e5):
            if arrangement == 'aligned':
                if (S_T/S_L > 0.7):
                    C = 0.27
                    m = 0.63
                else:
                    print('Warning: inefficient, S_T/S_L<0.7')
                
            elif arrangement == 'staggered':
                if (S_T/S_L < 2):
                    C = 0.35*(S_T/S_L)**(1./5.)
                    m = 0.6
                else:
                    C = 0.40
                    m = 0.6
            self.Nu = Corr*C*Re**m*Pr**(0.36)*(Pr/Pr_s)**(1./4.)
        elif (Re > 2e5) and (Re <= 2.e6):
            if arrangement == 'aligned':
                C = 0.021
                m = 0.84
            elif arrangement == 'staggered':
                C = 0.022
                m = 0.84
            self.Nu = Corr*C*Re**m*Pr**(0.36)*(Pr/Pr_s)**(1./4.)
        else:
            print('Warning: Re is out of bounds')


    def heat_rate(self,hbar,D,T_s,T_i,T_o):
        DT_lm = Delta_T_lm(T_s,T_i,T_o)
        self.q=self.N*hbar*np.pi*D*DT_lm

            
def temperature_outlet_tube_banks(T_s,T_i,D,N,N_T,hbar,rho,V_i,S_T,Cp):
    return T_s-(T_s-T_i)*np.exp(-np.pi*D*N*hbar/(rho*V_i*N_T*S_T*Cp))

def N_L_for_given_To(T_s,T_i,T_o,D,hbar,rho,V_i,S_T,Cp):
    return -np.log((T_s-T_o)/(T_s-T_i))/(np.pi*D*hbar)*(rho*V_i*S_T*Cp)

def Delta_T_lm(T_s,T_i,T_o):
    return ((T_s-T_i)-(T_s-T_o))/np.log((T_s-T_i)/(T_s-T_o))

