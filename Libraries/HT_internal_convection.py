""" 
Object name: PipeFlow
"""
import numpy as np
import scipy
import scipy.optimize

class PipeFlow(object):
    """ Determination of Nu, pressure drop, mean temperature for internal convection
        import HT_internal_convection.py as intconv
        
        pipe =intconv.PipeFlow(D, Re=0.0, Um = 0.0, mdot = 0.0, nu = 0.0, rho = 0.0)
        where 
        D is the only required input and one of the following combination (Re, nu) or (Um, nu) or (mdot, rho, nu)
        Hence the minimum calls for PipeFlow are
        pipe =intconv.PipeFlow(D, Re= Re_m, nu = nu_m) outputs pipe.Um
        pipe =intconv.PipeFlow(D, Re= Re_m, nu = nu_m, rho = rho_m) outputs pipe.Um (bulk velocity) 
                and pipe.mdot (mass flow)
        pipe =intconv.PipeFlow(D, Um = 0.0, nu = 0.0) outputs pipe.Re
        pipe =intconv.PipeFlow(D, Um = Um, nu = nu_m, rho = rho_m) outputs pipe.Re, pipe.mdot
        pipe =intconv.PipeFlow(D, mdot = 0.0, nu = 0.0, rho = 0.0) outputs pipe.Re, pipe.Um
        
        pipe.f_laminar(Re) outputs the friction factor for laminar flow pipe.f
        pipe.f_turbulent(Re,eps = 0.0, nu = 0.0) outputs the friction factor for turbulent flow pipe.f
        
        The following correlations output pipe.Nu
        pipe.laminar_isothermal for isothermal wall boundary condition
        pipe.laminar_isoflux for isoflux wall boundary condition
        pipe.Dittus_Boelter(mode, Pr, Re = 0.) for turbulent flow where mode is either "heating" or "cooling"
            The Re is optional if omitted, the Reynolds number calculated in the object PipeFlow will be used
        pipe.Sieder_Tate(Pr,mu,mu_s, Re = 0.0) mu and mu_s are the mean and wall dynamics viscosities
            The Re is optional if omitted, the Reynolds number calculated in the object PipeFlow will be used
        pipe.Gnielinski( Pr, f,Re = 0.0): where f is the friction factor
            The Re is optional if omitted, the Reynolds number calculated in the object PipeFlow will be used
        
        """
    def __init__(self,D, Re=0.0, Um = 0.0 , mdot = 0.0, nu = 0.0, rho = 0.0, L = 1.0 ):
        self.D = D
        self.L = L
            
        if Re == 0.0:
            if Um != 0.0 and nu != 0.0:
                Re = Um*D/nu
            elif mdot != 0 and rho != 0.0 and nu != 0.0:
                Um = mdot/(rho*np.pi*D**2/4.)
                Re = Um*D/nu
            else:
                print("Warning if Re == 0, Um, D and nu or mdot, rho and nu must be specified")
                
        self.Re = Re
        if Um == 0.:
            if Re != 0. and nu != 0.:
                Um = Re*nu/D
                if mdot == 0.0 and rho != 0.0:
                    mdot = rho*Um*np.pi*D**2/4.
            elif mdot !=0.0 and rho != 0.0:
                Um = mdot/(rho*np.pi*D**2/4.)
                 
                
        self.Um = Um
        if mdot == 0.0:
            if rho != 0.0:
                mdot = rho*Um*np.pi*D**2/4.
            else:
                self.rho = 1.0
                self.mdot = rho*Um*np.pi*D**2/4.
        self.mdot = mdot
        self.nu = nu
        if Re == 0. and nu != 0.:
            Re = Um*D/nu
        self.Re = Re
        
        if rho == 0.0:
            self.rho = 1.0
            
        else:
            self.rho = rho
            
    


    def f_laminar(self, Re = 0.0):
        if Re == 0. and self.Re !=0:
            Re = self.Re
        elif Re == 0 and self.Re == 0.0:
            print("Warning Reynolds number is not defined")
        self.f = 64./Re
        self.dPdx = self.f*(self.L/self.D)*(self.rho*self.Um**2)/2.

    def f_turbulent(self,Re = 0.0, eps = 0.0):
        if Re == 0. and self.Re !=0.0:
            Re = self.Re
        elif Re == 0 and self.Re == 0.0:
            print("Warning Reynolds number is not defined")
        if eps == 0.0:
            print("Pipe wall is assumed to be hydrodynamically smooth") 
        e = eps
     
        f_0 = (0.790*np.log(Re)- 1.64)**(-2.)
        if (e > 0.):
            f_1 = 1./(-2.0*np.log10(e/3.71))**2
        else:
            f_1 = f_0
        f_guess = np.max([f_0,f_1])
        #f_guess = 0.04
        def f_tmp(x):
            y = (-2*np.log10((2.51/(Re*np.sqrt(x))) + (e/(3.71))) - 1.0/np.sqrt(x))
            return y
        y = scipy.optimize.fsolve(f_tmp, f_guess)
        self.f = y[0]
        self.dPdx = self.f*(self.L/self.D)*(self.rho*self.Um**2)/2.
        
    def laminar_isothermal(self):
        self.Nu = 3.66

    def laminar_isoflux(self):
        self.Nu = 4.36

    def Dittus_Boelter(self,mode,Pr,Re = 0.0):
        if Re == 0. and self.Re !=0:
            Re = self.Re
        else:
            print("Warning Reynolds number is not defined")
        if (mode == 'heating'):
            n = 0.4
        elif (mode == 'cooling'):
            n = 0.3
        else:
            print("Warning you have to specify mode='heating' or 'cooling'")
        self.Nu = 0.023*Re**(4./5.)*Pr**n

    def Sieder_Tate(self,Pr,mu,mu_s, Re = 0.0):
        if Re == 0. and self.Re !=0:
            Re = self.Re
        else:
            print("Warning Reynolds number is not defined")
        self.Nu = 0.027*Re**(4/5)*Pr**(1/3)*(mu/mu_s)**0.14

    def Gnielinski(self, Pr, f,Re = 0.0):
        if Re == 0. and self.Re !=0:
            Re = self.Re
        else:
            print("Warning Reynolds number is not defined")
        self.Nu = (f/8.)*(Re-1000.)*Pr/(1+12.7*(f/8.)**0.5*(Pr**(2./3.)-1.))

    def Skupinski(self,Pr, Re = 0.0):
        if Re == 0. and self.Re !=0:
            Re = self.Re
        else:
            print("Warning Reynolds number is not defined")
        self.Nu = 4.82+0.0185*(Re*Pr)**0.827

    def Seban(self,Pr, Re = 0.0):
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

def T_mo_T_infty(T_infty,T_mi,mdot,Cp,R_tot):
    return T_infty-(Tinfty-T_mi)*np.exp(-1/(mdot*Cp*Rtot))

def L_given_other_params(T_infty,T_mo,T_mi,mdot,Cp,Rptot):
    return -mdot*Cp*Rptot*np.log((T_infty -T_mo)/(T_infty - T_mi))

