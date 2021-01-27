""" 
Object name: HorizontalCylinder
Functions: Gr(g,beta,DT,D,nu) gives the Grashoff number based on:
            gravity g, thermal expansion coefficient beta, Temperature difference DT, 
            length scale D, viscosity nu
           Ra(g,beta,DT,D,nu,alpha) gives the Rayleigh number where alpha is the thermal conductivity.
"""
import numpy as np
import scipy
import scipy.optimize

class HorizontalCylinder(object):
    """ Natural convection about a horizontal cylinder
        from NewLibraries import HT_natural_convection as natconv
        cyl = natconv.HorizontalCylinder(correlation, Ra, Pr = 0.0)
        where correlation is "Morgan" or "Churchill-Chu"
        cyl = natconv.HorizontalCylinder("Morgan", Ra)
        cyl = natconv.HorizontalCylinder("Churchill-Chu", Ra, Pr = xx)
    """

    def __init__(self,correlation="Morgan", Ra=0.0, Pr = 0.0):
        self.correlation = correlation
        self.Ra = Ra
        if correlation == "Morgan":
            if (Ra <= 1e-2):
                C=0.675
                n=0.058
            elif (Ra <= 1e2):
                C=1.02
                n=0.148
            elif (Ra <= 1e4):
                C=0.85
                n=0.188
            elif (Ra <= 1e7):
                C=0.480
                n=0.250
            elif (Ra <= 1e12):
                C=0.125
                n=0.333
            self.Nu = C*Ra**n
        elif correlation == "Churchill-Chu":
            if Pr == 0.:
                print("Warning you must specify Pr for Churchill and Chu correlation")
            else:
                self.Nu = (0.60+(0.387*Ra**(1./6.))/(1.+(0.559/Pr)**(9./16.))**(8./27.))**2
        else:
            print("Warning wrong correlation name")

class VerticalEnclosure(object):
    """ Natural convection about a horizontal cylinder
        from NewLibraries import HT_natural_convection as natconv
        cyl = natconv.HorizontalCylinder(correlation, Ra, Pr = 0.0)
        where correlation is "Morgan" or "Churchill-Chu"
        cyl = natconv.HorizontalCylinder("Morgan", Ra)
        cyl = natconv.HorizontalCylinder("Churchill-Chu", Ra, Pr = xx)
    """

    def __init__(self,Ra,Pr,H,L):
        self.Ra = Ra
        self.Pr = Pr
        self.H = H
        self.L = L
        if correlation == "Morgan":
            if (H/L) < 2.:
                if Ra*Pr/(0.2+Pr)> 1.e3:
                    self.Nu = 0.18*(Pr/(0.2+Pr)*Ra)**0.29
                else:
                    print('Ra is too low for this correlation')
                    self.Nu = np.inf
            elif H/L < 10:
                if Ra < 1e10:
                    self.Nu = 0.22*(Pr/(0.2+Pr)*Ra)**0.28*(H/L)**(-0.25)
                else:
                    print('Ra is too high for this correlation')
                    self.Nu = np.inf
            elif Ra < 1e4:
                print('Ra is too low for this correlation')
                self.Nu = np.inf
            elif Ra < 1e7:
                if Pr > 0.6 and Pr < 2e4:
                    print('ok')
                    self.Nu =0.42*Ra**0.25*Pr**0.012*(H/L)**(-0.3)
                else :
                    print('Pr is out of bounds for this correlation')
                    self.Nu = np.inf
            elif Ra < 1e9:
                if Pr > 0.6 and Pr < 20.:
                    self.Nu =0.46*Ra**(1./3.)
                else :
                    print('Pr is out of bounds for this correlation')
                    self.Nu = np.inf
            else:
                print('Ra is too high, got nothing for you')
                self.Nu = np.inf

def Gr(g=9.81,beta=0.0,DT=0.0,D=0.0,nu=1.0):
    return (g*beta*DT*D**3)/(nu**2)

def Ra(g=9.81,beta=0.0,DT=0.0,D=0.0,nu=1.0,alpha=1.0):
    return (g*beta*DT*D**3)/(nu*alpha)



