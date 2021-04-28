""" 
Object name: 
    - HorizontalCylinder
    - FlatPlate
    - VerticalEnclosure
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
class FlatPlate(object):
    """ Natural convection caused by a flat plate at temperature T_s
        in a fluid at given ambient temperature T_infty. 
        Inputs:
        - Ra based on the |T_s - T_infty| and L=A_s/P_s where 
          A_s and P_s are the area and the perimeter of the plate,
          respectively.
        - Pr
        - surface is either 'upper' or 'lower'
        - surfaceT is either 'cold' or 'hot'
        Output:
        - Average Nu_L
        Limits:
        from NewLibraries import HT_natural_convection as natconv
        plate = natconv.FlatPlate(Ra= Ra_L, Pr = Pr, 
                                 surface = 'upper', surfaceT = 'hot')
        plate = natconv.FlatPlate(Ra= Ra_L, Pr = Pr, 
                                 surface = 'lower', surfaceT = 'cold')
        - Ra > 1e4
        - If 1e4 <= Ra <= 1e7 then Pr>= 0.7
        - If 1e7 <= Ra <= 1e11 then all Pr
        plate = natconv.FlatPlate(Ra= Ra_L, Pr = Pr, 
                                 surface = 'lower', surfaceT = 'hot')
        plate = natconv.FlatPlate(Ra= Ra_L, Pr = Pr, 
                                 surface = 'upper', surfaceT = 'cold')
        - 1e4<= Ra <= 1e9, Pr>=0.7
       
    """

    def __init__(self,Ra,Pr,surface='none',surfaceT='none'):
        self.Ra = Ra
        self.Pr = Pr
        self.surface = surface
        self.surfaceT = surfaceT
#         print(self.surface,self.surfaceT)
        if (self.surface != 'upper') and (self.surface != 'lower'):
            print(self.surface)
            print("you must specify surface='upper' or 'lower'")
        if (self.surfaceT != 'hot') and (self.surfaceT != 'cold'):
            print("you must specify surfaceT='hot' or 'cold'")
        if ((self.surface == 'upper') and (self.surfaceT == 'hot')) or \
            ((self.surface == 'lower') and (self.surfaceT == 'cold')):
                if (self.Ra >= 1e4) and (self.Ra <= 1e7):
                    if self.Pr < 0.7:
                        print("Warning: For %s surface of %s plate, \
                              the correlation is only valid for Pr>= 0.7 " %(self.surface,self.surfaceT))
                    self.Nu = 0.54*self.Ra**(1/4)
                elif (self.Ra >= 1e7) : #and (self.Ra <= 1e11):
                    self.Nu = 0.15*self.Ra**(1/3)
                else:
                    print("Ra is too small")
                    self.Nu = 0
        if ((self.surface == 'lower') and (self.surfaceT == 'hot')) or \
            ((self.surface == 'upper') and (self.surfaceT == 'cold')):
            if (self.Ra >= 1e4) : # and (self.Ra <= 1e9):
                self.Nu = 0.15*self.Ra**(1/5)
                if self.Pr < 0.7:
                    print("Warning: the correlation is only valid for Pr>= 0.7")
            else:
                print("Ra is too small")
                self.Nu=0


def Gr(g=9.81,beta=0.0,DT=0.0,D=0.0,nu=1.0):
    return (g*beta*DT*D**3)/(nu**2)

def Ra(g=9.81,beta=0.0,DT=0.0,D=0.0,nu=1.0,alpha=1.0):
    return (g*beta*DT*D**3)/(nu*alpha)



