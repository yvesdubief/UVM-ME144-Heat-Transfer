""" Object name: Fluid"""
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
        
class Fluid(object):
    """ How to:
        from NewLibraries import thermodynamics as thermo
        
        air_inlet = thermo.Fluid(material) material can be air, water, argon and krypton (see below for ranges)
        air_inlet.get_properties(T) to get thermodynamics of the fluid at temperature T in Kelvin.
        
        Compute thermodynamics properties of air between -150 C and 400 C, 
        water between 274K and 373K, argon between 100 and 700K and
        krypton between 150 and 700 K under 1 atm. Argon, krypton and water were obtained 
        through http://webbook.nist.gov/chemistry/fluid/
        More fluids to be added in the future
        air_inlet.beta thermal expansion coefficient
        air_inlet.rho density
        air_inlet.Cp specific heat
        air_inlet.mu dynamic viscosity
        air_inlet.k thermal conductivity
        air_inlet.nu kinematic viscosity
        air_inlet.alpha thermal diffusivity
        air_inlet.Pr
        Outputs:
        
        """
    def __init__(self,name):
        self.name = name
        
    def get_properties(self,T_o):
        self.T = T_o
        if self.name == 'water':
            if T_o < 274 or T_o > 373:
                print("Temperature is out of bounds for liquid water")
                return 
            Ttab,ptab,rhotab,Cptab,mutab,ktab = \
            np.genfromtxt('Tables/water1atm.csv', delimiter=',', skip_header = 1, unpack=True, dtype=float)
            Ntab = len(Ttab)
            Cptab *= 1e3
            nutab = mutab/rhotab 
            alphatab = ktab/(rhotab*Cptab)
            Prtab = nutab/alphatab
            dTtab = Ttab[1] - Ttab[0]
            # compute beta from -rho(d rho/dT)
            betatab = -(1./rhotab)*np.gradient(rhotab)/dTtab
            i = int((T_o-Ttab[0])/dTtab)
            if (i == Ntab - 1):
                i == Ntab - 2
        elif self.name == 'argon':
            if T_o < 100 or T_o > 700:
                print("Temperature is out of bounds for argon")
                return 
            Ttab,ptab,rhotab,Cptab,mutab,ktab = \
            np.loadtxt('Tables/Argon1atm.csv', delimiter=',', skiprows = 1, unpack=True, dtype=float)
            Ntab = len(Ttab)
            Cptab *= 1e3
            nutab = mutab/rhotab 
            alphatab = ktab/(rhotab*Cptab)
            Prtab = nutab/alphatab
            dTtab = Ttab[1] - Ttab[0]
            # compute beta from -rho(d rho/dT)
            betatab = -(1./rhotab)*np.gradient(rhotab)/dTtab
            i = int((T_o-Ttab[0])/dTtab)
            if (i == Ntab - 1):
                i == Ntab - 2
        elif self.name == 'krypton':
            if T_o < 150 or T_o > 740:
                print("Temperature is out of bounds for krypton")
                return 
            Ttab,ptab,rhotab,Cptab,mutab,ktab = \
            np.loadtxt('Tables/Krypton1atm.csv', delimiter=',', skiprows = 1, unpack=True, dtype=float)
            Ntab = len(Ttab)
            Cptab *= 1e3
            nutab = mutab/rhotab 
            alphatab = ktab/(rhotab*Cptab)
            Prtab = nutab/alphatab
            dTtab = Ttab[1] - Ttab[0]
            # compute beta from -rho(d rho/dT)
            betatab = -(1./rhotab)*np.gradient(rhotab)/dTtab
            i = int((T_o-Ttab[0])/dTtab)
            if (i == Ntab - 1):
                i == Ntab - 2
        elif self.name == 'air':
            if T_o < C2K(-150.) or T_o > C2K(400.):
                print("Temperature is out of bounds of the table for air")
                return
            Ttab,rhotab,Cptab,ktab,nutab,betatab,Prtab = \
            np.genfromtxt('Tables/air1atm.csv', delimiter=',', skip_header = 1, unpack=True, dtype=float)
            Ntab = len(Ttab)
            Ttab = C2K(Ttab)
            Cptab *= 1e3
            nutab *= 1e-6
            mutab = rhotab*nutab
            alphatab = ktab/(rhotab*Cptab)
            Prtab = nutab/alphatab
            i = 0
            while (Ttab[i] < T_o) and (i<Ntab):
                i += 1
            i -=1
            if (i == Ntab - 1):
                i = Ntab - 2
            
        else:
            print("warning, no table available for", self.name)
            return
        
        self.rho = interpolate_table(T_o,i,Ttab,rhotab)
        self.Cp = interpolate_table(T_o,i,Ttab,Cptab)
        self.mu = interpolate_table(T_o,i,Ttab,mutab)
        self.k = interpolate_table(T_o,i,Ttab,ktab)
        self.nu = interpolate_table(T_o,i,Ttab,nutab)
        self.alpha = interpolate_table(T_o,i,Ttab,alphatab)
        self.Pr = interpolate_table(T_o,i,Ttab,Prtab)
        if (self.name == 'air'):
            self.beta = 1./T_o
        else:
            self.beta = interpolate_table(T_o,i,Ttab,betatab)
        