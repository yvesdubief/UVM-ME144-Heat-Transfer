""" Object name: Fluid"""
import sys  
sys.path.insert(0, 'Tables/') 
import numpy as np
import scipy
import scipy.optimize
from scipy.constants import convert_temperature
def C2K(T):
    return convert_temperature(T,'Celsius','Kelvin')
def C2F(T):
    return convert_temperature(T,'Celsius','Fahrenheit')
def F2K(T):
    return convert_temperature(T,'Fahrenheit','Kelvin')
def F2C(T):
    return convert_temperature(T,'Fahrenheit','Celsius')
def K2F(T):
    return convert_temperature(T,'Kelvin','Fahrenheit')
def K2C(T):
    return convert_temperature(T,'Kelvin','Celsius')
import scipy.constants as sc

def interpolate_table(target,index,xquantity,yquantity):
    return yquantity[index] + \
                (yquantity[index+1]-yquantity[index])* \
                (target-xquantity[index])/(xquantity[index+1]-xquantity[index])
        
class Fluid(object):
    """ How to:
        from NewLibraries import thermodynamics as thermo
        
        fluid_of_interest = thermo.Fluid(material,T) material can be air, water, argon and krypton (see below for ranges)
        and the temperature of the fluid T is in Kelvin.
        Outputs:
        The new object computes thermodynamic properties of air between -150 C and 400 C, 
        water between 274K and 373K, argon between 100 and 700K and
        krypton between 150 and 700 K under 1 atm. Argon, krypton and water were obtained 
        through http://webbook.nist.gov/chemistry/fluid/
        More fluids to be added in the future
        fluid_of_interest.beta thermal expansion coefficient
        fluid_of_interest.rho density
        fluid_of_interest.Cp specific heat
        fluid_of_interest.mu dynamic viscosity
        fluid_of_interest.k thermal conductivity
        fluid_of_interest.nu kinematic viscosity
        fluid_of_interest.alpha thermal diffusivity
        fluid_of_interest.Pr
        
        
        """
    def __init__(self,name,T,unit = "K",P = 101325.01):
        dirTables = 'https://raw.githubusercontent.com/yvesdubief/UVM-ME144-Heat-Transfer/master/Libraries/Tables/'
        self.name = name
        if unit == "C":
            T = C2K(T)
        elif unit == "F":
            T = F2K(T)
        self.T = T
        self.P = P
        if P != 101325.01:
            print("All available tables are for P=1ATM, reverting to P=101325.01Pa")
            self.P = 101325.01
        if self.name == 'water':
            if T < 274 or T > 373:
                print("Temperature is out of bounds for liquid water")
                return
            url = dirTables+'Water1atm.csv'
            Ttab,ptab,rhotab,Cptab,mutab,ktab = \
            np.genfromtxt(url, delimiter=',', skip_header = 1, unpack=True, dtype=float)
            Ntab = len(Ttab)
            Cptab *= 1e3
            nutab = mutab/rhotab 
            alphatab = ktab/(rhotab*Cptab)
            Prtab = nutab/alphatab
            dTtab = Ttab[1] - Ttab[0]
            # compute beta from -rho(d rho/dT)
            betatab = -(1./rhotab)*np.gradient(rhotab)/dTtab
            i = int((T-Ttab[0])/dTtab)
            if (i == Ntab - 1):
                i == Ntab - 2
        elif self.name == 'argon':
            if T < 100 or T > 700:
                print("Temperature is out of bounds for argon")
                return
            url = dirTables + 'Argon1atm.csv'
            Ttab,ptab,rhotab,Cptab,mutab,ktab = \
            np.loadtxt(url, delimiter=',', skiprows = 1, unpack=True, dtype=float)
            Ntab = len(Ttab)
            Cptab *= 1e3
            nutab = mutab/rhotab 
            alphatab = ktab/(rhotab*Cptab)
            Prtab = nutab/alphatab
            dTtab = Ttab[1] - Ttab[0]
            # compute beta from -rho(d rho/dT)
            betatab = -(1./rhotab)*np.gradient(rhotab)/dTtab
            i = int((T-Ttab[0])/dTtab)
            if (i == Ntab - 1):
                i == Ntab - 2
        elif self.name == 'krypton':
            if T < 150 or T > 740:
                print("Temperature is out of bounds for krypton")
                return
            url = dirTables + 'Krypton1atm.csv'
            Ttab,ptab,rhotab,Cptab,mutab,ktab = \
            np.loadtxt(url, delimiter=',', skiprows = 1, unpack=True, dtype=float)
            Ntab = len(Ttab)
            Cptab *= 1e3
            nutab = mutab/rhotab 
            alphatab = ktab/(rhotab*Cptab)
            Prtab = nutab/alphatab
            dTtab = Ttab[1] - Ttab[0]
            # compute beta from -rho(d rho/dT)
            betatab = -(1./rhotab)*np.gradient(rhotab)/dTtab
            i = int((T-Ttab[0])/dTtab)
            if (i == Ntab - 1):
                i == Ntab - 2
        elif self.name == 'air':
            if T < C2K(-150.) or T > C2K(400.):
                print("Temperature is out of bounds of the table for air")
                return
            url = dirTables + 'Air1atm.csv'
#             url = 'https://raw.githubusercontent.com/yvesdubief/UVM-ME144-Heat-Transfer/master/Libraries/Tables/Air1atm.csv'
#             url = 'Tables/Air1atm.csv'
            Ttab,rhotab,Cptab,ktab,nutab,betatab,Prtab = \
            np.genfromtxt(url, delimiter=',', skip_header = 1, unpack=True, dtype=float)
            Ntab = len(Ttab)
            Ttab = C2K(Ttab)
            Cptab *= 1e3
            nutab *= 1e-6
            mutab = rhotab*nutab
            alphatab = ktab/(rhotab*Cptab)
            Prtab = nutab/alphatab
            i = 0
            while (Ttab[i] < T) and (i<Ntab):
                i += 1
            i -=1
            if (i == Ntab - 1):
                i = Ntab - 2
            
        else:
            print("warning, no table available for", self.name)
            return
        
        self.rho = interpolate_table(T,i,Ttab,rhotab)
        self.Cp = interpolate_table(T,i,Ttab,Cptab)
        self.mu = interpolate_table(T,i,Ttab,mutab)
        self.k = interpolate_table(T,i,Ttab,ktab)
        self.nu = interpolate_table(T,i,Ttab,nutab)
        self.alpha = interpolate_table(T,i,Ttab,alphatab)
        self.Pr = interpolate_table(T,i,Ttab,Prtab)
        if (self.name == 'air'):
            self.beta = 1./T
        else:
            self.beta = interpolate_table(T,i,Ttab,betatab)
        
