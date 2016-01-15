import numpy as np
import scipy
import scipy.optimize

def Gr(g,beta,DT,L,nu):
    return (g*beta*DT*L**3)/(nu**2)

def Ra(g,beta,DT,L,nu,alpha):
    return (g*beta*DT*L**3)/(nu*alpha)

def Nu_vertical_enclosure(Ra,Pr,H,L):
    if (H/L) < 2.:
        if Ra*Pr/(0.2+Pr)> 1.e3:
            Nu = 0.18*(Pr/(0.2+Pr)*Ra)**0.29
        else:
            print('Ra is too low for this correlation')
            Nu = np.inf
    elif H/L < 10:
        if Ra < 1e10:
            Nu = 0.22*(Pr/(0.2+Pr)*Ra)**0.28*(H/L)**(-0.25)
        else:
            print('Ra is too high for this correlation')
            Nu = np.inf
    elif Ra < 1e4:
        print('Ra is too low for this correlation')
        Nu = np.inf
    elif Ra < 1e7:
        if Pr > 0.6 and Pr < 2e4:
            print('ok')
            Nu =0.42*Ra**0.25*Pr**0.012*(H/L)**(-0.3)
        else :
            print('Pr is out of bounds for this correlation')
            Nu = np.inf
    elif Ra < 1e9:
        if Pr > 0.6 and Pr < 20.:
            Nu =0.46*Ra**(1./3.)
        else :
            print('Pr is out of bounds for this correlation')
            Nu = np.inf
    else:
        print('Ra is too high, got nothing for you')
        Nu = np.inf
    return Nu
            