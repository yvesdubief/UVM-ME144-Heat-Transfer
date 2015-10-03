import numpy as np
import scipy
import scipy.optimize

def Gr(g,beta,DT,D,nu):
    return (g*beta*DT*D**3)/(nu**2)

def Ra(g,beta,DT,D,nu,alpha):
    return (g*beta*DT*D**3)/(nu*alpha)


def Nu_Morgan(Ra):
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
    return C*Ra**n

def Nu_Churchill_Chu(Ra,Pr):
    return (0.60+(0.387*Ra**(1./6.))/(1.+(0.559/Pr)**(9./16.))**(8./27.))**2 