import numpy as np
import scipy
import scipy.optimize

def linear_interpolation(x_t,x_1,x_2,y_1,y_2):
    return y_1+(y_2-y_1)*(x_t-x_1)/(x_2-x_1)

def pressure_drop_pipe(f,L,D,rho,u_m):
    return f*(L/D)*(rho*u_m**2)/2.

def f_pipe_laminar(Re_D):
    return 64./Re_D

def f_pipe_colebrook(Re_D,eps):
    Re = Re_D
    e = eps
     
    f_0 = (0.790*np.log(Re)- 1.64)**(-2.)
    if (e > 0.):
        f_1 = 1./(-2.0*np.log10(e/3.71))**2
    else:
        f_1 = f_0
    f_guess = np.max(f_0,f_1)
    #f_guess = 0.04
    def f_tmp(x):
        y = (-2*np.log10((2.51/(Re*np.sqrt(x))) + (e/(3.71))) - 1.0/np.sqrt(x))
        return y
    y = scipy.optimize.fsolve(f_tmp, f_guess)
    return y
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

def T_mo_T_infty(T_infty,T_mi,P,L,mdot,Cp,R_tot):
    return T_infty-(Tinfty-T_mi)*np.exp(-1/(mdot*Cp*Rtot))

def hbar_laminar_isothermal(k,D):
    return 3.66*k/D

def hbar_laminar_isoflux(k,D):
    return 4.36*k/D

def Nu_turbulent_Dittus_Boelter(Re,Pr,mode):
    if (mode == 'heating'):
        n = 0.4
    elif (mode == 'cooling'):
        n = 0.3
    return 0.023*Re**(4./5.)*Pr**n

def Nu_turbulent_Sieder_Tate(Re,Pr,mu,mu_s):
    return 0.027*Re**(4./5.)*Pr*(1./3.)*(mu/mu_s)**0.14

def Nu_turbulent_Gnielinski(Re,Pr,f):
    return (f/8.)*(Re-1000.)*Pr/(1+12.7*(f/8.)**0.5*(Pr**(2./3.)-1.))

def Nu_turbulent_Skupinski(Re,Pr):
    return 4.82+0.0185*(Re*Pr)**0.827

def Nu_turbulent_Seban(Re,Pr):
    return 5.0+0.025*(Re*Pr)**0.8