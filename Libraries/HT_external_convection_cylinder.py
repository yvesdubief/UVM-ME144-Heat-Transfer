def Nu_Hilbert(Re,Pr):
    if (Re < 0.4) or (Re > 400000.):
        Nu = 0.
    else:
        if (Re < 4.):
            C = 0.989
            m = 0.33
        elif (Re < 40.):
            C = 0.911
            m = 0.385
        elif (Re < 4000.):
            C = 0.683
            m = 0.466
        elif (Re < 40000.):
            C = 0.193
            m = 0.618
        else:
            C = 0.027
            m = 0.805
        Nu = C*Re**m*Pr**(1./3.)
    return Nu

def Nu_Churchill_Bernstein(Re,Pr):
    
    if (Re*Pr < 0.2):
        Nu = 0.
    
    else:
        Nu = 0.3+(0.62*Re**(0.5)*Pr**(1./3.)) \
          /(1.+(0.4/Pr)**(2./3.))**(1./4.) \
        *(1.+(Re/282000.)**(5./8.))**(4./5.)
    return Nu

def Nu_Zukauskas(Re,Pr,Pr_s):
    if (Pr <= 10):
        n = 0.37
    else:
        n = 0.36
    inbound = True
    if (Re < 1.) and (Re > 1.e6):
        Nu = 0.
    else:
        if (Re < 40.):
            C = 0.75
            m = 0.4
        elif (Re < 1000.):
            C = 0.51
            m = 0.5
        elif (Re < 2.e5):
            C = 0.26
            m = 0.6
        else:
            C = 0.076
            m = 0.7
        Nu = C*Re**m*Pr**n*(Pr/Pr_s)**(1./4.)
    return Nu