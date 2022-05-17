import numpy as np

def Kappaj(epsilon,rhograin,H,sigma1,sigma2,sigma3): # Opacity of a sphere
    return epsilon/rhograin*H/3.*(sigma1+sigma2+sigma3) #For spheres sigmais are equal

#MIE functions
def Hj(xj,e1,e2):
    return 1. + (xj**2.*((e1+2.)**2.+e2**2.))/90.

def xj(lambdak,aj,bj):
    return 2*np.pi/lambdak*(1./3.*aj+2./3.*bj)

def sigmajk(lambdak,e1,e2,ljk):
    return  2*np.pi*e2/(lambdak*ljk**2.)/((e1+1./ljk-1)**2.+e2**2)

# Dielectric function
def e1(nr,ni):
    return nr**2.-ni**2

def e2(nr,ni):
    return 2.*nr*ni

# POLLACK MRN
def mrn_pollack(r):
    P0=0.005e-4
    if(r>=5e-4):
        return 0.
    if(r<5e-4):
        if(r>=1e-4):
            return (1/P0)**2*(P0/r)**5.5
        if(r<1e-4):
            if(r>=P0):
                return (P0/r)**3.5
            if(r<P0):
                return 1.
        

def dic_comp(regime):
    #fraction
    if(regime=='NRM'):
        olivine = 2.64e-3
        iron= 1.26e-4
        pyroxene=7.70e-4
        Troilite=7.68e-4
        Refractory_organics=3.53e-3
        Volatile_organics=6.02e-4
        Water_ice=5.55e-3
    if(regime=='IRS'):
        olivine = 3.84e-3
        iron= 0.
        pyroxene=4.44E-5
        Troilite=3.8E-4
        Refractory_organics=3.53e-3
        Volatile_organics=6.02e-4
        Water_ice=5.55e-3
    if(regime=='IPS'):
        olivine = 6.3E-4
        iron= 7.97e-4
        pyroxene=1.91e-3
        Troilite=7.68e-4
        Refractory_organics=3.53e-3
        Volatile_organics=6.02e-4
        Water_ice=5.55e-3

    organics=Volatile_organics#+Refractory_organics

    total = olivine+iron+pyroxene+Troilite+organics+Water_ice
    olivine = olivine/total
    iron = iron/total
    pyroxene= pyroxene/total
    Troilite=Troilite/total
    organics=organics/total
    Water_ice=Water_ice/total
    
    if(regime=='NRM'):
        rho_olivine = 3.49
        rho_iron= 7.87
        rho_pyroxene= 3.4
        rho_Troilite=4.83
        rho_Refractory_organics= 1.5
        rho_Volatile_organics=1.
        rho_Water_ice= 0.92
    if(regime=='IRS'):
        rho_olivine = 3.59
        rho_iron= 0.
        rho_pyroxene= 3.42
        rho_Troilite=4.83
        rho_Refractory_organics= 1.5
        rho_Volatile_organics= 1.
        rho_Water_ice= 0.92
    if(regime=='IPS'):
        rho_olivine = 3.20
        rho_iron= 7.87
        rho_pyroxene= 3.2
        rho_Troilite=4.84
        rho_Refractory_organics= 1.5
        rho_Volatile_organics= 1.
        rho_Water_ice= 0.92
    rho_organics= Volatile_organics*rho_Volatile_organics#+rho_Refractory_organics*Refractory_organics

    rho_tot = rho_olivine*olivine+ rho_iron*iron + rho_pyroxene*pyroxene+rho_Troilite*Troilite+rho_organics*organics+rho_Water_ice*Water_ice
    dictcomp ={}
    dictcomp["frac_olivine"]=olivine
    dictcomp["frac_iron"]=iron
    dictcomp["frac_pyroxene"]=pyroxene
    dictcomp["frac_Troilite"]=Troilite
    dictcomp["frac_organics"]=Volatile_organics
    dictcomp["organics"]=organics
    dictcomp["frac_Water_ice"]=Water_ice
    dictcomp["dens_olivine"]=rho_olivine
    dictcomp["dens_iron"]=rho_iron
    dictcomp["dens_pyroxene"]=rho_pyroxene
    dictcomp["dens_Troilite"]=rho_Troilite
    dictcomp["dens_organics"]=rho_organics
    dictcomp["dens_Water_ice"]=rho_Water_ice
    dictcomp["dens_tot"]=rho_tot
    return dictcomp
