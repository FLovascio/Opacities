import functions_op as fop
import matplotlib.pyplot as plt
from scipy import optimize
import autograd.numpy as np
import scipy as sp
import mpmath as mpm
import autograd as ag


# To make the dust grid
def sizes_distr(smin,smax,ndust):
    sdust = np.zeros(ndust)
    sdb=np.logspace(np.log10(smin),np.log10(smax),ndust+1)
    for i in range(0,ndust):
        sdust[i]=np.sqrt(sdb[i]*sdb[i+1])
    return [sdust,sdb]

def SumFunc_Generic(delta,sigma,sigma_eff):
    n=len(delta)
    return np.sum(delta*(sigma-sigma_eff)/(sigma+((n-1)*sigma_eff)))

def makeProblem(delta_array,sigma_array):
    f=lambda sigma_e : SumFunc_Generic(delta_array,sigma_array,sigma_e)
    return f

smin  =  2.5e-7 # Min grain size
smax  =  5e-4    # Max grain size
ndust = 500    # Number of dust sizes
epsilon_0=0.013986 #Dust-to-gas ratio
sdust,sdb=sizes_distr(smin,smax,ndust) 

# Semenov 2003 is the reference for the opacity
# https://www2.mpia-hd.mpg.de/~semenov/Opacities/opacities.html
N_POL=np.zeros(ndust)  
for j in range(0,ndust):
    N_POL[j]=fop.mrn_pollack(sdust[j]) #To use the same distribution as in Semenov, it's not exactly an MRN
    
N_POL=N_POL/(sdust**3.) #Mass distribution
epsilondust=epsilon_0*N_POL/np.sum(N_POL)

#Semenov opacities
filename   = 'semenov.dat'
data       = np.loadtxt(filename,dtype = np.float128)
nusem      = data[:,0]
Kappasem   = data[:,1]


#Optical constants
dict_nrm=fop.dic_comp('NRM')
cons='new'
data=np.loadtxt(cons+"_cons/Normal_silicates/n_olivine.dat")
n_olivine=data[:,1]+1j*data[:,2]
data=np.loadtxt(cons+"_cons/Normal_silicates/n_orthopyroxene.dat")
n_pyroxene=data[:,1]+1j*data[:,2]
data =np.loadtxt(cons+"_cons/Normal_silicates/n_iron.dat")
n_iron =data[:,1]+1j*data[:,2]
data =np.loadtxt(cons+"_cons/Normal_silicates/n_troilite.dat")
n_troilite =data[:,1]+1j*data[:,2]
data =np.loadtxt(cons+"_cons/Normal_silicates/n_volatile_organics.dat")
n_organics =data[:,1]+1j*data[:,2]
data =np.loadtxt(cons+"_cons/Normal_silicates/n_water_ice.dat")
n_Water_ice=data[:,1]+1j*data[:,2]
lambda_k= data[:,0]*1e-4
nu=(sp.constants.speed_of_light*1e2)/lambda_k

dict_nrm=fop.dic_comp('NRM')
frac_olivine=dict_nrm["frac_olivine"]
frac_iron=dict_nrm["frac_iron"]
frac_pyroxene=dict_nrm["frac_pyroxene"]
frac_organics=dict_nrm["frac_organics"]
frac_Troilite=dict_nrm["frac_Troilite"]
frac_Water_ice=dict_nrm["frac_Water_ice"]

rho_olivine=dict_nrm["dens_olivine"]
rho_iron=dict_nrm["dens_iron"]
rho_pyroxene=dict_nrm["dens_pyroxene"]
rho_organics=dict_nrm["dens_organics"]
rho_Troilite=dict_nrm["dens_Troilite"]
rho_Water_ice=dict_nrm["dens_Water_ice"]



Frac_array=np.asarray([frac_olivine,frac_iron,frac_pyroxene,frac_organics,frac_Troilite,frac_Water_ice])
norm=np.sum(Frac_array)
Frac_array/=norm
sigma_array=np.asarray([n_olivine,n_iron,n_pyroxene,n_organics,n_troilite,n_Water_ice])  
rho_grain = (rho_olivine*frac_olivine+ rho_iron*frac_iron + rho_pyroxene*frac_pyroxene+rho_Troilite*frac_Troilite+rho_organics*frac_organics+rho_Water_ice*frac_Water_ice)/norm

#No ices
#Frac_array=np.asarray([frac_olivine,frac_iron,frac_pyroxene,frac_organics,frac_Troilite])
#norm=np.sum(Frac_array)
#Frac_array/=norm
#sigma_array=np.asarray([n_olivine,n_iron,n_pyroxene,n_organics,n_troilite])  
#rho_grain = (rho_olivine*frac_olivine+ rho_iron*frac_iron + rho_pyroxene*frac_pyroxene+rho_Troilite*frac_Troilite+rho_organics*frac_organics)/norm


#No organics
#Frac_array=np.asarray([frac_olivine,frac_iron,frac_pyroxene,frac_Troilite,frac_Water_ice])
#norm=np.sum(Frac_array)
#Frac_array/=norm
#sigma_array=np.asarray([n_olivine,n_iron,n_pyroxene,n_troilite,n_Water_ice])  
#rho_grain = (rho_olivine*frac_olivine+ rho_iron*frac_iron + rho_pyroxene*frac_pyroxene+rho_Troilite*frac_Troilite+rho_Water_ice*frac_Water_ice)/norm

#No pyroxen
#Frac_array=np.asarray([frac_olivine,frac_iron,frac_organics,frac_Troilite,frac_Water_ice])
#norm=np.sum(Frac_array)
#Frac_array/=norm
#sigma_array=np.asarray([n_olivine,n_iron,n_organics,n_troilite,n_Water_ice])  
#rho_grain = (rho_olivine*frac_olivine+ rho_iron*frac_iron +rho_Troilite*frac_Troilite+rho_organics*frac_organics+rho_Water_ice*frac_Water_ice)/norm
#n_mixture=np.zeros(len(n_Water_ice),dtype=complex)

#No troilite
Frac_array=np.asarray([frac_olivine,frac_iron,frac_pyroxene,frac_organics,frac_Water_ice])
norm=np.sum(Frac_array)
Frac_array/=norm
sigma_array=np.asarray([n_olivine,n_iron,n_pyroxene,n_organics,n_Water_ice])  
rho_grain = (rho_olivine*frac_olivine+ rho_iron*frac_iron + rho_pyroxene*frac_pyroxene+rho_organics*frac_organics+rho_Water_ice*frac_Water_ice)/norm


n_mixture=np.zeros(len(n_olivine),dtype=complex)
for i in range(len(n_mixture)):
    f=makeProblem(Frac_array,sigma_array[:,i])
    mixture=mpm.findroot(f,x0=1.0+1.0j)
    #print(mixture)
    n_mixture[i]=mixture

ee1 = fop.e1(np.real(n_mixture),np.imag(n_mixture))
ee2 = fop.e2(np.real(n_mixture),np.imag(n_mixture))
sigma= fop.sigmajk(lambda_k,ee1,ee2,0.333333333)

Kappadust=np.zeros(len(ee1))


for idust in range(ndust):
    Kappadust+=fop.Kappaj(epsilondust[idust],rho_grain,fop.Hj(fop.xj(lambda_k,sdust[idust],sdust[idust]),ee1,ee2),sigma,sigma,sigma)
plt.figure()
plt.loglog(nusem,Kappasem)
plt.loglog(nu,Kappadust)
plt.show()
