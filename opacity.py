import numpy as np
import functions_op as fop
import matplotlib.pyplot as plt
from scipy import optimize


# To make the dust grid
def sizes_distr(smin,smax,ndust):
    sdust = np.zeros(ndust)
    sdb=np.logspace(np.log10(smin),np.log10(smax),ndust+1)
    for i in range(0,ndust):
        sdust[i]=np.sqrt(sdb[i]*sdb[i+1])
    return [sdust,sdb]

def Bruggeman(eff,f,e1,e2):
    return f*(e1-eff)/(e1+2*eff)+(1.-f)*(e2-eff)/(e2+eff)
def SolvMix(fi,ei,em):
    f1=fi
    f2=1.0-fi
    b=(2.*f1-f2)*ei+(2.*f2-f1)*em
    root= 1./4.*(b+np.sqrt(8.*ei*em+b**2))
    
    return root

smin =2.5e-7 # Min grain size
smax=5e-4 # Max grain size
ndust=500 #Number of dust sizes
dict_nrm=fop.dic_comp('NRM')
epsilon_0=0.013986 #Dust-to-gas ratio
sdust,sdb=sizes_distr(smin,smax,ndust) 

# Semenov 2003 is the reference for the opacity
# https://www2.mpia-hd.mpg.de/~semenov/Opacities/opacities.html
N_POL=np.zeros(ndust)  
for j in range(0,ndust):
    N_POL[j]=fop.mrn_pollack(sdust[j]) #To use the same distribution as in Semenov, it's not exactly an MRN
    
N_POL=N_POL/(sdust**3.) #Mass distribution
epsilondust=epsilon_0*N_POL/np.sum(N_POL)

#Wavelenghts
lambda_sim=np.logspace(-1.,5.,100)

#Semenov opacities
filename = 'semenov.dat'
data       = np.loadtxt(filename,dtype = np.float128)
nusem  = data[:,0]
Kappasem = data[:,1]

#Optical constants
aaa =np.loadtxt("new_cons/Normal_silicates/n_olivine.dat")
nr_olivine =aaa[:,1]
ni_olivine=aaa[:,2]
aaa =np.loadtxt("new_cons/Normal_silicates/n_orthopyroxene.dat")
nr_pyroxene =aaa[:,1]
ni_pyroxene=aaa[:,2]
aaa =np.loadtxt("new_cons/Normal_silicates/n_iron.dat")
nr_iron =aaa[:,1]
ni_iron=aaa[:,2]
aaa =np.loadtxt("new_cons/Normal_silicates/n_troilite.dat")
nr_troilite =aaa[:,1]
ni_troilite=aaa[:,2]
aaa =np.loadtxt("new_cons/Normal_silicates/n_volatile_organics.dat")

nr_organics =aaa[:,1]
ni_organics=aaa[:,2]
aaa =np.loadtxt("new_cons/Normal_silicates/n_water_ice.dat")
lambdak = aaa[:,0]*1e-4
nr_Water_ice=aaa[:,1]
ni_Water_ice=aaa[:,2]



# Mass fractions of minerals
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

# Grain density
rho_grain = rho_olivine*frac_olivine+ rho_iron*frac_iron + rho_pyroxene*frac_pyroxene+rho_Troilite*frac_Troilite+rho_organics*frac_organics+rho_Water_ice*frac_Water_ice


ni = ni_olivine*frac_olivine+ni_iron*frac_iron+ni_pyroxene*frac_pyroxene+ni_troilite*frac_Troilite+ni_Water_ice*frac_Water_ice+ni_organics*frac_organics
nr = nr_olivine*frac_olivine+nr_iron*frac_iron+nr_pyroxene*frac_pyroxene+nr_troilite*frac_Troilite+nr_Water_ice*frac_Water_ice+nr_organics*frac_organics



Kappadust =np.zeros((ndust,len(lambdak)))
Kappad =np.zeros((len(lambdak)))
Kappadust_simple =np.zeros((ndust,len(lambdak)))
Kappad_simple =np.zeros((len(lambdak)))


clight = 3e10
nu = clight/lambdak
for j in range(0,ndust):
    for i in range(0,len(lambdak)):
        ee1_oliv = fop.e1(nr_olivine[i],ni_olivine[i])
        ee2_oliv = fop.e2(nr_olivine[i],ni_olivine[i])
        ee1_pyro = fop.e1(nr_pyroxene[i],ni_pyroxene[i])
        ee2_pyro = fop.e2(nr_pyroxene[i],ni_pyroxene[i])
        ee1_iron = fop.e1(nr_iron[i],ni_iron[i])
        ee2_iron = fop.e2(nr_iron[i],ni_iron[i])
        ee1_orga = fop.e1(nr_organics[i],ni_organics[i])
        ee2_orga = fop.e2(nr_organics[i],ni_organics[i])
        ee1_troi = fop.e1(nr_troilite[i],ni_troilite[i])
        ee2_troi = fop.e2(nr_troilite[i],ni_troilite[i])
        ee1_ice = fop.e1(nr_Water_ice[i],ni_Water_ice[i])
        ee2_ice = fop.e2(nr_Water_ice[i],ni_Water_ice[i])

        
        ee_oliv=complex(ee1_oliv,ee2_oliv)
        ee_pyro=complex(ee1_pyro,ee2_pyro)
        ee_iron=complex(ee1_iron,ee2_iron)
        ee_orga=complex(ee1_orga,ee2_orga)
        ee_troi=complex(ee1_troi,ee2_troi)
        ee_ice=complex(ee1_ice,ee2_ice)
        #H. Looyenga, "Dielectric const of heterogeneous mixtures",Physica 31, 401-406 [1965])
        folivine=(frac_olivine)
        fpyro=(frac_pyroxene)
        firon=(frac_iron)
        forga=(frac_organics)
        ftroi=(frac_Troilite)
        fice=(frac_Water_ice)
        
        
        f=folivine/(folivine+fpyro)
        eeff= SolvMix(f,ee_oliv,ee_pyro)

        f=(folivine+fpyro)/(folivine+fpyro+firon)
                
        eeff= SolvMix(f,eeff,ee_iron)

        f=(folivine+fpyro+firon)/(folivine+fpyro+firon+forga)
        eeff= SolvMix(f,eeff,ee_orga)

        f=(folivine+fpyro+firon+forga)/(folivine+fpyro+firon+forga+ftroi)
        eeff= SolvMix(f,eeff,ee_troi)

        f=(folivine+fpyro+firon+forga+ftroi)/(folivine+fpyro+firon+forga+ftroi+fice)
        eeff= SolvMix(f,eeff,ee_ice)
        
      

        eeff1=eeff.real
        eeff2=eeff.imag
        xdust=fop.xj(lambdak[i],sdust[j],sdust[j])
        Hdust=fop.Hj(xdust,eeff1,eeff2)
        sigma= fop.sigmajk(lambdak[i],eeff1,eeff2,1./3.)

        Kappadust_simple[j,i]=fop.Kappaj(epsilondust[j],rho_grain,Hdust,sigma,sigma,sigma)
        
        #olivine
        ee1 = fop.e1(nr_olivine[i],ni_olivine[i])
        ee2 = fop.e2(nr_olivine[i],ni_olivine[i])
        xdust=fop.xj(lambdak[i],sdust[j],sdust[j]) 
        Hdust=fop.Hj(xdust,ee1,ee2)
        sigma= fop.sigmajk(lambdak[i],ee1,ee2,1./3.)
        Kappadust[j,i]+=fop.Kappaj(epsilondust[j]*frac_olivine,rho_olivine,Hdust,sigma,sigma,sigma)

        #pyroxene
        ee1 = fop.e1(nr_pyroxene[i],ni_pyroxene[i])
        ee2 = fop.e2(nr_pyroxene[i],ni_pyroxene[i])
        Hdust=fop.Hj(xdust,ee1,ee2)
        sigma= fop.sigmajk(lambdak[i],ee1,ee2,1./3.)
        Kappadust[j,i]+=fop.Kappaj(epsilondust[j]*frac_pyroxene,rho_pyroxene,Hdust,sigma,sigma,sigma)

        #iron
        ee1 = fop.e1(nr_iron[i],ni_iron[i])
        ee2 = fop.e2(nr_iron[i],ni_iron[i])
        Hdust=fop.Hj(xdust,ee1,ee2)
        sigma= fop.sigmajk(lambdak[i],ee1,ee2,1./3.)
        Kappadust[j,i]+=fop.Kappaj(epsilondust[j]*frac_iron,rho_iron,Hdust,sigma,sigma,sigma)

        #organics
        ee1 = fop.e1(nr_organics[i],ni_organics[i])
        ee2 = fop.e2(nr_organics[i],ni_organics[i])
        Hdust=fop.Hj(xdust,ee1,ee2)
        sigma= fop.sigmajk(lambdak[i],ee1,ee2,1./3.)
        Kappadust[j,i]+=fop.Kappaj(epsilondust[j]*frac_organics,rho_organics,Hdust,sigma,sigma,sigma)

        #troilite
        ee1 = fop.e1(nr_troilite[i],ni_troilite[i])
        ee2 = fop.e2(nr_troilite[i],ni_troilite[i])
        Hdust=fop.Hj(xdust,ee1,ee2)
        sigma= fop.sigmajk(lambdak[i],ee1,ee2,1./3.)
        Kappadust[j,i]+=fop.Kappaj(epsilondust[j]*frac_Troilite,rho_Troilite,Hdust,sigma,sigma,sigma)

        #Water_ice
        ee1 = fop.e1(nr_Water_ice[i],ni_Water_ice[i])
        ee2 = fop.e2(nr_Water_ice[i],ni_Water_ice[i])
        Hdust=fop.Hj(xdust,ee1,ee2)
        sigma= fop.sigmajk(lambdak[i],ee1,ee2,1./3.)
        Kappadust[j,i]+=fop.Kappaj(epsilondust[j]*frac_Water_ice,rho_Water_ice,Hdust,sigma,sigma,sigma)

        
    Kappad+= Kappadust[j,:]
    Kappad_simple+= Kappadust_simple[j,:]
    
plt.loglog(nusem,Kappasem)

plt.loglog(nu,Kappad,label='small grains')
plt.loglog(nu,Kappad_simple,label='Bruggeman')

plt.xlim(np.amin(nu),np.amax(nu))
plt.legend()
plt.show()


