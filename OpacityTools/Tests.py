import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc

os.system("./bin/OpacityTests")
refValues=np.genfromtxt("./ReferenceValues/reference_data.csv",delimiter=",")
#resultsKappa=np.fromfile("./bin/KDust.bin",dtype="double")
resultsKappa = np.genfromtxt("./bin/KDust.csv",delimiter=" ")
resultsLambda = np.genfromtxt("./bin/Lambda.csv",delimiter=" ")
resultsDust = np.genfromtxt("./bin/Dust.csv",delimiter=" ")
tpe = np.dtype([('i','<i2'), ('q','<i2')])
resultsNu = np.genfromtxt("./bin/EffectiveConductivities.csv",delimiter=" ")
resultsNu = np.asarray([np.complex(arV[0],arV[1]) for arV in resultsNu])
refValues_nu = np.loadtxt("./ReferenceValues/reference_data_nu.csv").view(complex)

#print(resultsNu)
plt.figure()
#plt.loglog(refValues[:,0],refValues[:,1])
plt.loglog((spc.speed_of_light*1e2)/(resultsLambda),resultsKappa)
plt.show()

#print(resultsNu)
plt.figure()
#plt.loglog(refValues[:,0],refValues[:,1])
plt.loglog(resultsDust)
plt.show()
"""
plt.figure()
plt.title(r"Re($\nu$) comparison")
plt.xlabel(r"$\lambda$")
plt.ylabel(r"Re($\nu$)")
plt.loglog(refValues[:,0],refValues_nu.real,label="Python")
plt.loglog(refValues[:,0],resultsNu.real,label="OpacityTools")
plt.legend()
plt.show()
"""
"""
plt.figure()
plt.title(r"Im($\nu$) comparison")
plt.xlabel(r"$\lambda$")
plt.ylabel(r"Im($\nu$)")
plt.loglog(refValues[:,0],refValues_nu.imag,label="Python")
plt.loglog(refValues[:,0],resultsNu.imag,label="OpacityTools")
plt.legend()
plt.show() """