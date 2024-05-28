import aquila_borg as borg
import numpy as np
import matplotlib.pyplot as plt
import aquila as aq

N = 256
L=1000.
phi = np.fromfile("delta.bin", dtype=np.float64).reshape((N,N,N+2))
phi = phi[:,:,:N]


cpar = borg.cosmo.CosmologicalParameters()
cpar.omega_m =  0.3175 
cpar.omega_q = 0.6825   
cpar.omega_b = 0.049
cpar.h = 0.6711
cpar.sigma8 = 0.
cpar.n_s = 0.9624
#Redshift         127.       % Starting redshift
sigma8=0.834
#PrimordialIndex  0.9624       % may be used to tilt the primordial index, needed for nongaussian inital potential


print(cpar)
cc = borg.cosmo.ClassCosmo(cpar, 100,10., extra=dict(z_max_pk="200."))
cc.computeSigma8()
boost = (sigma8/cc.getCosmology()['sigma_8'])**2
cpar.A_s *= boost
cc = borg.cosmo.ClassCosmo(cpar, 100,10., extra=dict(z_max_pk="200."))
cc.computeSigma8()
print(cc.getCosmology())
cc.retrieve_Tk(127)
k, Pk= aq.clustering.compute_power_spectrum(phi, Lbox=L, Nk=256)

plt.loglog(k, Pk, label='IC')
plt.loglog(k, cc.get_Pk_matter(k),label='CLASS')
plt.legend()
plt.savefig("pk.png")

plt.clf()
plt.semilogx(k, Pk/cc.get_Pk_matter(k)) 
plt.axhline(1.0,lw=2.0,color='k')
plt.ylim(0.8,1.2)
plt.legend()
plt.savefig("ratio.png")
