import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt("energy_3x2.txt")
plt.scatter(data[:,0], data[:,1], c='black', label="ED")
plt.errorbar(data[:,0],data[:,2],yerr=data[:,3], fmt='x', c="green", label="METTS")
plt.legend(loc="best")
plt.xlabel(r"$\beta$")
plt.ylabel(r"$E$")
plt.savefig("energy_3x2.png", dpi=300)
