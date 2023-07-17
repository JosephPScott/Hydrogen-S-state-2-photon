# -*- coding: utf-8 -*-
"""
Last edited on Monday 17 July 2023

@author: joseph.p.scott@durham.ac.uk

A short script which creates a plot of the 1S and 2S polarisabilities in hydrogen across a specified range of wavelengths.
The intended range of wavelengths can be changed by editing "wlist"
The plots automatically use a symmetric-log scaling becasuse of the large difference between the polarisabilities of 1S and 2S in this regime

"""
# Import packages
import Main_calculation_functions as imp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import Locator

# Variables define the basis of radial Sturmian functions used in calculations
nmax = 300
k = 0.3

wlist = np.arange(395, 1000, 0.2) # Sets up the wavelength range to calculate across

###############################################################################
#The S_pol function from imp is designed to calculate the polarisability of a given state in a self contained way at a specific wavelength
#Since we are calculating the polarisability of the same states (with the same angular momentum) across multiple wavelengths, the matrix equations hardly change
#Thus, it is more efficient to represent these matricies once and to calculate the polarisability using the interal functions available in imp.
Es = imp.Schrodinger((0, 1), nmax, k)[0]
E1, V1 = imp.Schrodinger((1, 0), nmax, k)
E2, V2 = imp.Schrodinger((2, 0), nmax, k)

za = imp.Zplus(nmax, 0, k)
H = imp.Hamiltonian(nmax, 1, k)
T = imp.Overlap(nmax, 1, k)
zb = imp.Zminus(nmax, 1, k)

lhs1 = np.dot(V1, zb)
lhs2 = np.dot(V2, zb)
rhs1 = np.dot(za, V1)
rhs2 = np.dot(za, V2)

flist = [(imp.h*imp.c)/(w*10**(-9)*imp.Eh) for w in wlist]
pol1 = np.array([imp.Implicit_step(H, T, lhs1, rhs1, E1, f)/3 + imp.Implicit_step(H, T, lhs1, rhs1, E1, -1*f)/3 for f in flist])# Remembering the angular 1/3 factor
pol2 = np.array([imp.Implicit_step(H, T, lhs2, rhs2, E2, f)/3 + imp.Implicit_step(H, T, lhs2, rhs2, E2, -1*f)/3 for f in flist])

# There are poles on resonances, putting these poles in manually prevents the line which connects +infinty to -infinity across the pole.
# If changing the wavelength range, more may need to be considered
pol2[min(range(len(wlist)), key=lambda i : abs(wlist[i] - (imp.h*imp.c*10**(9))/((Es[1] - E2)*imp.Eh)))] = np.nan #2S - 3P
pol2[min(range(len(wlist)), key=lambda i : abs(wlist[i] - (imp.h*imp.c*10**(9))/((Es[2] - E2)*imp.Eh)))] = np.nan #2S - 4P
pol2[min(range(len(wlist)), key=lambda i : abs(wlist[i] - (imp.h*imp.c*10**(9))/((Es[3] - E2)*imp.Eh)))] = np.nan #2S - 5P
pol2[min(range(len(wlist)), key=lambda i : abs(wlist[i] - (imp.h*imp.c*10**(9))/((Es[4] - E2)*imp.Eh)))] = np.nan #2S - 6P
pol2[min(range(len(wlist)), key=lambda i : abs(wlist[i] - (imp.h*imp.c*10**(9))/((Es[5] - E2)*imp.Eh)))] = np.nan #2S - 7P

###############################################################################
#This class improves the ticks present for the symmetric log scaling of pyplot.
#It is taken from https://stackoverflow.com/a/20495928 and no ownership is claimed.
class MinorSymLogLocator(Locator):
    """
    Dynamically find minor tick positions based on the positions of
    major ticks for a symlog scaling.
    """
    def __init__(self, linthresh):
        """
        Ticks will be placed between the major ticks.
        The placement is linear for x between -linthresh and linthresh,
        otherwise its logarithmically
        """
        self.linthresh = linthresh

    def __call__(self):
        'Return the locations of the ticks'
        majorlocs = self.axis.get_majorticklocs()

        # iterate through minor locs
        minorlocs = []

        # handle the lowest part
        for i in range(1, len(majorlocs)):
            majorstep = majorlocs[i] - majorlocs[i-1]
            if abs(majorlocs[i-1] + majorstep/2) < self.linthresh:
                ndivs = 10
            else:
                ndivs = 9
            minorstep = majorstep / ndivs
            locs = np.arange(majorlocs[i-1], majorlocs[i], minorstep)[1:]
            minorlocs.extend(locs)

        return self.raise_if_exceeds(np.array(minorlocs))

    def tick_values(self, vmin, vmax):
        raise NotImplementedError('Cannot get tick locations for a '
                                  '%s type.' % type(self))


###############################################################################
#General plotting functions
FigDiffPol = plt.figure()
AxDiffPol = FigDiffPol.add_subplot()
AxDiffPol.set_yscale('symlog', linthreshy=1)
AxDiffPol.plot([395, 1000], [0,0], color="black", linewidth=1.5)
AxDiffPol.plot(wlist, pol1, linewidth=1.5, label="1S")
AxDiffPol.plot(wlist, pol2, linewidth=1.5, label="2S")

#Plots lines for the resonances
AxDiffPol.plot([656.3, 656.3], [-1000,1000], linestyle="dotted", color="black", linewidth=2)
AxDiffPol.plot([486.7, 486.7], [-1000,1000], linestyle="dotted", color="black", linewidth=2)
AxDiffPol.plot([434.8, 434.8], [-1000,1000], linestyle="dotted", color="black", linewidth=2)
AxDiffPol.plot([410.7, 410.7], [-1000,1000], linestyle="dotted", color="black", linewidth=2)
AxDiffPol.plot([397.9, 397.9], [-1000,1000], linestyle="dotted" ,color="black", linewidth=2)

#AxDiffPol.set_title("1S polarisability")
AxDiffPol.set_ylabel(r"$\alpha$[a.u.]")
AxDiffPol.set_xlabel(r"$\lambda$[nm]")
AxDiffPol.set_ylim(-1000, 1000)
AxDiffPol.set_xlim(395, 1000)

#Corrects the ticks associated with the symmetric log scaling
yaxis = AxDiffPol.yaxis
yaxis.set_minor_locator(MinorSymLogLocator(1))

AxDiffPol.legend()
plt.show()
