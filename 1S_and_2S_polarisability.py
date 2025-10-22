# -*- coding: utf-8 -*-
"""
Last edited on Monday 29 September 2025

@author: joseph.p.scott@durham.ac.uk

This code is made avaiable under the CC0 1.0 Universal license

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

wlist_1S = np.arange(395, 1000, 0.2) # Sets up the wavelength range to calculate across for the 1S state all values in nm
wlist_2S = np.arange(395, 1000, 0.2) # For the 2S state

w_min = min(min(wlist_1S), min(wlist_2S)) # finds the minimum wavelength value
w_max = max(max(wlist_1S), max(wlist_2S)) # finds the maximum wavelength value

plot_range = 500
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

flist_1S = [(imp.h*imp.c)/(w*10**(-9)*imp.Eh) for w in wlist_1S]
flist_2S = [(imp.h*imp.c)/(w*10**(-9)*imp.Eh) for w in wlist_2S]
pol1 = np.array([imp.Implicit_step(H, T, lhs1, rhs1, E1, f)/3 + imp.Implicit_step(H, T, lhs1, rhs1, E1, -1*f)/3 for f in flist_1S])# Remembering the angular 1/3 factor
pol2 = np.array([imp.Implicit_step(H, T, lhs2, rhs2, E2, f)/3 + imp.Implicit_step(H, T, lhs2, rhs2, E2, -1*f)/3 for f in flist_2S])

# There are poles on resonances, this prevents the line which connects +infinty to -infinity across the pole fpr the most notable.
for j in range(1, 7):
    pol2[min(range(len(wlist_2S)), key=lambda i : abs(wlist_2S[i] - (imp.h*imp.c*10**(9))/((Es[j] - E2)*imp.Eh)))] = np.nan #2S - (2+i)P
    pol1[min(range(len(wlist_1S)), key=lambda i : abs(wlist_1S[i] - (imp.h*imp.c*10**(9))/((Es[j-1] - E1)*imp.Eh)))] = np.nan #1S - (1+i)P

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
FigDiffPol = plt.figure(figsize=(9,3), dpi=100, tight_layout =True)
AxDiffPol = FigDiffPol.add_subplot()
AxDiffPol.set_yscale('symlog', linthresh=1)
AxDiffPol.plot([w_min, w_max], [0,0], color="black", linewidth=1.5)
AxDiffPol.plot(wlist_1S, pol1, linewidth=1.5, label="1S", linestyle="solid")
AxDiffPol.plot(wlist_2S, pol2, linewidth=1.5, label="2S")

#Plots lines for the resonances
for j in range(1, 7):
    resonance_2S = imp.h*imp.c*10**9/((Es[j]-E2)*imp.Eh)
    AxDiffPol.plot([resonance_2S, resonance_2S], [-plot_range, plot_range], linestyle="dotted", color="black", linewidth=2)
    resonance_1S = imp.h*imp.c*10**9/((Es[j-1]-E1)*imp.Eh)
    AxDiffPol.plot([resonance_1S, resonance_1S], [-plot_range, plot_range], linestyle="dotted", color="black", linewidth=2)

#AxDiffPol.set_title("1S polarisability")
AxDiffPol.set_ylabel(r"$\alpha$/a.u.")
AxDiffPol.set_xlabel(r"$\lambda$/nm")
AxDiffPol.set_ylim(-plot_range, plot_range)
AxDiffPol.set_xlim(w_min, w_max)

#Corrects the ticks associated with the symmetric log scaling
yaxis = AxDiffPol.yaxis
yaxis.set_minor_locator(MinorSymLogLocator(1))

AxDiffPol.legend()
plt.show()
