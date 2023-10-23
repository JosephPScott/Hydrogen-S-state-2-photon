# -*- coding: utf-8 -*-
"""
Last edited on Monday 17 July 2023

@author: joseph.p.scott@durham.ac.uk

A short script which plots the variation of the perturbalitve values of the inelastic scattering rate and two-photon ionisation rates out of the 2S state at the 1S--2S magic wavelengths.
"""

import numpy as np
import matplotlib.pyplot as plt
import Main_calculation_functions as imp
import matplotlib.font_manager as font_manager

nmax = 300
k = 0.1


#Plot the large scale log plot
Fig, ax = plt.subplots(1)

depths = np.arange(1, 100)
ion_rate = 31.6
scatter_rate = 61.49
#ion_rate = 28.33
#scatter_rate = 76.78

ax.plot(depths, [d**2*ion_rate for d in depths], label="Ionisation")
ax.plot(depths, [d*scatter_rate for d in depths], linestyle="dashed", label="Inelastic scattering rate")
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel("Depth", fontfamily="serif", fontsize=14)
ax.set_ylabel("Rate", fontfamily="serif", fontsize=14)

ax.tick_params(axis="x", labelsize=12)
ax.tick_params(axis="y", labelsize=12)

ax.legend()

plt.show()

##### This second plot is not used in the paper any more, but this is left in place for user interest.

#Set up the wavelengths and rates at D = 1Er (and at 1E8 W/cm^2 for ionisation rates), these have been calculated elsewhere
wavs = [514.646, 443.212, 414.483, 399.451]
ion_rates = [2.777E4, 1.442E4, 1.027E4, 8.132E3]
inelastic_rates = [61.5, 69.4, 74.0, 76.8]

depths = np.arange(1, 4, 0.02)

#Convert between depth and intensity for the ionisation rates. Returns scaling for the ionisation rate per unit depth.
def get_depth(wave):
    Erec = (imp.h**2)/(2*(imp.mp+imp.me)*(wave*10**(-9))**2*imp.Eh)
    pol = imp.S_pol(2, wave, nmax, k)
    I0 = 1.5536611486546777e-08
    return (2*np.pi*imp.fsc*abs(pol)*I0)/Erec

font = font_manager.FontProperties(family="serif", size=10)

Fig, ((ax_515, ax_443), (ax_414, ax_399)) = plt.subplots(2, 2, sharex=True, sharey=True, tight_layout =True)

ion_rate = np.zeros(4)

#Sets up the ionisation rate at 1Er depth at each magic wavelength
for i in range(0, 4):
    initial_depth = get_depth(wavs[i])
    ion_rate[i] = ion_rates[i]/(initial_depth**2)
    #dcrit = inelastic_rates[i]/ion_rate[i]
    #print(f"{wavs[i]} and {ion_rate[i]}, {dcrit}")

#Ioniation rates scale as depth^2 whilst inelastic rates scale as I
ax_515.plot(depths, [d**2*ion_rate[0] for d in depths], label="Ionisation")
ax_515.plot(depths, [d*inelastic_rates[0] for d in depths], label="Inelastic", linestyle="dashed")
ax_515.plot([1.948, 1.948], [10, 510], color="black", linestyle="dotted", linewidth=1.5)

ax_443.plot(depths, [d**2*ion_rate[1] for d in depths])
ax_443.plot(depths, [d*inelastic_rates[1] for d in depths], linestyle="dashed")
ax_443.plot([2.412, 2.412], [10, 510], color="black", linestyle="dotted", linewidth=1.5)

ax_414.plot(depths, [d**2*ion_rate[2] for d in depths])
ax_414.plot(depths, [d*inelastic_rates[2] for d in depths], linestyle="dashed")
ax_414.plot([2.819, 2.819], [10, 510], color="black", linestyle="dotted", linewidth=1.5)

ax_399.plot(depths, [d**2*ion_rate[3] for d in depths], label="Ionisation")
ax_399.plot(depths, [d*inelastic_rates[3] for d in depths], linestyle="dashed", label="Inelastic")
ax_399.plot([3.223, 3.223], [10, 510], color="black", linestyle="dotted", linewidth=1.5)

#ax_515.set_yscale("log")
ax_515.set_ylim(10, 510)
#ax_443.set_yscale("log")
ax_443.set_ylim(10, 510)
#ax_414.set_yscale("log")
ax_414.set_ylim(10, 510)
#ax_399.set_yscale("log")
ax_399.set_ylim(10, 510)

#ax_399.legend(prop=font)

ax_399.set_xlabel("Depth $E_{rec}$", fontfamily="serif", fontsize=14)
ax_414.set_xlabel("Depth $E_{rec}$", fontfamily="serif", fontsize=14)

ax_515.set_ylabel("Rate $[s^{-1}]$", fontfamily="serif", fontsize=14)
ax_414.set_ylabel("Rate $[s^{-1}]$", fontfamily="serif", fontsize=14)

ax_515.set_title(f"{wavs[0]} nm", fontfamily="serif", fontsize=14)
ax_443.set_title(f"{wavs[1]} nm", fontfamily="serif", fontsize=14)
ax_414.set_title(f"{wavs[2]} nm", fontfamily="serif", fontsize=14)
ax_399.set_title(f"{wavs[3]} nm", fontfamily="serif", fontsize=14)

ax_399.tick_params(axis="x", labelsize=12)
ax_414.tick_params(axis="x", labelsize=12)

ax_414.tick_params(axis="y", labelsize=12)
ax_515.tick_params(axis="y", labelsize=12)

plt.show()
