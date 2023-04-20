# -*- coding: utf-8 -*-
"""
Last edited on Thur Apr 20 2023

@author: joseph.p.scott@durham.ac.uk

A module to produce a comparative plot of 2S Rayleigh and Raman scattering rates across a given spectral range.
The default spectral range is between 395 and 1100 nm.
Rates are calculated using the module "Main_calculation_functions" in the same project.
Parameters are controlled directly in the module by altering the relevant variable assignments.
"""

import Main_calculation_functions as imp
import numpy as np
import matplotlib.pyplot as plt

#These define the basis of radial Sturmian functions used in calculations
nmax = 300
k = 0.3

def get_values(wav): # Function returns both the Rayleigh and total Raman scattering rates 
    return(imp.S_Rayleigh_Depth(2, wav, nmax, k, 1), get_Raman(wav))

def get_Raman(wav): # Function to return the total Raman scattering rate 
    Srates = [imp.Ram_D_dep(2, (i, 0), wav, nmax, k, 1) for i in range(1, 7) if i != 2] #The largest number 7 should be increased if final states of n=7 or higher are possible
    Drates = [imp.Ram_D_dep(2, (i, 2), wav, nmax, k, 1) for i in range(3, 7)]
    return sum(Srates) + sum(Drates)

wavs = np.append(np.append(np.arange(395, 450, 0.2), np.arange(450, 670, 1)),np.arange(670, 1100, 5)) #List of wavelengths to calculate at, requires higher density at higher energy
rates = np.array([get_values(wav) for wav in wavs]) # calculates the rates at each wavelength

def determine(i): # imp.Ram_D_dep often produces a tuple, but some exceptional cases it might produce a float, this prevents exeptions in these cases.
    if type(i)==tuple:
        return sum(i)
    else:
        return i

#Returns individual vectors of for Rayleigh and total Raman scattering 
ralrates = np.array([rates[i, 0] for i in range(0, len(wavs))])
ramrates = np.array([determine(rates[i, 1]) for i in range(0, len(wavs))])

E2 = imp.Schrodinger((2, 0), nmax, k)[0] # gets the energy of the 2S state (this is known analytically in Rydberg theroy, but to solve it with the Sturmian process is more internally consistent)

#General plotting functions
Fig1 = plt.figure()
ax1 = Fig1.add_subplot()
ax1.plot(wavs, ralrates, linestyle="dashed", label="Rayleigh")
ax1.plot(wavs, ramrates, label="Quenching")

for i in range(3, 8): # plots lines relating to resonances in the 2S state, changing the plotting range may requrie this to be updated 
    line = imp.h*imp.c/(imp.Eh*abs(imp.E2 - imp.Schrodinger((i, 0), nmax, k)[0]))*10**9
    ax1.plot([line, line], [0.01, 200], linestyle="dotted", color="black")

magicx = [514.646, 443.212, 414.483, 399.451] # plots lines relating to the known 1S-2S magic wavelengths, changing the plotting range may requrie this to be updated 
magicy = [determine(get_Raman(w))/3 for w in magicx]
ax1.scatter(magicx, magicy, marker="x", color="red")

ax1.set_yscale("log")
ax1.set_xlabel("Wavelength [nm]")
ax1.set_ylabel("R/D [per s]")
ax1.set_ylim(0.01, 100) # neccessary due to the presence of poles and points of 0 in the rate structure
ax1.set_xlim(395, 1100) # these should ideally match the range of "wavs"
ax1.legend()

plt.show()