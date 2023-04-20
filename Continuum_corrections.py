# -*- coding: utf-8 -*-
"""
Last edited on Wed Apr 19 2023

@author: joseph.p.scott@durham.ac.uk

A short script evauating the effect of neglecting continuum contributions when calculating the scattering rates of 2S hydrogen.
The "full" values are calculated using the module "Main_calculation_functions" in the same project.
Discrete terms are calculated by numerically integrating to obtain the dipole matrix elements, energies are given by the Rydberg theory.
The main function "Run_plot" produces a plot comparing the discrete sum calculations to the full value
The properties of these calculations can be controlled by changing the appropriate values in the main function "Run_plot"
"""
#Import packages
import Main_calculation_functions as imp
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as spec
from scipy import integrate

#Function to calculate a particular scattering rate via summation over a specified number of discrete states
def Disc_Rate(n, l, wav, sign, nmax, k, number):
    freq = (imp.h*imp.c)/(wav*10**(-9)*imp.Eh)
    E_a = -imp.redm*(1/8)
    E_b = -imp.redm*(1/(2*n**2))
    sumterms = np.array([Disc_n(p, freq, n, l, E_a, E_b, sign) for p in range(2, number+1)]) # References the calculation of the sum over radial matrix elements
    melms = sum(sumterms)
    pol = imp.S_pol(2, wav, nmax, k)/3
    freq_sc = -(E_b - E_a) + sign*freq
    if l==0: # Assigns angular components - calculated analytically
        A = 1/9
    else:
        A = 2/9
    return ((freq_sc**3)*(freq**2)*2*A*(imp.fsc**5)*(melms**2))/(3*pol*(imp.x+1))*(imp.Eh/imp.hbar) # Outputs a scattering rate in units of per second

#Calculates the sum across radial matrix elements
def Disc_n(nP, freq, n, l, E_a, E_b, sign):
    E_n = -imp.redm*(1/(2*nP**2))
    melm1 = integrate.quad(lambda r: wavfunc(r, 2, 0)*wavfunc(r, nP, 1)*r**3, 0, np.inf)[0] # matrix elements are calculated by numerical integration
    melm2 = integrate.quad(lambda r: wavfunc(r, n, l)*wavfunc(r, nP, 1)*r**3, 0, np.inf)[0]
    melms = melm1*melm2
    return melms/(E_n - E_a - sign*freq) + melms/(E_n - E_b + sign*freq) 

def NormL(n, l): # Normalising function for the wavelength
    return np.sqrt(np.math.factorial(n-l-1)/(2*n*np.math.factorial(n+l)))

def wavfunc(r, n, l): # Hydrogen bound state wavefunctions
    return NormL(n, l)*((2*imp.redm/n)**(3/2))*((2*r*imp.redm/n)**l)*np.exp(-r*imp.redm/n)*spec.assoc_laguerre(2*r*imp.redm/n, n-l-1, 2*l + 1)


#The main function of the module which produces the final plot.
def Run_plot():
    wavs = [514.646, 443.212, 414.483, 399.451]
    numbers = range(2, 21)
    Fig1 = plt.figure()
    ax1 = Fig1.add_subplot()
    markers = ["o", "s", "x", "v"]
    plt.xticks(numbers)
    for i in range(0, 4):
        lim = Overall(wavs[i], [40])
        print(lim)
        line1 = ax1.plot([numbers[0], numbers[-1]], [lim, lim], linestyle="dashed")
        ax1.plot(numbers, Overall(wavs[i], numbers), label = str(wavs[i]), marker=markers[i], color = line1[0].get_color())
    ax1.set_xlabel("Discrete inclusion")
    ax1.set_ylabel("Fractional correctness")
    plt.legend()
    plt.show()
    return


#Functions used to produce the lists of values to be plotted:
def Overall(wav, numbers): #Produces the total scattering rate and outputs overall ratios
    fulls = sum(imp.Ram_D_dep(2, (1, 0), wav, 300, 0.3, 1)) + imp.Ram_D_dep(2, (3, 0), wav, 300, 0.3, 1)[0] + imp.Ram_D_dep(2, (3, 2), wav, 300, 0.3, 1)[0] # there follows an inelegent solution to inclusing all relevant final states
    if wav < 514:
        fulls = fulls + imp.Ram_D_dep(2, (4, 0), wav, 300, 0.3, 1)[0] + imp.Ram_D_dep(2, (4, 2), wav, 300, 0.3, 1)[0]
    if wav < 443:
        fulls = fulls + imp.Ram_D_dep(2, (5, 0), wav, 300, 0.3, 1)[0] + imp.Ram_D_dep(2, (5, 2), wav, 300, 0.3, 1)[0]
    if wav < 414:
        fulls = fulls + imp.Ram_D_dep(2, (6, 0), wav, 300, 0.3, 1)[0] + imp.Ram_D_dep(2, (6, 2), wav, 300, 0.3, 1)[0]
    return np.array([Over_each(wav, n)/fulls for n in numbers]) # Generates a list of these ratios with the number of states in the discrete case increasing with each step

def Over_each(wav, number): #For a given maximum number of discrete states, calculates the total scattering rate at a particular wavelength
    discs = Disc_Rate(1, 0, wav, 1, 300, 0.3, number)+ Disc_Rate(1, 0, wav, -1, 300, 0.3, number)+ Disc_Rate(3, 0, wav, 1, 300, 0.3, number)+ Disc_Rate(3, 2, wav, 1, 300, 0.3, number)
    if wav < 514:
        discs = discs + Disc_Rate(4, 0, wav, 1, 300, 0.3, number) + Disc_Rate(4, 2, wav, 1, 300, 0.3, number)
    if wav < 443:
        discs = discs + Disc_Rate(5, 0, wav, 1, 300, 0.3, number) + Disc_Rate(5, 2, wav, 1, 300, 0.3, number)
    if wav < 414:
        discs = discs + Disc_Rate(6, 0, wav, 1, 300, 0.3, number) + Disc_Rate(6, 2, wav, 1, 300, 0.3, number)
    return discs