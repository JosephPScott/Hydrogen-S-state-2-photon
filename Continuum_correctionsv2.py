
# -*- coding: utf-8 -*-
"""
Last edited on Monday 17 July 2023

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


import os, psutil
import time

nmax = 300
k = 0.3

def Overall(wav, number, E_int, V_int, H, T, zb, rhs, lhs_pol): #Produces the total scattering rate and outputs overall ratios
    tic = time.perf_counter()
    freq = (imp.h*imp.c)/(wav*10**(-9)*imp.Eh)
    pol = abs((imp.Implicit_step(H, T, lhs_pol, rhs, E_int, freq) + imp.Implicit_step(H, T, lhs_pol, rhs, E_int, -1*freq))/3)
    Srates = [sum(get_Raman((i, 0), E_int, freq, H, T, rhs, zb, pol)) for i in range(1, 7) if i != 2] #The largest number 7 should be increased if final states of n=7 or higher are possible
    Drates = [sum(get_Raman((i, 2), E_int, freq, H, T, rhs, zb, pol)) for i in range(3, 7)]
    fulls = sum(Srates) + sum(Drates)
    Discs = Running_Disc(wav, pol, number[0])
    Lim = (Limit_Disc(wav, pol, number[1])/fulls)
    memory = psutil.Process().memory_info().rss
    print(f"memory at {memory*(1064)**(-2)} MB, or {memory*(1064)**(-3)} GB")
    toc = time.perf_counter()
    print(f"Time for Raman process: {toc - tic} s")
    print(wav)
    return np.array([D/fulls for D in Discs]), Lim
    
def get_Raman(final_st, E_int, freq, H, T, rhs, zb, pol):
    E_fin, V_fin = imp.Schrodinger(final_st, nmax, k)
    freq_scatt = -(E_fin - E_int) + freq
    if freq_scatt  <= 0:
        return 0, 0
    if final_st[1] == 0: # In the case where the final state is an S state
        angular = 1/9
        lhs = np.dot(V_fin, zb)
    elif final_st[1] == 2: # In the case where the final state is a D state
        angular = 2/9
        lhs = np.dot(V_fin, imp.Zplus(nmax, 1, k))
    else: # Handles requests for dipole forbidden final states
        print("Forbidden")
        return (0, 0)
    P1 = (2*(freq**2)*(freq_scatt**3)*(imp.fsc**5)*abs(imp.Implicit_step(H, T, lhs, rhs, E_int, freq) + imp.Implicit_step(H, T, lhs, rhs, E_fin, -1*freq))**2)/(3*pol*((imp.mp/imp.me) + 1))
    freq_scatt = -(E_fin - E_int) - freq
    if freq_scatt <= 0:
        return angular*P1*(imp.Eh/imp.hbar), 0
    else:
        P2 = (2*(freq**2)*(freq_scatt**3)*(imp.fsc**5)*abs(imp.Implicit_step(H, T, lhs, rhs, E_int, -1*freq) + imp.Implicit_step(H, T, lhs, rhs, E_fin, freq))**2)/(3*pol*((imp.mp/imp.me) + 1))
        return angular*P1*(imp.Eh/imp.hbar), angular*P2*(imp.Eh/imp.hbar) #Outputs a tuple of both rates in units of per second

def Running_Disc(wav, pol, number): #For a given maximum number of discrete states, calculates the total scattering rate at a particular wavelength
    E_a = -imp.redm*(1/8)
    freq = (imp.h*imp.c)/(wav*10**(-9)*imp.Eh)
    discs = Disc_Rate(1, 0, freq, E_a, pol, 1, number)+ Disc_Rate(1, 0, freq, E_a, pol, -1, number)+ Disc_Rate(3, 0, freq, E_a, pol, 1, number)+ Disc_Rate(3, 2, freq, E_a, pol, 1, number)
    if wav < 514:
        discs = discs + Disc_Rate(4, 0, freq, E_a, pol, 1, number) + Disc_Rate(4, 2,freq, E_a, pol, 1, number)
    if wav < 443:
        discs = discs + Disc_Rate(5, 0, freq, E_a, pol, 1, number) + Disc_Rate(5, 2, freq, E_a, pol, 1, number)
    if wav < 414:
        discs = discs + Disc_Rate(6, 0, freq, E_a, pol, 1, number) + Disc_Rate(6, 2, freq, E_a, pol, 1, number)
    return discs

def Disc_Rate(n, l, freq, E_a, pol, sign, number):
    E_b = -imp.redm*(1/(2*n**2))
    melms = np.cumsum(np.array([Disc_n(p, freq, n, l, E_a, E_b, sign) for p in range(2, number+1)]))
    freq_sc = -(E_b - E_a) + sign*freq
    if l==0: # Assigns angular components - calculated analytically
        A = 1/9
    else:
        A = 2/9
    return np.array([((freq_sc**3)*(freq**2)*2*A*(imp.fsc**5)*(ms**2))/(3*pol*((imp.mp/imp.me)+1))*(imp.Eh/imp.hbar) for ms in melms]) # Outputs a scattering rate in units of per second

def Limit_Disc(wav, pol, number): #For a given maximum number of discrete states, calculates the total scattering rate at a particular wavelength
    E_a = -imp.redm*(1/8)
    freq = (imp.h*imp.c)/(wav*10**(-9)*imp.Eh)
    discs = Disc_Lim(1, 0, freq, E_a, pol, 1, number)+ Disc_Lim(1, 0, freq, E_a, pol, -1, number)+ Disc_Lim(3, 0, freq, E_a, pol, 1, number)+ Disc_Lim(3, 2, freq, E_a, pol, 1, number)
    if wav < 514:
        discs = discs + Disc_Lim(4, 0, freq, E_a, pol, 1, number) + Disc_Lim(4, 2,freq, E_a, pol, 1, number)
    if wav < 443:
        discs = discs + Disc_Lim(5, 0, freq, E_a, pol, 1, number) + Disc_Lim(5, 2, freq, E_a, pol, 1, number)
    if wav < 414:
        discs = discs + Disc_Lim(6, 0, freq, E_a, pol, 1, number) + Disc_Lim(6, 2, freq, E_a, pol, 1, number)
    return discs

def Disc_Lim(n, l, freq, E_a, pol, sign, number):
    E_b = -imp.redm*(1/(2*n**2))
    melms = sum(np.array([Disc_n(p, freq, n, l, E_a, E_b, sign) for p in range(2, number+1)]))
    freq_sc = -(E_b - E_a) + sign*freq
    if l==0: # Assigns angular components - calculated analytically
        A = 1/9
    else:
        A = 2/9
    return ((freq_sc**3)*(freq**2)*2*A*(imp.fsc**5)*(melms**2))/(3*pol*((imp.mp/imp.me)+1))*(imp.Eh/imp.hbar) # Outputs a scattering rate in units of per second

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
    numbers = (20, 100)

    E2, V2 = imp.Schrodinger((2, 0), nmax, k)
    H = imp.Hamiltonian(nmax, 1, k)
    T = imp.Overlap(nmax, 1, k)
    zb = imp.Zminus(nmax, 1, k)
    rhs = np.dot(imp.Zplus(nmax, 0, k), V2)
    lhs_pol = np.dot(V2, zb)
    
    Fig1 = plt.figure()
    ax1 = Fig1.add_subplot()
    markers = ["o", "s", "x", "v"]
    plt.xticks(range(2, 21))
    for i in range(0, 4):
        Terms, lim = Overall(wavs[i], numbers, E2, V2, H, T, zb, rhs, lhs_pol)
        print(lim)
        line1 = ax1.plot([2, 20], [lim, lim], linestyle="dashed")
        ax1.plot(range(2, 21), Terms, label = str(wavs[i]), marker=markers[i], color = line1[0].get_color())
    ax1.set_xlabel("Discrete inclusion")
    ax1.set_ylabel("Fractional correctness")
    ax1.set_ylim(0, 1)
    plt.legend()
    plt.show()
    return



