# -*- coding: utf-8 -*-
"""
Last edited on Monday 17 July 2023

@author: joseph.p.scott@durham.ac.uk

This module contains a series of functions that are useful for evaluating the accuracy and precision of the Main_calculation_functions.
Accuracy is checked by a function that allows for the calculation of Raman scattering cross sections, comparison of the energy spectrum and eigenvectors with theoretical values, and the radial dipole matrix elements
Stability is assessed by plots of heatmaps showing the variation of calculated values with changes to the size of the Sturmian basis and the parameter k.
"""

import Main_calculation_functions as imp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import scipy.special as spec
import seaborn as sns
import pandas as pd

###############################################################################
#Function used to calculate the inelastic scattering cross sections for comparison to previous calculations.

def return_cross_section(n, fin_state, wav, nmax, k):
    energy= (imp.h*imp.c)/(wav*10**-9)
    if n==fin_state[0] and fin_state[1] == 0:
        Rate = imp.S_Rayleigh(n, wav, nmax, k, 1)
        return Rate*energy*10**4
    else:
        Rate = imp.S_Raman(n, fin_state, wav, nmax, k, 1)
        energy= (imp.h*imp.c)/(wav*10**-9)
        return (Rate[0]*energy*10**4, Rate[1]*energy*10**4)


###############################################################################
#Functions used to check the produced energy spectrum

def Check_EnergySpectrum(l, nmax, k): # Plots the correctness of the energy returned state by state up to n=20 and given l
    Es, Vs = imp.Schrodinger((0, l), nmax, k)
    Ediff = np.array([Es[n-l-1] + imp.redm/(2*n**2) for n in range(l+1, 21) if Es[n-l-1] < 0])
    Fig_Energy = plt.figure()
    ax_Energ = Fig_Energy.add_subplot()
    ax_Energ.plot(np.array([n for n in range(l+1, 21) if Es[n-l-1] < 0]), Ediff)
    ax_Energ.set_xlabel("n")
    ax_Energ.set_ylabel("Error")
    plt.xticks([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
    plt.show()
    return 


def Totalval(l, nmax, k): # Calculates the diferences between the analytic theory and the schrodinger equation for specific l 
    comps = imp.Schrodinger((0, l), nmax, k)[0]
    diffs = [abs(comps[i] + imp.redm/(2*(i+ l+ 1)**2)) for i in range(0, 10-l)] # Differences for specific l up to principle quantum number 10
    return diffs

def Avdiff(nmax, k): # Averages the differences between the analytic theory and the schrodinger equation for S, P and D states up to n=10 
    tot = 0
    num = 0
    for l in range(0, 3):
        diffs = Totalval(l, nmax, k)
        tot = tot + sum(diffs)
        num = num + len(diffs)
    return tot/num

def Spectrum_error_map(): # Plots a heatmap of the average correcness of the returned hydrogen spectrum with varying conditions on the sturmian basis
    ks = np.arange(0.1, 4.1, 0.2)
    ns = np.arange(200, 400, 10)
    
    vals = np.array([[Avdiff(nmax, k) for k in ks] for nmax in ns])
    df = pd.DataFrame(vals, columns=np.round(ks, 2), index=np.round(ns, 2))
    
    Fig_Spectrum_map = plt.figure()
    ax_Spectrum_map = Fig_Spectrum_map.add_subplot()
    # Plots the average difference as a heatmap with a log scale
    sns.heatmap(df, ax=ax_Spectrum_map, cmap="rocket_r", norm=LogNorm())
    ax_Spectrum_map.invert_yaxis()
    ax_Spectrum_map.set_xlabel("Free parameter")
    ax_Spectrum_map.set_ylabel("Number of Sturmian functions")

    plt.show()
    return df

###############################################################################
# Functions to check the vector representation of radial wavefunctions

def Check_Wavefunctions(st, nmax, k): # Plots the error between the analytic form of the radial wavefunction and that reconstructed from the Sturmian representation for given state "st"
    ESturm, VSturm = imp.Schrodinger(st, nmax, k)
    Rs = np.arange(0, 25*(st[0] - st[1]), 0.01) # The radius to plot over
    Diff = np.array([abs(r*Analytic_Wavefunction(r, st))**2 - abs(Sturmian_Wavefunction(r, st[1], k, VSturm))**2 for r in Rs])
    
    Fig_wavefunction = plt.figure()
    ax_wavefunction = Fig_wavefunction.add_subplot()
    ax_wavefunction.plot(Rs, Diff)
    ax_wavefunction.set_xlabel("r [a0]")
    ax_wavefunction.set_ylabel("Error")
    plt.show()
    return
    
def Sturmian_Wavefunction(r, l, k, Vrep): # Function calculates the radial wavefunction at position r from the represetnation across the Sturmian basis
    return sum(np.array([Vrep[i-1]*np.sqrt(np.math.factorial(i - 1)/np.math.factorial(i + 2*l))*(2*k*r)**(l+1)*np.exp(-k*r)*spec.assoc_laguerre(2*k*r, i - 1, 2*l + 1) for i in range(1, len(Vrep)+1)]))

def Analytic_Wavefunction(r, st):  # Function to calculare the he radial wavefunction at position r via the analytic solution
    return np.sqrt(np.math.factorial(st[0]-st[1]-1)/(2*st[0]*np.math.factorial(st[0]+st[1])))*((2*imp.redm)/st[0])**(3/2)*((2*r*imp.redm)/st[0])**st[1]*np.exp(-(r*imp.redm)/st[0])*spec.assoc_laguerre((2*r*imp.redm)/st[0], st[0]-st[1]-1, 2*st[1]+1)

###############################################################################
#Functions to check the dipole operators via the calculation of dipople matrix elements

def Calculate_DipoleElm(initial, final, nmax, k): # Function calculates the dipole matrix element between an intiial and final state, one must specify the initial l and its change
    l1 = initial[1]
    dl = final[1] - l1
    if dl == 1:
        z = imp.Zplus(nmax, l1, k)
        rhs = np.dot(z, initial)
        return np.dot(final, rhs)
    elif dl == -1:
        z = imp.Zminus(nmax, l1, k)
        rhs = np.dot(z, initial)
        return np.dot(final, rhs)
    else:
        return 0

def Dipole_stability(initial_st, final_st): # Plots a heatmap showing the variation of a particular dipole matrix element squared as k and nmax vary. Values are relative to that calculated at nmax = 300 and k = 0.2 by default
    ks = np.arange(0.1, 0.51, 0.025)
    ns = np.arange(290, 311, 1)
    #Compare to the value calculated at nmax = 300, k = 0.3
    base = abs(Calculate_DipoleElm(imp.Schrodinger(initial_st, 300, 0.3)[1], imp.Schrodinger(final_st, 300, 0.3)[1], 300, 0.3))**2
    #Calculates all the other values and puts into a data frame
    vals = np.array([[abs(Calculate_DipoleElm(imp.Schrodinger(initial_st, n, k)[1], imp.Schrodinger(final_st, n, k)[1], n, k))**2 - base for k in ks] for n in ns])
    df2 = pd.DataFrame(vals, columns=np.round(ks, 2), index=np.round(ns, 1))
    
    Fig_Magic_Wavelength = plt.figure()
    ax_Magic_Wavelength = Fig_Magic_Wavelength.add_subplot()
    # Plots a heatmap
    sns.heatmap(df2, ax=ax_Magic_Wavelength, cmap="rocket_r")
    ax_Magic_Wavelength.invert_yaxis()
    
    plt.show()

###############################################################################
#Functions to check the polarisability and values of magic wavelengths

def Polarisabilty_stability(n, wave): # Plots a heatmap showing the variation of given S state polarisation at a particular wavelength as k and nmax vary. Values are relative to that calculated at nmax = 300 and k = 0.2 by default
    ks = np.arange(0.1, 0.51, 0.025)
    ns = np.arange(290, 311, 1)
    #Compare to the value calculated at nmax = 300, k = 0.3
    base = imp.S_pol(n, wave, 300, 0.3)
    #Calculates all the other values and puts into a data frame
    vals = np.array([[abs(imp.S_pol(n, wave, nmax, k) - base) for k in ks] for nmax in ns])
    df = pd.DataFrame(vals, columns=np.round(ks, 2), index=np.round(ns, 1))
    
    Fig_Polarisability = plt.figure()
    ax_Polarisability = Fig_Polarisability.add_subplot()
    
    sns.heatmap(df, ax=ax_Polarisability, cmap="rocket_r") #for blue, change cmap="mako_r"
    ax_Polarisability.invert_yaxis()
    
    plt.show()

def Magic_wavelength_stability(wguess): # Plots a heatmap showing the variation of a particular 1S-2S magic wavelength as k and nmax vary. Values are relative to that calculated at nmax = 300 and k = 0.2 by default
    k2s = np.arange(0.1, 0.51, 0.025)
    n2s = np.arange(290, 311, 1)
    #Compare to the value calculated at nmax = 300, k = 0.3
    base = imp.Find_1S2S_magicwave(300, wguess, 0.3)
    #Calculates all the other values and puts into a data frame
    vals2 = np.array([[abs(imp.Find_1S2S_magicwave(nmax, 514, k) - base) for k in k2s] for nmax in n2s])
    df2 = pd.DataFrame(vals2, columns=np.round(k2s, 2), index=np.round(n2s, 1))
    
    Fig_Magic_Wavelength = plt.figure()
    ax_Magic_Wavelength = Fig_Magic_Wavelength.add_subplot()
    
    sns.heatmap(df2, ax=ax_Magic_Wavelength, cmap="rocket_r") #for blue, change cmap="mako_r"
    ax_Magic_Wavelength.invert_yaxis()
    
    plt.show()

###############################################################################
#Function to check the Raman scattering calculations

def Raman_scattering_stability(n, state, wave): # Plots a heatmap showing the variation of a particular Raman scattering rate as k and nmax vary. Values are relative to that calculated at nmax = 300 and k = 0.2 by default
    ks = np.arange(0.1, 0.51, 0.025)
    ns = np.arange(290, 311, 1)
    #Compare to the value calculated at nmax = 300, k = 0.3
    base = sum(imp.S_Raman_Depth(n, state, wave, 300, 0.3, 1))
    #Calculates all the other values and puts into a data frame
    vals = np.array([[abs(sum(imp.S_Raman_Depth(n, state, wave, nmax, k, 1)) - base) for k in ks] for nmax in ns])
    df = pd.DataFrame(vals, columns=np.round(ks, 2), index=np.round(ns, 1))
    
    Fig_Raman_Scattering = plt.figure()
    ax_Raman_Scattering = Fig_Raman_Scattering.add_subplot()
    
    sns.heatmap(df, ax=ax_Raman_Scattering, cmap="rocket_r") #for blue, change cmap="mako_r"
    ax_Raman_Scattering.invert_yaxis()
    
    plt.show()
