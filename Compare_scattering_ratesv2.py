# -*- coding: utf-8 -*-
"""
Last edited on Monday 17 July 2023

@author: joseph.p.scott@durham.ac.uk

A module to produce a comparative plot of 2S Rayleigh and Raman scattering rates across a given spectral range.
The default spectral range is between 395 and 1000 nm.
Rates are calculated using the module "Main_calculation_functions" in the same project.
Parameters are controlled directly in the module by altering the relevant variable assignments.
"""

import Main_calculation_functions as imp
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd

#Memory management
import os, psutil
import time

#These define the basis of radial Sturmian functions used in calculations
nmax = 300
k = 0.3

"""
Creating this plot requires the repeated calculation of a various rates at sucessive wavelengths. 
The calculation functions are of imp are designed to be self-contained and perform calculations from scratch.
As such, to call these for each wavelength would result in a redundant recalculations of the same term
Therefore, the calculation functions are cannibalised here to make for more efficient calculation
"""
#Scattering rates for a given trap depth (as measured in any specified state)
def get_values(wave, E_int, V_int, H, T, zb, rhs, lhs_rate, lhs_pol, rhs_pol, E_1):
    #tic = time.perf_counter()
    freq = (imp.h*imp.c)/(wave*10**(-9)*imp.Eh)
    #We want to be able to plot for depths defined for other states (the 1S), so we calculate two versions of polarisability
    pol_deep = abs((imp.Implicit_step(H, T, lhs_pol, rhs_pol, E_1, freq) + imp.Implicit_step(H, T, lhs_pol, rhs_pol, E_1, -1*freq))/3)
    pol_rate = abs((imp.Implicit_step(H, T, lhs_rate, rhs, E_int, freq) + imp.Implicit_step(H, T, lhs_rate, rhs, E_int, -1*freq))/3)
    Erec = (imp.h**2)/(2*(imp.mp+imp.me)*(wave*10**(-9))**2*imp.Eh)
    rayleigh_rate = (freq**3)*(imp.fsc**3)*(4/3)*((pol_rate)**2/abs(pol_deep))*Erec*(imp.Eh/imp.hbar)
    Srates = [sum(get_Raman((i, 0), E_int, freq, H, T, rhs, zb, pol_deep)) for i in range(1, 7) if i != 2] #The largest number 7 should be increased if final states of n=7 or higher are possible
    Drates = [sum(get_Raman((i, 2), E_int, freq, H, T, rhs, zb, pol_deep)) for i in range(3, 7)]
    #Memory tracking and benchmarking, can be removed if preferred.
    #memory = psutil.Process().memory_info().rss
    #print(f"memory at {memory*(1064)**(-2)} MB, or {memory*(1064)**(-3)} GB")
    #toc = time.perf_counter()
    #print(f"Time for Raman process: {toc - tic} s")
    #print(wave)
    return (rayleigh_rate, (sum(Srates) + sum(Drates)))

#Scattering rates for a given intensity
def get_values_intent(wave, E_int, V_int, H, T, zb, rhs, lhs_pol):
    #tic = time.perf_counter()
    freq = (imp.h*imp.c)/(wave*10**(-9)*imp.Eh)
    pol = abs((imp.Implicit_step(H, T, lhs_pol, rhs, E_int, freq) + imp.Implicit_step(H, T, lhs_pol, rhs, E_int, -1*freq))/3)
    I = 5.23820810495585e-8 # This intesnity relates to 337 MW/cm^2
    rayleigh_rate = (freq**3)*(imp.fsc**4)*(8*np.pi/3)*(pol**2)*I*(imp.Eh/imp.hbar) 
    Srates = [sum(get_Raman_intensity((i, 0), E_int, freq, H, T, rhs, zb, I)) for i in range(1, 7) if i != 2] #The largest number 7 should be increased if final states of n=7 or higher are possible
    Drates = [sum(get_Raman_intensity((i, 2), E_int, freq, H, T, rhs, zb, I)) for i in range(3, 7)]
    #Memory tracking and benchmarking, can be removed if preferred.
    #memory = psutil.Process().memory_info().rss
    #print(f"memory at {memory*(1064)**(-2)} MB, or {memory*(1064)**(-3)} GB")
    #toc = time.perf_counter()
    #print(f"Time for Raman process: {toc - tic} s")
    #print(wave)
    return (rayleigh_rate, (sum(Srates) + sum(Drates)))

###############################################################################
#Sub-functions to the above two, used to calculate individual ramamn scattering terms at a given depth (in any specified state)...
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
    
    
#...or at  a constant intensity.
def get_Raman_intensity(final_st, E_int, freq, H, T, rhs, zb, I):
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
    P1 = ((freq_scatt**3)*(imp.fsc**4)*abs(imp.Implicit_step(H, T, lhs, rhs, E_int, freq) + imp.Implicit_step(H, T, lhs, rhs, E_fin, -1*freq))**2)*I*(8*np.pi/3)
    freq_scatt = -(E_fin - E_int) - freq
    if freq_scatt <= 0:
        return angular*P1*(imp.Eh/imp.hbar), 0
    else:
        P2 = ((freq_scatt**3)*(imp.fsc**4)*abs(imp.Implicit_step(H, T, lhs, rhs, E_int, -1*freq) + imp.Implicit_step(H, T, lhs, rhs, E_fin, freq))**2)*I*(8*np.pi/3)
        return angular*P1*(imp.Eh/imp.hbar), angular*P2*(imp.Eh/imp.hbar) #Outputs a tuple of both rates in units of per second)

###############################################################################
#function to get the intensity from one unit recoil depth in either the 2S...
def rec_depth(wave):
    Erec = (imp.h**2)/(2*(imp.mp+imp.me)*(wave*10**(-9))**2*imp.Eh)
    pol = imp.S_pol(2, wave, nmax, k)
    I0 = 1.5536611486546777e-08
    I = Erec/(2*np.pi*imp.fsc*abs(pol))
    return (I/I0)**2
# or 1S states.
def S1rec_depth(wave):
    Erec = (imp.h**2)/(2*(imp.mp+imp.me)*(wave*10**(-9))**2*imp.Eh)
    pol = imp.S_pol(1, wave, nmax, k)
    I0 = 1.5536611486546777e-08
    I = Erec/(2*np.pi*imp.fsc*abs(pol))
    return (I/I0)**2

###############################################################################
"""
Core plotting terms. 
The choice of wavelength range is made to maximise resolution across features of interest whilst decreasing resolution (and so computational load) where it is not needed.
"""
wavs = np.append(np.append(np.append(np.append(np.append(np.append(np.append(np.append(np.append(np.arange(395, 401, 0.2), np.arange(401, 408, 0.5)), np.arange(408, 416, 0.2)), np.arange(416, 432, 0.5)), np.arange(432, 446, 0.5)), np.arange(446, 484, 1)), np.arange(484, 488, 0.5)), np.arange(488, 654, 1)), np.arange(654, 658.2, 0.2)), np.arange(660, 1005, 5))

#shared terms which are common to calculations at all wavelengths and should be kept as global variables
E1, V1 = imp.Schrodinger((1, 0), nmax, k)
E2, V2 = imp.Schrodinger((2, 0), nmax, k)
H = imp.Hamiltonian(nmax, 1, k)
T = imp.Overlap(nmax, 1, k)
zb = imp.Zminus(nmax, 1, k)
rhs_2S = np.dot(imp.Zplus(nmax, 0, k), V2)
lhs_2S = np.dot(V2, zb)
lhs_1S = np.dot(V1, zb)
rhs_1S = np.dot(imp.Zplus(nmax, 0, k), V1)


# calculates the scattering rates at each wavelength for constant depths in 2S, 1S and intensity
rates_2S = np.array([get_values(wav, E2, V2, H, T, zb, rhs_2S, lhs_2S, lhs_2S, rhs_2S, E2) for wav in wavs]) 
rates_1S = np.array([get_values(wav, E2, V2, H, T, zb, rhs_2S, lhs_2S, lhs_1S, rhs_1S, E1) for wav in wavs])
rates_i = np.array([get_values_intent(wav, E2, V2, H, T, zb, rhs_2S, lhs_2S) for wav in wavs])

#Returns individual vectors of for Rayleigh and total Raman scattering for each plot.
ralrates_1S = np.array([rates_1S[i, 0]*100 for i in range(0, len(wavs))])
ramrates_1S = np.array([rates_1S[i, 1]*100 for i in range(0, len(wavs))])

ralrates_2S = np.array([rates_2S[i, 0]*100 for i in range(0, len(wavs))])
ramrates_2S = np.array([rates_2S[i, 1]*100 for i in range(0, len(wavs))])

ralrates_intent = np.array([rates_i[i, 0] for i in range(0, len(wavs))])
ramrates_intent = np.array([rates_i[i, 1] for i in range(0, len(wavs))])


#Sets up the figure
Fig, (ax_intent, ax_2S, ax_1S) = plt.subplots(3, sharex=True, figsize=(9,9), dpi = 100, tight_layout =True)

#General p energy spectrum, needed for resonances.
Es = imp.Schrodinger((0, 1), nmax, k)[0]


# Used to add ionisation rates, initially calculated for an intensity of 1E8 W/cm^2 This is used to get the ionisation rates at specific depths or intensities
df = pd.read_excel('Ionisation_data.xlsx')
wavs_ion = [w for w in df.Wav]
rates_ion = [r for r in df.non_pert]
S1_rates_ion = [(100**2)*S1rec_depth(wavs_ion[i])*rates_ion[i] for i in range(0, len(rates_ion))]
S2_rates_ion = [(100**2)*rec_depth(wavs_ion[i])*rates_ion[i] for i in range(0, len(rates_ion))]
intent_rates_ion = [(((3.37E8)/((1E8)))**2)*rates_ion[i] for i in range(0, len(rates_ion))]

#General plotting functions
ax_1S.plot(wavs, ralrates_1S, linestyle="dotted", label="Elastic", linewidth=2)
ax_1S.plot(wavs, ramrates_1S, label="Inelastic", linewidth=2)
ax_1S.plot(wavs_ion, S1_rates_ion, linestyle="dashed", label="Ionisation", linewidth = 2)

ax_2S.plot(wavs, ralrates_2S, linestyle="dotted", label="Elastic", linewidth=2)
ax_2S.plot(wavs, ramrates_2S, label="Inelastic", linewidth=2)
ax_2S.plot(wavs_ion, S2_rates_ion, linestyle="dashed", label="Ionisation", linewidth = 2)

ax_intent.plot(wavs, ralrates_intent, linestyle="dotted", label="Elastic", linewidth=2)
ax_intent.plot(wavs, ramrates_intent, label="Inelastic", linewidth=2)
ax_intent.plot(wavs_ion, intent_rates_ion, linestyle="dashed", label="Ionisation", linewidth = 2)

for i in range(3, 8): # plots lines relating to resonances in the 2S state, changing the plotting range may requrie this to be updated 
    line = imp.h*imp.c/(imp.Eh*abs(E2 - Es[i - 2]))*10**9
    ax_1S.plot([line, line], [1, 10**8],  color="black", linewidth=1)
    ax_2S.plot([line, line], [1, 10000],  color="black", linewidth=1)
    ax_intent.plot([line, line], [1, 10**8],  color="black", linewidth=1)

magicx = [514.646, 443.212, 414.483, 399.451] # plots lines relating to the known 1S-2S magic wavelengths, changing the plotting range may requrie this to be updated 
for i in range(0, 4):
    ax_1S.plot([magicx[i], magicx[i]], [1, 10**8],  color ="#d62728", linewidth = 1)
    ax_2S.plot([magicx[i], magicx[i]], [1, 10000],  color ="#d62728", linewidth = 1)
    ax_intent.plot([magicx[i], magicx[i]], [1, 10**8] , color ="#d62728", linewidth = 1)
    
ax_1S.set_yscale("log")
ax_2S.set_yscale("log")
ax_intent.set_yscale("log")

ax_1S.set_xlabel(r'$\lambda$ [nm]')

ax_1S.set_ylabel(r'$R [s^{-1}]$')
ax_2S.set_ylabel(r'$R [s^{-1}]$')
ax_intent.set_ylabel(r'$R [s^{-1}]$')


ax_1S.set_ylim(10, 10**8) # neccessary due to the presence of poles and points of 0 in the rate structure
ax_2S.set_ylim(1, 10000)
ax_intent.set_ylim(10, 10**8)

ax_intent.set_xlim(395, 1000) # these should ideally match the range of "wavs"
ax_intent.legend()

plt.show()
