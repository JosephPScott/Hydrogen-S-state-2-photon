# -*- coding: utf-8 -*-
"""
Last edited on  Monday 17 July 2023

@author: joseph.p.scott@durham.ac.uk

A series of functions used to calculate certain second-order dipole effects in S-states of atomic hydrogen.
Including the A.C. polarisability, calculations of 1S-2S magic wavelengths, and off resonant scattering rates.
Calculation proceeds via implicit summation over a basis of radial Sturmian functions (see write up), angular parts are calculated analytically.
General user inputs to all calculations are the number of radial Sturmian functions, nmax, and the value of the free parameter, k.
Users should ensure that a sufficiently large basis is used to ensure convergence, nmax ~ 100's is usually sufficient for low lying states and small k.
"""

###############################################################################
# Import packages
import scipy.linalg as al
import scipy.optimize as opt
import scipy.misc as mic
import numpy as np

###############################################################################
#Define constants
fsc = 7.2973525693E-3 #Fine structure constant

#Atomic units
h = 6.62607015E-34 #plank constant
hbar = h/(2*np.pi) #reduced plank constant
e = 1.602176634E-19 #elementary charge
a0 = 5.29177210903E-11 #Bohr radius
me = 9.10938370E-31 #electron rest mass

Eh = (hbar**2)/(me*(a0**2)) #Hartree energy

#Additional constants
c = 299792458 #speed of light
ep0 = 8.854188E-12 #Vacuume permativity
Ry = 10973731.568160 #Rydberg constant
mp = 1.67262192369E-27 #Proton rest mass

#Reduced mass
redm = mp/(mp + me) #for hydrogen
#redm = (mp+mn)/(mp + mn + me) #for deuterium
#redm = (mp+2*mn)/(mp + 2*mn + me) #for tririum

###############################################################################
#Construct matrix representations of the operators in a basis of radial sturmian functions, require specified orbital angular momentum l

def Hamiltonian(nmax, l, k): # Representing the atomic Hamiltonian
    return np.array([[Hamelm(j, i, l, k) for i in range(1, nmax+1)]for j in range(1, nmax+1)])

def Hamelm(n1, n2, l, k): # Hamiltonian matrix elements
    dict={
        -1: 0.25*(k/redm)*np.sqrt((n2 - 1)*(n2 + 2*l)),
        0: (0.5*(k/redm)*(n2+l) - 1),
        1: 0.25*(k/redm)*np.sqrt(n2*(n2 + 2*l +1))
        }
    return dict.get(n1 - n2, 0)


def Overlap(nmax, l, k): # Representing the overlap matrix
    return np.array([[Overlapelm(i, j, l, k)/k for i in range(1, nmax+1)]for j in range(1, nmax+1)])

def Overlapelm(n1, n2, l, k): # Overlap matrix elements
    dict={
        -1: -0.5*np.sqrt(n2*(n2+2*l + 1)),
        0: (n2 + l),
        1: -0.5*np.sqrt((n2 - 1)*(n2 + 2*l))
        }
    return dict.get(n2 - n1, 0)


#The dipole matrix element can be broken into two parts, l raising...
def Zplus(nmax, l, k): # Representing the l raising part
    if l >= 0:
        angular = 1
        return np.array([[Zpluselm(i, j, l, k)*(angular/(4*k**2)) for i in range(1, nmax + 1)]for j in range(1, nmax + 1)])
    else:
        return np.zeros((nmax, nmax))
    
def Zpluselm(n2, n1, l, k): # l-raising matrix elements
    dict = {
        -1:-np.sqrt((n1 + 2*l + 2)*(n1 + 2*l + 1)*(n1 + 2*l)*(n1 -1)),
        0: 2*(2*n1 + l)*np.sqrt((n1 + 2*l + 2)*(n1 + 2*l + 1)),
        1: -6*(n1 + l + 1)*np.sqrt(n1*(n1 + 2*l + 2)),
        2: 2*(2*n1 + 3*l +4)*np.sqrt(n1*(n1 + 1)),
        3: -np.sqrt((n1 + 2*l + 3)*(n1 + 2)*(n1 + 1)*n1)
            }
    return dict.get(n2 - n1, 0)

#... and l lowering.
def Zminus(nmax, l, k): # Representing the l lowering part
    if l > 0:
        angular = 1
        return np.array([[Zminuselm(i, j, l, k)*(angular/(4*k**2)) for i in range(1, nmax+1)]for j in range(1, nmax+1)])
    else:
        return np.zeros((nmax, nmax))
    
def Zminuselm(n2, n1, l, k):  # l-lowering matrix elements
    dict = {
        -1:-np.sqrt((n2 + 2*l - 2)*(n2 + 2*l - 1)*(n2 + 2*l)*(n2 -1)),
        0: 2*(2*n2 + l - 1)*np.sqrt((n2 + 2*l)*(n2 + 2*l - 1)),
        1: -6*(n2 + l)*np.sqrt(n2*(n2 + 2*l)),
        2: 2*(2*n2 + 3*l + 1)*np.sqrt(n2*(n2 + 1)),
        3: -np.sqrt((n2 + 2*l + 1)*(n2 + 2)*(n2 + 1)*n2)
            }
    return dict.get(n1 - n2, 0)

###############################################################################
#General calculation functions

def Schrodinger(st, nmax, k): # Solves the Schrodinger equation, st = (n, l), if the full spectrum is desired, input (0, l)
    H = Hamiltonian(nmax, st[1], k)
    T = Overlap(nmax, st[1], k)
    w, v = al.eigh(H, T, type = 1)
    E, Vs = sorted(w.real), v[:, w.argsort()]
    if st[0] > 0:
        return E[(st[0] - st[1]) -1], Vs[:, (st[0] - st[1]) - 1] # Returns energy and vector representation for a given n.
    else:
        return E, Vs  # If n is given as 0, returns the entire spectrum


def Implicit_step(H, T, lhs, rhs, E, freq): # Calculates the implicit summation step via the relevant matrix element.
    H_adj = H -(E + freq)*T
    Psi = al.solve(H_adj, rhs)
    return np.dot(lhs, Psi)

###############################################################################
#S-state polarisability calculations

def S_pol(n, wave, nmax, k): # Calculates the polarisability of a given s-state
    E, V = Schrodinger((n, 0), nmax, k)
    za = Zplus(nmax, 0, k) # only a p intermediate state so l=1
    H = Hamiltonian(nmax, 1, k)
    T = Overlap(nmax, 1, k)
    zb = Zminus(nmax, 1, k)
    lhs = np.dot(V, zb)
    rhs = np.dot(za, V)
    freq = (h*c)/(wave*10**(-9)*Eh)
    angular = 1/3   # This term is calculated analytically
    return angular*(Implicit_step(H, T, lhs, rhs, E, freq) + Implicit_step(H, T, lhs, rhs, E, -1*freq)) # Polarisability is given in a.u.

###############################################################################
#Functions relating to the 1S-2S magic wavelengths: For calculating the differential polarisability, finding the magic wavelength according to an initial guess and calculating the slope in differential polarisability

#Functions for the differential polarisability of the 1S and 2S states for either a "wavelength" or "frequency" input.
def Differential_1S2Spol_wave(wave, H, T, lhs1, lhs2, rhs1, rhs2, E1, E2): # Takes wavelength input in nm
    freq = (h*c)/(wave*10**(-9)*Eh)
    S1 = Implicit_step(H, T, lhs1, rhs1, E1, freq) + Implicit_step(H, T, lhs1, rhs1, E1, -1*freq)
    S2 = Implicit_step(H, T, lhs2, rhs2, E2, freq) + Implicit_step(H, T, lhs2, rhs2, E2, -1*freq)
    return (S2 - S1)/3

def Differential_1S2Spol(freq, H, T, lhs1, lhs2, rhs1, rhs2, E1, E2): # Takes frequency input in a.u.
    S1 = Implicit_step(H, T, lhs1, rhs1, E1, freq) + Implicit_step(H, T, lhs1, rhs1, E1, -1*freq)
    S2 = Implicit_step(H, T, lhs2, rhs2, E2, freq) + Implicit_step(H, T, lhs2, rhs2, E2, -1*freq)
    return (S2 - S1)/3

def Find_1S2S_magicwave(nmax, wguess, k): # Calculate the magic wavelength from a given initial "guess" wavelength (in nm) by newton raphson
    E1, V1 = Schrodinger((1, 0), nmax, k)
    E2, V2 = Schrodinger((2, 0), nmax, k)
    za = Zplus(nmax, 0, k)
    H = Hamiltonian(nmax, 1, k)
    T = Overlap(nmax, 1, k)
    zb = Zminus(nmax, 1, k)
    lhs1 = np.dot(V1, zb)
    lhs2 = np.dot(V2, zb)
    rhs1 = np.dot(za, V1)    
    rhs2 = np.dot(za, V2)
    fguess = (h*c)/(wguess*10**(-9)*Eh)
    root = opt.newton(Differential_1S2Spol, fguess, args = (H, T, lhs1, lhs2, rhs1, rhs2, E1, E2)) #Newton-Raphson solve to find minimum in differential polarisability
    return (h*c*10**9)/(root*Eh) # Returns a wavelength in nm

def Magic_1S2S_stability(nmax, wave, k): # Calculate the slope on the differential 1S-2S light shift at a given wavelength
    E1, V1 = Schrodinger((1, 0), nmax, k)
    E2, V2 = Schrodinger((2, 0), nmax, k)
    za = Zplus(nmax, 0, k)
    H = Hamiltonian(nmax, 1, k)
    T = Overlap(nmax, 1, k)
    zb = Zminus(nmax, 1, k)
    lhs1 = np.dot(V1, zb)
    lhs2 = np.dot(V2, zb)
    rhs1 = np.dot(za, V1)    
    rhs2 = np.dot(za, V2)
    return mic.derivative(Differential_1S2Spol_wave, wave, dx = 0.00001, n=1, args=(H, T, lhs1, lhs2, rhs1, rhs2, E1, E2), order=5) # Returns in atomic units of polarisability per nm

#Then the same as above, but this time for general nS-n'S transitions.
def Differential_Spol_wave(wave, H, T, lhs1, lhs2, rhs1, rhs2, E1, E2): # Takes wavelength input in nm
    freq = (h*c)/(wave*10**(-9)*Eh)
    S1 = Implicit_step(H, T, lhs1, rhs1, E1, freq) + Implicit_step(H, T, lhs1, rhs1, E1, -1*freq)
    S2 = Implicit_step(H, T, lhs2, rhs2, E2, freq) + Implicit_step(H, T, lhs2, rhs2, E2, -1*freq)
    return (S2 - S1)/3

def Differential_Spol(freq, H, T, lhs1, lhs2, rhs1, rhs2, E1, E2):# Takes frequency input in a.u.
    S1 = Implicit_step(H, T, lhs1, rhs1, E1, freq) + Implicit_step(H, T, lhs1, rhs1, E1, -1*freq)
    S2 = Implicit_step(H, T, lhs2, rhs2, E2, freq) + Implicit_step(H, T, lhs2, rhs2, E2, -1*freq)
    return (S2 - S1)/3

def Find_magicwave(n1, n2, wguess, nmax, k): # Calculate the magic wavelength from a given initial "guess" wavelength (in nm) by newton raphson
    E1, V1 = Schrodinger((n1, 0), nmax, k)
    E2, V2 = Schrodinger((n2, 0), nmax, k)
    za = Zplus(nmax, 0, k)
    H = Hamiltonian(nmax, 1, k)
    T = Overlap(nmax, 1, k)
    zb = Zminus(nmax, 1, k)
    lhs1 = np.dot(V1, zb)
    lhs2 = np.dot(V2, zb)
    rhs1 = np.dot(za, V1)    
    rhs2 = np.dot(za, V2)
    fguess = (h*c)/(wguess*10**(-9)*Eh)
    root = opt.newton(Differential_Spol, fguess, args = (H, T, lhs1, lhs2, rhs1, rhs2, E1, E2)) #Newton-Raphson solve to find minimum in differential polarisability
    return (h*c*10**9)/(root*Eh) # Returns a wavelength in nm
    
def Magic_stability(n1, n2, nmax, wave, k): # Calculate the slope on the differential 1S-2S light shift at a given wavelength
    E1, V1 = Schrodinger((n1, 0), nmax, k)
    E2, V2 = Schrodinger((n2, 0), nmax, k)
    za = Zplus(nmax, 0, k)
    H = Hamiltonian(nmax, 1, k)
    T = Overlap(nmax, 1, k)
    zb = Zminus(nmax, 1, k)
    lhs1 = np.dot(V1, zb)
    lhs2 = np.dot(V2, zb)
    rhs1 = np.dot(za, V1)    
    rhs2 = np.dot(za, V2)
    return mic.derivative(Differential_Spol_wave, wave, dx = 0.00001, n=1, args=(H, T, lhs1, lhs2, rhs1, rhs2, E1, E2), order=5) # Returns in atomic units of polarisability per nm


###############################################################################
#Functions for calculating atom-photon scattering rates

# Rayleigh scattering rates for hydrogen S states 
# Rates are calculated for a specified wavelength and either a well defined intensity (input in SI) or at a given unitless depth - calculated for the specific (n, 0) state that was given as an input n.
def S_Rayleigh(n, wave, nmax, k, Inten): #Calculates the rayleigh scattering of a given S state at specified wavelength and intensity
    freq = (h*c)/(wave*10**(-9)*Eh)
    I = (Inten*hbar*a0**2)/Eh**2
    return (freq**3)*(fsc**4)*(8*np.pi/3)*(abs(S_pol(n, wave, nmax, k))**2)*I*(Eh/hbar) # Rate is given in per second

def S_Rayleigh_Depth(n, wave, nmax, k, D): #Same as the above, but calculates in terms of lattice depth (in the target state) rather than intensity
    Erec = (h**2)/(2*(mp+me)*(wave*10**(-9))**2*Eh)
    freq = (h*c)/(wave*10**(-9)*Eh)
    return (freq**3)*(fsc**3)*(4/3)*(abs(S_pol(n, wave, nmax, k)))*D*Erec*(Eh/hbar) # Rate is given in per second 

# Calculates the inelastic scattering rate out of an intial S state to specifed S or D state (functions are named "Raman" but allow for SSTPE as well)
# Rates are calculated for a specific wavelength and either a well defined intensity (input in SI) or at a given unitless depth - calculated for the specific (n, 0) state that was given as an input n.
def S_Raman(n, final_st, wave, nmax, k, Inten): #Calculates scattering rate for a given "n"S state and specified final state final_st = (n, l) at a given intensity "Inten" in S.I. units.
    freq = (h*c)/(wave*10**(-9)*Eh)
    #Represent atomic states
    E_int, V_int = Schrodinger((n, 0), nmax, k)
    E_fin, V_fin = Schrodinger(final_st, nmax, k) 
    # Check for energetically forbidden inputs
    freq_scatt = -(E_fin - E_int) + freq
    if  freq_scatt <= 0:
        return 0, 0
    I = (Inten*hbar*a0**2)/Eh**2
    za = Zplus(nmax, 0, k) # Only P intermetidate states, so l=1
    H = Hamiltonian(nmax, 1, k)
    T = Overlap(nmax, 1, k)
    # The angular components have been calculated analytically depending on if the final state is S or D
    if final_st[1] == 0: # In the case where the final state is an S state
        angular = 1/9
        zb = Zminus(nmax, 1, k)
    elif final_st[1] == 2: # In the case where the final state is a D state
        angular = 2/9
        zb = Zplus(nmax, 1, k)
    else: # Handles requests for dipole forbidden final states
        print("Forbidden")
        return (0, 0)
    rhs = np.dot(za, V_int)
    lhs = np.dot(V_fin, zb)
    #There are two types of scattering processes that can occur
    P1 = (freq_scatt**3)*(fsc**4)*(8*np.pi/3)*abs(Implicit_step(H, T, lhs, rhs, E_int, freq) + Implicit_step(H, T, lhs, rhs, E_fin, -1*freq))**2*I #Relates to Stokes/Anti-Stokes Raman scattering
    freq_scatt = -(E_fin - E_int) - freq
    if freq_scatt <= 0: # Conservation of energy constrains the scattering processes, can be seen as a positive definite scattered frequency
        return angular*P1*(Eh/hbar), 0 # When SSTPE is an energetically forbidden processs
    else:
        P2 = (freq_scatt**3)*(fsc**4)*(8*np.pi/3)*abs(Implicit_step(H, T, lhs, rhs, E_int, -1*freq) + Implicit_step(H, T, lhs, rhs, E_fin, freq))**2*I #Relates to SSTPE (when energetically allowed)
        return angular*P1*(Eh/hbar), angular*P2*(Eh/hbar) #Outputs a tuple of both rates (Stokes/anti-stokes, SSTPE) in units of per second


def Raman(H, T, lhs, rhs, Ea, Eb, freq, I): # The internal calculation step, called by the other functions
    freq_scatt = -(Eb - Ea) + freq
    if freq_scatt <= 0: # Conservation of energy constrains the scattering processes, can be seen as a positive definite scattered frequency
        return 0 # For energetically disallowed processes
    else:
        return (freq_scatt**3)*(fsc**4)*(8*np.pi/3)*abs(Implicit_step(H, T, lhs, rhs, Ea, freq) + Implicit_step(H, T, lhs, rhs, Eb, -1*freq))**2*I # For energetically allowed processes


def S_Raman_Depth(n, final_st, wave, nmax, k, D): #Calculates scattering rate for a given "n"S state and specified final state final_st = (n, l) at a given unitless depth (measured in the initial nS state) "D".
    freq = (h*c)/(wave*10**(-9)*Eh)
    #Represent atomic states
    E_int, V_int = Schrodinger((n, 0), nmax, k)
    E_fin, V_fin = Schrodinger(final_st, nmax, k)
    # Check for energetically forbidden inputs
    freq_scatt = -(E_fin - E_int) + freq
    if freq_scatt  <= 0:
        return 0, 0
    rhs = np.dot(Zplus(nmax, 0, k), V_int) # Only P intermetidate states, so l=1
    H = Hamiltonian(nmax, 1, k) 
    T = Overlap(nmax, 1, k)
    zb = Zminus(nmax, 1, k)
    lhs_pol = np.dot(V_int, zb)
    # The polarisability of the state is necessary to obtain it's depth.
    pol = abs((Implicit_step(H, T, lhs_pol, rhs, E_int, freq) + Implicit_step(H, T, lhs_pol, rhs, E_int, -1*freq))/3)
    # The angular components have been calculated analytically depending on if the final state is S or D
    if final_st[1] == 0: # In the case where the final state is an S state
        angular = 1/9
        lhs = np.dot(V_fin, zb)
    elif final_st[1] == 2: # In the case where the final state is a D state
        angular = 2/9
        lhs = np.dot(V_fin, Zplus(nmax, 1, k))
    else: # Handles requests for dipole forbidden final states
        print("Forbidden")
        return (0, 0)
    P1 = (2*D*(freq**2)*(freq_scatt**3)*(fsc**5)*abs(Implicit_step(H, T, lhs, rhs, E_int, freq) + Implicit_step(H, T, lhs, rhs, E_fin, -1*freq))**2)/(3*pol*((mp/me) + 1)) #Relates to Stokes/Anti-Stokes Raman scattering
    freq_scatt = -(E_fin - E_int) - freq
    if freq_scatt <= 0:  # Conservation of energy constrains the scattering processes, can be seen as a positive definite scattered frequency
        return angular*P1*(Eh/hbar), 0 # When SSTPE is an energetically forbidden processs
    else:
        P2 = (2*D*(freq**2)*(freq_scatt**3)*(fsc**5)*abs(Implicit_step(H, T, lhs, rhs, E_int, -1*freq) + Implicit_step(H, T, lhs, rhs, E_fin, freq))**2)/(3*pol*((mp/me) + 1)) #Relates to SSTPE (when energetically allowed)
        return angular*P1*(Eh/hbar), angular*P2*(Eh/hbar) #Outputs a tuple of both rates (Stokes/anti-stokes, SSTPE) in units of per second

