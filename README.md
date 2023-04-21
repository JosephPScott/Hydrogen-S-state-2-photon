# Hydrogen-S-state-2-photon
These modules contain functions relating to the calculations presented in [paper](). 
Specifically, the calculations of A.C. polarisability and atom-photon scattering rates in linearly polarised, off-resonance optical fields for S states in hydrogen. Calculations proceed via implicit summation over a basis of radial Sturmian functions, angular components take their analytical values.  The mathematical details of these calculations are described in the accompanying [supplementary material]().

## Installation

Modules are available to download directly and can be run in the users preffered environment.

## Usagae

### Core module
There is a single core module "Main_calculation_functions". This contains all of the functions related directly to calculation and is a requirement for the use of the additional modules. There are two parametes that are common to all functions in this module: "nmax" and "k". These parameters define the properties of the basis of radial Sturmian functions: "nmax" is the number of Sturmian functions and "k" is a free parameter. Together these are called the "basis parameters". 
The user should ensure that nmax is sufficiently large to ensure convergnence of calculated values. Whilst k is technically free, choosing k=0 will diverge and it takes convergence takes longer for large k. For low lying states, nmax ~ 100's and k~0.1 will usually be sufficient.

This module may be used for the following calculations:
- A.C. polarisability
```python
S_pol(n, wave, nmax, k)
``` 
Takes in the principal quantum number of the reference s state "n", the wavelegnth of the optical field in nm "wave", and the two basis parameters. Outputs the A.C. polarisability of the state at this wavelength in atomic units.
- Magic wavelengths for the 1s-2s transition
```python
Find_1S2S_magicwave(nmax, wguess, k)
Magic_1S2S_stability(nmax, wave, k)
```
The first function takes in a reasonable guess to a magic wavelength in nm "wguess" and the basis parameters. It outputs the nearest magic wavelength in nm.

The second function takes in the wavelength "wave" and the basis parameters. It outputs the "slope" of the 1s-2s differential A.C. polarisability at this wavelength in atomic units of polarisability per nm.
- Rayleigh scattering
```python
S_Rayleigh(n, wave, nmax, k, Inten)
S_Rayleigh_Depth(n, wave, nmax, k, D)
```
The first function takes in the principle quantum number of the state "n", the wavelength of the optical field in nm "wave", the intensity of the optical field in S.I. units "Inten", and the basis parameters. It outputs the Rayleigh scattering rate in units of per second. 

The second function is designed for calculations in optical lattices. It does the same as the first but instead of accepting an intensity, it takes in a lattice depth in recoil units.

- Raman scattering from an s state to final s or d state:
```python
S_Raman(n, final_st, wave, nmax, k, Inten)
S_Raman_Depth(n, final_st, wave, nmax, k, D)
```
The first function takes in the principle quantum number of the initial state "n", the final state "final_st" as a tuple (principle quantum number, orbital angular momentum number), the wavelength of the optical field in nm "wave", the intensity of the optical field in S.I. units "Inten", and the basis parameters. It outputs the Raman scattering rate in units of per second as a tuple, relating to the two possible scattering processes. If the scattering process between the specified intial and final state is energetically or dipole forbidden then the fnction returns 0 for this rate.

The second function is designed for calculations in optical lattices. It does the same as the first but instead of accepting an intensity, it takes in a lattice depth in recoil units.

### Additional modules
These modules are additional to the main module and are not strictly necessary to the general user, all depend upon the main module to operate. The main module is usually called with "as imp" to make usage simpler. The additional modules are as follows:

- "Convergence_and_correctness_tests.py"
Contains functions used to check the results of calculations made in the core module. The results of the schrodinger equation (the energy spectrum and the vector representation of the wavefunctions) are directly compared to the analytic results:
```python
Check_EnergySpectrum(l, nmax, k)
Check_Wavefunctions(st, nmax, k)
```
The module can also produce a number of heatmaps which can be used to assess the variation of the returned value or computational error as the basis parameters change:
```python
Spectrum_error_map()
Dipole_stability(initial_st, final_st)
Polarisabilty_stability(n, wave)
Magic_wavelength_stability(wguess)
Raman_scattering_stability(n, state, wave)
```
For the energy specturum, this is the average error compared to the Rydberg theory over a finite number of discrete states whilst for the others it is compared to a value calculated for nmax=300, k=0.3.

- "Continuum_corrections.py"
Produces a single plot to compare the Raman scattering rates calculated via implicit summation compared to that computed with by summing across a finite number of discrete states. This is used to invesitgate the importance of the continuum states to calculation. The plot gives the ratio of these results as we increase the number of discrete states and plots a large n limit. It currently runs on a single function:
```python
Run_plot()
```
Altering the parameters of this plot can be done directlt within the module.
**Note: this implemented in an inefficient way**

- "1S_and_2S_polarisability.py"
This module generates a comparative plot of the atomic polarisabilities of the 1S and 2S state across some spectral range, it also identifies and marks the 2S resonances in this region. It will generate the plot upon running the module and contains no governing function. Altering the parameters of this plot can be done directlt within the module.

- "Compare_scattering_rates.py"
This module generates a comparative plot of the Rayleigh and total Raman scattering rates of the 2s state across some spectral range, it also identifies and marks the 2S resonances in this region. It will generate the plot upon running the module and contains no governing function. Altering the parameters of this plot can be done directlt within the module.
**Note: this can take significant time to compute due to the design of the scattering rate functions in "Main_calculation_functions.py" and the wavelength resolution required**

## Dependencies
- "Main_calculation_functions.py"(imp) - numpy, scipy.
- "Convergence_and_correctness_tests.py" - imp, numpy, matplotlib, seaborn, pandas.
- "Continuum_corrections.py" - imp, numpy, scipy, matplotlib.
- "1S_and_2S_polarsiability.py" - imp, numpy, matplotlib.
- "Compare_scattering_rates.py" - imp, numpy, matplotlib.

## Scope for improvement

Currently, calculations are limited to s states in atomic hydrogen and linearly polarised optical fields. Additionally, it operates in the non-relativistic theory, however, relativisitc and QED corrections to these calculations are not usually large in hydrogen - especially for the S states. The most obvious extensions are to a general polarisation of light and initial states with different orbital angular momentum. 

Efficiency can also be improved upon - particularly when plotting. At the moment, calling the full calculation from "Main_calculation_functions.py" for plots of polarisation or scattering rates leads to many redundnat calculations. The calculation of the discrete states in "Convergence_and_correctness_tests.py" could also stand significant improvement.

## Authors
Joseph P. Scott

Department of Physics, Durham University, Durham, UK

joseph.p.scott@durham.ac.uk

## Date
Fri Apr 21 2023
