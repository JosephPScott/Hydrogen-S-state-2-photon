# Hydrogen-S-state-2-photon
These modules contain functions relating to the calculations and figures presented in ["Trap induced broadening in a potential hydrogen lattice clock"](https://arxiv.org/abs/2309.03753). 

The code covers the calculation of A.C. polarisability (including functions for finding S to S magic wavelenegths), and off-resonance atom-photon scattering rates of hydrogen S states in a linearly polarised, monochromatic, optical field. These terms are separated into radial and angular components: the radial parts are calculated via implicit summation over a basis of radial Sturmian functions, while angular parts are calculated analytically (see the appendices of [the paper](https://arxiv.org/abs/2309.03753) for details).

In addition to these calculations, this code contains a number of functions used to produce the plots seen in [the paper](https://arxiv.org/abs/2309.03753) and check the stability of calculations.

## Installation

Modules are available to download directly by a user. The most important module is "Main_calculaltion_functions.py", the others are subsiduary.

## General Overview

### Core module
There is a single core module "Main_calculation_functions". This contains all of the functions related directly to calculation and is a requirement for the use of the additional modules. There are two parametes that are common to all functions in this module: "nmax" and "k". These parameters define the properties of the basis of radial Sturmian functions: "nmax" is the number of Sturmian functions and "k" is a free parameter. We will call these the "basis parameters". 
The user should ensure that nmax is sufficiently large to ensure convergnence of calculated values. Whilst k is technically free, choosing k=0 will diverge and it takes convergence takes longer for large k. For low lying states, nmax ~ 100's and k ~ 0.1 will be sufficient.

This module contains functions for the following calculations:
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
The first function takes in the principle quantum number of the initial state "n", the wavelength of the optical field in nm "wave", the intensity of the optical field in S.I. units "Inten", and the basis parameters. It outputs the Rayleigh scattering rate in units of per second. 

The second function is designed for calculations in optical lattices. It does the same as the first but instead of accepting an intensity, it takes in a lattice depth (of the chosen initial state) in recoil units.

- Raman scattering from an S state to final S or D state:
```python
S_Raman(n, final_st, wave, nmax, k, Inten)
S_Raman_Depth(n, final_st, wave, nmax, k, D)
```
The first function takes in the principle quantum number of the initial state "n", the final state "final_st" as a tuple (principle quantum number, orbital angular momentum number), the wavelength of the optical field in nm "wave", the intensity of the optical field in S.I. units "Inten", and the basis parameters. It outputs the Raman scattering rate in units of per second as a tuple, relating to the two possible scattering processes. If the scattering process between the specified intial and final state is energetically or dipole forbidden then the fnction returns 0 for this rate.

The second function is designed for calculations in optical lattices. It does the same as the first but instead of accepting an intensity, it takes in a lattice depth (of the chosen initial state) in recoil units.

### Checking numerical accuracy
The file "Convergence_and_correctness.py" contains a number of functions that are useful for assessing the quality of these calculations. The accuracy of the calculations can be assessed by compating their results to those with known analytic solutions, e.g. comparing the spectrum of the Hamiltonian to the non-relativistic sstructure of hydrogen, or comparing calculated values of Raman scattering cross sections to literature values. The stability of these calculations is assessed by looking at the variation of the calculated value under small changes in the basis parameters nmax and k.

Functions are available to check the following:
- Calculate Raman scattering cross sections for comparison to existing literature such as [Klarsfeld](https://doi.org/10.1103/PhysRevA.6.506) or [Heno _et. al._](https://pubs.aip.org/aip/jap/article/51/1/11/12151/Raman-like-scattering-processes-in-metastable).
```python
obtain_cross_section(n, nmax, k)
```
- Check the correctness of the Hamiltonian energy spectrum and eigen-vector forms to analytical values of the hydrogen energy spectrum and wavefunctions.
```python
Check_EnergySpectrum(l, nmax, k)
Spectrum_error_map()
Check_Wavefunctions(st, nmax, k)
```
- Calculating the individual dipole matrix elements between specified hydogen states for comparison to analytical theory and a heatmap of the stability of such a value.
```python
Calculate_DipoleElm(initial, final, nmax, k)
Dipole_stability(initial_st, final_st)
```
- Heatmaps for the stability of the key calculated values - polarisability, magic wavelength, and Raman scattering rates.
```python
Polarisabilty_stability(n, wave)
Magic_wavelength_stability(wguess)
Raman_scattering_stability(n, state, wave)
```

### Plotting functions
These modules are used to generate the plots seen in ["Trap induced broadening in a potential hydrogen lattice clock"](https://arxiv.org/abs/2309.03753).

- "1S_and_2S_polarisability.py"
Relating to figure 1, this module generates a comparative plot of the atomic polarisabilities of the 1S and 2S state across some spectral range. It also identifies and marks the 2S resonances in this region. It will generate the plot upon running the module and contains no governing function. Altering the parameters of this plot can be done directlt within the module.

- "Variations_with_depth.py"
Relating to figure 2, this module plots the variation of the inelastic scattering rate and two-photon ionisation rate with increasing lattice at a chosen magic wavelength. Unlike the other plotting modules, it does not rely on data from other files in this repository, but uses specified values that have been calculated prior. The parameters of the plot can only be edited directly within the code.

- "Compare_scattering_rates.py"
Relating to figure 3, this module generates a comparative plot of the Rayleigh and total Raman scattering rates of the 2s state across some spectral range. It also marks the 2S resonances and 1S--2S magic wavelengths in this region. Three plots are generated, one with a constant intensity (100 Mw/cm^2), one at a constant depth in the 2S state, and one at a constant depth in the 1S state (see [the paper](https://arxiv.org/abs/2309.03753) for more detail). It will generate the plot upon running the module and contains no governing function. Altering the parameters of this plot can be done directlt within the module.
**Note: this can take a while to run and so contains code to monitor time and memory consumption. These can be activated by uncommenting them.**

- "Continuum_corrections.py"
Relating to figure 4, this module produces a single plot to compare the Raman scattering rates calculated via implicit summation compared to that computed with by summing across a finite number of discrete states. The plot gives the ratio of these results as we increase the number of discrete states and plots a large n limit. It currently runs on a single function:
```python
Run_plot()
```
Altering the parameters of this plot can be done directly within the module.

## Additional files
Also included are two xlsx data files. These contain the results of calculations done in the [STRFLO]() software and are used in "Compare_scattering_rates.py".
- Ionisation_data/xlsx - 2 photon ionisation rates, calculated perturbatively at 100 MW/cm^2.
- 3photon_ionisation_data.xlsx - 3 photon ionisation rates, calculated perturbatively at 337 MW/cm^2.

## Dependencies
- "Main_calculation_functions.py"(imp) - numpy, scipy.
- "Convergence_and_correctness_tests.py" - imp, numpy, matplotlib, seaborn, pandas.
- "1S_and_2S_polarsiability.py" - imp, numpy, matplotlib.
- "Variations_with_depth.py" - imp, numpy, matplotlib.
- "Compare_scattering_rates.py" - imp, numpy, matplotlib, "Ionisation_data.xlsx", "3photon_ionisation_data.py".
- "Continuum_corrections.py" - imp, numpy, scipy, matplotlib.

## Authors
Joseph P. Scott

Department of Physics, Durham University, Durham, UK

joseph.p.scott@durham.ac.uk

## Date of last update
Monday 23rd September 2023
