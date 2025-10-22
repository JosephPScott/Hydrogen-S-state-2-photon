# Hydrogen-S-state-2-photon
These modules contain functions relating to the calculations and figures presented first in ["Trap induced broadening in a potential hydrogen lattice clock"](https://doi.org/10.1088/1681-7575/ad1e37) and later in ["Towards precision measurement with trapped hydrogen atoms"](https://etheses.dur.ac.uk/15616/). 

The code covers the calculation of A.C. polarisability (including functions for finding S to S magic wavelenegths), and off-resonance atom-photon scattering rates of hydrogen S states in a linearly polarised, monochromatic, optical field. These terms are separated into radial and angular components: the radial parts are calculated via implicit summation over a basis of radial Sturmian functions, while angular parts are calculated analytically (see the appendices of [the metrologia paper](https://doi.org/10.1088/1681-7575/ad1e37) and [the thesis](https://etheses.dur.ac.uk/15616/) for details).

In addition to these calculations, this code contains a number of functions used to produce the specific plots seen in [the metrologia paper](https://doi.org/10.1088/1681-7575/ad1e37) and check the stability of calculations.

## Installation

Modules are available to download directly by a user. The most important module is "Main_calculaltion_functions.py", the others are secondary.

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

### Example inputs and outputs
There follows some example inputs and the expected outputs for the functions listed above. Note that when testing the results of [the metrologia paper](https://doi.org/10.1088/1681-7575/ad1e37) large gradients in the 2S polarisability near the magic wavelengths often mean that a higher precision value of these wavelengths should be used than is reported in the paper - this is discussed in [the thesis](https://etheses.dur.ac.uk/15616/).

#### Polarisability
```python
S_pol(n, wave, nmax, k)
``` 
<i>Inputs:</i> n = 1 , wave = 514.646, nmax = 300, k = 0.3  
<i>Outputs:</i> np.float64(4.727501984552037)

<i>Inputs:</i> n = 2 , wave = 514.646, nmax = 300, k = 0.3  
<i>Outputs:</i> np.float64(4.729787274477928)

<i>Inputs:</i> n = 2 , wave = 443.212, nmax = 300, k = 0.3  
<i>Outputs:</i> np.float64(4.813163820601365)

#### Magic wavelengths
```python
Find_1S2S_magicwave(nmax, wguess, k)
```
<i>Inputs:</i> nmax = 300, wguess = 510, k = 0.3  
<i>Outputs:</i> np.float64(514.646438310913)

<i>Inputs:</i> nmax = 300, wguess = 440, k = 0.3  
<i>Outputs:</i> np.float64(443.2122563829302)

```python
Magic_1S2S_stability(nmax, wave, k)
```
<i>Inputs:</i> nmax = 300, wave = 414.483, k = 0.3  
<i>Outputs:</i> np.float64(-24.885824335706936)

<i>Inputs:</i> nmax = 300, wave = 399.451, k = 0.3  
<i>Outputs:</i> np.float64(-42.09359271479248)

#### Rayleigh scattering rates
```python
S_Rayleigh(n, wave, nmax, k, Inten)
```
<i>Inputs:</i> n = 2, wave = 514.646, nmax = 300, k = 0.3, Inten = 3.372*10**10  
<i>Outputs:</i> np.float64(0.007987541387702496)

<i>Inputs:</i> n = 1, wave = 443.212, nmax = 300, k = 0.3, Inten = 4.460*10**10  
<i>Outputs:</i> np.float64(0.017105383822477796)

```python
S_Rayleigh_Depth(n, wave, nmax, k, D)
```
<i>Inputs:</i> n = 2, wave = 414.483, nmax = 300, k = 0.3, D = 1  
<i>Outputs:</i> np.float64(0.024233633316437422)

<i>Inputs:</i> n = 2, wave = 399.451, nmax = 300, k = 0.3, D = 1  
<i>Outputs:</i> np.float64(0.029289503005107607)

#### Raman scattering rates
```python
S_Raman(n, final_st, wave, nmax, k, Inten)
```
<i>Inputs:</i> n = 2, final_st = (3, 0), wave = 514.646, nmax = 300, k = 0.3, Inten = 3.372*10**10  
<i>Outputs:</i> (np.float64(3.8434991392405125), 0)

<i>Inputs:</i> n = 2, final_st = (3, 2), wave = 443.212, nmax = 300, k = 0.3, Inten = 4.460*10**10  
<i>Outputs:</i> (np.float64(0.18584871640794423), 0)

```python
S_Raman_Depth(n, final_st, wave, nmax, k, D)
```
<i>Inputs:</i> n = 2, final_st = (1, 0), wave = 414.483, nmax = 300, k = 0.3, D = 1  
<i>Outputs:</i> (np.float64(60.014366649296406), np.float64(5.735668500303974))

<i>Inputs:</i> n = 2, final_st = (4, 2) wave = 399.451, nmax = 300, k = 0.3, D = 1  
<i>Outputs:</i> (np.float64(0.46363495030601853), 0)

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
These modules are used to generate the plots seen in ["Trap induced broadening in a potential hydrogen lattice clock"](https://doi.org/10.1088/1681-7575/ad1e37). They are included both for the sake of transparency and to illustrate the potential use of the core module. Variations on these were also used to generate the plots produced in ["Towards precision measurement with trapped hydrogen atoms"](https://etheses.dur.ac.uk/15616/).

- "1S_and_2S_polarisability.py"
Relating to Figure 1, this module generates a comparative plot of the atomic polarisabilities of the 1S and 2S state across some spectral range. It also identifies and marks the 2S resonances in this region. It will generate the plot upon running the module and contains no governing function. Altering the parameters of this plot can be done directlt within the module.

- "Variations_with_depth.py"
Relating to Figure 2, this module plots the variation of the inelastic scattering rate and two-photon ionisation rate with increasing lattice at a chosen magic wavelength. Unlike the other plotting modules, it does not rely on data from other files in this repository, but uses specified values that have been calculated prior. The parameters of the plot can only be edited directly within the code.

- "Compare_scattering_rates.py"
Relating to Figure 3, this module generates a comparative plot of the Rayleigh and total Raman scattering rates of the 2s state across some spectral range. It also marks the 2S resonances and 1S--2S magic wavelengths in this region. Three plots are generated, one with a constant intensity (100 Mw/cm^2), one at a constant depth in the 2S state, and one at a constant depth in the 1S state (see [the paper](https://doi.org/10.1088/1681-7575/ad1e37) for more detail). It will generate the plot upon running the module and contains no governing function. Altering the parameters of this plot can be done directlt within the module.
**Note: this can take a while to run and so contains code to monitor time and memory consumption. These can be activated by uncommenting them.**

- "Continuum_corrections.py"
Relating to Figure 4, this module produces a single plot to compare the Raman scattering rates calculated via implicit summation compared to that computed with by summing across a finite number of discrete states. The plot gives the ratio of these results as we increase the number of discrete states and plots a large n limit. It currently runs on a single function:
```python
Run_plot()
```
Altering the parameters of this plot can be done directly within the module.

Running any of these modules will produce plots that can be directly checked against those in [the metrologia paper](https://doi.org/10.1088/1681-7575/ad1e37) (barring some additional markup).

## Additional files
Also included are two xlsx data files. These contain the results of calculations done in the [STRFLO](https://doi.org/10.1016/S0010-4655(98)00073-3) software and are used in "Compare_scattering_rates.py".
- Ionisation_data/xlsx - 2 photon ionisation rates, calculated perturbatively at 100 MW/cm^2.
- 3photon_ionisation_data.xlsx - 3 photon ionisation rates, calculated perturbatively at 337 MW/cm^2.

## Dependencies
The software was originally developed with Python 3.7, but its compatibilty has been confirmed with Python 3.13. The version numbers listed for the dependencies are the most recent for which compatibilty has been tested:

- "Main_calculation_functions.py"(imp) - numpy 2.3.3, scipy 1.16.2.
- "Convergence_and_correctness_tests.py" - imp, numpy 2.3.3, matplotlib 3.10.6, seaborn 0.13.2, pandas 1.3.0.
- "1S_and_2S_polarsiability.py" - imp, numpy 2.3.3, matplotlib 3.10.6.
- "Variations_with_depth.py" - imp, numpy 2.3.3, matplotlib 3.10.6.
- "Compare_scattering_rates.py" - imp, numpy 2.3.3, matplotlib 3.10.6, pandas 1.3.0, "Ionisation_data.xlsx", "3photon_ionisation_data.py".
- "Continuum_corrections.py" - imp, numpy 2.3.3, scipy 1.16.2, matplotlib 3.10.6.

Note that some functions were originally written using the scipy.misc.derivative() function which has since been removed. This has been replaced with appropriate use of scipy.optimize.approx_fprime(), the original lines have been commented out and can be re-enabled if desired.

## Authors
Joseph P. Scott

Department of Physics, Durham University, Durham, UK

joseph.p.scott@durham.ac.uk

## Date of last update
Wednesday 22nd October 2025
