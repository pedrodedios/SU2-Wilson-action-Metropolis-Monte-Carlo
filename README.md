# SU-2-Wilson-action-with-Monte-Carlo
Monte Carlo simulation of the SU(2) Wilson action in a 4D lattice. We compute plaquette and specific heat values

## 🌟 Highlights

This repository contains a Fortran implementation of a Metropolis Monte Carlo algorithm for SU(2) lattice gauge theory using the Wilson action.
The code is intended for numerical studies of pure gauge theories on a four-dimensional hypercubic lattice.

Main features:

-**Metropolis Monte Carlo updates of SU(2) link variables**

-**Periodic boundary conditions on a 4D lattice**

-**Hot starts using Haar-distributed SU(2) matrices**

-**Computation of fundamental thermodynamic observables: Average energy ⟨E⟩ and specific heat $\chi_{\beta}$** 

	
For a detailed theoretical overview of SU(2) lattice gauge theory with the Wilson action, see the accompanying documentation: **[Lattice_QCD.pdf](https://github.com/user-attachments/files/24371960/Lattice_QCD.pdf)**



## 🚀 Usage

1️⃣ Compile the utility module

This module contains lattice indexing routines and SU(2) matrix utilities.

gfortran -O3 -fopenmp -c gauge_utils.f90

2️⃣ Compile the Metropolis module

This module implements the Metropolis update algorithm and measurements.

gfortran -O3 -fopenmp -c metropolis_module.f90

3️⃣ Compile the main program

The main program sets simulation parameters and runs the Monte Carlo loop.

gfortran -O3 -fopenmp -c wilson_gauge.f90

4️⃣ Link all object files

gfortran -O3 -fopenmp gauge_utils.o metropolis_module.o wilson_gauge.o -o wilson_gauge

5️⃣ Run simulations

Simulation parameters (lattice size, number of sweeps, β values, etc.) are defined directly in the main source file.

./wilson_gauge
