# The `OpenPMX` project

This repository contains a copy of the latest stable version of `OpenPMX`, collection of numerical routines for estimation and simulation of mixed-effect pharmacokinetic and pharmacodynamic models.

- First-order conditional estimation results similar to industry standard [NONMEM](https://www.iconplc.com/solutions/technologies/nonmem)
- Code and simulate/estimate models using text files allows scripting and integration with other tools
- Analytic models for common comparmental models, ODE solver for complex models
- Utilize multi-core CPUs via [OpenMP](https://www.openmp.org/) or [pthreads](https://man7.org/linux/man-pages/man7/pthreads.7.html)
- Few dependencies [gcc](https://gcc.gnu.org/), [GSL](https://www.gnu.org/software/gsl/).
- Written in C, usable from C and C++

`OpenPMX` is free software, you can redistribute it and/or modify it under the terms of the GNU General Public License.
Unlike the licenses of proprietary pharmacometric software the license of `OpenPMX` does not restrict scientific cooperation. It allows you to share your programs freely with others.
The GNU General Public License does not permit this software to be redistributed in proprietary programs.

This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# Philosophy
- Do one thing and do it well.
- Fewer dependencies are better.
- Reproducability is good.

# Install
In the source directory run either:
- `./install` default build (using PThreads)
- `./install openmp` builds using OpenMP
- `./install singlethread` use only a single threads
This will build the library and a script `openpmx` with the installed paths.

# Examples
Some examples are available. Open a shell in the particular example directory and run:

### Propofol
In directory `schnider` the Schnider PK data reported in Schnider TW, Minto CF, Gambus PL,
Andresen C, Goodale DB, Shafer SL, Youngs EJ. The influence of method of
administration and covariates on the pharmacokinetics of propofol in adult
volunteers. Anesthesiology 1998;88:1170-82. 

In this example a 3-compartment allometric model is estimated.
- `../../openpmx schnider.gp`

The same example, but using a differential equation solver.
- `../../openpmx schnider_diffeqn.gp`

### Theophylline
In directory `theo` these are the Theophylline data from the NONMEM distribution.
- `../../openpmx theo.gp`

### Warfarin
In directory `warfarin` data obtained from [here](http://clinpharmacol.fmhs.auckland.ac.nz/docs/warfarin.csv) and [here](https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://holford.fmhs.auckland.ac.nz/docs/pkpd-workshop-nonmem7.pdf).
The PK analysis example is:
- `../../openpmx warfarin.gp`

