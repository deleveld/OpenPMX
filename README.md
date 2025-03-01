# The OpenPMX project

Hi everybody out there doing pharmacometrics,

I doing a (free) model estimation software as a side project (wont be big a professional like NONMEM) for PCs. This has been brewing for a number of years and is starting to get ready
for broader use. I'd like feedback on what people like/dislike with the control, data, and table file structures. It somwhat resembles NONMEM due to practical reasons.

https://www.reddit.com/r/agedlikemilk/comments/1ikfx3s/linux_is_now_the_most_widely_used_os_in_the_world/#lightbox

linux_is_now_the_most_widely_used_os_in_the_world


------------

OpenPMX is a collection of numerical routines for estimation and simulation of mixed-effect pharmacokinetic and pharmacodynamic models.
This repository contains a copy of the latest stable version.

- First-order conditional estimation results similar to industry standard [NONMEM](https://www.iconplc.com/solutions/technologies/nonmem)
- Code and simulate/estimate models using text files and C allows scripting and integration with other tools
- Analytic models for common comparmental models, ODE solver for complex models
- Utilize multi-core CPUs via [OpenMP](https://www.openmp.org/) or [pthreads](https://man7.org/linux/man-pages/man7/pthreads.7.html)
- Limited project dependencies: [gcc](https://gcc.gnu.org/), [GSL](https://www.gnu.org/software/gsl/)

OpenPMX is free software, you can redistribute it and/or modify it under the terms of the GNU General Public License.
Unlike the licenses of proprietary pharmacometric software the license of OpenPMX does not restrict scientific cooperation. It allows you to share your programs freely with others.
The GNU General Public License does not permit this software to be redistributed in proprietary programs.

This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# Install
* Windows:
	1. Install [MSYS2](https://www.msys2.org/)
	2. Install [gcc](https://gcc.gnu.org), I used `pacman -S mingw-w64-x86_64-gcc`
	3. Install [GSL](https://www.gnu.org/software/gsl/), I used `pacman -S mingw-w64-x86_64-gsl`
* Ubuntu:
	1.  Install [GSL](https://www.gnu.org/software/gsl/) with `sudo apt install libgsl-dev`

Then in the source directory choose one of:
 * `./install` default build using [pthreads](https://man7.org/linux/man-pages/man7/pthreads.7.html)
 * `./install openmp` builds using [OpenMP](https://www.openmp.org/)
 * `./install singlethread` use only a single threads

After install there will be a shell script `openpmx` with the installed paths.

# Examples
Some examples are available.

### Propofol
In the directory `examples/schnider` are the Schnider PK data reported in Schnider TW, Minto CF, Gambus PL,
Andresen C, Goodale DB, Shafer SL, Youngs EJ. The influence of method of
administration and covariates on the pharmacokinetics of propofol in adult
volunteers. Anesthesiology 1998;88:1170-82. 

In this example a 3-compartment allometric model is estimated.
- `../../openpmx schnider.gp`

The same example, but using a differential equation solver.
- `../../openpmx schnider_diffeqn.gp`

### Theophylline
In directory `examples/theo` these are the Theophylline data from the NONMEM distribution.
- `../../openpmx theo.gp`

### Warfarin
In directory `examples/warfarin` contains data obtained from [here](http://clinpharmacol.fmhs.auckland.ac.nz/docs/warfarin.csv) and [here](https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://holford.fmhs.auckland.ac.nz/docs/pkpd-workshop-nonmem7.pdf).
The PK analysis example is:
- `../../openpmx warfarin.gp`

# Contributors
- Your name here (your@email.com) Massively important thing elegantly coded.
### Version 0.01
- Douglas Eleveld (deleveld@dds.nl) Initial release. 

# Wish-list / TODO
- During checkout, calculate paramater gradients to identify paramaters that may be numerically unidentifiable
- Make estimation of caterorical variables possible.
- Post-estimation evaluation of the objective function in the spirit of NONMEM and $COV. Possibly the gradient can be calculated at the final estimate and the first and second derivatives calculated using splines. It is also possible to transform this back into the scale of the user paramaters?
- Should there be Bessel correction on CWRES?
- Calculate NPDE.
- Add code to generate data for VPCs.
- Implement more analytic models from Abuhelwa AY, Foster DJ, Upton RN. ADVAN-style analytical solutions for common pharmacokinetic models. J Pharmacol Toxicol Methods 2015; 73: 42-8

# Philosophy
- Do one thing and do it well.
- Fewer dependencies are better.
- Reproducability is good.

