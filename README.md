# The `OpenPMX` project

This is `OpenPMX`, a collection of numerical routines for estimation and simulation of mixed-effect pharmacokinetic and pharmacodynamic models.

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

# About this repository
This repository contains a copy of the latest stable version of `OpenPMX`.

# Installation
In the source directory run either:
- `./install` default build (using PThreads)
- `./install openmp` builds using OpenMP
- `./install pthreads` builds using PThreads
- `./install singlethread` use only a single threads

Once installed the library and a script `openpmx` will be built with the installed paths.

# Examples
Some examples are available. Open a shell in the particular example directory and run:

### Propofol
In directory `schnider` these are the Schnider PK data reported in Schnider TW, Minto CF, Gambus PL,
Andresen C, Goodale DB, Shafer SL, Youngs EJ. The influence of method of
administration and covariates on the pharmacokinetics of propofol in adult
volunteers. Anesthesiology 1998;88:1170-82. 
In this example a 3-compartment allometric model is estimated from the data.
- `../../openpmx schnider.gp`
- `../../openpmx schnider.c`
- `../../openpmx schnider.cpp`
- `../../openpmx schnider_diffeqn.gp`

### Theophylline
In directory `theo` these are the Theophylline data from the NONMEM distribution.
- `../../openpmx theo.gp`

### Warfarin
In directory `warfarin` - TODO: need to find and note cource
- `../../openpmx warfarin.gp`

# Philosophy
- Do one thing and do it well.
- Fewer dependencies are better.
- Reproducability is good.
- Fewer configuration parameters are better.
- Configuration parameters should as much as possible have a obvious value.
- Results should be as insensitive as possible to configuration parameters.

# TODO: Important
- should ext and phi files be removed if they arent produced? This can clash that they exist after a previous run 
- need consistent naming for options/config
- make a PKPDtools simulator and put it online!
	: TCI type advancer could be made, but where is its state being held?
	: probably the test DES solver can be reivsed to hold the state
	: but how to the covariates get sent to it?
- for symmetry we should have simulate_popmodel which we can comsider separating into a file
- remove yhat and yhatvar and replace with IPRED
	: but this will break NONMEM code
- rename POPDATA to IDATA or IPARAM
- OMEGA output should respect the blocks
- should do a gradient at posthoc and identify non-informative params. We should test the initial value again is OBJFN is the same then replace it!
	: this stops parameter accidental drift for non-informative values!
- warfarin example improves with allometric scaling! use this in example!
- make warfarin example following Holdord NONMEM 7 example in warfarin
	: it has RATE == -2 which is special NONMEM thing that I dont know how to program. I just flagged it as error in checkout
- first posthoc should run several times until individuals dont change so we have repeatable starts from finished models.
	: this is kind of what the posthoc runs were intended to do.
	: you can probably redo the refine step multiple times? Maybe onlly when all etas are zero.
	: you should probably re-run it at the end for posthoc
	: maybe need multiple runs at posthoc stage
- add screen feedback for nfunc so we can see the progress of very slow runs
- make comparison script which collects all performance metrics from a branch and compares to another branch and produces an output file which we can visually compare
- compare to straightforward cholesky encoding
- Work through Warfarin example from Holford
- fix compiler errors!!!
- it seems as ID in a dataset is case sensitive but all the rest are not?
- need repeatable seed for each run, so we can guarantee that each dataset is the same across runs
- add test as option so we can easily add it to the ananalysis and plots
- probably multiple runs should append to output files
	: maybe not, that the users responsibility, but what about multiple methods?
- test whether just minimizing cholesky makes eleveld example more stable with respect to restarts
- write C-code output to file so restarting is easy
- Make openpmx script work reading from a pipe instead of a control file
- I am not sure ALAG() is working robustly, make a regression test
- Probably the EVID and MDV does not match NONMEM style. Make these consistent to ease use of NONMEM compatible data files.
- Add diagnostics for estimation close to boundaries.
- Make test cases for OMEGA_SAME.
- From the unified regression tests, make a performance indicator so we can tweak default paramaters.
- Make sure examples pass tests with address sanitizer and with valgrind
- Re-add posthoc test for estimation
- add posthoc estimates of variances to resample output to show etas and errs are correct
- add regression tests to verify simulations
- remove C++ iterface by first release
- Unify fatal error and warning interface within openpmx.c
- warn about valid THETA variables after a THETA_INVALID. This catches if too few values used for THETA initilization.
- remove PHI file since its not correct, we have to save the eigenvalues and vectors
- Make full optimizing icov via cholesky decomposition as a comparator
- change Icov_resample naming to not conflict with resample in the sense of simulation
- Test NONMEM style of calculating Yhatvar at eta == 0 to see if it matches better
	: put this in the documentation
- We can do a single stage 1 with icov_resample and then calculate the population variability
  with weighted samples this could make a very fast EM algorithm, not unlike a unscented kalman filter
	: no, this could only work with mu-modelled models because otherwise there is no clear way to update
	theta values
 - Efficiency measures should count the function evaluations in the inner estimatio as well!
 - an encode object would contain extra information like fixed etc
	: have to maximum number of estimated params then, 100 maybe?
 - can you consider only printing non-default paramaters?
 - try piecewise linear encode transformation function
 - try NONMEM style transformation method
 - Control and output files and output are written in terms of $RECORDS
	: Appending output files allows repeat

# TODO: Maybe
- if during preprocess we would sort records to have all doses together and all observations together then advanace coud be simpler, no?
	: but this would be problems for tables I suppose because the records would be in wrong order
- consider using libprima for web-based solution, or maybe include in distribution? Performance is suggested to be better!
- when doing likelihood profiles, the besteta should not change from the estimated value, we dont want the start eta
	to differ depending on the process of runs along the ikelihood profile
- LOTS of error checking is necessary to be built in to the checkout function. The checkout function could check that derivatives are non-zero to identify paramaters that are not identifiable.
- Error limits via $COV step of likelihood profiles (which are probably most efficient).
- Make diagram of software components.
- Make prebuilt models for simulations Marsh, Schnider, Eleveld, Minto, etc.
- Estimation of probabilities directly as with LIKELIHOOD in NONMEM
- diagnostic in checkout for overlapping infusions
- Dont try to speed up estimation by reducing printing. It would only matter for very fast iterations and user can always use quiet or omit_* options.
- Maybe make kalman filter by just propogating the sigma points across iterations
- log should be shared across all runs?
- Separate openpmx settings like nthreads and filename from model results
- We could use a precondition matrix based on population covariance or the individual covariance matrix when
  doung the individual eta estimation. Then the initial steps make more sense (because they will be independant)
  and a size of 1 makes sense. Do this via cholesky decomposition? Can we make sure the covariance matrix wont collapse?
  This could probably be done more sensibly via eigenvalue decomposition just like the icov resample code
- If I make the evaluation in LUA then we can absolutely make a WEBassembly version!!!
- Should openpmxtran_data_preprocess address ID firled directly of should it calcculate it vis offsetof()?
- Try post-processing of estimate via likelihood profiles, just average the upper and lower +1 OBJFN and fix and estimate the rest
	this can be done for the THETA dnd OMEGA diagonals.
	the datasim analysis would be most interesting for this
 - I could recreate the LAPLACE as well using https://math.stackexchange.com/questions/2612795/4-point-like-central-finite-difference-for-second-partial-derivatives
 - https://github.com/saulwiggin/Numerical-Recipies-in-C/blob/master/README.md
 - https://www.astro.umd.edu/~ricotti/NEWWEB/teaching/ASTR415/InClassExamples/NR3/legacy/nr2/C_211/progs.htm
 - https://github.com/saulwiggin/Numerical-Recipies-in-C/tree/master
 
# TODO: Someday
- Consider piecewise linear paramater transformation functions.
- Consider SAEM and ITSB methods.
- Do exact root finding for individual covariance matrix because of the curse of high-dimensions importance sampling
- Consider removing libgsl dependency, probably only need matrix operations and Cholesky decomposition
- Make a web-version by compiling to Web-assembly with some simple fixed model types 1,2,3-compartment models. This would require for the C-code to call back into user javascript code?
- Make a R interface and R package?
- double check that no files are opened when we are omitting output tables
- Remove any log file, anything that goes to screen is log!
- Can I add the points tested for gradient calculation for the variance covariance matrix calculation! The OBJFN values will be quite low and the model is linear so probably it doesnt matter
- Use sampling around the final estimate to estimate a variance-covmatrix around the minimum. we can use the same encode-decode as for the popmodel covariance matrix, sampling at the
  sigma points, this can go through an optim, makbe even a nonlinear least squares method for speed
	: This allows model averaging for the final estimate and maybe this will correct for bias!!
- NONMEM RATE as -1 and -2 allows for estiumating duration etc, implement this?

# TODO: NOT
- Can we change the configuration between runs to make modifications and rerun? Or do we have to make a new PMXSTATE from this?
   -  We should skip this for the basis user. Expert users would use the full C interface.
- Output files should not contain any date or time information so that memoization and hashes can be used
- add include for hierarchial model development, this is a different problems and is better solved in another way
- should predict all be in ADVAN or table?
	: it cant be in table because we dont solve
	: it could be only in posthoc or predict run, no it cant because the user could rely on it for objective function
	: conclusion is that it must be in ADVAN



