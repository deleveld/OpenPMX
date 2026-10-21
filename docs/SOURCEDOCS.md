# 


# include/openpmx.h

This defines the version number and the various maximum size of
theta, omega, and sigma.

Also the maximum number of fields in a record, the number of state
variables etc.

If these are too low they can be adjusted here and the library must
be recompiled.

# 


# src/advan/advan.c

The user creates an ADVANFUNCS object via a method named something 
like: 
`ADVANFUNCS* pmx_advan_FOO(const DATACONFIG* const dataconfig, const ADVANCONFIG* const advanconfig)`
This function allocates and constructs an object containing function 
pointers to construct, advance across records, and destruct an 
advancer object. This is used in ievaluate.c.

This file is the base class ADVAN which contains information that 
all of the advancer objects have in common, for example: time, 
state, lag, bioavailability, and the infusions running at that 
moment. Probably the most important function is `advan_advance()` 
which handles the infusion AMT, RATE, EVID and the various special 
flags like CMT. 

### Advancing through an individuals records

While advancing the init() function is called (which calls the users
`IMODEL()` code in openpmxtran) when:

+ First record of an individual
+ EVID is 3 (reset event) or EVID is 4 (reset-and-dose event)
+ For dose records EVID is 1 or EVID is 4 
+ At the time given by a call to `pmx_advan_inittime()` which is 
`INITTIME()` in openpmxtran

The moment a dose record is encountered it is saved in a buffer with
the time of the RECORD and the lag valid at the moment of 
processing. Subsequent changes in lag do not change the moment when 
the dose is actually applied.

Advancing is done in "steps" to the earliest of:

+ The next RECORD time is achieved
+ A previous infusion starts or stopped
+ A bolus dose is given
+ If `pmx_advan_inittime()` has been called, which is `INITTIME()` in openpmxtran
Calling `pmx_advan_amtlag()` (in openpmxtran this is `ALAG()`) 
delays the application of subsequent dosing. It does not change the 
moment of application of doses that are already in the buffer.


# src/advan/diffeqn_libgsl.c

This file implements the ODE solvers inplemented in libGSL.

The default setting are:
- ODE stepping function `gsl_odeiv2_step_rk8pd`

# src/advan/diffeqn_test.c

This file implements a very simple fixed-step size Runge-Kutta 
fourth order ODE solver. This is for testing purposes and shouldnt
be used for serious models.


# src/advan/eigen.c

This advancer uses eigendecomposition to obtain exact analytical
(no step size error) solutions for any linear compartmental model,
regardless of the number of compartments or topology.

Instead of hardcoding the analytical solution for a specific model
structure or using numerical ODE integration, this advancer:

  1. User must supply system matrix
  2. Eigendecomposes the matrix: A = V * diag(eigvals) * V^-1
  3. Advances state analytically in modal coordinates

This is equivalent to NONMEM's ADVAN5 (linear general compartmental
model). The technique is described in:

- Joachim J, "TCI Eigendecomposition" (2026)
- Minto C, Schnider T, "From Eigenvalues to Plasma Coefficients" (2026)

The modal state update for constant input over interval dt:

  x_modal(t+dt) = x_modal(t) * exp(eigval*dt)
                 + (V^-1 * J) / (-eigval) * (1 - exp(eigval*dt))

where eigval < 0 for stable systems, V is the eigenvector matrix,
and J is the input vector (infusion rates into each compartment).

Uses libGSL for eigendecomposition (gsl_eigen_nonsymmv) and matrix
inversion (LU decomposition).

For the eigensystem advancer the system matrix must be provided
in row-major order. The user can write this using
`pmx_advan_eigen_sysmat()` or via openpmxtran as `SYSMAT()`.

This file also implements a three compartment mammilary model
via the underlying eigensystem solver. It is available as
`pmx_advan_eigen_threecomp()` of via openpmxtran as
`$ADVAN(eigen_threecomp)`.

This file also implements a two compartment mammilary model via the
underlying eigensystem solver. It is available as
`pmx_advan_eigen_twocomp()` of via openpmxtran as
`$ADVAN(eigen_twocomp)`.

This file also implements a one compartment model with depot
via the underlying eigensystem solver. It is available as
`pmx_advan_eigen_onecomp_depot()` or via openpmxtran as
`$ADVAN(eigen_onecomp_depot)`.

This file also implements a two compartment model with depot
via the underlying eigensystem solver. It is available as
`pmx_advan_eigen_twocomp_depot()` or via openpmxtran as
`$ADVAN(eigen_twocomp_depot)`.


# src/advan/onecomp.c

This file implements an advancer for a one compartment model.


# src/advan/onecomp_depot.c

This file implements an advancer for a one compartment model with a 
depot compartment.

The equations were obtained from: Abuhelwa AY, Foster DJ, Upton RN.
ADVAN-style analytical solutions for common pharmacokinetic imodels.
Journal of pharmacological and toxicological methods. 2015 Jun 30;73:42-8.

This advancer does not yet support infusions. See the GitHub issue. 
The workaround is to code the model as a eigensystem or differential
equation.


# src/advan/pred.c

This file implements a simple predictor without state.


# src/advan/tci.c

This file implements at Target-Controlled-Infusion dose controller.

Jona Joachim (jona@joachim.cc) Contributed the code in `src/advan/ministan` which is 
a modernized version of [Stanpump](https://opentci.org/code/stanpump).


# src/advan/threecomp.c

This file implements an advancer for a three compartment mammilary 
model.

The equations were obtained from: Abuhelwa AY, Foster DJ, Upton RN.
ADVAN-style analytical solutions for common pharmacokinetic imodels.
Journal of pharmacological and toxicological methods. 2015 Jun 30;73:42-8.

Applied some AI guided optimizations for speed.


# src/advan/twocomp.c

This file implements an advancer for a two compartment mammilary 
model.

The equations were obtained from: Abuhelwa AY, Foster DJ, Upton RN.
ADVAN-style analytical solutions for common pharmacokinetic imodels.
Journal of pharmacological and toxicological methods. 2015 Jun 30;73:42-8.

Applied some AI guided optimizations for speed.


# src/bobyqa/bobyqa.c

Implementation of Mike Powell's BOBYQA algorithm for minimizing a function
of many variables.  The method is "derivatives free" (only the function
values are needed) and accounts for bound constraints on the variables.  The
algorithm is described in:

M.J.D. Powell, "The BOBYQA Algorithm for Bound Constrained Optimization
Without Derivatives."  Technical report, Department of Applied Mathematics
and Theoretical Physics, University of Cambridge (2009).

The present code is based on the original FORTRAN version written by Mike
Powell who kindly provides his code on demand (at mjdp@cam.ac.uk) and has
been converted to C by É. Thiébaut.
Copyright (c) 2009, Mike Powell (FORTRAN version).
Copyright (c) 2015, Éric Thiébaut (C version).
Read the accompanying `LICENSE` file for details.

Some small changes made to bobyqa.c by Douglas Eleveld 
deleveld@dds.nl to avoid compiler warnings about possibly unused 
unitialized values


# src/checkout.c

This file implements the high-level checkout functions. The actual
checkout function is individual_checkout() which is implemented in
ievaluate.c. It advances the individual records with ETA zet to 0
and has lots of checks to detect errors. This is important because
the other advancing of the individual records for objective function
calculation, prediction, etc do not have any check to maximize
speed.


# src/covariance.c

This file calculates the covariance matrix of the estimate to provide
information about upper and lower bounds of paramater estimates.

For off-diagonal values it is not really clear how to transform these from the
unbounded encoded space to the bounded paramater space, so no limits
are calculated for them.

All work is done in the encoded (transformed) parameter space using
central finite differences of the objective function (-2LL).
Two matrices are relevant:

- **R**: the Hessian of the total objfn (NONMEM's R matrix).
- **S**: the OPG/BHHH cross-product of per-individual gradients (S = J'J).

They are combined into the sandwich estimator `Cov = Rinv * S * Rinv`, or
used alone (`Sinv` or `Rinv`, scaled because the objfn is -2LL). 
Quadratic smoothing can then refine the result. The covariance type is
chosen with `COVARIANCE_S_MATRIX`, `COVARIANCE_R_MATRIX`,
`COVARIANCE_SANDWICH` which is the default.

## Reporting

- Correlation pairs sorted by absolute value are reported.
- Eigen-analysis of S and/or R finds all parameters confined to the
  near-null subspace (eigenvalues below `reltol * lambda_max`), also
  reporting negative eigenvalues.
- When diagnostic samples are requested the predicted vs observed
  dobjfn is written to the diagnostics file and summarised
  (log-ratio mean and sd).

# src/dataconfig/dataconfig.c

This file initializes a RECORDINFO object from a DATACONFIG object.
In this way the offsets of some important RECORD fields are cached
for rapid access, 


# src/defines.h

This file define the file extenstions and numeric formats used.


# src/encode.c

This file performs the encoding and decoding of a POPMODEL to and 
from an unbounded vector. This is a part of the Stage 2 of 
optimization.

Only non-fixed parameters are encoded, in the order THETA, SIGMA,
then the lower triangle of the non-fixed OMEGA. Encoding is relative
to an offset vector, so a zero vector reproduces the POPMODEL given
to encode_transform(). THETA is mapped from its (lower, upper) bounds
(compile-time selectable: ENCODE_SPECIAL, ATANH, LINEAR, NONMEM).
SIGMA and the OMEGA diagonal are log-transformed (as standard
deviations if ESTIMATE_SD_NOT_VAR is defined). OMEGA off-diagonals
are encoded via the Cholesky factor of the correlation matrix and
canonical partial correlations (atanh transformed). The file also
handles OMEGA_SAME block filling, and the creation of parameter
labels such as THETA(1) and OMEGA(2,1).

Calling encode_transform() makes such that a zero-initialized vector 
will reproduce the given POPMODEL. 

Different transformations are possible to address theta bounds. 

Sigma is log-transformed. 

The diagonal of omega is log-transformed.

The off-diagonals of omega are encoded as a specialized Cholesky decomposition of a 
correlation matrix. Only the lower triangular is encoded. The 
approach is described in: <https://mc-stan.org/docs/reference-manual/transforms.html#cholesky-factors-of-correlation-matrices>
and also: <https://mc-stan.org/docs/reference-manual/transforms.html#correlation-matrix-transform.section>

If an OMEGA_SAME block is used then (OMEGASAME() in openpmxtran)
then omega is filled by block transposed by its size. Basically it
reproduces a block of omega by a lower block. This allows two eta
values with the same variances to be estimated.


# src/estimate.c

This file does the outer (stage 2) estimation.

OpenPMX estimates population models in two nested stages. Stage 1
(see stage1.c) finds, for fixed population parameters, the individual
random effects (etas) that minimize each individual's objective
function contribution, and evaluates the first-order conditional
estimation with interaction (FOCEI) approximation of the likelihood.
Stage 2, implemented here, searches over the population parameters
(theta, omega and sigma) to minimize the total objective function
summed over all individuals. Every objective function evaluation of
stage 2 therefore runs stage 1 for all individuals, in parallel over
the available threads.

The public entry points are pmx_estimate() and pmx_evaluate().
Evaluation is estimation with a single objective function evaluation:
the model is not changed, but the individual etas, the objective
function, the predictions and the output files are produced exactly
as they would be for an estimation.

Parameter encoding

The population parameters are not optimized directly. They are encoded
(see encode.h) as an unconstrained vector of offsets relative to the
current model, so an all-zero vector reproduces the current model.
This puts parameters of very different scales on a common footing.
The step sizes used by the optimizer are therefore in encoded units,
not in the units of the parameters. Only the parameters that are
actually estimated are part of the vector. If there are none, no
optimization is done after the initial evaluation.

Optimization

The objective function is minimized with BOBYQA, a derivative-free
bounded quadratic-model trust region optimizer. Estimation runs in
phases, controlled by the estimate options (maxeval, step_initial,
step_refine, nsig, dobjfn):

 1. Initial evaluation. The model as given is evaluated once, which
    sets the starting objective function and the starting etas. A
    warning is given if it is not finite.
 2. Initial optimization. BOBYQA runs with a trust region radius
    from step_initial down to step_refine. This makes large changes
    to the parameters and does the bulk of the work.
 3. Final step size. The smallest trust region radius that can still
    resolve nsig significant digits in the parameters is found by
    halving a radius, starting at half of step_refine, until
    perturbing each encoded parameter by plus or minus that radius
    leaves the decoded parameters agreeing to more than nsig digits.
    The search stops at 1e-8.
 4. Refinement. BOBYQA is restarted from the best model so far with a
    radius from step_refine down to the final step size. After each
    restart the number of significant digits (nsig) shared by the
    parameters before and after the restart is reported, as the
    minimum over all theta, sigma and the lower triangle of omega.
    Restarts continue until nsig exceeds the requested nsig
    (converged), until the change in the objective function between
    restarts is smaller than dobjfn (stalled), until BOBYQA reports
    an error, or until the evaluations allowed by maxeval are used up.

The result is marked as final only if the refinement converged,
otherwise it is the best model found so far. BOBYQA errors are
reported as warnings and do not discard the progress made.

Best model so far

The best model is tracked separately from the model under test and is
returned to the caller, so the result is always the best model seen
and not the last model evaluated. Only a finite objective function
is ever accepted as best: NaN and infinite values are rejected, and
a non-finite starting objective function is replaced by the first
finite one. When the best model improves, the etas of every individual
are saved and used as starting values for the next evaluations.

Output

Progress is written to the screen and, if the model has a filename, to
the out file. By default a single line is updated in place on a
terminal, the details and verbose options print every improvement or
every evaluation instead. At the start a header with the version,
the parallel configuration and the data is written, followed by a
data checkout (see checkout.c and ievaluate.c) that detects errors
before any time is spent on estimation. During estimation each
improvement is appended to the ext file. At the end the phi file
(individual parameters and etas) and the yhat file (pred, yhat,
yhatvar and the state) are written and the ext file gets its trailer.

Each time the objective function improves during estimation the eta
values are saved. They will be used as initial values for subsequent
optimization.

After initial optimization, a smaller search space is used from 
rho_refine to rho_final. This is repeated until the change in 
parameters between restarts is greater than nsig or objective 
function between restarts is less than dobjfn.

Model estimation begins with an intial optimization with large 
changes to paramater values with BOBYQA rho values from 
step_initial, stopping at step_refine. 

At start of estimation the header of the ext file is written.

Before estimation a data checkout is done (see checkout.c and 
ievaluate.c) to detect various errors. 

At the end of estimation the phi file is written and a trailer is 
put onto the ext file. Also the yhat file is written with prediction
varaiables including pred, yhat, yhatvar, and the state.

Evaluation is the same as estimation but with maxeval=1. Note that 
maxeval=0 would trigger the default number of evaluations.


# src/idata.c

This file handles the individualized data of a population. This 
include the pointers to the individual data records, to the 
initialized models (IMDOEL), state, etas, predictions (PREDICTVARS),
It also keeps track of 4 of the 5 terms in the objective function.

Here the objective function for each individual is summed from its
components. The terms match those from the individual objective
function equation in the manuscript. 
The omega_nonzero_lndet term is the log(det(omega)) term and is 
calculated in omegainfo.c.


# src/ievaluate.c

This file uses the advan to evaluate each individial, advancing the
model over the data records. This is done for various reasons:
during optimization to calculate the objective function, for 
prediction of PREDICTVARS after Stage 1, to checkout the data before
estimation to detect errors, and to perform simulations.

individual_fasteval() is the optimized core used by the objective
function: it advances the ADVAN over the records and accumulates
separate sums for continuous observations (squared error divided by
the prediction variance, and its log) and for discrete observations
given as a log-likelihood. individual_evaluate() does the same but
also stores per-record predictions, variances and states.
The prediction variance (yhatvar) is obtained by propagating the
SIGMA error through the model with central differences around zero
error. loglik_llq() gives the log-likelihood for observations below
the limit of quantification. individual_simulate() advances the
model while applying supplied residual errors. All functions work
on a single individual and must be thread safe.

The variance of y prediction (yhatvar) is estimated by central 
differences around err=0.

### Checkout

Before estimation a data checkout is done to detect various errors.

+ If RATE exists but AMT does not it is an error 
+ Non-integer ID is probably an error.
+ Non-integer CMT is an error.
+ Time must increase except for reset events.
+ Some checking is done that state is not accessed outside of its limits.
+ Observations are given when EVID is 0.
+ A warning is given if DV is 0 for an observation.
+ If any state values become non-finite is an error.

# src/linalg.c

This file contains linear algebra functions for Cholesky 
decomposition. This is necessary for calculating log(det()) of a 
covariance matrix, and calculating the sample likelihood.


# src/models/propofol_eleveld.c

This file implements the Eleveld propofol model.

Eleveld DJ, Colin P, Absalom AR, Struys MM. Pharmacokinetic–pharmacodynamic
model for propofol for broad application in anaesthesia and sedation.
British journal of anaesthesia. 2018 May 1;120(5):942-59.



# src/models/propofol_schnider.c

This file implements the Schnider propofol model.

Schnider T, Minto C, Gambus P, Andresen C, Goodale D, Shafer S, Youngs E:
The influence of method of administration and covariates on the pharmacokinetics
of propofol in adult volunteers. Anesthesiology 1998; 88:1170–82 PMID: 9605675



# src/models/remifentanil_eleveld.c

This file implements the Eleveld remifentanil model.

Eleveld DJ, Proost JH, Vereecke H, Absalom AR, Olofsen E, Vuyk J,
Struys MM. An allometric model of remifentanil pharmacokinetics and
pharmacodynamics. Anesthesiology. 2017 Jun 1;126(6):1005-18.



# src/models/remifentanil_minto.c

This file implements the Minto remifentanil model.

Minto CF, Schnider TW, Egan TD, Youngs E, Lemmens HJ, Gambus PL,
Billard V, Hoke JF, Moore KH, Hermann DJ, Muir KT. Influence of age
and gender on the pharmacokinetics and pharmacodynamics of
remifentanil: I. Model development. Anesthesiology. 1997 Jan 1;86(1):10-23.

Minto CF, Schnider TW, Shafer SL. Pharmacokinetics and pharmacodynamics
of remifentanil: II. Model application. Anesthesiology. 1997 Jan
1;86(1):24-33.



# src/models/remimazolam_eleveld.c

This file implements the Eleveld remimazolam model.

Eleveld DJ, Colin PJ, Van den Berg JP, Koomen JV, Stoehr T, Struys MM. 
Development and analysis of a remimazolam pharmacokinetics and 
pharmacodynamics model with proposed dosing and concentrations for 
anaesthesia and sedation. British Journal of Anaesthesia. 
2025 Jul 1;135(1):206-17.



# src/models/sufentanil_gepts.c

This file implements the Gepts sufentanil model

Gepts E, Shafer SL, Camu F, Stanski DR, Woestenborghs R, Van Peer A,
Heykants JJ. Linearity of pharmacokinetics and model estimation of
sufentanil. Anesthesiology. 1995 Dec 1;83(6):1194-204. 



# src/openpmxtran.c

openpmxtran: translator from OpenPMX control files to C source.

Reads a NONMEM-like control file made of `$BLOCK` sections ($DATA,
$ADVAN, $IMODEL, $PREDICT, $DIFFEQN, $THETA, $OMEGA, $OMEGABLOCK,
$OMEGASAME, $SIGMA, $FUNCTIONS, $MAIN, $FILE), strips C-style
comments, splits the text into blocks and dispatches each to a
handler that fills a RESULT structure. The header of the data file
is read to determine the RECORD fields, names are checked against
the model, and everything is expanded into the embedded C template
(`openpmxtran_template`). The generated program defines the model
callbacks (imodel_init, imodel_diffeqn, imodel_predict), the data
preprocessing, the initial THETA/OMEGA/SIGMA values and a main()
that runs the user's $MAIN code (estimate, evaluate, simulate,
covariance, table, ...). Any $FILE sections are written out as
auxiliary files alongside the generated source.

Any text preceeding the blocks defined by $... is ignored and can
be used for a description of the control file.

The `$DATA(...)`	 block defines the data filename and any code for
manipulation of the data file before analysis.

The `$ADVAN(...)` block defines the advancer type and any options
for the advancer.

The `$IMODEL(...)` block defines the model paramaters and the code
to initialize them.

The `$PREDICT(...)` block defines the prediction paramaters and the
code to initialize them. The predition for model fitting is returned
as `Y`.

The `$DIFFEQN` block defines the differential equation for the
advan solvers that require this.

The `$THETA` block defines the fixed effects. These can be either
given as `{ lower, initial, upper, ESTIMATE },` or in NONMEM style
as `(lower, initial, upper)` or as `(value FIXED)`.

The `$OMEGA()` block defines a diagonal omega matrix block of random
effects. The values are separated by `,` or whitespace. Negative
variances are treated as fixed. 

The `$OMEGABLOCK()` block defines a block omega matrix block of
random effects. The values are separated by `,` or whitespace.
Negative variances along the diagonal are treated as fixed. 

The `$OMEGASAME(...)` block defines a block which copies the values
from lower indicies in the omega matrix. This is used for
inter-occasion variability, i.e. multiple etas with variances and
covariances matching those of lower indicies.

The `$SIGMA()` block defines a diagonal omega matrix block of
residual variability. The values are separated by `,` or whitespace.
Negative variances are treated as fixed. 

The `$FUNCTIONS` block defines C code which will run after problem
initialization. 

The `$MAIN` block defines C code which will run after problem
initialization. 

Sections defined by `$FILE(...)` are written to the control filename
appended by '.' and then the name given. For the `graphics.R` script
this can be used to provide a configuration file for modifying the
graphics produced.

Within the datafile a `.` as a data entry is replaced by NAN.

Within the `$ADVAN(...)` block an advancer must be indicated.


# src/profile.c

This file implements a function that allows likelihood profiling.
It is available in openpmxtran as `profile()`.

If .dobjfn is 0. then it will be set to the chi-squared distribution 
for alpha 0.01 and 1 degree of freedom. This is approximately 6.63.

The default tolerance is 1 objfn units. This is usually enough that a
smoothed line through the profile points is close to the confidence
limits.

The default number of root finding evaluations is 10.

If .value is 0 and .lower < .upper then two profiles are run, one to
the .lower bound and one to the .upper bound. On return .lower is the
crossing point from the lower profile and .upper is the crossing point
from the upper profile. The .value is set to the source estimate. 
Otherwise a single profile is run to .value. If a bound does not 
cross the .dobjfn threshold then the limit on that side is open and 
is returned as -INFINITY (lower) or +INFINITY (upper).

pmx_profile_add() evaluates the profile at the single value in args->param.value and
append the point to the profile file. The profiled parameter is fixed
and the others are re-estimated starting from the source fit. Only
args->param and args->estimate are used.


# src/reload.c

This file implements a function that modifies an OPENPMX object
from information from a .ext file. It is available in openpmxtran
as `reload()`. This loads the population parameters, setting their
structure to that in the file.

Each line of the .ext file is read and saved. This allows
continuation of terminated estimation runs from the population
parameters from the last full iteration.

The default is that reloading of the population parameters fails if: 

+ Any theta bounds dont match
+ Any theta FIXED/ESTIMATE dont match
+ The structure of any the omega blocks differ
+ Any omega values fixed or estimated dont match
+ The number of sigma values differ
+ Any sigma values fixed or estimated dont match

For the configuration object there are various settings.

- `.filename="...",` The filename is used as an .ext file to
load the population parameters (theta, omega, and sigma values)
and save them in the OPENPMX object. The default behavior is to load
the population parameters from the filename of the destination
OPENPMX object with a .ext extension added.
- `.force=true,` The loaded population parameters are copied over
those of the destinantion ignoring any mismatch in structure. The
default is to exit the program if the structures differ.
- `.optional=true,` If reloading fails then the destination OPENPMX
object is not modified and the program continues. The default
behavior is to consider this a fatal error and the program exits.
- `.silent=true,` Setting this supresses a message showing the
structure of the newly loaded population parameters.

A successful reload invalidates the objective function result.

A successful reload invalidates any existing individual level data
in the destination OPENPMX object because the number of etas could
have been changed. 


# src/scatter.c

This file implements a threadpool of waiting processes to process
an individual. The task must be of a THREADTASK type and should only
touch the data of the individual passed. 

The threadpool is implemented by pthreads or OpenMP. A
single-threaded scatter is also possible.

The number of threads is limited to the number of individuals.
Otherwise threads will always be waiting and do nothing.
The main thread does work as well so we make nthread-1 worker
threads.

When scattering tasks, we sort individuals based on thier expected
runtime. We do the slowest first so we are the most efficient
because of lower risk of one long process slowing down finishing the
last task. This avoids the "long-tail" inefficiency.


# src/server_queue.c

This file implements the server queue. 

It compiles to nothing when not compiling with the install server
option.


# src/simulate.c

This file implements the top-level function for simulation.

For simulation random observation residuals are chosen according to
sigma. 

For simulation random ETA values are chosen according to the omega
matrix.

Upon simulation the individualized values the objective function components are set to zero.

Upon simulation the individualized eta values are randomly selected from 
from the distribution of the omega matrix.

After simulation a prediction occurs and the DV in the dataset is
replaced by observation with noise. PRED is set to zero and yhat is
the observation without noise. 


# src/stage1.c

This file does the inner (stage 1, conditional, posthoc) estimation.

All of the functions in this file only touch individual data so they
can be called in parallel, i.e. individual can be processed in seperate
threads.

If all of the eta values are zero then an initial estimation run is done
with step sizes from step_initial to step_refine. 

After this a refinement estimation run is done with step sizes from
step_refine to step_final. 

After minimization the individual covariance matrix is calculated
using the derivative of predictions with respect to eta i.e.
the jacobian. This is the third term in the objective function.
The approach does take into account the changes in V for different
eta values, so it properly handles eta-sigma interaction.

The inner (Stage 1) optimization only optimizes the first, second, 
and third terms in the objective function. The fourth term is not
dependant on the individual and the fifth term is calculated at the
minimum of the first three terms.


# src/table.c

This file implements writing tables. 

If a .filename is defined in the TABLECONFIG then the output is
written to that filename.  

If a .name is defined in the TABLECONFIG then the output is
written to a filename which is the name in the PMX object extended
with the given name.   

if a .stream (FILE*) is provided then table will be output to it.
For example it could be stdout or stderr. If given in this way then
fclose() is not called on the stream when the table is closed.

If no TABLECONFIG, then see if we have a filename in PMX then make
extension with numbered tables.

If no TABLECONFIG and no PMX filename, then use stdout

The table fields are indicated by a string which will be tokenized
with whitespace or comma.

Table generation is single threaded because no advancing takes 
place, only the saved state is used at each step. Thus it is I/O
limited and not cpu limited.

If firstonly flag in the TABLECONFIG is set then only the
first record of each individual in output to the table.

All of the fields of the data record are accessible.

All of the fields of $IMODEL(...) are accessible.

All of the fields of $PREDICT(...) are accessible.

The follwing fields are accesible: YHAT, YHATVAR, PRED, OBJ, INEVAL,
EVID, MDV, EVID, AMT, RATE, and CMT.

