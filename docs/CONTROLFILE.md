# OpenPMX Control File Reference

OpenPMX control files describe a pharmacokinetic/pharmacodynamic model and the analysis to perform. 

Control files can have any extention. The examples use `.gr`. 

The `openpmx` script calls `openpmxtran` which reads the control file and generates C-source code to perform the
analysis. The C-source code is then compiled, linked with `libopenpmx.a`, and the program is executed
to perform the analysis and output the results.

## General info

* The control file is not fully parsed, but portions are passed 
through to the C compiler. Thus many usage constructs come directly from the C language. 
OpenPMX does not check for all illegal syntax and relies on 
the C compiler to provide error and warning messages.
* Arrays such as `THETA()`, `ETA()` and `CMT` are 1-based when accessed within code generate
by `openpmxtran`. This eases compatability with NONMEM code. Internally, OpenPMX
uses 0-based arrays native to the C language.
* If variances in `$OMEGA*` and `$SIGMA` along a diagonal are given as 0 or negative then these are treated as fixed
and not optimized.

---

## Block overview

Control files consist of blocks which are defined by a `$` after a newline and followed by a block name.

Any text before the first `$` marker is treated as a preamble and ignored. This is a good
place to document the control file. 

Comments use C syntax: `/* ... */` and `//`. These are stripped from the control file 
when read in by `openpmxtran`.

| Block name | Arguments | Required |
|---|---|---|
| [`$DATA`](#data) | `$DATA("file.csv")` | No |
| [`$ADVAN`](#advan) | `$ADVAN(method)` | Yes |
| [`$IMODEL`](#imodel) | `$IMODEL(param, ...)` | Yes |
| [`$DIFFEQN`](#diffeqn) | `$DIFFEQN` | Only with ODE advans |
| [`$PREDICT`](#predict) | `$PREDICT(var, ...)` | No (warning if absent) |
| [`$THETA`](#theta) | `$THETA` | No |
| [`$OMEGA`](#omega) | `$OMEGA(val, ...)` | No |
| [`$OMEGABLOCK`](#omegablock) | `$OMEGABLOCK(val, ...)` | No |
| [`$OMEGASAME`](#omegasame) | `$OMEGASAME(n)` | No |
| [`$SIGMA`](#sigma) | `$SIGMA(val, ...)` | No |
| [`$FUNCTIONS`](#functions) | `$FUNCTIONS` | No |
| [`$MAIN`](#main) | `$MAIN` | Yes |
| [`$FILE`](#file) | `$FILE("name")` | No |

Multiple `$OMEGA*` blocks are allowed and define successive blocks of the omega matrix. All other blocks may appear at most once.

---

## $DATA

Specifies the dataset and optional per-record preprocessing.

The data file is read by `openpmxtran` and translated to a C language array.

The data file can be whitespace or comma separated and must start with a header.

The header names defines the variable names so they must be valid C language identifers.

All data values are treated as `double` variables. Values given as `.` are translated to `NAN`.

### Arguments

```
$DATA("data.csv")
```

Alternatively the `$DATA` block can be omitted and the data file passed as the second
argument to `openpmx` script.

### Content

Optional C code executed once per data row, before any analysis. 

All CSV column names are available as writable `double` variables.

Any modifications to variables are written back to the data array in-place at startup. 

```c
$DATA("warfarin.csv")
    RATE = 0.;
    if (isnan(AMT) && MDV == 1)
        AMT = 0.;
```

### Recognized CSV columns

Columns are discovered automatically from the CSV header row. Standard columns:

| Column | Description | Required |
|---|---|---|
| `ID` | Subject identifier | Yes
| `TIME` | Time point | Yes
| `DV` | Dependent variable / observation | Yes
| `EVID` | Event type: `0` = observation, `1` = dose, `2` = other | No
| `MDV` | Missing DV flag: `0` = present, `1` = missing | No
| `AMT` | Dose amount | No
| `RATE` | Infusion rate | No
| `CMT` | Compartment number | No
| Any other | Covariates (e.g. `WGT`, `AGE`)  | No

All data record fields are accessible in `$IMODEL` and `$PREDICT` and `table(...);` within `$MAIN`.

---

## $ADVAN

Selects the method to advance the compartmental model over time and configures its options.

### Arguments

```
$ADVAN(method)
```

### Supported methods

| Method | Description |
|---|---|
| `pred` | Prediction only; no compartmental state |
| `onecomp` | One-compartment analytical solution |
| `onecomp_depot` | One-compartment with depot (note: infusions not supported) |
| `twocomp` | Two-compartment mammillary analytical solution |
| `threecomp` | Three-compartment mammillary analytical solution |
| `eigen` | Linear eigensystem (any topology); requires `SYSMAT()` in `$IMODEL` |
| `eigen_twocomp` | Specialized two-compartment eigensystem |
| `eigen_threecomp` | Specialized three-compartment eigensystem |
| `eigen_onecomp_depot` | Specialized one-compartment with absorption eigensystem |
| `eigen_twocomp_depot` | Specialized two-compartment with absorption eigensystem |
| `diffeqn_libgsl` | Numerical ODE via GSL RK45; requires `$DIFFEQN` |


### Content 

The block content is a C language designated initializer, thus comma separated  
key-value pairs. These define the options of the advancer:

| Option | Type | Default | Description |
|---|---|---|---|
| `.firstonly` | `bool` | `false` | Call `$IMODEL` code only for the first ecord per individual |
| `.predictall` | `bool` | `false` | Compute predictions for all records, not just observations |
| `.nstate` | `int` | — | Number of state variables (required for `diffeqn_*` or `eigen_*` advancers) |
| `.args.diffeqn.abstol` | `double` | 1e-9 | Absolute tolerance for ODE solver |
| `.args.diffeqn.reltol` | `double` | 1e-9 | Relative tolerance for ODE solver |
| `.args.diffeqn.hstart` | `double` | 0.1 | Initial step size for ODE solver |
| `.args.diffeqn.steptype` | `string` |  | ODE stepping method (e.g. `"rkf45"`) |

```c
$ADVAN(threecomp)
    .firstonly = true,
```

```c
$ADVAN(diffeqn_libgsl)
    .nstate = 3,
    .args.diffeqn.abstol = 1e-9,
    .args.diffeqn.reltol = 1e-9,
```

### Eigensystem solver

This `eigen_*` advancers use eigendecomposition to obtain exact analytical
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

#### Eigensystem advancer for a three-compartment mammilary model

```c
double k10 = CL / V1;
double k12 = Q2 / V1;
double k21 = Q2 / V2;
double k13 = Q3 / V1;
double k31 = Q3 / V3;
SYSMAT(-k10-k12-k13,	 k21,	 k31,
	    k12,			-k21,	 0.,
	    k13, 			 0., 	-k31); 
```


---

## $IMODEL(param, ...)

Declares the individual model parameters and the C code to initialize them per individual.

The model paramaters are used as C language struct members 
so they must be valid C language identifers.

### Arguments

Parameter names become fields of the `IMODEL` struct. Names may be comma- or space-separated. Reserved names (e.g. `THETA`, `ETA`, `A`) are rejected.

### Content

C code executed once for every record (or optionally on for the first record) 
while advancing over the records of an individual. All of the model parameters should be initialized.

Withing the function the following are available read-only:

| Symbol | Description |
|---|---|
| `THETA(i)` | Population parameter i (1-indexed) |
| `ETA(i)` | Individual random effect i (1-indexed) |
| `A(i)` | Compartment amount i at current time (1-indexed) |
| CSV column names | All data columns (e.g. `WT`, `AGE`) |
| Declared `$IMODEL` params | Parameters defined on earlier lines |

All `$IMODEL(...)` paramaters are accessible in `$PREDICT` and `table(...);` within `$MAIN`.

Special macros available in the body:

| Macro | Description |
|---|---|
| `ALAG(cmt, t)` | Set absorption lag time for doses to compartment `cmt` |
| `BIOAVAIL(cmt, f)` | Set dose bioavailability fraction for compartment `cmt` |
| `A_0(cmt, v)` | Initialize compartment `cmt` to value `v` at time zero |
| `SYSMAT(...)` | Define the system matrix for `eigen` advans (row-major) |
| `INITTIME(t)` | Set model start time |
| `INITCOUNT` | Number of times initialization has been called for this individual |

#### $IMODEL definition for a three-compartment mammilary model
```c
$IMODEL(V1, V2, V3, CL, Q2, Q3)
	const double SIZE = WT/70.;
	V1 = THETA(1) * SIZE * exp(ETA(1));
	V2 = THETA(2) * SIZE * exp(ETA(2));
	V3 = THETA(3) * SIZE * exp(ETA(3));
	CL = THETA(4) * pow(SIZE, 0.75) * exp(ETA(4);
	Q2 = THETA(5) * pow(SIZE, 0.75) * exp(ETA(5));
	Q3 = THETA(6) * pow(SIZE, 0.75) * exp(ETA(6));
```

### Target-controlled-infusion (TCI) dose controller

Dosing can be driven by a TCI controller. Within the `$IMODEL(...)` call `TCINIT(...)`
with a `TCICONFIG` object containing the controller options. The first call does setup
for the controller. Subsequent calls return immediately. The initial target is
0 so no dosing is given.

If it is computationally expensive to initialze the controller then `TCISTARTED()`
is available which returns if the controller has been started or not. This allows
the call to `TCIINIT()`, and the constrcution of the `TCICONFIG` object to be skipped.

Calling `TCITARGET()` sets the target for the TCI controller. This function passes
doses to the advancer. Setting a negative target turns off the TCI controller. The
function return the cumulative dose given so far since initialisation of the controller.

Doses given by the TCI controller are independent of any in the data records
`AMT` and `RATE` coloumns.

Some models for the TCI controller are available in `openpmx_model.h`.

```c
PROPOFOL_ELEVELD model = pmx_model_propofol_eleveld(&(PROPOFOL_ELEVELD_COVARIATES) {
	.age = age,
	.weight = wgt,
	.height = hgt,
	.female = (m1f2 == 2) ? true : false,
	.opiates = true,
});
TCIINIT(.k10 = model.k10,
		.k12 = model.k12,
		.k21 = model.k21,
		.k13 = model.k13,
		.k31 = model.k31,
		.ke0 = model.KE0,
		.vc = V1,
		.target_effect = true);
TOTAMT = TCITARGET(modelE50);
```
### TCI models available

Some models for the TCI controller are available in `openpmx_model.h`.

#### `pmx_model_propofol_schnider();`
	
Schnider T, Minto C, Gambus P, Andresen C, Goodale D, Shafer S, Youngs E:
The influence of method of administration and covariates on the pharmacokinetics
of propofol in adult volunteers. Anesthesiology 1998; 88:1170–82 PMID: 9605675

#### `pmx_model_propofol_eleveld();`

Eleveld DJ, Colin P, Absalom AR, Struys MM. Pharmacokinetic–pharmacodynamic
model for propofol for broad application in anaesthesia and sedation.
British journal of anaesthesia. 2018 May 1;120(5):942-59.

#### `pmx_model_remifentanil_eleveld();`

Eleveld DJ, Proost JH, Vereecke H, Absalom AR, Olofsen E, Vuyk J,
Struys MM. An allometric model of remifentanil pharmacokinetics and
pharmacodynamics. Anesthesiology. 2017 Jun 1;126(6):1005-18.

#### `pmx_model_remifentanil_minto();`

Minto CF, Schnider TW, Egan TD, Youngs E, Lemmens HJ, Gambus PL,
Billard V, Hoke JF, Moore KH, Hermann DJ, Muir KT. Influence of age
and gender on the pharmacokinetics and pharmacodynamics of
remifentanil: I. Model development. Anesthesiology. 1997 Jan 1;86(1):10-23.

Minto CF, Schnider TW, Shafer SL. Pharmacokinetics and pharmacodynamics
of remifentanil: II. Model application. Anesthesiology. 1997 Jan
1;86(1):24-33.

#### `pmx_model_sufentanil_gepts();`

Gepts E, Shafer SL, Camu F, Stanski DR, Woestenborghs R, Van Peer A,
Heykants JJ. Linearity of pharmacokinetics and model estimation of
sufentanil. Anesthesiology. 1995 Dec 1;83(6):1194-204. 

#### `pmx_model_remimazolam_eleveld();`

Eleveld DJ, Colin PJ, Van den Berg JP, Koomen JV, Stoehr T, Struys MM. 
Development and analysis of a remimazolam pharmacokinetics and 
pharmacodynamics model with proposed dosing and concentrations for 
anaesthesia and sedation. British Journal of Anaesthesia. 
2025 Jul 1;135(1):206-17. 

---

## $DIFFEQN

Defines the differential equations for ODE-based advans. Only used with `diffeqn_libgsl` or `diffeqn_test`.

C code called at each ODE integration step. The following are available:

| Symbol | Description |
|---|---|
| `A(i)` / `STATE(i)` | Current compartment amount i (read-only, 1-indexed) |
| `DADT(i)` | Rate of change of compartment i (write, 1-indexed) |
| `$DATA` parameters | Values from each data record |
| `$IMODEL` parameters | Values computed in `$IMODEL` for this individual |
| CSV column names | Data columns |

All `DADT(i)` from 1 to `.nstate` must be assigned.

```c
$DIFFEQN
    DADT(1) = -A(1) * KA;
    DADT(2) =  A(1) * KA - A(2) * (CL / V);
    DADT(3) =  E0 * KOUT * (1 - EMAX * (A(2)/V) / (C50 + A(2)/V)) - KOUT * A(3);
```

---

## $PREDICT

Declares prediction output variables and the C code to compute them per observation record.

The prediction output variables are used as C language struct members 
so they must be valid C language identifers.

### Arguments

```
$PREDICT(var1, var2, ...)
```

Variable names become fields of the generated `PREDICTVARS` struct and appear as columns in the `.yhat` output file. Names may be comma- or space-separated. `Y` does not need to be declared here — it is always available.

All `$PREDICT(...)` variables are accessible in `table(...);` within `$MAIN`.

### Content

C code executed once per observation record. The following are available read-only:

| Symbol | Description |
|---|---|
| `THETA(i)` | Population parameter i (1-indexed) |
| `ETA(i)` | Individual random effect i (1-indexed) |
| `A(i)` | Compartment amount i at observation time (1-indexed) |
| `ERR(i)` / `EPS(i)` | Residual error i (1-indexed); zero during prediction |
| `$IMODEL` parameters | Values computed per individual |
| CSV column names | Data columns (e.g. `DVID`, `TIME`) |

`Y` must be assigned to the predicted observation. It may be set conditionally (e.g. for multiple DVIDs). Set `Y = NAN` to indicate no prediction for a record.

```c
$PREDICT(IPRED, BPRED)
    IPRED = A(2) / V;
    BPRED = A(3);
    if (DVID == 1) 	
		Y = IPRED * (1 + ERR(1)) + ERR(2);
    if (DVID == 2) 
		Y = BPRED + ERR(3);
```

---

## $THETA

Define values for population-level fixed effects. All parameter
estimates are either bounded or fixed. Two syntaxes are accepted:

**C struct form:**

```c
{ lower, initial, upper, type }
```

| Field | Description |
|---|---|
| `lower` | Lower bound for optimization |
| `initial` | Initial value |
| `upper` | Upper bound for optimization |
| `type` | `ESTIMATE` or `FIXED` |

```c
$THETA
    {   1,   8,  20, ESTIMATE },   // C struct form
    {   0, 0.1,   1, ESTIMATE },
    {   0,   0,   0, FIXED    },
```

**NONMEM style:**

```
(lower, initial, upper)    // all three bounds
(initial FIXED)            // or (initial FIX)
```

Parameters marked as `FIXED` or `FIX` are held constant and not optimized.


```c
$THETA
    (1 8 20)     // NONMEM style: lower=1, initial=8, upper=20
    (0.1 FIXED)  // fixed at 0.1
```

---

## $OMEGA

Specifies a diagonal block of the inter-individual variability (IIV) covariance matrix.

**Zero or negative values are treated as fixed** (not optimized).

### Arguments

```
$OMEGA(val1, val2, ...)
```

### Content

Comma- or space-separated list of variance values. 
The number of values equals the dimension of this diagonal block. 
Off-diagonal covariances are zero (thus fixed). 


```c
$OMEGA(0.1, 0.1, 0.1, 0., 0.)
```

**Alternative: C struct form.** When the section body is not followed 
by `(`, the text is passed as a raw C struct initializer to the 
generated file. 

```c
$OMEGA
    { OMEGA_DIAGONAL, 3, { 0.1, 0.1, 0.1 } }
```

---

## $OMEGABLOCK

Specifies a full lower-triangular block of the omega matrix, allowing off-diagonal covariances.

**Zero or negative values on the diagonal are treated as fixed** (not optimized).

### Arguments

```
$OMEGABLOCK(val1, val2, val3, ...)
```

### Content

Lower-triangular elements listed row by row. For a block of dimension N, N*(N+1)/2 values are required. **Negative diagonal values are treated as fixed.**

```
Row 1: val_11
Row 2: val_21  val_22
Row 3: val_31  val_32  val_33
...
```

```c
$OMEGABLOCK(
    0.1,
    0.01, 0.1,
    0.01, 0.01, 0.1)
```
**Alternative: C struct form.** When the section body is not followed 
by `(`, the text is passed as a raw C struct initializer to the generated 
file. 

```c
$OMEGA
    { OMEGA_BLOCK, 3, { 0.1,
                        0.01, 0.1,
                        0.01, 0.01, 0.1 } }
```

---

## $OMEGASAME()

Duplicates a previously defined omega block by reproducing its structure.

### Arguments

```
$OMEGASAME(n)
```

Where `n` is the index of the previously defined omega block to replicate.

```c
$OMEGASAME(1)
```

---

## $SIGMA

Specifies initial variances for the residual error terms.

**Zero or negative values are treated as fixed.**

### Arguments

```
$SIGMA(val1, val2, ...)
```

### Content

Comma- or space-separated list of variance values. 
Each value corresponds to one `ERR(i)` or `EPS(i)` term in `$PREDICT`. 

The parentheses in the header are optional; values may follow `$SIGMA` directly:

```c
$SIGMA(1, 1, 100)   // parenthesised form
$SIGMA 1, 1         // bare form (same result)
```

---

## $FUNCTIONS

The body is arbitrary text written verbatim to the C-file just 
before the main function. 

The primary use is to define utility functions that may be 
called in `$MAIN`.

```
$FUNCTIONS
static void print_message(void)
{
	fputs("hello world\n", stdout);
}
```
---

## $MAIN

Contains arbitrary C code specifying what analysis to run. 

A global object of the name openpmx of the type `OPENPMX`. 

It is possible to change the number of parallel threads/processes
are used by setting `openpmx.nthread`. If set to `0` then the 
default is used which is equal to the number of CPU cores minus 1.
Setting a negative value reserves that number of CPU cores and
will use all others. So setting `openpmx.nthread = -1;` and `openpmx.nthread = 0;`
will both use the same amount of CPU cores.

Within the block content some functions and C-macros are available:

#### `simulate(...);`

Monte Carlo simulation: samples ETAs from OMEGA and residuals from SIGMA, 
computes predictions, and overwrites DV in the in-memory dataset. 

Configuration is passed in a `SIMCONFIG` object:

| Option | Type | Default | Description |
|---|---|---|---|
| `.seed` | `unsigned long` | `200501041406` | RNG seed |

```c
$MAIN
	simulate(.seed=1234);
```

#### `evaluate(...);`

Evaluate the objective function at the current parameters without optimizing the outer loop.
Only inner ETA optimization is run.

Configuration is passed in a `STAGE1CONFIG` object:

| Option | Type | Default | Description |
|---|---|---|---|
| `.gradient_step` | `double` | `1e-4` | Step size for gradient on first iteration |
| `.step_initial` | `double` | `0.2` | Initial trust-region radius |
| `.step_refine` | `double` | `0.1` | Refined trust-region radius |
| `.step_final` | `double` | ~`1.3e-5` | Final trust-region radius |
| `.maxeval` | `int` | `1000` | Maximum function evaluations |

```c
$MAIN
	evaluate(.maxeval=100);
```

#### `estimate(...);`

Run population parameter estimation with an outer- and inner-optimization.

Configuration is passed in a `ESTIMCONFIG` object:

| Option | Type | Default | Description |
|---|---|---|---|
| `.step_initial` | `double` | `0.2` | Initial BOBYQA trust-region radius |
| `.step_refine` | `double` | `0.05` | Refined trust-region radius |
| `.nsig` | `double` | `3.0` | Convergence criterion: require this many significant digits of stability in parameters |
| `.dobjfn` | `double` | `1e-3` | Convergence criterion on OFV change |
| `.maxeval` | `int` | `10000` | Maximum outer function evaluations |
| `.details` | `bool` | `false` | Print per-improvement model information |
| `.verbose` | `bool` | `false` | Print per-evaluation model information |
| `.stage1.*` | — | see `STAGE1CONFIG` | Inner (Stage 1) optimizer options |

Estimation optimization terminates when:

1. Parameter estimates to precision greater than `nsig`.
This is considerd successful convergence.
2. Change in objective function less than `dobjfn`. 
The estimation is considered to have stalled and terminates.

After estimation terminates:

1. The achieved `nsig` is estimated by modifying the encoded model by the calculated `step_final`.
2. Final model information is printed
3. Standard estimation tables are written.


```c
$MAIN
	estimate(.nsig=5);
```

#### `predict();`

Compute IPRED predictions (individual, using current ETA) at the 
current population parameters without estimation or simulation.

#### `predict_pred();`

Like `predict()` but computes PRED (population predictions with all ETA fixed to zero). 

Used to evaluate the structural model contribution.

#### `reload(...);`

Load population parameter estimates from a previous run's `.ext` file into the current model. 
Individual level data (`ETA` values) are not loaded.

Useful for warm-starting or chaining runs.

| Option | Type | Default | Description |
|---|---|---|---|
| `.filename` | `string` | `<controlfile>.ext` | Path to a `.ext` file from a previous run |
| `.force` | `bool` | `false` | Accept the file even if the parameter structure differs from the current model (truncates or pads as needed) |
| `.optional` | `bool` | `false` | If the file is absent or loading fails, continue silently rather than exiting |
| `.silent` | `bool` | `false` | Suppress the message showing the newly loaded parameter values |

```c
$MAIN
    reload();
    estimate();
```

#### `table("field, field ...", ...)`

Write a table of output columns to a file or stream.


| Option | Type | Default | Description |
|---|---|---|---|
| first argument | `string` | — | Comma-separated list of column names (required) |
| `.stream` | `FILE*` | — | Output FILE* stream (e.g. `stdout`) |
| `.filename` | `string` | — | Output file path |
| `.name` | `string` | — | Output to `<controlfile>.<name>`|
| `.firstonly` | `bool` | `false` | Write only the first record per individual |

If `.filename`, `.name`, `.stream` are missing then output is written to
a file `<controlfile>.table.XXX.txt` where `XXX` is incremented at each table generation.

```c
$MAIN
	estimate();
	table("ID,TIME,DV,IPRED", .filename = "output.csv")
	table("ID,TIME,IPRED", .name = "mytable")
```

---

## $FILE

Writes the body to a verbatim side-file generated by `openpmxtran`. 

The primary use is to provide a graphics configuration file for `graphics.R`.

### Arguments

```
$FILE("name")
```

The `name` is used as the filename suffix: the file is written as
`<controlfile>.<name>`. For example, with a control file named 
`schnider.gr` and `$FILE("graphics")`, the output is `schnider.gr.graphics`.

```
$FILE("graphics")
NAME    VALUE  TYPE     SETTING
DVID    1      name     Concentration
DVID    1      scale    log
DVID    2      name     Effect
```
