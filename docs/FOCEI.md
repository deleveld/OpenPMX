# FOCEI Implementation

## Objective function

The FOCEI objective function is −2 log-likelihood, assembled from five terms per
individual as defined in Bae and Yim (referenced in source comments). The full objective
is the sum over all individuals, assembled in `idata_objfn()` (`src/idata.c`), which
receives the population-level Term 4 (`omegainfo->omega_nonzero_lndet`) as an argument.

```c
iobjfn = term1 + term2 + term2_5 + term3 + term4 + term5
```

`term2_5` is `obs_logp`, the `−2·loglik` of discrete observations (see below), which
stands in for Terms 1+2 for those records. Each individual's `iobjfn` is stored, and the
terms are also summed by type across individuals before being combined, to add numbers of
similar magnitude and reduce floating-point error. For an individual with no
observations (`nobs == 0`) all terms are zero, including Term 4: the population term
does not count when an individual has no observations. Here `nobs` counts EVID=0 records
with a non-NaN DV.

### Observation types

Each EVID=0 record is either:

- a **continuous observation** — the predict function returns a non-NaN `Y`. It
  contributes to Terms 1 and 2 through `yhatvar`.
- a **discrete observation** — `Y` is NaN and the predict function returns a
  log-likelihood. It contributes `−2·loglik` (`obs_logp`) instead. `individual_fasteval()`
  returns `obs_min2ll + obs_lndet + obs_logp`, and this is what Stage 1 minimizes
  together with Term 3.

### Term 1 — observation variance log-determinant

`obs_lndet = Σ log(yhatvar)`

Summed over the continuous observation records (EVID=0, `Y` not NaN) of one individual.

### Term 2 — weighted squared residuals

`obs_min2ll = Σ (DV − yhat)² / yhatvar`

Summed over the same records as Term 1. Terms 1 and 2 are computed together in
`individual_fasteval()` and `individual_evaluate()` in `src/ievaluate.c`.
`individual_fasteval()` is used inside the Stage 1 optimization callback;
`individual_evaluate()` is used once at the final ETAs (which also fills `yhat`,
`yhatvar` and `obs_logp`) and for the Jacobian evaluations in Term 5. Both advance the
ADVAN over the individual's records and accumulate each sum separately, to avoid loss of
precision when their magnitudes differ strongly.

### Prediction variance (yhatvar)

`yhatvar` is computed in `evaluate_yhatvar()` (`src/ievaluate.c`) by central differences
around ERR = 0 for each SIGMA — the delta-method propagation of residual variance through
the user's predict function:

```
yhatvar = Σ_j  (f(ERR_j = +√σⱼ) − f(ERR_j = −√σⱼ))² / 4
```

The sum skips any `σⱼ = 0`. This simplifies from the standard delta-method formula
because the step size is `√σⱼ` and `(2√σⱼ)² = 4σⱼ` cancels with the `σⱼ` factor.
`yhatvar` is evaluated at the current ETA, and is set to zero in `individual_evaluate()`
for non-observation and discrete-observation records.

### Term 3 — ETA prior penalty

`eta_min2ll = η' Ω⁻¹ η`

The Mahalanobis distance of the individual's ETAs from zero under the population
distribution. Computed by `sample_min2ll()` (`src/nonzero.c`) on the reduced ETA vector,
by solving with the Cholesky factor `L` of the reduced Omega and taking `|L⁻¹η|²`
(`sample_min2ll_from_cholesky()` in `src/linalg.c`; an inverse-based alternative is
compiled out). Only the non-zero elements of Ω are included: the `NONZERO` structure
lists the rows/columns of Ω whose row/column sum is non-zero, and `reduce_eta()` /
`unreduce_eta()` convert between the reduced and the full ETA vector.

### Term 4 — population omega log-determinant

`omega_nonzero_lndet = log|Ω|`

The log-determinant of the reduced (non-zero) Omega, computed from its Cholesky
decomposition by `omegainfo_update_inverse_lndet()` (`src/omegainfo.c`) as
`2 Σ log(L_ii)` (`matrix_lndet_from_cholesky()`). The same call stores the Cholesky
factor and the inverse used in Terms 3 and 5. It is called by `encode_transform()` and
`encode_untransform()`, i.e. once per Stage 2 candidate. `idata_objfn()` adds it once for
every individual that has observations.

### Term 5 — individual information matrix log-determinant

`icov_lndet = log|J'J + Ω⁻¹|`

The Laplace approximation term. After Stage 1 converges to the conditional mode η̂, it is
computed by `stage1_individcov()` (`src/stage1.c`). There is a single implementation;
`stage1.c` has no options to select a different Hessian approximation. Despite the
name, `icov_lndet` is the log-determinant of the information matrix, i.e. of the
inverse of the individual covariance matrix `icov`.

1. Compute a numerical Jacobian `J` of the standardized residual `(yhat − DV) / √yhatvar`
   with respect to each active (non-zero Omega) ETA by central differences. Each ETA `j`
   is perturbed by `±step`, where `step = max(√icov[j,j], gradient_step)` and `icov` is
   the individual's stored covariance matrix from the previous evaluation (if its
   diagonal element is zero, `gradient_step` is used). The derivative is divided by
   `above − below` (= 2·step). Rows of `J` are the EVID=0 records of the individual
   (`nobs` rows). `yhatvar` is evaluated at each perturbed ETA, so the dependence of the
   residual variance on η (eta-sigma interaction) enters through `J`.
   For records where `yhatvar` is zero at either perturbed point (discrete
   observations), the derivative of the per-record `LOGP` (`−2·loglik`) is used instead.
2. Form `J'J + Ω⁻¹` (the Gauss-Newton approximation to the Hessian) using the reduced
   inverse Omega stored in `nonzero->inversedata`.
3. Cholesky-decompose `J'J + Ω⁻¹` (`gsl_linalg_cholesky_decomp1`). A copy of the matrix
   is kept first because a failed decomposition overwrites it. If the decomposition
   fails, `1e-6` is added to the diagonal of the copy and the decomposition is tried
   once more (the second attempt is not checked). `icov_lndet` is `2 Σ log(L_ii)` of
   the resulting Cholesky factor, i.e. the log-determinant of the information matrix
   before it is inverted.
4. Invert from the Cholesky factor to give the individual covariance matrix and write it
   back to the full-size `icov` of the individual. This `icov` is the one used for the
   step sizes in step 1 on the next evaluation.

In `stage1_thread()` this calculation is repeated up to 4 times when Stage 1 started from
all-zero ETAs (the step sizes depend on `icov`, which is not yet available), and once
otherwise. It stops early when `icov_lndet` changes by less than 0.01 between
repetitions. Individuals with no observations or no active ETAs have their ETAs and `icov`
set to zero and `icov_lndet = 0`; if they have observations, Terms 1 and 2 are still
evaluated at zero ETA.

---

## Hessian approximation

The Hessian used in Term 5 is `J'J + Ω⁻¹` as described above. `J` is the derivative of
`(f − DV) / √V` with `V` re-evaluated at each perturbed ETA, so the actual residuals
(through DV) and the η-dependence of V both enter `J`. It is a Gauss-Newton
approximation: second-derivative terms such as `(DV − f)·∂²f/∂η²` are dropped. 

### Role of the Hessian in optimization

The Hessian is not fed to the BOBYQA optimizer — BOBYQA is derivative-free and builds its
own internal quadratic surrogate from function evaluations. The Hessian enters only
through Term 5 (`icov_lndet`) of the objective value, and through the `icov` diagonal
that sets the finite-difference step sizes in the next Term 5 calculation. After Stage 1
has converged to η̂, `icov_lndet` is added to the FOCEI objective, so the choice of
Hessian approximation changes the objective landscape and therefore the optimum BOBYQA
converges to, but not the mechanism by which it gets there.

### FOCE vs FOCEI in OpenPMX

OpenPMX always computes FOCEI — there is no FOCE mode (the out file reports
`optim FOCEI BOBYQA`). `yhatvar` is always evaluated at the current η, so the interaction
between random effects and residual error is always included in the objective function.

---

## Optimization structure

Estimation uses two nested BOBYQA optimization loops. Both use `npt = 2n+1`
interpolation points and bounds of `±1e6` (the `±DBL_MAX` bounds are commented out in
the source).

### Stage 1 — inner loop (per individual)

**File:** `src/stage1.c` — `estimate_individual_posthoc_eta()`, called from
`stage1_thread()`

Optimizes the individual ETA values (reduced to the non-zero Omega elements) for fixed
population parameters. The callback `stage1_evaluate_individual_iobjfn()` evaluates
terms 1+2 (and the discrete-observation term) via `individual_fasteval()`, plus term 3.
Terms 4 and 5 are not part of Stage 1 — term 4 is population-level, term 5 is computed
after Stage 1 has converged.

Two passes are possible, controlled by the stage 1 options:

1. **Initial pass** — `rho` from `step_initial` to `step_refine`. Done when the starting
   ETAs are all zero, or when `stage1_evaluate` is set in the scatter options (the very
   first objective function evaluation of an estimation or evaluation, and the
   evaluations done for the covariance step).
2. **Refinement pass** — `rho` from `step_refine` to `step_final`. Always done, starting
   from the result of the initial pass or from the previous ETAs when warm-starting.

The evaluation limit per pass is the stage 1 `maxeval`. A BOBYQA failure is reported as
a warning with the individual's ID and does not abort the estimation.

### Parallelisation (`src/scatter.c`)

Stage 1 runs independently per individual, so `scatter_threads()` distributes
`stage1_thread()` over individuals using pthreads, OpenMP, or a single thread depending
on the build.

- The number of threads is limited to the number of individuals. The main thread also
  does work, so `nthread − 1` workers are created.
- With no workers, individuals are processed serially in index order.
- Otherwise individuals are sorted by their previous run time (`stage1_msec` when
  `stage1_order` is set, otherwise `eval_msec`) and the slowest are processed first, to
  avoid a long tail at the end.
- The pthreads pool persists between calls and is cleaned up at program exit; OpenMP uses
  a dynamic schedule over the sorted order.
- Output is serialized with `print_serialize()` during the parallel section.

### Stage 2 — outer loop (population)

**File:** `src/estimate.c` — `focei()`, with the objective callback `stage2_objfn()`

Optimizes the population parameters (THETA, OMEGA, SIGMA). The parameters are held in an
*encoded* vector (see below). The callback `stage2_objfn()`:

1. Decodes the vector into a candidate `POPMODEL` via `encode_untransform()`, which also
   updates the Omega Cholesky factor, inverse and log-determinant (Term 4)
2. Resets all individual ETAs and `icov` matrices to the best known values
   (`idata_etas_write()`)
3. Calls `encode_evaluate()`, which dispatches `stage1_thread()` across all individuals
   via `scatter_threads()` and assembles the full objective function with
   `idata_objfn()`
4. Updates the best model so far (`update_best_imodel()`)

Only a finite objective function is accepted as the best model. The best ETAs and `icov`
matrices are saved each time the best model improves and used as starting values for
later evaluations.

`focei()` runs the following phases:

1. **Initial evaluation** (`test_initial_objfn()`) — the model as given is evaluated once
   with `stage1_evaluate` set, giving the starting objective function and ETAs. A
   warning is given if it is not finite.
2. **Coarse search** — BOBYQA with `rho` from `step_initial` to `step_refine`.
3. **Final step size** — the final `rho` is not an option. Starting at `step_refine / 2`,
   it is halved until perturbing each encoded parameter by `±rho` changes the decoded
   parameters in more than `nsig` significant digits (smaller change than `nsig`
   digits), or until it falls below `1e-8`.
4. **Fine refinement** (`refine()`) — BOBYQA restarted from the best model with `rho`
   from `step_refine` to the final step size. The individual ETAs and `icov` matrices are
   reset to zero at each restart (so Stage 1 does a full initial pass again). After each restart the number of significant digits shared between the
   parameters before and after the restart (minimum over THETA, SIGMA and the lower
   triangle of OMEGA) is reported. Restarts continue until one of:
   - the shared significant digits exceed the `nsig` option (converged),
   - `|Δobjfn| < dobjfn` between restarts (stalled),
   - BOBYQA reports an error, or
   - the evaluations allowed by `maxeval` are used up.

The result is marked final (`OBJFN_FINAL`) only if the refinement converged; otherwise
the best model found so far is returned. `pmx_evaluate()` is the same procedure with
`maxeval = 1`, so only the initial evaluation is done.

### Parameter encoding (`src/encode.c`)

Population parameters are not optimized directly. Only non-fixed parameters are encoded,
in the order THETA, SIGMA, then the lower triangle of the non-fixed OMEGA (rows and
columns that are non-zero and not `OMEGA_SAME`). The vector `x` is unbounded and is
relative to an offset vector set by `encode_transform()`: the parameter is decoded from
`x + offset`, so an all-zero `x` reproduces the current model. Step sizes (`rho`) are in
encoded units.

- **THETA**: mapped from its `[lower, upper]` bounds. The mapping is chosen at compile
  time: `ENCODE_SPECIAL` (default), `ENCODE_ATANH`, `ENCODE_LINEAR` or `ENCODE_NONMEM`
  (`log(θ − lower) − log(upper − θ)`). The default scales θ to `[−1, 1]` and is linear
  over the central 90% of the range, with atanh-shaped tails near the bounds so the
  encoded value is unbounded.
- **SIGMA**: log-transformed. `ESTIMATE_SD_NOT_VAR` is defined, so the log of the
  standard deviation is encoded (`σ = exp(x + offset)²`).
- **OMEGA diagonal**: log of the standard deviation, as for SIGMA.
- **OMEGA off-diagonals**: encoded through the Cholesky factor of the correlation
  matrix as canonical partial correlations (atanh-transformed), following the Stan
  correlation-matrix transform. OMEGA is rebuilt as the correlation matrix scaled to the
  decoded variances, so it stays positive definite.
- `OMEGA_SAME` blocks are filled in by copying the earlier block after each update.

If there are no estimated parameters, no optimization is done after the initial
evaluation.

---

## BOBYQA algorithm

BOBYQA (Bound Optimization BY Quadratic Approximation, Powell 2009) is a derivative-free
trust-region method. It maintains a quadratic surrogate model built from
`npt = 2n+1` interpolation points (the recommended setting; `n+2` and
`(n+1)(n+2)/2` are the minimum and maximum). At each iteration it
minimizes the surrogate within a trust region of radius `rho`, evaluates the true
objective at the suggested point, and updates the model. `rho` shrinks from `rhobeg` to
`rhoend` as the search converges.

Because it requires only function values, it is well-suited to the noisy, non-smooth
FOCEI objective where numerical derivatives would be unreliable.

