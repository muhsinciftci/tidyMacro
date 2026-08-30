# LP-DiD review

## Scope and conclusion

Reviewed the LP-DiD Stata replications under `inst/LPDID/codes/`, their
called preparation scripts, `src/fLPDID.cpp`, `src/fLPDID.h`,
`src/panel_shift.h`, and `R/fLPDID.R`.  There are no MATLAB (`.m`) sources
under `inst/LPDID/`, nor is there a `src/fLPDID.cpp.R`; the relevant R
counterpart is `R/fLPDID.R`.

The C++ engine gets the core construction right: gap-aware shifts, the
absorbing and DDCG transition samples, the all-past PMD baseline used in the
main simulation, weighted time demeaning, and the CR1 score formula are
coherent.  A direct check of the DDCG ANRR/no-CCC impact regression against
`fixest::feols(D0y ~ tdemoc + lag1y + ... + lag4y | year, cluster = ~wbcode2,
ssc = fixest::ssc(fixef.K = "full"))` produced the same coefficient and standard error to
about `9e-15` (C++: `0.3603594`, `0.6504521`, `N = 3029`, `G = 126`).

There is one high-priority interval-construction bug, a smaller
`reghdfe`-compatibility discrepancy, and input-validation problems that can
silently produce a wrong analysis.

## Findings

### High — clustered confidence intervals use a normal, not cluster-t, critical value

The C++ covariance uses a CR1 factor with the number of clusters, but the R
wrapper builds intervals with `qnorm` (`R/fLPDID.R:367-371`).  The documented
and Stata-compatible cluster inference policy is a t critical value with
`G - 1` degrees of freedom (or the appropriate finite-sample residual df when
that is the adopted convention).

For a 20-cluster test, the wrapper multiplier was `1.959964`, while
`qt(.975, 19)` is `2.093024`.  Point estimates and SEs are unaffected, but the
published `conf_low` and `conf_high` are too narrow.  This is the highest
priority finding because those bands are what users read from `fPlotLPDID()`.
The normalized `-1` point can remain exactly zero.

### Medium — singleton time-FE observations are retained

`reghdfe` drops singleton absorbed-FE observations by default.  The C++ code
selects all usable rows, demeans them by time, and keeps a time cell containing
only one row (`src/fLPDID.cpp:210-251`, `321-336`).  Such a row becomes zero
after demeaning, so it cannot change the coefficient, cluster score, or
`denom = m - (Pk + nY)`: adding a singleton increases both `m` and `nY` by
one.  It does, however, change the `(m - 1)` term in the CR1 adjustment and
the reported `nobs` (`391-400`).

This is observable in the Figure 3 banking data, impact horizon, only-FE
interstate-banking specification:

| estimator | estimate | clustered SE | N |
| --- | ---: | ---: | ---: |
| `fLPDID` | -0.004862917 | 0.001785323 | 753 |
| `fixest`, `fixef.K = "full"` | -0.004862917 | 0.001784135 | 752 |

`fixest` reports one singleton time FE removed.  The SE difference is small in
this example and is bounded by the CR1 numerator adjustment: with `s`
retained singleton rows, the C++/singleton-dropped SE ratio is
`sqrt((m - 1) / (m - 1 - s))`.  The issue is therefore a medium-severity
inference and reporting discrepancy, not a point-estimate problem.  It is
also confined to `reweight = FALSE`, because a singleton time cell has
`p_t` equal to zero or one and is already excluded by the undefined implicit
weight.  Drop singleton selected time cells to match default `reghdfe` and
report the post-drop `nobs`.  With only time FEs, one pass is sufficient.

### Medium — transformed LHS terms are silently ignored

The wrapper derives the response name with `all.vars()` and then passes the
raw column to C++ (`R/fLPDID.R:278-280`, `335-339`).  Consequently,
`fLPDID(log(y) ~ treat, ...)` estimates the model for `y`, not `log(y)`, with
no warning.  A deterministic test returned identical estimates and SEs for
`y ~ tr` and `log(y) ~ tr`.

The right fix is to evaluate the LHS expression in `data`/the formula
environment (for example, `eval(main_formula[[2L]], data, env)`) and pass the
resulting numeric vector to C++.  RHS transformations already work through
`model.matrix`, so rejecting transformed LHS terms would leave an unnecessary
and inconsistent interface restriction.

### Medium — treatment, panel-id, time, and cluster inputs lack required checks

The documentation says treatment is binary (`R/fLPDID.R:98-100`), but the
wrapper only verifies numeric type.  The engine treats only `dD == 1` as an
entry and `dD == 0` as a candidate control (`src/fLPDID.cpp:217-239`).  A
treatment coded `0/2` therefore completes without an error and returns `NaN`
for every estimable horizon rather than explaining the invalid coding.

Missing `id` and `cluster` values are also passed through
`as.integer(as.factor(...))` (`R/fLPDID.R:340-342`).  The C++ hash maps then
treat each missing integer code as an ordinary shared unit/cluster.  The
duplicate check happens to catch two missing ids at the same time, but one
missing id (or missing ids with non-overlapping times) passes and creates a
synthetic unit.  In a 20-cluster test, setting three units' cluster id to `NA`
reduced the returned count to 18 and changed the SE from `0.16872` to
`0.16645`.

Missing time is worse: `sort(unique(t_raw))` drops `NA`, the remaining time
step is computed normally, and `as.integer(...)` passes an `NA_integer_` to
C++.  The hash maps treat it as an ordinary time period.  A direct test with
one missing time returned a finite, wrong result rather than an error.  Reject
or explicitly remove missing id/time/cluster values before time normalization
and factor conversion.

Validate before any factor conversion:

- `treat[!is.na(treat)] %in% c(0, 1)`;
- finite, non-missing unit id, time, and cluster values (or explicitly remove
  incomplete identifier rows and say so);
- finite numeric time with a clear diagnostic before computing the step.

### Medium — the implementation does not reproduce every bundled replication variant

These are mostly well-known feature gaps rather than accidental errors, but
they matter when describing the implementation as a replication of the whole
`inst/LPDID` folder.

- The main-text simulations use the all-past PMD baseline, which C++ matches
  (`first_sim_pmd_LPDiD_estimation.do:14-25`; `src/fLPDID.cpp:97-135`).  The
  appendix instead uses a moving average over the prior `L` outcomes and
  removes baselines contaminated by treatment changes
  (`appendix_sim_pmd_LPDiD_estimation.do:18-26`).  The C++ `pmd = TRUE` option
  cannot reproduce that appendix estimator.
- The appendix non-absorbing scripts use `D.treat == 0` as their comparison
  state, which includes stable treated spells
  (`appendix_sim_LPDiD_estimation.do:57`).  The C++ non-absorbing controls
  require both `dD == 0` and treatment level zero
  (`src/fLPDID.cpp:235`).  The latter is correct for the DDCG `tdemoc`
  construction (`ddcg_prepare_dataset.do:12-14`) and Figure 4's CCC1/CCC2
  regressions, but it is not a literal port of the appendix sample.
- DDCG and the reweighted appendix scripts use regression adjustment/ATET
  (`teffects ra`, or interacted regression followed by `margins`), not the
  inverse-implicit-weight WLS used by C++ (`figure_4.do:157-166`,
  `appendix_sim_rwLPDiD_estimation.do:68-91`).  The R warning correctly
  discloses this, so this is a limitation rather than a hidden mismatch.
- Figure 3 specification B has a horizon-varying set of leads of the other
  treatment (`figure_3.do:88-107`).  A single `fLPDID()` call builds one fixed
  control matrix, so it cannot reproduce that specification exactly.

The DDCG CCC definitions themselves agree with the current C++ design:
CCC1 uses `CCS_0` for both entry and controls; CCC2 uses `CCS_0` for entries
and `CCS_h` for post-horizon controls; both use `CCS_mh` on the pre side.

### Low — diagnostics are discarded or incomplete

The C++ result records `ndrop` for collinear controls, but the R result drops
it (`src/fLPDID.cpp:430-446`; `R/fLPDID.R:354-373`).  It also returns all-NaN
horizons with `nobs = nclust = 0` without a wrapper warning.  Retaining
`ndrop` and warning with horizon/sample information would make failed or weak
cells much easier to diagnose.

## Performance review

The important architecture is already good: LP-DiD precomputes shifts and
clean-control sets, avoids R calls inside OpenMP, compacts time/cluster ids,
and parallelizes independent horizons.  For typical small control sets, the
highest-value remaining improvements are:

1. **LP-DiD: avoid the redundant cross-product and full covariance inverse.**
   `A = Zs'Zs` is formed for rank detection, then `Ak = Zsk'Zsk` is formed
   again (`src/fLPDID.cpp:344`, `376-382`).  Use `A` directly when no columns
   were dropped, or its kept-column submatrix otherwise.  Only the treatment
   coefficient and its variance are returned, so solve for the treatment
   column of `Ak^{-1}` and compute `q' meat q`; do not form `Ainv` or the full
   `V` matrix.  This is exact and matters as the control count rises.
2. **LP-DiD: reuse per-thread workspace.**  Each horizon reallocates
   remapping vectors, sample vectors, design matrices, and cluster-score
   matrices.  Thread-local reusable buffers would reduce allocator traffic on
   long panels and many horizons.  Preserve the current no-R-API-in-OpenMP
   rule.
3. **Panel LP: remove repeated hash lookups and avoid unconditional SVDs.**
   `fLPPanel` shifts every `sX` column separately for each lag
   (`src/fLPPanel.cpp:203-211`), builds a new time hash map every horizon
   (`312-319`), and computes `pinv(X'X)` unconditionally (`321-322`).
   Cache row indices for every required lag and gather all `sX` columns in one
   pass; use the existing global compact-time/scratch-remap pattern from
   LP-DiD; use Cholesky/QR solves with `pinv` only as a rank-deficient
   fallback.  These are the clearest general Panel-LP speedups.
4. **Panel LP: control peak memory.**  `sX_lags` is an
   `n_obs × n_s × p_max` cube and each horizon also materializes `W_lag`, its
   complete-case mask, and selected matrices.  This can dominate runtime via
   memory bandwidth.  A compact lag-row index cache plus horizon-local gathers
   trades some recomputation for far lower peak memory, which is preferable
   for wide interactions or large panels.
5. **Classical LP and LP-IV: batch HAC work across response columns.**
   Both engines already use QR and `fLP` has a good shock-only fast path.
   Their remaining hot loop is scalar Bartlett autocovariance accumulation by
   response and lag.  For multi-outcome fits, form the projected score matrix
   once and use BLAS-level cross-products/columnwise diagonals for each lag.
   For LP-IV with controls, residualize `[Y, Z, D]` together against the same
   thin-Q matrix rather than multiplying by the residual maker three times.
   These changes matter much more for many outcomes or long horizons than for
   the usual small VAR-sized use case.

Do not optimize by forming dense time or fixed-effect dummies, or by allowing
Armadillo to create nested OpenMP work inside the horizon loop; both would
work against the current efficient design.

## Verification performed

No package code was changed and no Stata/MATLAB execution was performed.
Read-only checks used the installed build corresponding to this clean source
tree plus the bundled processed data:

- DDCG impact, `ccc = 0`, compared with `fixest`: coefficient and clustered
  SE match to machine precision.
- Banking impact, only time FEs, compared with `fixest`: coefficient matches,
  while the singleton treatment above exposes the SE/N discrepancy.
- Synthetic wrapper tests confirmed transformed-LHS ignoring, normal rather
  than cluster-t intervals, silent non-binary-treatment failure, and silent
  missing-cluster inclusion.

Recommended regression tests after fixes: the two empirical comparisons above,
a one-row time cell, `log(y)` evaluation, `0/2` treatment rejection,
and missing id/time/cluster rejection.
