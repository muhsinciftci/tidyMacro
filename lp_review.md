# LP and LP-IV Review

## Scope and conclusion

Reviewed only the classical time-series LP paths:

- JT2025 replication driver: `inst/Bianchi/Replic/JT2025/GO_JT2025.m`
- Toolbox implementation called by that driver: `inst/Bianchi/VAR/LPmodel.m`,
  `LPoption.m`, `OLSmodel.m`, `VARmakexy.m`, and `VARmakelags.m`
- tidyMacro interfaces and engines: `R/fLP.R`, `R/fLPIV.R`, `src/fLP.cpp`,
  `src/fLPIV.cpp`, and their headers

I did not review panel LP, LP-DID, or generated Rcpp export files.

The numerical kernels are sound for the JT2025 estimators. With correctly
aligned inputs, the C++ kernels reproduce a direct implementation of the
toolbox equations to machine precision. The important defect is in the
high-level R wrappers: their global complete-case filtering removes the
pre-sample outcome required by a cumulative/long-difference LP. As a result,
the documented R equivalents do not replicate either JT2025 figure.

## Findings

### High: cumulative LP preprocessing drops one valid observation

`fLP()` expands lag terms and then removes every row containing an `NA`
(`R/fLP.R:451-477`). For the JT2025 OLS specification, `l(..., 1:4)` therefore
removes the first four rows. The C++ cumulative branch then independently
drops its first input row in order to form `y_t - y_{t-1}`
(`src/fLP.cpp:112-123`).

This is one drop too many. The toolbox first creates the control lags, retains
the outcome at the preceding date, and uses it as the long-difference base:
`endo_lag1 = ENDO(lag_trim:end-1)` (`inst/Bianchi/VAR/LPmodel.m:333-335`),
then pairs it with the first usable RHS date (`LPmodel.m:371-387`).

Concrete JT2025 result at horizon zero:

| Estimator | Raw sample | Correct h=0 rows | Current h=0 rows | Coefficient |
| --- | ---: | ---: | ---: | ---: |
| LP-OLS, CPI response | 92 | 88 | 87 | toolbox-aligned `0.039260334371`; current `fLP()` `0.025611683034` |
| LP-IV, unemployment response | 181 | 175 | 174 | toolbox-aligned `-0.100596626234`; current `fLPIV()` `-0.081772319501` |

The current R results exactly equal a reference calculated after deliberately
discarding that first valid observation, so this is an alignment defect rather
than a tolerance or HAC issue. It applies at every horizon whenever
`cumulative = TRUE` and preprocessing has removed leading rows, including the
JT2025 LP-OLS and LP-IV calls.

Recommended resolution: retain one valid LHS pre-sample row separately from
the estimable RHS rows and pass it to the cumulative engine as its first `Y`
row. Its corresponding `X`/`D`/`Z`/`C` row is not used by the C++ cumulative
loop. Do not try to fix this by changing the C++ alignment; that alignment is
correct for a kernel input that includes the required base outcome.

### High: global NA removal compresses time series with interior gaps

Both wrappers delete all incomplete rows before invoking C++:
`R/fLP.R:455-463` and `R/fLPIV.R:166-171`. The C++ engines subsequently treat
the retained rows as consecutive time periods (`src/fLP.cpp:117-122` and
`src/fLPIV.cpp:107-116`). An interior missing value therefore makes, for
example, observations at `t-1` and `t+1` adjacent in the projection. Leads,
long differences, and HAC autocovariances then use the wrong time spacing.

This does not affect the supplied JT2025 samples after their intended leading
lag trim, but it is unsafe for ordinary time-series data with internal missing
observations. The wrappers should either reject internal gaps, retain an
explicit time index and select contiguous horizon-specific samples, or split
the data into contiguous segments. A blanket `na.omit`-style compaction is not
valid for a time-series LP.

### Medium: rank-deficient designs silently produce inference that is not defined

The C++ engines fall back to a Moore-Penrose inverse for singular `X'X`,
`C'C`, or `Z'Z` matrices (`src/fLP.cpp:142-146`; `src/fLPIV.cpp:135-149`).
That avoids a crash, but no rank diagnostic is returned. In LP-IV, the
reported first-stage F still uses the raw column counts `kc` and `nz`
(`src/fLPIV.cpp:153-167`), which are wrong when controls or instruments are
collinear. In LP-OLS, a pseudoinverse can return a coefficient for a
nonidentified shock while its HAC standard error has no usual interpretation.

For a user-facing estimator, detect deficient rank and stop with an actionable
error (or explicitly report effective rank and suppress affected inference).
The JT2025 designs are full rank, so this does not explain the replication
discrepancy.

### Low: two non-estimation mismatches

- The toolbox's LP-IV F statistic uses `n_h - k_z` degrees of freedom after
  FWL residualisation (`inst/Bianchi/VAR/LPmodel.m:474-480`). The C++ engine
  uses `T_h - k_c - k_z` (`src/fLPIV.cpp:153-164`). The C++ denominator is the
  correct residual degrees of freedom for the original full first-stage
  regression. Thus this is an intentional, statistically correct divergence
  from the MATLAB *reported F statistic*, not an IRF or confidence-band error.
- `GO_JT2025.m:280` closes `\textbf{...}` before the LP-IV shock label, unlike
  the OLS title at line 261. The rendered title is valid but only its first
  phrase is bold. It has no effect on estimation.

## Replication mapping

| JT2025 path | Required tidyMacro specification | Notes |
| --- | --- | --- |
| LP-OLS (`GO_JT2025.m:85-104`) | `fLP(lcpi ~ rr_shock + l(dlrgdp,1:4) + l(dlcpi,1:4) + l(dstir,1:4), horizons = 17, shock = "rr_shock", conf = c(.68,.95), nw_lags = 0, nw_offset = 0, cumulative = TRUE)` | MATLAB's 18 `nsteps` are horizons 0 through 17; `OLSmodel(..., hh-1)` means HAC bandwidth `h`. |
| LP-IV (`GO_JT2025.m:179-202`) | `fLPIV(urate ~ ffr + l(urate,1:6) + l(infl,1:6) + l(ffr,1:6), instruments = ~ l(RRCGShock,0:6), endog = "ffr", horizons = 48, conf = c(.68,.95), nw_lags_iv = 6, cumulative = TRUE)` | MATLAB's 49 `nsteps` are horizons 0 through 48; instrument lags and fixed HAC bandwidth are 6. |

The first specification must not use the default `fLP()` HAC settings. The
default uses a data-dependent base bandwidth and an offset of one
(`R/fLP.R:482-489`); the code correctly documents that `nw_lags = 0` and
`nw_offset = 0` are needed for the toolbox convention (`R/fLP.R:292-300`).

## C++ validation

I read the toolbox call chain as follows:

`GO_JT2025.m` -> `LPmodel.m` -> `VARmakexy.m` / `VARmakelags.m` ->
`OLSmodel.m`; in IV mode `LPmodel.m` instead does horizon-specific FWL,
first-stage projection, and Newey-West delta-method inference
(`LPmodel.m:427-522`).

For the supplied workbook, I passed each C++ kernel the same inputs the
toolbox uses, including a padded pre-sample outcome row. The largest absolute
differences over all horizons were:

| Engine | IRF | HAC SE | First-stage F |
| --- | ---: | ---: | ---: |
| `fLP_cpp`, h = 0..17 | `4.16e-13` | `9.51e-14` | n/a |
| `fLPIV_cpp`, h = 0..48 | `1.57e-14` | `4.58e-12` | `7.99e-14` |

This validates the following C++ choices for this replication:

- Long-difference alignment in the kernels is correct when their input begins
  with the prior outcome (`src/fLP.cpp:115-118`, `src/fLPIV.cpp:107-111`).
- The OLS Bartlett/Newey-West sandwich matches `OLSmodel.m`
  (`src/fLP.cpp:178-223`; `inst/Bianchi/VAR/OLSmodel.m:123-148`).
- LP-IV FWL, 2SLS point estimates, and the HAC delta-method standard error
  match the toolbox construction (`src/fLPIV.cpp:120-217`).

## Performance opportunities

These are improvements, not correctness blockers. Profile before implementing
them on the expected workload.

1. **Make the `fLP` fast path actually allocation-light.** With
   `store_full = FALSE`, the code still constructs `Xreg`, `XU`, and `S` for
   every outcome and horizon (`src/fLP.cpp:164-191`). The shock-only variance
   can calculate `a_t = u_t * (v_0 + x_t' v_x)` directly from `Xh`, `u`, and
   the shock row of `(X'X)^{-1}`. This removes the `Xreg`, `XU`, and
   transposed-score `S` full-size intermediates from the common path without changing the HAC
   estimator.

2. **Batch LP-IV residualisation.** For `nc > 0`, `fLPIV` performs three
   separate FWL residualisations of `Yh`, `Zh`, and `Dh`
   (`src/fLPIV.cpp:132-140`). Residualising a single concatenated
   `[Yh | Zh | Dh]` matrix shares the two large matrix multiplications and is
   especially valuable for multiple outcome columns.

3. **Avoid horizon-slice copies where possible.** Each horizon allocates
   copies of aligned `Y`, `X` (LP) or `Y`, `D`, `Z`, `C` (LP-IV)
   (`src/fLP.cpp:113-123`; `src/fLPIV.cpp:101-117`). Standard LP can use
   Armadillo subviews directly; cumulative LP necessarily materialises the
   difference but can still avoid copying the RHS blocks. This lowers memory
   bandwidth and reduces per-thread peak memory.

4. **Prefer factorizations/solves to explicit normal-equation inverses.**
   `inv_sympd(X'X)` and `inv_sympd(Z'Z)` are fast for well-conditioned small
   designs, but square the condition number and force pseudoinverse fallbacks.
   QR-based solves or a cached Cholesky factor would improve numerical
   robustness and can be faster for larger regressor/instrument sets. This is
   a stability improvement first; benchmark it before replacing the current
   short path.

5. **Clarify thread behavior before tuning it.** Documentation says zero
   threads means "all cores minus one" (`R/fLP.R:307-308`, `R/fLPIV.R:55-56`),
   while both engines use all OpenMP threads available
   (`src/fLP.cpp:89-94`, `src/fLPIV.cpp:80-84`). Correct the documentation or
   implementation. On multi-threaded BLAS builds, all-core outer OpenMP can
   oversubscribe the machine, so benchmarks should test both `n_threads = 1`
   and a capped thread count.

## Verification and test gaps

- `R CMD INSTALL` succeeded in a temporary library.
- Actual-data comparisons above used `JT2025_Data.xlsx` and direct reference
  calculations; no source code was changed.
- Existing LP smoke coverage verifies baseline `fLP()` and the corrected
  full-model LP-IV F degrees of freedom (`tests/smoke_lp.R:14-76`). It does
  not exercise `cumulative = TRUE` together with generated lags, so it cannot
  catch the JT2025 alignment failure.

Add focused tests for: (1) OLS and IV long-difference calls with generated
lags, including the first valid observation; (2) internal missing-time gaps;
and (3) collinear controls/instruments.
