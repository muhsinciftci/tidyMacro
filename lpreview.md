# Local Projections C++ Review

Scope: `src/fLP.cpp`, `src/fLPIV.cpp`, `src/fLPPanel.cpp`, `src/fLPDID.cpp`, their headers, the handwritten wrappers `R/fLP.R`, `R/fLPIV.R`, `R/fLPPanel.R`, `R/fLPDID.R`, and the stated translation sources `inst/PanelLP/R_Code/panel_LP.R` and selected LP-DiD Stata scripts under `inst/LPDID/codes`.

I did not inspect unrelated VAR C++ files or generated documentation.

## Findings

### 1. Panel formula lags/leads do not preserve time gaps, while the C++ engines and source code do

Files:

- `R/fLPPanel.R:47-79`
- `R/fLPDID.R:40-60`
- `R/fLPDID.R:234-236`
- `src/fLPPanel.cpp:32-50`
- `inst/PanelLP/R_Code/panel_LP.R:15-35`

`src/fLPPanel.cpp` and the source `panel_LP.R` both implement shifts by looking for the same unit at exact time `t + L`. That means if a unit has observations at times `1` and `3` but not `2`, `l(x, 1)` at time `3` should be missing.

The R-side formula expander used by `fLPPanel()` sorts by `(id, time)` and shifts by row position within unit. It only checks that the shifted row is in the same unit (`R/fLPPanel.R:66-77`), not that time differs by exactly `k`. As a result, user-written formula controls such as `l(y, 1)` and `f(x, 1)` can silently bridge gaps in unbalanced panels, even though the C++ engine's internally generated `p_max` lags do not.

`fLPDID()` reuses this expander for controls and also has a separate `.flpdid_shift()` for `d(x)`. `.flpdid_shift()` checks raw time increments of exactly `k`, but it runs before time normalization. That breaks for `Date`, `POSIXct`, or numeric time indexes whose natural period step is not `1`. The C++ LP-DiD engine uses normalized integer time, so the wrapper and engine can disagree.

Impact:

- Incorrect control lags/leads in unbalanced panels.
- `d(x)` and `l(d(x), ...)` can become mostly missing for Date/POSIX time indexes.
- Translation mismatch against both `panel_LP.R:28-32` and the C++ `panel_time_shift()` implementation.

Suggested fix:

- Normalize time before all panel-aware R formula expansion.
- Use a gap-aware shift helper for `fLPPanel()` and `fLPDID()` that looks up `(id, normalized_time + k)`, not adjacent rows.
- Add a duplicate `(id, time)` check before building shifts.

### 2. `fLPPanel` calls R math and throws Rcpp errors inside an OpenMP loop

Files:

- `src/fLPPanel.cpp:263-265`
- `src/fLPPanel.cpp:294`
- `src/fLPPanel.cpp:404-416`

The horizon loop is parallelized with OpenMP, but it calls `Rcpp::stop()` when a horizon has no complete observations and calls `R::qt()` / `R::pt()` for confidence intervals and p-values inside the parallel region.

This violates the project profile's rule to avoid R API calls inside parallel regions. Even when it appears to work, R's C API and error machinery are not a safe dependency inside OpenMP workers.

Suggested fix:

- In the parallel loop, set per-horizon status flags instead of calling `Rcpp::stop()`.
- Compute only estimates, SEs, df, and status in C++.
- Move `qt`, `pt`, and confidence-band construction after the parallel loop or leave it entirely in the R wrapper, as `R/fLPPanel.R:472-500` already rebuilds user-requested bands.

### 3. LP-IV first-stage F statistic uses the wrong residual degrees of freedom after FWL residualization

Files:

- `src/fLPIV.cpp:115-129`
- `src/fLPIV.cpp:141-149`
- `src/fLPIV.cpp:56`
- `R/fLPIV.R:199-203`

The IV engine residualizes `D` and `Z` on `[1, C]`, then computes the first-stage F statistic as:

```cpp
((rss_res - rss_unr) / nz) / (rss_unr / (Th - nz))
```

After partialling out `kc = 1 + n_controls` controls, the unrestricted first-stage residual degrees of freedom should account for both the controls and the instruments, i.e. roughly `Th - kc - nz` under full rank. Using `Th - nz` overstates the denominator degrees of freedom and can inflate the reported first-stage F statistic, especially with many controls.

The C++ largest-horizon observation check also only tests against `kc + 1`, not `kc + nz`, while the R wrapper's `min_obs` does include `ncol(Z)`.

Impact:

- Point estimates are not affected.
- `Fstat_fs` can be materially wrong and too optimistic.
- Direct low-level calls can pass weakly identified or underdetermined first stages that the wrapper would reject.

Suggested fix:

- Use `df2 = Th - rank([1, C]) - rank(Zr)` or at least `Th - kc - nz`.
- Return `NA` or `NaN` when the denominator degrees of freedom are nonpositive.
- Align the C++ input check with the wrapper's stricter dimension requirement.

### 4. LP-DiD PMD baseline differs from the Stata source on unbalanced or missing panels

Files:

- `src/fLPDID.cpp:91-101`
- `inst/LPDID/codes/main_text_simulation/first_sim_pmd_LPDiD_estimation.do:14-25`
- `inst/LPDID/codes/main_text_simulation/first_sim_rw_pmd_LPDiD_estimation.do:14-25`

The Stata PMD baseline is:

```stata
bysort stateid (year) : gen cumulative_y = sum(ln_y)
gen time = year - start_year + 1
gen aveLY = L.cumulative_y/(time-1)
```

The C++ baseline is a running mean over finite past observed outcomes within each unit:

```cpp
base(k) = (c > 0) ? s / double(c) : datum::nan;
if (fin(y(k))) { s += y(k); ++c; }
```

These coincide on balanced panels with no missing outcomes. They differ when a unit has gaps, late entry, or missing outcomes. The Stata denominator is based on calendar time from the sample start, while the C++ denominator is based on observed finite past rows for that unit.

Impact:

- PMD LP-DiD estimates can differ from the Stata replication code outside balanced/no-missing settings.
- The difference may be desirable for robustness, but then it should be documented as an intentional deviation.

Suggested fix:

- Decide whether the package should match Stata exactly or use observed-past PMD.
- If exact translation is the goal, compute the denominator from normalized time as in the Stata code.
- If observed-past PMD is intended, document it and add tests showing the expected divergence.

### 5. Duplicate unit-time rows silently corrupt panel shifts

Files:

- `src/fLPPanel.cpp:79-83`
- `src/fLPDID.cpp:66-70`
- `R/fLPPanel.R:406-441`
- `R/fLPDID.R:292-300`

Both C++ panel engines build `time -> row` maps with `unordered_map::emplace`. If a unit has duplicate time rows, the first row wins and later duplicates are ignored. There is no validation in the wrappers.

Impact:

- Leads/lags, treatment changes, and event outcomes can use the wrong row with no warning.

Suggested fix:

- In `fLPPanel()` and `fLPDID()`, reject duplicate `(panel_id[1], panel_id[2])` pairs before calling C++.
- In the low-level C++ engines, validate map insertion success and fail cleanly for direct calls.

### 6. `fLPDID()` does not validate several arguments that can change semantics or cause out-of-bounds access

Files:

- `R/fLPDID.R:163-170`
- `R/fLPDID.R:193-200`
- `src/fLPDID.cpp:118-138`
- `src/fLPDID.cpp:172-181`

The wrapper validates `conf`, `panel_id`, and `treat`, but it does not fully validate:

- `pre >= 0` and `post >= 0`
- `ccc` in `{0, 1, 2}`
- `L` as a nonnegative integer when required
- scalar logical values for `nonabsorbing`, `pmd`, and `reweight`

In C++, `ccs0`, `ccsF`, and `ccsm` are only allocated when `nonabsorbing && ccc > 0`. A negative `ccc` with `nonabsorbing = TRUE` can fall into branches that index empty clean-control vectors.

Impact:

- Invalid user input can silently map to unintended clean-control behavior.
- Some invalid inputs can lead to undefined behavior in C++.

Suggested fix:

- Validate these arguments in `R/fLPDID.R` before the `.Call`.
- Add defensive C++ checks in `fLPDID_cpp()` because it is exported and callable directly.

### 7. `fLPPanel` HDFE residualization is not numerically equivalent to the reference implementation

Files:

- `src/fLPPanel.cpp:89-153`
- `inst/PanelLP/R_Code/panel_LP.R:64-98`

The source code uses `fixest::demean(..., iter = 10000, tol = 1e-8)` and then `pracma::mldivide()`. The C++ port implements alternating demeaning manually with `max_iter = 500`, an absolute Frobenius-norm stopping rule, and `arma::solve(..., solve_opts::fast)` without a rank fallback.

Risks:

- Multiple high-dimensional fixed effects may not converge within 500 iterations.
- The absolute tolerance is scale-dependent.
- Collinear or nearly collinear designs can throw or produce unstable coefficients.
- The reference code has explicit ill-conditioning warnings; the C++ port silently proceeds.

Suggested fix:

- Return or record a convergence flag per horizon.
- Use a relative convergence criterion.
- Add a rank-aware fallback to `solve()` or explicit column dropping, at least for controls beyond the target shock terms.
- Add tests comparing `fLPPanel_cpp()` to `panel_LP()` on designs with two fixed effects and near-collinear controls.

### 8. Low-level exported engines need stronger direct-call validation

Files:

- `src/fLP.cpp:278-299`
- `src/fLPIV.cpp:243-250`
- `src/fLPDID.cpp:33-49`

The high-level R wrappers do most validation, but the C++ functions are exported. Direct calls can still pass non-finite values, inconsistent lengths, duplicate panel keys, or invalid clean-control settings.

Specific gaps:

- `fLP_cpp()` and `fLPIV_cpp()` do not reject `Inf` values, while the R wrappers only drop `NA` via `complete.cases()`.
- `fLPDID_cpp()` does not check that `y`, `treat`, `i_index`, `t_index`, and `cl_index` have the same length, or that `X.n_rows == y.n_elem`.
- `fLPDID_cpp()` does not validate `pre_window`, `post_window`, `ccc`, `Lwin`, or thread count.

Suggested fix:

- Add low-level input checks at the exported C++ boundary.
- Prefer returning clear `Rcpp::stop()` messages before any OpenMP region.

### 9. Reweighted LP-DiD with covariates is not the same as the DDCG application code

Files:

- `src/fLPDID.cpp:208-220`
- `R/fLPDID.R:136-139`
- `inst/LPDID/codes/main_text_simulation/first_sim_rwLPDiD_estimation.do:24-55`
- `inst/LPDID/codes/ddcg_app/figure_4.do:152-165`

For the main simulation, the C++ reweighting is broadly consistent with inverse implicit weights up to a common scale. But the DDCG empirical application uses `teffects ra` with covariates for the reweighted specification. The wrapper documentation already says the covariate regression-adjustment version is not implemented.

Impact:

- This is not a bug if the package intentionally implements only the simpler reweighting path.
- It is a translation gap if the goal is to replicate every LP-DiD Stata application result.

Suggested fix:

- Keep the current warning/documentation explicit.
- Consider making `reweight = TRUE` with nonempty controls emit a message or warning that it is not the DDCG `teffects ra` estimator.

## Performance And Memory

### `fLP.cpp`

Hot spots:

- `src/fLP.cpp:108-118`: copies `Yh` and `Xh` for every horizon.
- `src/fLP.cpp:161`: allocates `Xreg` for every horizon.
- `src/fLP.cpp:179-195`: allocates `XU`, `S`, `G`, `Gamma_a`, and full `V` for every equation.

Recommendations:

- When `store_full = FALSE`, compute only the shock coefficient variance `V(sc, sc)` instead of the full covariance matrix for every equation.
- Avoid constructing `V = invXpX * G * invXpX` when only one diagonal element is needed. Use quadratic forms with row/column `sc`.
- Reuse per-horizon buffers where possible.
- For non-cumulative mode, consider Armadillo subviews rather than materializing `Yh` and `Xh`.
- Use `solve()` for beta where possible and reserve explicit inverse/pseudoinverse for the sandwich calculation.

### `fLPIV.cpp`

Hot spots:

- `src/fLPIV.cpp:115-129`: residualizes `Y`, `Z`, and `D` by allocating full residual matrices each horizon.
- `src/fLPIV.cpp:132-139`: recomputes first-stage projection at every horizon, even though in non-cumulative mode only the sample length changes.

Recommendations:

- Add a fast path for `nc == 0` that only demeans by an intercept or directly centers variables.
- For fixed samples or repeated horizon windows, cache reusable cross-products where sample truncation permits.
- Use a tolerance for weak or zero denominators instead of checking only exact zero at `src/fLPIV.cpp:168` and `src/fLPIV.cpp:187`.

### `fLPPanel.cpp`

Hot spots:

- `src/fLPPanel.cpp:140-148`: copies full residualized `y` and `X` every HDFE iteration.
- `src/fLPPanel.cpp:226-235`: precomputes a potentially large `n_obs x n_s x p_max` cube.
- `src/fLPPanel.cpp:360-380`: builds dense `T x T` matrices for the small-sample correction.

Recommendations:

- Replace full old-matrix copies in HDFE convergence checks with a running update norm.
- Make precomputing `sX_lags` optional or chunked for large panels.
- In the IK block, avoid materializing dense `T x T` matrices where only diagonal or quadratic forms are needed.
- Move fixed 90/95/99 confidence interval construction out of C++ if the R wrapper will rebuild bands anyway.

### `fLPDID.cpp`

Hot spots:

- `src/fLPDID.cpp:151-192`: every horizon scans all rows and builds temporary vectors.
- `src/fLPDID.cpp:197-245`: every horizon rebuilds compact year and cluster maps.
- `src/fLPDID.cpp:247-274`: every horizon builds and demeans a full design matrix.

Recommendations:

- Precompute `row_at(k, offset)` for all offsets used by the task list.
- Precompute raw `dtr`, `base`, and horizon-specific `dy` arrays when memory allows.
- Replace per-horizon unordered maps with sorted unique time/cluster codes built once, then remapped through integer arrays.
- Consider a separate no-control fast path because many LP-DiD specifications are treatment plus time fixed effects only.

### Threading And Console Output

Files:

- `src/fLP.cpp:87-99`
- `src/fLPIV.cpp:77-88`
- `src/fLPPanel.cpp:247-259`
- `src/fLPDID.cpp:355-359`

The engines print thread counts with `std::printf()` and several use `omp_set_num_threads()`, which changes OpenMP process state. For an R package, this creates noisy output and can surprise callers.

Recommendations:

- Add a `verbose` flag at the R wrapper level if user-facing messages are desired.
- Prefer `#pragma omp parallel for num_threads(actual_threads)` over global `omp_set_num_threads()`.
- Document possible BLAS/OpenMP oversubscription when R is linked to a multithreaded BLAS.

## Readability And API Consistency

- `src/fLP.h:61` has a default argument in the header. Project convention says defaults should live in `.cpp` definitions, not headers.
- Thread defaults are inconsistent: `fLP`, `fLPIV`, and `fLPPanel` use `0 = max_threads - 1`; `fLPDID` uses `n_threads <= 0 = max_threads` and the R wrapper defaults to `-1`.
- Panel index and gap-aware shift logic is duplicated in `fLPPanel.cpp` and `fLPDID.cpp`. A small shared internal helper would reduce translation drift.
- Several files mix algorithm comments, translation notes, and user-facing behavior in long blocks. The algorithms would be easier to audit if each engine separated: input validation, alignment/sample construction, estimation, variance, output assembly.

## Targeted Tests To Add

1. `fLPPanel()` formula lags on a panel with unit `A` observed at times `1` and `3`; assert `l(x, 1)` at time `3` is missing, not `x[time == 1]`.
2. `fLPDID()` with a `Date` time column and `d(y)` / `l(d(y), 1)` controls; assert the wrapper produces the same controls as normalized integer time.
3. Duplicate `(id, time)` rows for `fLPPanel()` and `fLPDID()`; assert a clear error.
4. `fLPIV()` first-stage F against an R reference that residualizes `D` and `Z` on controls and uses `df2 = T - rank([1, C]) - rank(Zr)`.
5. `fLPPanel_cpp()` versus `inst/PanelLP/R_Code/panel_LP.R` on a balanced toy panel for:
   - `small_sample = FALSE`
   - `small_sample = TRUE`
   - `cumulative = TRUE`
   - two fixed effects
6. `fLPDID()` absorbing treatment against the Stata formulas:
   - post: `D_j y = F_j y - L.y`, sample `D.treat == 1 | F_j.treat == 0`
   - pre: `Dm_j y = L_j y - L.y`, sample `D.treat == 1 | treat == 0`
7. `fLPDID()` nonabsorbing `ccc = 0, 1, 2` on a manually constructed switching panel with known clean-control sets.
8. Direct `.Call`/low-level calls with invalid dimensions and non-finite values; assert clean errors before computation.

## Overall Assessment

The LP engines have the right broad structure: horizon-level parallelism, formula wrappers that keep user syntax out of C++, and translation-aware sample construction. The highest-priority improvements are correctness guardrails rather than broad rewrites:

1. Fix panel-aware R formula shifts so wrapper-generated controls obey the same gap logic as the C++ engines.
2. Remove R API calls from OpenMP regions.
3. Correct LP-IV first-stage F degrees of freedom.
4. Decide and test the intended PMD baseline semantics for LP-DiD.
5. Add duplicate panel-key and direct-call validation.

After those, the largest performance gains should come from reducing per-horizon/per-equation allocation in `fLP.cpp`, moving confidence-band work out of `fLPPanel.cpp`, and replacing repeated unordered-map construction in `fLPDID.cpp` with precomputed integer indexes.
