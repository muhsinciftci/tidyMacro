# Local Projections C++ Review — Claude Verification

Scope: I re-verified the 9 findings in [lpreview.md](lpreview.md) against the actual source under [src/fLP.cpp](src/fLP.cpp), [src/fLPIV.cpp](src/fLPIV.cpp), [src/fLPPanel.cpp](src/fLPPanel.cpp), [src/fLPDID.cpp](src/fLPDID.cpp) and their R wrappers, plus the translation sources [inst/PanelLP/R_Code/panel_LP.R](inst/PanelLP/R_Code/panel_LP.R) and the LP-DiD Stata scripts under [inst/LPDID/codes/](inst/LPDID/codes/). VAR files and generated `docs/`/`man/` were skipped per project rules.

Verdicts below flag which items are real, which need nuance, and which I would deprioritize when we implement.

---

## Findings — Verdicts

### 1. Panel formula lags/leads do not preserve time gaps — **CORRECT (P0)**

Confirmed. [R/fLPPanel.R:66-77](R/fLPPanel.R#L66-L77) shifts by row position within unit and only checks `unit_id[idx] == unit_id[idx - k]`. The C++ engine ([src/fLPPanel.cpp:32-50](src/fLPPanel.cpp#L32-L50) `panel_time_shift`) and the reference ([inst/PanelLP/R_Code/panel_LP.R:15-36](inst/PanelLP/R_Code/panel_LP.R#L15-L36)) both look up exact `(unit, t+L)`. The wrapper and engine therefore disagree on unbalanced panels.

[R/fLPDID.R:40-60](R/fLPDID.R#L40-L60) `.flpdid_shift()` does check raw time increments equal to `k`, but the wrapper never normalizes the time index first, so `Date`/`POSIXct`/step-≠-1 numeric columns break `d(x)` and `l(d(x), …)`.

Implementation plan: normalize time once at the top of both wrappers, then reuse a single gap-aware shift helper. Add a duplicate `(id, time)` guard at the same entry point (see Finding 5).

### 2. `fLPPanel` calls R math and throws `Rcpp::stop` inside an OpenMP loop — **CORRECT (P0)**

Confirmed. The horizon loop at [src/fLPPanel.cpp:263-265](src/fLPPanel.cpp#L263-L265) is `#pragma omp parallel for`, and inside it we call `Rcpp::stop()` at line 294, plus `R::qt()` / `R::pt()` at lines 404-416. This violates the project rule against R API calls in parallel regions.

Implementation plan: replace `Rcpp::stop()` with per-horizon status flags in a plain struct; drop `R::qt`/`R::pt` inside the loop and either (a) do band construction serially after the join, or (b) let [R/fLPPanel.R:472-500](R/fLPPanel.R#L472-L500) rebuild bands (it already does for user-requested confidence levels).

### 3. LP-IV first-stage F uses wrong denominator df after FWL — **CORRECT (P0)**

Confirmed at [src/fLPIV.cpp:141-148](src/fLPIV.cpp#L141-L148). The denominator is `rss_unr / (Th - nz)`. But `Dr` and `Zr` have already been partialled out on `Cext = [1, C]` (see [src/fLPIV.cpp:126-129](src/fLPIV.cpp#L126-L129)), so under FWL the correct residual df for the full regression `D ~ [1, C, Z]` is `Th - kc - nz`, not `Th - nz`. This *biases the reported first-stage F upward* whenever there are controls, and the bias grows with `nc`. Point estimates are unaffected.

Also [src/fLPIV.cpp:56](src/fLPIV.cpp#L56) input-check uses `kc + 1` while [R/fLPIV.R:199-203](R/fLPIV.R#L199-L203) uses `ncol(Z) + kc`; the C++ boundary is looser than the R wrapper.

Implementation plan: use `df2 = Th - kc - nz` (or ideally `Th - rank([1,C]) - rank(Zr)`); return `NaN` for `Fstat_fs` when df2 ≤ 0; tighten the C++ boundary check to match the wrapper.

### 4. LP-DiD PMD baseline differs from the Stata source — **CORRECT (P1, needs a decision)**

Confirmed. Stata ([inst/LPDID/codes/main_text_simulation/first_sim_pmd_LPDiD_estimation.do:14-25](inst/LPDID/codes/main_text_simulation/first_sim_pmd_LPDiD_estimation.do)) divides `L.cumulative_y` by *calendar* `time - 1`. C++ ([src/fLPDID.cpp:91-101](src/fLPDID.cpp#L91-L101)) divides by the *count of finite observed* past outcomes. On balanced, complete panels they coincide; on gaps / late entry / NA outcomes they diverge.

This is not obviously a bug — the observed-past denominator is arguably more defensible on unbalanced panels — but it is a translation deviation. We should decide policy before implementing.

Implementation plan: pick one (I recommend keeping observed-past as the default and adding a `pmd_denominator = c("observed", "calendar")` toggle) and add a test that pins the chosen semantics against the Stata numbers on a balanced toy panel.

### 5. Duplicate `(id, time)` rows silently corrupt panel shifts — **CORRECT (P0, cheap fix)**

Confirmed. Both maps use `emplace`, which is no-op on collision: [src/fLPPanel.cpp:82](src/fLPPanel.cpp#L82) and [src/fLPDID.cpp:69](src/fLPDID.cpp#L69). Neither wrapper rejects duplicates.

Implementation plan: add a duplicate-key check in [R/fLPPanel.R](R/fLPPanel.R) and [R/fLPDID.R](R/fLPDID.R) before `.Call`, and have the C++ boundary at least assert or return an error status for direct callers.

### 6. `fLPDID()` does not validate several arguments — **CORRECT (P1)**

Confirmed. [R/fLPDID.R:163-200](R/fLPDID.R#L163-L200) checks `conf`, `panel_id`, `treat`, but not `pre >= 0`, `post >= 0`, `ccc ∈ {0,1,2}`, `L` nonnegative integer when required, or scalar-logical shape of `nonabsorbing` / `pmd` / `reweight`. In C++, `ccs0`/`ccsF`/`ccsm` are only allocated when `nonabsorbing && ccc > 0` ([src/fLPDID.cpp:118-138](src/fLPDID.cpp#L118-L138)), so a negative `ccc` with `nonabsorbing = TRUE` can reach branches that index empty vectors.

Implementation plan: add wrapper-level checks, plus a small validation block at the top of `fLPDID_cpp` since it is exported.

### 7. `fLPPanel` HDFE is not numerically equivalent to `panel_LP.R` — **CORRECT (P2)**

Confirmed. [src/fLPPanel.cpp:89-153](src/fLPPanel.cpp#L89-L153) uses alternating demeaning with `max_iter = 500`, absolute Frobenius tolerance `1e-8`, and `arma::solve(..., solve_opts::fast)` — no rank fallback. Reference ([inst/PanelLP/R_Code/panel_LP.R:64-98](inst/PanelLP/R_Code/panel_LP.R#L64-L98)) uses `fixest::demean(iter = 10000, tol = 1e-8)` and `pracma::mldivide()` which is rank-aware.

For one FE this is fine; for two high-cardinality FEs, 500 iterations of alternating demeaning may not converge and the absolute norm is scale-sensitive. Also silent failures on rank-deficient controls.

Implementation plan: (a) bump `max_iter`, switch to relative norm; (b) fall back to `arma::solve(..., solve_opts::equilibrate)` or `pinv` on failure and record it; (c) surface a per-horizon convergence flag. Not urgent unless someone reports numerical drift.

### 8. Low-level exported engines need direct-call validation — **CORRECT (P1)**

Confirmed. The exported `fLP_cpp` / `fLPIV_cpp` / `fLPDID_cpp` accept whatever a direct `.Call` gives them: no `Inf` guard ([src/fLP.cpp:278-299](src/fLP.cpp#L278-L299), [src/fLPIV.cpp:243-250](src/fLPIV.cpp#L243-L250)), no length-consistency check for `y/treat/i_index/t_index/cl_index`, no `pre_window/post_window/ccc/Lwin/thread_count` validation ([src/fLPDID.cpp:33-49](src/fLPDID.cpp#L33-L49)).

Implementation plan: one small validation block per engine at the C++ boundary, before any OpenMP region. This is defensive; the R wrappers already cover most legitimate users.

### 9. Reweighted LP-DiD with covariates ≠ DDCG `teffects ra` — **CORRECT but not a bug (P3)**

Confirmed. [src/fLPDID.cpp:208-220](src/fLPDID.cpp#L208-L220) implements simple inverse-implicit-weight reweighting; the DDCG application ([inst/LPDID/codes/ddcg_app/figure_4.do:152-165](inst/LPDID/codes/ddcg_app/figure_4.do)) uses `teffects ra` (regression adjustment). The wrapper docstring already flags that the RA version isn't implemented.

Implementation plan: strengthen the message — emit a `warning()` (not just doc) whenever `reweight = TRUE` and `!is.null(controls)`, telling the user this is not `teffects ra`. Leave the estimator as is.

---

## Performance And Memory — Verdicts

All hot-spot lines in the review reproduce. Concrete confirmations:

- [src/fLP.cpp:108-118](src/fLP.cpp#L108-L118), [161](src/fLP.cpp#L161), [179-195](src/fLP.cpp#L179-L195): per-horizon `Yh`/`Xh` copies, per-horizon `Xreg`, per-equation full `V`. When `store_full = FALSE` we only need `V(sc, sc)`.
- [src/fLPIV.cpp:115-129](src/fLPIV.cpp#L115-L129), [132-139](src/fLPIV.cpp#L132-L139): FWL residuals recomputed each horizon; add an `nc == 0` fast path (just centering).
- [src/fLPPanel.cpp:140-148](src/fLPPanel.cpp#L140-L148), [226-235](src/fLPPanel.cpp#L226-L235), [360-380](src/fLPPanel.cpp#L360-L380): full old-matrix copy for HDFE convergence check, large `n_obs × n_s × p_max` cube, dense `T × T` matrices in the small-sample correction block.
- [src/fLPDID.cpp:151-192](src/fLPDID.cpp#L151-L192), [197-245](src/fLPDID.cpp#L197-L245), [247-274](src/fLPDID.cpp#L247-L274): every horizon re-scans rows and rebuilds year/cluster maps and the demeaned design matrix.

Priority ordering when we implement:
1. `fLP.cpp` — drop full-`V` allocation when `store_full = FALSE`; large win, small diff.
2. `fLPDID.cpp` — precompute row-offset tables and integer-remapped year/cluster codes once; also add a no-control fast path (this is the majority of LP-DiD calls).
3. `fLPPanel.cpp` — running-norm HDFE and optional/chunked `sX_lags` precompute.
4. `fLPIV.cpp` — `nc == 0` fast path.

## Threading And Console Output — **CORRECT**

Confirmed at [src/fLP.cpp:87-99](src/fLP.cpp#L87-L99), [src/fLPIV.cpp:77-88](src/fLPIV.cpp#L77-L88), [src/fLPPanel.cpp:247-259](src/fLPPanel.cpp#L247-L259), [src/fLPDID.cpp:355-359](src/fLPDID.cpp#L355-L359). Uses `std::printf` for thread messages and `omp_set_num_threads()` for global state.

Implementation plan: gate prints behind a `verbose` flag, prefer `num_threads(t)` clause over `omp_set_num_threads()`.

## Readability And API Consistency — **CORRECT**

- Header default at [src/fLP.h:61](src/fLP.h#L61) violates the project convention (defaults in `.cpp`).
- Thread convention drift: `fLP`/`fLPIV`/`fLPPanel` treat `n_threads == 0` as `max_threads - 1`; `fLPDID` treats `n_threads <= 0` as `max_threads`; the `fLPDID` R wrapper defaults to `-1` while the others default to `0`. Pick one and unify.
- Panel index + gap-shift logic is duplicated in `fLPPanel.cpp` and `fLPDID.cpp` — a shared internal header would prevent drift.

---

## Overall Assessment

The upstream review is accurate. Every numbered finding reproduces against the current tree, and the file/line pointers are correct.

**Suggested implementation order** (bundle logically, not one PR each):

1. **Correctness bundle (P0):** Findings 1, 2, 3, 5. These are user-visible bugs and the fixes are localized. Ship with tests from the review's "Targeted Tests To Add" list (items 1, 2, 3, 4, 6).
2. **Policy decision + fix (P1):** Finding 4 — pick observed-past vs calendar-time PMD, document it, add a test that pins the choice.
3. **Defensive validation (P1):** Findings 6 and 8 together — one pass across `R/fLPDID.R` and all three exported `_cpp` engines. Also add duplicate-key checks from Finding 5.
4. **Numerical robustness (P2):** Finding 7 HDFE hardening — convergence flag, relative tolerance, rank-aware solve fallback.
5. **Documentation nudge (P3):** Finding 9 — upgrade the reweight-with-controls note from doc to runtime warning.
6. **Performance pass:** in the priority order above, ideally after (1)-(3) so we aren't optimizing code we're about to rewrite.
7. **Cleanup:** header default, thread-count unification, shared panel-shift helper, gate `printf` behind `verbose`.
