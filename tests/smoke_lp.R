## Smoke + targeted tests for the local-projection engines after the
## lpreview_claude.md implementation and the follow-up performance pass.

suppressPackageStartupMessages({
  library(tidyMacro)
})

pass <- function(msg) cat(sprintf("  PASS %s\n", msg))
fail <- function(msg) stop(sprintf("FAIL %s\n", msg))

set.seed(2024)


# ---------------------------------------------------------------------
# 1. fLP: baseline smoke
# ---------------------------------------------------------------------
cat("[1] fLP baseline\n")
T <- 200L
df <- data.frame(
  y = arima.sim(list(ar = 0.5), n = T),
  x = rnorm(T),
  z = rnorm(T)
)
m <- fLP(y ~ x + z, data = df, horizons = 6L, shock = "x")
stopifnot(inherits(m, "fLP"), nrow(m$irfs) == 7L)
pass("fLP returns fLP object with H+1 IRF rows")


# ---------------------------------------------------------------------
# 2. fLPIV: F-stat df fix — numeric match to a hand-computed reference
# ---------------------------------------------------------------------
cat("[2] fLPIV first-stage F df\n")
n <- 250L
z_iv <- rnorm(n)
d_end <- 0.4 * z_iv + rnorm(n, sd = 0.5)
c1 <- rnorm(n)
c2 <- rnorm(n)
c3 <- rnorm(n)
c4 <- rnorm(n)
c5 <- rnorm(n)
y_iv <- 0.2 * d_end + 0.1 * c1 + rnorm(n, sd = 0.5)
iv_df <- data.frame(y_iv, d_end, z_iv, c1, c2, c3, c4, c5)
m_iv <- fLPIV(
  formula = y_iv ~ d_end + c1 + c2 + c3 + c4 + c5,
  instruments = ~z_iv,
  data = iv_df,
  endog = "d_end",
  horizons = 4L,
  conf = 0.90,
  nw_lags_iv = 2L,
  cumulative = FALSE
)
stopifnot(
  is.numeric(m_iv$Fstat_fs),
  length(m_iv$Fstat_fs) == 5L,
  all(is.finite(m_iv$Fstat_fs))
)
Th_ref <- n
kc_ref <- ncol(iv_df) - 3L
nz_ref <- 1L
lm1 <- lm(d_end ~ c1 + c2 + c3 + c4 + c5 + z_iv, data = iv_df)
lm0 <- lm(d_end ~ c1 + c2 + c3 + c4 + c5, data = iv_df)
rss1 <- sum(residuals(lm1)^2)
rss0 <- sum(residuals(lm0)^2)
F_ref <- ((rss0 - rss1) / nz_ref) / (rss1 / (Th_ref - kc_ref - 1L - nz_ref))
if (abs(m_iv$Fstat_fs[1] - F_ref) / max(1, F_ref) > 1e-6) {
  fail(sprintf(
    "F-stat mismatch: engine=%.6f, reference=%.6f",
    m_iv$Fstat_fs[1],
    F_ref
  ))
}
pass(sprintf(
  "fLPIV F(h=0) matches Th - kc - nz reference: %.4f",
  m_iv$Fstat_fs[1]
))


# ---------------------------------------------------------------------
# 3. fLPPanel: gap-aware lag runs to completion on a gap panel
# ---------------------------------------------------------------------
cat("[3] fLPPanel gap-aware lag (smoke)\n")
gap_df <- data.frame(
  id = rep(1:4, each = 20),
  time = rep(1:20, times = 4),
  y = rnorm(80),
  x = rnorm(80)
)
gap_df <- gap_df[!(gap_df$id == 1 & gap_df$time == 5), ]
res_gap <- tryCatch(
  fLPPanel(
    y ~ x + l(x, 1) | id + time,
    data = gap_df,
    panel_id = c("id", "time"),
    shock = "x",
    horizons = 3L
  ),
  error = function(e) e
)
if (inherits(res_gap, "error")) {
  fail(sprintf("fLPPanel gap ran: %s", conditionMessage(res_gap)))
}
pass("fLPPanel completed on gap-containing panel")


# ---------------------------------------------------------------------
# 4. fLPPanel: duplicate (id, time) rejection
# ---------------------------------------------------------------------
cat("[4] fLPPanel duplicate (id, time) rejection\n")
dup_df <- rbind(gap_df, gap_df[1, , drop = FALSE])
err <- tryCatch(
  fLPPanel(
    y ~ x | id + time,
    data = dup_df,
    panel_id = c("id", "time"),
    shock = "x",
    horizons = 2L
  ),
  error = function(e) conditionMessage(e)
)
if (!grepl("duplicate", err)) {
  fail(paste("expected duplicate error, got:", err))
}
pass("fLPPanel rejects duplicate (id, time) rows")


# ---------------------------------------------------------------------
# 5. fLPDID: absorbing DiD smoke
# ---------------------------------------------------------------------
cat("[5] fLPDID absorbing DiD\n")
n_i <- 30L
n_t <- 30L
dgp <- expand.grid(id = 1:n_i, year = 1:n_t)
dgp$treat <- as.integer(dgp$id > n_i / 2 & dgp$year >= 15)
dgp$y <- 0.5 * dgp$treat + rnorm(nrow(dgp))
m_did <- fLPDID(
  y ~ treat,
  data = dgp,
  panel_id = c("id", "year"),
  treat = "treat",
  post = 5L,
  pre = 3L
)
stopifnot(
  inherits(m_did, "fLPDID"),
  all(c("event_time", "estimate", "se") %in% names(m_did))
)
pass("fLPDID returns fLPDID tibble with expected columns")


# ---------------------------------------------------------------------
# 6. fLPDID: PMD baseline runs on gap panel (Stata calendar-time formula)
# ---------------------------------------------------------------------
cat("[6] fLPDID PMD baseline on gap panel\n")
pmd_df <- expand.grid(id = 1L:3L, year = 1L:8L)
pmd_df$y <- as.numeric(seq_len(nrow(pmd_df)))
pmd_df$treat <- as.integer(pmd_df$id == 3L & pmd_df$year >= 6L)
pmd_df <- pmd_df[!(pmd_df$id == 2L & pmd_df$year == 4L), ]
m_pmd <- tryCatch(
  fLPDID(
    y ~ treat,
    data = pmd_df,
    panel_id = c("id", "year"),
    treat = "treat",
    post = 2L,
    pre = 1L,
    pmd = TRUE
  ),
  error = function(e) e
)
if (inherits(m_pmd, "error")) {
  fail(paste("PMD baseline engine call failed:", conditionMessage(m_pmd)))
}
pass("fLPDID pmd = TRUE runs on gap panel (Stata calendar-time baseline)")


# ---------------------------------------------------------------------
# 7. fLPDID argument validation
# ---------------------------------------------------------------------
cat("[7] fLPDID argument validation\n")
err1 <- tryCatch(
  fLPDID(
    y ~ treat,
    data = dgp,
    panel_id = c("id", "year"),
    treat = "treat",
    post = 3L,
    pre = 3L,
    ccc = 5L
  ),
  error = function(e) conditionMessage(e)
)
if (!grepl("ccc", err1)) {
  fail(paste("expected ccc error, got:", err1))
}
err2 <- tryCatch(
  fLPDID(
    y ~ treat,
    data = dgp,
    panel_id = c("id", "year"),
    treat = "treat",
    post = -1L,
    pre = 3L
  ),
  error = function(e) conditionMessage(e)
)
if (!grepl("post", err2)) {
  fail(paste("expected post error, got:", err2))
}
pass("fLPDID rejects invalid ccc and post values")


# ---------------------------------------------------------------------
# 8. fLP_cpp: Inf inputs rejected at the C++ boundary
# ---------------------------------------------------------------------
cat("[8] fLP_cpp Inf rejection\n")
bad_Y <- matrix(c(1, Inf, 2, 3), ncol = 1L)
bad_X <- matrix(rnorm(4L), ncol = 1L)
err_inf <- tryCatch(
  tidyMacro:::fLP_cpp(
    Y = bad_Y,
    X = bad_X,
    H = 1L,
    shock_col = 0L,
    conf_level = 0.9,
    nw_lags_base = 1L
  ),
  error = function(e) conditionMessage(e)
)
if (!grepl("finite", err_inf)) {
  fail(paste("expected finite error, got:", err_inf))
}
pass("fLP_cpp rejects Inf inputs at the C++ boundary")


# =====================================================================
# TARGETED TESTS from lpreview_claude.md §"Targeted Tests To Add"
# =====================================================================

# ---------------------------------------------------------------------
# T1. fLPPanel gap-aware l() value: unit observed at times 1 and 3 —
#     l(x, 1) at t=3 must be NA (not x at t=1). We verify this via the
#     C++ engine directly to inspect the constructed control column.
# ---------------------------------------------------------------------
cat("[T1] fLPPanel l(x, 1) at t=3 is NA when t=2 missing\n")
mini <- data.frame(
  id = rep(1L, 2L),
  time = c(1L, 3L),
  x = c(10, 20)
)
# Replicate the wrapper's helper: normalize + gap-aware shift.
tn <- tidyMacro:::.fLPPanel_normalize_time(mini$time)
# gap panel: normalized times should be 1 and 3 (step = min(diff) = 2 -> 1, 2?)
# with step = 2 the normalization gives (1 - 1)/2 + 1 = 1 and (3 - 1)/2 + 1 = 2.
# That collapses the gap. Use a bigger panel so step = 1 and the gap survives.
mini <- data.frame(
  id = rep(c(1L, 2L), each = 3L),
  time = c(1L, 3L, 4L, 1L, 2L, 3L),
  x = c(10, 20, 30, 40, 50, 60)
)
tn <- tidyMacro:::.fLPPanel_normalize_time(mini$time)
stopifnot(all(tn$on_grid))
# Gap-aware helper directly (bypass formula parsing):
res_ex <- tidyMacro:::.fLPPanel_expand_lag_terms(
  formula_chr = "y ~ l(x, 1)",
  data = mini,
  id_vec = mini$id,
  t_norm = tn$t_norm
)
x_l1 <- res_ex$data$x_l1
# Row for id=1, time=3 (index 2 in mini): l(x,1) requires t=2 for id=1, missing -> NA.
# Row for id=1, time=4 (index 3): requires t=3 for id=1, present -> 20.
if (!is.na(x_l1[2])) {
  fail(sprintf("expected NA at (id=1, t=3), got %s", x_l1[2]))
}
if (!isTRUE(all.equal(x_l1[3], 20))) {
  fail(sprintf("expected 20 at (id=1, t=4), got %s", x_l1[3]))
}
if (!is.na(x_l1[4])) {
  fail(sprintf("expected NA at (id=2, t=1), got %s", x_l1[4]))
}
if (!isTRUE(all.equal(x_l1[5], 40))) {
  fail(sprintf("expected 40 at (id=2, t=2), got %s", x_l1[5]))
}
pass("fLPPanel l(x, 1) preserves calendar-time gap for id=1 at t=3")


# ---------------------------------------------------------------------
# T2. fLPDID with a Date time column and d(y): the d() column produced
#     from the Date-typed wrapper must equal the one from a numeric
#     integer-time version of the same data.
# ---------------------------------------------------------------------
cat("[T2] fLPDID d() with Date time == integer-time equivalent\n")
n_i <- 8L
n_t <- 20L
grid <- expand.grid(id = 1L:n_i, k = 0L:(n_t - 1L))
grid$date_time <- as.Date("2000-01-01") + grid$k * 30L # monthly-ish
grid$int_time <- grid$k + 1L
grid$y <- rnorm(nrow(grid))
grid$treat <- as.integer(grid$id > n_i / 2 & grid$k >= 10L)
# Drop a row per unit to create gaps, so d() is a real test.
drop_ix <- which(grid$id == 1L & grid$k == 5L)
grid <- grid[-drop_ix, ]

# Wrapper does not expose the built column, but the shift helper does.
tn_date <- tidyMacro:::.fLPPanel_normalize_time(grid$date_time)
tn_int <- tidyMacro:::.fLPPanel_normalize_time(grid$int_time)
# Normalized time indexes should match modulo the step
if (!identical(tn_date$t_norm, tn_int$t_norm)) {
  fail("normalized time differs between Date and integer parameterizations")
}

d_date <- grid$y -
  tidyMacro:::.flpdid_shift(grid$y, grid$id, tn_date$t_norm, 1L)
d_int <- grid$y - tidyMacro:::.flpdid_shift(grid$y, grid$id, tn_int$t_norm, 1L)
if (
  !identical(is.na(d_date), is.na(d_int)) ||
    !isTRUE(all.equal(d_date[!is.na(d_date)], d_int[!is.na(d_int)]))
) {
  fail("d(y) differs between Date-typed and integer-typed time columns")
}
pass("fLPDID d(y) with Date time equals integer-time equivalent")


# ---------------------------------------------------------------------
# T3. fLPPanel_cpp on a balanced toy panel — verify the point estimate
#     against a hand-coded within-year OLS on s*x (matches the panel_LP
#     asymptotic LAHR regression at horizon 0 with no FE and no lags).
#     This is a numeric sanity check against R's lm(), not a rerun of
#     the reference R port (which requires fixest + pracma to install).
# ---------------------------------------------------------------------
cat("[T3] fLPPanel_cpp horizon 0 estimate matches lm() reference\n")
set.seed(11)
N <- 60L
T <- 20L
pdf <- expand.grid(unit = 1:N, tt = 1:T)
pdf$shock <- rep(rnorm(T), each = N)
pdf$size <- rep(rnorm(N), times = T)
pdf$y <- pdf$size * pdf$shock + rnorm(nrow(pdf))
m_panel <- fLPPanel(
  y ~ shock:size,
  data = pdf,
  panel_id = c("unit", "tt"),
  shock = "shock:size",
  horizons = 0L,
  small_sample = FALSE,
  cumulative = FALSE,
  p_max = 0L,
  n_threads = 1L
)
# No `|` FE and no intercept term ⇒ engine solves `y = b * (shock*size)`.
ref <- unname(coef(lm(y ~ I(shock * size) - 1, data = pdf))[1])
if (abs(m_panel$irfs[1, 1] - ref) / max(1, abs(ref)) > 1e-8) {
  fail(sprintf(
    "fLPPanel horizon-0 estimate mismatch: engine=%.6f, lm=%.6f",
    m_panel$irfs[1, 1],
    ref
  ))
}
pass(sprintf("fLPPanel_cpp horizon 0 matches lm(): %.6f", m_panel$irfs[1, 1]))


# ---------------------------------------------------------------------
# T4. fLPDID absorbing — check against the direct Stata formula:
#     post h = 1 estimate should equal E[F1.y - L.y | D.treat==1]
#                                    - E[F1.y - L.y | F1.treat==0, D.treat==0]
#     on a tiny hand-constructed panel with no controls and 2-cohort data.
# ---------------------------------------------------------------------
cat("[T4] fLPDID absorbing post h=1 matches direct DiD contrast\n")
set.seed(2)
tiny <- expand.grid(id = 1L:20L, year = 1L:12L)
tiny$treat <- as.integer(tiny$id >= 11L & tiny$year >= 6L) # cohort 6, half treated
tiny$y <- 0.3 * tiny$treat + 0.05 * tiny$year + rnorm(nrow(tiny), sd = 0.5)

m <- fLPDID(
  y ~ treat,
  data = tiny,
  panel_id = c("id", "year"),
  treat = "treat",
  post = 3L,
  pre = 3L
)
b_engine <- m$estimate[m$event_time == 1L]

tiny <- tiny[order(tiny$id, tiny$year), ]
tiny$y_L1 <- ave(tiny$y, tiny$id, FUN = function(v) c(NA, head(v, -1)))
tiny$y_F1 <- ave(tiny$y, tiny$id, FUN = function(v) c(tail(v, -1), NA))
tiny$t_F1 <- ave(tiny$treat, tiny$id, FUN = function(v) c(tail(v, -1), NA))
tiny$D <- tiny$treat -
  tiny$y_L1 * 0 -
  ave(tiny$treat, tiny$id, FUN = function(v) c(NA, head(v, -1)))
# clean-control sample: D.treat == 1  |  F1.treat == 0
sample <- tiny$D == 1 | tiny$t_F1 == 0
d <- tiny[which(sample), ]
d$dy <- d$y_F1 - d$y_L1
# Add year fixed effects and regress dy on D.treat.
lm_ref <- lm(dy ~ factor(year) + D, data = d, na.action = na.exclude)
b_ref <- unname(coef(lm_ref)["D"])
if (abs(b_engine - b_ref) / max(1, abs(b_ref)) > 1e-6) {
  fail(sprintf(
    "fLPDID absorbing post h=1 mismatch: engine=%.6f, ref=%.6f",
    b_engine,
    b_ref
  ))
}
pass(sprintf(
  "fLPDID absorbing post h=1 matches direct DiD contrast: %.6f",
  b_engine
))


# ---------------------------------------------------------------------
# T5. fLPDID nonabsorbing ccc = 0, 1, 2 — all three flavors should run
#     on a hand-constructed switching panel and produce numeric estimates.
#     ccc = 2 should filter more rows than ccc = 1 which filters more
#     than ccc = 0, so the reported nobs weakly decreases.
# ---------------------------------------------------------------------
cat("[T5] fLPDID nonabsorbing ccc = 0,1,2 on switching panel\n")
set.seed(3)
sw <- expand.grid(id = 1L:15L, year = 1L:25L)
sw <- sw[order(sw$id, sw$year), ]
# Randomly switch treatment on/off with 0.15 prob per year, per unit.
sw$treat <- 0L
for (uid in unique(sw$id)) {
  ix <- which(sw$id == uid)
  toggles <- rbinom(length(ix), 1L, 0.15)
  state <- 0L
  for (j in seq_along(ix)) {
    if (toggles[j] == 1L) {
      state <- 1L - state
    }
    sw$treat[ix[j]] <- state
  }
}
sw$y <- 0.4 * sw$treat + 0.02 * sw$year + rnorm(nrow(sw))
nobs_by_ccc <- integer(3)
for (cc in 0L:2L) {
  m_cc <- fLPDID(
    y ~ treat,
    data = sw,
    panel_id = c("id", "year"),
    treat = "treat",
    post = 3L,
    pre = 3L,
    nonabsorbing = TRUE,
    L = 5L,
    ccc = cc
  )
  if (!all(is.finite(m_cc$estimate[m_cc$event_time >= 0]))) {
    fail(sprintf("ccc=%d produced non-finite post estimates", cc))
  }
  # Use the observation count at the first post horizon as the proxy.
  nobs_by_ccc[cc + 1L] <- m_cc$nobs[m_cc$event_time == 0L]
}
if (!(nobs_by_ccc[1] >= nobs_by_ccc[2] && nobs_by_ccc[2] >= nobs_by_ccc[3])) {
  fail(sprintf(
    "expected monotone nobs across ccc: got (%d, %d, %d)",
    nobs_by_ccc[1],
    nobs_by_ccc[2],
    nobs_by_ccc[3]
  ))
}
pass(sprintf(
  "fLPDID nonabsorbing ccc={0,1,2} monotone in nobs: (%d >= %d >= %d)",
  nobs_by_ccc[1],
  nobs_by_ccc[2],
  nobs_by_ccc[3]
))


# ---------------------------------------------------------------------
# T6. Direct C++ boundary: length-mismatch inputs must fail cleanly.
# ---------------------------------------------------------------------
cat("[T6] fLPDID_cpp rejects length-mismatched inputs\n")
n_bad <- 20L
y_bad <- rnorm(n_bad)
tr_bad <- rbinom(n_bad, 1L, 0.5)
X_bad <- matrix(0, n_bad, 0L)
i_bad <- rep(1L, n_bad)
t_bad <- seq_len(n_bad)
cl_bad <- rep(1L, n_bad - 1L) # length mismatch
err_len <- tryCatch(
  tidyMacro:::fLPDID_cpp(
    y = y_bad,
    treat = tr_bad,
    X = X_bad,
    i_index = i_bad,
    t_index = t_bad,
    cl_index = cl_bad,
    pre_window = 2L,
    post_window = 2L,
    nonabsorbing = FALSE,
    Lwin = 0L,
    ccc = 0L,
    pmd = FALSE,
    reweight = FALSE,
    n_threads = 1L
  ),
  error = function(e) conditionMessage(e)
)
if (!grepl("length|rows", err_len, ignore.case = TRUE)) {
  fail(paste("expected length-mismatch error, got:", err_len))
}
pass("fLPDID_cpp rejects length-mismatched inputs")
