#ifndef PANEL_SHIFT_H
#define PANEL_SHIFT_H

#include <RcppArmadillo.h>
#include <unordered_map>
#include <vector>

// -----------------------------------------------------------------------
// Shared panel index and gap-aware shift helpers for the panel LP engines
// (fLPPanel.cpp and fLPDID.cpp). Previously duplicated in both files.
//
// Design:
//   * `unit_of[k]`      : integer unit index of row k (0-based, dense)
//   * `unit_rows[u]`    : rows belonging to unit u, in original row order
//   * `t_to_row[u]`     : map from time value to row index, one map per unit
//
// All shifts look up `(unit_of[k], t_index[k] + L)`; if the target row does
// not exist for that unit, the shifted value is NaN. Duplicate (unit, time)
// pairs are rejected up-front by the R wrappers, so `unordered_map::emplace`
// is safe (first-wins is a valid invariant here).
// -----------------------------------------------------------------------

// Build unit_of, unit_rows, and t_to_row from (i_index, t_index). Called
// once per engine invocation and shared across horizons.
inline void panel_build_index(
    const arma::ivec& i_index,
    const arma::ivec& t_index,
    std::vector<arma::uword>&                         unit_of,
    std::vector<std::vector<arma::uword>>&            unit_rows,
    std::vector<std::unordered_map<int, arma::uword>>& t_to_row
) {
  const arma::uword n = i_index.n_elem;

  std::unordered_map<int, arma::uword> uix;
  uix.reserve(512);
  arma::uword next = 0;

  unit_of.assign(n, 0);
  for (arma::uword k = 0; k < n; ++k) {
    auto it = uix.find(i_index(k));
    if (it == uix.end()) {
      uix.emplace(i_index(k), next);
      unit_of[k] = next++;
    } else {
      unit_of[k] = it->second;
    }
  }

  unit_rows.assign(next, {});
  for (arma::uword k = 0; k < n; ++k) unit_rows[unit_of[k]].push_back(k);

  t_to_row.assign(next, {});
  for (arma::uword u = 0; u < next; ++u) {
    t_to_row[u].reserve(unit_rows[u].size() * 2);
    for (arma::uword k : unit_rows[u]) t_to_row[u].emplace(t_index(k), k);
  }
}

// Row lookup: for row k, return the row index of the same unit's
// observation at time t_index[k] + off, or -1 if it does not exist.
// Signed return type mirrors the fLPDID convention (uses `sword`).
inline arma::sword panel_row_at(
    arma::uword k,
    int off,
    const std::vector<arma::uword>&                         unit_of,
    const std::vector<std::unordered_map<int, arma::uword>>& t_to_row,
    const arma::ivec& t_index
) {
  const auto& m = t_to_row[unit_of[k]];
  auto it = m.find(t_index(k) + off);
  return (it == m.end()) ? arma::sword(-1) : arma::sword(it->second);
}

// Vector shift: for each row k, produce y1[k] = y0[row'] where row' is
// the same unit's observation at t_index[k] + L, else NaN. L > 0 = lead;
// L < 0 = lag (opposite to `panel_row_at`'s `off`, matching the
// convention used inside fLPPanel_internal).
inline arma::vec panel_time_shift(
    const arma::vec& y0,
    const std::vector<std::vector<arma::uword>>&            unit_rows,
    const std::vector<std::unordered_map<int, arma::uword>>& t_to_row,
    const arma::ivec& t_index,
    int L
) {
  const arma::uword n = y0.n_elem;
  arma::vec y1(n, arma::fill::value(arma::datum::nan));
  for (arma::uword u = 0; u < unit_rows.size(); ++u) {
    const auto& rows = unit_rows[u];
    const auto& map  = t_to_row[u];
    for (arma::uword k : rows) {
      auto it = map.find(t_index(k) + L);
      if (it != map.end()) y1(k) = y0(it->second);
    }
  }
  return y1;
}

#endif // PANEL_SHIFT_H
