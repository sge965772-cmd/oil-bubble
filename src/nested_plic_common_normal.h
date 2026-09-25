#ifndef NESTED_PLIC_COMMON_NORMAL_H
#define NESTED_PLIC_COMMON_NORMAL_H

#include <float.h>
#include <math.h>
#include <stdbool.h>

typedef enum {
  NESTED_PLIC_ALIGNED_BLEND = 0,
  NESTED_PLIC_NONFINITE_INPUT = 1,
  NESTED_PLIC_UNRESOLVED_NORMALS = 2,
  NESTED_PLIC_OPPOSING_NORMALS = 3,
  NESTED_PLIC_DEGENERATE_BLEND = 4,
  NESTED_PLIC_SLIVER_FALLBACK = 5,
  NESTED_PLIC_NEIGHBOR_FALLBACK = 6,
  NESTED_PLIC_PAIRED_NEIGHBOR_CONSENSUS = 7,
  NESTED_PLIC_NEGLIGIBLE_ENDPOINT_FALLBACK = 8,
  NESTED_PLIC_CENTER_PAIRED_NEIGHBOR_CONSENSUS = 9,
  NESTED_PLIC_ENDPOINT_GRADIENT_FALLBACK = 10,
  NESTED_PLIC_ENDPOINT_CENTER_BLEND = 11,
  NESTED_PLIC_COINCIDENT_INVENTORY_BLEND = 12
} NestedPlicReason;

#define NESTED_PLIC_NEIGHBOR_COUNT 4

typedef struct {
  double inner_fraction;
  double outer_fraction;
  double inner_x;
  double inner_y;
  double outer_x;
  double outer_y;
} NestedPlicNeighbor;

typedef struct {
  double x;
  double y;
  double inner_x;
  double inner_y;
  double outer_x;
  double outer_y;
  double inner_norm_l1;
  double outer_norm_l1;
  double inner_raw_weight;
  double outer_raw_weight;
  double inner_weight;
  double outer_weight;
  double normal_dot;
  int valid;
  int used_sliver_fallback;
  int fallback_neighbor;
  int fallback_uses_outer;
  int fallback_support_count;
  double fallback_quality;
  NestedPlicReason reason;
} NestedPlicCommonNormal;

#ifndef NESTED_PLIC_SLIVER_QUALITY_MAX
# define NESTED_PLIC_SLIVER_QUALITY_MAX 1.e-8
#endif
#ifndef NESTED_PLIC_SLIVER_DOMINANCE_RATIO
# define NESTED_PLIC_SLIVER_DOMINANCE_RATIO 16.
#endif
#ifndef NESTED_PLIC_NEIGHBOR_QUALITY_MIN
# define NESTED_PLIC_NEIGHBOR_QUALITY_MIN 1.e-4
#endif
#ifndef NESTED_PLIC_CENTER_FRACTION_TOLERANCE
# define NESTED_PLIC_CENTER_FRACTION_TOLERANCE 1.e-6
#endif
#ifndef NESTED_PLIC_NEGLIGIBLE_ENDPOINT_TOLERANCE
# define NESTED_PLIC_NEGLIGIBLE_ENDPOINT_TOLERANCE 1.e-12
#endif
#ifndef NESTED_PLIC_NEGLIGIBLE_ENDPOINT_DOMINANCE_RATIO
# define NESTED_PLIC_NEGLIGIBLE_ENDPOINT_DOMINANCE_RATIO 2.
#endif
#ifndef NESTED_PLIC_PAIRED_NORMAL_COSINE_MIN
# define NESTED_PLIC_PAIRED_NORMAL_COSINE_MIN 0.999
#endif
#ifndef NESTED_PLIC_NEIGHBOR_CONSENSUS_COSINE_MIN
# define NESTED_PLIC_NEIGHBOR_CONSENSUS_COSINE_MIN 0.96
#endif
#ifndef NESTED_PLIC_NEIGHBOR_CONSENSUS_COUNT_MIN
# define NESTED_PLIC_NEIGHBOR_CONSENSUS_COUNT_MIN 2
#endif
#ifndef NESTED_PLIC_CENTER_PAIRED_DOMINANCE_RATIO_MIN
# define NESTED_PLIC_CENTER_PAIRED_DOMINANCE_RATIO_MIN 2.
#endif
#ifndef NESTED_PLIC_CENTER_PAIRED_COSINE_MIN
# define NESTED_PLIC_CENTER_PAIRED_COSINE_MIN \
  NESTED_PLIC_NEIGHBOR_CONSENSUS_COSINE_MIN
#endif
#ifndef NESTED_PLIC_ENDPOINT_GRADIENT_DOMINANCE_RATIO
# define NESTED_PLIC_ENDPOINT_GRADIENT_DOMINANCE_RATIO 2.
#endif

/*
 * Cuts at the same pure endpoint within roundoff, and two unresolved cuts at
 * opposite endpoints of an almost-pure intermediate-material cell, do not
 * define one shared physical orientation. Their stored fractions and normals
 * remain untouched; the caller skips only common-normal reconciliation and
 * still applies the strict face-partition gate. Non-finite, out-of-range,
 * one-sided or resolved cuts continue through strict reconstruction.
 */
static inline bool nested_plic_requires_common_normal (
  double inner_fraction, double outer_fraction)
{
  if (!isfinite(inner_fraction) || !isfinite(outer_fraction) ||
      inner_fraction < 0. || outer_fraction > 1. ||
      inner_fraction > outer_fraction + 64.*DBL_EPSILON)
    return true;
  const double tolerance = 64.*DBL_EPSILON;
  const bool both_empty =
    inner_fraction <= tolerance && outer_fraction <= tolerance;
  const bool both_full =
    inner_fraction >= 1. - tolerance && outer_fraction >= 1. - tolerance;
  const double inner_quality =
    4.*inner_fraction*(1. - inner_fraction);
  const double outer_quality =
    4.*outer_fraction*(1. - outer_fraction);
  const bool separated_unresolved_cuts =
    inner_fraction < 0.5 && outer_fraction > 0.5 &&
    inner_quality < (double)NESTED_PLIC_NEIGHBOR_QUALITY_MIN &&
    outer_quality < (double)NESTED_PLIC_NEIGHBOR_QUALITY_MIN;
  return !(both_empty || both_full || separated_unresolved_cuts);
}

/*
 * Return true when the original inner PLIC half-plane is already a subset of
 * the original outer PLIC half-plane inside the unit cell.  In that case every
 * swept subregion also preserves inner <= outer, so forcing one common normal
 * would discard a valid pair of spatially separated cuts.  Non-finite or
 * degenerate inputs fail closed and continue through common-normal recovery.
 */
static inline bool nested_plic_independent_cuts_are_nested (
  double inner_x, double inner_y, double inner_alpha,
  double outer_x, double outer_y, double outer_alpha)
{
  const double values[6] = {
    inner_x, inner_y, inner_alpha,
    outer_x, outer_y, outer_alpha
  };
  for (int index = 0; index < 6; index++)
    if (!isfinite(values[index]))
      return false;
  const double inner_norm = fabs(inner_x) + fabs(inner_y);
  const double outer_norm = fabs(outer_x) + fabs(outer_y);
  if (inner_norm <= 64.*DBL_EPSILON ||
      outer_norm <= 64.*DBL_EPSILON)
    return false;

  const double tolerance = 256.*DBL_EPSILON*(
    1. + inner_norm + outer_norm +
    fabs(inner_alpha) + fabs(outer_alpha));
  const double corners[4][2] = {
    {-0.5, -0.5}, {0.5, -0.5}, {0.5, 0.5}, {-0.5, 0.5}
  };
  const int edge_start[4] = {0, 1, 2, 3};
  const int edge_end[4] = {1, 2, 3, 0};
  int candidate_count = 0;

  for (int corner = 0; corner < 4; corner++) {
    const double x = corners[corner][0], y = corners[corner][1];
    const double inner_signed = inner_x*x + inner_y*y - inner_alpha;
    if (inner_signed <= tolerance) {
      candidate_count++;
      if (outer_x*x + outer_y*y - outer_alpha > tolerance)
        return false;
    }
  }

  for (int edge = 0; edge < 4; edge++) {
    const double x0 = corners[edge_start[edge]][0];
    const double y0 = corners[edge_start[edge]][1];
    const double x1 = corners[edge_end[edge]][0];
    const double y1 = corners[edge_end[edge]][1];
    const double signed0 = inner_x*x0 + inner_y*y0 - inner_alpha;
    const double signed1 = inner_x*x1 + inner_y*y1 - inner_alpha;
    const double denominator = signed0 - signed1;
    if (fabs(denominator) <= tolerance ||
        !((signed0 <= tolerance && signed1 >= -tolerance) ||
          (signed1 <= tolerance && signed0 >= -tolerance)))
      continue;
    const double position = signed0/denominator;
    if (position < -tolerance || position > 1. + tolerance)
      continue;
    const double x = x0 + position*(x1 - x0);
    const double y = y0 + position*(y1 - y0);
    candidate_count++;
    if (outer_x*x + outer_y*y - outer_alpha > tolerance)
      return false;
  }
  return candidate_count > 0;
}

static inline double nested_plic_direction_cosine (
  double first_x, double first_y, double second_x, double second_y)
{
  double first_norm = hypot(first_x, first_y);
  double second_norm = hypot(second_x, second_y);
  if (!isfinite(first_norm) || !isfinite(second_norm) ||
      first_norm <= 64.*DBL_EPSILON || second_norm <= 64.*DBL_EPSILON)
    return -1.;
  return (first_x*second_x + first_y*second_y)/
         (first_norm*second_norm);
}

/*
 * Select one L1-normalized PLIC normal for two nested cuts in the same cell.
 * Cuts nearest one-half receive the largest weight because their normals are
 * less sensitive to volume-fraction noise than sliver cuts.
 */
static inline NestedPlicCommonNormal nested_plic_common_normal (
  double inner_fraction, double outer_fraction,
  double inner_x, double inner_y,
  double outer_x, double outer_y)
{
  NestedPlicCommonNormal common = {0};
  common.fallback_neighbor = -1;
  common.reason = NESTED_PLIC_NONFINITE_INPUT;
  double values[6] = {
    inner_fraction, outer_fraction,
    inner_x, inner_y, outer_x, outer_y
  };
  for (int index = 0; index < 6; index++)
    if (!isfinite(values[index]))
      return common;

  common.inner_norm_l1 = fabs(inner_x) + fabs(inner_y);
  common.outer_norm_l1 = fabs(outer_x) + fabs(outer_y);
  common.inner_raw_weight = 4.*inner_fraction*(1. - inner_fraction);
  common.outer_raw_weight = 4.*outer_fraction*(1. - outer_fraction);
  common.inner_weight = fmax(common.inner_raw_weight, 64.*DBL_EPSILON);
  common.outer_weight = fmax(common.outer_raw_weight, 64.*DBL_EPSILON);
  bool inner_resolved = common.inner_norm_l1 > 64.*DBL_EPSILON;
  bool outer_resolved = common.outer_norm_l1 > 64.*DBL_EPSILON;
  if (!inner_resolved && !outer_resolved) {
    common.reason = NESTED_PLIC_UNRESOLVED_NORMALS;
    return common;
  }
  if (inner_resolved)
    inner_x /= common.inner_norm_l1, inner_y /= common.inner_norm_l1;
  if (outer_resolved)
    outer_x /= common.outer_norm_l1, outer_y /= common.outer_norm_l1;
  common.inner_x = inner_x;
  common.inner_y = inner_y;
  common.outer_x = outer_x;
  common.outer_y = outer_y;
  common.normal_dot = inner_x*outer_x + inner_y*outer_y;

  bool use_outer = outer_resolved &&
    common.inner_raw_weight <= (double)NESTED_PLIC_SLIVER_QUALITY_MAX &&
    common.outer_raw_weight >
      (double)NESTED_PLIC_SLIVER_DOMINANCE_RATIO*common.inner_raw_weight;
  bool use_inner = inner_resolved &&
    common.outer_raw_weight <= (double)NESTED_PLIC_SLIVER_QUALITY_MAX &&
    common.inner_raw_weight >
      (double)NESTED_PLIC_SLIVER_DOMINANCE_RATIO*common.outer_raw_weight;
  if (use_outer || use_inner) {
    common.x = use_outer ? outer_x : inner_x;
    common.y = use_outer ? outer_y : inner_y;
    common.valid = 1;
    common.used_sliver_fallback = 1;
    common.fallback_uses_outer = use_outer;
    common.fallback_support_count = 1;
    common.fallback_quality = use_outer ?
      common.outer_raw_weight : common.inner_raw_weight;
    common.reason = NESTED_PLIC_SLIVER_FALLBACK;
    return common;
  }
  if (!inner_resolved || !outer_resolved) {
    common.reason = NESTED_PLIC_UNRESOLVED_NORMALS;
    return common;
  }
  /* Equal, finite endpoint slivers carry no intermediate-material inventory
     in this cell.  Reconcile only the narrow quality band in which both PLIC
     cuts are geometrically present but still below resolved-neighbor quality.
     Machine-level slivers and resolved equal-volume cuts retain the existing
     fail-closed policy.  The caller still applies the strict face-partition
     gate to the selected non-cancelling center direction. */
  bool coincident_inventory = fabs(outer_fraction - inner_fraction) <=
    64.*DBL_EPSILON;
  bool coincident_endpoint = coincident_inventory &&
    ((inner_fraction <= (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE &&
      outer_fraction <= (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE) ||
     (inner_fraction >=
        1. - (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE &&
      outer_fraction >=
        1. - (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE));
  bool finite_endpoint_sliver = coincident_endpoint &&
    common.inner_raw_weight > (double)NESTED_PLIC_SLIVER_QUALITY_MAX &&
    common.outer_raw_weight > (double)NESTED_PLIC_SLIVER_QUALITY_MAX &&
    common.inner_raw_weight < (double)NESTED_PLIC_NEIGHBOR_QUALITY_MIN &&
    common.outer_raw_weight < (double)NESTED_PLIC_NEIGHBOR_QUALITY_MIN;
  if (finite_endpoint_sliver && common.normal_dot <= 0.) {
    double candidate_x = common.inner_weight*inner_x +
      common.outer_weight*outer_x;
    double candidate_y = common.inner_weight*inner_y +
      common.outer_weight*outer_y;
    double candidate_norm = fabs(candidate_x) + fabs(candidate_y);
    double required_norm = fmax(common.inner_weight, common.outer_weight);
    if (isfinite(candidate_norm) &&
        candidate_norm >= required_norm*(1. - 64.*DBL_EPSILON)) {
      common.x = candidate_x/candidate_norm;
      common.y = candidate_y/candidate_norm;
      common.valid = 1;
      common.used_sliver_fallback = 1;
      common.fallback_uses_outer = -1;
      common.fallback_support_count = 2;
      common.fallback_quality = candidate_norm/
        (common.inner_weight + common.outer_weight);
      common.reason = NESTED_PLIC_COINCIDENT_INVENTORY_BLEND;
      return common;
    }
  }
  if (common.normal_dot <= 0.) {
    common.reason = NESTED_PLIC_OPPOSING_NORMALS;
    return common;
  }

  common.x = common.inner_weight*inner_x + common.outer_weight*outer_x;
  common.y = common.inner_weight*inner_y + common.outer_weight*outer_y;
  double common_norm = fabs(common.x) + fabs(common.y);
  if (common_norm <= 64.*DBL_EPSILON || !isfinite(common_norm)) {
    common.x = common.y = 0.;
    common.reason = NESTED_PLIC_DEGENERATE_BLEND;
    return common;
  }
  common.x /= common_norm;
  common.y /= common_norm;
  common.valid = 1;
  common.reason = NESTED_PLIC_ALIGNED_BLEND;
  return common;
}

/*
 * A cell in which a nested cut is a sliver has no reliable center-cell
 * orientation. Prefer a finite, properly nested, genuinely mixed one-ring
 * neighbor. If both center cuts and neighboring cuts are machine-level
 * slivers, require paired inner/outer agreement in at least two directionally
 * consistent neighbors. The caller still recomputes both plane constants and
 * applies the strict face-partition gate.
 */
static inline NestedPlicCommonNormal
nested_plic_common_normal_with_neighbors (
  double inner_fraction, double outer_fraction,
  double inner_x, double inner_y,
  double outer_x, double outer_y,
  const NestedPlicNeighbor neighbors[NESTED_PLIC_NEIGHBOR_COUNT])
{
  NestedPlicCommonNormal common = nested_plic_common_normal(
    inner_fraction, outer_fraction,
    inner_x, inner_y, outer_x, outer_y);
  if (common.valid || common.reason == NESTED_PLIC_NONFINITE_INPUT)
    return common;
  if (common.reason != NESTED_PLIC_OPPOSING_NORMALS &&
      common.reason != NESTED_PLIC_UNRESOLVED_NORMALS)
    return common;
  if (inner_fraction < 0. || outer_fraction > 1. ||
      inner_fraction > outer_fraction + 64.*DBL_EPSILON)
    return common;
  double weak_quality = fmin(
    common.inner_raw_weight, common.outer_raw_weight);
  double strong_quality = fmax(
    common.inner_raw_weight, common.outer_raw_weight);
  bool inner_center_weak =
    inner_fraction <= (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE ||
    inner_fraction >= 1. - (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE;
  bool outer_center_weak =
    outer_fraction <= (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE ||
    outer_fraction >= 1. - (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE;
  bool both_weak_slivers = inner_center_weak && outer_center_weak;
  bool weak_uses_outer =
    common.outer_raw_weight < common.inner_raw_weight;
  bool weak_center_is_sliver =
    weak_uses_outer ? outer_center_weak : inner_center_weak;
  bool weak_side_dominated =
    weak_center_is_sliver &&
    strong_quality >
      (double)NESTED_PLIC_SLIVER_DOMINANCE_RATIO*weak_quality;
  bool weak_side_below_neighbor_resolution =
    weak_quality < (double)NESTED_PLIC_NEIGHBOR_QUALITY_MIN;
  bool strong_side_dominates = strong_quality >
    (double)NESTED_PLIC_SLIVER_DOMINANCE_RATIO*weak_quality;
  bool needs_neighbor_confirmation =
    (weak_center_is_sliver && strong_quality > weak_quality) ||
    weak_side_below_neighbor_resolution || strong_side_dominates;
  /* Opposing center cuts cannot select an orientation by local quality alone.
     Whenever one side is resolved, require an aligned two-sample cluster from
     the center and one-ring data. This also covers similarly resolved cuts for
     which neither side satisfies the sliver-dominance ratio. */
  bool require_resolved_neighbor_consensus =
    common.reason == NESTED_PLIC_OPPOSING_NORMALS &&
    !both_weak_slivers && !weak_side_dominated &&
    (needs_neighbor_confirmation ||
     strong_quality >= (double)NESTED_PLIC_NEIGHBOR_QUALITY_MIN);
  if (!both_weak_slivers && !weak_side_dominated &&
      !require_resolved_neighbor_consensus)
    return common;
  int strong_uses_outer =
    common.outer_raw_weight > common.inner_raw_weight;

  double best_quality = 0., best_x = 0., best_y = 0.;
  int best_neighbor = -1, best_uses_outer = 0;
  double supported_x[NESTED_PLIC_NEIGHBOR_COUNT] = {0};
  double supported_y[NESTED_PLIC_NEIGHBOR_COUNT] = {0};
  double supported_quality[NESTED_PLIC_NEIGHBOR_COUNT] = {0};
  bool supported_valid[NESTED_PLIC_NEIGHBOR_COUNT] = {false};
  for (int neighbor = 0; neighbor < NESTED_PLIC_NEIGHBOR_COUNT; neighbor++) {
    const NestedPlicNeighbor * sample = &neighbors[neighbor];
    double values[6] = {
      sample->inner_fraction, sample->outer_fraction,
      sample->inner_x, sample->inner_y,
      sample->outer_x, sample->outer_y
    };
    bool finite = true;
    for (int index = 0; index < 6; index++)
      finite = finite && isfinite(values[index]);
    if (!finite || sample->inner_fraction < 0. ||
        sample->outer_fraction > 1. ||
        sample->inner_fraction >
          sample->outer_fraction + 64.*DBL_EPSILON)
      continue;

    for (int uses_outer = 0; uses_outer <= 1; uses_outer++) {
      if (!both_weak_slivers && uses_outer != strong_uses_outer)
        continue;
      double fraction = uses_outer ?
        sample->outer_fraction : sample->inner_fraction;
      double candidate_x = uses_outer ? sample->outer_x : sample->inner_x;
      double candidate_y = uses_outer ? sample->outer_y : sample->inner_y;
      double quality = 4.*fraction*(1. - fraction);
      double norm = fabs(candidate_x) + fabs(candidate_y);
      if (quality < (double)NESTED_PLIC_NEIGHBOR_QUALITY_MIN ||
          norm <= 64.*DBL_EPSILON)
        continue;
      candidate_x /= norm;
      candidate_y /= norm;
      double center_x = uses_outer ? common.outer_x : common.inner_x;
      double center_y = uses_outer ? common.outer_y : common.inner_y;
      double center_norm = fabs(center_x) + fabs(center_y);
      if (center_norm > 64.*DBL_EPSILON &&
          candidate_x*center_x + candidate_y*center_y <= 0.)
        continue;
      if (require_resolved_neighbor_consensus) {
        supported_valid[neighbor] = true;
        supported_x[neighbor] = candidate_x;
        supported_y[neighbor] = candidate_y;
        supported_quality[neighbor] = quality;
      }
      if (quality > best_quality) {
        best_quality = quality;
        best_x = candidate_x;
        best_y = candidate_y;
        best_neighbor = neighbor;
        best_uses_outer = uses_outer;
      }
    }
  }
  int resolved_support_count = best_neighbor >= 0 ? 1 : 0;
  if (best_neighbor >= 0 && require_resolved_neighbor_consensus) {
    double center_x = strong_uses_outer ? common.outer_x : common.inner_x;
    double center_y = strong_uses_outer ? common.outer_y : common.inner_y;
    double center_norm = fabs(center_x) + fabs(center_y);
    bool center_valid =
      strong_quality >= (double)NESTED_PLIC_NEIGHBOR_QUALITY_MIN &&
      center_norm > 64.*DBL_EPSILON;
    if (center_valid) {
      center_x /= center_norm;
      center_y /= center_norm;
    }

    /* Select the largest anchor-consistent support cluster. Anchoring the
       consensus to the single highest-quality sample can discard a coherent
       center-plus-neighbor pair when an isolated mixed neighbor has higher
       fraction quality. Support count is therefore primary and total quality
       is only the deterministic tie-breaker. */
    int selected_count = 0, selected_neighbor = -1;
    double selected_weight = -1.;
    double selected_x = 0., selected_y = 0.;
    double selected_minimum_cosine = 1.;
    for (int anchor = -1; anchor < NESTED_PLIC_NEIGHBOR_COUNT; anchor++) {
      if ((anchor < 0 && !center_valid) ||
          (anchor >= 0 && !supported_valid[anchor]))
        continue;
      double anchor_x = anchor < 0 ? center_x : supported_x[anchor];
      double anchor_y = anchor < 0 ? center_y : supported_y[anchor];
      int cluster_count = 0, cluster_neighbor = -1;
      double cluster_weight = 0., cluster_x = 0., cluster_y = 0.;
      double cluster_minimum_cosine = 1.;
      if (center_valid) {
        double cosine = nested_plic_direction_cosine(
          anchor_x, anchor_y, center_x, center_y);
        if (cosine >=
            (double)NESTED_PLIC_NEIGHBOR_CONSENSUS_COSINE_MIN) {
          cluster_count++;
          cluster_minimum_cosine = fmin(cluster_minimum_cosine, cosine);
          cluster_x += strong_quality*center_x;
          cluster_y += strong_quality*center_y;
          cluster_weight += strong_quality;
        }
      }
      for (int neighbor = 0; neighbor < NESTED_PLIC_NEIGHBOR_COUNT;
           neighbor++) {
        if (!supported_valid[neighbor])
          continue;
        double cosine = nested_plic_direction_cosine(
          anchor_x, anchor_y,
          supported_x[neighbor], supported_y[neighbor]);
        if (cosine < (double)NESTED_PLIC_NEIGHBOR_CONSENSUS_COSINE_MIN)
          continue;
        cluster_count++;
        cluster_minimum_cosine = fmin(cluster_minimum_cosine, cosine);
        cluster_x += supported_quality[neighbor]*supported_x[neighbor];
        cluster_y += supported_quality[neighbor]*supported_y[neighbor];
        cluster_weight += supported_quality[neighbor];
        if (cluster_neighbor < 0)
          cluster_neighbor = neighbor;
      }
      if (cluster_count > selected_count ||
          (cluster_count == selected_count &&
           cluster_weight > selected_weight)) {
        selected_count = cluster_count;
        selected_neighbor = cluster_neighbor;
        selected_weight = cluster_weight;
        selected_x = cluster_x;
        selected_y = cluster_y;
        selected_minimum_cosine = cluster_minimum_cosine;
      }
    }
    resolved_support_count = selected_count;
    if (resolved_support_count >=
          (int)NESTED_PLIC_NEIGHBOR_CONSENSUS_COUNT_MIN &&
        selected_neighbor >= 0 && selected_weight > 64.*DBL_EPSILON) {
      double consensus_norm = fabs(selected_x) + fabs(selected_y);
      if (!isfinite(consensus_norm) || consensus_norm <= 64.*DBL_EPSILON)
        return common;
      best_x = selected_x/consensus_norm;
      best_y = selected_y/consensus_norm;
      best_quality = selected_minimum_cosine;
      best_neighbor = selected_neighbor;
      best_uses_outer = strong_uses_outer;
    }
    else {
      /* The subordinate cut is below the declared normal-resolution floor,
         while the other cut is resolved by the existing dominance test.
         Preserve the resolved side using its best finite one-ring sample;
         both reconstructed cuts still pass the caller's strict face gate. */
      bool dominant_resolved_side =
        weak_side_below_neighbor_resolution &&
        strong_quality >
          (double)NESTED_PLIC_SLIVER_DOMINANCE_RATIO*weak_quality;
      if (!dominant_resolved_side)
        return common;
    }
  }
  if (best_neighbor < 0 && both_weak_slivers &&
      strong_quality <= (double)NESTED_PLIC_SLIVER_QUALITY_MAX) {
    bool center_low =
      inner_fraction <= (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE &&
      outer_fraction <= (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE;
    bool center_high =
      inner_fraction >= 1. -
        (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE &&
      outer_fraction >= 1. -
        (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE;
    if (center_low || center_high) {
      double paired_x[NESTED_PLIC_NEIGHBOR_COUNT] = {0};
      double paired_y[NESTED_PLIC_NEIGHBOR_COUNT] = {0};
      double paired_quality[NESTED_PLIC_NEIGHBOR_COUNT] = {0};
      bool paired_valid[NESTED_PLIC_NEIGHBOR_COUNT] = {false};
      int paired_best_neighbor = -1;
      for (int neighbor = 0; neighbor < NESTED_PLIC_NEIGHBOR_COUNT;
           neighbor++) {
        const NestedPlicNeighbor * sample = &neighbors[neighbor];
        double values[6] = {
          sample->inner_fraction, sample->outer_fraction,
          sample->inner_x, sample->inner_y,
          sample->outer_x, sample->outer_y
        };
        bool finite = true;
        for (int index = 0; index < 6; index++)
          finite = finite && isfinite(values[index]);
        bool nested = sample->inner_fraction >= 0. &&
          sample->outer_fraction <= 1. &&
          sample->inner_fraction <=
            sample->outer_fraction + 64.*DBL_EPSILON;
        bool same_low_state = center_low &&
          sample->inner_fraction > 0. && sample->outer_fraction > 0. &&
          sample->inner_fraction <=
            (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE &&
          sample->outer_fraction <=
            (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE;
        bool same_high_state = center_high &&
          sample->inner_fraction < 1. && sample->outer_fraction < 1. &&
          sample->inner_fraction >= 1. -
            (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE &&
          sample->outer_fraction >= 1. -
            (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE;
        if (!finite || !nested || (!same_low_state && !same_high_state))
          continue;

        double inner_norm = fabs(sample->inner_x) + fabs(sample->inner_y);
        double outer_norm = fabs(sample->outer_x) + fabs(sample->outer_y);
        double paired_inner_quality = 4.*sample->inner_fraction*
          (1. - sample->inner_fraction);
        double paired_outer_quality = 4.*sample->outer_fraction*
          (1. - sample->outer_fraction);
        if (inner_norm <= 64.*DBL_EPSILON ||
            outer_norm <= 64.*DBL_EPSILON ||
            fmax(paired_inner_quality, paired_outer_quality) >
              (double)NESTED_PLIC_SLIVER_QUALITY_MAX)
          continue;
        double paired_inner_x = sample->inner_x/inner_norm;
        double paired_inner_y = sample->inner_y/inner_norm;
        double paired_outer_x = sample->outer_x/outer_norm;
        double paired_outer_y = sample->outer_y/outer_norm;
        if (nested_plic_direction_cosine(
              paired_inner_x, paired_inner_y,
              paired_outer_x, paired_outer_y) <
            (double)NESTED_PLIC_PAIRED_NORMAL_COSINE_MIN)
          continue;
        double candidate_x = paired_inner_x + paired_outer_x;
        double candidate_y = paired_inner_y + paired_outer_y;
        double candidate_norm = fabs(candidate_x) + fabs(candidate_y);
        if (!isfinite(candidate_norm) ||
            candidate_norm <= 64.*DBL_EPSILON)
          continue;
        paired_x[neighbor] = candidate_x/candidate_norm;
        paired_y[neighbor] = candidate_y/candidate_norm;
        paired_quality[neighbor] = fmax(
          paired_inner_quality, paired_outer_quality);
        paired_valid[neighbor] = true;
        if (paired_best_neighbor < 0 ||
            paired_quality[neighbor] >
              paired_quality[paired_best_neighbor])
          paired_best_neighbor = neighbor;
      }
      if (paired_best_neighbor >= 0) {
        int support_count = 0;
        double minimum_cosine = 1.;
        for (int neighbor = 0; neighbor < NESTED_PLIC_NEIGHBOR_COUNT;
             neighbor++) {
          if (!paired_valid[neighbor])
            continue;
          double cosine = nested_plic_direction_cosine(
            paired_x[paired_best_neighbor], paired_y[paired_best_neighbor],
            paired_x[neighbor], paired_y[neighbor]);
          if (cosine >=
              (double)NESTED_PLIC_NEIGHBOR_CONSENSUS_COSINE_MIN) {
            support_count++;
            minimum_cosine = fmin(minimum_cosine, cosine);
          }
        }
        if (support_count >=
            (int)NESTED_PLIC_NEIGHBOR_CONSENSUS_COUNT_MIN) {
          common.x = paired_x[paired_best_neighbor];
          common.y = paired_y[paired_best_neighbor];
          common.valid = 1;
          common.used_sliver_fallback = 1;
          common.fallback_neighbor = paired_best_neighbor;
          common.fallback_uses_outer = -1;
          common.fallback_support_count = support_count;
          common.fallback_quality = minimum_cosine;
          common.reason = NESTED_PLIC_PAIRED_NEIGHBOR_CONSENSUS;
          return common;
        }
      }
    }
  }
  if (best_neighbor < 0) {
    /* A center cut can lie between the machine-sliver and resolved-neighbor
       quality bands. Accept its better-resolved side only when a genuinely
       nested, same-endpoint neighbor supplies an independently paired normal
       in the same direction. The center and paired neighbor are the two
       required supports; unsupported endpoint cuts remain hard failures. */
    bool center_low =
      inner_fraction <= (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE &&
      outer_fraction <= (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE;
    bool center_high =
      inner_fraction >= 1. -
        (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE &&
      outer_fraction >= 1. -
        (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE;
    bool center_side_is_dominant = weak_quality > 0. &&
      strong_quality >
        (double)NESTED_PLIC_CENTER_PAIRED_DOMINANCE_RATIO_MIN*weak_quality;
    double center_x = strong_uses_outer ? common.outer_x : common.inner_x;
    double center_y = strong_uses_outer ? common.outer_y : common.inner_y;
    double center_norm = fabs(center_x) + fabs(center_y);
    int paired_neighbor = -1;
    double paired_x = 0., paired_y = 0., paired_quality = 0.;
    double paired_alignment = -1.;
    if ((center_low || center_high) && center_side_is_dominant &&
        center_norm > 64.*DBL_EPSILON) {
      center_x /= center_norm;
      center_y /= center_norm;
      for (int neighbor = 0; neighbor < NESTED_PLIC_NEIGHBOR_COUNT;
           neighbor++) {
        const NestedPlicNeighbor * sample = &neighbors[neighbor];
        double values[6] = {
          sample->inner_fraction, sample->outer_fraction,
          sample->inner_x, sample->inner_y,
          sample->outer_x, sample->outer_y
        };
        bool finite = true;
        for (int index = 0; index < 6; index++)
          finite = finite && isfinite(values[index]);
        bool nested = sample->inner_fraction >= 0. &&
          sample->outer_fraction <= 1. &&
          sample->inner_fraction <=
            sample->outer_fraction + 64.*DBL_EPSILON;
        bool same_low_state = center_low &&
          sample->inner_fraction > 0. && sample->outer_fraction > 0. &&
          sample->inner_fraction <=
            (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE &&
          sample->outer_fraction <=
            (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE;
        bool same_high_state = center_high &&
          sample->inner_fraction < 1. && sample->outer_fraction < 1. &&
          sample->inner_fraction >= 1. -
            (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE &&
          sample->outer_fraction >= 1. -
            (double)NESTED_PLIC_CENTER_FRACTION_TOLERANCE;
        if (!finite || !nested || (!same_low_state && !same_high_state))
          continue;
        double inner_norm = fabs(sample->inner_x) + fabs(sample->inner_y);
        double outer_norm = fabs(sample->outer_x) + fabs(sample->outer_y);
        if (inner_norm <= 64.*DBL_EPSILON ||
            outer_norm <= 64.*DBL_EPSILON)
          continue;
        double inner_x = sample->inner_x/inner_norm;
        double inner_y = sample->inner_y/inner_norm;
        double outer_x = sample->outer_x/outer_norm;
        double outer_y = sample->outer_y/outer_norm;
        if (nested_plic_direction_cosine(
              inner_x, inner_y, outer_x, outer_y) <
            (double)NESTED_PLIC_PAIRED_NORMAL_COSINE_MIN)
          continue;
        double candidate_x = inner_x + outer_x;
        double candidate_y = inner_y + outer_y;
        double candidate_norm = fabs(candidate_x) + fabs(candidate_y);
        if (!isfinite(candidate_norm) || candidate_norm <= 64.*DBL_EPSILON)
          continue;
        candidate_x /= candidate_norm;
        candidate_y /= candidate_norm;
        double alignment = nested_plic_direction_cosine(
          center_x, center_y, candidate_x, candidate_y);
        if (alignment < (double)NESTED_PLIC_CENTER_PAIRED_COSINE_MIN)
          continue;
        double quality = fmax(
          4.*sample->inner_fraction*(1. - sample->inner_fraction),
          4.*sample->outer_fraction*(1. - sample->outer_fraction));
        if (paired_neighbor < 0 || quality > paired_quality) {
          paired_neighbor = neighbor;
          paired_x = candidate_x;
          paired_y = candidate_y;
          paired_quality = quality;
          paired_alignment = alignment;
        }
      }
    }
    if (paired_neighbor >= 0) {
      common.x = strong_quality*center_x + paired_quality*paired_x;
      common.y = strong_quality*center_y + paired_quality*paired_y;
      double consensus_norm = fabs(common.x) + fabs(common.y);
      if (isfinite(consensus_norm) && consensus_norm > 64.*DBL_EPSILON) {
        common.x /= consensus_norm;
        common.y /= consensus_norm;
        common.valid = 1;
        common.used_sliver_fallback = 1;
        common.fallback_neighbor = paired_neighbor;
        common.fallback_uses_outer = strong_uses_outer;
        common.fallback_support_count = 2;
        common.fallback_quality = paired_alignment;
        common.reason = NESTED_PLIC_CENTER_PAIRED_NEIGHBOR_CONSENSUS;
        return common;
      }
      common.x = common.y = 0.;
    }

    /* At a common endpoint, two PLIC normals within the declared endpoint
       fraction band may be dominated by fraction noise and point in opposite
       directions even when their weights exceed the machine-sliver band.
       After all resolved-normal consensus paths fail, use the centered
       gradient of the mean nested occupancy only when its one-ring signal
       dominates the center endpoint deficit. The caller still reconstructs
       both cuts and applies the strict face-partition gate. */
    if ((center_low || center_high) && both_weak_slivers) {
      double mean_fraction[NESTED_PLIC_NEIGHBOR_COUNT] = {0};
      bool gradient_inputs_valid = true;
      for (int neighbor = 0; neighbor < NESTED_PLIC_NEIGHBOR_COUNT;
           neighbor++) {
        const NestedPlicNeighbor * sample = &neighbors[neighbor];
        gradient_inputs_valid = gradient_inputs_valid &&
          isfinite(sample->inner_fraction) &&
          isfinite(sample->outer_fraction) &&
          sample->inner_fraction >= 0. &&
          sample->outer_fraction <= 1. &&
          sample->inner_fraction <=
            sample->outer_fraction + 64.*DBL_EPSILON;
        mean_fraction[neighbor] =
          .5*(sample->inner_fraction + sample->outer_fraction);
      }
      double gradient_x = mean_fraction[1] - mean_fraction[0];
      double gradient_y = mean_fraction[3] - mean_fraction[2];
      double gradient_norm = fabs(gradient_x) + fabs(gradient_y);
      double center_mean = .5*(inner_fraction + outer_fraction);
      double endpoint_deficit = center_low ? center_mean : 1. - center_mean;
      double signal_floor = fmax(64.*DBL_EPSILON,
        (double)NESTED_PLIC_ENDPOINT_GRADIENT_DOMINANCE_RATIO*
        endpoint_deficit);
      if (gradient_inputs_valid && isfinite(gradient_norm) &&
          gradient_norm > signal_floor) {
        common.x = gradient_x/gradient_norm;
        common.y = gradient_y/gradient_norm;
        common.valid = 1;
        common.used_sliver_fallback = 1;
        common.fallback_neighbor = -1;
        common.fallback_uses_outer = -1;
        common.fallback_support_count = 2;
        common.fallback_quality = gradient_norm;
        common.reason = NESTED_PLIC_ENDPOINT_GRADIENT_FALLBACK;
        return common;
      }
    }

    /* If no resolved one-ring or paired consensus exists, negligible or
       machine-level same-endpoint inventory still needs one transport
       orientation. Use the better-resolved center side only when its PLIC
       support dominates the other side; both cell fractions are retained and
       the caller still applies the strict face-partition gate. */
    const double endpoint_tolerance =
      (double)NESTED_PLIC_NEGLIGIBLE_ENDPOINT_TOLERANCE;
    bool machine_endpoint = (center_low || center_high) &&
      both_weak_slivers &&
      strong_quality <= (double)NESTED_PLIC_SLIVER_QUALITY_MAX;
    bool same_low_endpoint =
      center_low &&
      ((inner_fraction <= endpoint_tolerance &&
        outer_fraction <= endpoint_tolerance) || machine_endpoint);
    bool same_high_endpoint =
      center_high &&
      ((inner_fraction >= 1. - endpoint_tolerance &&
        outer_fraction >= 1. - endpoint_tolerance) || machine_endpoint);
    double endpoint_weak_quality = fmin(
      common.inner_raw_weight, common.outer_raw_weight);
    double endpoint_strong_quality = fmax(
      common.inner_raw_weight, common.outer_raw_weight);
    bool endpoint_side_is_dominant = endpoint_strong_quality >
      (double)NESTED_PLIC_NEGLIGIBLE_ENDPOINT_DOMINANCE_RATIO*
      endpoint_weak_quality;
    if ((same_low_endpoint || same_high_endpoint) &&
        endpoint_side_is_dominant) {
      bool use_outer_endpoint = common.outer_norm_l1 > 64.*DBL_EPSILON &&
        (common.inner_norm_l1 <= 64.*DBL_EPSILON ||
         common.outer_raw_weight > common.inner_raw_weight ||
         (common.outer_raw_weight == common.inner_raw_weight &&
          same_low_endpoint));
      common.x = use_outer_endpoint ? common.outer_x : common.inner_x;
      common.y = use_outer_endpoint ? common.outer_y : common.inner_y;
      common.valid = 1;
      common.used_sliver_fallback = 1;
      common.fallback_uses_outer = use_outer_endpoint;
      common.fallback_support_count = 1;
      common.fallback_quality = use_outer_endpoint ?
        common.outer_raw_weight : common.inner_raw_weight;
      common.reason = NESTED_PLIC_NEGLIGIBLE_ENDPOINT_FALLBACK;
      return common;
    }

    /* At a machine-level common endpoint, the PLIC normals can be mutually
       transverse even though both cuts carry negligible interfacial support.
       If every stronger neighborhood policy above is unavailable, retain the
       two center directions only when their fraction-weighted resultant does
       not cancel either input's resolved magnitude. This preserves both cell
       inventories and leaves the strict face-partition gate unchanged; truly
       opposing or underdetermined endpoint normals remain hard failures. */
    bool both_center_normals_resolved =
      common.inner_norm_l1 > 64.*DBL_EPSILON &&
      common.outer_norm_l1 > 64.*DBL_EPSILON;
    bool both_center_fractions_resolved =
      common.inner_raw_weight > 64.*DBL_EPSILON &&
      common.outer_raw_weight > 64.*DBL_EPSILON;
    if (machine_endpoint && both_center_normals_resolved &&
        both_center_fractions_resolved) {
      double candidate_x = common.inner_weight*common.inner_x +
        common.outer_weight*common.outer_x;
      double candidate_y = common.inner_weight*common.inner_y +
        common.outer_weight*common.outer_y;
      double candidate_norm = fabs(candidate_x) + fabs(candidate_y);
      double required_norm = fmax(
        common.inner_weight, common.outer_weight);
      if (isfinite(candidate_norm) &&
          candidate_norm >= required_norm*(1. - 64.*DBL_EPSILON)) {
        common.x = candidate_x/candidate_norm;
        common.y = candidate_y/candidate_norm;
        common.valid = 1;
        common.used_sliver_fallback = 1;
        common.fallback_neighbor = -1;
        common.fallback_uses_outer = -1;
        common.fallback_support_count = 2;
        common.fallback_quality = candidate_norm/
          (common.inner_weight + common.outer_weight);
        common.reason = NESTED_PLIC_ENDPOINT_CENTER_BLEND;
        return common;
      }
    }
  }
  if (best_neighbor < 0)
    return common;

  common.x = best_x;
  common.y = best_y;
  common.valid = 1;
  common.used_sliver_fallback = 1;
  common.fallback_neighbor = best_neighbor;
  common.fallback_uses_outer = best_uses_outer;
  common.fallback_support_count = resolved_support_count;
  common.fallback_quality = best_quality;
  common.reason = NESTED_PLIC_NEIGHBOR_FALLBACK;
  return common;
}

#endif
