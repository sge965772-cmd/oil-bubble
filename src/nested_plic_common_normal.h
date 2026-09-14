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
  NESTED_PLIC_CENTER_PAIRED_NEIGHBOR_CONSENSUS = 9
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
  bool needs_neighbor_confirmation =
    (weak_center_is_sliver && strong_quality > weak_quality) ||
    weak_side_below_neighbor_resolution;
  /* A dominant center cut is insufficient when the two center normals oppose.
     If the weak cut is below the declared resolved-neighbor quality, require
     two resolved one-ring samples to confirm the strong-side direction. */
  bool require_resolved_neighbor_consensus =
    common.reason == NESTED_PLIC_OPPOSING_NORMALS &&
    !both_weak_slivers && !weak_side_dominated &&
    needs_neighbor_confirmation;
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
    double consensus_x = 0., consensus_y = 0.;
    double consensus_weight = 0., minimum_cosine = 1.;
    resolved_support_count = 0;
    double center_x = strong_uses_outer ? common.outer_x : common.inner_x;
    double center_y = strong_uses_outer ? common.outer_y : common.inner_y;
    double center_norm = fabs(center_x) + fabs(center_y);
    if (strong_quality >= (double)NESTED_PLIC_NEIGHBOR_QUALITY_MIN &&
        center_norm > 64.*DBL_EPSILON) {
      center_x /= center_norm;
      center_y /= center_norm;
      double cosine = nested_plic_direction_cosine(
        best_x, best_y, center_x, center_y);
      if (cosine >=
          (double)NESTED_PLIC_NEIGHBOR_CONSENSUS_COSINE_MIN) {
        resolved_support_count++;
        minimum_cosine = fmin(minimum_cosine, cosine);
        consensus_x += strong_quality*center_x;
        consensus_y += strong_quality*center_y;
        consensus_weight += strong_quality;
      }
    }
    for (int neighbor = 0; neighbor < NESTED_PLIC_NEIGHBOR_COUNT;
         neighbor++) {
      if (!supported_valid[neighbor])
        continue;
      double cosine = nested_plic_direction_cosine(
        best_x, best_y, supported_x[neighbor], supported_y[neighbor]);
      if (cosine < (double)NESTED_PLIC_NEIGHBOR_CONSENSUS_COSINE_MIN)
        continue;
      resolved_support_count++;
      minimum_cosine = fmin(minimum_cosine, cosine);
      consensus_x += supported_quality[neighbor]*supported_x[neighbor];
      consensus_y += supported_quality[neighbor]*supported_y[neighbor];
      consensus_weight += supported_quality[neighbor];
    }
    if (resolved_support_count <
          (int)NESTED_PLIC_NEIGHBOR_CONSENSUS_COUNT_MIN ||
        consensus_weight <= 64.*DBL_EPSILON)
      return common;
    double consensus_norm = fabs(consensus_x) + fabs(consensus_y);
    if (!isfinite(consensus_norm) || consensus_norm <= 64.*DBL_EPSILON)
      return common;
    best_x = consensus_x/consensus_norm;
    best_y = consensus_y/consensus_norm;
    best_quality = minimum_cosine;
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

    /* If no resolved one-ring or paired consensus exists, negligible
       same-endpoint inventory still needs one transport orientation: using
       the better-resolved center side prevents O(1) non-nested swept
       fractions without changing either cell fraction. */
    const double endpoint_tolerance =
      (double)NESTED_PLIC_NEGLIGIBLE_ENDPOINT_TOLERANCE;
    bool same_low_endpoint =
      inner_fraction <= endpoint_tolerance &&
      outer_fraction <= endpoint_tolerance;
    bool same_high_endpoint =
      inner_fraction >= 1. - endpoint_tolerance &&
      outer_fraction >= 1. - endpoint_tolerance;
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
