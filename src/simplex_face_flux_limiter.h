#ifndef SIMPLEX_FACE_FLUX_LIMITER_H
#define SIMPLEX_FACE_FLUX_LIMITER_H

#include <float.h>
#include <math.h>

#define SIMPLEX_FLUX_MATERIALS 5

#ifndef SIMPLEX_LOW_ORDER_BUDGET_REGULARIZATION
# define SIMPLEX_LOW_ORDER_BUDGET_REGULARIZATION 1.e-30
#endif

typedef enum {
  SIMPLEX_FACE_FLUX_OK = 0,
  SIMPLEX_FACE_FLUX_INVALID_SCALAR,
  SIMPLEX_FACE_FLUX_NONFINITE_INPUT,
  SIMPLEX_FACE_FLUX_HIGH_OUT_OF_BOUNDS,
  SIMPLEX_FACE_FLUX_LOW_OUT_OF_BOUNDS,
  SIMPLEX_FACE_FLUX_LEFT_RATIO_OUT_OF_BOUNDS,
  SIMPLEX_FACE_FLUX_RIGHT_RATIO_OUT_OF_BOUNDS,
  SIMPLEX_FACE_FLUX_HIGH_CLOSURE,
  SIMPLEX_FACE_FLUX_LOW_CLOSURE,
  SIMPLEX_FACE_FLUX_LIMITED_OUT_OF_BOUNDS,
  SIMPLEX_FACE_FLUX_LIMITED_CLOSURE
} SimplexFaceFluxFailure;

typedef struct {
  double fraction[SIMPLEX_FLUX_MATERIALS];
  double alpha;
  double correction_linf;
  double closure_defect;
  int limited;
  int valid;
  SimplexFaceFluxFailure failure_reason;
  int failure_material;
} SimplexFaceFlux;

typedef enum {
  SIMPLEX_LOW_ORDER_BUDGET_OK = 0,
  SIMPLEX_LOW_ORDER_BUDGET_NONFINITE_INPUT,
  SIMPLEX_LOW_ORDER_BUDGET_INVALID_COMPRESSION_SHARE,
  SIMPLEX_LOW_ORDER_BUDGET_NEGATIVE
} SimplexLowOrderBudgetFailure;

typedef struct {
  double value;
  double donor_flux_contribution;
  double compression_contribution;
  double left_antidiffusive_contribution;
  double right_antidiffusive_contribution;
  double depletion_ratio;
  int valid;
  SimplexLowOrderBudgetFailure failure_reason;
} SimplexLowOrderBudget;

typedef struct {
  SimplexLowOrderBudget material[SIMPLEX_FLUX_MATERIALS];
  int used_fallback;
  int valid;
} SimplexLowOrderBudgetSet;

/*
 * Evaluate the monotone material budget used by the face-flux limiter.
 * compression_share is the selected material share associated with the
 * directional velocity divergence.  Keeping this calculation behind a
 * pure interface makes a failed Basilisk cell directly replayable in C.
 */
static inline SimplexLowOrderBudget simplex_low_order_budget (
  double initial,
  double left_low_flux, double right_low_flux,
  double left_antidiffusive_flux, double right_antidiffusive_flux,
  double compression_share, double divergence, double scale,
  double tolerance)
{
  SimplexLowOrderBudget budget = {0};
  const double inputs[] = {
    initial, left_low_flux, right_low_flux,
    left_antidiffusive_flux, right_antidiffusive_flux,
    compression_share, divergence, scale, tolerance
  };
  for (unsigned int index = 0;
       index < sizeof(inputs)/sizeof(inputs[0]); index++)
    if (!isfinite(inputs[index])) {
      budget.failure_reason = SIMPLEX_LOW_ORDER_BUDGET_NONFINITE_INPUT;
      return budget;
    }
  if (tolerance < 0. || compression_share < -tolerance ||
      compression_share > 1. + tolerance) {
    budget.failure_reason = SIMPLEX_LOW_ORDER_BUDGET_INVALID_COMPRESSION_SHARE;
    return budget;
  }

  budget.donor_flux_contribution =
    scale*(left_low_flux - right_low_flux);
  budget.compression_contribution =
    scale*compression_share*divergence;
  budget.left_antidiffusive_contribution =
    scale*left_antidiffusive_flux;
  budget.right_antidiffusive_contribution =
    -scale*right_antidiffusive_flux;
  budget.value = initial + budget.donor_flux_contribution +
    budget.compression_contribution;
  if (!isfinite(budget.value) || budget.value < -tolerance) {
    budget.failure_reason = SIMPLEX_LOW_ORDER_BUDGET_NEGATIVE;
    return budget;
  }

  const double negative_antidiffusive =
    fmin(0., budget.left_antidiffusive_contribution) +
    fmin(0., budget.right_antidiffusive_contribution);
  budget.depletion_ratio = negative_antidiffusive < 0. ?
    fmin(1., fmax(0., budget.value)/
         (-negative_antidiffusive +
          (double)SIMPLEX_LOW_ORDER_BUDGET_REGULARIZATION)) :
    1.;
  budget.failure_reason = SIMPLEX_LOW_ORDER_BUDGET_OK;
  budget.valid = 1;
  return budget;
}

/*
 * Evaluate one five-material low-order cell budget.  The requested
 * compression partition is retained bit-for-bit when it is admissible.  If
 * any material would become negative, retry all materials with one supplied
 * simplex-preserving fallback partition.  A caller may then use the returned
 * per-material depletion ratios and the same selected compression partition.
 */
static inline SimplexLowOrderBudgetSet simplex_low_order_budget_set (
  const double initial[SIMPLEX_FLUX_MATERIALS],
  const double left_low_flux[SIMPLEX_FLUX_MATERIALS],
  const double right_low_flux[SIMPLEX_FLUX_MATERIALS],
  const double left_antidiffusive_flux[SIMPLEX_FLUX_MATERIALS],
  const double right_antidiffusive_flux[SIMPLEX_FLUX_MATERIALS],
  const double requested_compression_share[SIMPLEX_FLUX_MATERIALS],
  const double fallback_compression_share[SIMPLEX_FLUX_MATERIALS],
  double divergence, double scale, double tolerance)
{
  SimplexLowOrderBudgetSet set = {0};
  double requested_sum = 0., fallback_sum = 0.;
  int requested_partition_valid = isfinite(tolerance) && tolerance >= 0.;
  int fallback_partition_valid = requested_partition_valid;
  for (int material = 0; material < SIMPLEX_FLUX_MATERIALS; material++) {
    requested_sum += requested_compression_share[material];
    fallback_sum += fallback_compression_share[material];
    if (!isfinite(requested_compression_share[material]) ||
        requested_compression_share[material] < -tolerance ||
        requested_compression_share[material] > 1. + tolerance)
      requested_partition_valid = 0;
    if (!isfinite(fallback_compression_share[material]) ||
        fallback_compression_share[material] < -tolerance ||
        fallback_compression_share[material] > 1. + tolerance)
      fallback_partition_valid = 0;
  }
  if (!isfinite(requested_sum) || fabs(requested_sum - 1.) > tolerance)
    requested_partition_valid = 0;
  if (!isfinite(fallback_sum) || fabs(fallback_sum - 1.) > tolerance)
    fallback_partition_valid = 0;

  set.used_fallback = !requested_partition_valid;
  if (requested_partition_valid)
    for (int material = 0; material < SIMPLEX_FLUX_MATERIALS; material++) {
      set.material[material] = simplex_low_order_budget(
        initial[material],
        left_low_flux[material], right_low_flux[material],
        left_antidiffusive_flux[material],
        right_antidiffusive_flux[material],
        requested_compression_share[material],
        divergence, scale, tolerance);
      if (!set.material[material].valid)
        set.used_fallback = 1;
    }
  if (set.used_fallback)
    for (int material = 0; material < SIMPLEX_FLUX_MATERIALS; material++)
      set.material[material] = simplex_low_order_budget(
        initial[material],
        left_low_flux[material], right_low_flux[material],
        left_antidiffusive_flux[material],
        right_antidiffusive_flux[material],
        fallback_compression_share[material],
        divergence, scale, tolerance);
  set.valid = set.used_fallback ? fallback_partition_valid :
    requested_partition_valid;
  for (int material = 0; material < SIMPLEX_FLUX_MATERIALS; material++)
    if (!set.material[material].valid)
      set.valid = 0;
  return set;
}

static inline double simplex_face_flux_sum (const SimplexFaceFlux * flux)
{
  double sum = 0.;
  for (int material = 0; material < SIMPLEX_FLUX_MATERIALS; material++)
    sum += flux->fraction[material];
  return sum;
}

/*
 * Blend a geometric material face fraction toward a monotone upwind
 * fraction.  One alpha is shared by all materials, so the simplex closure
 * is retained.  The sign of (high - low)*face_velocity identifies the cell
 * depleted by each antidiffusive material flux.
 */
static inline SimplexFaceFlux simplex_face_flux_limit (
  const double high[SIMPLEX_FLUX_MATERIALS],
  const double low[SIMPLEX_FLUX_MATERIALS],
  double face_velocity,
  const double left_depletion_ratio[SIMPLEX_FLUX_MATERIALS],
  const double right_depletion_ratio[SIMPLEX_FLUX_MATERIALS],
  double tolerance)
{
  SimplexFaceFlux limited = {0};
  limited.alpha = 1.;
  limited.failure_material = -1;
  if (!isfinite(face_velocity) || !isfinite(tolerance) || tolerance < 0.) {
    limited.failure_reason = SIMPLEX_FACE_FLUX_INVALID_SCALAR;
    return limited;
  }

  double high_sum = 0., low_sum = 0.;
  for (int material = 0; material < SIMPLEX_FLUX_MATERIALS; material++) {
    double values[] = {
      high[material], low[material],
      left_depletion_ratio[material], right_depletion_ratio[material]
    };
    for (int index = 0; index < 4; index++)
      if (!isfinite(values[index])) {
        limited.failure_reason = SIMPLEX_FACE_FLUX_NONFINITE_INPUT;
        limited.failure_material = material;
        return limited;
      }
    if (high[material] < -tolerance || high[material] > 1. + tolerance) {
      limited.failure_reason = SIMPLEX_FACE_FLUX_HIGH_OUT_OF_BOUNDS;
      limited.failure_material = material;
      return limited;
    }
    if (low[material] < -tolerance || low[material] > 1. + tolerance) {
      limited.failure_reason = SIMPLEX_FACE_FLUX_LOW_OUT_OF_BOUNDS;
      limited.failure_material = material;
      return limited;
    }
    if (left_depletion_ratio[material] < 0. ||
        left_depletion_ratio[material] > 1.) {
      limited.failure_reason = SIMPLEX_FACE_FLUX_LEFT_RATIO_OUT_OF_BOUNDS;
      limited.failure_material = material;
      return limited;
    }
    if (right_depletion_ratio[material] < 0. ||
        right_depletion_ratio[material] > 1.) {
      limited.failure_reason = SIMPLEX_FACE_FLUX_RIGHT_RATIO_OUT_OF_BOUNDS;
      limited.failure_material = material;
      return limited;
    }
    high_sum += high[material];
    low_sum += low[material];
  }
  if (fabs(high_sum - 1.) > tolerance) {
    limited.failure_reason = SIMPLEX_FACE_FLUX_HIGH_CLOSURE;
    return limited;
  }
  if (fabs(low_sum - 1.) > tolerance) {
    limited.failure_reason = SIMPLEX_FACE_FLUX_LOW_CLOSURE;
    return limited;
  }

  for (int material = 0; material < SIMPLEX_FLUX_MATERIALS; material++) {
    double antidiffusive_flux =
      (high[material] - low[material])*face_velocity;
    double ratio = antidiffusive_flux > 0. ?
      left_depletion_ratio[material] : antidiffusive_flux < 0. ?
      right_depletion_ratio[material] : 1.;
    limited.alpha = fmin(limited.alpha, ratio);
  }

  double sum = 0.;
  for (int material = 0; material < SIMPLEX_FLUX_MATERIALS; material++) {
    limited.fraction[material] = low[material] + limited.alpha*
      (high[material] - low[material]);
    limited.correction_linf = fmax(limited.correction_linf,
      fabs(limited.fraction[material] - high[material]));
    sum += limited.fraction[material];
    if (limited.fraction[material] < -tolerance ||
        limited.fraction[material] > 1. + tolerance) {
      limited.failure_reason = SIMPLEX_FACE_FLUX_LIMITED_OUT_OF_BOUNDS;
      limited.failure_material = material;
      return limited;
    }
  }
  limited.closure_defect = fabs(sum - 1.);
  if (limited.closure_defect > tolerance) {
    limited.failure_reason = SIMPLEX_FACE_FLUX_LIMITED_CLOSURE;
    return limited;
  }
  limited.limited = limited.alpha < 1. - 64.*DBL_EPSILON;
  limited.failure_reason = SIMPLEX_FACE_FLUX_OK;
  limited.valid = 1;
  return limited;
}

#endif
