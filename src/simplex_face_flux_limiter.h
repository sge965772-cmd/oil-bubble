#ifndef SIMPLEX_FACE_FLUX_LIMITER_H
#define SIMPLEX_FACE_FLUX_LIMITER_H

#include <float.h>
#include <math.h>

#define SIMPLEX_FLUX_MATERIALS 5

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
