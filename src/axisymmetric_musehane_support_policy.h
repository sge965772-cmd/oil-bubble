#ifndef AXISYMMETRIC_MUSEHANE_SUPPORT_POLICY_H
#define AXISYMMETRIC_MUSEHANE_SUPPORT_POLICY_H

#include <math.h>
#include <stdbool.h>

typedef struct {
  bool valid;
  bool geometric_fallback;
  double support;
} AxisymmetricMusehaneSupportSelection;

/* Preserve the phase-weighted support whenever it exists.  A required
   component may use the same compact geometric kernel without identity
   weighting only when discrete side assignment leaves it with no support.
   Observation and deposition must use this selection symmetrically. */
static inline AxisymmetricMusehaneSupportSelection
axisymmetric_musehane_select_component_support (
  double phase_weighted_support,
  double geometric_support,
  bool required)
{
  AxisymmetricMusehaneSupportSelection result = {0};
  if (!isfinite(phase_weighted_support) || phase_weighted_support < 0. ||
      !isfinite(geometric_support) || geometric_support < 0.)
    return result;
  if (!required || phase_weighted_support > 0.) {
    result.valid = true;
    result.support = phase_weighted_support;
    return result;
  }
  if (geometric_support > 0.) {
    result.valid = true;
    result.geometric_fallback = true;
    result.support = geometric_support;
  }
  return result;
}

#endif
