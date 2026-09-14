#ifndef AGGREGATE_MASS_WEIGHTING_H
#define AGGREGATE_MASS_WEIGHTING_H

#include <stdbool.h>
#include <math.h>

typedef struct {
  double value;
  bool valid;
  bool zero_measure;
} AggregateMassWeightedValue;

static inline AggregateMassWeightedValue aggregate_mass_weighted_value (
  double weighted_total, double total_mass)
{
  AggregateMassWeightedValue result = {0., false, false};
  if (!isfinite(weighted_total) || !isfinite(total_mass) || total_mass < 0.)
    return result;
  if (total_mass == 0.) {
    result.zero_measure = true;
    result.valid = weighted_total == 0.;
    return result;
  }
  result.value = weighted_total/total_mass;
  result.valid = isfinite(result.value);
  if (!result.valid)
    result.value = 0.;
  return result;
}

#endif
