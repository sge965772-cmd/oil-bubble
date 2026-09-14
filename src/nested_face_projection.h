#ifndef NESTED_FACE_PROJECTION_H
#define NESTED_FACE_PROJECTION_H

#include <float.h>
#include <math.h>

#include "positive_material_partition.h"

typedef struct {
  double raw_lower_gas;
  double raw_lower_envelope;
  double raw_upper_gas;
  double raw_upper_envelope;
  double vof_lower_gas;
  double vof_lower_envelope;
  double vof_upper_gas;
  double vof_upper_envelope;
  double lower_gas;
  double lower_envelope;
  double upper_gas;
  double upper_envelope;
  PositiveMaterialPartition materials;
  double raw_bound_defect;
  double raw_identity_defect;
  double correction_linf;
  int projected;
  int valid;
} NestedFaceProjection;

/*
 * Project the five raw identity-material fractions onto the probability
 * simplex. Reconstructing the four indicators from the projected materials
 * enforces gas <= envelope and disjoint envelopes on every transported face.
 */
static inline NestedFaceProjection nested_face_projection (
  double lower_gas, double lower_envelope,
  double upper_gas, double upper_envelope)
{
  NestedFaceProjection projection = {0};
  projection.raw_lower_gas = lower_gas;
  projection.raw_lower_envelope = lower_envelope;
  projection.raw_upper_gas = upper_gas;
  projection.raw_upper_envelope = upper_envelope;
  double indicators[4] = {
    lower_gas, lower_envelope, upper_gas, upper_envelope
  };
  for (int index = 0; index < 4; index++)
    if (!isfinite(indicators[index])) {
      projection.raw_bound_defect = INFINITY;
      projection.raw_identity_defect = INFINITY;
      return projection;
    }

  projection.vof_lower_gas = fmin(fmax(lower_gas, 0.), 1.);
  projection.vof_lower_envelope = fmin(fmax(lower_envelope, 0.), 1.);
  projection.vof_upper_gas = fmin(fmax(upper_gas, 0.), 1.);
  projection.vof_upper_envelope = fmin(fmax(upper_envelope, 0.), 1.);

  double raw[POSITIVE_MATERIAL_COUNT] = {
    1. - projection.vof_lower_envelope - projection.vof_upper_envelope,
    projection.vof_lower_envelope - projection.vof_lower_gas,
    projection.vof_lower_gas,
    projection.vof_upper_envelope - projection.vof_upper_gas,
    projection.vof_upper_gas
  };
  for (int index = 0; index < 4; index++) {
    projection.raw_bound_defect = fmax(
      projection.raw_bound_defect, -indicators[index]);
    projection.raw_bound_defect = fmax(
      projection.raw_bound_defect, indicators[index] - 1.);
  }
  for (int material = 0; material < POSITIVE_MATERIAL_COUNT; material++)
    projection.raw_identity_defect = fmax(
      projection.raw_identity_defect, -raw[material]);

  double sorted[POSITIVE_MATERIAL_COUNT];
  for (int material = 0; material < POSITIVE_MATERIAL_COUNT; material++)
    sorted[material] = raw[material];
  for (int left = 0; left < POSITIVE_MATERIAL_COUNT - 1; left++)
    for (int right = left + 1; right < POSITIVE_MATERIAL_COUNT; right++)
      if (sorted[right] > sorted[left]) {
        double swap = sorted[left];
        sorted[left] = sorted[right];
        sorted[right] = swap;
      }

  int active = -1;
  double cumulative = 0., theta = 0.;
  for (int material = 0; material < POSITIVE_MATERIAL_COUNT; material++) {
    cumulative += sorted[material];
    double candidate = (cumulative - 1.)/(material + 1.);
    if (sorted[material] - candidate > 0.) {
      active = material;
      theta = candidate;
    }
  }
  if (active < 0 || !isfinite(theta)) {
    projection.raw_bound_defect = INFINITY;
    projection.raw_identity_defect = INFINITY;
    return projection;
  }

  double corrected[POSITIVE_MATERIAL_COUNT], corrected_sum = 0.;
  int largest = 0;
  for (int material = 0; material < POSITIVE_MATERIAL_COUNT; material++) {
    corrected[material] = fmax(raw[material] - theta, 0.);
    corrected_sum += corrected[material];
    if (corrected[material] > corrected[largest])
      largest = material;
  }
  corrected[largest] += 1. - corrected_sum;
  if (corrected[largest] < 0. || !isfinite(corrected[largest]))
    return projection;

  double correction_l1 = 0.;
  for (int material = 0; material < POSITIVE_MATERIAL_COUNT; material++) {
    double correction = fabs(corrected[material] - raw[material]);
    correction_l1 += correction;
    projection.correction_linf = fmax(
      projection.correction_linf, correction);
  }

  projection.materials.water = corrected[POSITIVE_MATERIAL_WATER];
  projection.materials.lower_oil =
    corrected[POSITIVE_MATERIAL_LOWER_OIL];
  projection.materials.lower_gas =
    corrected[POSITIVE_MATERIAL_LOWER_GAS];
  projection.materials.upper_oil =
    corrected[POSITIVE_MATERIAL_UPPER_OIL];
  projection.materials.upper_gas =
    corrected[POSITIVE_MATERIAL_UPPER_GAS];
  projection.materials.bound_defect = projection.raw_bound_defect;
  projection.materials.identity_defect = projection.raw_identity_defect;
  projection.materials.water_defect = fmax(0., -raw[0]);
  projection.materials.lower_oil_defect = fmax(0., -raw[1]);
  projection.materials.upper_oil_defect = fmax(0., -raw[3]);
  projection.materials.defect = fmax(
    projection.raw_bound_defect, projection.raw_identity_defect);
  projection.materials.correction_l1 = correction_l1;
  projection.materials.valid = 1;

  projection.lower_gas = projection.materials.lower_gas;
  projection.lower_envelope = projection.materials.lower_oil +
                              projection.materials.lower_gas;
  projection.upper_gas = projection.materials.upper_gas;
  projection.upper_envelope = projection.materials.upper_oil +
                              projection.materials.upper_gas;
  projection.projected = projection.correction_linf > 64.*DBL_EPSILON;
  projection.valid = 1;
  return projection;
}

#endif
