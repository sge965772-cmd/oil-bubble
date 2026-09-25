#ifndef FILM_MUSEHANE_MATCHED_TRACTION_H
#define FILM_MUSEHANE_MATCHED_TRACTION_H

#include <stddef.h>

typedef enum {
  MUSEHANE_MATCHED_TRACTION_OK = 0,
  MUSEHANE_MATCHED_TRACTION_INVALID_INPUT = 1,
  MUSEHANE_MATCHED_TRACTION_ALLOCATION_FAILED = 2
} MusehaneMatchedTractionStatus;

typedef struct {
  size_t cell_count;
  double water_viscosity;
  double matching_cells;
  const unsigned char *active;
  const double *thickness;
  const double *local_delta;
  /* Physical axisymmetric annular area, including 2*pi. */
  const double *area_weight;
  const double *pressure_gradient_tangent;
  const double *lower_tangential_velocity;
  const double *upper_tangential_velocity;
} MusehaneMatchedTractionInput;

typedef struct {
  double *lower_correction_traction;
  double *upper_correction_traction;
  size_t corrected_cells;
  double minimum_reference_thickness;
  double maximum_resolved_region_correction;
  double analytical_film_power;
  double resolved_reference_power;
  /* Stable constitutive power uses the squared relative velocity. The direct
     force-dot-velocity sum is retained only to audit cancellation. */
  double correction_power;
  double direct_correction_power;
  double power_identity_residual;
  double analytical_relative_shear_dissipation;
  double resolved_relative_shear_dissipation;
  double additional_relative_shear_dissipation;
} MusehaneMatchedTractionOutput;

/* Viscous part of equation-(27) minus the same viscous traction evaluated at
 * h_ref=max(h, matching_cells*Delta). The pressure-gradient contribution is
 * diagnostic only because DNS already resolves pressure. This is a stress
 * replacement seam, not an added collision force: the correction is exactly
 * zero in the resolved region and no passivity clipping is applied. */
MusehaneMatchedTractionStatus musehane_evaluate_matched_traction (
  const MusehaneMatchedTractionInput *input,
  MusehaneMatchedTractionOutput *output);

#endif
