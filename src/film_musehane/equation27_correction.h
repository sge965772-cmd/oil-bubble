#ifndef FILM_MUSEHANE_EQUATION27_CORRECTION_H
#define FILM_MUSEHANE_EQUATION27_CORRECTION_H

#include <stddef.h>

typedef enum {
  MUSEHANE_EQUATION27_CORRECTION_OK = 0,
  MUSEHANE_EQUATION27_CORRECTION_INVALID_INPUT = 1,
  MUSEHANE_EQUATION27_CORRECTION_ALLOCATION_FAILED = 2
} MusehaneEquation27CorrectionStatus;

typedef struct {
  size_t cell_count;
  double film_viscosity;
  double matching_cells;
  const unsigned char *active;
  const double *thickness;
  const double *local_delta;
  const double *area_weight;
  const double *pressure_gradient_tangent;
  const double *lower_tangential_velocity;
  const double *upper_tangential_velocity;
} MusehaneEquation27CorrectionInput;

typedef struct {
  double *lower_correction_traction;
  double *upper_correction_traction;
  /* Equal and opposite force assigned to the resolved film control volume.
     It closes the nonzero pair sum of the pressure-driven Eq. (27) term. */
  double *resolved_reaction_force;
  size_t corrected_cells;
  double minimum_reference_thickness;
  double maximum_resolved_region_correction;
  double maximum_pressure_correction_traction;
  double maximum_shear_correction_traction;
  double total_resolved_reaction_force;
  double maximum_pair_closure_residual;
  double pressure_correction_power;
  double shear_correction_power;
  double total_correction_power;
  double direct_correction_power;
  double power_identity_residual;
} MusehaneEquation27CorrectionOutput;

/* Full equation-(27) traction minus the same constitutive traction at
   h_ref=max(h, matching_cells*Delta). Pressure-transfer work is signed;
   the shear contribution must remain dissipative. No clipping is applied. */
MusehaneEquation27CorrectionStatus musehane_evaluate_equation27_correction (
  const MusehaneEquation27CorrectionInput *input,
  MusehaneEquation27CorrectionOutput *output);

/* Validate the exact pressure-plus-shear work identity at a scale which stays
   well-conditioned when the two constituent powers nearly cancel. */
int musehane_equation27_power_audit_is_valid (
  const MusehaneEquation27CorrectionOutput *output);

#endif
