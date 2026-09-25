#include "film_model/musehane_traction_deposition.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int finite_ring_input (
  const FilmRingObservation *ring, const FilmRingTraction *traction)
{
  if (!ring || !traction || !ring->geometry_valid ||
      !ring->lower_identity_valid || !ring->upper_identity_valid ||
      !isfinite(ring->radius) || ring->radius < 0. ||
      !isfinite(ring->lower_position) ||
      !isfinite(ring->upper_position) ||
      !isfinite(ring->lower_tangential_velocity) ||
      !isfinite(ring->upper_tangential_velocity) ||
      !isfinite(ring->lower_normal_x) ||
      !isfinite(ring->lower_normal_r) ||
      !isfinite(ring->upper_normal_x) ||
      !isfinite(ring->upper_normal_r) ||
      !isfinite(traction->lower_normal_force) ||
      !isfinite(traction->upper_normal_force) ||
      !isfinite(traction->lower_tangential_force) ||
      !isfinite(traction->upper_tangential_force) ||
      !isfinite(traction->resolved_tangential_reaction_force) ||
      traction->lower_normal_force != 0. ||
      traction->upper_normal_force != 0.)
    return 0;
  const double lower_norm = hypot(ring->lower_normal_x,
    ring->lower_normal_r);
  const double upper_norm = hypot(ring->upper_normal_x,
    ring->upper_normal_r);
  return isfinite(lower_norm) && lower_norm > DBL_MIN &&
    isfinite(upper_norm) && upper_norm > DBL_MIN;
}

MusehaneTractionDepositionStatus musehane_plan_traction_deposition (
  const MusehaneTractionDepositionInput *input,
  MusehaneTractionDepositionOutput *output)
{
  if (!input || !output || input->ring_count == 0 ||
      input->ring_count > SIZE_MAX/(4*sizeof(double)) ||
      !input->rings || !input->tractions || !output->lower_force_x ||
      !output->lower_force_r || !output->upper_force_x ||
      !output->upper_force_r)
    return MUSEHANE_TRACTION_DEPOSITION_INVALID_INPUT;

  double *candidate = calloc(4*input->ring_count, sizeof *candidate);
  if (!candidate)
    return MUSEHANE_TRACTION_DEPOSITION_ALLOCATION_FAILED;
  double *lower_x = candidate;
  double *lower_r = lower_x + input->ring_count;
  double *upper_x = lower_r + input->ring_count;
  double *upper_r = upper_x + input->ring_count;
  MusehaneTractionDepositionOutput result = {
    .lower_force_x = lower_x,
    .lower_force_r = lower_r,
    .upper_force_x = upper_x,
    .upper_force_r = upper_r
  };

  for (size_t ring_index = 0; ring_index < input->ring_count;
       ring_index++) {
    const FilmRingObservation *ring = &input->rings[ring_index];
    const FilmRingTraction *traction = &input->tractions[ring_index];
    if (!finite_ring_input(ring, traction)) {
      free(candidate);
      return MUSEHANE_TRACTION_DEPOSITION_INVALID_INPUT;
    }
    const double pair_residual = traction->lower_tangential_force +
      traction->upper_tangential_force +
      traction->resolved_tangential_reaction_force;
    const double pair_scale = fmax(DBL_MIN,
      fmax(fabs(traction->lower_tangential_force),
        fmax(fabs(traction->upper_tangential_force),
          fabs(traction->resolved_tangential_reaction_force))));
    result.maximum_scalar_pair_residual = fmax(
      result.maximum_scalar_pair_residual, fabs(pair_residual));
    if (fabs(pair_residual) > 512.*DBL_EPSILON*pair_scale) {
      free(candidate);
      return MUSEHANE_TRACTION_DEPOSITION_UNBALANCED_PAIR;
    }

    const double lower_norm = hypot(ring->lower_normal_x,
      ring->lower_normal_r);
    const double upper_norm = hypot(ring->upper_normal_x,
      ring->upper_normal_r);
    const double lower_tangent_x = -ring->lower_normal_r/lower_norm;
    const double lower_tangent_r = ring->lower_normal_x/lower_norm;
    const double upper_tangent_x = -ring->upper_normal_r/upper_norm;
    const double upper_tangent_r = ring->upper_normal_x/upper_norm;
    lower_x[ring_index] = traction->lower_tangential_force*lower_tangent_x;
    lower_r[ring_index] = traction->lower_tangential_force*lower_tangent_r;
    upper_x[ring_index] = traction->upper_tangential_force*upper_tangent_x;
    upper_r[ring_index] = traction->upper_tangential_force*upper_tangent_r;
    result.requested_lower_power += traction->lower_tangential_force*
      ring->lower_tangential_velocity;
    result.requested_upper_power += traction->upper_tangential_force*
      ring->upper_tangential_velocity;
    result.requested_resolved_reaction_force +=
      traction->resolved_tangential_reaction_force;
    result.requested_axial_force += lower_x[ring_index] +
      upper_x[ring_index];
    result.requested_meridional_radial_force += lower_r[ring_index] +
      upper_r[ring_index];
    result.requested_meridional_moment +=
      ring->lower_position*lower_r[ring_index] -
        ring->radius*lower_x[ring_index] +
      ring->upper_position*upper_r[ring_index] -
        ring->radius*upper_x[ring_index];
  }
  result.requested_total_power = result.requested_lower_power +
    result.requested_upper_power;
  if (!isfinite(result.requested_lower_power) ||
      !isfinite(result.requested_upper_power) ||
      !isfinite(result.requested_total_power) ||
      !isfinite(result.requested_resolved_reaction_force) ||
      !isfinite(result.requested_axial_force) ||
      !isfinite(result.requested_meridional_radial_force) ||
      !isfinite(result.requested_meridional_moment)) {
    free(candidate);
    return MUSEHANE_TRACTION_DEPOSITION_INVALID_INPUT;
  }

  double *const output_lower_x = output->lower_force_x;
  double *const output_lower_r = output->lower_force_r;
  double *const output_upper_x = output->upper_force_x;
  double *const output_upper_r = output->upper_force_r;
  memcpy(output_lower_x, lower_x, input->ring_count*sizeof(double));
  memcpy(output_lower_r, lower_r, input->ring_count*sizeof(double));
  memcpy(output_upper_x, upper_x, input->ring_count*sizeof(double));
  memcpy(output_upper_r, upper_r, input->ring_count*sizeof(double));
  result.lower_force_x = output_lower_x;
  result.lower_force_r = output_lower_r;
  result.upper_force_x = output_upper_x;
  result.upper_force_r = output_upper_r;
  *output = result;
  free(candidate);
  return MUSEHANE_TRACTION_DEPOSITION_OK;
}
