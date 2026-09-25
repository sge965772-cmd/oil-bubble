#include "film_musehane/equation27_correction.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int equation27_correction_input_is_valid (
  const MusehaneEquation27CorrectionInput *input,
  const MusehaneEquation27CorrectionOutput *output)
{
  if (!input || !output || input->cell_count == 0 ||
      input->cell_count > SIZE_MAX/(3*sizeof(double)) ||
      !isfinite(input->film_viscosity) || input->film_viscosity <= 0. ||
      !isfinite(input->matching_cells) || input->matching_cells <= 0. ||
      !input->active || !input->thickness || !input->local_delta ||
      !input->area_weight || !input->pressure_gradient_tangent ||
      !input->lower_tangential_velocity ||
      !input->upper_tangential_velocity ||
      !output->lower_correction_traction ||
      !output->upper_correction_traction ||
      !output->resolved_reaction_force ||
      output->lower_correction_traction ==
        output->upper_correction_traction ||
      output->lower_correction_traction == output->resolved_reaction_force ||
      output->upper_correction_traction == output->resolved_reaction_force)
    return 0;
  for (size_t cell = 0; cell < input->cell_count; cell++)
    if (input->active[cell] > 1 ||
        !isfinite(input->thickness[cell]) || input->thickness[cell] <= 0. ||
        !isfinite(input->local_delta[cell]) || input->local_delta[cell] <= 0. ||
        !isfinite(input->area_weight[cell]) || input->area_weight[cell] <= 0. ||
        !isfinite(input->pressure_gradient_tangent[cell]) ||
        !isfinite(input->lower_tangential_velocity[cell]) ||
        !isfinite(input->upper_tangential_velocity[cell]))
      return 0;
  return 1;
}

int musehane_equation27_power_audit_is_valid (
  const MusehaneEquation27CorrectionOutput *output)
{
  if (!output ||
      !isfinite(output->pressure_correction_power) ||
      !isfinite(output->shear_correction_power) ||
      !isfinite(output->direct_correction_power) ||
      !isfinite(output->total_correction_power) ||
      output->shear_correction_power > 0.)
    return 0;
  const double scale = fmax(DBL_MIN,
    fmax(fabs(output->pressure_correction_power),
      fmax(fabs(output->shear_correction_power),
        fmax(fabs(output->direct_correction_power),
          fabs(output->total_correction_power)))));
  return fabs(output->direct_correction_power -
    output->total_correction_power) <= 512.*DBL_EPSILON*scale;
}

MusehaneEquation27CorrectionStatus musehane_evaluate_equation27_correction (
  const MusehaneEquation27CorrectionInput *input,
  MusehaneEquation27CorrectionOutput *output)
{
  if (!equation27_correction_input_is_valid(input, output))
    return MUSEHANE_EQUATION27_CORRECTION_INVALID_INPUT;

  double *storage = calloc(3*input->cell_count, sizeof *storage);
  if (!storage)
    return MUSEHANE_EQUATION27_CORRECTION_ALLOCATION_FAILED;
  double *lower = storage;
  double *upper = lower + input->cell_count;
  double *reaction = upper + input->cell_count;
  MusehaneEquation27CorrectionOutput candidate = {
    .lower_correction_traction = lower,
    .upper_correction_traction = upper,
    .resolved_reaction_force = reaction,
    .minimum_reference_thickness = INFINITY
  };

  for (size_t cell = 0; cell < input->cell_count; cell++) {
    if (!input->active[cell])
      continue;
    const double h = input->thickness[cell];
    const double h_reference = fmax(h,
      input->matching_cells*input->local_delta[cell]);
    const double gradient = input->pressure_gradient_tangent[cell];
    const double lower_velocity = input->lower_tangential_velocity[cell];
    const double upper_velocity = input->upper_tangential_velocity[cell];
    const double relative_velocity = lower_velocity - upper_velocity;
    const double pressure_correction = -.5*(h - h_reference)*gradient;
    const double shear_correction = input->film_viscosity*
      (1./h - 1./h_reference)*relative_velocity;
    lower[cell] = pressure_correction - shear_correction;
    upper[cell] = pressure_correction + shear_correction;
    const double area = input->area_weight[cell];
    const double lower_force = area*lower[cell];
    const double upper_force = area*upper[cell];
    reaction[cell] = -(lower_force + upper_force);
    const double closure = lower_force + upper_force + reaction[cell];

    candidate.minimum_reference_thickness = fmin(
      candidate.minimum_reference_thickness, h_reference);
    if (h_reference > h)
      candidate.corrected_cells++;
    else
      candidate.maximum_resolved_region_correction = fmax(
        candidate.maximum_resolved_region_correction,
        fmax(fabs(lower[cell]), fabs(upper[cell])));
    candidate.maximum_pressure_correction_traction = fmax(
      candidate.maximum_pressure_correction_traction,
      fabs(pressure_correction));
    candidate.maximum_shear_correction_traction = fmax(
      candidate.maximum_shear_correction_traction,
      fabs(shear_correction));
    candidate.total_resolved_reaction_force += reaction[cell];
    candidate.maximum_pair_closure_residual = fmax(
      candidate.maximum_pair_closure_residual, fabs(closure));
    candidate.pressure_correction_power +=
      area*pressure_correction*(lower_velocity + upper_velocity);
    candidate.shear_correction_power -=
      area*shear_correction*relative_velocity;
    candidate.direct_correction_power += area*(
      lower[cell]*lower_velocity + upper[cell]*upper_velocity);
  }
  candidate.total_correction_power =
    candidate.pressure_correction_power + candidate.shear_correction_power;
  candidate.power_identity_residual = candidate.direct_correction_power -
    candidate.total_correction_power;
  if (!isfinite(candidate.minimum_reference_thickness))
    candidate.minimum_reference_thickness = 0.;
  if (!isfinite(candidate.maximum_resolved_region_correction) ||
      !isfinite(candidate.maximum_pressure_correction_traction) ||
      !isfinite(candidate.maximum_shear_correction_traction) ||
      !isfinite(candidate.total_resolved_reaction_force) ||
      !isfinite(candidate.maximum_pair_closure_residual) ||
      !isfinite(candidate.pressure_correction_power) ||
      !isfinite(candidate.shear_correction_power) ||
      !isfinite(candidate.total_correction_power) ||
      !isfinite(candidate.direct_correction_power) ||
      !isfinite(candidate.power_identity_residual) ||
      candidate.shear_correction_power > 0.) {
    free(storage);
    return MUSEHANE_EQUATION27_CORRECTION_INVALID_INPUT;
  }

  double *const output_lower = output->lower_correction_traction;
  double *const output_upper = output->upper_correction_traction;
  double *const output_reaction = output->resolved_reaction_force;
  memcpy(output_lower, lower, input->cell_count*sizeof *lower);
  memcpy(output_upper, upper, input->cell_count*sizeof *upper);
  memcpy(output_reaction, reaction, input->cell_count*sizeof *reaction);
  candidate.lower_correction_traction = output_lower;
  candidate.upper_correction_traction = output_upper;
  candidate.resolved_reaction_force = output_reaction;
  *output = candidate;
  free(storage);
  return MUSEHANE_EQUATION27_CORRECTION_OK;
}
