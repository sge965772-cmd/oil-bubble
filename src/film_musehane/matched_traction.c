#include "film_musehane/matched_traction.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int matched_traction_input_is_valid (
  const MusehaneMatchedTractionInput *input,
  const MusehaneMatchedTractionOutput *output)
{
  if (!input || !output || input->cell_count == 0 ||
      input->cell_count > SIZE_MAX/sizeof(double) ||
      !isfinite(input->water_viscosity) || input->water_viscosity <= 0. ||
      !isfinite(input->matching_cells) || input->matching_cells <= 0. ||
      !input->active || !input->thickness || !input->local_delta ||
      !input->area_weight || !input->pressure_gradient_tangent ||
      !input->lower_tangential_velocity ||
      !input->upper_tangential_velocity ||
      !output->lower_correction_traction ||
      !output->upper_correction_traction ||
      output->lower_correction_traction ==
        output->upper_correction_traction)
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

MusehaneMatchedTractionStatus musehane_evaluate_matched_traction (
  const MusehaneMatchedTractionInput *input,
  MusehaneMatchedTractionOutput *output)
{
  if (!matched_traction_input_is_valid(input, output))
    return MUSEHANE_MATCHED_TRACTION_INVALID_INPUT;

  double *lower = calloc(input->cell_count, sizeof *lower);
  double *upper = calloc(input->cell_count, sizeof *upper);
  if (!lower || !upper) {
    free(lower);
    free(upper);
    return MUSEHANE_MATCHED_TRACTION_ALLOCATION_FAILED;
  }

  MusehaneMatchedTractionOutput candidate = {
    .lower_correction_traction = lower,
    .upper_correction_traction = upper,
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
    const double velocity_difference = lower_velocity - upper_velocity;
    const double film_pressure_part = -.5*h*gradient;
    const double reference_pressure_part = -.5*h_reference*gradient;
    const double film_shear_part =
      input->water_viscosity/h*velocity_difference;
    const double reference_shear_part =
      input->water_viscosity/h_reference*velocity_difference;
    const double film_lower = film_pressure_part - film_shear_part;
    const double film_upper = film_pressure_part + film_shear_part;
    const double reference_lower =
      reference_pressure_part - reference_shear_part;
    const double reference_upper =
      reference_pressure_part + reference_shear_part;
    lower[cell] = -(film_shear_part - reference_shear_part);
    upper[cell] = film_shear_part - reference_shear_part;
    if (!isfinite(h_reference) || !isfinite(lower[cell]) ||
        !isfinite(upper[cell])) {
      free(lower);
      free(upper);
      return MUSEHANE_MATCHED_TRACTION_INVALID_INPUT;
    }
    candidate.minimum_reference_thickness = fmin(
      candidate.minimum_reference_thickness, h_reference);
    if (h_reference > h)
      candidate.corrected_cells++;
    else
      candidate.maximum_resolved_region_correction = fmax(
        candidate.maximum_resolved_region_correction,
        fmax(fabs(lower[cell]), fabs(upper[cell])));

    const double area = input->area_weight[cell];
    const double film_power = area*(film_lower*lower_velocity +
      film_upper*upper_velocity);
    const double reference_power = area*(reference_lower*lower_velocity +
      reference_upper*upper_velocity);
    const double film_dissipation = area*input->water_viscosity/h*
      velocity_difference*velocity_difference;
    const double reference_dissipation = area*input->water_viscosity/
      h_reference*velocity_difference*velocity_difference;
    candidate.analytical_film_power += film_power;
    candidate.resolved_reference_power += reference_power;
    candidate.direct_correction_power += area*(
      lower[cell]*lower_velocity + upper[cell]*upper_velocity);
    candidate.analytical_relative_shear_dissipation += film_dissipation;
    candidate.resolved_relative_shear_dissipation += reference_dissipation;
    candidate.additional_relative_shear_dissipation +=
      film_dissipation - reference_dissipation;
  }
  candidate.correction_power =
    -candidate.additional_relative_shear_dissipation;
  candidate.power_identity_residual = candidate.direct_correction_power -
    candidate.correction_power;
  if (!isfinite(candidate.minimum_reference_thickness))
    candidate.minimum_reference_thickness = 0.;
  if (!isfinite(candidate.analytical_film_power) ||
      !isfinite(candidate.resolved_reference_power) ||
      !isfinite(candidate.correction_power) ||
      !isfinite(candidate.direct_correction_power) ||
      !isfinite(candidate.power_identity_residual) ||
      !isfinite(candidate.analytical_relative_shear_dissipation) ||
      !isfinite(candidate.resolved_relative_shear_dissipation) ||
      !isfinite(candidate.additional_relative_shear_dissipation) ||
      candidate.additional_relative_shear_dissipation < 0.) {
    free(lower);
    free(upper);
    return MUSEHANE_MATCHED_TRACTION_INVALID_INPUT;
  }

  double *const output_lower = output->lower_correction_traction;
  double *const output_upper = output->upper_correction_traction;
  memcpy(output_lower, lower, input->cell_count*sizeof *lower);
  memcpy(output_upper, upper, input->cell_count*sizeof *upper);
  candidate.lower_correction_traction = output_lower;
  candidate.upper_correction_traction = output_upper;
  *output = candidate;
  free(lower);
  free(upper);
  return MUSEHANE_MATCHED_TRACTION_OK;
}
