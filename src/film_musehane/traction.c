#include "film_musehane/traction.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int traction_input_is_valid (
  const MusehaneTractionInput *input, const MusehaneTractionOutput *output)
{
  if (!input || !output || input->cell_count == 0 ||
      input->cell_count > SIZE_MAX/sizeof(double) ||
      !isfinite(input->water_viscosity) || input->water_viscosity <= 0. ||
      !input->active || !input->thickness ||
      !input->pressure_gradient_tangent ||
      !input->lower_tangential_velocity ||
      !input->upper_tangential_velocity ||
      !output->lower_tangential_traction ||
      !output->upper_tangential_traction ||
      output->lower_tangential_traction == output->upper_tangential_traction)
    return 0;
  for (size_t cell = 0; cell < input->cell_count; cell++)
    if (input->active[cell] > 1 ||
        !isfinite(input->thickness[cell]) || input->thickness[cell] <= 0. ||
        !isfinite(input->pressure_gradient_tangent[cell]) ||
        !isfinite(input->lower_tangential_velocity[cell]) ||
        !isfinite(input->upper_tangential_velocity[cell]))
      return 0;
  return 1;
}

MusehaneTractionStatus musehane_evaluate_traction (
  const MusehaneTractionInput *input, MusehaneTractionOutput *output)
{
  if (!traction_input_is_valid(input, output))
    return MUSEHANE_TRACTION_INVALID_INPUT;

  double *lower = malloc(input->cell_count*sizeof *lower);
  double *upper = malloc(input->cell_count*sizeof *upper);
  if (!lower || !upper) {
    free(lower);
    free(upper);
    return MUSEHANE_TRACTION_ALLOCATION_FAILED;
  }

  double maximum_residual = 0.;
  for (size_t cell = 0; cell < input->cell_count; cell++) {
    if (!input->active[cell]) {
      lower[cell] = 0.;
      upper[cell] = 0.;
      continue;
    }
    const double h = input->thickness[cell];
    const double gradient = input->pressure_gradient_tangent[cell];
    const double pressure_part = -.5*h*gradient;
    const double velocity_difference =
      input->lower_tangential_velocity[cell] -
      input->upper_tangential_velocity[cell];
    const double shear_part = input->water_viscosity/h*velocity_difference;
    lower[cell] = pressure_part - shear_part;
    upper[cell] = pressure_part + shear_part;
    if (!isfinite(lower[cell]) || !isfinite(upper[cell])) {
      free(lower);
      free(upper);
      return MUSEHANE_TRACTION_INVALID_INPUT;
    }
    const double pressure_residual =
      lower[cell] + upper[cell] + h*gradient;
    const double shear_residual =
      upper[cell] - lower[cell] - 2.*shear_part;
    maximum_residual = fmax(maximum_residual,
      fmax(fabs(pressure_residual), fabs(shear_residual)));
  }

  double *const output_lower = output->lower_tangential_traction;
  double *const output_upper = output->upper_tangential_traction;
  memcpy(output_lower, lower, input->cell_count*sizeof *lower);
  memcpy(output_upper, upper, input->cell_count*sizeof *upper);
  *output = (MusehaneTractionOutput) {
    .lower_tangential_traction = output_lower,
    .upper_tangential_traction = output_upper,
    .maximum_constitutive_residual = maximum_residual
  };
  free(lower);
  free(upper);
  return MUSEHANE_TRACTION_OK;
}
