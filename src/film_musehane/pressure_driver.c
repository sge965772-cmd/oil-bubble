#include "film_musehane/pressure_driver.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int pressure_driver_input_is_valid (
  const MusehanePressureDriverConfig *config,
  const MusehanePressureDriverInput *input,
  const MusehanePressureDriverOutput *output)
{
  if (!config || !input || !output || !input->mesh ||
      !isfinite(config->smoothing_length_squared) ||
      config->smoothing_length_squared < 0. ||
      !isfinite(config->linear_residual_relative_tolerance) ||
      config->linear_residual_relative_tolerance <= 0. ||
      input->mesh->cell_count == 0 ||
      input->mesh->cell_count > SIZE_MAX/(5*sizeof(double)) ||
      !input->mesh->cell_radius || !input->mesh->cell_width ||
      !input->mesh->cell_measure || !input->mesh->face_measure ||
      !input->raw_pressure_gradient_tangent_cell ||
      !input->capillary_force_density_tangent_cell ||
      !output->decontaminated_gradient_cell ||
      !output->smoothed_gradient_cell ||
      !output->smoothed_gradient_face ||
      output->decontaminated_gradient_cell == output->smoothed_gradient_cell ||
      output->decontaminated_gradient_cell == output->smoothed_gradient_face ||
      output->smoothed_gradient_cell == output->smoothed_gradient_face)
    return 0;

  const size_t cells = input->mesh->cell_count;
  for (size_t face = 0; face <= cells; face++)
    if (!isfinite(input->mesh->face_measure[face]) ||
        input->mesh->face_measure[face] < 0.)
      return 0;
  for (size_t cell = 0; cell < cells; cell++)
    if (!isfinite(input->mesh->cell_radius[cell]) ||
        input->mesh->cell_radius[cell] < 0. ||
        !isfinite(input->mesh->cell_width[cell]) ||
        input->mesh->cell_width[cell] <= 0. ||
        !isfinite(input->mesh->cell_measure[cell]) ||
        input->mesh->cell_measure[cell] <= 0. ||
        !isfinite(input->raw_pressure_gradient_tangent_cell[cell]) ||
        !isfinite(input->capillary_force_density_tangent_cell[cell]) ||
        (cell > 0 && input->mesh->cell_radius[cell] <=
          input->mesh->cell_radius[cell - 1]))
      return 0;
  return 1;
}

static double left_coefficient (
  const MusehanePressureDriverConfig *config,
  const MusehaneSurfaceMesh *mesh, size_t cell)
{
  if (cell == 0)
    return 0.;
  const double distance = mesh->cell_radius[cell] -
    mesh->cell_radius[cell - 1];
  return config->smoothing_length_squared*mesh->face_measure[cell]/distance;
}

static double right_coefficient (
  const MusehanePressureDriverConfig *config,
  const MusehaneSurfaceMesh *mesh, size_t cell)
{
  if (cell + 1 == mesh->cell_count)
    return 0.;
  const double distance = mesh->cell_radius[cell + 1] -
    mesh->cell_radius[cell];
  return config->smoothing_length_squared*mesh->face_measure[cell + 1]/
    distance;
}

MusehanePressureDriverStatus musehane_prepare_pressure_driver (
  const MusehanePressureDriverConfig *config,
  const MusehanePressureDriverInput *input,
  MusehanePressureDriverOutput *output)
{
  if (!pressure_driver_input_is_valid(config, input, output))
    return MUSEHANE_PRESSURE_DRIVER_INVALID_INPUT;

  const size_t cells = input->mesh->cell_count;
  double *storage = malloc(5*cells*sizeof *storage);
  double *face_gradient = malloc((cells + 1)*sizeof *face_gradient);
  if (!storage || !face_gradient) {
    free(storage);
    free(face_gradient);
    return MUSEHANE_PRESSURE_DRIVER_ALLOCATION_FAILED;
  }
  double *corrected = storage;
  double *lower = corrected + cells;
  double *diagonal = lower + cells;
  double *upper = diagonal + cells;
  double *rhs = upper + cells;

  for (size_t cell = 0; cell < cells; cell++) {
    corrected[cell] = input->raw_pressure_gradient_tangent_cell[cell] -
      input->capillary_force_density_tangent_cell[cell];
    if (!isfinite(corrected[cell])) {
      free(storage);
      free(face_gradient);
      return MUSEHANE_PRESSURE_DRIVER_INVALID_INPUT;
    }
    const double left = left_coefficient(config, input->mesh, cell);
    const double right = right_coefficient(config, input->mesh, cell);
    lower[cell] = -left;
    diagonal[cell] = input->mesh->cell_measure[cell] + left + right;
    upper[cell] = -right;
    rhs[cell] = input->mesh->cell_measure[cell]*corrected[cell];
    if (!isfinite(diagonal[cell]) || diagonal[cell] <= 0. ||
        !isfinite(rhs[cell])) {
      free(storage);
      free(face_gradient);
      return MUSEHANE_PRESSURE_DRIVER_INVALID_INPUT;
    }
  }

  for (size_t cell = 1; cell < cells; cell++) {
    if (!isfinite(diagonal[cell - 1]) ||
        fabs(diagonal[cell - 1]) <= DBL_MIN) {
      free(storage);
      free(face_gradient);
      return MUSEHANE_PRESSURE_DRIVER_SINGULAR;
    }
    const double factor = lower[cell]/diagonal[cell - 1];
    diagonal[cell] -= factor*upper[cell - 1];
    rhs[cell] -= factor*rhs[cell - 1];
  }
  if (!isfinite(diagonal[cells - 1]) ||
      fabs(diagonal[cells - 1]) <= DBL_MIN) {
    free(storage);
    free(face_gradient);
    return MUSEHANE_PRESSURE_DRIVER_SINGULAR;
  }
  rhs[cells - 1] /= diagonal[cells - 1];
  for (size_t cell = cells - 1; cell > 0; cell--) {
    const size_t previous = cell - 1;
    if (!isfinite(diagonal[previous]) ||
        fabs(diagonal[previous]) <= DBL_MIN) {
      free(storage);
      free(face_gradient);
      return MUSEHANE_PRESSURE_DRIVER_SINGULAR;
    }
    rhs[previous] = (rhs[previous] - upper[previous]*rhs[cell])/
      diagonal[previous];
  }

  double maximum_relative_residual = 0.;
  for (size_t cell = 0; cell < cells; cell++) {
    if (!isfinite(rhs[cell])) {
      free(storage);
      free(face_gradient);
      return MUSEHANE_PRESSURE_DRIVER_SINGULAR;
    }
    const double left = left_coefficient(config, input->mesh, cell);
    const double right = right_coefficient(config, input->mesh, cell);
    const double diagonal_term =
      (input->mesh->cell_measure[cell] + left + right)*rhs[cell];
    const double left_term = cell > 0 ? left*rhs[cell - 1] : 0.;
    const double right_term = cell + 1 < cells ? right*rhs[cell + 1] : 0.;
    const double lhs = diagonal_term - left_term - right_term;
    const double target = input->mesh->cell_measure[cell]*corrected[cell];
    const double scale = fmax(DBL_MIN,
      fabs(diagonal_term) + fabs(left_term) + fabs(right_term) +
        fabs(target));
    maximum_relative_residual = fmax(maximum_relative_residual,
      fabs(lhs - target)/scale);
  }
  if (!isfinite(maximum_relative_residual) ||
      maximum_relative_residual >
        config->linear_residual_relative_tolerance) {
    free(storage);
    free(face_gradient);
    return MUSEHANE_PRESSURE_DRIVER_RESIDUAL_MISMATCH;
  }

  face_gradient[0] = rhs[0];
  for (size_t face = 1; face < cells; face++) {
    const double left_radius = input->mesh->cell_radius[face - 1];
    const double right_radius = input->mesh->cell_radius[face];
    const double face_radius = left_radius +
      .5*input->mesh->cell_width[face - 1];
    const double fraction = (face_radius - left_radius)/
      (right_radius - left_radius);
    face_gradient[face] = (1. - fraction)*rhs[face - 1] +
      fraction*rhs[face];
  }
  face_gradient[cells] = rhs[cells - 1];

  double *const output_corrected = output->decontaminated_gradient_cell;
  double *const output_smoothed = output->smoothed_gradient_cell;
  double *const output_smoothed_face = output->smoothed_gradient_face;
  memcpy(output_corrected, corrected, cells*sizeof *corrected);
  memcpy(output_smoothed, rhs, cells*sizeof *rhs);
  memcpy(output_smoothed_face, face_gradient,
    (cells + 1)*sizeof *face_gradient);
  *output = (MusehanePressureDriverOutput) {
    .decontaminated_gradient_cell = output_corrected,
    .smoothed_gradient_cell = output_smoothed,
    .smoothed_gradient_face = output_smoothed_face,
    .maximum_linear_relative_residual = maximum_relative_residual
  };
  free(storage);
  free(face_gradient);
  return MUSEHANE_PRESSURE_DRIVER_OK;
}
