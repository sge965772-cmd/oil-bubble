#include "film_musehane/published_equation23_operator.h"

#include <math.h>

typedef struct {
  double value;
  double left_derivative;
  double right_derivative;
} FaceThickness;

static int arrays_are_valid (
  const MusehanePublishedEquation23Input *input,
  const MusehanePublishedEquation23Output *output)
{
  if (!input || !output || !input->mesh || !input->active ||
      !input->thickness || !input->pressure_gradient_cell ||
      !input->pressure_gradient_face ||
      !output->pressure_inventory_rate ||
      !output->pressure_flux_integral_face || !output->lower_derivative ||
      !output->diagonal_derivative || !output->upper_derivative ||
      !isfinite(input->viscosity) || input->viscosity <= 0. ||
      input->mesh->cell_count == 0 || !input->mesh->cell_radius ||
      !input->mesh->cell_width || !input->mesh->cell_measure ||
      !input->mesh->face_measure)
    return 0;
  const size_t cells = input->mesh->cell_count;
  for (size_t cell = 0; cell < cells; cell++)
    if (!isfinite(input->mesh->cell_radius[cell]) ||
        !isfinite(input->mesh->cell_width[cell]) ||
        input->mesh->cell_width[cell] <= 0. ||
        !isfinite(input->mesh->cell_measure[cell]) ||
        input->mesh->cell_measure[cell] <= 0. ||
        !isfinite(input->thickness[cell]) ||
        input->thickness[cell] <= 0. ||
        !isfinite(input->pressure_gradient_cell[cell]))
      return 0;
  for (size_t face = 0; face <= cells; face++)
    if (!isfinite(input->mesh->face_measure[face]) ||
        input->mesh->face_measure[face] < 0. ||
        !isfinite(input->pressure_gradient_face[face]))
      return 0;
  return 1;
}

static int face_thickness (
  const MusehanePublishedEquation23Input *input, size_t face,
  FaceThickness *result)
{
  const size_t cells = input->mesh->cell_count;
  *result = (FaceThickness) {0};
  if (face == 0) {
    if (input->active[0]) {
      result->value = input->thickness[0];
      result->right_derivative = 1.;
    }
    return 1;
  }
  if (face == cells) {
    if (input->active[cells - 1]) {
      result->value = input->thickness[cells - 1];
      result->left_derivative = 1.;
    }
    return 1;
  }
  const size_t left = face - 1, right = face;
  const int left_active = input->active[left] != 0;
  const int right_active = input->active[right] != 0;
  if (left_active && right_active) {
    const double left_radius = input->mesh->cell_radius[left];
    const double right_radius = input->mesh->cell_radius[right];
    const double face_radius = input->mesh->face_measure[face];
    const double span = right_radius - left_radius;
    if (!isfinite(span) || span <= 0. || face_radius < left_radius ||
        face_radius > right_radius)
      return 0;
    result->right_derivative = (face_radius - left_radius)/span;
    result->left_derivative = 1. - result->right_derivative;
    result->value = result->left_derivative*input->thickness[left] +
      result->right_derivative*input->thickness[right];
  }
  else if (left_active) {
    result->value = input->thickness[left];
    result->left_derivative = 1.;
  }
  else if (right_active) {
    result->value = input->thickness[right];
    result->right_derivative = 1.;
  }
  return isfinite(result->value);
}

MusehaneStatus musehane_published_equation23_pressure (
  const MusehanePublishedEquation23Input *input,
  MusehanePublishedEquation23Output *output)
{
  if (!arrays_are_valid(input, output))
    return MUSEHANE_INVALID_INPUT;
  const size_t cells = input->mesh->cell_count;
  for (size_t face = 0; face <= cells; face++) {
    FaceThickness value;
    if (!face_thickness(input, face, &value))
      return MUSEHANE_INVALID_INPUT;
    output->pressure_flux_integral_face[face] =
      input->mesh->face_measure[face]*value.value*value.value*value.value/
      (12.*input->viscosity)*input->pressure_gradient_face[face];
    if (!isfinite(output->pressure_flux_integral_face[face]))
      return MUSEHANE_INVALID_INPUT;
  }
  for (size_t cell = 0; cell < cells; cell++) {
    output->pressure_inventory_rate[cell] = 0.;
    output->lower_derivative[cell] = 0.;
    output->diagonal_derivative[cell] = 0.;
    output->upper_derivative[cell] = 0.;
    if (!input->active[cell])
      continue;
    FaceThickness inward, outward;
    if (!face_thickness(input, cell, &inward) ||
        !face_thickness(input, cell + 1, &outward))
      return MUSEHANE_INVALID_INPUT;
    const double inward_area = input->mesh->face_measure[cell];
    const double outward_area = input->mesh->face_measure[cell + 1];
    const double pressure_laplacian_integral =
      outward_area*input->pressure_gradient_face[cell + 1] -
      inward_area*input->pressure_gradient_face[cell];
    const double thickness_gradient_integral =
      outward_area*outward.value - inward_area*inward.value;
    const double h = input->thickness[cell];
    const double cell_gradient = input->pressure_gradient_cell[cell];
    const double inverse_mu_12 = 1./(12.*input->viscosity);
    const double inverse_mu_4 = 1./(4.*input->viscosity);
    output->pressure_inventory_rate[cell] =
      h*h*h*inverse_mu_12*pressure_laplacian_integral +
      cell_gradient*h*h*inverse_mu_4*thickness_gradient_integral;

    const double inward_self = inward.right_derivative;
    const double outward_self = outward.left_derivative;
    const double thickness_gradient_self =
      outward_area*outward_self - inward_area*inward_self;
    output->diagonal_derivative[cell] =
      3.*h*h*inverse_mu_12*pressure_laplacian_integral +
      cell_gradient*inverse_mu_4*(
        2.*h*thickness_gradient_integral +
        h*h*thickness_gradient_self);
    if (cell > 0)
      output->lower_derivative[cell] = cell_gradient*h*h*inverse_mu_4*
        (-inward_area*inward.left_derivative);
    if (cell + 1 < cells)
      output->upper_derivative[cell] = cell_gradient*h*h*inverse_mu_4*
        (outward_area*outward.right_derivative);
    if (!isfinite(output->pressure_inventory_rate[cell]) ||
        !isfinite(output->lower_derivative[cell]) ||
        !isfinite(output->diagonal_derivative[cell]) ||
        !isfinite(output->upper_derivative[cell]))
      return MUSEHANE_INVALID_INPUT;
  }
  return MUSEHANE_OK;
}
