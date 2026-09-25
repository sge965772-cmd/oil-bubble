#ifndef FILM_MUSEHANE_PUBLISHED_EQUATION23_OPERATOR_H
#define FILM_MUSEHANE_PUBLISHED_EQUATION23_OPERATOR_H

#include "film_musehane/core.h"

typedef struct {
  const MusehaneSurfaceMesh *mesh;
  const unsigned char *active;
  const double *thickness;
  const double *pressure_gradient_cell;
  const double *pressure_gradient_face;
  double viscosity;
} MusehanePublishedEquation23Input;

typedef struct {
  double *pressure_inventory_rate;
  double *pressure_flux_integral_face;
  double *lower_derivative;
  double *diagonal_derivative;
  double *upper_derivative;
} MusehanePublishedEquation23Output;

/* Pressure terms and their tridiagonal thickness Jacobian from Musehane
   et al. (2018), equation (23). Face thickness is linearly interpolated
   between active cells and uses the active-cell value at an active-region
   edge, implementing the paper's homogeneous thickness-gradient condition. */
MusehaneStatus musehane_published_equation23_pressure (
  const MusehanePublishedEquation23Input *input,
  MusehanePublishedEquation23Output *output);

#endif
