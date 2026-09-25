#ifndef FILM_MUSEHANE_PRESSURE_DRIVER_H
#define FILM_MUSEHANE_PRESSURE_DRIVER_H

#include "film_musehane/core.h"

typedef enum {
  MUSEHANE_PRESSURE_DRIVER_OK = 0,
  MUSEHANE_PRESSURE_DRIVER_INVALID_INPUT = 1,
  MUSEHANE_PRESSURE_DRIVER_SINGULAR = 2,
  MUSEHANE_PRESSURE_DRIVER_RESIDUAL_MISMATCH = 3,
  MUSEHANE_PRESSURE_DRIVER_ALLOCATION_FAILED = 4
} MusehanePressureDriverStatus;

typedef struct {
  /* eta in (I - div(eta grad)) grad(p')_smooth = grad(p'). */
  double smoothing_length_squared;
  double linear_residual_relative_tolerance;
} MusehanePressureDriverConfig;

typedef struct {
  const MusehaneSurfaceMesh *mesh;
  const double *raw_pressure_gradient_tangent_cell;
  const double *capillary_force_density_tangent_cell;
} MusehanePressureDriverInput;

typedef struct {
  double *decontaminated_gradient_cell;
  double *smoothed_gradient_cell;
  double *smoothed_gradient_face;
  double maximum_linear_relative_residual;
} MusehanePressureDriverOutput;

/* Implements Musehane et al. (2018), equations (24)-(25), on the supplied
   axisymmetric surface mesh with homogeneous Neumann edge conditions. */
MusehanePressureDriverStatus musehane_prepare_pressure_driver (
  const MusehanePressureDriverConfig *config,
  const MusehanePressureDriverInput *input,
  MusehanePressureDriverOutput *output);

#endif
