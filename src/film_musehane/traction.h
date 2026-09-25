#ifndef FILM_MUSEHANE_TRACTION_H
#define FILM_MUSEHANE_TRACTION_H

#include <stddef.h>

typedef enum {
  MUSEHANE_TRACTION_OK = 0,
  MUSEHANE_TRACTION_INVALID_INPUT = 1,
  MUSEHANE_TRACTION_ALLOCATION_FAILED = 2
} MusehaneTractionStatus;

typedef struct {
  size_t cell_count;
  double water_viscosity;
  const unsigned char *active;
  const double *thickness;
  const double *pressure_gradient_tangent;
  const double *lower_tangential_velocity;
  const double *upper_tangential_velocity;
} MusehaneTractionInput;

typedef struct {
  double *lower_tangential_traction;
  double *upper_tangential_traction;
  double maximum_constitutive_residual;
} MusehaneTractionOutput;

/* Musehane et al. (2018), equation (27). The lower envelope is droplet 1
   and the upper envelope is droplet 2. Values are traction densities, not
   ring-integrated forces. Output is committed only after full validation. */
MusehaneTractionStatus musehane_evaluate_traction (
  const MusehaneTractionInput *input, MusehaneTractionOutput *output);

#endif
