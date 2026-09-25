#ifndef FILM_MUSEHANE_BETA_VISCOSITY_H
#define FILM_MUSEHANE_BETA_VISCOSITY_H

#include <stddef.h>

#include "film_model/contract.h"

typedef struct {
  size_t ring_count;
  const FilmRingObservation *rings;
  const unsigned char *active;
  double axial_padding_cells;
} MusehaneBetaViscosityInput;

/* Binary beta support for the equation-(28) replacement. The radial support
 * is the sampled annulus; the axial support spans the two interfaces plus a
 * grid-scaled allowance for the smoothed-marker tails used in equation (22). */
double musehane_beta_viscosity_at_point (
  const MusehaneBetaViscosityInput *input,
  double axial_position, double radius, double local_delta);

#endif
