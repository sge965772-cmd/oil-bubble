#ifndef AXISYMMETRIC_FILM_RUNTIME_H
#define AXISYMMETRIC_FILM_RUNTIME_H

#include "axisymmetric_film_state.h"

/* The implementation is compiled as ordinary C so Basilisk's dimensional
   interpreter does not enumerate the 32-bin state arrays. Callers must first
   validate both the active state and the observation. */
double axisymmetric_film_runtime_timestep_limit (
  const AxisymmetricFilmState * state,
  const AxisymmetricFilmObservation * observation);

#endif
