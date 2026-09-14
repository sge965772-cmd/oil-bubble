#include "axisymmetric_film_runtime.h"

double axisymmetric_film_runtime_timestep_limit (
  const AxisymmetricFilmState * state,
  const AxisymmetricFilmObservation * observation)
{
  double limit = INFINITY;
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++) {
    const double closing = observation->lower_normal_velocity[bin] -
      observation->upper_normal_velocity[bin];
    if (closing > 0.) {
      const double candidate = AXISYMMETRIC_FILM_STATE_CFL*
        state->thickness[bin]/closing;
      limit = fmin(limit, candidate);
    }
  }
  return limit;
}
