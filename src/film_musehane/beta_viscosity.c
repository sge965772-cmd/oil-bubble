#include "film_musehane/beta_viscosity.h"

#include <math.h>

double musehane_beta_viscosity_at_point (
  const MusehaneBetaViscosityInput *input,
  double axial_position, double radius, double local_delta)
{
  if (!input || input->ring_count == 0 || !input->rings || !input->active ||
      !isfinite(input->axial_padding_cells) ||
      input->axial_padding_cells < 0. || !isfinite(axial_position) ||
      !isfinite(radius) || radius < 0. || !isfinite(local_delta) ||
      local_delta <= 0.)
    return 0.;

  for (size_t i = 0; i < input->ring_count; i++) {
    const FilmRingObservation *ring = &input->rings[i];
    if (!input->active[i] || !ring->geometry_valid ||
        !isfinite(ring->radius) || ring->radius < 0. ||
        !isfinite(ring->width) || ring->width <= 0. ||
        !isfinite(ring->lower_position) ||
        !isfinite(ring->upper_position) ||
        ring->upper_position <= ring->lower_position)
      continue;
    const double inner = fmax(0., ring->radius - .5*ring->width);
    const double outer = ring->radius + .5*ring->width;
    const double padding = input->axial_padding_cells*local_delta;
    if (radius >= inner && radius <= outer &&
        axial_position >= ring->lower_position - padding &&
        axial_position <= ring->upper_position + padding)
      return 1.;
  }
  return 0.;
}
