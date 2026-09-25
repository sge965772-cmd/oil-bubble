#include "film_musehane/pressure_projection.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int distinct_outputs (
  const double *raw, const double *capillary, const double *decontaminated)
{
  return raw && capillary && decontaminated &&
    raw != capillary && raw != decontaminated &&
    capillary != decontaminated;
}

MusehanePressureProjectionStatus musehane_project_pressure_driver (
  size_t count,
  const MusehanePlaneVector *unit_tangent,
  const MusehanePlaneVector *pressure_gradient,
  const MusehanePlaneVector *capillary_force_density,
  double *raw_pressure_gradient_tangent,
  double *capillary_force_density_tangent,
  double *decontaminated_pressure_gradient_tangent,
  MusehanePressureProjectionAudit *audit)
{
  if (count == 0 || count > SIZE_MAX/(3*sizeof(double)) ||
      !unit_tangent || !pressure_gradient || !capillary_force_density ||
      !distinct_outputs(raw_pressure_gradient_tangent,
        capillary_force_density_tangent,
        decontaminated_pressure_gradient_tangent) || !audit)
    return MUSEHANE_PRESSURE_PROJECTION_INVALID_INPUT;

  double *storage = malloc(3*count*sizeof *storage);
  if (!storage)
    return MUSEHANE_PRESSURE_PROJECTION_ALLOCATION_FAILED;
  double *raw = storage;
  double *capillary = raw + count;
  double *decontaminated = capillary + count;
  double maximum_norm_error = 0.;
  double maximum_static_residual = 0.;

  for (size_t sample_index = 0; sample_index < count; sample_index++) {
    const MusehanePlaneVector tangent = unit_tangent[sample_index];
    const MusehanePlaneVector pressure = pressure_gradient[sample_index];
    const MusehanePlaneVector force = capillary_force_density[sample_index];
    if (!isfinite(tangent.x) || !isfinite(tangent.r) ||
        !isfinite(pressure.x) || !isfinite(pressure.r) ||
        !isfinite(force.x) || !isfinite(force.r)) {
      free(storage);
      return MUSEHANE_PRESSURE_PROJECTION_INVALID_INPUT;
    }
    const double norm = hypot(tangent.x, tangent.r);
    const double norm_error = fabs(norm - 1.);
    if (!isfinite(norm) || norm_error > 1024.*DBL_EPSILON) {
      free(storage);
      return MUSEHANE_PRESSURE_PROJECTION_INVALID_INPUT;
    }
    raw[sample_index] = pressure.x*tangent.x + pressure.r*tangent.r;
    capillary[sample_index] = force.x*tangent.x + force.r*tangent.r;
    decontaminated[sample_index] =
      raw[sample_index] - capillary[sample_index];
    if (!isfinite(raw[sample_index]) ||
        !isfinite(capillary[sample_index]) ||
        !isfinite(decontaminated[sample_index])) {
      free(storage);
      return MUSEHANE_PRESSURE_PROJECTION_INVALID_INPUT;
    }
    maximum_norm_error = fmax(maximum_norm_error, norm_error);
    maximum_static_residual = fmax(maximum_static_residual,
      fabs(decontaminated[sample_index]));
  }

  memcpy(raw_pressure_gradient_tangent, raw, count*sizeof *raw);
  memcpy(capillary_force_density_tangent, capillary,
    count*sizeof *capillary);
  memcpy(decontaminated_pressure_gradient_tangent, decontaminated,
    count*sizeof *decontaminated);
  *audit = (MusehanePressureProjectionAudit) {
    .maximum_tangent_norm_error = maximum_norm_error,
    .maximum_static_residual = maximum_static_residual
  };
  free(storage);
  return MUSEHANE_PRESSURE_PROJECTION_OK;
}
