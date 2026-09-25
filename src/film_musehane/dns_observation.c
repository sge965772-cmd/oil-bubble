#include "film_musehane/dns_observation.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int finite_vector (MusehanePlaneVector vector)
{
  return isfinite(vector.x) && isfinite(vector.r);
}

MusehaneDNSObservationStatus musehane_assemble_dns_observations (
  size_t count,
  const MusehaneDNSRingSample *samples,
  MusehaneDNSRingObservation *observations,
  MusehaneDNSObservationAudit *audit)
{
  if (count == 0 ||
      count > SIZE_MAX/sizeof(MusehaneDNSRingObservation) ||
      count > SIZE_MAX/(6*sizeof(double)) || !samples || !observations ||
      !audit)
    return MUSEHANE_DNS_OBSERVATION_INVALID_INPUT;

  MusehaneDNSRingObservation *candidate =
    calloc(count, sizeof *candidate);
  MusehanePlaneVector *vectors = malloc(3*count*sizeof *vectors);
  double *projections = malloc(3*count*sizeof *projections);
  if (!candidate || !vectors || !projections) {
    free(candidate);
    free(vectors);
    free(projections);
    return MUSEHANE_DNS_OBSERVATION_ALLOCATION_FAILED;
  }
  MusehanePlaneVector *tangent = vectors;
  MusehanePlaneVector *pressure = tangent + count;
  MusehanePlaneVector *capillary = pressure + count;
  double *raw_tangent = projections;
  double *capillary_tangent = raw_tangent + count;
  double *decontaminated = capillary_tangent + count;
  double maximum_normal_norm_error = 0.;

  for (size_t ring = 0; ring < count; ring++) {
    const MusehaneDNSRingSample sample = samples[ring];
    const double gap = sample.upper_position - sample.lower_position;
    const double normal_norm = hypot(sample.lower_to_upper_normal.x,
      sample.lower_to_upper_normal.r);
    const double normal_error = fabs(normal_norm - 1.);
    if (!isfinite(sample.radius) || sample.radius < 0. ||
        !isfinite(sample.width) || sample.width <= 0. ||
        !isfinite(sample.area_weight) || sample.area_weight <= 0. ||
        !isfinite(sample.lower_position) ||
        !isfinite(sample.upper_position) || !isfinite(gap) || gap <= 0. ||
        !isfinite(sample.local_delta) || sample.local_delta <= 0. ||
        !finite_vector(sample.lower_velocity) ||
        !finite_vector(sample.upper_velocity) ||
        !isfinite(sample.water_pressure) ||
        !finite_vector(sample.pressure_gradient) ||
        !finite_vector(sample.capillary_force_density) ||
        !finite_vector(sample.lower_to_upper_normal) ||
        sample.lower_to_upper_normal.x <= 0. ||
        !isfinite(normal_norm) || normal_error > 1024.*DBL_EPSILON ||
        !sample.geometry_valid || !sample.lower_identity_valid ||
        !sample.upper_identity_valid) {
      free(candidate);
      free(vectors);
      free(projections);
      return MUSEHANE_DNS_OBSERVATION_INVALID_INPUT;
    }
    tangent[ring] = (MusehanePlaneVector) {
      .x = -sample.lower_to_upper_normal.r,
      .r = sample.lower_to_upper_normal.x
    };
    pressure[ring] = sample.pressure_gradient;
    capillary[ring] = sample.capillary_force_density;
    maximum_normal_norm_error = fmax(maximum_normal_norm_error, normal_error);
  }

  MusehanePressureProjectionAudit projection_audit = {0};
  MusehanePressureProjectionStatus projection_status =
    musehane_project_pressure_driver(count, tangent, pressure, capillary,
      raw_tangent, capillary_tangent, decontaminated, &projection_audit);
  if (projection_status != MUSEHANE_PRESSURE_PROJECTION_OK) {
    free(candidate);
    free(vectors);
    free(projections);
    return projection_status == MUSEHANE_PRESSURE_PROJECTION_ALLOCATION_FAILED ?
      MUSEHANE_DNS_OBSERVATION_ALLOCATION_FAILED :
      MUSEHANE_DNS_OBSERVATION_INVALID_INPUT;
  }

  for (size_t ring = 0; ring < count; ring++) {
    const MusehaneDNSRingSample sample = samples[ring];
    const MusehanePlaneVector normal = sample.lower_to_upper_normal;
    candidate[ring] = (MusehaneDNSRingObservation) {
      .radius = sample.radius,
      .width = sample.width,
      .area_weight = sample.area_weight,
      .lower_position = sample.lower_position,
      .upper_position = sample.upper_position,
      .gap = sample.upper_position - sample.lower_position,
      .midpoint = .5*(sample.lower_position + sample.upper_position),
      .local_delta = sample.local_delta,
      .lower_normal_velocity = sample.lower_velocity.x*normal.x +
        sample.lower_velocity.r*normal.r,
      .upper_normal_velocity = sample.upper_velocity.x*normal.x +
        sample.upper_velocity.r*normal.r,
      .lower_tangential_velocity = sample.lower_velocity.x*tangent[ring].x +
        sample.lower_velocity.r*tangent[ring].r,
      .upper_tangential_velocity = sample.upper_velocity.x*tangent[ring].x +
        sample.upper_velocity.r*tangent[ring].r,
      .water_pressure = sample.water_pressure,
      .pressure_gradient_tangent = raw_tangent[ring],
      .capillary_force_density_tangent = capillary_tangent[ring],
      .lower_normal_x = normal.x,
      .lower_normal_r = normal.r,
      .upper_normal_x = normal.x,
      .upper_normal_r = normal.r,
      .geometry_valid = true,
      .lower_identity_valid = true,
      .upper_identity_valid = true
    };
  }
  memcpy(observations, candidate, count*sizeof *candidate);
  *audit = (MusehaneDNSObservationAudit) {
    .maximum_normal_norm_error = maximum_normal_norm_error,
    .maximum_decontaminated_pressure_gradient =
      projection_audit.maximum_static_residual
  };
  free(candidate);
  free(vectors);
  free(projections);
  return MUSEHANE_DNS_OBSERVATION_OK;
}
