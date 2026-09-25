#ifndef FILM_MUSEHANE_DNS_OBSERVATION_H
#define FILM_MUSEHANE_DNS_OBSERVATION_H

#include <stdbool.h>
#include <stddef.h>

#include "film_musehane/pressure_projection.h"

typedef struct {
  double radius;
  double width;
  double area_weight;
  double lower_position;
  double upper_position;
  double local_delta;
  MusehanePlaneVector lower_velocity;
  MusehanePlaneVector upper_velocity;
  double water_pressure;
  MusehanePlaneVector pressure_gradient;
  MusehanePlaneVector capillary_force_density;
  MusehanePlaneVector lower_to_upper_normal;
  bool geometry_valid;
  bool lower_identity_valid;
  bool upper_identity_valid;
} MusehaneDNSRingSample;

/* Track A owns this representation.  Conversion to the repository-wide
 * FilmRingObservation contract belongs in the outer adapter layer. */
typedef struct {
  double radius;
  double width;
  double area_weight;
  double lower_position;
  double upper_position;
  double gap;
  double midpoint;
  double local_delta;
  double lower_normal_velocity;
  double upper_normal_velocity;
  double lower_tangential_velocity;
  double upper_tangential_velocity;
  double water_pressure;
  double pressure_gradient_tangent;
  double capillary_force_density_tangent;
  double lower_normal_x;
  double lower_normal_r;
  double upper_normal_x;
  double upper_normal_r;
  bool geometry_valid;
  bool lower_identity_valid;
  bool upper_identity_valid;
} MusehaneDNSRingObservation;

typedef enum {
  MUSEHANE_DNS_OBSERVATION_OK = 0,
  MUSEHANE_DNS_OBSERVATION_INVALID_INPUT = 1,
  MUSEHANE_DNS_OBSERVATION_ALLOCATION_FAILED = 2
} MusehaneDNSObservationStatus;

typedef struct {
  double maximum_normal_norm_error;
  double maximum_decontaminated_pressure_gradient;
} MusehaneDNSObservationAudit;

MusehaneDNSObservationStatus musehane_assemble_dns_observations (
  size_t count,
  const MusehaneDNSRingSample *samples,
  MusehaneDNSRingObservation *observations,
  MusehaneDNSObservationAudit *audit);

#endif
