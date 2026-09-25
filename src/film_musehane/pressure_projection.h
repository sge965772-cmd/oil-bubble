#ifndef FILM_MUSEHANE_PRESSURE_PROJECTION_H
#define FILM_MUSEHANE_PRESSURE_PROJECTION_H

#include <stddef.h>

typedef struct {
  double x;
  double r;
} MusehanePlaneVector;

typedef enum {
  MUSEHANE_PRESSURE_PROJECTION_OK = 0,
  MUSEHANE_PRESSURE_PROJECTION_INVALID_INPUT = 1,
  MUSEHANE_PRESSURE_PROJECTION_ALLOCATION_FAILED = 2
} MusehanePressureProjectionStatus;

typedef struct {
  double maximum_tangent_norm_error;
  double maximum_static_residual;
} MusehanePressureProjectionAudit;

/* Projects equation (24) into the supplied lower-to-upper film tangent.
   Outputs are transactional: rejection leaves every caller buffer unchanged. */
MusehanePressureProjectionStatus musehane_project_pressure_driver (
  size_t count,
  const MusehanePlaneVector *unit_tangent,
  const MusehanePlaneVector *pressure_gradient,
  const MusehanePlaneVector *capillary_force_density,
  double *raw_pressure_gradient_tangent,
  double *capillary_force_density_tangent,
  double *decontaminated_pressure_gradient_tangent,
  MusehanePressureProjectionAudit *audit);

#endif
