#ifndef AXISYMMETRIC_PROJECTION_AUDIT_H
#define AXISYMMETRIC_PROJECTION_AUDIT_H

#include "axisymmetric_projection_metric.h"

static inline double axisymmetric_projected_cell_volume_change (void)
{
  double maximum_change = 0.;
  foreach(reduction(max:maximum_change)) {
    double metric_flux_difference = 0.;
    foreach_dimension()
      metric_flux_difference += uf.x[1] - uf.x[];
    double local_change = axisymmetric_fractional_cell_volume_change(
      metric_flux_difference, dt, Delta, cm[]);
    maximum_change = max(maximum_change, local_change);
  }
  return maximum_change;
}

#endif
