#ifndef AXISYMMETRIC_PROJECTION_METRIC_H
#define AXISYMMETRIC_PROJECTION_METRIC_H

#include <math.h>

static inline double axisymmetric_fractional_cell_volume_change (
  double metric_flux_difference,
  double timestep,
  double cell_size,
  double cell_metric)
{
  return fabs(metric_flux_difference)*timestep/
         (cell_size*fmax(cell_metric, 1.e-300));
}

#endif
