#ifndef CONTACT_RESOLUTION_METRIC_H
#define CONTACT_RESOLUTION_METRIC_H

#include <math.h>

static inline double contact_resolution_cells (
  double gap, double lower_interface_delta, double upper_interface_delta)
{
  if (!isfinite(gap) || !isfinite(lower_interface_delta) ||
      !isfinite(upper_interface_delta) || lower_interface_delta <= 0. ||
      upper_interface_delta <= 0.)
    return NAN;
  double controlling_delta = fmax(lower_interface_delta, upper_interface_delta);
  return gap/controlling_delta;
}

#endif
