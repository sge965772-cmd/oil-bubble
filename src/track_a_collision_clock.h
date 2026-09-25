#ifndef TRACK_A_COLLISION_CLOCK_H
#define TRACK_A_COLLISION_CLOCK_H

#include <float.h>
#include <math.h>
#include <stdbool.h>

/* The collision observer runs after a completed Basilisk DNS step.  Its t is
 * the interval end and dt is the completed interval length, while the film
 * contract labels an update by its interval start. */
static inline bool track_a_completed_step_interval_start (
  double completed_time, double completed_dt, double *interval_start)
{
  if (!interval_start || !isfinite(completed_time) || completed_time < 0. ||
      !isfinite(completed_dt) || completed_dt <= 0.)
    return false;
  double start = completed_time - completed_dt;
  const double scale = fmax(DBL_MIN,
    fmax(fabs(completed_time), fabs(completed_dt)));
  const double roundoff = 128.*DBL_EPSILON*scale;
  if (!isfinite(start) || start < -roundoff)
    return false;
  if (fabs(start) <= roundoff)
    start = 0.;
  *interval_start = start;
  return true;
}

#endif
