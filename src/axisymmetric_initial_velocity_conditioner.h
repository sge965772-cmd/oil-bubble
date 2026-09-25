#ifndef AXISYMMETRIC_INITIAL_VELOCITY_CONDITIONER_H
#define AXISYMMETRIC_INITIAL_VELOCITY_CONDITIONER_H

#include <math.h>

static inline double axisymmetric_conditioning_clamp01 (double value)
{
  return value < 0. ? 0. : (value > 1. ? 1. : value);
}

static inline double axisymmetric_conditioning_smootherstep (double value)
{
  double s = axisymmetric_conditioning_clamp01(value);
  return s*s*s*(10. + s*(-15. + 6.*s));
}

static inline double axisymmetric_conditioning_window (double distance,
                                                       double core_radius,
                                                       double outer_radius)
{
  if (distance <= core_radius)
    return 1.;
  if (distance >= outer_radius || outer_radius <= core_radius)
    return 0.;
  return 1. - axisymmetric_conditioning_smootherstep(
    (distance - core_radius)/(outer_radius - core_radius));
}

// Stokes streamfunction: q_z = d(psi)/dr and q_r = -d(psi)/dz are
// metric face fluxes for an axisymmetric divergence-free velocity field.
static inline double axisymmetric_conditioning_streamfunction (
  double z, double r, double center_z, double axial_speed,
  double core_radius, double outer_radius)
{
  double dz = z - center_z;
  double distance = sqrt(dz*dz + r*r);
  double window = axisymmetric_conditioning_window(
    distance, core_radius, outer_radius);
  return 0.5*axial_speed*r*r*window;
}

static inline int axisymmetric_conditioning_solve_2x2 (
  const double response[2][2], const double target_change[2],
  double coefficient[2], double determinant_tolerance)
{
  double determinant = response[0][0]*response[1][1] -
                       response[0][1]*response[1][0];
  double scale = fmax(fabs(response[0][0]*response[1][1]),
                      fabs(response[0][1]*response[1][0]));
  if (!isfinite(determinant) ||
      fabs(determinant) <= determinant_tolerance*fmax(1., scale))
    return 0;
  coefficient[0] = (target_change[0]*response[1][1] -
                    response[0][1]*target_change[1])/determinant;
  coefficient[1] = (response[0][0]*target_change[1] -
                    target_change[0]*response[1][0])/determinant;
  return isfinite(coefficient[0]) && isfinite(coefficient[1]);
}

#endif
