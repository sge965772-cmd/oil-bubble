#ifndef AXISYMMETRIC_COMPOUND_PROFILE_H
#define AXISYMMETRIC_COMPOUND_PROFILE_H

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#ifndef AXISYMMETRIC_PROFILE_MAX_POINTS
# define AXISYMMETRIC_PROFILE_MAX_POINTS 4096
#endif

#define AXISYMMETRIC_PROFILE_PI 3.141592653589793238462643383279502884

typedef enum {
  AXISYMMETRIC_PROFILE_OUTER = 0,
  AXISYMMETRIC_PROFILE_GAS = 1
} AxisymmetricProfilePhase;

typedef struct {
  int count;
  double z[AXISYMMETRIC_PROFILE_MAX_POINTS];
  double outer_radius[AXISYMMETRIC_PROFILE_MAX_POINTS];
  double gas_radius[AXISYMMETRIC_PROFILE_MAX_POINTS];
  double outer_volume;
  double gas_volume;
  double oil_volume;
  double oil_fraction;
  double outer_centroid;
  double outer_support_min;
  double outer_support_max;
  double gas_support_min;
  double gas_support_max;
  double minimum_radial_film;
} AxisymmetricCompoundProfile;

static inline int axisymmetric_profile_error (
  char * error, size_t capacity, const char * message)
{
  if (error && capacity) {
    snprintf(error, capacity, "%s", message);
    error[capacity - 1] = '\0';
  }
  return 0;
}

static inline double axisymmetric_profile_segment_volume (
  double z0, double z1, double r0, double r1)
{
  double dz = z1 - z0;
  return AXISYMMETRIC_PROFILE_PI*dz*(r0*r0 + r0*r1 + r1*r1)/3.;
}

static inline double axisymmetric_profile_segment_first_moment (
  double z0, double z1, double r0, double r1)
{
  double dz = z1 - z0;
  double slope = r1 - r0;
  double radius_integral = r0*r0 + r0*slope + slope*slope/3.;
  double weighted_integral =
    r0*r0/2. + 2.*r0*slope/3. + slope*slope/4.;
  return AXISYMMETRIC_PROFILE_PI*dz*
         (z0*radius_integral + dz*weighted_integral);
}

static inline int axisymmetric_profile_load (
  AxisymmetricCompoundProfile * profile,
  const char * path,
  char * error,
  size_t error_capacity)
{
  memset(profile, 0, sizeof(*profile));
  FILE * input = fopen(path, "r");
  if (!input)
    return axisymmetric_profile_error(
      error, error_capacity, "cannot open profile file");

  char line[512];
  while (fgets(line, sizeof(line), input)) {
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\r')
      continue;
    for (char * cursor = line; *cursor; cursor++)
      if (*cursor == ',' || *cursor == ';')
        *cursor = ' ';
    double z, outer, gas;
    if (sscanf(line, "%lf %lf %lf", &z, &outer, &gas) != 3) {
      fclose(input);
      return axisymmetric_profile_error(
        error, error_capacity, "profile contains a malformed data row");
    }
    if (profile->count >= AXISYMMETRIC_PROFILE_MAX_POINTS) {
      fclose(input);
      return axisymmetric_profile_error(
        error, error_capacity, "profile exceeds point capacity");
    }
    if (!isfinite(z) || !isfinite(outer) || !isfinite(gas)) {
      fclose(input);
      return axisymmetric_profile_error(
        error, error_capacity, "profile contains a non-finite value");
    }
    if (profile->count && z <= profile->z[profile->count - 1]) {
      fclose(input);
      return axisymmetric_profile_error(
        error, error_capacity, "profile z coordinates are not strictly increasing");
    }
    gas = gas < 0. ? 0. : gas;
    if (outer < 0. || gas < 0. || gas > outer) {
      fclose(input);
      return axisymmetric_profile_error(
        error, error_capacity, "profile violates 0 <= gas radius <= outer radius");
    }
    int index = profile->count++;
    profile->z[index] = z;
    profile->outer_radius[index] = outer;
    profile->gas_radius[index] = gas;
  }
  fclose(input);

  if (profile->count < 4)
    return axisymmetric_profile_error(
      error, error_capacity, "profile contains fewer than four points");

  double outer_moment = 0.;
  profile->minimum_radial_film = HUGE_VAL;
  int positive_outer = 0, positive_gas = 0;
  int outer_support_ended = 0, gas_support_ended = 0;
  for (int index = 0; index < profile->count; index++) {
    double outer = profile->outer_radius[index];
    double gas = profile->gas_radius[index];
    if (outer > 0.) {
      if (!positive_outer)
        profile->outer_support_min = profile->z[index];
      profile->outer_support_max = profile->z[index];
      if (outer_support_ended)
        return axisymmetric_profile_error(
          error, error_capacity, "outer profile support is disconnected");
      positive_outer++;
    }
    else if (positive_outer)
      outer_support_ended = 1;
    if (gas > 0.) {
      if (gas_support_ended)
        return axisymmetric_profile_error(
          error, error_capacity, "gas profile support is disconnected");
      if (!positive_gas)
        profile->gas_support_min = profile->z[index];
      profile->gas_support_max = profile->z[index];
      positive_gas++;
      profile->minimum_radial_film = fmin(
        profile->minimum_radial_film, outer - gas);
    }
    else if (positive_gas)
      gas_support_ended = 1;
    if (index + 1 < profile->count) {
      double next_z = profile->z[index + 1];
      double next_outer = profile->outer_radius[index + 1];
      double next_gas = profile->gas_radius[index + 1];
      profile->outer_volume += axisymmetric_profile_segment_volume(
        profile->z[index], next_z, outer, next_outer);
      profile->gas_volume += axisymmetric_profile_segment_volume(
        profile->z[index], next_z, gas, next_gas);
      outer_moment += axisymmetric_profile_segment_first_moment(
        profile->z[index], next_z, outer, next_outer);
    }
  }
  if (!positive_outer || !positive_gas || profile->outer_volume <= 0. ||
      profile->gas_volume <= 0. || profile->gas_volume >= profile->outer_volume)
    return axisymmetric_profile_error(
      error, error_capacity, "profile has invalid gas or outer volume");
  if (profile->gas_support_min <= profile->outer_support_min ||
      profile->gas_support_max >= profile->outer_support_max ||
      profile->minimum_radial_film <= 0.)
    return axisymmetric_profile_error(
      error, error_capacity, "profile does not contain a positive closed oil coating");

  profile->oil_volume = profile->outer_volume - profile->gas_volume;
  profile->oil_fraction = profile->oil_volume/profile->outer_volume;
  profile->outer_centroid = outer_moment/profile->outer_volume;
  if (!isfinite(profile->minimum_radial_film))
    profile->minimum_radial_film = 0.;
  if (error && error_capacity)
    error[0] = '\0';
  return 1;
}

static inline double axisymmetric_profile_radius (
  const AxisymmetricCompoundProfile * profile,
  double local_z,
  AxisymmetricProfilePhase phase)
{
  if (profile->count < 2 || local_z < profile->z[0] ||
      local_z > profile->z[profile->count - 1])
    return -1.;
  int lower = 0, upper = profile->count - 1;
  while (upper - lower > 1) {
    int middle = (lower + upper)/2;
    if (profile->z[middle] <= local_z)
      lower = middle;
    else
      upper = middle;
  }
  double span = profile->z[upper] - profile->z[lower];
  double weight = span > 0. ? (local_z - profile->z[lower])/span : 0.;
  const double * radius = phase == AXISYMMETRIC_PROFILE_GAS ?
                          profile->gas_radius : profile->outer_radius;
  return (1. - weight)*radius[lower] + weight*radius[upper];
}

static inline double axisymmetric_profile_levelset (
  const AxisymmetricCompoundProfile * profile,
  double axial_coordinate,
  double radial_coordinate,
  double axial_offset,
  AxisymmetricProfilePhase phase)
{
  double radius = axisymmetric_profile_radius(
    profile, axial_coordinate - axial_offset, phase);
  return radius < 0. ? -1. : radius - fabs(radial_coordinate);
}

#endif
