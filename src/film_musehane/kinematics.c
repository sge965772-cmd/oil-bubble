#include "film_musehane/kinematics.h"

#include <math.h>

static int vector_is_finite (const double vector[2])
{
  return vector && isfinite(vector[0]) && isfinite(vector[1]);
}

static double dot (const double first[2], const double second[2])
{
  return first[0]*second[0] + first[1]*second[1];
}

MusehaneKinematicsStatus musehane_project_pair_kinematics (
  const MusehanePairKinematicsInput *input,
  MusehanePairKinematics *output)
{
  if (!input || !output ||
      !vector_is_finite(input->lower_to_upper_normal) ||
      !vector_is_finite(input->lower_velocity) ||
      !vector_is_finite(input->upper_velocity))
    return MUSEHANE_KINEMATICS_INVALID_INPUT;
  const double magnitude = hypot(input->lower_to_upper_normal[0],
    input->lower_to_upper_normal[1]);
  if (!isfinite(magnitude) || magnitude <= 0.)
    return MUSEHANE_KINEMATICS_INVALID_INPUT;

  MusehanePairKinematics candidate = {.valid = true};
  candidate.normal[0] = input->lower_to_upper_normal[0]/magnitude;
  candidate.normal[1] = input->lower_to_upper_normal[1]/magnitude;
  candidate.tangent[0] = -candidate.normal[1];
  candidate.tangent[1] = candidate.normal[0];
  candidate.lower_normal_velocity = dot(input->lower_velocity,
    candidate.normal);
  candidate.upper_normal_velocity = dot(input->upper_velocity,
    candidate.normal);
  candidate.lower_tangential_velocity = dot(input->lower_velocity,
    candidate.tangent);
  candidate.upper_tangential_velocity = dot(input->upper_velocity,
    candidate.tangent);
  candidate.closing_velocity = candidate.lower_normal_velocity -
    candidate.upper_normal_velocity;
  candidate.relative_tangential_velocity =
    candidate.lower_tangential_velocity -
    candidate.upper_tangential_velocity;
  if (!isfinite(candidate.normal[0]) || !isfinite(candidate.normal[1]) ||
      !isfinite(candidate.tangent[0]) || !isfinite(candidate.tangent[1]) ||
      !isfinite(candidate.lower_normal_velocity) ||
      !isfinite(candidate.upper_normal_velocity) ||
      !isfinite(candidate.lower_tangential_velocity) ||
      !isfinite(candidate.upper_tangential_velocity) ||
      !isfinite(candidate.closing_velocity) ||
      !isfinite(candidate.relative_tangential_velocity))
    return MUSEHANE_KINEMATICS_NONFINITE_RESULT;
  *output = candidate;
  return MUSEHANE_KINEMATICS_OK;
}
