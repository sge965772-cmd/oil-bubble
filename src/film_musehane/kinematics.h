#ifndef FILM_MUSEHANE_KINEMATICS_H
#define FILM_MUSEHANE_KINEMATICS_H

#include <stdbool.h>

typedef enum {
  MUSEHANE_KINEMATICS_OK = 0,
  MUSEHANE_KINEMATICS_INVALID_INPUT = 1,
  MUSEHANE_KINEMATICS_NONFINITE_RESULT = 2
} MusehaneKinematicsStatus;

typedef struct {
  /* The normal orientation is part of the interface: lower envelope to upper
     envelope. Callers must not independently orient the two interfaces. */
  double lower_to_upper_normal[2];
  double lower_velocity[2];
  double upper_velocity[2];
} MusehanePairKinematicsInput;

typedef struct {
  bool valid;
  double normal[2];
  double tangent[2];
  double lower_normal_velocity;
  double upper_normal_velocity;
  double lower_tangential_velocity;
  double upper_tangential_velocity;
  double closing_velocity;
  double relative_tangential_velocity;
} MusehanePairKinematics;

/* Projects both interface velocities into one common orthonormal frame.
   Pair-relative outputs are invariant under rigid translation and simultaneous
   rotation of geometry and velocity. Output is transactional. */
MusehaneKinematicsStatus musehane_project_pair_kinematics (
  const MusehanePairKinematicsInput *input,
  MusehanePairKinematics *output);

#endif
