#ifndef FILM_MUSEHANE_TRACTION_DEPOSITION_H
#define FILM_MUSEHANE_TRACTION_DEPOSITION_H

#include <stddef.h>

#include "film_model/contract.h"

typedef enum {
  MUSEHANE_TRACTION_DEPOSITION_OK = 0,
  MUSEHANE_TRACTION_DEPOSITION_INVALID_INPUT = 1,
  MUSEHANE_TRACTION_DEPOSITION_UNBALANCED_PAIR = 2,
  MUSEHANE_TRACTION_DEPOSITION_ALLOCATION_FAILED = 3
} MusehaneTractionDepositionStatus;

typedef struct {
  size_t ring_count;
  const FilmRingObservation *rings;
  const FilmRingTraction *tractions;
} MusehaneTractionDepositionInput;

typedef struct {
  double *lower_force_x;
  double *lower_force_r;
  double *upper_force_x;
  double *upper_force_r;
  double requested_lower_power;
  double requested_upper_power;
  double requested_total_power;
  double requested_resolved_reaction_force;
  double requested_axial_force;
  /* Signed meridional radial sum; a complete axisymmetric ring has zero
     Cartesian radial resultant after azimuthal integration. */
  double requested_meridional_radial_force;
  double requested_meridional_moment;
  double maximum_scalar_pair_residual;
} MusehaneTractionDepositionOutput;

/* Convert ring-integrated scalar tangential forces into meridional x/r
   components. The two interface forces plus the declared resolved-film
   reaction must close ring by ring. Output is transactional and no Eulerian
   field is modified. */
MusehaneTractionDepositionStatus musehane_plan_traction_deposition (
  const MusehaneTractionDepositionInput *input,
  MusehaneTractionDepositionOutput *output);

#endif
