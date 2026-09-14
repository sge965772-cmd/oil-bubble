#ifndef AXISYMMETRIC_FILM_STATE_H
#define AXISYMMETRIC_FILM_STATE_H

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>

#include "axisymmetric_matched_film.h"

#ifndef AXISYMMETRIC_FILM_STATE_BINS
# define AXISYMMETRIC_FILM_STATE_BINS AXISYMMETRIC_MATCHED_FILM_BINS
#endif

#define AXISYMMETRIC_FILM_STATE_HANDOFF_CELLS 4.
#define AXISYMMETRIC_FILM_STATE_CFL .20
#define AXISYMMETRIC_FILM_STATE_RELEASE_FACTOR 1.125
/* Raw area is reconstructed from finite fractions: this accepts only roundoff
   below one, not a user-tuned resolved-area threshold. */
#define AXISYMMETRIC_FILM_STATE_FULL_RAW_AREA_TOLERANCE (256.*DBL_EPSILON)
#define AXISYMMETRIC_FILM_STATE_REACQUISITION_OBSERVATIONS 4UL

enum {
  AXISYMMETRIC_FILM_STATE_INACTIVE = 0,
  AXISYMMETRIC_FILM_STATE_HANDOFF = 1,
  AXISYMMETRIC_FILM_STATE_ACTIVE = 2,
  AXISYMMETRIC_FILM_STATE_REACQUIRED = 3,
  AXISYMMETRIC_FILM_STATE_INVALID_CONFIG = 4,
  AXISYMMETRIC_FILM_STATE_INVALID_OBSERVATION = 5,
  AXISYMMETRIC_FILM_STATE_NONPOSITIVE_STATE = 6,
  AXISYMMETRIC_FILM_STATE_FLUX_SOLVE_FAILED = 7,
  AXISYMMETRIC_FILM_STATE_RESTART_MISMATCH = 8,
  AXISYMMETRIC_FILM_STATE_INVENTORY_MISMATCH = 9
};
typedef int AxisymmetricFilmStateStatus;

enum {
  AXISYMMETRIC_FILM_OBSERVATION_NONE = 0,
  AXISYMMETRIC_FILM_OBSERVATION_RAW_UNAVAILABLE = 1,
  AXISYMMETRIC_FILM_OBSERVATION_RAW_MALFORMED = 2,
  AXISYMMETRIC_FILM_OBSERVATION_LOW_QUALITY_ENTERING = 3,
  AXISYMMETRIC_FILM_OBSERVATION_CROSSED_ENTERING = 4,
  AXISYMMETRIC_FILM_OBSERVATION_UNDER_REFINED_ENTERING = 5,
  AXISYMMETRIC_FILM_OBSERVATION_RAW_GEOMETRY_UNRESOLVED = 6
};
typedef int AxisymmetricFilmObservationReason;

typedef struct {
  AxisymmetricMatchedFilmConfig reynolds;
  double minimum_raw_valid_area_fraction;
  double flux_balance_tolerance;
} AxisymmetricFilmStateConfig;

typedef struct {
  bool raw_input_available;
  double raw_valid_area_fraction;
  int raw_crossed_bins;
  bool contact_corridor_at_maxlevel;
  double time;
  AxisymmetricFilmInput raw;
  double lower_normal_velocity[AXISYMMETRIC_FILM_STATE_BINS];
  double upper_normal_velocity[AXISYMMETRIC_FILM_STATE_BINS];
  double lower_tangential_velocity[AXISYMMETRIC_FILM_STATE_BINS];
  double upper_tangential_velocity[AXISYMMETRIC_FILM_STATE_BINS];
} AxisymmetricFilmObservation;

typedef struct {
  bool active;
  double thickness[AXISYMMETRIC_FILM_STATE_BINS];
  double midpoint[AXISYMMETRIC_FILM_STATE_BINS];
  double radius[AXISYMMETRIC_FILM_STATE_BINS];
  double width[AXISYMMETRIC_FILM_STATE_BINS];
  double delta[AXISYMMETRIC_FILM_STATE_BINS];
  double face_flux[AXISYMMETRIC_FILM_STATE_BINS + 1];
  unsigned long activation_count;
  unsigned long update_count;
  unsigned long reacquisition_consecutive_observations;
  unsigned long cumulative_unresolved_observations;
  unsigned long cumulative_raw_crossing_observations;
  double source_time;
  double handoff_inventory;
  double inventory;
  double cumulative_edge_outflow;
  double cumulative_work;
} AxisymmetricFilmState;

typedef struct {
  bool valid;
  bool active;
  bool handed_off;
  bool reacquired;
  bool raw_available;
  bool raw_crossed;
  bool raw_geometry_unresolved;
  AxisymmetricFilmStateStatus status;
  AxisymmetricFilmObservationReason observation_reason;
  double inventory_before;
  double inventory_after;
  double expected_inventory_after;
  double inventory_closure_error;
  int nonpositive_bin;
  double nonpositive_thickness;
  double minimum_state_thickness;
  double minimum_raw_gap;
  double maximum_state_raw_gap_discrepancy;
  double raw_inventory;
  double state_raw_inventory_discrepancy;
  double state_raw_inventory_tolerance;
  unsigned long reacquisition_consecutive_observations;
  AxisymmetricReynoldsSolution reynolds_solution;
} AxisymmetricFilmStateResult;

static inline double axisymmetric_film_state_max (double first, double second)
{
  return first > second ? first : second;
}

static inline bool axisymmetric_film_state_config_is_valid (
  const AxisymmetricFilmStateConfig * config)
{
  return config && config->reynolds.bins == AXISYMMETRIC_FILM_STATE_BINS &&
    config->reynolds.bins == AXISYMMETRIC_MATCHED_FILM_BINS &&
    isfinite(config->reynolds.water_viscosity) &&
    config->reynolds.water_viscosity > 0. &&
    config->reynolds.match_cells == AXISYMMETRIC_FILM_STATE_HANDOFF_CELLS &&
    config->reynolds.release_factor ==
      AXISYMMETRIC_FILM_STATE_RELEASE_FACTOR &&
    isfinite(config->minimum_raw_valid_area_fraction) &&
    config->minimum_raw_valid_area_fraction >= .90 &&
    config->minimum_raw_valid_area_fraction <= 1. &&
    isfinite(config->flux_balance_tolerance) &&
    config->flux_balance_tolerance > 0.;
}

static inline double axisymmetric_film_state_ring_area (
  double radius, double width)
{
  const double inner_radius = radius - .5*width;
  const double outer_radius = radius + .5*width;
  return .5*(outer_radius*outer_radius - inner_radius*inner_radius);
}

static inline double axisymmetric_film_state_scale (double first, double second)
{
  return axisymmetric_film_state_max(DBL_MIN, axisymmetric_film_state_max(
    fabs(first), fabs(second)));
}

static inline bool axisymmetric_film_state_raw_geometry_is_well_formed (
  const AxisymmetricFilmObservation * observation)
{
  if (!observation || !observation->raw_input_available ||
      !isfinite(observation->raw_valid_area_fraction) ||
      observation->raw_valid_area_fraction < 0. ||
      observation->raw_valid_area_fraction > 1. ||
      observation->raw_crossed_bins < 0)
    return false;

  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++)
    if (!isfinite(observation->raw.radius[bin]) ||
        observation->raw.radius[bin] < 0. ||
        !isfinite(observation->raw.width[bin]) ||
        observation->raw.width[bin] <= 0. ||
        !isfinite(observation->raw.gap[bin]) ||
        !isfinite(observation->raw.delta[bin]) ||
        observation->raw.delta[bin] <= 0. ||
        !isfinite(observation->raw.lower_position[bin]) ||
        !isfinite(observation->raw.upper_position[bin]))
      return false;

  const double reference_width = observation->raw.width[0];
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++) {
    const double radius = observation->raw.radius[bin];
    const double width = observation->raw.width[bin];
    const double inner_radius = radius - .5*width;
    const double expected_inner_radius = bin == 0 ? 0. :
      observation->raw.radius[bin - 1] +
      .5*observation->raw.width[bin - 1];
    const double scale = axisymmetric_film_state_max(width, reference_width);
    const double gap_from_positions = observation->raw.upper_position[bin] -
      observation->raw.lower_position[bin];
    const double gap_subtraction_scale = axisymmetric_film_state_max(
      axisymmetric_film_state_scale(gap_from_positions,
        observation->raw.gap[bin]),
      axisymmetric_film_state_max(
        fabs(observation->raw.lower_position[bin]),
        fabs(observation->raw.upper_position[bin])));
    if (!isfinite(gap_from_positions) ||
        fabs(gap_from_positions - observation->raw.gap[bin]) >
          256.*DBL_EPSILON*gap_subtraction_scale ||
        inner_radius < 0. ||
        fabs(inner_radius - expected_inner_radius) >
          128.*DBL_EPSILON*scale ||
        fabs(width - reference_width) > 128.*DBL_EPSILON*scale)
      return false;
  }
  return true;
}

static inline bool axisymmetric_film_state_raw_observation_is_well_formed (
  const AxisymmetricFilmStateConfig * config,
  const AxisymmetricFilmObservation * observation, double dt)
{
  if (!config || !observation || !isfinite(observation->time) ||
      !isfinite(dt) || dt <= 0. ||
      !axisymmetric_film_state_raw_geometry_is_well_formed(observation))
    return false;
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++)
    if (observation->raw.gap[bin] <= 0. ||
        !isfinite(observation->lower_normal_velocity[bin]) ||
        !isfinite(observation->upper_normal_velocity[bin]) ||
        !isfinite(observation->lower_tangential_velocity[bin]) ||
        !isfinite(observation->upper_tangential_velocity[bin]))
      return false;
  return true;
}

static inline double axisymmetric_film_state_raw_minimum_gap (
  const AxisymmetricFilmObservation * observation)
{
  if (!axisymmetric_film_state_raw_geometry_is_well_formed(observation))
    return NAN;
  double minimum_gap = observation->raw.gap[0];
  for (int bin = 1; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++)
    minimum_gap = fmin(minimum_gap, observation->raw.gap[bin]);
  return minimum_gap;
}

static inline bool axisymmetric_film_state_raw_is_crossed (
  const AxisymmetricFilmObservation * observation)
{
  return axisymmetric_film_state_raw_geometry_is_well_formed(observation) &&
    (observation->raw_crossed_bins != 0 ||
     axisymmetric_film_state_raw_minimum_gap(observation) <= 0.);
}

static inline bool axisymmetric_film_state_raw_matches_state_mesh (
  const AxisymmetricFilmState * state,
  const AxisymmetricFilmObservation * observation)
{
  if (!state || !axisymmetric_film_state_raw_geometry_is_well_formed(
        observation))
    return false;
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++) {
    const double scale = axisymmetric_film_state_scale(state->radius[bin],
      observation->raw.radius[bin]);
    const double width_scale = axisymmetric_film_state_scale(state->width[bin],
      observation->raw.width[bin]);
    if (fabs(state->radius[bin] - observation->raw.radius[bin]) >
          128.*DBL_EPSILON*scale ||
        fabs(state->width[bin] - observation->raw.width[bin]) >
          128.*DBL_EPSILON*width_scale)
      return false;
  }
  return true;
}

static inline double axisymmetric_film_state_raw_inventory (
  const AxisymmetricFilmObservation * observation)
{
  if (!axisymmetric_film_state_raw_geometry_is_well_formed(observation))
    return NAN;
  const long double two_pi = 2.L*acosl(-1.L);
  long double inventory = 0.L;
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++) {
    const long double inner_radius = (long double) observation->raw.radius[bin]
      - .5L*(long double) observation->raw.width[bin];
    const long double outer_radius = (long double) observation->raw.radius[bin]
      + .5L*(long double) observation->raw.width[bin];
    const long double area = .5L*(outer_radius*outer_radius -
      inner_radius*inner_radius);
    inventory += two_pi*area*(long double) observation->raw.gap[bin];
  }
  return isfinite(inventory) && fabsl(inventory) <= (long double) DBL_MAX ?
    (double) inventory : NAN;
}

static inline double axisymmetric_film_state_raw_inventory_tolerance (
  const AxisymmetricFilmState * state)
{
  if (!state)
    return NAN;
  const long double two_pi = 2.L*acosl(-1.L);
  long double tolerance = 0.L;
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++) {
    const long double inner_radius = (long double) state->radius[bin] -
      .5L*(long double) state->width[bin];
    const long double outer_radius = (long double) state->radius[bin] +
      .5L*(long double) state->width[bin];
    const long double area = .5L*(outer_radius*outer_radius -
      inner_radius*inner_radius);
    tolerance += two_pi*area*.5L*(long double) state->delta[bin];
  }
  return isfinite(tolerance) && tolerance <= (long double) DBL_MAX ?
    (double) tolerance : NAN;
}

static inline double axisymmetric_film_state_inventory (
  const AxisymmetricFilmState * state)
{
  if (!state)
    return NAN;
  const double two_pi = 2.*acos(-1.);
  double inventory = 0.;
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++) {
    const double area = axisymmetric_film_state_ring_area(
      state->radius[bin], state->width[bin]);
    if (!isfinite(area) || area <= 0. ||
        !isfinite(state->thickness[bin]))
      return NAN;
    inventory += two_pi*area*state->thickness[bin];
  }
  return inventory;
}

static inline double axisymmetric_film_state_minimum_thickness (
  const AxisymmetricFilmState * state)
{
  if (!state)
    return NAN;
  double minimum = state->thickness[0];
  if (!isfinite(minimum))
    return NAN;
  for (int bin = 1; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++) {
    if (!isfinite(state->thickness[bin]))
      return NAN;
    minimum = fmin(minimum, state->thickness[bin]);
  }
  return minimum;
}

static inline bool axisymmetric_film_state_kinematics_are_valid (
  const AxisymmetricFilmObservation * observation)
{
  if (!observation)
    return false;
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++)
    if (!isfinite(observation->lower_normal_velocity[bin]) ||
        !isfinite(observation->upper_normal_velocity[bin]) ||
        !isfinite(observation->lower_tangential_velocity[bin]) ||
        !isfinite(observation->upper_tangential_velocity[bin]))
      return false;
  return true;
}

static inline bool axisymmetric_film_state_active_state_is_valid (
  const AxisymmetricFilmState * state)
{
  if (!state || !state->active || !isfinite(state->source_time) ||
      !isfinite(state->handoff_inventory) || state->handoff_inventory <= 0. ||
      !isfinite(state->inventory) || state->inventory <= 0. ||
      !isfinite(state->cumulative_edge_outflow) ||
      !isfinite(state->cumulative_work))
    return false;
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++) {
    const double inner_radius = state->radius[bin] - .5*state->width[bin];
    const double expected_inner_radius = bin == 0 ? 0. :
      state->radius[bin - 1] + .5*state->width[bin - 1];
    const double width_scale = axisymmetric_film_state_max(
      state->width[bin], state->width[0]);
    if (!isfinite(state->thickness[bin]) || state->thickness[bin] <= 0. ||
        !isfinite(state->midpoint[bin]) || !isfinite(state->radius[bin]) ||
        !isfinite(state->width[bin]) || state->width[bin] <= 0. ||
        !isfinite(state->delta[bin]) || state->delta[bin] <= 0. ||
        inner_radius < 0. ||
        fabs(inner_radius - expected_inner_radius) >
          128.*DBL_EPSILON*width_scale ||
        fabs(state->width[bin] - state->width[0]) >
          128.*DBL_EPSILON*width_scale)
      return false;
  }
  for (int face = 0; face <= AXISYMMETRIC_FILM_STATE_BINS; face++)
    if (!isfinite(state->face_flux[face]))
      return false;
  return true;
}

static inline bool axisymmetric_film_state_active_ledger_is_valid (
  const AxisymmetricFilmState * state)
{
  if (!state || !isfinite(state->source_time) ||
      state->activation_count == 0 ||
      state->activation_count == ULONG_MAX ||
      state->update_count == ULONG_MAX ||
      state->reacquisition_consecutive_observations >=
        AXISYMMETRIC_FILM_STATE_REACQUISITION_OBSERVATIONS ||
      state->cumulative_unresolved_observations == ULONG_MAX ||
      state->cumulative_raw_crossing_observations == ULONG_MAX ||
      !isfinite(state->handoff_inventory) || state->handoff_inventory <= 0. ||
      !isfinite(state->inventory) || state->inventory <= 0. ||
      !isfinite(state->cumulative_edge_outflow) ||
      fabs(state->cumulative_edge_outflow) > .5*DBL_MAX ||
      !isfinite(state->cumulative_work))
    return false;
  for (int face = 0; face <= AXISYMMETRIC_FILM_STATE_BINS; face++)
    if (!isfinite(state->face_flux[face]))
      return false;
  return true;
}

static inline bool axisymmetric_film_state_close (double first, double second)
{
  return isfinite(first) && isfinite(second) && fabs(first - second) <=
    1.e-12*axisymmetric_film_state_scale(first, second);
}

static inline bool axisymmetric_film_state_active_observation_is_valid (
  const AxisymmetricFilmStateConfig * config,
  const AxisymmetricFilmObservation * observation, double dt)
{
  return config && observation && isfinite(observation->time) &&
    isfinite(dt) && dt > 0. &&
    axisymmetric_film_state_kinematics_are_valid(observation);
}

static inline void axisymmetric_film_state_set_raw_diagnostics (
  AxisymmetricFilmStateResult * result,
  const AxisymmetricFilmState * state,
  const AxisymmetricFilmObservation * observation)
{
  if (!result)
    return;
  result->minimum_raw_gap = NAN;
  result->maximum_state_raw_gap_discrepancy = NAN;
  result->raw_inventory = NAN;
  result->state_raw_inventory_discrepancy = NAN;
  result->state_raw_inventory_tolerance = NAN;
  result->raw_available = observation && observation->raw_input_available;
  result->raw_crossed = observation && observation->raw_crossed_bins > 0;
  if (!observation || !observation->raw_input_available) {
    result->raw_geometry_unresolved = true;
    result->observation_reason =
      AXISYMMETRIC_FILM_OBSERVATION_RAW_GEOMETRY_UNRESOLVED;
    return;
  }
  if (!axisymmetric_film_state_raw_geometry_is_well_formed(observation)) {
    result->raw_geometry_unresolved = true;
    result->observation_reason =
      AXISYMMETRIC_FILM_OBSERVATION_RAW_GEOMETRY_UNRESOLVED;
    return;
  }

  result->minimum_raw_gap = axisymmetric_film_state_raw_minimum_gap(
    observation);
  result->raw_crossed |= axisymmetric_film_state_raw_is_crossed(observation);
  result->raw_geometry_unresolved = result->raw_crossed;
  if (result->raw_geometry_unresolved)
    result->observation_reason =
      AXISYMMETRIC_FILM_OBSERVATION_RAW_GEOMETRY_UNRESOLVED;
  if (!state || !axisymmetric_film_state_raw_matches_state_mesh(state,
        observation))
    return;

  result->raw_inventory = axisymmetric_film_state_raw_inventory(observation);
  result->state_raw_inventory_tolerance =
    axisymmetric_film_state_raw_inventory_tolerance(state);
  const double state_inventory = axisymmetric_film_state_inventory(state);
  if (isfinite(result->raw_inventory) && isfinite(state_inventory))
    result->state_raw_inventory_discrepancy = fabs(result->raw_inventory -
      state_inventory);
  double maximum_discrepancy = 0.;
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++)
    maximum_discrepancy = fmax(maximum_discrepancy,
      fabs(observation->raw.gap[bin] - state->thickness[bin]));
  result->maximum_state_raw_gap_discrepancy = maximum_discrepancy;
}

static inline bool axisymmetric_film_state_can_reacquire (
  const AxisymmetricFilmStateConfig * config,
  const AxisymmetricFilmState * state,
  const AxisymmetricFilmObservation * observation,
  const AxisymmetricFilmStateResult * diagnostics)
{
  if (!config || !state || !observation || !diagnostics ||
      !axisymmetric_film_state_raw_geometry_is_well_formed(observation) ||
      !axisymmetric_film_state_raw_matches_state_mesh(state, observation) ||
      observation->raw_valid_area_fraction < 1. -
        AXISYMMETRIC_FILM_STATE_FULL_RAW_AREA_TOLERANCE ||
      observation->raw_crossed_bins != 0 || diagnostics->raw_crossed ||
      !isfinite(diagnostics->state_raw_inventory_discrepancy) ||
      !isfinite(diagnostics->state_raw_inventory_tolerance))
    return false;
  /* The physical bound is the sum of per-bin 0.5 Delta discrepancies.  This
     allowance only covers independently rounded inventory reductions. */
  const double inventory_roundoff = 256.*DBL_EPSILON*
    axisymmetric_film_state_max(DBL_MIN, axisymmetric_film_state_max(
      fabs(diagnostics->raw_inventory), fabs(diagnostics->state_raw_inventory_tolerance)));
  if (diagnostics->state_raw_inventory_discrepancy >
      diagnostics->state_raw_inventory_tolerance + inventory_roundoff)
    return false;
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++) {
    const double resolved_gap = AXISYMMETRIC_FILM_STATE_RELEASE_FACTOR*
      AXISYMMETRIC_FILM_STATE_HANDOFF_CELLS*state->delta[bin];
    const double gap_agreement_tolerance = .5*state->delta[bin];
    const double gap_agreement_roundoff = 128.*DBL_EPSILON*
      axisymmetric_film_state_scale(observation->raw.gap[bin],
        state->thickness[bin]);
    const double closing = observation->lower_normal_velocity[bin] -
      observation->upper_normal_velocity[bin];
    if (!isfinite(resolved_gap) || !isfinite(closing) ||
        !(state->thickness[bin] > resolved_gap) ||
        !(observation->raw.gap[bin] > resolved_gap) ||
        fabs(observation->raw.gap[bin] - state->thickness[bin]) >
          gap_agreement_tolerance + gap_agreement_roundoff || !(closing < 0.))
      return false;
  }
  return true;
}

static inline double axisymmetric_film_state_timestep_limit_for_valid_state (
  const AxisymmetricFilmState * state,
  const AxisymmetricFilmObservation * observation)
{
  double limit = INFINITY;
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++) {
    const double closing = observation->lower_normal_velocity[bin] -
      observation->upper_normal_velocity[bin];
    const double candidate = AXISYMMETRIC_FILM_STATE_CFL*
      state->thickness[bin]/fmax(closing, 0.);
    limit = fmin(limit, candidate);
  }
  return limit;
}

static inline double axisymmetric_film_state_timestep_limit (
  const AxisymmetricFilmStateConfig * config,
  const AxisymmetricFilmState * state,
  const AxisymmetricFilmObservation * observation)
{
  if (!axisymmetric_film_state_config_is_valid(config) ||
      !axisymmetric_film_state_active_state_is_valid(state) ||
      !axisymmetric_film_state_active_ledger_is_valid(state) ||
      !observation || !isfinite(observation->time) ||
      !axisymmetric_film_state_kinematics_are_valid(observation))
    return NAN;

  double limit = INFINITY;
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++) {
    const double closing = observation->lower_normal_velocity[bin] -
      observation->upper_normal_velocity[bin];
    if (!isfinite(closing))
      return NAN;
    if (closing > 0.) {
      const long double candidate_long =
        (long double) AXISYMMETRIC_FILM_STATE_CFL*
        (long double) state->thickness[bin]/(long double) closing;
      double candidate;
      if (candidate_long > (long double) DBL_MAX)
        candidate = DBL_MAX;
      else {
        if (!isfinite(candidate_long) || !(candidate_long > 0.L))
          return NAN;
        candidate = (double) candidate_long;
        if (!isfinite(candidate) || candidate <= 0.)
          return NAN;
      }
      limit = fmin(limit, candidate);
    }
  }
  return limit;
}

static inline AxisymmetricFilmInput axisymmetric_film_state_reynolds_input (
  const AxisymmetricFilmState * state,
  const AxisymmetricFilmObservation * observation)
{
  AxisymmetricFilmInput input = {0};
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++) {
    input.radius[bin] = state->radius[bin];
    input.width[bin] = state->width[bin];
    input.gap[bin] = state->thickness[bin];
    input.delta[bin] = state->delta[bin];
    input.lower_position[bin] = state->midpoint[bin] -
      .5*state->thickness[bin];
    input.upper_position[bin] = state->midpoint[bin] +
      .5*state->thickness[bin];
    input.relative_normal_velocity[bin] =
      observation->lower_normal_velocity[bin] -
      observation->upper_normal_velocity[bin];
    input.lower_tangential_velocity[bin] =
      observation->lower_tangential_velocity[bin];
    input.upper_tangential_velocity[bin] =
      observation->upper_tangential_velocity[bin];
  }
  return input;
}

static inline AxisymmetricFilmStateResult axisymmetric_film_state_update (
  const AxisymmetricFilmStateConfig * config,
  AxisymmetricFilmState * state,
  const AxisymmetricFilmObservation * observation, double dt)
{
  AxisymmetricFilmStateResult result = {
    .status = AXISYMMETRIC_FILM_STATE_INVALID_CONFIG,
    .nonpositive_bin = -1,
    .nonpositive_thickness = NAN,
    .minimum_state_thickness = NAN,
    .minimum_raw_gap = NAN,
    .maximum_state_raw_gap_discrepancy = NAN,
    .raw_inventory = NAN,
    .state_raw_inventory_discrepancy = NAN,
    .state_raw_inventory_tolerance = NAN
  };
  if (!axisymmetric_film_state_config_is_valid(config))
    return result;
  if (!state || !observation || !isfinite(dt) || dt <= 0.) {
    result.status = AXISYMMETRIC_FILM_STATE_INVALID_OBSERVATION;
    return result;
  }

  if (state->active) {
    if (!axisymmetric_film_state_active_ledger_is_valid(state)) {
      result.status = AXISYMMETRIC_FILM_STATE_RESTART_MISMATCH;
      return result;
    }
    if (!axisymmetric_film_state_active_state_is_valid(state)) {
      result.status = AXISYMMETRIC_FILM_STATE_NONPOSITIVE_STATE;
      return result;
    }
    axisymmetric_film_state_set_raw_diagnostics(&result, state, observation);
    if (!axisymmetric_film_state_active_observation_is_valid(
          config, observation, dt)) {
      result.status = AXISYMMETRIC_FILM_STATE_INVALID_OBSERVATION;
      return result;
    }

    const bool can_reacquire = axisymmetric_film_state_can_reacquire(config,
      state, observation, &result);
    if (can_reacquire && state->reacquisition_consecutive_observations ==
        AXISYMMETRIC_FILM_STATE_REACQUISITION_OBSERVATIONS - 1UL) {
      AxisymmetricFilmState candidate = *state;
      candidate.active = false;
      candidate.reacquisition_consecutive_observations =
        AXISYMMETRIC_FILM_STATE_REACQUISITION_OBSERVATIONS;
      *state = candidate;
      result.valid = true;
      result.active = false;
      result.reacquired = true;
      result.status = AXISYMMETRIC_FILM_STATE_REACQUIRED;
      result.inventory_before = state->inventory;
      result.inventory_after = state->inventory;
      result.expected_inventory_after = state->inventory;
      result.minimum_state_thickness =
        axisymmetric_film_state_minimum_thickness(state);
      result.reacquisition_consecutive_observations =
        candidate.reacquisition_consecutive_observations;
      return result;
    }
    if (state->update_count >= ULONG_MAX - 1UL ||
        (result.raw_geometry_unresolved &&
         state->cumulative_unresolved_observations >= ULONG_MAX - 1UL) ||
        (result.raw_crossed && state->cumulative_raw_crossing_observations >=
          ULONG_MAX - 1UL)) {
      result.status = AXISYMMETRIC_FILM_STATE_RESTART_MISMATCH;
      return result;
    }

    AxisymmetricFilmInput input = axisymmetric_film_state_reynolds_input(
      state, observation);
    const double inventory_before = axisymmetric_film_state_inventory(state);
    const double ledger_inventory = state->handoff_inventory -
      state->cumulative_edge_outflow;
    if (!axisymmetric_film_state_close(inventory_before, state->inventory) ||
        !axisymmetric_film_state_close(state->inventory, ledger_inventory)) {
      result.status = AXISYMMETRIC_FILM_STATE_INVENTORY_MISMATCH;
      result.inventory_before = inventory_before;
      result.inventory_after = state->inventory;
      result.expected_inventory_after = ledger_inventory;
      result.inventory_closure_error = fabs(state->inventory -
        ledger_inventory)/axisymmetric_film_state_scale(state->inventory,
          ledger_inventory);
      return result;
    }
    AxisymmetricReynoldsSolution solution = axisymmetric_reynolds_solution(
      &config->reynolds, &input, state->thickness);
    if (!solution.valid ||
        !isfinite(solution.maximum_flux_balance_error) ||
        solution.maximum_flux_balance_error > config->flux_balance_tolerance) {
      result.reynolds_solution = solution;
      result.status = AXISYMMETRIC_FILM_STATE_FLUX_SOLVE_FAILED;
      return result;
    }

    AxisymmetricFilmState candidate = *state;
    candidate.reacquisition_consecutive_observations = can_reacquire ?
      state->reacquisition_consecutive_observations + 1UL : 0UL;
    if (result.raw_geometry_unresolved)
      candidate.cumulative_unresolved_observations++;
    if (result.raw_crossed)
      candidate.cumulative_raw_crossing_observations++;
    for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++) {
      const double area = axisymmetric_film_state_ring_area(
        state->radius[bin], state->width[bin]);
      const double flux_difference = solution.face_flux[bin + 1] -
        solution.face_flux[bin];
      const double candidate_thickness = state->thickness[bin] -
        dt*flux_difference/area;
      if (!isfinite(candidate_thickness) || candidate_thickness <= 0.) {
        result.status = AXISYMMETRIC_FILM_STATE_NONPOSITIVE_STATE;
        result.nonpositive_bin = bin;
        result.nonpositive_thickness = candidate_thickness;
        result.inventory_before = inventory_before;
        result.reynolds_solution = solution;
        return result;
      }
      candidate.thickness[bin] = candidate_thickness;
      candidate.midpoint[bin] = state->midpoint[bin] + dt*.5*(
        observation->lower_normal_velocity[bin] +
        observation->upper_normal_velocity[bin]);
      candidate.face_flux[bin] = solution.face_flux[bin];
    }
    candidate.face_flux[AXISYMMETRIC_FILM_STATE_BINS] =
      solution.face_flux[AXISYMMETRIC_FILM_STATE_BINS];
    candidate.inventory = axisymmetric_film_state_inventory(&candidate);
    const double two_pi = 2.*acos(-1.);
    const double expected_inventory_after = inventory_before - two_pi*dt*
      solution.face_flux[AXISYMMETRIC_FILM_STATE_BINS];
    const double closure_scale = axisymmetric_film_state_max(1.e-30,
      axisymmetric_film_state_max(fabs(inventory_before),
        fabs(expected_inventory_after)));
    const double closure_error = fabs(candidate.inventory -
      expected_inventory_after)/closure_scale;
    if (!isfinite(candidate.inventory) || !isfinite(closure_error) ||
        closure_error > 1.e-12) {
      result.status = AXISYMMETRIC_FILM_STATE_INVENTORY_MISMATCH;
      result.inventory_before = inventory_before;
      result.inventory_after = candidate.inventory;
      result.expected_inventory_after = expected_inventory_after;
      result.inventory_closure_error = closure_error;
      result.reynolds_solution = solution;
      return result;
    }
    const double edge_increment = two_pi*dt*
      solution.face_flux[AXISYMMETRIC_FILM_STATE_BINS];
    candidate.cumulative_edge_outflow += edge_increment;
    if (!isfinite(edge_increment) ||
        !isfinite(candidate.cumulative_edge_outflow)) {
      result.status = AXISYMMETRIC_FILM_STATE_RESTART_MISMATCH;
      return result;
    }
    if (!axisymmetric_film_state_close(candidate.inventory,
          candidate.handoff_inventory - candidate.cumulative_edge_outflow)) {
      result.status = AXISYMMETRIC_FILM_STATE_INVENTORY_MISMATCH;
      return result;
    }
    candidate.update_count++;
    *state = candidate;
    result.valid = true;
    result.status = AXISYMMETRIC_FILM_STATE_ACTIVE;
    result.active = true;
    result.inventory_before = inventory_before;
    result.inventory_after = candidate.inventory;
    result.expected_inventory_after = expected_inventory_after;
    result.inventory_closure_error = closure_error;
    result.minimum_state_thickness =
      axisymmetric_film_state_minimum_thickness(&candidate);
    result.reacquisition_consecutive_observations =
      candidate.reacquisition_consecutive_observations;
    result.reynolds_solution = solution;
    return result;
  }

  if (!observation->raw_input_available) {
    result.status = AXISYMMETRIC_FILM_STATE_INVALID_OBSERVATION;
    result.observation_reason = AXISYMMETRIC_FILM_OBSERVATION_RAW_UNAVAILABLE;
    return result;
  }
  if (!axisymmetric_film_state_raw_observation_is_well_formed(
        config, observation, dt)) {
    result.status = AXISYMMETRIC_FILM_STATE_INVALID_OBSERVATION;
    result.observation_reason = AXISYMMETRIC_FILM_OBSERVATION_RAW_MALFORMED;
    return result;
  }

  bool enters_handoff = false;
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++)
    enters_handoff |= observation->raw.gap[bin] <
      AXISYMMETRIC_FILM_STATE_HANDOFF_CELLS*observation->raw.delta[bin];
  if (!enters_handoff) {
    result.valid = true;
    result.status = AXISYMMETRIC_FILM_STATE_INACTIVE;
    return result;
  }
  if (observation->raw_crossed_bins != 0 ||
      observation->raw_valid_area_fraction <
        config->minimum_raw_valid_area_fraction ||
      !observation->contact_corridor_at_maxlevel) {
    result.status = AXISYMMETRIC_FILM_STATE_INVALID_OBSERVATION;
    result.observation_reason = observation->raw_crossed_bins != 0 ?
      AXISYMMETRIC_FILM_OBSERVATION_CROSSED_ENTERING :
      observation->raw_valid_area_fraction <
        config->minimum_raw_valid_area_fraction ?
      AXISYMMETRIC_FILM_OBSERVATION_LOW_QUALITY_ENTERING :
      AXISYMMETRIC_FILM_OBSERVATION_UNDER_REFINED_ENTERING;
    return result;
  }
  if (state->activation_count == ULONG_MAX) {
    result.status = AXISYMMETRIC_FILM_STATE_RESTART_MISMATCH;
    return result;
  }

  AxisymmetricFilmState candidate = {0};
  candidate.active = true;
  candidate.activation_count = state->activation_count + 1;
  candidate.source_time = observation->time;
  candidate.update_count = 0;
  candidate.cumulative_edge_outflow = 0.;
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++) {
    candidate.thickness[bin] = observation->raw.gap[bin];
    candidate.midpoint[bin] = .5*(observation->raw.lower_position[bin] +
      observation->raw.upper_position[bin]);
    candidate.radius[bin] = observation->raw.radius[bin];
    candidate.width[bin] = observation->raw.width[bin];
    candidate.delta[bin] = observation->raw.delta[bin];
    candidate.face_flux[bin] = 0.;
  }
  candidate.face_flux[AXISYMMETRIC_FILM_STATE_BINS] = 0.;
  candidate.inventory = axisymmetric_film_state_inventory(&candidate);
  if (!isfinite(candidate.inventory) || candidate.inventory <= 0.) {
    result.status = AXISYMMETRIC_FILM_STATE_INVALID_OBSERVATION;
    return result;
  }
  candidate.handoff_inventory = candidate.inventory;
  if (!isfinite(candidate.source_time) ||
      !isfinite(candidate.handoff_inventory) ||
      !isfinite(candidate.inventory) ||
      !isfinite(candidate.cumulative_edge_outflow) ||
      !isfinite(candidate.cumulative_work)) {
    result.status = AXISYMMETRIC_FILM_STATE_INVALID_OBSERVATION;
    return result;
  }

  *state = candidate;
  result.valid = true;
  result.active = true;
  result.handed_off = true;
  result.status = AXISYMMETRIC_FILM_STATE_HANDOFF;
  result.inventory_before = candidate.inventory;
  result.inventory_after = candidate.inventory;
  result.expected_inventory_after = candidate.inventory;
  result.minimum_state_thickness =
    axisymmetric_film_state_minimum_thickness(&candidate);
  return result;
}

#endif
