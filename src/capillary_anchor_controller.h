#ifndef CAPILLARY_ANCHOR_CONTROLLER_H
#define CAPILLARY_ANCHOR_CONTROLLER_H

#include <stdbool.h>
#include <math.h>

typedef int CapillaryAnchorStatus;

enum {
  CAPILLARY_ANCHOR_HOLDING = 0,
  CAPILLARY_ANCHOR_RELEASED = 1
};

typedef struct {
  double release_gap;
  double proportional_gain;
  double integral_gain;
  double acceleration_limit;
  double feedforward_acceleration;
  bool require_approach;
} CapillaryAnchorConfig;

typedef struct {
  CapillaryAnchorStatus status;
  double integral_error;
  double commanded_acceleration;
  double release_time;
  double release_gap;
  double last_reaction_force;
  double last_power;
  double cumulative_work;
  unsigned long saturation_count;
} CapillaryAnchorState;

typedef struct {
  double time;
  double dt;
  double outer_gap;
  double closing_speed;
  double upper_center_velocity;
} CapillaryAnchorObservation;

typedef struct {
  bool active;
  bool released_this_step;
  double acceleration;
} CapillaryAnchorCommand;

typedef struct {
  CapillaryAnchorStatus status;
  double integral_error;
  double commanded_acceleration;
  double release_time;
  double release_gap;
  double last_reaction_force;
  double last_power;
  double cumulative_work;
  unsigned long saturation_count;
} CapillaryAnchorRestartSnapshot;

static inline CapillaryAnchorConfig capillary_anchor_config (
  double release_gap,
  double proportional_gain,
  double integral_gain,
  double acceleration_limit,
  double feedforward_acceleration,
  bool require_approach)
{
  CapillaryAnchorConfig config = {
    .release_gap = release_gap,
    .proportional_gain = proportional_gain,
    .integral_gain = integral_gain,
    .acceleration_limit = acceleration_limit,
    .feedforward_acceleration = feedforward_acceleration,
    .require_approach = require_approach
  };
  return config;
}

static inline CapillaryAnchorState capillary_anchor_initial_state (void)
{
  CapillaryAnchorState state;
  state.status = CAPILLARY_ANCHOR_HOLDING;
  state.integral_error = 0.;
  state.commanded_acceleration = 0.;
  state.release_time = NAN;
  state.release_gap = NAN;
  state.last_reaction_force = 0.;
  state.last_power = 0.;
  state.cumulative_work = 0.;
  state.saturation_count = 0;
  return state;
}

static inline CapillaryAnchorState capillary_anchor_released_state (void)
{
  CapillaryAnchorState state = capillary_anchor_initial_state();
  state.status = CAPILLARY_ANCHOR_RELEASED;
  return state;
}

static inline bool capillary_anchor_restore (
  CapillaryAnchorState * state,
  CapillaryAnchorRestartSnapshot snapshot,
  double restart_time)
{
  if (!state || !isfinite(restart_time) || restart_time < 0. ||
      !isfinite(snapshot.integral_error) ||
      !isfinite(snapshot.commanded_acceleration) ||
      !isfinite(snapshot.last_reaction_force) ||
      !isfinite(snapshot.last_power) || !isfinite(snapshot.cumulative_work))
    return false;
  if (snapshot.status == CAPILLARY_ANCHOR_RELEASED &&
      (!isfinite(snapshot.release_time) ||
       !isfinite(snapshot.release_gap) || snapshot.release_gap <= 0. ||
       snapshot.release_time > restart_time))
    return false;
  if (snapshot.status != CAPILLARY_ANCHOR_HOLDING &&
      snapshot.status != CAPILLARY_ANCHOR_RELEASED)
    return false;
  state->status = snapshot.status;
  state->integral_error = snapshot.integral_error;
  state->commanded_acceleration = snapshot.commanded_acceleration;
  state->release_time = snapshot.status == CAPILLARY_ANCHOR_RELEASED ?
    snapshot.release_time : NAN;
  state->release_gap = snapshot.status == CAPILLARY_ANCHOR_RELEASED ?
    snapshot.release_gap : NAN;
  state->last_reaction_force = snapshot.last_reaction_force;
  state->last_power = snapshot.last_power;
  state->cumulative_work = snapshot.cumulative_work;
  state->saturation_count = snapshot.saturation_count;
  return true;
}

static inline double capillary_anchor_clamp (
  double value, double magnitude_limit)
{
  if (value > magnitude_limit)
    return magnitude_limit;
  if (value < -magnitude_limit)
    return -magnitude_limit;
  return value;
}

static inline CapillaryAnchorCommand capillary_anchor_update (
  CapillaryAnchorState * state,
  const CapillaryAnchorConfig * config,
  CapillaryAnchorObservation observation)
{
  CapillaryAnchorCommand command = {0};
  if (state->status == CAPILLARY_ANCHOR_RELEASED) {
    state->commanded_acceleration = 0.;
    return command;
  }

  bool release = observation.outer_gap <= config->release_gap &&
                 (!config->require_approach ||
                  observation.closing_speed > 0.);
  if (release) {
    state->status = CAPILLARY_ANCHOR_RELEASED;
    state->commanded_acceleration = 0.;
    state->release_time = observation.time;
    state->release_gap = observation.outer_gap;
    command.released_this_step = true;
    return command;
  }

  command.active = true;
  double error = -observation.upper_center_velocity;
  double candidate_integral = state->integral_error +
                              error*observation.dt;
  double unconstrained = config->feedforward_acceleration +
                         config->proportional_gain*error +
                         config->integral_gain*candidate_integral;
  double constrained = capillary_anchor_clamp(
    unconstrained, config->acceleration_limit);

  if (constrained == unconstrained)
    state->integral_error = candidate_integral;
  else
    state->saturation_count++;

  state->commanded_acceleration = constrained;
  command.acceleration = constrained;
  return command;
}

static inline void capillary_anchor_record_reaction (
  CapillaryAnchorState * state,
  double reaction_force,
  double upper_center_velocity,
  double dt)
{
  state->last_reaction_force = reaction_force;
  state->last_power = reaction_force*upper_center_velocity;
  state->cumulative_work += state->last_power*dt;
}

#endif
