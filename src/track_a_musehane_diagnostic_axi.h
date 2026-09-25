#ifndef TRACK_A_MUSEHANE_DIAGNOSTIC_AXI_H
#define TRACK_A_MUSEHANE_DIAGNOSTIC_AXI_H

#include <stdbool.h>

#include "axisymmetric_musehane_native_fields_axi.h"
#include "axisymmetric_musehane_observation_axi.h"
#include "film_model/musehane.h"

typedef enum {
  TRACK_A_MUSEHANE_DIAGNOSTIC_UNINITIALIZED = 0,
  TRACK_A_MUSEHANE_DIAGNOSTIC_READY = 1,
  TRACK_A_MUSEHANE_DIAGNOSTIC_OK = 2,
  TRACK_A_MUSEHANE_DIAGNOSTIC_INVALID_CONFIG = 3,
  TRACK_A_MUSEHANE_DIAGNOSTIC_INVALID_STATE = 4,
  TRACK_A_MUSEHANE_DIAGNOSTIC_NATIVE_FIELD_FAILED = 5,
  TRACK_A_MUSEHANE_DIAGNOSTIC_OBSERVATION_FAILED = 6,
  TRACK_A_MUSEHANE_DIAGNOSTIC_MODEL_FAILED = 7,
  TRACK_A_MUSEHANE_DIAGNOSTIC_NONZERO_DNS_TRACTION = 8
} TrackAMusehaneDiagnosticStatus;

typedef struct {
  TrackAMusehaneDiagnosticStatus status;
  FilmModelConfig config;
  FilmModelState model_state;
  bool model_initialized;
  unsigned long accepted_steps;
  AxisymmetricMusehaneNativeFieldStatus latest_native_status;
  AxisymmetricMusehaneNativeFieldAudit latest_native_audit;
  AxisymmetricMusehaneObservationStatus latest_observation_status;
  AxisymmetricMusehaneObservation latest_observation;
  FilmModelResult latest_result;
  MusehaneFilmDiagnostics latest_diagnostics;
} TrackAMusehaneDiagnosticState;

#define TRACK_A_MUSEHANE_DIAGNOSTIC_STATE_INITIALIZER { \
  .status = TRACK_A_MUSEHANE_DIAGNOSTIC_UNINITIALIZED, \
  .model_state = FILM_MODEL_STATE_INITIALIZER \
}

static inline TrackAMusehaneDiagnosticStatus
track_a_musehane_diagnostic_fail (
  TrackAMusehaneDiagnosticState *state,
  TrackAMusehaneDiagnosticStatus status)
{
  if (state)
    state->status = status;
  return status;
}

static inline TrackAMusehaneDiagnosticStatus
track_a_musehane_diagnostic_initialize (
  TrackAMusehaneDiagnosticState *state,
  const FilmModelConfig *config)
{
  if (!state || !config ||
      state->status != TRACK_A_MUSEHANE_DIAGNOSTIC_UNINITIALIZED ||
      state->model_state.guard != 0U || state->model_state.implementation ||
      config->model_kind != FILM_MODEL_MUSEHANE_INVERSE ||
      config->contract_version != FILM_MODEL_CONTRACT_VERSION ||
      config->traction_policy != MUSEHANE_TRACTION_DIAGNOSTIC_ONLY)
    return track_a_musehane_diagnostic_fail(
      state, TRACK_A_MUSEHANE_DIAGNOSTIC_INVALID_CONFIG);
  state->config = *config;
  state->status = TRACK_A_MUSEHANE_DIAGNOSTIC_READY;
  return state->status;
}

static inline bool track_a_musehane_result_has_zero_dns_traction (
  const FilmModelResult *result)
{
  if (!result || result->ring_count == 0 || !result->tractions)
    return false;
  for (size_t ring = 0; ring < result->ring_count; ring++) {
    const FilmRingTraction traction = result->tractions[ring];
    if (!isfinite(traction.lower_normal_force) ||
        !isfinite(traction.upper_normal_force) ||
        !isfinite(traction.lower_tangential_force) ||
        !isfinite(traction.upper_tangential_force) ||
        traction.lower_normal_force != 0. ||
        traction.upper_normal_force != 0. ||
        traction.lower_tangential_force != 0. ||
        traction.upper_tangential_force != 0.)
      return false;
  }
  return true;
}

static inline TrackAMusehaneDiagnosticStatus
track_a_musehane_diagnostic_step (
  TrackAMusehaneDiagnosticState *state,
  scalar water_pressure,
  vector pressure_gradient,
  vector capillary_force_density,
  double support_half_width,
  double sample_time,
  double sample_dt,
  unsigned long sample_iteration,
  bool anchor_released,
  uint32_t lower_identity,
  uint32_t upper_identity)
{
  if (!state ||
      (state->status != TRACK_A_MUSEHANE_DIAGNOSTIC_READY &&
       state->status != TRACK_A_MUSEHANE_DIAGNOSTIC_OK))
    return track_a_musehane_diagnostic_fail(
      state, TRACK_A_MUSEHANE_DIAGNOSTIC_INVALID_STATE);

  state->latest_native_status =
    build_axisymmetric_musehane_native_fields(
      water_pressure, pressure_gradient, capillary_force_density,
      &state->latest_native_audit);
  if (state->latest_native_status !=
      AXISYMMETRIC_MUSEHANE_NATIVE_FIELD_OK)
    return track_a_musehane_diagnostic_fail(
      state, TRACK_A_MUSEHANE_DIAGNOSTIC_NATIVE_FIELD_FAILED);

  state->latest_observation_status =
    sample_axisymmetric_musehane_observation(
      support_half_width, water_pressure, pressure_gradient,
      capillary_force_density, sample_time, sample_dt, sample_iteration,
      anchor_released, lower_identity, upper_identity,
      NULL,
      &state->latest_observation);
  if (state->latest_observation_status !=
      AXISYMMETRIC_MUSEHANE_OBSERVATION_OK)
    return track_a_musehane_diagnostic_fail(
      state, TRACK_A_MUSEHANE_DIAGNOSTIC_OBSERVATION_FAILED);

  if (!state->model_initialized) {
    const FilmModelStatus initialize_status = film_model_initialize(
      &state->model_state, &state->config,
      &state->latest_observation.view);
    if (initialize_status != FILM_MODEL_OK)
      return track_a_musehane_diagnostic_fail(
        state, TRACK_A_MUSEHANE_DIAGNOSTIC_MODEL_FAILED);
    state->model_initialized = true;
  }

  state->latest_result = film_model_advance(
    &state->model_state, &state->config,
    &state->latest_observation.view, sample_dt);
  if (state->latest_result.status != FILM_MODEL_OK ||
      !state->latest_result.commit_allowed ||
      state->latest_result.ring_count !=
        state->latest_observation.view.ring_count)
    return track_a_musehane_diagnostic_fail(
      state, TRACK_A_MUSEHANE_DIAGNOSTIC_MODEL_FAILED);
  if (!track_a_musehane_result_has_zero_dns_traction(
        &state->latest_result))
    return track_a_musehane_diagnostic_fail(
      state, TRACK_A_MUSEHANE_DIAGNOSTIC_NONZERO_DNS_TRACTION);
  if (musehane_film_model_diagnostics(
        &state->model_state, &state->latest_diagnostics) != FILM_MODEL_OK ||
      !state->latest_diagnostics.valid)
    return track_a_musehane_diagnostic_fail(
      state, TRACK_A_MUSEHANE_DIAGNOSTIC_MODEL_FAILED);

  state->accepted_steps++;
  state->status = TRACK_A_MUSEHANE_DIAGNOSTIC_OK;
  return state->status;
}

static inline void track_a_musehane_diagnostic_destroy (
  TrackAMusehaneDiagnosticState *state)
{
  if (!state)
    return;
  if (state->model_initialized)
    film_model_destroy(&state->model_state);
  *state = (TrackAMusehaneDiagnosticState)
    TRACK_A_MUSEHANE_DIAGNOSTIC_STATE_INITIALIZER;
}

#endif
