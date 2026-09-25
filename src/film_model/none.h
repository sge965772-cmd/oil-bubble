#ifndef FILM_MODEL_NONE_H
#define FILM_MODEL_NONE_H

#include <string.h>

#include "contract.h"

typedef struct FilmModelConfig {
  FilmModelKind model_kind;
  uint32_t contract_version;
} FilmModelConfig;

typedef struct FilmModelState {
  FilmModelKind model_kind;
  uint32_t contract_version;
  unsigned long completed_steps;
} FilmModelState;

static inline FilmModelConfig film_model_none_config (void)
{
  return (FilmModelConfig) {
    .model_kind = FILM_MODEL_NONE,
    .contract_version = FILM_MODEL_CONTRACT_VERSION
  };
}

static inline FilmModelStatus film_model_initialize (
  FilmModelState *state, const FilmModelConfig *config, const FilmDNSView *dns)
{
  if (!state || !config || config->model_kind != FILM_MODEL_NONE ||
      config->contract_version != FILM_MODEL_CONTRACT_VERSION)
    return FILM_MODEL_INVALID_CONFIG;
  if (!film_dns_view_is_well_formed(dns))
    return FILM_MODEL_INVALID_DNS_VIEW;
  *state = (FilmModelState) {
    .model_kind = FILM_MODEL_NONE,
    .contract_version = FILM_MODEL_CONTRACT_VERSION,
    .completed_steps = 0
  };
  return FILM_MODEL_OK;
}

static inline FilmModelResult film_model_advance (
  FilmModelState *state, const FilmModelConfig *config, const FilmDNSView *dns,
  double dt)
{
  if (!config || config->model_kind != FILM_MODEL_NONE ||
      config->contract_version != FILM_MODEL_CONTRACT_VERSION)
    return (FilmModelResult) {
      .status = FILM_MODEL_INVALID_CONFIG,
      .model_kind = FILM_MODEL_NONE,
      .reason = "invalid no-subgrid config"
    };
  if (!state || state->model_kind != FILM_MODEL_NONE ||
      state->contract_version != FILM_MODEL_CONTRACT_VERSION)
    return (FilmModelResult) {
      .status = FILM_MODEL_INVALID_STATE,
      .model_kind = FILM_MODEL_NONE,
      .reason = "invalid no-subgrid state"
    };
  if (!dns || !isfinite(dt) || dt != dns->dt ||
      !film_dns_view_is_well_formed(dns))
    return (FilmModelResult) {
      .status = FILM_MODEL_INVALID_DNS_VIEW,
      .model_kind = FILM_MODEL_NONE,
      .reason = "invalid no-subgrid DNS view"
    };
  state->completed_steps++;
  return (FilmModelResult) {
    .status = FILM_MODEL_OK,
    .model_kind = FILM_MODEL_NONE,
    .active = false,
    .commit_allowed = true,
    .ring_count = 0,
    .tractions = NULL,
    .minimum_thickness = 0.,
    .inventory = 0.,
    .paired_force_residual = 0.,
    .paired_moment_residual = 0.,
    .step_work = 0.,
    .cumulative_work = 0.,
    .completed_stage = FILM_STAGE_EVALUATE_TRACTION_N,
    .reason = "no subgrid film correction"
  };
}

static inline FilmModelStatus film_model_apply_traction (
  const FilmModelResult *result, FilmDNSTractionTarget *target)
{
  if (!result || !target || result->status != FILM_MODEL_OK ||
      result->model_kind != FILM_MODEL_NONE || result->active ||
      result->ring_count != 0 || result->tractions != NULL)
    return FILM_MODEL_TRACTION_REJECTED;
  return FILM_MODEL_OK;
}

static inline FilmModelStatus film_model_checkpoint (
  const FilmModelState *state, FilmRestartRecord *record)
{
  if (!state || !record || !record->data ||
      record->capacity < sizeof *state || state->model_kind != FILM_MODEL_NONE)
    return FILM_MODEL_RESTART_MISMATCH;
  memcpy(record->data, state, sizeof *state);
  record->schema_version = FILM_MODEL_RESTART_SCHEMA_VERSION;
  record->model_kind = FILM_MODEL_NONE;
  record->size = sizeof *state;
  return FILM_MODEL_OK;
}

static inline FilmModelStatus film_model_restore (
  FilmModelState *state, const FilmModelConfig *config,
  const FilmRestartRecord *record)
{
  if (!state || !config || !record || !record->data ||
      config->model_kind != FILM_MODEL_NONE ||
      record->model_kind != FILM_MODEL_NONE ||
      record->schema_version != FILM_MODEL_RESTART_SCHEMA_VERSION ||
      record->size != sizeof *state || record->capacity < record->size)
    return FILM_MODEL_RESTART_MISMATCH;
  FilmModelState candidate;
  memcpy(&candidate, record->data, sizeof candidate);
  if (candidate.model_kind != FILM_MODEL_NONE ||
      candidate.contract_version != FILM_MODEL_CONTRACT_VERSION)
    return FILM_MODEL_RESTART_MISMATCH;
  *state = candidate;
  return FILM_MODEL_OK;
}

#endif
