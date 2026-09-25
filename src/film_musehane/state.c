#include "film_musehane/state.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Retry the complete parent step on each refinement level. Local-only
   splitting is path-dependent for the nonlinear backward-Euler update. */
#define MUSEHANE_STATE_MAX_SUBSTEPS 1024U

struct MusehaneState {
  size_t cell_count;
  uint64_t mesh_signature;
  double time;
  unsigned long update_count;
  double cumulative_poiseuille_dissipation;
  double *estimated_thickness;
  double *thickness;
  unsigned char *active;
};

static int mesh_is_valid (const MusehaneSurfaceMesh *mesh)
{
  if (!mesh || mesh->cell_count == 0 || !mesh->cell_radius ||
      !mesh->cell_width || !mesh->cell_measure || !mesh->face_measure)
    return 0;
  for (size_t face = 0; face <= mesh->cell_count; face++)
    if (!isfinite(mesh->face_measure[face]) ||
        mesh->face_measure[face] < 0. ||
        (face > 0 && mesh->face_measure[face] <
          mesh->face_measure[face - 1]))
      return 0;
  for (size_t cell = 0; cell < mesh->cell_count; cell++)
    if (!isfinite(mesh->cell_radius[cell]) ||
        mesh->cell_radius[cell] < 0. ||
        !isfinite(mesh->cell_width[cell]) ||
        mesh->cell_width[cell] <= 0. ||
        !isfinite(mesh->cell_measure[cell]) ||
        mesh->cell_measure[cell] <= 0.)
      return 0;
  return 1;
}

static uint64_t hash_bytes (uint64_t hash, const void *value, size_t bytes)
{
  const unsigned char *data = value;
  for (size_t index = 0; index < bytes; index++) {
    hash ^= data[index];
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

static uint64_t mesh_signature (const MusehaneSurfaceMesh *mesh)
{
  uint64_t hash = UINT64_C(1469598103934665603);
  hash = hash_bytes(hash, &mesh->cell_count, sizeof mesh->cell_count);
  for (size_t cell = 0; cell < mesh->cell_count; cell++) {
    hash = hash_bytes(hash, &mesh->cell_radius[cell], sizeof(double));
    hash = hash_bytes(hash, &mesh->cell_width[cell], sizeof(double));
    hash = hash_bytes(hash, &mesh->cell_measure[cell], sizeof(double));
  }
  for (size_t face = 0; face <= mesh->cell_count; face++)
    hash = hash_bytes(hash, &mesh->face_measure[face], sizeof(double));
  return hash;
}

static int thickness_is_valid (const double *thickness, size_t cells)
{
  if (!thickness)
    return 0;
  for (size_t cell = 0; cell < cells; cell++)
    if (!isfinite(thickness[cell]) || thickness[cell] <= 0.)
      return 0;
  return 1;
}

static MusehaneStatus advance_core_adaptively (
  const MusehaneConfig *config, const MusehaneStepInput *parent,
  unsigned char *working_active, double *working_old, double *next,
  double *flux, const double *resolved_separation_threshold,
  MusehaneStepOutput *output, MusehaneCoreFailure *failure,
  unsigned *accepted_substeps, unsigned *resolved_release_count,
  double *resolved_release_inventory_change,
  double *dissipation_integral)
{
  const size_t cells = parent->mesh->cell_count;
  MusehaneStatus status = MUSEHANE_NONLINEAR_NOT_CONVERGED;
  *failure = (MusehaneCoreFailure) {0};
  *accepted_substeps = 0;
  *resolved_release_count = 0;
  *resolved_release_inventory_change = 0.;
  *dissipation_integral = 0.;

  for (unsigned count = 1; count <= MUSEHANE_STATE_MAX_SUBSTEPS;
       count *= 2) {
    memcpy(working_old, parent->old_thickness, cells*sizeof *working_old);
    memcpy(working_active, parent->active, cells*sizeof *working_active);
    double pressure_change = 0.;
    double kinematic_change = 0.;
    double release_inventory_change = 0.;
    double dissipation = 0.;
    double last_dissipation_rate = 0.;
    double inventory_before = 0.;
    double inventory_after = 0.;
    double minimum_thickness = INFINITY;
    unsigned total_iterations = 0;
    unsigned releases = 0;
    unsigned completed = 0;
    const double substep_dt = parent->dt/count;
    MusehaneCoreFailure last_failure = {0};

    for (unsigned substep = 0; substep < count; substep++) {
      MusehaneStepInput input = *parent;
      input.time = parent->time + substep*substep_dt;
      input.dt = substep_dt;
      input.active = working_active;
      input.old_thickness = working_old;
      MusehaneStepOutput candidate = {
        .thickness = next,
        .pressure_flux_integral_face = flux
      };
      last_failure = (MusehaneCoreFailure) {0};
      status = musehane_surface_step_diagnosed(config, &input,
        &candidate, &last_failure);
      if (status != MUSEHANE_OK)
        break;
      if (substep == 0)
        inventory_before = candidate.inventory_before;
      inventory_after = candidate.inventory_after;
      pressure_change += candidate.pressure_inventory_change;
      kinematic_change += candidate.kinematic_inventory_change;
      dissipation += substep_dt*candidate.poiseuille_dissipation_rate;
      last_dissipation_rate = candidate.poiseuille_dissipation_rate;
      total_iterations += candidate.iterations;
      double substep_release_inventory_change = 0.;
      if (resolved_separation_threshold)
        for (size_t cell = 0; cell < cells; cell++)
          if (working_active[cell] &&
              next[cell] >= resolved_separation_threshold[cell]) {
            substep_release_inventory_change -=
              parent->mesh->cell_measure[cell]*next[cell];
            working_active[cell] = 0;
            next[cell] = parent->estimated_thickness[cell];
            releases++;
          }
      release_inventory_change += substep_release_inventory_change;
      inventory_after = candidate.inventory_after +
        substep_release_inventory_change;
      minimum_thickness = INFINITY;
      int has_active = 0;
      for (size_t cell = 0; cell < cells; cell++) {
        if (!working_active[cell])
          continue;
        has_active = 1;
        minimum_thickness = fmin(minimum_thickness, next[cell]);
      }
      if (!has_active)
        for (size_t cell = 0; cell < cells; cell++)
          minimum_thickness = fmin(minimum_thickness,
            parent->estimated_thickness[cell]);
      memcpy(working_old, next, cells*sizeof *working_old);
      completed++;
      if (!has_active) {
        completed = count;
        break;
      }
    }

    if (completed == count) {
      const double expected_after = inventory_before + pressure_change +
        kinematic_change + release_inventory_change;
      const double closure_error = inventory_after - expected_after;
      const double inventory_scale = fmax(DBL_MIN,
        fmax(fabs(inventory_before), fabs(expected_after)));
      if (!isfinite(dissipation) || dissipation < 0. ||
          !isfinite(closure_error) || fabs(closure_error) >
            config->inventory_relative_tolerance*inventory_scale) {
        *failure = (MusehaneCoreFailure) {
          .inventory_diagnostic_valid = true,
          .inventory_before = inventory_before,
          .inventory_after = inventory_after,
          .expected_inventory_after = expected_after,
          .inventory_closure_error = closure_error,
          .inventory_relative_error = fabs(closure_error)/inventory_scale,
          .inventory_relative_tolerance =
            config->inventory_relative_tolerance
        };
        return MUSEHANE_INVENTORY_MISMATCH;
      }
      *output = (MusehaneStepOutput) {
        .thickness = next,
        .pressure_flux_integral_face = flux,
        .converged = true,
        .iterations = total_iterations,
        .minimum_thickness = minimum_thickness,
        .inventory_before = inventory_before,
        .inventory_after = inventory_after,
        .expected_inventory_after = expected_after,
        .pressure_inventory_change = pressure_change,
        .kinematic_inventory_change = kinematic_change,
        .inventory_closure_error = closure_error,
        .poiseuille_dissipation_rate = count == 1 ?
          last_dissipation_rate : dissipation/parent->dt
      };
      *accepted_substeps = count;
      *resolved_release_count = releases;
      *resolved_release_inventory_change = release_inventory_change;
      *dissipation_integral = dissipation;
      return MUSEHANE_OK;
    }

    *failure = last_failure;
    if (status != MUSEHANE_NONLINEAR_NOT_CONVERGED ||
        !last_failure.valid ||
        last_failure.kind != MUSEHANE_CORE_FAILURE_LINE_SEARCH_EXHAUSTED ||
        count == MUSEHANE_STATE_MAX_SUBSTEPS)
      return status;
  }
  return status;
}

static MusehaneState *allocate_state (size_t cells)
{
  if (cells == 0 ||
      cells > SIZE_MAX/(2*sizeof(double) + sizeof(unsigned char)))
    return NULL;
  MusehaneState *state = calloc(1, sizeof *state);
  if (!state)
    return NULL;
  state->estimated_thickness = malloc(cells*sizeof(double));
  state->thickness = malloc(cells*sizeof(double));
  state->active = calloc(cells, sizeof(unsigned char));
  if (!state->estimated_thickness || !state->thickness || !state->active) {
    musehane_state_destroy(state);
    return NULL;
  }
  state->cell_count = cells;
  return state;
}

MusehaneStateStatus musehane_state_create (
  const MusehaneSurfaceMesh *mesh,
  const double *initial_thickness,
  double initial_time,
  MusehaneState **output)
{
  if (!output || !mesh_is_valid(mesh) ||
      !thickness_is_valid(initial_thickness, mesh->cell_count) ||
      !isfinite(initial_time) || initial_time < 0.)
    return MUSEHANE_STATE_INVALID_INPUT;
  MusehaneState *candidate = allocate_state(mesh->cell_count);
  if (!candidate)
    return MUSEHANE_STATE_ALLOCATION_FAILED;
  memcpy(candidate->estimated_thickness, initial_thickness,
    mesh->cell_count*sizeof(double));
  memcpy(candidate->thickness, initial_thickness,
    mesh->cell_count*sizeof(double));
  candidate->mesh_signature = mesh_signature(mesh);
  candidate->time = initial_time;
  *output = candidate;
  return MUSEHANE_STATE_OK;
}

void musehane_state_destroy (MusehaneState *state)
{
  if (!state)
    return;
  free(state->estimated_thickness);
  free(state->thickness);
  free(state->active);
  free(state);
}

MusehaneStateStatus musehane_state_advance (
  const MusehaneConfig *config,
  MusehaneState *state,
  const MusehaneStateAdvanceInput *input,
  MusehaneStateAdvanceResult *output)
{
  if (!state || !input || !output || !mesh_is_valid(input->mesh) ||
      input->mesh->cell_count != state->cell_count ||
      mesh_signature(input->mesh) != state->mesh_signature ||
      !input->active || !input->estimated_thickness ||
      !thickness_is_valid(input->estimated_thickness,
        input->mesh->cell_count) ||
      !input->mean_surface_velocity_divergence ||
      !input->decontaminated_pressure_gradient_face ||
      !isfinite(input->dt) || input->dt <= 0.)
    return MUSEHANE_STATE_INVALID_INPUT;
  if (state->update_count == ULONG_MAX)
    return MUSEHANE_STATE_COUNTER_OVERFLOW;

  const size_t cells = state->cell_count;
  int has_active = 0;
  double minimum_estimated = INFINITY;
  for (size_t cell = 0; cell < cells; cell++) {
    if (input->active[cell] > 1 ||
        !isfinite(input->mean_surface_velocity_divergence[cell]))
      return MUSEHANE_STATE_INVALID_INPUT;
    if (input->resolved_separation_threshold &&
        (!isfinite(input->resolved_separation_threshold[cell]) ||
         input->resolved_separation_threshold[cell] <= 0.))
      return MUSEHANE_STATE_INVALID_INPUT;
    has_active |= input->active[cell] != 0;
    minimum_estimated = fmin(minimum_estimated,
      input->estimated_thickness[cell]);
  }
  for (size_t face = 0; face <= cells; face++)
    if (!isfinite(input->decontaminated_pressure_gradient_face[face]))
      return MUSEHANE_STATE_INVALID_INPUT;

  const double time_after = state->time + input->dt;
  if (!isfinite(time_after) || time_after < state->time)
    return MUSEHANE_STATE_CORE_REJECTED;
  if (!has_active) {
    MusehaneStateAdvanceResult candidate = {
      .valid = true,
      .core_status = MUSEHANE_NO_ACTIVE_REGION,
      .time_before = state->time,
      .time_after = time_after,
      .update_count = state->update_count + 1UL,
      .substep_count = 0,
      .minimum_thickness = minimum_estimated,
      .cumulative_poiseuille_dissipation =
        state->cumulative_poiseuille_dissipation
    };
    memcpy(state->estimated_thickness, input->estimated_thickness,
      cells*sizeof(double));
    memcpy(state->thickness, input->estimated_thickness,
      cells*sizeof(double));
    memcpy(state->active, input->active, cells*sizeof(unsigned char));
    state->time = time_after;
    state->update_count = candidate.update_count;
    *output = candidate;
    return MUSEHANE_STATE_OK;
  }

  double *old_for_step = malloc(cells*sizeof(double));
  unsigned char *working_active = malloc(cells*sizeof(unsigned char));
  double *working_old = malloc(cells*sizeof(double));
  double *next = malloc(cells*sizeof(double));
  double *flux = malloc((cells + 1)*sizeof(double));
  if (!old_for_step || !working_active || !working_old || !next || !flux) {
    free(old_for_step);
    free(working_active);
    free(working_old);
    free(next);
    free(flux);
    return MUSEHANE_STATE_ALLOCATION_FAILED;
  }
  for (size_t cell = 0; cell < cells; cell++) {
    old_for_step[cell] = input->active[cell] && state->active[cell] ?
      state->thickness[cell] : input->estimated_thickness[cell];
  }
  MusehaneStepInput core_input = {
    .time = state->time,
    .dt = input->dt,
    .mesh = input->mesh,
    .active = input->active,
    .estimated_thickness = input->estimated_thickness,
    .old_thickness = old_for_step,
    .mean_surface_velocity_divergence =
      input->mean_surface_velocity_divergence,
    .decontaminated_pressure_gradient_cell =
      input->decontaminated_pressure_gradient_cell,
    .decontaminated_pressure_gradient_face =
      input->decontaminated_pressure_gradient_face
  };
  MusehaneStepOutput core_output = {
    .thickness = next,
    .pressure_flux_integral_face = flux
  };
  MusehaneCoreFailure core_failure = {0};
  unsigned substep_count = 0;
  unsigned resolved_release_count = 0;
  double resolved_release_inventory_change = 0.;
  double dissipation_integral = 0.;
  const MusehaneStatus core_status = advance_core_adaptively(config,
    &core_input, working_active, working_old, next, flux,
    input->resolved_separation_threshold, &core_output, &core_failure,
    &substep_count, &resolved_release_count,
    &resolved_release_inventory_change, &dissipation_integral);
  if (core_status != MUSEHANE_OK) {
    free(old_for_step);
    free(working_active);
    free(working_old);
    free(next);
    free(flux);
    if (core_status == MUSEHANE_INVALID_CONFIG)
      return MUSEHANE_STATE_INVALID_CONFIG;
    if (core_status == MUSEHANE_ALLOCATION_FAILED)
      return MUSEHANE_STATE_ALLOCATION_FAILED;
    *output = (MusehaneStateAdvanceResult) {
      .valid = false,
      .core_status = core_status,
      .core_failure = core_failure,
      .time_before = state->time,
      .time_after = time_after,
      .update_count = state->update_count
    };
    return MUSEHANE_STATE_CORE_REJECTED;
  }

  const double cumulative_dissipation =
    state->cumulative_poiseuille_dissipation + dissipation_integral;
  if (!isfinite(time_after) || time_after < state->time ||
      !isfinite(cumulative_dissipation) || cumulative_dissipation < 0.) {
    free(old_for_step);
    free(working_active);
    free(working_old);
    free(next);
    free(flux);
    return MUSEHANE_STATE_CORE_REJECTED;
  }

  MusehaneStateAdvanceResult candidate = {
    .valid = true,
    .core_status = core_status,
    .time_before = state->time,
    .time_after = time_after,
    .update_count = state->update_count + 1UL,
    .substep_count = substep_count,
    .resolved_release_count = resolved_release_count,
    .minimum_thickness = core_output.minimum_thickness,
    .inventory_before = core_output.inventory_before,
    .inventory_after = core_output.inventory_after,
    .inventory_closure_error = core_output.inventory_closure_error,
    .pressure_inventory_change = core_output.pressure_inventory_change,
    .kinematic_inventory_change = core_output.kinematic_inventory_change,
    .resolved_release_inventory_change =
      resolved_release_inventory_change,
    .poiseuille_dissipation_rate =
      core_output.poiseuille_dissipation_rate,
    .cumulative_poiseuille_dissipation = cumulative_dissipation
  };
  memcpy(state->estimated_thickness, input->estimated_thickness,
    cells*sizeof(double));
  memcpy(state->thickness, next, cells*sizeof(double));
  memcpy(state->active, working_active, cells*sizeof(unsigned char));
  state->time = time_after;
  state->update_count = candidate.update_count;
  state->cumulative_poiseuille_dissipation = cumulative_dissipation;
  *output = candidate;
  free(old_for_step);
  free(working_active);
  free(working_old);
  free(next);
  free(flux);
  return MUSEHANE_STATE_OK;
}

MusehaneStateStatus musehane_state_snapshot (
  const MusehaneState *state,
  MusehaneStateSnapshot *output)
{
  if (!state || !output || output->capacity < state->cell_count ||
      !output->estimated_thickness || !output->thickness ||
      !output->active ||
      output->estimated_thickness == output->thickness ||
      !thickness_is_valid(state->estimated_thickness, state->cell_count) ||
      !thickness_is_valid(state->thickness, state->cell_count) ||
      !isfinite(state->time) || state->time < 0. ||
      !isfinite(state->cumulative_poiseuille_dissipation) ||
      state->cumulative_poiseuille_dissipation < 0.)
    return MUSEHANE_STATE_INVALID_INPUT;
  double *const estimated_buffer = output->estimated_thickness;
  double *const thickness_buffer = output->thickness;
  unsigned char *const active_buffer = output->active;
  const size_t capacity = output->capacity;
  memcpy(estimated_buffer, state->estimated_thickness,
    state->cell_count*sizeof(double));
  memcpy(thickness_buffer, state->thickness,
    state->cell_count*sizeof(double));
  memcpy(active_buffer, state->active,
    state->cell_count*sizeof(unsigned char));
  *output = (MusehaneStateSnapshot) {
    .schema_version = MUSEHANE_STATE_SNAPSHOT_SCHEMA,
    .cell_count = state->cell_count,
    .capacity = capacity,
    .mesh_signature = state->mesh_signature,
    .time = state->time,
    .update_count = state->update_count,
    .cumulative_poiseuille_dissipation =
      state->cumulative_poiseuille_dissipation,
    .estimated_thickness = estimated_buffer,
    .thickness = thickness_buffer,
    .active = active_buffer
  };
  return MUSEHANE_STATE_OK;
}

MusehaneStateStatus musehane_state_restore (
  const MusehaneSurfaceMesh *mesh,
  const MusehaneStateSnapshot *snapshot,
  MusehaneState **output)
{
  if (!output || !mesh_is_valid(mesh) || !snapshot ||
      snapshot->schema_version != MUSEHANE_STATE_SNAPSHOT_SCHEMA ||
      snapshot->cell_count != mesh->cell_count ||
      snapshot->capacity < snapshot->cell_count ||
      snapshot->mesh_signature != mesh_signature(mesh) ||
      !thickness_is_valid(snapshot->estimated_thickness,
        snapshot->cell_count) ||
      !thickness_is_valid(snapshot->thickness, snapshot->cell_count) ||
      !snapshot->active ||
      !isfinite(snapshot->time) || snapshot->time < 0. ||
      !isfinite(snapshot->cumulative_poiseuille_dissipation) ||
      snapshot->cumulative_poiseuille_dissipation < 0.)
    return MUSEHANE_STATE_RESTART_MISMATCH;
  for (size_t cell = 0; cell < snapshot->cell_count; cell++)
    if (snapshot->active[cell] > 1)
      return MUSEHANE_STATE_RESTART_MISMATCH;
  MusehaneState *candidate = allocate_state(snapshot->cell_count);
  if (!candidate)
    return MUSEHANE_STATE_ALLOCATION_FAILED;
  memcpy(candidate->estimated_thickness, snapshot->estimated_thickness,
    snapshot->cell_count*sizeof(double));
  memcpy(candidate->thickness, snapshot->thickness,
    snapshot->cell_count*sizeof(double));
  memcpy(candidate->active, snapshot->active,
    snapshot->cell_count*sizeof(unsigned char));
  candidate->mesh_signature = snapshot->mesh_signature;
  candidate->time = snapshot->time;
  candidate->update_count = snapshot->update_count;
  candidate->cumulative_poiseuille_dissipation =
    snapshot->cumulative_poiseuille_dissipation;
  *output = candidate;
  return MUSEHANE_STATE_OK;
}
