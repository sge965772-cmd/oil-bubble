#include "film_model/musehane.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "film_musehane/activation.h"
#include "film_musehane/equation27_correction.h"
#include "film_musehane/matched_traction.h"
#include "film_musehane/state.h"
#include "film_musehane/traction.h"

#define MUSEHANE_MODEL_GUARD UINT64_C(0x4d55534548414e45)
#define MUSEHANE_RESTART_MAGIC UINT64_C(0x4d55534552535431)
#define TWO_PI 6.283185307179586476925286766559005768

typedef struct {
  FilmModelConfig config;
  size_t ring_count;
  MusehaneSurfaceMesh mesh;
  double *radius;
  double *width;
  double *measure;
  double *face_measure;
  MusehaneState *evolution;
  unsigned long completed_steps;
  unsigned long last_iteration;
  unsigned substep_count;
  unsigned resolved_release_count;
  bool diagnostics_valid;
  FilmRingTraction *tractions;
  unsigned char *active_n;
  unsigned char *active_np1;
  unsigned char *traction_active;
  double *estimated;
  double *thickness_n;
  double *thickness_np1;
  double *divergence;
  double *decontaminated_gradient;
  double *smoothed_gradient;
  double *lower_analytic_traction;
  double *upper_analytic_traction;
  double *lower_matched_correction;
  double *upper_matched_correction;
  double *resolved_reaction_force;
  size_t matched_corrected_ring_count;
  double minimum_matched_reference_thickness;
  double maximum_resolved_region_correction;
  double analytical_film_power;
  double resolved_reference_power;
  double matched_correction_power;
  double additional_relative_shear_dissipation_rate;
  double step_matched_correction_work;
  double cumulative_matched_correction_work;
  double equation27_pressure_correction_power;
  double equation27_shear_correction_power;
  double equation27_correction_power;
  double equation27_total_resolved_reaction_force;
  double maximum_equation27_pair_closure_residual;
  double time_n;
  double time_np1;
  double pressure_residual;
  double traction_residual;
  double pressure_inventory_change;
  double kinematic_inventory_change;
  double resolved_release_inventory_change;
  double dissipation_rate;
  double cumulative_dissipation;
  bool restart_clock_rebased;
  bool rejected_core_input_valid;
  double rejected_core_time;
  double rejected_core_dt;
  unsigned char *rejected_core_active;
  double *rejected_core_estimated;
  double *rejected_core_old;
  double *rejected_core_divergence;
  double *rejected_core_gradient_cell;
  double *rejected_core_gradient_face;
} MusehaneModelImplementation;

typedef struct {
  uint64_t magic;
  uint32_t schema_version;
  uint32_t contract_version;
  uint32_t model_kind;
  uint32_t traction_policy;
  uint64_t ring_count;
  uint64_t completed_steps;
  uint64_t last_iteration;
  uint64_t mesh_signature;
  double time;
  double cumulative_dissipation;
  double cumulative_matched_correction_work;
  double water_density;
  double water_viscosity;
  double nonlinear_tolerance;
  double inventory_tolerance;
  uint32_t maximum_iterations;
  uint32_t discretization;
  double smoothing_length_squared;
  double pressure_residual_tolerance;
  double activation_cells;
} MusehaneRestartHeader;

typedef struct {
  unsigned char *active_n;
  unsigned char *active_np1;
  unsigned char *traction_active;
  double *estimated_n;
  double *estimated;
  double *thickness_n;
  double *thickness_np1;
  double *divergence;
  double *raw_gradient;
  double *capillary_gradient;
  double *decontaminated_gradient;
  double *smoothed_gradient;
  double *smoothed_gradient_face;
  double *lower_velocity;
  double *upper_velocity;
  double *lower_analytic_traction;
  double *upper_analytic_traction;
  double *lower_matched_correction;
  double *upper_matched_correction;
  double *resolved_reaction_force;
  double *area_weight;
  FilmRingTraction *tractions;
  double *local_delta;
  double *resolved_separation_threshold;
  double *lower_identity_support;
  double *upper_identity_support;
  double *surface_velocity;
  double *surface_velocity_face;
} MusehaneStepWorkspace;

static FilmModelResult rejected_result (
  FilmModelStatus status, const char *reason)
{
  return (FilmModelResult) {
    .status = status,
    .model_kind = FILM_MODEL_MUSEHANE_INVERSE,
    .commit_allowed = false,
    .completed_stage = FILM_STAGE_OBSERVE_N,
    .reason = reason
  };
}

static FilmModelResult traction_rejected_result (
  const char *reason,
  int traction_status,
  size_t failure_cell,
  const MusehaneEquation27CorrectionOutput *equation27,
  double step_work,
  double cumulative_work)
{
  FilmModelResult result = rejected_result(
    FILM_MODEL_TRACTION_REJECTED, reason);
  result.failure_traction_diagnostic_valid = true;
  result.failure_traction_status = traction_status;
  result.failure_traction_cell = failure_cell;
  result.failure_traction_pressure_power = equation27 ?
    equation27->pressure_correction_power : NAN;
  result.failure_traction_shear_power = equation27 ?
    equation27->shear_correction_power : NAN;
  result.failure_traction_direct_power = equation27 ?
    equation27->direct_correction_power : NAN;
  result.failure_traction_total_power = equation27 ?
    equation27->total_correction_power : NAN;
  result.failure_traction_power_identity_residual = equation27 ?
    equation27->power_identity_residual : NAN;
  result.failure_traction_step_work = step_work;
  result.failure_traction_cumulative_work = cumulative_work;
  return result;
}

static int close_enough (double first, double second)
{
  const double scale = fmax(DBL_MIN, fmax(fabs(first), fabs(second)));
  return fabs(first - second) <= 512.*DBL_EPSILON*scale;
}

static int close_enough_at_scale (
  double first, double second, double reference_scale)
{
  const double scale = fmax(DBL_MIN, fmax(fabs(reference_scale),
    fmax(fabs(first), fabs(second))));
  return fabs(first - second) <= 512.*DBL_EPSILON*scale;
}

static int config_is_valid (const FilmModelConfig *config)
{
  return config && config->model_kind == FILM_MODEL_MUSEHANE_INVERSE &&
    config->contract_version == FILM_MODEL_CONTRACT_VERSION &&
    isfinite(config->equation.water_density) &&
    config->equation.water_density > 0. &&
    isfinite(config->equation.water_viscosity) &&
    config->equation.water_viscosity > 0. &&
    isfinite(config->equation.nonlinear_relative_tolerance) &&
    config->equation.nonlinear_relative_tolerance > 0. &&
    isfinite(config->equation.inventory_relative_tolerance) &&
    config->equation.inventory_relative_tolerance > 0. &&
    config->equation.maximum_iterations > 0 &&
    isfinite(config->pressure_driver.smoothing_length_squared) &&
    config->pressure_driver.smoothing_length_squared >= 0. &&
    isfinite(config->pressure_driver.linear_residual_relative_tolerance) &&
    config->pressure_driver.linear_residual_relative_tolerance > 0. &&
    (config->equation.discretization ==
       MUSEHANE_DISCRETIZATION_CONSERVATIVE_HARMONIC ||
     config->equation.discretization ==
       MUSEHANE_DISCRETIZATION_PUBLISHED_EQUATION23) &&
    isfinite(config->activation_cells) && config->activation_cells > 0. &&
    (config->traction_policy == MUSEHANE_TRACTION_DIAGNOSTIC_ONLY ||
     config->traction_policy == MUSEHANE_MATCHED_VISCOUS_CORRECTION ||
     config->traction_policy == MUSEHANE_MATCHED_EQUATION27_CORRECTION ||
     config->traction_policy == MUSEHANE_PUBLISHED_EQUATION27_TRACTION);
}

static int configs_match (
  const FilmModelConfig *first, const FilmModelConfig *second)
{
  return first && second && first->model_kind == second->model_kind &&
    first->contract_version == second->contract_version &&
    first->equation.water_density == second->equation.water_density &&
    first->equation.water_viscosity == second->equation.water_viscosity &&
    first->equation.nonlinear_relative_tolerance ==
      second->equation.nonlinear_relative_tolerance &&
    first->equation.inventory_relative_tolerance ==
      second->equation.inventory_relative_tolerance &&
    first->equation.maximum_iterations ==
      second->equation.maximum_iterations &&
    first->equation.discretization == second->equation.discretization &&
    first->pressure_driver.smoothing_length_squared ==
      second->pressure_driver.smoothing_length_squared &&
    first->pressure_driver.linear_residual_relative_tolerance ==
      second->pressure_driver.linear_residual_relative_tolerance &&
    first->activation_cells == second->activation_cells &&
    first->traction_policy == second->traction_policy;
}

FilmModelConfig film_model_musehane_config (void)
{
  return (FilmModelConfig) {
    .model_kind = FILM_MODEL_MUSEHANE_INVERSE,
    .contract_version = FILM_MODEL_CONTRACT_VERSION,
    .equation = {
      .water_density = 998.2,
      .water_viscosity = .89e-3,
      .discretization = MUSEHANE_DISCRETIZATION_CONSERVATIVE_HARMONIC,
      .nonlinear_relative_tolerance = 1.e-10,
      .inventory_relative_tolerance = 1.e-10,
      .maximum_iterations = 100
    },
    .pressure_driver = {
      .smoothing_length_squared = 0.,
      .linear_residual_relative_tolerance = 1.e-12
    },
    .activation_cells = 7.,
    .traction_policy = MUSEHANE_TRACTION_DIAGNOSTIC_ONLY
  };
}

static void free_implementation (MusehaneModelImplementation *implementation)
{
  if (!implementation)
    return;
  musehane_state_destroy(implementation->evolution);
  free(implementation->radius);
  free(implementation->width);
  free(implementation->measure);
  free(implementation->face_measure);
  free(implementation->tractions);
  free(implementation->active_n);
  free(implementation->active_np1);
  free(implementation->traction_active);
  free(implementation->estimated);
  free(implementation->thickness_n);
  free(implementation->thickness_np1);
  free(implementation->divergence);
  free(implementation->decontaminated_gradient);
  free(implementation->smoothed_gradient);
  free(implementation->lower_analytic_traction);
  free(implementation->upper_analytic_traction);
  free(implementation->lower_matched_correction);
  free(implementation->upper_matched_correction);
  free(implementation->resolved_reaction_force);
  free(implementation->rejected_core_active);
  free(implementation->rejected_core_estimated);
  free(implementation->rejected_core_old);
  free(implementation->rejected_core_divergence);
  free(implementation->rejected_core_gradient_cell);
  free(implementation->rejected_core_gradient_face);
  free(implementation);
}

static MusehaneModelImplementation *allocate_implementation (size_t cells)
{
  if (cells == 0 || cells > SIZE_MAX/sizeof(double) - 1)
    return NULL;
  MusehaneModelImplementation *result = calloc(1, sizeof *result);
  if (!result)
    return NULL;
#define ALLOCATE_MEMBER(member, count) \
  do { \
    result->member = calloc((count), sizeof *result->member); \
    if (!result->member) { \
      free_implementation(result); \
      return NULL; \
    } \
  } while (0)
  ALLOCATE_MEMBER(radius, cells);
  ALLOCATE_MEMBER(width, cells);
  ALLOCATE_MEMBER(measure, cells);
  ALLOCATE_MEMBER(face_measure, cells + 1);
  ALLOCATE_MEMBER(tractions, cells);
  ALLOCATE_MEMBER(active_n, cells);
  ALLOCATE_MEMBER(active_np1, cells);
  ALLOCATE_MEMBER(traction_active, cells);
  ALLOCATE_MEMBER(estimated, cells);
  ALLOCATE_MEMBER(thickness_n, cells);
  ALLOCATE_MEMBER(thickness_np1, cells);
  ALLOCATE_MEMBER(divergence, cells);
  ALLOCATE_MEMBER(decontaminated_gradient, cells);
  ALLOCATE_MEMBER(smoothed_gradient, cells);
  ALLOCATE_MEMBER(lower_analytic_traction, cells);
  ALLOCATE_MEMBER(upper_analytic_traction, cells);
  ALLOCATE_MEMBER(lower_matched_correction, cells);
  ALLOCATE_MEMBER(upper_matched_correction, cells);
  ALLOCATE_MEMBER(resolved_reaction_force, cells);
  ALLOCATE_MEMBER(rejected_core_active, cells);
  ALLOCATE_MEMBER(rejected_core_estimated, cells);
  ALLOCATE_MEMBER(rejected_core_old, cells);
  ALLOCATE_MEMBER(rejected_core_divergence, cells);
  ALLOCATE_MEMBER(rejected_core_gradient_cell, cells);
  ALLOCATE_MEMBER(rejected_core_gradient_face, cells + 1);
#undef ALLOCATE_MEMBER
  result->ring_count = cells;
  result->mesh = (MusehaneSurfaceMesh) {
    .cell_count = cells,
    .cell_radius = result->radius,
    .cell_width = result->width,
    .cell_measure = result->measure,
    .face_measure = result->face_measure
  };
  return result;
}

static int build_mesh (
  MusehaneModelImplementation *implementation, const FilmDNSView *dns)
{
  if (!implementation || !dns || dns->ring_count != implementation->ring_count)
    return 0;
  double previous_right = 0.;
  for (size_t cell = 0; cell < dns->ring_count; cell++) {
    const FilmRingObservation *ring = &dns->rings[cell];
    const double left = ring->radius - .5*ring->width;
    const double right = ring->radius + .5*ring->width;
    if (!isfinite(left) || !isfinite(right) || left < 0. || right <= left ||
        (cell == 0 && !close_enough_at_scale(left, 0., right)) ||
        (cell > 0 && !close_enough_at_scale(left, previous_right,
          fmax(ring->width, right))))
      return 0;
    implementation->radius[cell] = ring->radius;
    implementation->width[cell] = ring->width;
    implementation->measure[cell] = .5*(right*right - left*left);
    implementation->face_measure[cell] = left;
    previous_right = right;
  }
  implementation->face_measure[dns->ring_count] = previous_right;
  return 1;
}

static int mesh_matches_dns (
  const MusehaneModelImplementation *implementation, const FilmDNSView *dns)
{
  if (!implementation || !dns || dns->ring_count != implementation->ring_count)
    return 0;
  for (size_t cell = 0; cell < dns->ring_count; cell++)
    if (!close_enough(implementation->radius[cell], dns->rings[cell].radius) ||
        !close_enough(implementation->width[cell], dns->rings[cell].width))
      return 0;
  return 1;
}

FilmModelStatus film_model_initialize (
  FilmModelState *state, const FilmModelConfig *config,
  const FilmDNSView *dns)
{
  if (!config_is_valid(config))
    return FILM_MODEL_INVALID_CONFIG;
  if (!state || state->guard != 0U || state->implementation != NULL)
    return FILM_MODEL_INVALID_STATE;
  if (!film_dns_view_has_usable_geometry(dns))
    return FILM_MODEL_INVALID_DNS_VIEW;
  MusehaneModelImplementation *candidate =
    allocate_implementation(dns->ring_count);
  if (!candidate)
    return FILM_MODEL_INVALID_STATE;
  candidate->config = *config;
  candidate->last_iteration = dns->iteration;
  if (!build_mesh(candidate, dns)) {
    free_implementation(candidate);
    return FILM_MODEL_INVALID_DNS_VIEW;
  }
  for (size_t cell = 0; cell < dns->ring_count; cell++)
    candidate->estimated[cell] = dns->rings[cell].gap;
  if (musehane_state_create(&candidate->mesh, candidate->estimated,
      dns->time, &candidate->evolution) != MUSEHANE_STATE_OK) {
    free_implementation(candidate);
    return FILM_MODEL_INVALID_STATE;
  }
  state->guard = MUSEHANE_MODEL_GUARD;
  state->implementation = candidate;
  return FILM_MODEL_OK;
}

void film_model_destroy (FilmModelState *state)
{
  if (!state)
    return;
  if (state->guard == MUSEHANE_MODEL_GUARD)
    free_implementation(state->implementation);
  state->guard = 0U;
  state->implementation = NULL;
}

static MusehaneModelImplementation *implementation_of (FilmModelState *state)
{
  return state && state->guard == MUSEHANE_MODEL_GUARD ?
    state->implementation : NULL;
}

static const MusehaneModelImplementation *constant_implementation_of (
  const FilmModelState *state)
{
  return state && state->guard == MUSEHANE_MODEL_GUARD ?
    state->implementation : NULL;
}

static void free_workspace (MusehaneStepWorkspace *workspace)
{
  if (!workspace)
    return;
  free(workspace->active_n);
  free(workspace->active_np1);
  free(workspace->traction_active);
  free(workspace->estimated_n);
  free(workspace->estimated);
  free(workspace->thickness_n);
  free(workspace->thickness_np1);
  free(workspace->divergence);
  free(workspace->raw_gradient);
  free(workspace->capillary_gradient);
  free(workspace->decontaminated_gradient);
  free(workspace->smoothed_gradient);
  free(workspace->smoothed_gradient_face);
  free(workspace->lower_velocity);
  free(workspace->upper_velocity);
  free(workspace->lower_analytic_traction);
  free(workspace->upper_analytic_traction);
  free(workspace->lower_matched_correction);
  free(workspace->upper_matched_correction);
  free(workspace->resolved_reaction_force);
  free(workspace->area_weight);
  free(workspace->tractions);
  free(workspace->local_delta);
  free(workspace->resolved_separation_threshold);
  free(workspace->lower_identity_support);
  free(workspace->upper_identity_support);
  free(workspace->surface_velocity);
  free(workspace->surface_velocity_face);
  memset(workspace, 0, sizeof *workspace);
}

static int allocate_workspace (MusehaneStepWorkspace *workspace, size_t cells)
{
  if (!workspace || cells == 0 || cells > SIZE_MAX/sizeof(double) - 1)
    return 0;
  memset(workspace, 0, sizeof *workspace);
#define ALLOCATE_WORK(member, count) \
  do { \
    workspace->member = calloc((count), sizeof *workspace->member); \
    if (!workspace->member) { \
      free_workspace(workspace); \
      return 0; \
    } \
  } while (0)
  ALLOCATE_WORK(active_n, cells);
  ALLOCATE_WORK(active_np1, cells);
  ALLOCATE_WORK(traction_active, cells);
  ALLOCATE_WORK(estimated_n, cells);
  ALLOCATE_WORK(estimated, cells);
  ALLOCATE_WORK(thickness_n, cells);
  ALLOCATE_WORK(thickness_np1, cells);
  ALLOCATE_WORK(divergence, cells);
  ALLOCATE_WORK(raw_gradient, cells);
  ALLOCATE_WORK(capillary_gradient, cells);
  ALLOCATE_WORK(decontaminated_gradient, cells);
  ALLOCATE_WORK(smoothed_gradient, cells);
  ALLOCATE_WORK(smoothed_gradient_face, cells + 1);
  ALLOCATE_WORK(lower_velocity, cells);
  ALLOCATE_WORK(upper_velocity, cells);
  ALLOCATE_WORK(lower_analytic_traction, cells);
  ALLOCATE_WORK(upper_analytic_traction, cells);
  ALLOCATE_WORK(lower_matched_correction, cells);
  ALLOCATE_WORK(upper_matched_correction, cells);
  ALLOCATE_WORK(resolved_reaction_force, cells);
  ALLOCATE_WORK(area_weight, cells);
  ALLOCATE_WORK(tractions, cells);
  ALLOCATE_WORK(local_delta, cells);
  ALLOCATE_WORK(resolved_separation_threshold, cells);
  ALLOCATE_WORK(lower_identity_support, cells);
  ALLOCATE_WORK(upper_identity_support, cells);
  ALLOCATE_WORK(surface_velocity, cells);
  ALLOCATE_WORK(surface_velocity_face, cells + 1);
#undef ALLOCATE_WORK
  return 1;
}

static void capture_rejected_core_input (
  MusehaneModelImplementation *implementation,
  const MusehaneStateSnapshot *before,
  const MusehaneStepWorkspace *workspace,
  double dt)
{
  implementation->rejected_core_input_valid = false;
  if (!before || !workspace || before->cell_count != implementation->ring_count)
    return;
  const size_t cells = implementation->ring_count;
  for (size_t cell = 0; cell < cells; cell++) {
    implementation->rejected_core_active[cell] = workspace->active_np1[cell];
    implementation->rejected_core_estimated[cell] =
      workspace->estimated[cell];
    implementation->rejected_core_old[cell] =
      workspace->active_np1[cell] && before->active[cell] ?
        before->thickness[cell] : workspace->estimated[cell];
    implementation->rejected_core_divergence[cell] =
      workspace->divergence[cell];
    implementation->rejected_core_gradient_cell[cell] =
      workspace->smoothed_gradient[cell];
  }
  memcpy(implementation->rejected_core_gradient_face,
    workspace->smoothed_gradient_face, (cells + 1)*sizeof(double));
  implementation->rejected_core_time = before->time;
  implementation->rejected_core_dt = dt;
  implementation->rejected_core_input_valid = true;
}

static int prepare_observation (
  const MusehaneModelImplementation *implementation, const FilmDNSView *dns,
  MusehaneStepWorkspace *workspace)
{
  const size_t cells = dns->ring_count;
  for (size_t cell = 0; cell < cells; cell++) {
    const FilmRingObservation *ring = &dns->rings[cell];
    const double lower_magnitude = hypot(ring->lower_normal_x,
      ring->lower_normal_r);
    const double upper_magnitude = hypot(ring->upper_normal_x,
      ring->upper_normal_r);
    const double normal_alignment =
      (ring->lower_normal_x*ring->upper_normal_x +
       ring->lower_normal_r*ring->upper_normal_r)/
      (lower_magnitude*upper_magnitude);
    if (!isfinite(normal_alignment) || normal_alignment <= 0.)
      return 0;
    workspace->estimated[cell] = ring->gap;
    workspace->raw_gradient[cell] = ring->pressure_gradient_tangent;
    workspace->capillary_gradient[cell] =
      ring->capillary_force_density_tangent;
    workspace->lower_velocity[cell] = ring->lower_tangential_velocity;
    workspace->upper_velocity[cell] = ring->upper_tangential_velocity;
    workspace->area_weight[cell] = ring->area_weight;
    workspace->surface_velocity[cell] = .5*(ring->lower_tangential_velocity +
      ring->upper_tangential_velocity);
    workspace->local_delta[cell] = ring->local_delta;
    workspace->lower_identity_support[cell] =
      ring->lower_identity_valid ? 1. : 0.;
    workspace->upper_identity_support[cell] =
      ring->upper_identity_valid ? 1. : 0.;
  }

  workspace->surface_velocity_face[0] = 0.;
  for (size_t face = 1; face < cells; face++) {
    const double left_radius = implementation->radius[face - 1];
    const double right_radius = implementation->radius[face];
    const double face_radius = implementation->face_measure[face];
    const double fraction = (face_radius - left_radius)/
      (right_radius - left_radius);
    workspace->surface_velocity_face[face] =
      workspace->surface_velocity[face - 1] + fraction*
      (workspace->surface_velocity[face] -
       workspace->surface_velocity[face - 1]);
  }
  const size_t last = cells - 1;
  if (cells == 1) {
    const double center_radius = implementation->radius[0];
    workspace->surface_velocity_face[1] = center_radius > 0. ?
      workspace->surface_velocity[0]*implementation->face_measure[1]/
        center_radius : 0.;
  }
  else {
    const double distance = implementation->radius[last] -
      implementation->radius[last - 1];
    const double gradient = (workspace->surface_velocity[last] -
      workspace->surface_velocity[last - 1])/distance;
    workspace->surface_velocity_face[cells] =
      workspace->surface_velocity[last] + gradient*
      (implementation->face_measure[cells] - implementation->radius[last]);
  }
  for (size_t cell = 0; cell < cells; cell++) {
    workspace->divergence[cell] =
      (implementation->face_measure[cell + 1]*
         workspace->surface_velocity_face[cell + 1] -
       implementation->face_measure[cell]*
         workspace->surface_velocity_face[cell])/
      implementation->measure[cell];
    if (!isfinite(workspace->divergence[cell]))
      return 0;
  }
  return 1;
}

static FilmModelStatus map_state_status (MusehaneStateStatus status)
{
  if (status == MUSEHANE_STATE_INVALID_CONFIG)
    return FILM_MODEL_INVALID_CONFIG;
  if (status == MUSEHANE_STATE_INVALID_INPUT ||
      status == MUSEHANE_STATE_RESTART_MISMATCH)
    return FILM_MODEL_INVALID_STATE;
  return FILM_MODEL_STEP_REJECTED;
}

FilmModelResult film_model_advance (
  FilmModelState *state, const FilmModelConfig *config,
  const FilmDNSView *dns, double dt)
{
  if (!config_is_valid(config))
    return rejected_result(FILM_MODEL_INVALID_CONFIG,
      "invalid Musehane model config");
  MusehaneModelImplementation *implementation = implementation_of(state);
  if (!implementation || !configs_match(config, &implementation->config))
    return rejected_result(FILM_MODEL_INVALID_STATE,
      "invalid or mismatched Musehane model state");
  implementation->rejected_core_input_valid = false;
  if (!film_dns_view_has_usable_geometry(dns) || dt != dns->dt ||
      !mesh_matches_dns(implementation, dns))
    return rejected_result(FILM_MODEL_INVALID_DNS_VIEW,
      "invalid Musehane DNS observation");

  MusehaneStepWorkspace workspace;
  if (!allocate_workspace(&workspace, dns->ring_count))
    return rejected_result(FILM_MODEL_STEP_REJECTED,
      "Musehane step workspace allocation failed");
  MusehaneStateSnapshot before = {
    .capacity = dns->ring_count,
    .estimated_thickness = workspace.estimated_n,
    .thickness = workspace.thickness_n,
    .active = workspace.active_n
  };
  const MusehaneStateStatus snapshot_status = musehane_state_snapshot(
    implementation->evolution, &before);
  if (snapshot_status != MUSEHANE_STATE_OK ||
      !close_enough_at_scale(before.time, dns->time, dns->dt) ||
      (implementation->completed_steps == 0 &&
       dns->iteration != implementation->last_iteration) ||
      (implementation->completed_steps > 0 &&
       dns->iteration <= implementation->last_iteration)) {
    free_workspace(&workspace);
    return rejected_result(FILM_MODEL_INVALID_DNS_VIEW,
      "Musehane DNS time or iteration is not continuous");
  }
  if (!prepare_observation(implementation, dns, &workspace)) {
    free_workspace(&workspace);
    return rejected_result(FILM_MODEL_INVALID_DNS_VIEW,
      "Musehane surface kinematics reconstruction failed");
  }

  const MusehaneActivationConfig activation_config = {
    .activation_cells = config->activation_cells,
    .minimum_smoothed_fraction = 1.e-3
  };
  const MusehaneActivationInput activation_input = {
    .cell_count = dns->ring_count,
    .estimated_separation = workspace.estimated,
    .local_delta = workspace.local_delta,
    .lower_smoothed_fraction = workspace.lower_identity_support,
    .upper_smoothed_fraction = workspace.upper_identity_support
  };
  MusehaneActivationOutput activation_output = {
    .active = workspace.active_np1
  };
  if (musehane_evaluate_activation(&activation_config, &activation_input,
      &activation_output) != MUSEHANE_ACTIVATION_OK) {
    free_workspace(&workspace);
    return rejected_result(FILM_MODEL_STEP_REJECTED,
      "Musehane activation evaluation failed");
  }
  for (size_t cell = 0; cell < dns->ring_count; cell++) {
    workspace.resolved_separation_threshold[cell] =
      config->activation_cells*workspace.local_delta[cell];
    if (dns->rings[cell].resolved_release_latched)
      workspace.active_np1[cell] = 0;
  }

  const MusehanePressureDriverInput pressure_input = {
    .mesh = &implementation->mesh,
    .raw_pressure_gradient_tangent_cell = workspace.raw_gradient,
    .capillary_force_density_tangent_cell = workspace.capillary_gradient
  };
  MusehanePressureDriverOutput pressure_output = {
    .decontaminated_gradient_cell = workspace.decontaminated_gradient,
    .smoothed_gradient_cell = workspace.smoothed_gradient,
    .smoothed_gradient_face = workspace.smoothed_gradient_face
  };
  const MusehanePressureDriverStatus pressure_status =
    musehane_prepare_pressure_driver(&config->pressure_driver,
      &pressure_input, &pressure_output);
  if (pressure_status != MUSEHANE_PRESSURE_DRIVER_OK) {
    const char *reason = "Musehane pressure driver rejected";
    if (pressure_status == MUSEHANE_PRESSURE_DRIVER_RESIDUAL_MISMATCH)
      reason = "Musehane pressure driver residual mismatch";
    else if (pressure_status == MUSEHANE_PRESSURE_DRIVER_SINGULAR)
      reason = "Musehane pressure driver singular system";
    else if (pressure_status == MUSEHANE_PRESSURE_DRIVER_INVALID_INPUT)
      reason = "Musehane pressure driver invalid input";
    else if (pressure_status == MUSEHANE_PRESSURE_DRIVER_ALLOCATION_FAILED)
      reason = "Musehane pressure driver allocation failed";
    free_workspace(&workspace);
    return rejected_result(FILM_MODEL_STEP_REJECTED, reason);
  }

  MusehaneState *candidate_state = NULL;
  if (musehane_state_restore(&implementation->mesh, &before,
      &candidate_state) != MUSEHANE_STATE_OK) {
    free_workspace(&workspace);
    return rejected_result(FILM_MODEL_STEP_REJECTED,
      "Musehane candidate state restore failed");
  }
  const MusehaneStateAdvanceInput state_input = {
    .mesh = &implementation->mesh,
    .active = workspace.active_np1,
    .estimated_thickness = workspace.estimated,
    .resolved_separation_threshold =
      workspace.resolved_separation_threshold,
    .mean_surface_velocity_divergence = workspace.divergence,
    .decontaminated_pressure_gradient_cell = workspace.smoothed_gradient,
    .decontaminated_pressure_gradient_face =
      workspace.smoothed_gradient_face,
    .dt = dt
  };
  MusehaneStateAdvanceResult state_result = {0};
  const MusehaneStateStatus advance_status = musehane_state_advance(
    &config->equation, candidate_state, &state_input, &state_result);
  if (advance_status != MUSEHANE_STATE_OK) {
    const char *reason = "Musehane thickness update rejected";
    if (advance_status == MUSEHANE_STATE_CORE_REJECTED) {
      switch (state_result.core_status) {
        case MUSEHANE_NONPOSITIVE_THICKNESS:
          reason = "Musehane nonpositive thickness or step denominator";
          break;
        case MUSEHANE_NONLINEAR_NOT_CONVERGED:
          reason = "Musehane thickness nonlinear iteration not converged";
          break;
        case MUSEHANE_INVENTORY_MISMATCH:
          reason = "Musehane thickness inventory mismatch";
          break;
        default:
          break;
      }
    }
    if (advance_status == MUSEHANE_STATE_CORE_REJECTED)
      capture_rejected_core_input(implementation, &before, &workspace, dt);
    musehane_state_destroy(candidate_state);
    free_workspace(&workspace);
    FilmModelResult rejected = rejected_result(
      map_state_status(advance_status), reason);
    rejected.failure_diagnostic_valid = state_result.core_failure.valid;
    rejected.failure_kind = state_result.core_failure.kind;
    rejected.failure_cell = state_result.core_failure.cell_index;
    rejected.failure_iteration = state_result.core_failure.iteration;
    rejected.failure_line_search_backtracks =
      state_result.core_failure.line_search_backtracks;
    rejected.failure_initial_merit =
      state_result.core_failure.initial_merit;
    rejected.failure_final_merit = state_result.core_failure.final_merit;
    rejected.failure_numerator = state_result.core_failure.numerator;
    rejected.failure_denominator = state_result.core_failure.denominator;
    rejected.failure_old_thickness =
      state_result.core_failure.old_thickness;
    rejected.failure_current_thickness =
      state_result.core_failure.current_thickness;
    rejected.failure_flux_in = state_result.core_failure.flux_in;
    rejected.failure_flux_out = state_result.core_failure.flux_out;
    rejected.failure_surface_velocity_divergence =
      state_result.core_failure.surface_velocity_divergence;
    rejected.failure_equation_relative_residual =
      state_result.core_failure.equation_relative_residual;
    rejected.failure_inventory_diagnostic_valid =
      state_result.core_failure.inventory_diagnostic_valid;
    rejected.failure_inventory_before =
      state_result.core_failure.inventory_before;
    rejected.failure_inventory_after =
      state_result.core_failure.inventory_after;
    rejected.failure_expected_inventory_after =
      state_result.core_failure.expected_inventory_after;
    rejected.failure_inventory_closure_error =
      state_result.core_failure.inventory_closure_error;
    rejected.failure_inventory_relative_error =
      state_result.core_failure.inventory_relative_error;
    rejected.failure_inventory_relative_tolerance =
      state_result.core_failure.inventory_relative_tolerance;
    return rejected;
  }

  MusehaneStateSnapshot after = {
    .capacity = dns->ring_count,
    .estimated_thickness = workspace.estimated,
    .thickness = workspace.thickness_np1,
    .active = workspace.active_np1
  };
  if (musehane_state_snapshot(candidate_state, &after) != MUSEHANE_STATE_OK) {
    musehane_state_destroy(candidate_state);
    free_workspace(&workspace);
    return rejected_result(FILM_MODEL_STEP_REJECTED,
      "Musehane candidate snapshot rejected");
  }
  size_t final_active_cells = 0;
  for (size_t cell = 0; cell < dns->ring_count; cell++)
    final_active_cells += workspace.active_np1[cell] != 0;

  for (size_t cell = 0; cell < dns->ring_count; cell++)
    workspace.traction_active[cell] = workspace.active_n[cell] &&
      workspace.active_np1[cell];
  const MusehaneTractionInput traction_input = {
    .cell_count = dns->ring_count,
    .water_viscosity = config->equation.water_viscosity,
    .active = workspace.traction_active,
    .thickness = workspace.thickness_n,
    .pressure_gradient_tangent = workspace.smoothed_gradient,
    .lower_tangential_velocity = workspace.lower_velocity,
    .upper_tangential_velocity = workspace.upper_velocity
  };
  MusehaneTractionOutput traction_output = {
    .lower_tangential_traction = workspace.lower_analytic_traction,
    .upper_tangential_traction = workspace.upper_analytic_traction
  };
  const MusehaneTractionStatus traction_status =
    musehane_evaluate_traction(&traction_input, &traction_output);
  if (traction_status != MUSEHANE_TRACTION_OK) {
    musehane_state_destroy(candidate_state);
    free_workspace(&workspace);
    return traction_rejected_result(
      "Musehane analytical traction rejected", (int)traction_status,
      SIZE_MAX, NULL, NAN, NAN);
  }

  MusehaneMatchedTractionOutput matched_output = {
    .lower_correction_traction = workspace.lower_matched_correction,
    .upper_correction_traction = workspace.upper_matched_correction
  };
  MusehaneEquation27CorrectionOutput equation27_output = {
    .lower_correction_traction = workspace.lower_matched_correction,
    .upper_correction_traction = workspace.upper_matched_correction,
    .resolved_reaction_force = workspace.resolved_reaction_force
  };
  if (config->traction_policy == MUSEHANE_MATCHED_VISCOUS_CORRECTION) {
    const MusehaneMatchedTractionInput matched_input = {
      .cell_count = dns->ring_count,
      .water_viscosity = config->equation.water_viscosity,
      .matching_cells = config->activation_cells,
      .active = workspace.traction_active,
      .thickness = workspace.thickness_n,
      .local_delta = workspace.local_delta,
      .area_weight = workspace.area_weight,
      .pressure_gradient_tangent = workspace.smoothed_gradient,
      .lower_tangential_velocity = workspace.lower_velocity,
      .upper_tangential_velocity = workspace.upper_velocity
    };
    if (musehane_evaluate_matched_traction(&matched_input, &matched_output) !=
        MUSEHANE_MATCHED_TRACTION_OK) {
      musehane_state_destroy(candidate_state);
      free_workspace(&workspace);
      return traction_rejected_result(
        "Musehane matched viscous correction rejected", -1,
        SIZE_MAX, NULL, NAN, NAN);
    }
    if (matched_output.correction_power > 0. ||
        !close_enough_at_scale(matched_output.correction_power,
          -matched_output.additional_relative_shear_dissipation,
          matched_output.additional_relative_shear_dissipation)) {
      musehane_state_destroy(candidate_state);
      free_workspace(&workspace);
      return traction_rejected_result(
        "Musehane matched viscous power audit failed", 0,
        SIZE_MAX, NULL, NAN, NAN);
    }
  }
  else if (config->traction_policy ==
      MUSEHANE_MATCHED_EQUATION27_CORRECTION) {
    const MusehaneEquation27CorrectionInput equation27_input = {
      .cell_count = dns->ring_count,
      .film_viscosity = config->equation.water_viscosity,
      .matching_cells = config->activation_cells,
      .active = workspace.traction_active,
      .thickness = workspace.thickness_n,
      .local_delta = workspace.local_delta,
      .area_weight = workspace.area_weight,
      .pressure_gradient_tangent = workspace.smoothed_gradient,
      .lower_tangential_velocity = workspace.lower_velocity,
      .upper_tangential_velocity = workspace.upper_velocity
    };
    const MusehaneEquation27CorrectionStatus equation27_status =
      musehane_evaluate_equation27_correction(&equation27_input,
        &equation27_output);
    if (equation27_status != MUSEHANE_EQUATION27_CORRECTION_OK) {
      musehane_state_destroy(candidate_state);
      free_workspace(&workspace);
      return traction_rejected_result(
        "Musehane matched equation-27 correction rejected",
        (int)equation27_status, SIZE_MAX, &equation27_output, NAN, NAN);
    }
    if (!musehane_equation27_power_audit_is_valid(&equation27_output)) {
      musehane_state_destroy(candidate_state);
      free_workspace(&workspace);
      return traction_rejected_result(
        "Musehane matched equation-27 power audit failed", 0,
        SIZE_MAX, &equation27_output, NAN, NAN);
    }
  }
  else if (config->traction_policy ==
      MUSEHANE_PUBLISHED_EQUATION27_TRACTION) {
    equation27_output.minimum_reference_thickness = INFINITY;
    for (size_t cell = 0; cell < dns->ring_count; cell++) {
      if (!workspace.traction_active[cell])
        continue;
      const double lower = workspace.lower_analytic_traction[cell];
      const double upper = workspace.upper_analytic_traction[cell];
      const double pressure = .5*(lower + upper);
      const double shear = .5*(upper - lower);
      const double area = workspace.area_weight[cell];
      const double lower_force = area*lower;
      const double upper_force = area*upper;
      workspace.lower_matched_correction[cell] = lower;
      workspace.upper_matched_correction[cell] = upper;
      workspace.resolved_reaction_force[cell] =
        -(lower_force + upper_force);
      equation27_output.corrected_cells++;
      equation27_output.minimum_reference_thickness = fmin(
        equation27_output.minimum_reference_thickness,
        workspace.thickness_n[cell]);
      equation27_output.maximum_pressure_correction_traction = fmax(
        equation27_output.maximum_pressure_correction_traction,
        fabs(pressure));
      equation27_output.maximum_shear_correction_traction = fmax(
        equation27_output.maximum_shear_correction_traction,
        fabs(shear));
      equation27_output.total_resolved_reaction_force +=
        workspace.resolved_reaction_force[cell];
      equation27_output.maximum_pair_closure_residual = fmax(
        equation27_output.maximum_pair_closure_residual,
        fabs(lower_force + upper_force +
          workspace.resolved_reaction_force[cell]));
      equation27_output.pressure_correction_power += area*pressure*
        (workspace.lower_velocity[cell] + workspace.upper_velocity[cell]);
      equation27_output.shear_correction_power -= area*shear*
        (workspace.lower_velocity[cell] - workspace.upper_velocity[cell]);
      equation27_output.direct_correction_power += area*(
        lower*workspace.lower_velocity[cell] +
        upper*workspace.upper_velocity[cell]);
    }
    if (!isfinite(equation27_output.minimum_reference_thickness))
      equation27_output.minimum_reference_thickness = 0.;
    equation27_output.total_correction_power =
      equation27_output.pressure_correction_power +
      equation27_output.shear_correction_power;
    equation27_output.power_identity_residual =
      equation27_output.direct_correction_power -
      equation27_output.total_correction_power;
    if (!musehane_equation27_power_audit_is_valid(&equation27_output)) {
      musehane_state_destroy(candidate_state);
      free_workspace(&workspace);
      return traction_rejected_result(
        "Musehane published equation-27 power audit failed", 0,
        SIZE_MAX, &equation27_output, NAN, NAN);
    }
  }

  double paired_force_residual = 0.;
  for (size_t cell = 0; cell < dns->ring_count; cell++) {
    FilmRingTraction *ring_traction = &workspace.tractions[cell];
    ring_traction->lower_tangential_force =
      workspace.lower_matched_correction[cell]*workspace.area_weight[cell];
    ring_traction->upper_tangential_force =
      workspace.upper_matched_correction[cell]*workspace.area_weight[cell];
    ring_traction->resolved_tangential_reaction_force =
      workspace.resolved_reaction_force[cell];
    if (!isfinite(ring_traction->lower_tangential_force) ||
        !isfinite(ring_traction->upper_tangential_force) ||
        !isfinite(ring_traction->resolved_tangential_reaction_force)) {
      musehane_state_destroy(candidate_state);
      free_workspace(&workspace);
      return traction_rejected_result(
        "Musehane ring-integrated correction is nonfinite", 0,
        cell, &equation27_output, NAN, NAN);
    }
    paired_force_residual = fmax(paired_force_residual, fabs(
      ring_traction->lower_tangential_force +
      ring_traction->upper_tangential_force +
      ring_traction->resolved_tangential_reaction_force));
  }
  const double applied_correction_power = config->traction_policy ==
    MUSEHANE_MATCHED_VISCOUS_CORRECTION ? matched_output.correction_power :
    (config->traction_policy == MUSEHANE_MATCHED_EQUATION27_CORRECTION ||
     config->traction_policy == MUSEHANE_PUBLISHED_EQUATION27_TRACTION) ?
      equation27_output.total_correction_power : 0.;
  const double step_matched_work = applied_correction_power*dt;
  const double cumulative_matched_work =
    implementation->cumulative_matched_correction_work + step_matched_work;
  if (!isfinite(step_matched_work) || !isfinite(cumulative_matched_work)) {
    musehane_state_destroy(candidate_state);
    free_workspace(&workspace);
    return traction_rejected_result(
      "Musehane matched correction work is nonfinite", 0,
      SIZE_MAX, &equation27_output, step_matched_work,
      cumulative_matched_work);
  }

  double minimum = INFINITY;
  double inventory = 0.;
  for (size_t cell = 0; cell < dns->ring_count; cell++) {
    if (workspace.active_np1[cell]) {
      minimum = fmin(minimum, workspace.thickness_np1[cell]);
      inventory += TWO_PI*implementation->measure[cell]*
        workspace.thickness_np1[cell];
    }
  }
  if (final_active_cells == 0) {
    for (size_t cell = 0; cell < dns->ring_count; cell++)
      minimum = fmin(minimum, workspace.estimated[cell]);
  }
  if (!isfinite(minimum) || !isfinite(inventory) || inventory < 0.) {
    musehane_state_destroy(candidate_state);
    free_workspace(&workspace);
    return rejected_result(FILM_MODEL_STEP_REJECTED,
      "Musehane diagnostics are nonfinite");
  }

  musehane_state_destroy(implementation->evolution);
  implementation->evolution = candidate_state;
  implementation->completed_steps++;
  implementation->last_iteration = dns->iteration;
  implementation->substep_count = state_result.substep_count;
  implementation->resolved_release_count =
    state_result.resolved_release_count;
  implementation->diagnostics_valid = true;
  implementation->time_n = before.time;
  implementation->time_np1 = after.time;
  implementation->pressure_residual =
    pressure_output.maximum_linear_relative_residual;
  implementation->traction_residual =
    traction_output.maximum_constitutive_residual;
  implementation->pressure_inventory_change =
    state_result.pressure_inventory_change;
  implementation->kinematic_inventory_change =
    state_result.kinematic_inventory_change;
  implementation->resolved_release_inventory_change =
    state_result.resolved_release_inventory_change;
  implementation->dissipation_rate = state_result.poiseuille_dissipation_rate;
  implementation->cumulative_dissipation =
    state_result.cumulative_poiseuille_dissipation;
  implementation->matched_corrected_ring_count =
    (config->traction_policy == MUSEHANE_MATCHED_EQUATION27_CORRECTION ||
     config->traction_policy == MUSEHANE_PUBLISHED_EQUATION27_TRACTION) ?
      equation27_output.corrected_cells : matched_output.corrected_cells;
  implementation->minimum_matched_reference_thickness =
    (config->traction_policy == MUSEHANE_MATCHED_EQUATION27_CORRECTION ||
     config->traction_policy == MUSEHANE_PUBLISHED_EQUATION27_TRACTION) ?
      equation27_output.minimum_reference_thickness :
      matched_output.minimum_reference_thickness;
  implementation->maximum_resolved_region_correction =
    (config->traction_policy == MUSEHANE_MATCHED_EQUATION27_CORRECTION ||
     config->traction_policy == MUSEHANE_PUBLISHED_EQUATION27_TRACTION) ?
      equation27_output.maximum_resolved_region_correction :
      matched_output.maximum_resolved_region_correction;
  implementation->analytical_film_power =
    matched_output.analytical_film_power;
  implementation->resolved_reference_power =
    matched_output.resolved_reference_power;
  implementation->matched_correction_power = matched_output.correction_power;
  implementation->additional_relative_shear_dissipation_rate =
    matched_output.additional_relative_shear_dissipation;
  implementation->step_matched_correction_work = step_matched_work;
  implementation->cumulative_matched_correction_work =
    cumulative_matched_work;
  implementation->equation27_pressure_correction_power =
    equation27_output.pressure_correction_power;
  implementation->equation27_shear_correction_power =
    equation27_output.shear_correction_power;
  implementation->equation27_correction_power =
    equation27_output.total_correction_power;
  implementation->equation27_total_resolved_reaction_force =
    equation27_output.total_resolved_reaction_force;
  implementation->maximum_equation27_pair_closure_residual =
    equation27_output.maximum_pair_closure_residual;
  memcpy(implementation->active_n, workspace.active_n,
    dns->ring_count*sizeof *implementation->active_n);
  memcpy(implementation->active_np1, workspace.active_np1,
    dns->ring_count*sizeof *implementation->active_np1);
  memcpy(implementation->traction_active, workspace.traction_active,
    dns->ring_count*sizeof *implementation->traction_active);
  memcpy(implementation->estimated, workspace.estimated,
    dns->ring_count*sizeof *implementation->estimated);
  memcpy(implementation->thickness_n, workspace.thickness_n,
    dns->ring_count*sizeof *implementation->thickness_n);
  memcpy(implementation->thickness_np1, workspace.thickness_np1,
    dns->ring_count*sizeof *implementation->thickness_np1);
  memcpy(implementation->divergence, workspace.divergence,
    dns->ring_count*sizeof *implementation->divergence);
  memcpy(implementation->decontaminated_gradient,
    workspace.decontaminated_gradient,
    dns->ring_count*sizeof *implementation->decontaminated_gradient);
  memcpy(implementation->smoothed_gradient, workspace.smoothed_gradient,
    dns->ring_count*sizeof *implementation->smoothed_gradient);
  memcpy(implementation->lower_analytic_traction,
    workspace.lower_analytic_traction,
    dns->ring_count*sizeof *implementation->lower_analytic_traction);
  memcpy(implementation->upper_analytic_traction,
    workspace.upper_analytic_traction,
    dns->ring_count*sizeof *implementation->upper_analytic_traction);
  memcpy(implementation->lower_matched_correction,
    workspace.lower_matched_correction,
    dns->ring_count*sizeof *implementation->lower_matched_correction);
  memcpy(implementation->upper_matched_correction,
    workspace.upper_matched_correction,
    dns->ring_count*sizeof *implementation->upper_matched_correction);
  memcpy(implementation->resolved_reaction_force,
    workspace.resolved_reaction_force,
    dns->ring_count*sizeof *implementation->resolved_reaction_force);
  memcpy(implementation->tractions, workspace.tractions,
    dns->ring_count*sizeof *implementation->tractions);
  free_workspace(&workspace);

  return (FilmModelResult) {
    .status = FILM_MODEL_OK,
    .model_kind = FILM_MODEL_MUSEHANE_INVERSE,
    .active = final_active_cells > 0,
    .commit_allowed = true,
    .ring_count = dns->ring_count,
    .tractions = implementation->tractions,
    .minimum_thickness = minimum,
    .inventory = inventory,
    .paired_force_residual = paired_force_residual,
    .paired_moment_residual = 0.,
    .step_work = step_matched_work,
    .cumulative_work = cumulative_matched_work,
    .completed_stage = FILM_STAGE_EVALUATE_TRACTION_N,
    .reason = config->traction_policy == MUSEHANE_MATCHED_VISCOUS_CORRECTION ?
      "Musehane matched viscous correction evaluated" :
      config->traction_policy == MUSEHANE_MATCHED_EQUATION27_CORRECTION ?
        "Musehane matched equation-27 correction evaluated" :
      config->traction_policy == MUSEHANE_PUBLISHED_EQUATION27_TRACTION ?
        "Musehane published equation-27 traction evaluated" :
        "Musehane diagnostic-only step"
  };
}

FilmModelStatus film_model_apply_traction (
  const FilmModelResult *result, FilmDNSTractionTarget *target)
{
  if (!result || !target || result->status != FILM_MODEL_OK ||
      result->model_kind != FILM_MODEL_MUSEHANE_INVERSE ||
      !result->commit_allowed || result->ring_count == 0 ||
      !result->tractions || target->ring_count != result->ring_count ||
      !target->lower_normal_force || !target->upper_normal_force ||
      !target->lower_tangential_force || !target->upper_tangential_force)
    return FILM_MODEL_TRACTION_REJECTED;
  for (size_t cell = 0; cell < result->ring_count; cell++) {
    const FilmRingTraction *traction = &result->tractions[cell];
    if (traction->lower_normal_force != 0. ||
        traction->upper_normal_force != 0. ||
        !isfinite(traction->lower_tangential_force) ||
        !isfinite(traction->upper_tangential_force) ||
        !isfinite(traction->resolved_tangential_reaction_force) ||
        !isfinite(target->lower_normal_force[cell]) ||
        !isfinite(target->upper_normal_force[cell]) ||
        !isfinite(target->lower_tangential_force[cell]) ||
        !isfinite(target->upper_tangential_force[cell]) ||
        !isfinite(target->lower_tangential_force[cell] +
          traction->lower_tangential_force) ||
        !isfinite(target->upper_tangential_force[cell] +
          traction->upper_tangential_force))
      return FILM_MODEL_TRACTION_REJECTED;
  }
  for (size_t cell = 0; cell < result->ring_count; cell++) {
    target->lower_tangential_force[cell] +=
      result->tractions[cell].lower_tangential_force;
    target->upper_tangential_force[cell] +=
      result->tractions[cell].upper_tangential_force;
  }
  return FILM_MODEL_OK;
}

static int checkpoint_size_for_cells (size_t cells, size_t *size)
{
  if (!size || cells == 0 || cells > (SIZE_MAX - 1)/6)
    return 0;
  const size_t doubles = 6*cells + 1;
  if (doubles > (SIZE_MAX - sizeof(MusehaneRestartHeader) - cells)/
      sizeof(double))
    return 0;
  *size = sizeof(MusehaneRestartHeader) + doubles*sizeof(double) + cells;
  return 1;
}

size_t film_model_checkpoint_size (const FilmModelState *state)
{
  const MusehaneModelImplementation *implementation =
    constant_implementation_of(state);
  size_t size = 0;
  return implementation && checkpoint_size_for_cells(
    implementation->ring_count, &size) ? size : 0;
}

static unsigned char *write_bytes (
  unsigned char *cursor, const void *source, size_t size)
{
  memcpy(cursor, source, size);
  return cursor + size;
}

static const unsigned char *read_bytes (
  const unsigned char *cursor, void *target, size_t size)
{
  memcpy(target, cursor, size);
  return cursor + size;
}

FilmModelStatus film_model_checkpoint (
  const FilmModelState *state, FilmRestartRecord *record)
{
  const MusehaneModelImplementation *implementation =
    constant_implementation_of(state);
  const size_t required = film_model_checkpoint_size(state);
  if (!implementation || !record || !record->data || required == 0 ||
      record->capacity < required)
    return FILM_MODEL_RESTART_MISMATCH;
  unsigned char *candidate = malloc(required);
  double *estimated = malloc(implementation->ring_count*sizeof(double));
  double *thickness = malloc(implementation->ring_count*sizeof(double));
  unsigned char *active = malloc(implementation->ring_count);
  if (!candidate || !estimated || !thickness || !active) {
    free(candidate);
    free(estimated);
    free(thickness);
    free(active);
    return FILM_MODEL_RESTART_MISMATCH;
  }
  MusehaneStateSnapshot snapshot = {
    .capacity = implementation->ring_count,
    .estimated_thickness = estimated,
    .thickness = thickness,
    .active = active
  };
  if (musehane_state_snapshot(implementation->evolution, &snapshot) !=
      MUSEHANE_STATE_OK ||
      snapshot.update_count != implementation->completed_steps) {
    free(candidate);
    free(estimated);
    free(thickness);
    free(active);
    return FILM_MODEL_RESTART_MISMATCH;
  }
  MusehaneRestartHeader header;
  memset(&header, 0, sizeof header);
  header.magic = MUSEHANE_RESTART_MAGIC;
  header.schema_version = MUSEHANE_FILM_MODEL_RESTART_SCHEMA;
  header.contract_version = FILM_MODEL_CONTRACT_VERSION;
  header.model_kind = FILM_MODEL_MUSEHANE_INVERSE;
  header.traction_policy = implementation->config.traction_policy;
  header.ring_count = implementation->ring_count;
  header.completed_steps = implementation->completed_steps;
  header.last_iteration = implementation->last_iteration;
  header.mesh_signature = snapshot.mesh_signature;
  header.time = snapshot.time;
  header.cumulative_dissipation =
    snapshot.cumulative_poiseuille_dissipation;
  header.cumulative_matched_correction_work =
    implementation->cumulative_matched_correction_work;
  header.water_density = implementation->config.equation.water_density;
  header.water_viscosity = implementation->config.equation.water_viscosity;
  header.nonlinear_tolerance =
    implementation->config.equation.nonlinear_relative_tolerance;
  header.inventory_tolerance =
    implementation->config.equation.inventory_relative_tolerance;
  header.maximum_iterations =
    implementation->config.equation.maximum_iterations;
  header.discretization = (uint32_t)
    implementation->config.equation.discretization;
  header.smoothing_length_squared =
    implementation->config.pressure_driver.smoothing_length_squared;
  header.pressure_residual_tolerance = implementation->config.pressure_driver.
    linear_residual_relative_tolerance;
  header.activation_cells = implementation->config.activation_cells;
  unsigned char *cursor = candidate;
  cursor = write_bytes(cursor, &header, sizeof header);
  cursor = write_bytes(cursor, implementation->radius,
    implementation->ring_count*sizeof(double));
  cursor = write_bytes(cursor, implementation->width,
    implementation->ring_count*sizeof(double));
  cursor = write_bytes(cursor, implementation->measure,
    implementation->ring_count*sizeof(double));
  cursor = write_bytes(cursor, implementation->face_measure,
    (implementation->ring_count + 1)*sizeof(double));
  cursor = write_bytes(cursor, estimated,
    implementation->ring_count*sizeof(double));
  cursor = write_bytes(cursor, thickness,
    implementation->ring_count*sizeof(double));
  cursor = write_bytes(cursor, active, implementation->ring_count);
  if ((size_t)(cursor - candidate) != required) {
    free(candidate);
    free(estimated);
    free(thickness);
    free(active);
    return FILM_MODEL_RESTART_MISMATCH;
  }
  memcpy(record->data, candidate, required);
  record->schema_version = FILM_MODEL_RESTART_SCHEMA_VERSION;
  record->model_kind = FILM_MODEL_MUSEHANE_INVERSE;
  record->size = required;
  free(candidate);
  free(estimated);
  free(thickness);
  free(active);
  return FILM_MODEL_OK;
}

static int restart_header_matches_config (
  const MusehaneRestartHeader *header, const FilmModelConfig *config)
{
  return header && config && header->magic == MUSEHANE_RESTART_MAGIC &&
    header->schema_version == MUSEHANE_FILM_MODEL_RESTART_SCHEMA &&
    header->contract_version == config->contract_version &&
    header->model_kind == (uint32_t) config->model_kind &&
    header->traction_policy == (uint32_t) config->traction_policy &&
    header->water_density == config->equation.water_density &&
    header->water_viscosity == config->equation.water_viscosity &&
    header->nonlinear_tolerance ==
      config->equation.nonlinear_relative_tolerance &&
    header->inventory_tolerance ==
      config->equation.inventory_relative_tolerance &&
    header->maximum_iterations == config->equation.maximum_iterations &&
    header->discretization == (uint32_t)config->equation.discretization &&
    header->smoothing_length_squared ==
      config->pressure_driver.smoothing_length_squared &&
    header->pressure_residual_tolerance ==
      config->pressure_driver.linear_residual_relative_tolerance &&
    header->activation_cells == config->activation_cells &&
    isfinite(header->cumulative_matched_correction_work) &&
    header->ring_count > 0 && header->ring_count <= SIZE_MAX;
}

FilmModelStatus film_model_restore (
  FilmModelState *state, const FilmModelConfig *config,
  const FilmRestartRecord *record)
{
  if (!config_is_valid(config))
    return FILM_MODEL_INVALID_CONFIG;
  if (!state || state->guard != 0U || state->implementation != NULL)
    return FILM_MODEL_INVALID_STATE;
  if (!record || !record->data ||
      record->schema_version != FILM_MODEL_RESTART_SCHEMA_VERSION ||
      record->model_kind != FILM_MODEL_MUSEHANE_INVERSE ||
      record->size < sizeof(MusehaneRestartHeader) ||
      record->capacity < record->size)
    return FILM_MODEL_RESTART_MISMATCH;
  MusehaneRestartHeader header;
  const unsigned char *cursor = read_bytes(record->data, &header,
    sizeof header);
  if (!restart_header_matches_config(&header, config))
    return FILM_MODEL_RESTART_MISMATCH;
  if (header.completed_steps > ULONG_MAX || header.last_iteration > ULONG_MAX)
    return FILM_MODEL_RESTART_MISMATCH;
  size_t expected = 0;
  if (!checkpoint_size_for_cells((size_t) header.ring_count, &expected) ||
      expected != record->size)
    return FILM_MODEL_RESTART_MISMATCH;
  MusehaneModelImplementation *candidate =
    allocate_implementation((size_t) header.ring_count);
  if (!candidate)
    return FILM_MODEL_RESTART_MISMATCH;
  candidate->config = *config;
  candidate->completed_steps = (unsigned long) header.completed_steps;
  candidate->last_iteration = (unsigned long) header.last_iteration;
  candidate->cumulative_matched_correction_work =
    header.cumulative_matched_correction_work;
  cursor = read_bytes(cursor, candidate->radius,
    candidate->ring_count*sizeof(double));
  cursor = read_bytes(cursor, candidate->width,
    candidate->ring_count*sizeof(double));
  cursor = read_bytes(cursor, candidate->measure,
    candidate->ring_count*sizeof(double));
  cursor = read_bytes(cursor, candidate->face_measure,
    (candidate->ring_count + 1)*sizeof(double));
  double *estimated = malloc(candidate->ring_count*sizeof(double));
  double *thickness = malloc(candidate->ring_count*sizeof(double));
  unsigned char *active = malloc(candidate->ring_count);
  if (!estimated || !thickness || !active) {
    free(estimated);
    free(thickness);
    free(active);
    free_implementation(candidate);
    return FILM_MODEL_RESTART_MISMATCH;
  }
  cursor = read_bytes(cursor, estimated, candidate->ring_count*sizeof(double));
  cursor = read_bytes(cursor, thickness, candidate->ring_count*sizeof(double));
  cursor = read_bytes(cursor, active, candidate->ring_count);
  if ((size_t)(cursor - record->data) != record->size) {
    free(estimated);
    free(thickness);
    free(active);
    free_implementation(candidate);
    return FILM_MODEL_RESTART_MISMATCH;
  }
  const MusehaneStateSnapshot snapshot = {
    .schema_version = MUSEHANE_STATE_SNAPSHOT_SCHEMA,
    .cell_count = candidate->ring_count,
    .capacity = candidate->ring_count,
    .mesh_signature = header.mesh_signature,
    .time = header.time,
    .update_count = candidate->completed_steps,
    .cumulative_poiseuille_dissipation = header.cumulative_dissipation,
    .estimated_thickness = estimated,
    .thickness = thickness,
    .active = active
  };
  const MusehaneStateStatus restored = musehane_state_restore(
    &candidate->mesh, &snapshot, &candidate->evolution);
  if (restored == MUSEHANE_STATE_OK) {
    memcpy(candidate->estimated, estimated,
      candidate->ring_count*sizeof *candidate->estimated);
    memcpy(candidate->thickness_np1, thickness,
      candidate->ring_count*sizeof *candidate->thickness_np1);
    memcpy(candidate->active_np1, active,
      candidate->ring_count*sizeof *candidate->active_np1);
    candidate->time_np1 = header.time;
  }
  free(estimated);
  free(thickness);
  free(active);
  if (restored != MUSEHANE_STATE_OK) {
    free_implementation(candidate);
    return FILM_MODEL_RESTART_MISMATCH;
  }
  state->guard = MUSEHANE_MODEL_GUARD;
  state->implementation = candidate;
  return FILM_MODEL_OK;
}

FilmModelStatus musehane_film_model_restored_state_view (
  const FilmModelState *state, MusehaneFilmRestoredStateView *output)
{
  const MusehaneModelImplementation *implementation =
    constant_implementation_of(state);
  if (!implementation || !output || implementation->completed_steps == 0 ||
      implementation->diagnostics_valid)
    return FILM_MODEL_INVALID_STATE;
  *output = (MusehaneFilmRestoredStateView) {
    .valid = true,
    .ring_count = implementation->ring_count,
    .completed_steps = implementation->completed_steps,
    .time = implementation->time_np1,
    .active = implementation->active_np1,
    .estimated_thickness = implementation->estimated,
    .thickness = implementation->thickness_np1
  };
  return FILM_MODEL_OK;
}

FilmModelStatus film_model_rebase_restart_clock (
  FilmModelState *state, double new_time, double maximum_gap)
{
  MusehaneModelImplementation *implementation = implementation_of(state);
  if (!implementation || !isfinite(new_time) || new_time < 0. ||
      !isfinite(maximum_gap) || maximum_gap <= 0. ||
      implementation->completed_steps == 0 ||
      implementation->diagnostics_valid ||
      implementation->restart_clock_rebased)
    return FILM_MODEL_INVALID_STATE;

  const size_t cells = implementation->ring_count;
  double *estimated = malloc(cells*sizeof(double));
  double *thickness = malloc(cells*sizeof(double));
  unsigned char *active = malloc(cells);
  if (!estimated || !thickness || !active) {
    free(estimated);
    free(thickness);
    free(active);
    return FILM_MODEL_INVALID_STATE;
  }

  MusehaneStateSnapshot snapshot = {
    .capacity = cells,
    .estimated_thickness = estimated,
    .thickness = thickness,
    .active = active
  };
  const MusehaneStateStatus snapshot_status = musehane_state_snapshot(
    implementation->evolution, &snapshot);
  const double clock_gap = new_time - snapshot.time;
  if (snapshot_status != MUSEHANE_STATE_OK || !isfinite(clock_gap) ||
      clock_gap < 0. || clock_gap > maximum_gap ||
      snapshot.update_count != implementation->completed_steps) {
    free(estimated);
    free(thickness);
    free(active);
    return FILM_MODEL_INVALID_STATE;
  }

  const unsigned long update_count = snapshot.update_count;
  const double cumulative_dissipation =
    snapshot.cumulative_poiseuille_dissipation;
  snapshot.time = new_time;
  MusehaneState *candidate = NULL;
  const MusehaneStateStatus restore_status = musehane_state_restore(
    &implementation->mesh, &snapshot, &candidate);
  free(estimated);
  free(thickness);
  free(active);
  if (restore_status != MUSEHANE_STATE_OK || !candidate ||
      snapshot.update_count != update_count ||
      snapshot.cumulative_poiseuille_dissipation !=
        cumulative_dissipation) {
    musehane_state_destroy(candidate);
    return FILM_MODEL_INVALID_STATE;
  }

  MusehaneState *previous = implementation->evolution;
  implementation->evolution = candidate;
  implementation->restart_clock_rebased = true;
  musehane_state_destroy(previous);
  return FILM_MODEL_OK;
}

FilmModelStatus musehane_film_model_diagnostics (
  const FilmModelState *state, MusehaneFilmDiagnostics *output)
{
  const MusehaneModelImplementation *implementation =
    constant_implementation_of(state);
  if (!implementation || !output || !implementation->diagnostics_valid)
    return FILM_MODEL_INVALID_STATE;
  *output = (MusehaneFilmDiagnostics) {
    .valid = true,
    .one_step_lag = true,
    .viscous_correction_applied = implementation->config.traction_policy ==
      MUSEHANE_MATCHED_VISCOUS_CORRECTION,
    .equation27_correction_applied =
      implementation->config.traction_policy ==
        MUSEHANE_MATCHED_EQUATION27_CORRECTION,
    .published_equation27_applied =
      implementation->config.traction_policy ==
        MUSEHANE_PUBLISHED_EQUATION27_TRACTION,
    .ring_count = implementation->ring_count,
    .completed_steps = implementation->completed_steps,
    .substep_count = implementation->substep_count,
    .resolved_release_count = implementation->resolved_release_count,
    .time_n = implementation->time_n,
    .time_np1 = implementation->time_np1,
    .active_n = implementation->active_n,
    .active_np1 = implementation->active_np1,
    .traction_active = implementation->traction_active,
    .estimated_thickness = implementation->estimated,
    .thickness_n = implementation->thickness_n,
    .thickness_np1 = implementation->thickness_np1,
    .mean_surface_velocity_divergence = implementation->divergence,
    .decontaminated_pressure_gradient =
      implementation->decontaminated_gradient,
    .smoothed_pressure_gradient = implementation->smoothed_gradient,
    .lower_analytic_tangential_traction =
      implementation->lower_analytic_traction,
    .upper_analytic_tangential_traction =
      implementation->upper_analytic_traction,
    .lower_matched_correction_traction =
      implementation->lower_matched_correction,
    .upper_matched_correction_traction =
      implementation->upper_matched_correction,
    .matched_corrected_ring_count =
      implementation->matched_corrected_ring_count,
    .minimum_matched_reference_thickness =
      implementation->minimum_matched_reference_thickness,
    .maximum_resolved_region_correction =
      implementation->maximum_resolved_region_correction,
    .analytical_film_power = implementation->analytical_film_power,
    .resolved_reference_power = implementation->resolved_reference_power,
    .matched_correction_power = implementation->matched_correction_power,
    .additional_relative_shear_dissipation_rate =
      implementation->additional_relative_shear_dissipation_rate,
    .step_matched_correction_work =
      implementation->step_matched_correction_work,
    .cumulative_matched_correction_work =
      implementation->cumulative_matched_correction_work,
    .equation27_pressure_correction_power =
      implementation->equation27_pressure_correction_power,
    .equation27_shear_correction_power =
      implementation->equation27_shear_correction_power,
    .equation27_correction_power =
      implementation->equation27_correction_power,
    .equation27_total_resolved_reaction_force =
      implementation->equation27_total_resolved_reaction_force,
    .maximum_equation27_pair_closure_residual =
      implementation->maximum_equation27_pair_closure_residual,
    .maximum_pressure_driver_residual = implementation->pressure_residual,
    .maximum_traction_constitutive_residual =
      implementation->traction_residual,
    .pressure_inventory_change =
      implementation->pressure_inventory_change,
    .kinematic_inventory_change =
      implementation->kinematic_inventory_change,
    .resolved_release_inventory_change =
      implementation->resolved_release_inventory_change,
    .poiseuille_dissipation_rate = implementation->dissipation_rate,
    .cumulative_poiseuille_dissipation =
      implementation->cumulative_dissipation
  };
  return FILM_MODEL_OK;
}

FilmModelStatus musehane_film_model_rejected_core_input_view (
  const FilmModelState *state, MusehaneRejectedCoreInputView *output)
{
  const MusehaneModelImplementation *implementation =
    constant_implementation_of(state);
  if (!implementation || !output ||
      !implementation->rejected_core_input_valid)
    return FILM_MODEL_INVALID_STATE;
  *output = (MusehaneRejectedCoreInputView) {
    .valid = true,
    .config = implementation->config.equation,
    .mesh = implementation->mesh,
    .time = implementation->rejected_core_time,
    .dt = implementation->rejected_core_dt,
    .active = implementation->rejected_core_active,
    .estimated_thickness = implementation->rejected_core_estimated,
    .old_thickness = implementation->rejected_core_old,
    .mean_surface_velocity_divergence =
      implementation->rejected_core_divergence,
    .decontaminated_pressure_gradient_cell =
      implementation->rejected_core_gradient_cell,
    .decontaminated_pressure_gradient_face =
      implementation->rejected_core_gradient_face
  };
  return FILM_MODEL_OK;
}
