#ifndef TRACK_A_MUSEHANE_COUPLED_AXI_H
#define TRACK_A_MUSEHANE_COUPLED_AXI_H

#include <stdint.h>
#include <string.h>

#include "axisymmetric_musehane_native_fields_axi.h"
#include "axisymmetric_musehane_observation_axi.h"
#include "axisymmetric_musehane_traction_axi.h"
#include "film_model/musehane.h"
#include "film_model/musehane_checkpoint_file.h"
#include "film_model/musehane_transaction.h"
#include "film_musehane/beta_viscosity.h"

#ifndef TRACK_A_MUSEHANE_ACTIVATION_CELLS
# define TRACK_A_MUSEHANE_ACTIVATION_CELLS 7.
#endif
#ifndef TRACK_A_MUSEHANE_SMOOTHING_LENGTH_SQUARED
# define TRACK_A_MUSEHANE_SMOOTHING_LENGTH_SQUARED 0.
#endif
#ifndef TRACK_A_MUSEHANE_RESTART_CLOCK_GAP_FACTOR
# define TRACK_A_MUSEHANE_RESTART_CLOCK_GAP_FACTOR 0.25
#endif
#ifndef TRACK_A_MUSEHANE_TRACTION_POLICY
# define TRACK_A_MUSEHANE_TRACTION_POLICY \
  MUSEHANE_MATCHED_VISCOUS_CORRECTION
#endif
#ifndef TRACK_A_MUSEHANE_REPLACE_NATIVE_VISCOSITY
# define TRACK_A_MUSEHANE_REPLACE_NATIVE_VISCOSITY 0
#endif
#ifndef TRACK_A_MUSEHANE_EQ28_AXIAL_PADDING_CELLS
# define TRACK_A_MUSEHANE_EQ28_AXIAL_PADDING_CELLS 2.
#endif
#ifndef TRACK_A_PHASE2C_AUDIT
# define TRACK_A_PHASE2C_AUDIT 0
#endif

typedef enum {
  TRACK_A_COUPLED_UNINITIALIZED = 0,
  TRACK_A_COUPLED_READY = 1,
  TRACK_A_COUPLED_OK = 2,
  TRACK_A_COUPLED_NATIVE_FIELD_FAILED = 3,
  TRACK_A_COUPLED_OBSERVATION_FAILED = 4,
  TRACK_A_COUPLED_MODEL_FAILED = 5,
  TRACK_A_COUPLED_DEPOSITION_FAILED = 6,
  TRACK_A_COUPLED_VELOCITY_FAILED = 7,
  TRACK_A_COUPLED_CHECKPOINT_FAILED = 8,
  TRACK_A_COUPLED_ROLLBACK_FAILED = 9,
  TRACK_A_COUPLED_COMMIT_FAILED = 10,
  TRACK_A_COUPLED_CHECKPOINT_FILE_FAILED = 11,
  TRACK_A_COUPLED_RESTORE_FILE_FAILED = 12,
  TRACK_A_COUPLED_RESTART_CLOCK_FAILED = 13
} TrackAMusehaneCoupledStatus;

typedef struct {
  TrackAMusehaneCoupledStatus status;
  bool initialized;
  unsigned long accepted_steps;
  FilmModelConfig config;
  FilmModelState model_state;
  MusehaneModelTransaction transaction;
  AxisymmetricMusehaneNativeFieldAudit native_audit;
  AxisymmetricMusehaneObservationStatus observation_status;
  AxisymmetricMusehaneObservation observation;
  AxisymmetricMusehaneVelocityLedger velocity_observation;
  FilmModelResult result;
  MusehaneFilmDiagnostics diagnostics;
  AxisymmetricMusehaneTractionLedger deposition;
  MusehaneCheckpointFileAudit checkpoint_audit;
  bool restart_clock_alignment_pending;
  bool restart_clock_rebased;
  double restart_clock_source_time;
  double restart_clock_gap;
  double restart_clock_maximum_gap;
  uint64_t observation_release_latch_mask;
} TrackAMusehaneCoupledState;

static vector track_a_coupled_pressure_gradient[];
static vector track_a_coupled_capillary_force_density[];
static TrackAMusehaneCoupledState track_a_coupled = {
  .status = TRACK_A_COUPLED_UNINITIALIZED,
  .model_state = FILM_MODEL_STATE_INITIALIZER,
  .transaction = MUSEHANE_MODEL_TRANSACTION_INITIALIZER
};
static bool track_a_eq28_reported = false;
static long track_a_phase2c_eq28_masked_faces = 0;
static double track_a_phase2c_eq28_mask_time = NAN;

#if TRACK_A_PHASE2C_AUDIT
typedef struct {
  bool valid;
  size_t ring_count;
  double time_np1;
  unsigned char active_n[AXISYMMETRIC_MATCHED_FILM_BINS];
  unsigned char active_np1[AXISYMMETRIC_MATCHED_FILM_BINS];
  unsigned char traction_active[AXISYMMETRIC_MATCHED_FILM_BINS];
  double estimated_thickness[AXISYMMETRIC_MATCHED_FILM_BINS];
  double thickness_n[AXISYMMETRIC_MATCHED_FILM_BINS];
  double thickness_np1[AXISYMMETRIC_MATCHED_FILM_BINS];
} TrackAPhase2CFailureSnapshot;

static TrackAPhase2CFailureSnapshot track_a_phase2c_failure_snapshot = {0};
#endif

static void track_a_coupled_replace_native_viscosity (void)
{
  track_a_phase2c_eq28_masked_faces = 0;
  track_a_phase2c_eq28_mask_time = t;
  if (!(int)TRACK_A_MUSEHANE_REPLACE_NATIVE_VISCOSITY ||
      !track_a_coupled.initialized || !track_a_coupled.diagnostics.valid)
    return;
  const MusehaneBetaViscosityInput beta_input = {
    .ring_count = track_a_coupled.observation.view.ring_count,
    .rings = track_a_coupled.observation.rings,
    .active = track_a_coupled.diagnostics.traction_active,
    .axial_padding_cells =
      (double)TRACK_A_MUSEHANE_EQ28_AXIAL_PADDING_CELLS
  };
  face vector viscosity = mu;
  long masked_faces = 0;
  foreach_face(reduction(+:masked_faces)) {
    const double beta = musehane_beta_viscosity_at_point(
      &beta_input, x, y, Delta);
    viscosity.x[] *= 1. - beta;
    if (beta > 0.)
      masked_faces++;
  }
  track_a_phase2c_eq28_masked_faces = masked_faces;
  if (masked_faces > 0 && !track_a_eq28_reported && pid() == 0) {
    fprintf(stderr,
      "# TRACK_A_EQ28_NATIVE_VISCOSITY_REPLACED t=%.17g faces=%ld\n",
      t, masked_faces);
    fflush(stderr);
    track_a_eq28_reported = true;
  }
}

static bool track_a_coupled_start (void);

#define COLLISION_TRACTION_EXTENSION_HAS_COMMON_STATE_PREPARE 1
static void collision_traction_extension_prepare_common_state_checkpoint (void)
{
  track_a_coupled_pressure_gradient.x.nodump = true;
  track_a_coupled_pressure_gradient.y.nodump = true;
  track_a_coupled_capillary_force_density.x.nodump = true;
  track_a_coupled_capillary_force_density.y.nodump = true;
}

static bool track_a_coupled_rollback (
  TrackAMusehaneCoupledStatus failure_status)
{
  if (musehane_model_transaction_rollback(&track_a_coupled.transaction,
      &track_a_coupled.model_state, &track_a_coupled.config) !=
      MUSEHANE_MODEL_TRANSACTION_OK) {
    track_a_coupled.status = TRACK_A_COUPLED_ROLLBACK_FAILED;
    return false;
  }
  track_a_coupled.status = failure_status;
  return false;
}

static void track_a_coupled_write_header (void)
{
  if (pid() != 0)
    return;
  FILE *output = fopen("data/track_a_musehane_two_way.dat", "w");
  if (!output) {
    perror("data/track_a_musehane_two_way.dat");
    exit(6);
  }
  fprintf(output,
    "# t i dt status accepted active corrected_rings min_gap min_h "
    "step_work cumulative_work requested_power delivered_power "
    "equation27_pressure_power equation27_shear_power "
    "poiseuille_dissipation_rate cumulative_poiseuille_dissipation "
    "resolved_reaction_force "
    "power_residual requested_lower_power requested_upper_power "
    "delivered_lower_power delivered_upper_power max_ring_power_residual "
    "max_interface_velocity_mismatch max_component_force_residual "
    "max_bin_force_residual pressure_residual traction_residual "
    "velocity_status plic_velocity_adjustment "
    "velocity_geometric_fallback_components "
    "velocity_first_geometric_fallback_bin "
    "velocity_first_geometric_fallback_mask "
    "deposition_geometric_fallback_components restart_rebased "
    "restart_clock_gap restart_clock_maximum_gap "
    "resolved_release_count resolved_release_inventory_change "
    "observation_release_latch_mask observation_release_latch_count\n");
  fclose(output);
  output = fopen("data/track_a_musehane_checkpoints.dat", "w");
  if (!output) {
    perror("data/track_a_musehane_checkpoints.dat");
    exit(6);
  }
  fprintf(output,
    "# time iteration accepted_steps payload_size payload_checksum path\n");
  fclose(output);
  output = fopen("data/track_a_musehane_substeps.dat", "w");
  if (!output) {
    perror("data/track_a_musehane_substeps.dat");
    exit(6);
  }
  fprintf(output, "# t i dt accepted_steps substep_count\n");
  fclose(output);
}

#if TRACK_A_PHASE2C_AUDIT
static void track_a_phase2c_write_headers (void)
{
  if (pid() != 0)
    return;
  FILE *output = fopen("data/track_a_phase2c_steps.dat", "w");
  if (!output) {
    perror("data/track_a_phase2c_steps.dat");
    exit(6);
  }
  fprintf(output,
    "# t i dt model_time_n model_time_np1 accepted_steps model_active "
    "active_rings traction_active_rings lower_center upper_center "
    "center_separation lower_velocity upper_velocity relative_velocity "
    "eq28_mask_time eq28_masked_faces invalid_identity_faces "
    "requested_lower_x_force requested_upper_x_force "
    "requested_lower_r_force requested_upper_r_force "
    "delivered_lower_x_force delivered_upper_x_force "
    "delivered_lower_r_force delivered_upper_r_force "
    "requested_power delivered_power power_residual\n");
  fclose(output);

  output = fopen("data/track_a_phase2c_rings.dat", "w");
  if (!output) {
    perror("data/track_a_phase2c_rings.dat");
    exit(6);
  }
  fprintf(output,
    "# t i dt ring radius width h_native_plic h0 h_track_a_n "
    "h_track_a_np1 active_n active_np1 traction_active "
    "lower_identity_valid upper_identity_valid pressure_gradient_used "
    "eq27_pressure_traction eq27_shear_traction "
    "lower_total_traction upper_total_traction "
    "requested_lower_x_force requested_upper_x_force "
    "requested_lower_r_force requested_upper_r_force "
    "delivered_lower_x_force delivered_upper_x_force "
    "delivered_lower_r_force delivered_upper_r_force\n");
  fclose(output);

  output = fopen("data/track_a_phase2c_failure_rings.dat", "w");
  if (!output) {
    perror("data/track_a_phase2c_failure_rings.dat");
    exit(6);
  }
  fprintf(output,
    "# failure_time i dt model_state_time ring radius native_valid "
    "lower_samples upper_samples lower_position upper_position "
    "h_native_plic h0 h_track_a_n h_track_a_np1 active_n active_np1 "
    "traction_active h_track_a_positive\n");
  fclose(output);
}

static void track_a_phase2c_append_accepted (
  double sample_time, int sample_iteration, double sample_dt)
{
  const CollisionState state = measure_collision_state();
  const AxisymmetricFilmSamples native_geometry =
    sample_axisymmetric_water_film(
      (double)CONTACT_SUPPORT_HALF_WIDTH);
  if (pid() != 0)
    return;
  const size_t rings = track_a_coupled.diagnostics.ring_count;
  size_t active_rings = 0, traction_active_rings = 0;
  for (size_t ring = 0; ring < rings; ring++) {
    active_rings += track_a_coupled.diagnostics.active_np1[ring] ? 1U : 0U;
    traction_active_rings +=
      track_a_coupled.diagnostics.traction_active[ring] ? 1U : 0U;
  }
  FILE *output = fopen("data/track_a_phase2c_steps.dat", "a");
  if (!output) {
    perror("data/track_a_phase2c_steps.dat");
    exit(6);
  }
  fprintf(output, "%.17g %d %.17g %.17g %.17g %lu %d %zu %zu ",
    sample_time, sample_iteration, sample_dt,
    track_a_coupled.diagnostics.time_n,
    track_a_coupled.diagnostics.time_np1,
    track_a_coupled.accepted_steps, track_a_coupled.result.active,
    active_rings, traction_active_rings);
  fprintf(output, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g ",
    state.centers[0], state.centers[1],
    state.centers[1] - state.centers[0],
    state.center_velocities[0], state.center_velocities[1],
    state.center_velocities[0] - state.center_velocities[1],
    track_a_phase2c_eq28_mask_time);
  fprintf(output, "%ld %ld %.17g %.17g %.17g %.17g ",
    track_a_phase2c_eq28_masked_faces,
    track_a_coupled.deposition.invalid_identity_faces,
    track_a_coupled.deposition.requested_lower_x_force,
    track_a_coupled.deposition.requested_upper_x_force,
    track_a_coupled.deposition.requested_lower_r_force,
    track_a_coupled.deposition.requested_upper_r_force);
  fprintf(output, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
    track_a_coupled.deposition.delivered_lower_x_force,
    track_a_coupled.deposition.delivered_upper_x_force,
    track_a_coupled.deposition.delivered_lower_r_force,
    track_a_coupled.deposition.delivered_upper_r_force,
    track_a_coupled.deposition.requested_power,
    track_a_coupled.deposition.delivered_power,
    track_a_coupled.deposition.power_residual);
  fclose(output);

  output = fopen("data/track_a_phase2c_rings.dat", "a");
  if (!output) {
    perror("data/track_a_phase2c_rings.dat");
    exit(6);
  }
  for (size_t ring = 0; ring < rings; ring++) {
    const AxisymmetricFilmSample *native = &native_geometry.bin[ring];
    const FilmRingObservation *coupling =
      &track_a_coupled.observation.rings[ring];
    const bool native_positions_available =
      native->lower_sample_count > 0 && native->upper_sample_count > 0 &&
      isfinite(native->lower_point.x) && isfinite(native->upper_point.x);
    const double h_native_plic = native_positions_available ?
      native->upper_point.x - native->lower_point.x : NAN;
    const double lower_total =
      track_a_coupled.diagnostics.lower_matched_correction_traction[ring];
    const double upper_total =
      track_a_coupled.diagnostics.upper_matched_correction_traction[ring];
    const double eq27_pressure_traction = .5*(lower_total + upper_total);
    const double eq27_shear_traction = .5*(upper_total - lower_total);
    fprintf(output, "%.17g %d %.17g %zu %.17g %.17g ",
      sample_time, sample_iteration, sample_dt, ring,
      native->radius, native->width);
    fprintf(output, "%.17g %.17g %.17g %.17g %u %u %u %d %d ",
      h_native_plic,
      track_a_coupled.diagnostics.estimated_thickness[ring],
      track_a_coupled.diagnostics.thickness_n[ring],
      track_a_coupled.diagnostics.thickness_np1[ring],
      (unsigned)track_a_coupled.diagnostics.active_n[ring],
      (unsigned)track_a_coupled.diagnostics.active_np1[ring],
      (unsigned)track_a_coupled.diagnostics.traction_active[ring],
      coupling->lower_identity_valid, coupling->upper_identity_valid);
    fprintf(output, "%.17g %.17g %.17g %.17g %.17g ",
      track_a_coupled.diagnostics.smoothed_pressure_gradient[ring],
      eq27_pressure_traction, eq27_shear_traction,
      lower_total, upper_total);
    fprintf(output, "%.17g %.17g %.17g %.17g ",
      track_a_coupled.deposition.requested_lower_x_force_by_ring[ring],
      track_a_coupled.deposition.requested_upper_x_force_by_ring[ring],
      track_a_coupled.deposition.requested_lower_r_force_by_ring[ring],
      track_a_coupled.deposition.requested_upper_r_force_by_ring[ring]);
    fprintf(output, "%.17g %.17g %.17g %.17g\n",
      track_a_coupled.deposition.delivered_lower_x_force_by_ring[ring],
      track_a_coupled.deposition.delivered_upper_x_force_by_ring[ring],
      track_a_coupled.deposition.delivered_lower_r_force_by_ring[ring],
      track_a_coupled.deposition.delivered_upper_r_force_by_ring[ring]);
  }
  fclose(output);
}

static void track_a_phase2c_append_failure (
  double sample_time, int sample_iteration, double sample_dt,
  const AxisymmetricFilmSamples *geometry)
{
  if (pid() != 0 || !geometry)
    return;
  FILE *output = fopen("data/track_a_phase2c_failure_rings.dat", "a");
  if (!output) {
    perror("data/track_a_phase2c_failure_rings.dat");
    exit(6);
  }
  const bool model_valid = track_a_phase2c_failure_snapshot.valid &&
    track_a_phase2c_failure_snapshot.ring_count ==
      AXISYMMETRIC_MATCHED_FILM_BINS;
  for (int ring = 0; ring < AXISYMMETRIC_MATCHED_FILM_BINS; ring++) {
    const AxisymmetricFilmSample *native = &geometry->bin[ring];
    double lower_position = native->valid ? native->lower_point.x : NAN;
    double upper_position = native->valid ? native->upper_point.x : NAN;
    if (ring == geometry->first_nonpositive_gap_bin) {
      lower_position = geometry->first_nonpositive_lower_position;
      upper_position = geometry->first_nonpositive_upper_position;
    }
    const double h_native_plic =
      isfinite(lower_position) && isfinite(upper_position) ?
        upper_position - lower_position : NAN;
    const double h0 = model_valid ?
      track_a_phase2c_failure_snapshot.estimated_thickness[ring] : NAN;
    const double h_track_a_n = model_valid ?
      track_a_phase2c_failure_snapshot.thickness_n[ring] : NAN;
    const double h_track_a_np1 = model_valid ?
      track_a_phase2c_failure_snapshot.thickness_np1[ring] : NAN;
    fprintf(output, "%.17g %d %.17g %.17g %d %.17g %d %d %d ",
      sample_time, sample_iteration, sample_dt,
      model_valid ? track_a_phase2c_failure_snapshot.time_np1 : NAN,
      ring, native->radius, native->valid,
      native->lower_sample_count, native->upper_sample_count);
    fprintf(output, "%.17g %.17g %.17g %.17g %.17g %.17g ",
      lower_position, upper_position, h_native_plic,
      h0, h_track_a_n, h_track_a_np1);
    fprintf(output, "%u %u %u %d\n",
      model_valid ?
        (unsigned)track_a_phase2c_failure_snapshot.active_n[ring] : 0U,
      model_valid ?
        (unsigned)track_a_phase2c_failure_snapshot.active_np1[ring] : 0U,
      model_valid ?
        (unsigned)track_a_phase2c_failure_snapshot.traction_active[ring] : 0U,
      isfinite(h_track_a_np1) && h_track_a_np1 > 0.);
  }
  fclose(output);
}

static void track_a_phase2c_capture_failure_state (void)
{
  track_a_phase2c_failure_snapshot.valid = false;
  const MusehaneFilmDiagnostics *source = &track_a_coupled.diagnostics;
  if (!source->valid ||
      source->ring_count != AXISYMMETRIC_MATCHED_FILM_BINS)
    return;
  track_a_phase2c_failure_snapshot.ring_count = source->ring_count;
  track_a_phase2c_failure_snapshot.time_np1 = source->time_np1;
  for (size_t ring = 0; ring < source->ring_count; ring++) {
    track_a_phase2c_failure_snapshot.active_n[ring] = source->active_n[ring];
    track_a_phase2c_failure_snapshot.active_np1[ring] =
      source->active_np1[ring];
    track_a_phase2c_failure_snapshot.traction_active[ring] =
      source->traction_active[ring];
    track_a_phase2c_failure_snapshot.estimated_thickness[ring] =
      source->estimated_thickness[ring];
    track_a_phase2c_failure_snapshot.thickness_n[ring] =
      source->thickness_n[ring];
    track_a_phase2c_failure_snapshot.thickness_np1[ring] =
      source->thickness_np1[ring];
  }
  track_a_phase2c_failure_snapshot.valid = true;
}

static void track_a_phase2c_write_core_failure_fixture (void)
{
  if (pid() != 0)
    return;
  MusehaneRejectedCoreInputView fixture = {0};
  if (musehane_film_model_rejected_core_input_view(
        &track_a_coupled.model_state, &fixture) != FILM_MODEL_OK ||
      !fixture.valid)
    return;
  FILE *output = fopen(
    "data/track_a_musehane_core_failure_fixture.dat", "w");
  if (!output) {
    perror("data/track_a_musehane_core_failure_fixture.dat");
    exit(6);
  }
  fprintf(output, "# musehane_core_failure_fixture_schema 1\n");
  fprintf(output,
    "# config %.17g %.17g %d %.17g %.17g %u\n",
    fixture.config.water_density, fixture.config.water_viscosity,
    (int)fixture.config.discretization,
    fixture.config.nonlinear_relative_tolerance,
    fixture.config.inventory_relative_tolerance,
    fixture.config.maximum_iterations);
  fprintf(output, "# step %.17g %.17g %zu\n",
    fixture.time, fixture.dt, fixture.mesh.cell_count);
  fprintf(output,
    "# cell radius width measure active estimated old divergence gradient_cell\n");
  for (size_t ring = 0; ring < fixture.mesh.cell_count; ring++)
    fprintf(output, "cell %zu %.17g %.17g %.17g %u %.17g %.17g %.17g %.17g\n",
      ring, fixture.mesh.cell_radius[ring], fixture.mesh.cell_width[ring],
      fixture.mesh.cell_measure[ring], (unsigned)fixture.active[ring],
      fixture.estimated_thickness[ring], fixture.old_thickness[ring],
      fixture.mean_surface_velocity_divergence[ring],
      fixture.decontaminated_pressure_gradient_cell[ring]);
  fprintf(output, "# face measure gradient_face\n");
  for (size_t face = 0; face <= fixture.mesh.cell_count; face++)
    fprintf(output, "face %zu %.17g %.17g\n", face,
      fixture.mesh.face_measure[face],
      fixture.decontaminated_pressure_gradient_face[face]);
  fclose(output);
}
#endif

static bool track_a_coupled_sidecar_path (
  const char *dump_path, char *sidecar_path, size_t capacity)
{
  if (!dump_path || !dump_path[0] || !sidecar_path || capacity == 0)
    return false;
  const int count = snprintf(sidecar_path, capacity, "%s.track-a.bin",
    dump_path);
  return count > 0 && (size_t)count < capacity;
}

static bool track_a_coupled_write_restart_metadata (
  const char *dump_path, double dump_time, int dump_iteration)
{
  char metadata_path[512], temporary_path[516];
  int count = snprintf(metadata_path, sizeof metadata_path,
    "%s.track-a.json", dump_path);
  if (count <= 0 || (size_t)count >= sizeof metadata_path)
    return false;
  count = snprintf(temporary_path, sizeof temporary_path, "%s.tmp",
    metadata_path);
  if (count <= 0 || (size_t)count >= sizeof temporary_path)
    return false;
  FILE *output = fopen(temporary_path, "w");
  if (!output)
    return false;
  fprintf(output,
    "{\n"
    "  \"schema\": \"track-a-collision-restart-v2\",\n"
    "  \"time\": %.17g,\n"
    "  \"iteration\": %d,\n"
    "  \"accepted_steps\": %lu,\n"
    "  \"model_payload_size\": %zu,\n"
    "  \"model_payload_checksum\": \"%llu\",\n"
    "  \"observation_release_latch_mask\": \"%llu\",\n"
    "  \"initial_lower_gas_volume\": %.17g,\n"
    "  \"initial_lower_oil_volume\": %.17g,\n"
    "  \"initial_upper_gas_volume\": %.17g,\n"
    "  \"initial_upper_oil_volume\": %.17g,\n"
    "  \"anchor_status\": %d,\n"
    "  \"anchor_integral_error\": %.17g,\n"
    "  \"anchor_commanded_acceleration\": %.17g,\n"
    "  \"anchor_release_time\": %.17g,\n"
    "  \"anchor_release_gap\": %.17g,\n"
    "  \"anchor_last_reaction_force\": %.17g,\n"
    "  \"anchor_last_power\": %.17g,\n"
    "  \"anchor_cumulative_work\": %.17g,\n"
    "  \"anchor_saturation_count\": %lu,\n"
    "  \"anchor_release_threshold\": %.17g,\n"
    "  \"anchor_proportional_gain\": %.17g,\n"
    "  \"anchor_integral_gain\": %.17g,\n"
    "  \"anchor_acceleration_limit\": %.17g,\n"
    "  \"anchor_feedforward_acceleration\": %.17g,\n"
    "  \"anchor_require_approach\": %d\n"
    "}\n",
    dump_time, dump_iteration, track_a_coupled.accepted_steps,
    track_a_coupled.checkpoint_audit.payload_size,
    (unsigned long long)track_a_coupled.checkpoint_audit.payload_checksum,
    (unsigned long long)track_a_coupled.observation_release_latch_mask,
    initial_gas_volumes[0], initial_oil_volumes[0],
    initial_gas_volumes[1], initial_oil_volumes[1],
    anchor_state.status, anchor_state.integral_error,
    anchor_state.commanded_acceleration, anchor_state.release_time,
    anchor_state.release_gap, anchor_state.last_reaction_force,
    anchor_state.last_power, anchor_state.cumulative_work,
    anchor_state.saturation_count, anchor_config.release_gap,
    anchor_config.proportional_gain, anchor_config.integral_gain,
    anchor_config.acceleration_limit, anchor_config.feedforward_acceleration,
    anchor_config.require_approach);
  const bool closed = fclose(output) == 0;
  if (!closed || rename(temporary_path, metadata_path) != 0) {
    remove(temporary_path);
    return false;
  }
  return true;
}

static bool track_a_coupled_read_restart_release_latch (
  const char *dump_path, uint64_t *latch_mask)
{
  if (!dump_path || !dump_path[0] || !latch_mask)
    return false;
  char metadata_path[512];
  const int count = snprintf(metadata_path, sizeof metadata_path,
    "%s.track-a.json", dump_path);
  if (count <= 0 || (size_t)count >= sizeof metadata_path)
    return false;
  FILE *input = fopen(metadata_path, "r");
  if (!input)
    return false;
  bool schema_v1 = false, schema_v2 = false, found_latch = false;
  unsigned long long parsed_latch = 0ULL;
  char line[512];
  while (fgets(line, sizeof line, input)) {
    if (strstr(line, "\"schema\": \"track-a-collision-restart-v1\""))
      schema_v1 = true;
    if (strstr(line, "\"schema\": \"track-a-collision-restart-v2\""))
      schema_v2 = true;
    if (sscanf(line,
        " \"observation_release_latch_mask\": \"%llu\"",
        &parsed_latch) == 1)
      found_latch = true;
  }
  const bool closed = fclose(input) == 0;
  if (!closed || schema_v1 == schema_v2 || (schema_v2 && !found_latch))
    return false;
  *latch_mask = found_latch ? (uint64_t)parsed_latch : UINT64_C(0);
  return true;
}

static bool collision_traction_extension_checkpoint (
  const char *dump_path, double dump_time, int dump_iteration)
{
  if (!track_a_coupled.initialized)
    return true;
  if (track_a_coupled.transaction.open || dump_iteration < 0 ||
      !isfinite(dump_time) || dump_time < 0.) {
    track_a_coupled.status = TRACK_A_COUPLED_CHECKPOINT_FILE_FAILED;
    return false;
  }
  char sidecar_path[512];
  if (!track_a_coupled_sidecar_path(dump_path, sidecar_path,
      sizeof sidecar_path)) {
    track_a_coupled.status = TRACK_A_COUPLED_CHECKPOINT_FILE_FAILED;
    return false;
  }
  int local_success = 1;
  if (pid() == 0) {
    local_success = musehane_checkpoint_file_write(
      sidecar_path, &track_a_coupled.model_state, dump_time,
      (unsigned long)dump_iteration, track_a_coupled.accepted_steps,
      &track_a_coupled.checkpoint_audit) == MUSEHANE_CHECKPOINT_FILE_OK;
    if (local_success)
      local_success = track_a_coupled_write_restart_metadata(
        dump_path, dump_time, dump_iteration);
    if (local_success) {
      FILE *output = fopen("data/track_a_musehane_checkpoints.dat", "a");
      if (!output)
        local_success = 0;
      else {
        fprintf(output, "%.17g %d %lu %zu %llu %s\n",
          dump_time, dump_iteration, track_a_coupled.accepted_steps,
          track_a_coupled.checkpoint_audit.payload_size,
          (unsigned long long)
            track_a_coupled.checkpoint_audit.payload_checksum,
          sidecar_path);
        if (fclose(output) != 0)
          local_success = 0;
      }
    }
  }
#if _MPI
  MPI_Bcast(&local_success, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  if (!local_success) {
    track_a_coupled.status = TRACK_A_COUPLED_CHECKPOINT_FILE_FAILED;
    return false;
  }
  return true;
}

#if COLLISION_FULL_STATE_RESTART
static bool collision_traction_extension_restore (
  const char *dump_path, double dump_time, int dump_iteration)
{
  if (dump_iteration < 0 || !isfinite(dump_time) || dump_time < 0. ||
      track_a_coupled.initialized ||
      track_a_coupled.status != TRACK_A_COUPLED_UNINITIALIZED)
    return false;
  if (!track_a_coupled_start())
    return false;
  char sidecar_path[512];
  if (!track_a_coupled_sidecar_path(dump_path, sidecar_path,
      sizeof sidecar_path)) {
    track_a_coupled.status = TRACK_A_COUPLED_RESTORE_FILE_FAILED;
    return false;
  }
  unsigned long accepted_steps = 0;
  uint64_t release_latch_mask = 0U;
  const MusehaneCheckpointFileStatus restored =
    musehane_checkpoint_file_restore(
      sidecar_path, &track_a_coupled.config, dump_time,
      (unsigned long)dump_iteration, &track_a_coupled.model_state,
      &accepted_steps, &track_a_coupled.checkpoint_audit);
  int local_success = restored == MUSEHANE_CHECKPOINT_FILE_OK &&
    track_a_coupled_read_restart_release_latch(
      dump_path, &release_latch_mask);
#if _MPI
  int global_success = 0;
  MPI_Allreduce(&local_success, &global_success, 1, MPI_INT, MPI_MIN,
    MPI_COMM_WORLD);
  local_success = global_success;
#endif
  if (!local_success) {
    if (restored == MUSEHANE_CHECKPOINT_FILE_OK)
      film_model_destroy(&track_a_coupled.model_state);
    track_a_coupled.status = TRACK_A_COUPLED_RESTORE_FILE_FAILED;
    return false;
  }
  track_a_coupled.accepted_steps = accepted_steps;
  track_a_coupled.observation_release_latch_mask = release_latch_mask;
  track_a_coupled.restart_clock_alignment_pending = true;
  track_a_coupled.restart_clock_rebased = false;
  track_a_coupled.restart_clock_source_time = dump_time;
  track_a_coupled.restart_clock_gap = NAN;
  track_a_coupled.restart_clock_maximum_gap =
    (double)TRACK_A_MUSEHANE_RESTART_CLOCK_GAP_FACTOR*(double)DT_MAX;
  track_a_coupled.initialized = true;
  track_a_coupled.status = TRACK_A_COUPLED_READY;
  return true;
}
#endif

static bool track_a_coupled_start (void)
{
  collision_traction_extension_prepare_common_state_checkpoint();
  track_a_coupled.config = film_model_musehane_config();
  track_a_coupled.config.activation_cells =
    (double)TRACK_A_MUSEHANE_ACTIVATION_CELLS;
  track_a_coupled.config.pressure_driver.smoothing_length_squared =
    (double)TRACK_A_MUSEHANE_SMOOTHING_LENGTH_SQUARED;
  track_a_coupled.config.traction_policy =
    (MusehaneTractionPolicy)TRACK_A_MUSEHANE_TRACTION_POLICY;
  track_a_coupled.status = TRACK_A_COUPLED_READY;
  track_a_coupled_write_header();
#if TRACK_A_PHASE2C_AUDIT
  track_a_phase2c_write_headers();
#endif
  return true;
}

static void track_a_coupled_append (double sample_time, int iteration,
  double sample_dt)
{
  if (pid() != 0)
    return;
  double minimum_gap = HUGE;
  for (size_t ring = 0; ring < track_a_coupled.observation.view.ring_count;
       ring++)
    minimum_gap = min(minimum_gap,
      track_a_coupled.observation.rings[ring].gap);
  FILE *output = fopen("data/track_a_musehane_two_way.dat", "a");
  if (!output) {
    perror("data/track_a_musehane_two_way.dat");
    exit(6);
  }
  fprintf(output,
    "%.17g %d %.17g %d %lu %d %zu %.17g %.17g %.17g %.17g %.17g "
    "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g "
    "%.17g %.17g %.17g %.17g %.17g "
    "%.17g %.17g %d %.17g %d %d %u %d %d %.17g %.17g %u %.17g "
    "%llu %u\n",
    sample_time, iteration, sample_dt, track_a_coupled.status,
    track_a_coupled.accepted_steps, track_a_coupled.result.active,
    track_a_coupled.diagnostics.matched_corrected_ring_count,
    minimum_gap, track_a_coupled.result.minimum_thickness,
    track_a_coupled.result.step_work,
    track_a_coupled.result.cumulative_work,
    track_a_coupled.deposition.requested_power,
    track_a_coupled.deposition.delivered_power,
    track_a_coupled.diagnostics.equation27_pressure_correction_power,
    track_a_coupled.diagnostics.equation27_shear_correction_power,
    track_a_coupled.diagnostics.poiseuille_dissipation_rate,
    track_a_coupled.diagnostics.cumulative_poiseuille_dissipation,
    track_a_coupled.diagnostics.equation27_total_resolved_reaction_force,
    track_a_coupled.deposition.power_residual,
    track_a_coupled.deposition.requested_lower_power,
    track_a_coupled.deposition.requested_upper_power,
    track_a_coupled.deposition.delivered_lower_power,
    track_a_coupled.deposition.delivered_upper_power,
    track_a_coupled.deposition.maximum_ring_power_residual,
    track_a_coupled.deposition.maximum_interface_velocity_mismatch,
    track_a_coupled.deposition.maximum_component_force_residual,
    track_a_coupled.deposition.maximum_bin_force_residual,
    track_a_coupled.diagnostics.maximum_pressure_driver_residual,
    track_a_coupled.diagnostics.maximum_traction_constitutive_residual,
    track_a_coupled.velocity_observation.status,
    track_a_coupled.velocity_observation.maximum_plic_velocity_adjustment,
    track_a_coupled.velocity_observation.geometric_fallback_components,
    track_a_coupled.velocity_observation.first_geometric_fallback_bin,
    track_a_coupled.velocity_observation.first_geometric_fallback_mask,
    track_a_coupled.deposition.geometric_fallback_components,
    track_a_coupled.restart_clock_rebased,
    track_a_coupled.restart_clock_gap,
    track_a_coupled.restart_clock_maximum_gap,
    track_a_coupled.diagnostics.resolved_release_count,
    track_a_coupled.diagnostics.resolved_release_inventory_change,
    (unsigned long long)track_a_coupled.observation_release_latch_mask,
    track_a_observation_release_latch_count(
      track_a_coupled.observation_release_latch_mask));
  fclose(output);
  output = fopen("data/track_a_musehane_substeps.dat", "a");
  if (!output) {
    perror("data/track_a_musehane_substeps.dat");
    exit(6);
  }
  fprintf(output, "%.17g %d %.17g %lu %u\n",
    sample_time, iteration, sample_dt, track_a_coupled.accepted_steps,
    track_a_coupled.diagnostics.substep_count);
  fclose(output);
}

static bool collision_traction_extension_apply (
  face vector acceleration_field,
  double sample_time,
  double sample_dt,
  int sample_iteration)
{
#if TRACK_A_PHASE2C_AUDIT
  track_a_phase2c_failure_snapshot.valid = false;
#endif
  if (collision_terminal_event_time(sample_time, (double)END_TIME))
    return true;
  if (track_a_coupled.status == TRACK_A_COUPLED_UNINITIALIZED &&
      !track_a_coupled_start())
    return false;
  const AxisymmetricMusehaneNativeFieldStatus native_status =
    build_axisymmetric_musehane_native_fields(
      p, track_a_coupled_pressure_gradient,
      track_a_coupled_capillary_force_density,
      &track_a_coupled.native_audit);
  if (native_status != AXISYMMETRIC_MUSEHANE_NATIVE_FIELD_OK) {
    track_a_coupled.status = TRACK_A_COUPLED_NATIVE_FIELD_FAILED;
    return false;
  }
  TrackASubgridObservationState subgrid_observation = {0};
  if (track_a_coupled.initialized) {
    const FilmModelStatus diagnostics_status =
      musehane_film_model_diagnostics(&track_a_coupled.model_state,
        &track_a_coupled.diagnostics);
    if (diagnostics_status == FILM_MODEL_OK &&
        track_a_coupled.diagnostics.valid) {
      subgrid_observation = (TrackASubgridObservationState) {
        .valid = true,
        .ring_count = (int)track_a_coupled.diagnostics.ring_count,
        .active_n = track_a_coupled.diagnostics.active_n,
        .active_np1 = track_a_coupled.diagnostics.active_np1,
        .thickness = track_a_coupled.diagnostics.thickness_np1,
        .release_latch_mask =
          track_a_coupled.observation_release_latch_mask
      };
    }
    else if (track_a_coupled.restart_clock_alignment_pending) {
      MusehaneFilmRestoredStateView restored_state = {0};
      if (musehane_film_model_restored_state_view(
          &track_a_coupled.model_state, &restored_state) != FILM_MODEL_OK ||
          !restored_state.valid) {
        track_a_coupled.status = TRACK_A_COUPLED_MODEL_FAILED;
        return false;
      }
      subgrid_observation = (TrackASubgridObservationState) {
        .valid = true,
        .ring_count = (int)restored_state.ring_count,
        .active_n = restored_state.active,
        .active_np1 = restored_state.active,
        .thickness = restored_state.thickness,
        .release_latch_mask =
          track_a_coupled.observation_release_latch_mask
      };
    }
    else {
      track_a_coupled.status = TRACK_A_COUPLED_MODEL_FAILED;
      return false;
    }
  }
  track_a_coupled.observation_status =
    sample_axisymmetric_musehane_observation(
      (double)CONTACT_SUPPORT_HALF_WIDTH, p,
      track_a_coupled_pressure_gradient,
      track_a_coupled_capillary_force_density,
      sample_time, sample_dt, (unsigned long)sample_iteration,
      anchor_state.status == CAPILLARY_ANCHOR_RELEASED, 1U, 2U,
      subgrid_observation.valid ? &subgrid_observation : NULL,
      &track_a_coupled.observation);
  if (track_a_coupled.observation_status !=
      AXISYMMETRIC_MUSEHANE_OBSERVATION_OK) {
    track_a_coupled.status = TRACK_A_COUPLED_OBSERVATION_FAILED;
    return false;
  }
  track_a_coupled.velocity_observation =
    axisymmetric_musehane_observe_support_velocity_axi(
      u, lower_envelope_fraction, upper_envelope_fraction,
      track_a_coupled.observation.rings,
      track_a_coupled.observation.view.ring_count,
      (double)CONTACT_SUPPORT_HALF_WIDTH);
  if (!track_a_coupled.velocity_observation.observed ||
      track_a_coupled.velocity_observation.status !=
        AXISYMMETRIC_MUSEHANE_VELOCITY_OK) {
    track_a_coupled.status = TRACK_A_COUPLED_VELOCITY_FAILED;
    return false;
  }
  if (!track_a_coupled.initialized) {
    if (film_model_initialize(&track_a_coupled.model_state,
        &track_a_coupled.config, &track_a_coupled.observation.view) !=
        FILM_MODEL_OK) {
      track_a_coupled.status = TRACK_A_COUPLED_MODEL_FAILED;
      return false;
    }
    track_a_coupled.initialized = true;
  }
  if (track_a_coupled.restart_clock_alignment_pending) {
    const double clock_gap = sample_time -
      track_a_coupled.restart_clock_source_time;
    if (!isfinite(clock_gap) || clock_gap < 0. ||
        film_model_rebase_restart_clock(&track_a_coupled.model_state,
          sample_time, track_a_coupled.restart_clock_maximum_gap) !=
          FILM_MODEL_OK) {
      track_a_coupled.status = TRACK_A_COUPLED_RESTART_CLOCK_FAILED;
      return false;
    }
    track_a_coupled.restart_clock_gap = clock_gap;
    track_a_coupled.restart_clock_rebased = true;
    track_a_coupled.restart_clock_alignment_pending = false;
  }
  if (musehane_model_transaction_begin(&track_a_coupled.transaction,
      &track_a_coupled.model_state) != MUSEHANE_MODEL_TRANSACTION_OK) {
    track_a_coupled.status = TRACK_A_COUPLED_CHECKPOINT_FAILED;
    return false;
  }
  track_a_coupled.result = film_model_advance(
    &track_a_coupled.model_state, &track_a_coupled.config,
    &track_a_coupled.observation.view, sample_dt);
  if (track_a_coupled.result.status != FILM_MODEL_OK ||
      !track_a_coupled.result.commit_allowed ||
      musehane_film_model_diagnostics(&track_a_coupled.model_state,
        &track_a_coupled.diagnostics) != FILM_MODEL_OK) {
#if TRACK_A_PHASE2C_AUDIT
    track_a_phase2c_write_core_failure_fixture();
    track_a_phase2c_capture_failure_state();
#endif
    return track_a_coupled_rollback(TRACK_A_COUPLED_MODEL_FAILED);
  }
  track_a_coupled.deposition = axisymmetric_musehane_apply_traction_axi(
    acceleration_field, dual_compound_alpha, u,
    lower_envelope_fraction, upper_envelope_fraction,
    &track_a_coupled.result, &track_a_coupled.observation.view,
    (double)CONTACT_SUPPORT_HALF_WIDTH);
  if (!track_a_coupled.deposition.applied ||
      track_a_coupled.deposition.status !=
        AXISYMMETRIC_MUSEHANE_TRACTION_OK) {
#if TRACK_A_PHASE2C_AUDIT
    track_a_phase2c_capture_failure_state();
#endif
    return track_a_coupled_rollback(TRACK_A_COUPLED_DEPOSITION_FAILED);
  }
  if (musehane_model_transaction_commit(&track_a_coupled.transaction) !=
      MUSEHANE_MODEL_TRANSACTION_OK) {
    track_a_coupled.status = TRACK_A_COUPLED_COMMIT_FAILED;
    return false;
  }
  for (size_t ring = 0; ring < track_a_coupled.diagnostics.ring_count;
       ring++) {
    const FilmRingObservation *observation =
      &track_a_coupled.observation.rings[ring];
    const bool native_resolved =
      track_a_coupled.observation.native_valid[ring] != 0U &&
      observation->gap >= track_a_coupled.config.activation_cells*
        observation->local_delta;
    track_a_coupled.observation_release_latch_mask =
      track_a_update_observation_release_latch(
        track_a_coupled.observation_release_latch_mask, ring,
        native_resolved,
        track_a_coupled.diagnostics.active_n[ring] != 0U,
        track_a_coupled.diagnostics.active_np1[ring] != 0U);
  }
  track_a_coupled.accepted_steps++;
  track_a_coupled.status = TRACK_A_COUPLED_OK;
#if TRACK_A_PHASE2C_AUDIT
  track_a_phase2c_append_accepted(
    sample_time, sample_iteration, sample_dt);
#endif
  track_a_coupled_append(sample_time, sample_iteration, sample_dt);
  return true;
}

static void collision_traction_extension_report_failure (
  double sample_time, double sample_dt, int sample_iteration)
{
  /* Geometry sampling contains MPI reductions and interpolation, so every
     rank must enter it before rank-zero-only formatting.  This is a read-only
     failure snapshot; it does not alter the accepted Track A state or force. */
  const bool geometry_sampled = track_a_coupled.status ==
    TRACK_A_COUPLED_OBSERVATION_FAILED;
  AxisymmetricFilmSamples failed_geometry = {0};
  axisymmetric_film_samples_set_widths(
    &failed_geometry, (double)CONTACT_SUPPORT_HALF_WIDTH);
  if (geometry_sampled)
    failed_geometry = sample_axisymmetric_water_film(
      (double)CONTACT_SUPPORT_HALF_WIDTH);
#if TRACK_A_PHASE2C_AUDIT
  if (!track_a_phase2c_failure_snapshot.valid)
    track_a_phase2c_capture_failure_state();
  track_a_phase2c_append_failure(
    sample_time, sample_iteration, sample_dt, &failed_geometry);
#endif
  if (pid() != 0)
    return;
  char geometry_validity[AXISYMMETRIC_MATCHED_FILM_BINS + 1];
  for (int bin = 0; bin < AXISYMMETRIC_MATCHED_FILM_BINS; bin++)
    geometry_validity[bin] = geometry_sampled ?
      (failed_geometry.bin[bin].valid ? '1' : '0') : '?';
  geometry_validity[AXISYMMETRIC_MATCHED_FILM_BINS] = '\0';
  const double first_nonpositive_gap =
    isfinite(failed_geometry.first_nonpositive_lower_position) &&
    isfinite(failed_geometry.first_nonpositive_upper_position) ?
      failed_geometry.first_nonpositive_upper_position -
      failed_geometry.first_nonpositive_lower_position : NAN;
  fprintf(stderr,
    "[TRACK-A-COUPLING-REJECT] t=%.17g i=%d dt=%.17g status=%d "
    "initialized=%d accepted_steps=%lu result_status=%d "
    "result_stage=%d velocity_status=%d "
    "first_missing_support_bin=%d first_missing_support_mask=%u "
    "first_missing_lower_x_support=%.17g "
    "first_missing_lower_r_support=%.17g "
    "first_missing_upper_x_support=%.17g "
    "first_missing_upper_r_support=%.17g "
    "first_missing_lower_x_geometric_support=%.17g "
    "first_missing_lower_r_geometric_support=%.17g "
    "first_missing_upper_x_geometric_support=%.17g "
    "first_missing_upper_r_geometric_support=%.17g "
    "first_missing_lower_x_max_identity=%.17g "
    "first_missing_lower_r_max_identity=%.17g "
    "first_missing_upper_x_max_identity=%.17g "
    "first_missing_upper_r_max_identity=%.17g deposition_status=%d "
    "plan_status=%d transaction_open=%d result_reason=\"%s\" "
    "failure_diagnostic_valid=%d failure_kind=%d "
    "failure_iteration=%u failure_line_search_backtracks=%u "
    "failure_initial_merit=%.17g failure_final_merit=%.17g "
    "failure_cell=%zu failure_numerator=%.17g failure_denominator=%.17g "
    "failure_old_thickness=%.17g failure_current_thickness=%.17g "
    "failure_flux_in=%.17g failure_flux_out=%.17g "
    "failure_surface_velocity_divergence=%.17g "
    "failure_equation_relative_residual=%.17g "
    "observation_status=%d observation_bins=%d "
    "observation_valid_bins=%d observation_interpolated_bins=%d "
    "observation_valid_area_fraction=%.17g "
    "geometry_sampled=%d geometry_input_ready=%d "
    "geometry_valid_bins=%d geometry_valid_area_fraction=%.17g "
    "geometry_nonpositive_gap_bins=%d "
    "geometry_invalid_plic_segments=%d "
    "geometry_degenerate_plic_facets=%d "
    "geometry_first_nonpositive_bin=%d "
    "geometry_first_nonpositive_lower_samples=%d "
    "geometry_first_nonpositive_upper_samples=%d "
    "geometry_first_nonpositive_lower_position=%.17g "
    "geometry_first_nonpositive_upper_position=%.17g "
    "geometry_first_nonpositive_gap=%.17g "
    "geometry_first_nonpositive_lower_delta=%.17g "
    "geometry_first_nonpositive_upper_delta=%.17g "
    "geometry_validity=%s "
    "traction_audit_valid=%d traction_status=%d traction_cell=%zu "
    "pressure_power=%.17g shear_power=%.17g direct_power=%.17g "
    "total_power=%.17g power_identity_residual=%.17g "
    "step_work=%.17g cumulative_work=%.17g\n",
    sample_time, sample_iteration, sample_dt, track_a_coupled.status,
    track_a_coupled.initialized, track_a_coupled.accepted_steps,
    track_a_coupled.result.status, track_a_coupled.result.completed_stage,
    track_a_coupled.velocity_observation.status,
    track_a_coupled.velocity_observation.first_missing_support_bin,
    track_a_coupled.velocity_observation.first_missing_support_mask,
    track_a_coupled.velocity_observation.first_missing_lower_x_support,
    track_a_coupled.velocity_observation.first_missing_lower_r_support,
    track_a_coupled.velocity_observation.first_missing_upper_x_support,
    track_a_coupled.velocity_observation.first_missing_upper_r_support,
    track_a_coupled.velocity_observation.
      first_missing_lower_x_geometric_support,
    track_a_coupled.velocity_observation.
      first_missing_lower_r_geometric_support,
    track_a_coupled.velocity_observation.
      first_missing_upper_x_geometric_support,
    track_a_coupled.velocity_observation.
      first_missing_upper_r_geometric_support,
    track_a_coupled.velocity_observation.first_missing_lower_x_max_identity,
    track_a_coupled.velocity_observation.first_missing_lower_r_max_identity,
    track_a_coupled.velocity_observation.first_missing_upper_x_max_identity,
    track_a_coupled.velocity_observation.first_missing_upper_r_max_identity,
    track_a_coupled.deposition.status,
    track_a_coupled.deposition.plan_status,
    track_a_coupled.transaction.open,
    track_a_coupled.result.reason ? track_a_coupled.result.reason : "(none)",
    track_a_coupled.result.failure_diagnostic_valid,
    track_a_coupled.result.failure_kind,
    track_a_coupled.result.failure_iteration,
    track_a_coupled.result.failure_line_search_backtracks,
    track_a_coupled.result.failure_initial_merit,
    track_a_coupled.result.failure_final_merit,
    track_a_coupled.result.failure_cell,
    track_a_coupled.result.failure_numerator,
    track_a_coupled.result.failure_denominator,
    track_a_coupled.result.failure_old_thickness,
    track_a_coupled.result.failure_current_thickness,
    track_a_coupled.result.failure_flux_in,
    track_a_coupled.result.failure_flux_out,
    track_a_coupled.result.failure_surface_velocity_divergence,
    track_a_coupled.result.failure_equation_relative_residual,
    track_a_coupled.observation_status,
    track_a_coupled.observation.sampled_bins,
    track_a_coupled.observation.source_valid_bins,
    track_a_coupled.observation.interpolated_bins,
    track_a_coupled.observation.source_valid_area_fraction,
    geometry_sampled, failed_geometry.input_ready,
    failed_geometry.valid_bins, failed_geometry.valid_area_fraction,
    failed_geometry.nonpositive_gap_bins,
    failed_geometry.invalid_plic_segments,
    failed_geometry.skipped_degenerate_plic_facets,
    failed_geometry.first_nonpositive_gap_bin,
    failed_geometry.first_nonpositive_lower_samples,
    failed_geometry.first_nonpositive_upper_samples,
    failed_geometry.first_nonpositive_lower_position,
    failed_geometry.first_nonpositive_upper_position,
    first_nonpositive_gap,
    failed_geometry.first_nonpositive_lower_delta,
    failed_geometry.first_nonpositive_upper_delta,
    geometry_validity,
    track_a_coupled.result.failure_traction_diagnostic_valid,
    track_a_coupled.result.failure_traction_status,
    track_a_coupled.result.failure_traction_cell,
    track_a_coupled.result.failure_traction_pressure_power,
    track_a_coupled.result.failure_traction_shear_power,
    track_a_coupled.result.failure_traction_direct_power,
    track_a_coupled.result.failure_traction_total_power,
    track_a_coupled.result.failure_traction_power_identity_residual,
    track_a_coupled.result.failure_traction_step_work,
    track_a_coupled.result.failure_traction_cumulative_work);
  fflush(stderr);
}

static void collision_traction_extension_finish (
  double finish_time, int finish_iteration)
{
  if (pid() == 0)
    fprintf(stderr,
      "# TRACK_A_TWO_WAY_FINISH t=%.17g i=%d accepted_steps=%lu\n",
      finish_time, finish_iteration, track_a_coupled.accepted_steps);
  if (track_a_coupled.initialized)
    film_model_destroy(&track_a_coupled.model_state);
  musehane_model_transaction_destroy(&track_a_coupled.transaction);
  track_a_coupled.initialized = false;
}

#endif
