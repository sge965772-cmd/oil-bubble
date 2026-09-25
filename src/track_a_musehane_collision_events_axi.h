#ifndef TRACK_A_MUSEHANE_COLLISION_EVENTS_AXI_H
#define TRACK_A_MUSEHANE_COLLISION_EVENTS_AXI_H

#include "track_a_musehane_diagnostic_axi.h"
#include "track_a_collision_clock.h"

#ifndef TRACK_A_MUSEHANE_ACTIVATION_CELLS
# define TRACK_A_MUSEHANE_ACTIVATION_CELLS 7.
#endif
#ifndef TRACK_A_MUSEHANE_SMOOTHING_LENGTH_SQUARED
# define TRACK_A_MUSEHANE_SMOOTHING_LENGTH_SQUARED 0.
#endif

static vector track_a_musehane_pressure_gradient[];
static vector track_a_musehane_capillary_force_density[];
static TrackAMusehaneDiagnosticState track_a_musehane_state =
  TRACK_A_MUSEHANE_DIAGNOSTIC_STATE_INITIALIZER;

static void track_a_musehane_write_header (void)
{
  if (pid() != 0)
    return;
  FILE *output = fopen("data/track_a_musehane_diagnostic.dat", "w");
  if (!output) {
    perror("data/track_a_musehane_diagnostic.dat");
    exit(6);
  }
  fprintf(output,
    "# t i dt status accepted_steps anchor_released active sampled_bins "
    "source_valid_bins interpolated_bins valid_area_fraction "
    "minimum_geometric_gap minimum_model_thickness "
    "model_time_n model_time_np1 "
    "max_gradp max_csf max_decontaminated_gradient pressure_residual "
    "traction_residual poiseuille_dissipation cumulative_dissipation\n");
  fclose(output);
}

static void track_a_musehane_append_diagnostic (
  double output_time, int output_iteration, double output_dt)
{
  if (pid() != 0)
    return;
  double minimum_geometric_gap = HUGE;
  for (size_t ring = 0;
       ring < track_a_musehane_state.latest_observation.view.ring_count;
       ring++)
    minimum_geometric_gap = min(minimum_geometric_gap,
      track_a_musehane_state.latest_observation.rings[ring].gap);
  FILE *output = fopen("data/track_a_musehane_diagnostic.dat", "a");
  if (!output) {
    perror("data/track_a_musehane_diagnostic.dat");
    exit(6);
  }
  fprintf(output,
    "%.17g %d %.17g %d %lu %d %d %d %d %d %.17g %.17g %.17g %.17g %.17g "
    "%.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
    output_time, output_iteration, output_dt,
    track_a_musehane_state.status,
    track_a_musehane_state.accepted_steps,
    anchor_state.status == CAPILLARY_ANCHOR_RELEASED,
    track_a_musehane_state.latest_result.active,
    track_a_musehane_state.latest_observation.sampled_bins,
    track_a_musehane_state.latest_observation.source_valid_bins,
    track_a_musehane_state.latest_observation.interpolated_bins,
    track_a_musehane_state.latest_observation.source_valid_area_fraction,
    minimum_geometric_gap,
    track_a_musehane_state.latest_result.minimum_thickness,
    track_a_musehane_state.latest_diagnostics.time_n,
    track_a_musehane_state.latest_diagnostics.time_np1,
    track_a_musehane_state.latest_native_audit.maximum_pressure_gradient,
    track_a_musehane_state.latest_native_audit.
      maximum_capillary_force_density,
    track_a_musehane_state.latest_observation.assembly_audit.
      maximum_decontaminated_pressure_gradient,
    track_a_musehane_state.latest_diagnostics.
      maximum_pressure_driver_residual,
    track_a_musehane_state.latest_diagnostics.
      maximum_traction_constitutive_residual,
    track_a_musehane_state.latest_diagnostics.
      poiseuille_dissipation_rate,
    track_a_musehane_state.latest_diagnostics.
      cumulative_poiseuille_dissipation);
  fclose(output);
}

static bool track_a_musehane_start (void)
{
  track_a_musehane_pressure_gradient.x.nodump = true;
  track_a_musehane_pressure_gradient.y.nodump = true;
  track_a_musehane_capillary_force_density.x.nodump = true;
  track_a_musehane_capillary_force_density.y.nodump = true;
  FilmModelConfig config = film_model_musehane_config();
  config.activation_cells = (double)TRACK_A_MUSEHANE_ACTIVATION_CELLS;
  config.pressure_driver.smoothing_length_squared =
    (double)TRACK_A_MUSEHANE_SMOOTHING_LENGTH_SQUARED;
  const TrackAMusehaneDiagnosticStatus status =
    track_a_musehane_diagnostic_initialize(
      &track_a_musehane_state, &config);
  if (status != TRACK_A_MUSEHANE_DIAGNOSTIC_READY) {
    if (pid() == 0)
      fprintf(stderr, "[TRACK-A-REJECT] initialization status=%d\n", status);
    return false;
  }
  track_a_musehane_write_header();
  return true;
}

/* Declared after projection/snapshot events and before adapt in the collision
 * source.  This observer does not write acceleration, velocity, pressure or
 * phase fields. */
event track_a_musehane_observe (i++)
{
  if (collision_terminal_event_time(t, (double)END_TIME))
    return 0;
  double interval_start = 0.;
  if (!track_a_completed_step_interval_start(t, dt, &interval_start)) {
    if (pid() == 0)
      fprintf(stderr,
        "[TRACK-A-REJECT] invalid completed-step clock t=%.17g i=%d "
        "dt=%.17g\n", t, i, dt);
    exit(6);
  }
  if (track_a_musehane_state.status ==
        TRACK_A_MUSEHANE_DIAGNOSTIC_UNINITIALIZED &&
      !track_a_musehane_start())
    exit(6);
  const TrackAMusehaneDiagnosticStatus status =
    track_a_musehane_diagnostic_step(
      &track_a_musehane_state, p,
      track_a_musehane_pressure_gradient,
      track_a_musehane_capillary_force_density,
      (double)CONTACT_SUPPORT_HALF_WIDTH, interval_start, dt, i,
      anchor_state.status == CAPILLARY_ANCHOR_RELEASED, 1U, 2U);
  if (status != TRACK_A_MUSEHANE_DIAGNOSTIC_OK) {
    AxisymmetricFilmSamples failed_geometry =
      sample_axisymmetric_water_film(
        (double)CONTACT_SUPPORT_HALF_WIDTH);
    if (pid() == 0) {
      char validity[AXISYMMETRIC_MATCHED_FILM_BINS + 1];
      for (int bin = 0; bin < AXISYMMETRIC_MATCHED_FILM_BINS; bin++)
        validity[bin] = failed_geometry.bin[bin].valid ? '1' : '0';
      validity[AXISYMMETRIC_MATCHED_FILM_BINS] = '\0';
      fprintf(stderr,
        "[TRACK-A-REJECT] t=%.17g i=%d dt=%.17g status=%d native=%d "
        "observation=%d model=%d valid_bins=%d bins=%d valid_area=%.17g "
        "nonpositive=%d invalid_plic=%d degenerate_plic=%d reason=%s "
        "validity=%s\n",
        t, i, dt, status,
        track_a_musehane_state.latest_native_status,
        track_a_musehane_state.latest_observation_status,
        track_a_musehane_state.latest_result.status,
        failed_geometry.valid_bins, failed_geometry.bins,
        failed_geometry.valid_area_fraction,
        failed_geometry.nonpositive_gap_bins,
        failed_geometry.invalid_plic_segments,
        failed_geometry.skipped_degenerate_plic_facets,
        track_a_musehane_state.latest_result.reason ?
          track_a_musehane_state.latest_result.reason : "unavailable",
        validity);
    }
    exit(6);
  }
  track_a_musehane_append_diagnostic(t, i, dt);
}

event track_a_musehane_finish (t = END_TIME)
{
  if (pid() == 0)
    fprintf(stderr, "# TRACK_A_FINISH t=%.17g i=%d accepted_steps=%lu\n",
      t, i, track_a_musehane_state.accepted_steps);
  track_a_musehane_diagnostic_destroy(&track_a_musehane_state);
}

#endif
