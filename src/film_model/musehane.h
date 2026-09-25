#ifndef FILM_MODEL_MUSEHANE_H
#define FILM_MODEL_MUSEHANE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "film_model/contract.h"
#include "film_musehane/core.h"
#include "film_musehane/pressure_driver.h"

#define MUSEHANE_FILM_MODEL_RESTART_SCHEMA 3U
#define FILM_MODEL_STATE_INITIALIZER {0U, NULL}

typedef enum {
  MUSEHANE_TRACTION_DIAGNOSTIC_ONLY = 0,
  MUSEHANE_MATCHED_VISCOUS_CORRECTION = 1,
  MUSEHANE_MATCHED_EQUATION27_CORRECTION = 2,
  MUSEHANE_PUBLISHED_EQUATION27_TRACTION = 3
} MusehaneTractionPolicy;

typedef struct FilmModelConfig {
  FilmModelKind model_kind;
  uint32_t contract_version;
  MusehaneConfig equation;
  MusehanePressureDriverConfig pressure_driver;
  double activation_cells;
  MusehaneTractionPolicy traction_policy;
} FilmModelConfig;

/* The implementation pointer is opaque to callers. A state must start with
   FILM_MODEL_STATE_INITIALIZER and must be destroyed before reuse. */
typedef struct FilmModelState {
  uint64_t guard;
  void *implementation;
} FilmModelState;

typedef struct {
  bool valid;
  bool one_step_lag;
  bool viscous_correction_applied;
  bool equation27_correction_applied;
  bool published_equation27_applied;
  size_t ring_count;
  unsigned long completed_steps;
  unsigned substep_count;
  unsigned resolved_release_count;
  double time_n;
  double time_np1;
  const unsigned char *active_n;
  const unsigned char *active_np1;
  const unsigned char *traction_active;
  const double *estimated_thickness;
  const double *thickness_n;
  const double *thickness_np1;
  const double *mean_surface_velocity_divergence;
  const double *decontaminated_pressure_gradient;
  const double *smoothed_pressure_gradient;
  const double *lower_analytic_tangential_traction;
  const double *upper_analytic_tangential_traction;
  const double *lower_matched_correction_traction;
  const double *upper_matched_correction_traction;
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
  double maximum_pressure_driver_residual;
  double maximum_traction_constitutive_residual;
  double pressure_inventory_change;
  double kinematic_inventory_change;
  double resolved_release_inventory_change;
  double poiseuille_dissipation_rate;
  double cumulative_poiseuille_dissipation;
} MusehaneFilmDiagnostics;

/* Minimal accepted state exposed only between checkpoint restore and the
   first successful post-restart advance.  This is not a reconstructed
   diagnostics record: it contains only the state serialized in the
   checkpoint and is therefore safe to use for restart-time geometry
   coupling.  Returned arrays are state-owned. */
typedef struct {
  bool valid;
  size_t ring_count;
  unsigned long completed_steps;
  double time;
  const unsigned char *active;
  const double *estimated_thickness;
  const double *thickness;
} MusehaneFilmRestoredStateView;

/* Exact input to the most recent rejected Musehane core solve. The arrays are
   model-owned and remain valid only until the next advance or destruction.
   This diagnostics-only view lets callers preserve a deterministic pure-C
   failure fixture without reconstructing inputs from post-rollback state. */
typedef struct {
  bool valid;
  MusehaneConfig config;
  MusehaneSurfaceMesh mesh;
  double time;
  double dt;
  const unsigned char *active;
  const double *estimated_thickness;
  const double *old_thickness;
  const double *mean_surface_velocity_divergence;
  const double *decontaminated_pressure_gradient_cell;
  const double *decontaminated_pressure_gradient_face;
} MusehaneRejectedCoreInputView;

FilmModelConfig film_model_musehane_config (void);

FilmModelStatus film_model_initialize (
  FilmModelState *state, const FilmModelConfig *config,
  const FilmDNSView *dns);

void film_model_destroy (FilmModelState *state);

FilmModelResult film_model_advance (
  FilmModelState *state, const FilmModelConfig *config,
  const FilmDNSView *dns, double dt);

FilmModelStatus film_model_apply_traction (
  const FilmModelResult *result, FilmDNSTractionTarget *target);

size_t film_model_checkpoint_size (const FilmModelState *state);

FilmModelStatus film_model_checkpoint (
  const FilmModelState *state, FilmRestartRecord *record);

FilmModelStatus film_model_restore (
  FilmModelState *state, const FilmModelConfig *config,
  const FilmRestartRecord *record);

/* Basilisk restarts skip events at the dump clock and enter through a small
   startup step. This one-shot operation aligns only the restored model clock;
   it does not evolve thickness, inventory, dissipation or correction work. */
FilmModelStatus film_model_rebase_restart_clock (
  FilmModelState *state, double new_time, double maximum_gap);

FilmModelStatus musehane_film_model_restored_state_view (
  const FilmModelState *state, MusehaneFilmRestoredStateView *output);

FilmModelStatus musehane_film_model_rejected_core_input_view (
  const FilmModelState *state, MusehaneRejectedCoreInputView *output);

/* Returned arrays are state-owned and remain valid only until the next
   successful advance or destruction. Analytical tractions remain diagnostic.
   Matched correction arrays contain the applied stress replacement for either
   matched policy. The Eq. (27) policy also reports the resolved-film reaction
   needed to close its pressure-driven interface-pair force. */
FilmModelStatus musehane_film_model_diagnostics (
  const FilmModelState *state, MusehaneFilmDiagnostics *output);

#endif
