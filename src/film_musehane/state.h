#ifndef FILM_MUSEHANE_STATE_H
#define FILM_MUSEHANE_STATE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "film_musehane/core.h"

#define MUSEHANE_STATE_SNAPSHOT_SCHEMA 2U

typedef enum {
  MUSEHANE_STATE_OK = 0,
  MUSEHANE_STATE_INVALID_CONFIG = 1,
  MUSEHANE_STATE_INVALID_INPUT = 2,
  MUSEHANE_STATE_ALLOCATION_FAILED = 3,
  MUSEHANE_STATE_CORE_REJECTED = 4,
  MUSEHANE_STATE_RESTART_MISMATCH = 5,
  MUSEHANE_STATE_COUNTER_OVERFLOW = 6
} MusehaneStateStatus;

typedef struct MusehaneState MusehaneState;

typedef struct {
  const MusehaneSurfaceMesh *mesh;
  const unsigned char *active;
  const double *estimated_thickness;
  /* Optional per-cell handoff threshold. An active analytical-film cell is
     released after an accepted substep once its evolved thickness reaches
     this resolved-separation scale. */
  const double *resolved_separation_threshold;
  const double *mean_surface_velocity_divergence;
  const double *decontaminated_pressure_gradient_cell;
  const double *decontaminated_pressure_gradient_face;
  double dt;
} MusehaneStateAdvanceInput;

typedef struct {
  bool valid;
  MusehaneStatus core_status;
  MusehaneCoreFailure core_failure;
  double time_before;
  double time_after;
  unsigned long update_count;
  unsigned substep_count;
  unsigned resolved_release_count;
  double minimum_thickness;
  double inventory_before;
  double inventory_after;
  double inventory_closure_error;
  double pressure_inventory_change;
  double kinematic_inventory_change;
  double resolved_release_inventory_change;
  double poiseuille_dissipation_rate;
  double cumulative_poiseuille_dissipation;
} MusehaneStateAdvanceResult;

/* Caller supplies arrays and capacity for export. The schema and mesh
   signature bind a restart to the same radial surface discretisation. */
typedef struct {
  unsigned schema_version;
  size_t cell_count;
  size_t capacity;
  uint64_t mesh_signature;
  double time;
  unsigned long update_count;
  double cumulative_poiseuille_dissipation;
  double *estimated_thickness;
  double *thickness;
  unsigned char *active;
} MusehaneStateSnapshot;

MusehaneStateStatus musehane_state_create (
  const MusehaneSurfaceMesh *mesh,
  const double *initial_thickness,
  double initial_time,
  MusehaneState **output);

void musehane_state_destroy (MusehaneState *state);

MusehaneStateStatus musehane_state_advance (
  const MusehaneConfig *config,
  MusehaneState *state,
  const MusehaneStateAdvanceInput *input,
  MusehaneStateAdvanceResult *output);

MusehaneStateStatus musehane_state_snapshot (
  const MusehaneState *state,
  MusehaneStateSnapshot *output);

MusehaneStateStatus musehane_state_restore (
  const MusehaneSurfaceMesh *mesh,
  const MusehaneStateSnapshot *snapshot,
  MusehaneState **output);

#endif
