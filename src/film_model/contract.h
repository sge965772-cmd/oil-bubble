#ifndef FILM_MODEL_CONTRACT_H
#define FILM_MODEL_CONTRACT_H

#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

#define FILM_MODEL_CONTRACT_VERSION 1U
#define FILM_MODEL_RESTART_SCHEMA_VERSION 1U

typedef enum {
  FILM_MODEL_NONE = 0,
  FILM_MODEL_MUSEHANE_INVERSE = 1,
  FILM_MODEL_REYNOLDS_V2 = 2
} FilmModelKind;

typedef enum {
  FILM_MODEL_OK = 0,
  FILM_MODEL_INVALID_CONFIG = 1,
  FILM_MODEL_INVALID_DNS_VIEW = 2,
  FILM_MODEL_INVALID_STATE = 3,
  FILM_MODEL_STEP_REJECTED = 4,
  FILM_MODEL_TRACTION_REJECTED = 5,
  FILM_MODEL_RESTART_MISMATCH = 6
} FilmModelStatus;

typedef enum {
  FILM_STAGE_OBSERVE_N = 0,
  FILM_STAGE_UPDATE_N_TO_NP1 = 1,
  FILM_STAGE_EVALUATE_TRACTION_N = 2,
  FILM_STAGE_APPLY_TRACTION_N = 3,
  FILM_STAGE_PROJECT_DNS_N_TO_NP1 = 4
} FilmStepStage;

/* Axisymmetric convention: x is upward along the collision axis, r >= 0 is
   radial, and the film normal points from the lower to the upper envelope.
   A positive closing_rate is lower_normal_velocity-upper_normal_velocity. */
typedef struct {
  double radius;
  double width;
  /* Physical axisymmetric annular surface area, including the 2*pi factor. */
  double area_weight;
  double lower_position;
  double upper_position;
  double gap;
  double midpoint;
  double local_delta;
  double lower_normal_velocity;
  double upper_normal_velocity;
  double lower_tangential_velocity;
  double upper_tangential_velocity;
  double water_pressure;
  double pressure_gradient_tangent;
  double capillary_force_density_tangent;
  double lower_normal_x;
  double lower_normal_r;
  double upper_normal_x;
  double upper_normal_r;
  bool geometry_valid;
  /* Native PLIC validity is distinct from a usable subgrid-backed geometry.
     The release latch is caller-owned restart state and prevents analytical
     film reactivation until native geometry has recovered. */
  bool native_geometry_valid;
  bool resolved_release_latched;
  bool lower_identity_valid;
  bool upper_identity_valid;
} FilmRingObservation;

typedef struct {
  double time;
  double dt;
  unsigned long iteration;
  size_t ring_count;
  const FilmRingObservation *rings;
  bool anchor_released;
  uint32_t lower_identity;
  uint32_t upper_identity;
} FilmDNSView;

typedef struct {
  double lower_normal_force;
  double upper_normal_force;
  double lower_tangential_force;
  double upper_tangential_force;
  /* Reaction carried by the resolved film control volume. This is zero for
     an internally balanced pair and closes the pressure-driven Eq. (27)
     matched correction when that correction is selected. */
  double resolved_tangential_reaction_force;
} FilmRingTraction;

typedef struct {
  FilmModelStatus status;
  FilmModelKind model_kind;
  bool active;
  bool commit_allowed;
  size_t ring_count;
  const FilmRingTraction *tractions;
  double minimum_thickness;
  double inventory;
  double paired_force_residual;
  double paired_moment_residual;
  double step_work;
  double cumulative_work;
  FilmStepStage completed_stage;
  const char *reason;
  bool failure_diagnostic_valid;
  int failure_kind;
  size_t failure_cell;
  unsigned failure_iteration;
  unsigned failure_line_search_backtracks;
  double failure_initial_merit;
  double failure_final_merit;
  double failure_numerator;
  double failure_denominator;
  double failure_old_thickness;
  double failure_current_thickness;
  double failure_flux_in;
  double failure_flux_out;
  double failure_surface_velocity_divergence;
  double failure_equation_relative_residual;
  bool failure_inventory_diagnostic_valid;
  double failure_inventory_before;
  double failure_inventory_after;
  double failure_expected_inventory_after;
  double failure_inventory_closure_error;
  double failure_inventory_relative_error;
  double failure_inventory_relative_tolerance;
  bool failure_traction_diagnostic_valid;
  int failure_traction_status;
  size_t failure_traction_cell;
  double failure_traction_pressure_power;
  double failure_traction_shear_power;
  double failure_traction_direct_power;
  double failure_traction_total_power;
  double failure_traction_power_identity_residual;
  double failure_traction_step_work;
  double failure_traction_cumulative_work;
} FilmModelResult;

typedef struct {
  size_t ring_count;
  double *lower_normal_force;
  double *upper_normal_force;
  double *lower_tangential_force;
  double *upper_tangential_force;
} FilmDNSTractionTarget;

typedef struct {
  uint32_t schema_version;
  FilmModelKind model_kind;
  size_t size;
  size_t capacity;
  unsigned char *data;
} FilmRestartRecord;

static inline bool film_dns_view_is_well_formed (const FilmDNSView *dns)
{
  if (!dns || !isfinite(dns->time) || dns->time < 0. ||
      !isfinite(dns->dt) || !(dns->dt > 0.) ||
      dns->ring_count == 0 || !dns->rings ||
      dns->lower_identity == 0U || dns->upper_identity == 0U ||
      dns->lower_identity == dns->upper_identity)
    return false;
  for (size_t i = 0; i < dns->ring_count; i++) {
    const FilmRingObservation *ring = &dns->rings[i];
    const double position_gap = ring->upper_position - ring->lower_position;
    const double position_scale = fmax(DBL_MIN, fmax(fabs(ring->gap),
      fmax(fabs(ring->lower_position), fabs(ring->upper_position))));
    const double midpoint = .5*(ring->lower_position + ring->upper_position);
    const double midpoint_scale = fmax(DBL_MIN, fmax(fabs(midpoint),
      fabs(ring->midpoint)));
    const double lower_normal_norm = hypot(ring->lower_normal_x,
      ring->lower_normal_r);
    const double upper_normal_norm = hypot(ring->upper_normal_x,
      ring->upper_normal_r);
    if (!isfinite(ring->radius) || !(ring->radius >= 0.) ||
        !isfinite(ring->width) || !(ring->width > 0.) ||
        !isfinite(ring->area_weight) || !(ring->area_weight > 0.) ||
        !isfinite(ring->lower_position) ||
        !isfinite(ring->upper_position) ||
        !isfinite(ring->gap) || !(ring->gap > 0.) ||
        fabs(position_gap - ring->gap) > 256.*DBL_EPSILON*position_scale ||
        !isfinite(ring->midpoint) ||
        fabs(midpoint - ring->midpoint) > 256.*DBL_EPSILON*midpoint_scale ||
        !isfinite(ring->local_delta) || !(ring->local_delta > 0.) ||
        !isfinite(ring->lower_normal_velocity) ||
        !isfinite(ring->upper_normal_velocity) ||
        !isfinite(ring->lower_tangential_velocity) ||
        !isfinite(ring->upper_tangential_velocity) ||
        !isfinite(ring->water_pressure) ||
        !isfinite(ring->pressure_gradient_tangent) ||
        !isfinite(ring->capillary_force_density_tangent) ||
        !isfinite(lower_normal_norm) || !(lower_normal_norm > DBL_MIN) ||
        !isfinite(upper_normal_norm) || !(upper_normal_norm > DBL_MIN))
      return false;
  }
  return true;
}

static inline bool film_dns_view_has_usable_geometry (const FilmDNSView *dns)
{
  if (!film_dns_view_is_well_formed(dns))
    return false;
  for (size_t i = 0; i < dns->ring_count; i++)
    if (!dns->rings[i].geometry_valid ||
        !dns->rings[i].lower_identity_valid ||
        !dns->rings[i].upper_identity_valid)
      return false;
  return true;
}

#endif
