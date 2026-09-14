#ifdef COLLISION_CASE_CONFIG_HEADER
# define COLLISION_CONFIG_STRINGIFY_INNER(value) #value
# define COLLISION_CONFIG_STRINGIFY(value) \
    COLLISION_CONFIG_STRINGIFY_INNER(value)
# include COLLISION_CONFIG_STRINGIFY(COLLISION_CASE_CONFIG_HEADER)
# undef COLLISION_CONFIG_STRINGIFY
# undef COLLISION_CONFIG_STRINGIFY_INNER
#endif

#include "grid/quadtree.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "axisymmetric_numerical_policy.h"
#include "axisymmetric_compound_profile.h"
#include "dual_compound_phases.h"
#ifndef MOMENTUM_TRANSPORT_VERSION
# define MOMENTUM_TRANSPORT_VERSION 3
#endif
#include "positive_material_momentum_policy.h"
#include "tension.h"
#include "separated_solver_tolerance_policy.h"
#include "collision_state.h"
#include "axisymmetric_contact_geometry.h"
#include "contact_resolution_metric.h"
#include "contact_amr_policy.h"
#include "axisymmetric_film_distribution.h"
#include "axisymmetric_film_state_axi.h"
#include "axisymmetric_film_runtime.h"
#include "axisymmetric_collision_observables.h"
#include "axisymmetric_velocity_audit.h"
#include "capillary_anchor_controller.h"
#include "capillary_anchor_axi.h"
#include <float.h>
#include <string.h>

#ifndef MINLEVEL
# define MINLEVEL 7
#endif
#ifndef MAXLEVEL
# define MAXLEVEL 10
#endif
#ifndef DOMAIN_LENGTH
# define DOMAIN_LENGTH 0.080
#endif
#ifndef AXIAL_ORIGIN
# define AXIAL_ORIGIN -0.030
#endif
#ifndef OUTER_RADIUS
# define OUTER_RADIUS 0.002
#endif
#ifndef OIL_VOLUME_FRACTION
# define OIL_VOLUME_FRACTION 0.20
#endif
#ifndef INITIAL_GAP
# define INITIAL_GAP 8.e-4
#endif
#ifndef SIGMA_OW
# define SIGMA_OW 40.9e-3
#endif
#ifndef SIGMA_OG
# define SIGMA_OG 19.1e-3
#endif
#ifndef GRAVITY_MAGNITUDE
# define GRAVITY_MAGNITUDE 9.81
#endif
#ifndef END_TIME
# define END_TIME 1.e-3
#endif
#ifndef DT_MAX
# define DT_MAX 5.e-7
#endif
#ifndef OUTPUT_INTERVAL
# define OUTPUT_INTERVAL 1.e-5
#endif
#ifndef DUMP_INTERVAL
# define DUMP_INTERVAL 1.e-4
#endif
#ifndef FRACTION_TOLERANCE
# define FRACTION_TOLERANCE 1.e-6
#endif
#ifndef VELOCITY_TOLERANCE
# define VELOCITY_TOLERANCE 1.e-4
#endif
#ifndef ENABLE_CAPILLARY_ANCHOR
# define ENABLE_CAPILLARY_ANCHOR 1
#endif
#ifndef ANCHOR_RELEASE_GAP
# define ANCHOR_RELEASE_GAP 2.e-4
#endif
#ifndef ANCHOR_REQUIRE_APPROACH
# define ANCHOR_REQUIRE_APPROACH 1
#endif
#ifndef ANCHOR_PROPORTIONAL_GAIN
# define ANCHOR_PROPORTIONAL_GAIN 4000.
#endif
#ifndef ANCHOR_INTEGRAL_GAIN
# define ANCHOR_INTEGRAL_GAIN 4.e6
#endif
#ifndef ANCHOR_ACCELERATION_LIMIT
# define ANCHOR_ACCELERATION_LIMIT 100.
#endif
#ifndef CONTACT_SUPPORT_HALF_WIDTH
# define CONTACT_SUPPORT_HALF_WIDTH 5.e-4
#endif
#ifndef UPPER_PROFILE_MODE
# define UPPER_PROFILE_MODE 0
#endif
#ifndef UPPER_PROFILE_DIAMETER_RELATIVE_TOLERANCE
# define UPPER_PROFILE_DIAMETER_RELATIVE_TOLERANCE 0.03
#endif
#ifndef UPPER_PROFILE_OIL_FRACTION_TOLERANCE
# define UPPER_PROFILE_OIL_FRACTION_TOLERANCE 0.01
#endif
#ifndef INITIAL_GAP_TOLERANCE_CELLS
# define INITIAL_GAP_TOLERANCE_CELLS 0.5
#endif
#ifndef LOG_GAP_ALIGNMENT_SEARCH
# define LOG_GAP_ALIGNMENT_SEARCH 0
#endif
#ifndef ENABLE_STATEFUL_FILM
# define ENABLE_STATEFUL_FILM 0
#endif
#ifndef MATCHED_FILM_WATER_VISCOSITY
# define MATCHED_FILM_WATER_VISCOSITY 0.89e-3
#endif
#ifndef MATCHED_FILM_MATCH_CELLS
# define MATCHED_FILM_MATCH_CELLS 4.
#endif
#ifndef MATCHED_FILM_RELEASE_FACTOR
# define MATCHED_FILM_RELEASE_FACTOR 1.125
#endif
#ifndef FILM_STATE_FLUX_BALANCE_TOLERANCE
# define FILM_STATE_FLUX_BALANCE_TOLERANCE 1.e-12
#endif
#ifndef CONTACT_AMR_ACTIVATION_CELLS
# define CONTACT_AMR_ACTIVATION_CELLS 8.
#endif
#ifndef CONTACT_AMR_AXIAL_PADDING_CELLS
# define CONTACT_AMR_AXIAL_PADDING_CELLS 2.
#endif
#ifndef CONTACT_AMR_RADIAL_PADDING_CELLS
# define CONTACT_AMR_RADIAL_PADDING_CELLS 1.
#endif
#ifndef COLLISION_FULL_STATE_RESTART
# define COLLISION_FULL_STATE_RESTART 0
#endif
#ifndef RESTART_FIELD_AUDIT
# define RESTART_FIELD_AUDIT 0
#endif
#ifndef RESTART_FORCE_CONTACT_AMR
# define RESTART_FORCE_CONTACT_AMR 0
#endif
#if COLLISION_FULL_STATE_RESTART
# ifndef RESTART_FILM_STATE_INCLUDE
#  error "full-state restart requires RESTART_FILM_STATE_INCLUDE"
# endif
# include RESTART_FILM_STATE_INCLUDE
# ifndef COLLISION_RESTART_FILM_STATE_V4_INCLUDED
#  error "full-state restart requires collision-restart-state-v4 include"
# endif
# ifndef RESTART_INITIAL_LOWER_GAS_VOLUME
#  error "full-state restart requires RESTART_INITIAL_LOWER_GAS_VOLUME"
# endif
# ifndef RESTART_INITIAL_LOWER_OIL_VOLUME
#  error "full-state restart requires RESTART_INITIAL_LOWER_OIL_VOLUME"
# endif
# ifndef RESTART_INITIAL_UPPER_GAS_VOLUME
#  error "full-state restart requires RESTART_INITIAL_UPPER_GAS_VOLUME"
# endif
# ifndef RESTART_INITIAL_UPPER_OIL_VOLUME
#  error "full-state restart requires RESTART_INITIAL_UPPER_OIL_VOLUME"
# endif
# ifndef RESTART_RELEASE_TIME
#  error "full-state restart requires RESTART_RELEASE_TIME"
# endif
# ifndef RESTART_RELEASE_GAP
#  error "full-state restart requires RESTART_RELEASE_GAP"
# endif
# ifndef RESTART_CONSTRAINT_WORK
#  error "full-state restart requires RESTART_CONSTRAINT_WORK"
# endif
# ifndef RESTART_SOURCE_TIME
#  error "full-state restart requires RESTART_SOURCE_TIME"
# endif
# ifndef RESTART_EXPECTED_FILM_STATE_SHA256
#  error "full-state restart requires RESTART_EXPECTED_FILM_STATE_SHA256"
# endif
# ifndef RESTART_EXPECTED_CONTRACT_SHA256
#  error "full-state restart requires RESTART_EXPECTED_CONTRACT_SHA256"
# endif
#endif

static const char * restore_path;
#if UPPER_PROFILE_MODE == 1
static const char * upper_profile_path;
#endif
static double inner_radius;
#if !COLLISION_FULL_STATE_RESTART
static double inserted_upper_center;
#endif
#if UPPER_PROFILE_MODE == 1
static double upper_profile_offset;
#endif
static double upper_outer_curvature_radius;
static double upper_gas_curvature_radius;
static double initial_gas_volumes[2];
static double initial_oil_volumes[2];
static CapillaryAnchorConfig anchor_config;
static CapillaryAnchorState anchor_state;
static CapillaryAnchorCommand latest_anchor_command;
static CapillaryAnchorLoad latest_anchor_load;
static AxisymmetricContactGeometry latest_contact_geometry;
static AxisymmetricFilmSamples latest_film_samples;
static AxisymmetricFilmStateConfig film_state_config;
#if ENABLE_STATEFUL_FILM
static AxisymmetricFilmState film_state;
static AxisymmetricFilmStateResult latest_film_state_result;
static AxisymmetricFilmStateAxiObservation latest_film_observation;
static AxisymmetricFilmStateAxiPairLedger latest_film_pair_ledger;
static AxisymmetricMatchedFilmResult latest_film_traction;
static face vector film_owner_mask[];
static double latest_film_timestep_limit = HUGE_VAL;
static bool latest_contact_corridor_at_maxlevel;
static bool positive_film_work;
static bool film_restart_restored;
static double film_restart_continuity_error;
static double film_handoff_time = NAN;
static double film_handoff_raw_valid_area = NAN;
static int film_handoff_raw_crossed_bins = -1;
static bool film_handoff_corridor_at_maxlevel;
static double last_film_primitive_time = NAN;
#endif
static bool film_state_failed;
static long latest_contact_amr_forced_cells;
static double estimated_buoyancy_force;
#if COLLISION_FULL_STATE_RESTART
static bool restart_pressure_boundary_initialized;
static bool restart_projection_input_failed;
static double restart_dump_time = NAN;
static int restart_dump_iteration = -1;
#if RESTART_FIELD_AUDIT
static int restart_projection_audit_count;
#endif
#endif

static inline bool collision_terminal_event_time (double current_time,
                                                   double end_time)
{
  const double tolerance = 128.*DBL_EPSILON*max(1., fabs(end_time));
  return current_time >= end_time - tolerance;
}

static void append_collision_diagnostics (double output_time,
                                          int output_iteration);
#if ENABLE_STATEFUL_FILM
static void append_film_state_primitive (double output_time,
                                         int output_iteration);
#endif

#if COLLISION_FULL_STATE_RESTART
static bool collision_read_restart_dump_clock (const char * path,
                                                double * dump_time,
                                                int * dump_iteration)
{
  if (!path || !dump_time || !dump_iteration)
    return false;
  FILE * fp = fopen(path, "rb");
  if (!fp)
    return false;
  struct DumpHeader header = {0};
  bool valid = fread(&header, sizeof(header), 1, fp) == 1;
  fclose(fp);
  valid = valid && isfinite(header.t) && header.t >= 0. &&
    header.i >= 0 && header.depth >= 0 && header.npe > 0;
  if (!valid)
    return false;
  *dump_time = header.t;
  *dump_iteration = header.i;
  return true;
}
#endif

static ContactAmrPolicy collision_contact_amr_policy (void)
{
  return (ContactAmrPolicy){
    .finest_delta = (double)DOMAIN_LENGTH/(1 << MAXLEVEL),
    .activation_cells = (double)CONTACT_AMR_ACTIVATION_CELLS,
    .axial_padding_cells = (double)CONTACT_AMR_AXIAL_PADDING_CELLS,
    .radial_support = (double)CONTACT_SUPPORT_HALF_WIDTH,
    .radial_padding_cells = (double)CONTACT_AMR_RADIAL_PADDING_CELLS
  };
}

static ContactAmrGeometry collision_raw_contact_geometry (
  AxisymmetricContactGeometry contact)
{
  return (ContactAmrGeometry){
    .valid = contact.valid,
    .lower_top = contact.lower_top,
    .upper_bottom = contact.upper_bottom
  };
}

static bool collision_contact_corridor_at_maxlevel (
  const ContactAmrGeometry * geometry)
{
  ContactAmrPolicy policy = collision_contact_amr_policy();
  long corridor_cells = 0, coarse_cells = 0;
  foreach(reduction(+:corridor_cells) reduction(+:coarse_cells))
    if (contact_amr_cell_intersects(&policy, geometry, x, y, Delta)) {
      corridor_cells++;
      if (level < MAXLEVEL)
        coarse_cells++;
    }
  return corridor_cells > 0 && coarse_cells == 0;
}

static long collision_refine_contact_corridor (
  const ContactAmrGeometry * geometry)
{
#if TREE
  ContactAmrPolicy policy = collision_contact_amr_policy();
  long forced_cells = 0;
  foreach(reduction(+:forced_cells))
    if (level < MAXLEVEL &&
        contact_amr_cell_intersects(&policy, geometry, x, y, Delta))
      forced_cells++;
  for (int pass = 0; pass < MAXLEVEL - MINLEVEL; pass++)
    refine(level < MAXLEVEL &&
           contact_amr_cell_intersects(&policy, geometry, x, y, Delta));
  return forced_cells;
#else
  (void) geometry;
  return 0;
#endif
}

#if COLLISION_FULL_STATE_RESTART && ENABLE_STATEFUL_FILM
static bool axisymmetric_film_state_restore (AxisymmetricFilmState * state)
{
  const double restart_length_unit = 1. [1];
  const double restart_time_unit = 1. [0,1];
  const double restart_volume_unit = 1. [3];
  const double restart_volume_flux_unit = 1. [3,-1];
  const double restart_work_unit = 1. [2,-2,1];
  if (!state || strcmp(COLLISION_RESTART_FILM_STATE_SCHEMA,
        "collision-restart-state-v4") ||
      strcmp(COLLISION_RESTART_FILM_MODEL_IDENTIFIER,
        "axisymmetric-stateful-reynolds-film-v1") ||
      strcmp(COLLISION_RESTART_FILM_STATE_SHA256,
        RESTART_EXPECTED_FILM_STATE_SHA256) ||
      strcmp(COLLISION_RESTART_CONTRACT_SHA256,
        RESTART_EXPECTED_CONTRACT_SHA256) ||
      fabs((double)COLLISION_RESTART_FILM_STATE_TIME -
        (double)RESTART_SOURCE_TIME) >
        128.*DBL_EPSILON*max(1., fabs((double)RESTART_SOURCE_TIME)))
    return false;
  AxisymmetricFilmState restored = {0};
  restored.active = (bool)COLLISION_RESTART_FILM_STATE_ACTIVE;
  restored.activation_count =
    COLLISION_RESTART_FILM_STATE_ACTIVATION_COUNT;
  restored.update_count = COLLISION_RESTART_FILM_STATE_UPDATE_COUNT;
  restored.reacquisition_consecutive_observations =
    COLLISION_RESTART_FILM_STATE_REACQUISITION_CONSECUTIVE_OBSERVATIONS;
  restored.cumulative_unresolved_observations =
    COLLISION_RESTART_FILM_STATE_CUMULATIVE_UNRESOLVED_OBSERVATIONS;
  restored.cumulative_raw_crossing_observations =
    COLLISION_RESTART_FILM_STATE_CUMULATIVE_RAW_CROSSING_OBSERVATIONS;
  restored.source_time =
    COLLISION_RESTART_FILM_STATE_SOURCE_TIME*restart_time_unit;
  restored.handoff_inventory =
    COLLISION_RESTART_FILM_STATE_HANDOFF_INVENTORY*restart_volume_unit;
  restored.inventory =
    COLLISION_RESTART_FILM_STATE_INVENTORY*restart_volume_unit;
  restored.cumulative_edge_outflow =
    COLLISION_RESTART_FILM_STATE_CUMULATIVE_EDGE_OUTFLOW*restart_volume_unit;
  restored.cumulative_work =
    COLLISION_RESTART_FILM_STATE_CUMULATIVE_WORK*restart_work_unit;
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++) {
    restored.thickness[bin] =
      collision_restart_film_thickness[bin]*restart_length_unit;
    restored.midpoint[bin] =
      collision_restart_film_midpoint[bin]*restart_length_unit;
    restored.radius[bin] =
      collision_restart_film_radius[bin]*restart_length_unit;
    restored.width[bin] =
      collision_restart_film_width[bin]*restart_length_unit;
    restored.delta[bin] =
      collision_restart_film_delta[bin]*restart_length_unit;
    restored.face_flux[bin] =
      collision_restart_film_face_flux[bin]*restart_volume_flux_unit;
  }
  restored.face_flux[AXISYMMETRIC_FILM_STATE_BINS] =
    collision_restart_film_face_flux[AXISYMMETRIC_FILM_STATE_BINS]*
    restart_volume_flux_unit;
  if (restored.active &&
      (!axisymmetric_film_state_active_state_is_valid(&restored) ||
       !axisymmetric_film_state_active_ledger_is_valid(&restored)))
    return false;
  if (!restored.active && restored.activation_count == 0UL) {
    if (restored.update_count || restored.handoff_inventory != 0. ||
        restored.inventory != 0. || restored.cumulative_work != 0.)
      return false;
  }
  *state = restored;
  return true;
}
#endif

#if ENABLE_STATEFUL_FILM
static bool collision_refresh_film_observation (void)
{
  latest_film_samples = (AxisymmetricFilmSamples){0};
  latest_film_state_result = (AxisymmetricFilmStateResult){
    .status = AXISYMMETRIC_FILM_STATE_INACTIVE,
    .minimum_state_thickness = NAN,
    .minimum_raw_gap = NAN,
    .maximum_state_raw_gap_discrepancy = NAN
  };
  latest_film_observation = (AxisymmetricFilmStateAxiObservation){0};
  latest_film_pair_ledger = (AxisymmetricFilmStateAxiPairLedger){0};
  latest_film_traction = (AxisymmetricMatchedFilmResult){0};
  positive_film_work = false;

  const double finest_delta = (double)DOMAIN_LENGTH/(1 << MAXLEVEL);
  const bool sample_required = film_state.active ||
    (latest_contact_geometry.valid && latest_contact_geometry.gap <=
      AXISYMMETRIC_FILM_STATE_HANDOFF_CELLS*finest_delta);
  if (!sample_required) {
    latest_film_state_result.valid = true;
    latest_film_timestep_limit = HUGE_VAL;
    return true;
  }

  latest_film_samples = sample_axisymmetric_water_film(
    (double)CONTACT_SUPPORT_HALF_WIDTH);
  ContactAmrGeometry observation_geometry = film_state.active ?
    axisymmetric_film_state_contact_geometry(&film_state) :
    collision_raw_contact_geometry(latest_contact_geometry);
  latest_contact_corridor_at_maxlevel =
    collision_contact_corridor_at_maxlevel(&observation_geometry);
  latest_film_observation = axisymmetric_film_state_observe_axi(
    &film_state, &latest_film_samples,
    latest_contact_corridor_at_maxlevel);
  if (!latest_film_observation.valid) {
    latest_film_state_result.status =
      AXISYMMETRIC_FILM_STATE_INVALID_OBSERVATION;
    return false;
  }
  return true;
}

static void collision_abort_film_transaction (
  const AxisymmetricFilmState * state_before_update)
{
  if (state_before_update)
    film_state = *state_before_update;
  film_state_failed = true;
}
#endif
#if UPPER_PROFILE_MODE == 1
static AxisymmetricCompoundProfile upper_profile;
#endif

#if !COLLISION_FULL_STATE_RESTART
static void set_upper_geometry (void)
{
#if UPPER_PROFILE_MODE == 1
  vertex scalar outer_levelset[], gas_levelset[];
  foreach_vertex() {
    outer_levelset[] = axisymmetric_profile_levelset(
      &upper_profile, x, y, upper_profile_offset,
      AXISYMMETRIC_PROFILE_OUTER);
    gas_levelset[] = axisymmetric_profile_levelset(
      &upper_profile, x, y, upper_profile_offset,
      AXISYMMETRIC_PROFILE_GAS);
  }
  fractions(outer_levelset, upper_envelope_fraction);
  fractions(gas_levelset, upper_gas_fraction);
  foreach()
    upper_gas_fraction[] = min(
      upper_gas_fraction[], upper_envelope_fraction[]);
  boundary((scalar *){upper_envelope_fraction, upper_gas_fraction});
#else
  fraction(upper_envelope_fraction,
           sq((double)OUTER_RADIUS) -
           sq(x - inserted_upper_center) - sq(y));
  fraction(upper_gas_fraction,
           sq(inner_radius) - sq(x - inserted_upper_center) - sq(y));
#endif
}

static void align_upper_geometry_to_initial_gap (bool enforce_tolerance)
{
  AxisymmetricContactGeometry initial =
    measure_axisymmetric_contact_geometry((double)CONTACT_SUPPORT_HALF_WIDTH);
  if (!initial.valid) {
    fprintf(stderr, "upper profile has no valid PLIC contact support\n");
    exit(2);
  }

  const int samples = 33;
  double finest_dx = (double)DOMAIN_LENGTH/(1 << MAXLEVEL);
  double tolerance = (double)INITIAL_GAP_TOLERANCE_CELLS*finest_dx;
  double original_center = inserted_upper_center;
#if UPPER_PROFILE_MODE == 1
  double original_offset = upper_profile_offset;
#endif
  double predicted_shift = (double)INITIAL_GAP - initial.gap;
  double best_shift = 0., best_error = HUGE_VAL, best_gap = NAN;
  for (int sample = 0; sample < samples; sample++) {
    double search_coordinate = -1. + 2.*sample/(samples - 1.);
    double shift = predicted_shift + search_coordinate*finest_dx;
    inserted_upper_center = original_center + shift;
#if UPPER_PROFILE_MODE == 1
    upper_profile_offset = original_offset + shift;
#endif
    set_upper_geometry();
    AxisymmetricContactGeometry candidate =
      measure_axisymmetric_contact_geometry((double)CONTACT_SUPPORT_HALF_WIDTH);
#if LOG_GAP_ALIGNMENT_SEARCH
    if (pid() == 0)
      fprintf(stderr, "# GAP_SEARCH sample=%d shift=%g gap=%g valid=%d\n",
              sample, shift, candidate.gap, candidate.valid);
#endif
    if (candidate.valid) {
      double error = fabs(candidate.gap - (double)INITIAL_GAP);
      if (error < best_error) {
        best_error = error;
        best_shift = shift;
        best_gap = candidate.gap;
      }
    }
  }
  inserted_upper_center = original_center + best_shift;
#if UPPER_PROFILE_MODE == 1
  upper_profile_offset = original_offset + best_shift;
#endif
  set_upper_geometry();
  if (pid() == 0)
    fprintf(stderr,
      "# GAP_ALIGNMENT requested=%g measured=%g error_cells=%g samples=%d\n",
      (double)INITIAL_GAP, best_gap, best_error/finest_dx, samples);
  if (!isfinite(best_gap) || (enforce_tolerance && best_error > tolerance)) {
    fprintf(stderr,
      "initial PLIC gap alignment failed: requested=%g measured=%g "
      "error_cells=%g tolerance_cells=%g\n",
      (double)INITIAL_GAP, best_gap, best_error/finest_dx,
      (double)INITIAL_GAP_TOLERANCE_CELLS);
    exit(2);
  }
}

static bool collision_prepare_fresh_contact_corridor (void)
{
#if TREE
  const int maximum_passes = MAXLEVEL - MINLEVEL + 2;
  for (int pass = 0; pass < maximum_passes; pass++) {
    AxisymmetricContactGeometry contact =
      measure_axisymmetric_contact_geometry(
        (double)CONTACT_SUPPORT_HALF_WIDTH);
    ContactAmrGeometry geometry = collision_raw_contact_geometry(contact);
    ContactAmrPolicy policy = collision_contact_amr_policy();
    if (!contact_amr_is_active(&policy, &geometry))
      return true;
    if (collision_contact_corridor_at_maxlevel(&geometry))
      return true;
    const long forced_cells = collision_refine_contact_corridor(&geometry);
    latest_contact_amr_forced_cells += forced_cells;
    if (forced_cells <= 0)
      return false;
    set_upper_geometry();
    align_upper_geometry_to_initial_gap(true);
  }
  AxisymmetricContactGeometry contact =
    measure_axisymmetric_contact_geometry(
      (double)CONTACT_SUPPORT_HALF_WIDTH);
  ContactAmrGeometry geometry = collision_raw_contact_geometry(contact);
  return collision_contact_corridor_at_maxlevel(&geometry);
#else
  return true;
#endif
}
#endif

int main (int argc, char ** argv)
{
#if UPPER_PROFILE_MODE == 1
  if (argc != 3) {
    fprintf(stderr,
      "usage: %s NORMALIZED_PRECURSOR_DUMP UPPER_PROFILE_CSV\n", argv[0]);
    return 2;
  }
  upper_profile_path = argv[2];
  char profile_error[256];
  if (!axisymmetric_profile_load(
        &upper_profile, upper_profile_path,
        profile_error, sizeof(profile_error))) {
    fprintf(stderr, "invalid upper profile %s: %s\n",
            upper_profile_path, profile_error);
    return 2;
  }
  double equivalent_diameter = cbrt(
    6.*upper_profile.outer_volume/AXISYMMETRIC_PROFILE_PI);
  double diameter_error = fabs(
    equivalent_diameter - 2.*(double)OUTER_RADIUS)/
    (2.*(double)OUTER_RADIUS);
  double oil_fraction_error = fabs(
    upper_profile.oil_fraction - (double)OIL_VOLUME_FRACTION);
  if (diameter_error > (double)UPPER_PROFILE_DIAMETER_RELATIVE_TOLERANCE ||
      oil_fraction_error > (double)UPPER_PROFILE_OIL_FRACTION_TOLERANCE) {
    fprintf(stderr,
      "upper profile mismatch: diameter_error=%g oil_fraction_error=%g\n",
      diameter_error, oil_fraction_error);
    return 2;
  }
#else
  if (argc != 2) {
    fprintf(stderr, "usage: %s NORMALIZED_PRECURSOR_DUMP\n", argv[0]);
    return 2;
  }
#endif
  restore_path = argv[1];
  inner_radius = (double)OUTER_RADIUS*
    pow(1. - (double)OIL_VOLUME_FRACTION, 1./3.);
#if UPPER_PROFILE_MODE == 1
  upper_outer_curvature_radius = cbrt(
    3.*upper_profile.outer_volume/(4.*AXISYMMETRIC_PROFILE_PI));
  upper_gas_curvature_radius = cbrt(
    3.*upper_profile.gas_volume/(4.*AXISYMMETRIC_PROFILE_PI));
#else
  upper_outer_curvature_radius = (double)OUTER_RADIUS;
  upper_gas_curvature_radius = inner_radius;
#endif
  anchor_config = capillary_anchor_config(
    (double)ANCHOR_RELEASE_GAP,
    (double)ANCHOR_PROPORTIONAL_GAIN,
    (double)ANCHOR_INTEGRAL_GAIN,
    (double)ANCHOR_ACCELERATION_LIMIT,
    0.,
    (bool)ANCHOR_REQUIRE_APPROACH);
  anchor_state = ENABLE_CAPILLARY_ANCHOR ?
    capillary_anchor_initial_state() : capillary_anchor_released_state();
  film_state_config = (AxisymmetricFilmStateConfig){
    .reynolds = {
      .water_viscosity = (double)MATCHED_FILM_WATER_VISCOSITY,
      .match_cells = (double)MATCHED_FILM_MATCH_CELLS,
      .release_factor = (double)MATCHED_FILM_RELEASE_FACTOR,
      .bins = AXISYMMETRIC_FILM_STATE_BINS
    },
    .minimum_raw_valid_area_fraction = 0.90,
    .flux_balance_tolerance = (double)FILM_STATE_FLUX_BALANCE_TOLERANCE
  };

  size((double)DOMAIN_LENGTH);
  origin((double)AXIAL_ORIGIN, 0.);
  init_grid(1 << MINLEVEL);
  DT = (double)DT_MAX;
  CFL = 0.20;
  TOLERANCE = (double)PRESSURE_TOLERANCE;
  NITERMAX = 10000;
  lower_gas_fraction.sigma = (double)SIGMA_OG;
  lower_envelope_fraction.sigma = (double)SIGMA_OW;
  upper_gas_fraction.sigma = (double)SIGMA_OG;
  upper_envelope_fraction.sigma = (double)SIGMA_OW;
  run();
#if COLLISION_FULL_STATE_RESTART
  if (restart_projection_input_failed)
    return 5;
#endif
#if ENABLE_STATEFUL_FILM
  if (film_state_failed)
    return 4;
#endif
  return dual_momentum_transport_failed ? 3 : 0;
}

u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);

event init (i = 0)
{
  p.nodump = false;
  pf.nodump = false;
#if COLLISION_FULL_STATE_RESTART
  if (!collision_read_restart_dump_clock(
        restore_path, &restart_dump_time, &restart_dump_iteration)) {
    fprintf(stderr, "cannot read restart dump clock: %s\n", restore_path);
    exit(5);
  }
#endif
  if (!restore(file = restore_path)) {
    fprintf(stderr, "restore failed: %s\n", restore_path);
    exit(2);
  }
#if COLLISION_FULL_STATE_RESTART
  boundary((scalar *){
    lower_gas_fraction, lower_envelope_fraction,
    upper_gas_fraction, upper_envelope_fraction, u, g
  });
  dual_compound_update_properties();
  boundary((scalar *){dual_compound_alpha, mu});
  // The normal projection refreshes pressure ghosts after face coefficients
  // have been restricted across every level of the restored tree.
  CollisionState initial = measure_collision_state();
  initial_gas_volumes[0] = (double)RESTART_INITIAL_LOWER_GAS_VOLUME;
  initial_oil_volumes[0] = (double)RESTART_INITIAL_LOWER_OIL_VOLUME;
  initial_gas_volumes[1] = (double)RESTART_INITIAL_UPPER_GAS_VOLUME;
  initial_oil_volumes[1] = (double)RESTART_INITIAL_UPPER_OIL_VOLUME;
  anchor_state = capillary_anchor_released_state();
  anchor_state.release_time = (double)RESTART_RELEASE_TIME;
  anchor_state.release_gap = (double)RESTART_RELEASE_GAP;
  anchor_state.cumulative_work = (double)RESTART_CONSTRAINT_WORK;
#if ENABLE_STATEFUL_FILM
  film_restart_restored = axisymmetric_film_state_restore(&film_state);
  film_restart_continuity_error = fabs(restart_dump_time -
    (double)COLLISION_RESTART_FILM_STATE_TIME);
  if (!film_restart_restored || film_restart_continuity_error >
      128.*DBL_EPSILON*max(1., fabs(restart_dump_time))) {
    fprintf(stderr,
      "stateful film restart mismatch: restored=%d dump_t=%.17g "
      "dump_i=%d init_event_t=%.17g init_event_i=%d state_t=%.17g "
      "error=%.17g\n",
      film_restart_restored, restart_dump_time, restart_dump_iteration, t, i,
      (double)COLLISION_RESTART_FILM_STATE_TIME,
      film_restart_continuity_error);
    film_state_failed = true;
    exit(5);
  }
#endif
  if (pid() == 0)
    fprintf(stderr,
      "# COLLISION_FULL_STATE_RESTART source=%s source_t=%.17g "
      "dump_t=%.17g dump_i=%d init_event_t=%.17g init_event_i=%d "
      "target_t=%.17g "
      "restart_lower_gas=%.17g restart_lower_oil=%.17g "
      "restart_upper_gas=%.17g restart_upper_oil=%.17g "
      "film_state_sha256=%s restart_contract_sha256=%s\n",
      restore_path, (double)RESTART_SOURCE_TIME,
      restart_dump_time, restart_dump_iteration, t, i, (double)END_TIME,
      initial.gas_volumes[0], initial.oil_volumes[0],
      initial.gas_volumes[1], initial.oil_volumes[1],
      COLLISION_RESTART_FILM_STATE_SHA256,
      COLLISION_RESTART_CONTRACT_SHA256);
#else
  dual_compound_update_properties();
  CollisionState before = measure_collision_state();
  inserted_upper_center = before.centers[0] +
                          2.*(double)OUTER_RADIUS + (double)INITIAL_GAP;
#if UPPER_PROFILE_MODE == 1
  upper_profile_offset = inserted_upper_center - upper_profile.outer_centroid;
#endif
  set_upper_geometry();
  align_upper_geometry_to_initial_gap(false);
#if TREE
  for (int pass = 0; pass < MAXLEVEL - MINLEVEL + 1; pass++) {
    adapt_wavelet((scalar *){
        lower_gas_fraction, lower_envelope_fraction,
        upper_gas_fraction, upper_envelope_fraction, u
      },
      (double[]){
        (double)FRACTION_TOLERANCE, (double)FRACTION_TOLERANCE,
        (double)FRACTION_TOLERANCE, (double)FRACTION_TOLERANCE,
        (double)VELOCITY_TOLERANCE, (double)VELOCITY_TOLERANCE
      }, MAXLEVEL, MINLEVEL);
    set_upper_geometry();
  }
#endif
  align_upper_geometry_to_initial_gap(true);
  if (!collision_prepare_fresh_contact_corridor()) {
    fprintf(stderr,
      "fresh collision contact corridor did not reach MAXLEVEL before "
      "the first film observation\n");
    exit(2);
  }
  dual_compound_update_properties();
  foreach() {
    double laplace_jump =
      2.*(double)SIGMA_OW/upper_outer_curvature_radius*
        upper_envelope_fraction[] +
      2.*(double)SIGMA_OG/upper_gas_curvature_radius*
        upper_gas_fraction[];
    p[] += laplace_jump;
    pf[] += laplace_jump;
  }
  boundary((scalar *){
    lower_gas_fraction, lower_envelope_fraction,
    upper_gas_fraction, upper_envelope_fraction, u, p, pf
  });
  CollisionState initial = measure_collision_state();
  for (int identity = 0; identity < 2; identity++) {
    initial_gas_volumes[identity] = initial.gas_volumes[identity];
    initial_oil_volumes[identity] = initial.oil_volumes[identity];
  }
#endif
  latest_contact_geometry = measure_axisymmetric_contact_geometry(
    (double)CONTACT_SUPPORT_HALF_WIDTH);
#if ENABLE_STATEFUL_FILM
  ContactAmrGeometry initial_film_amr_geometry = film_state.active ?
    axisymmetric_film_state_contact_geometry(&film_state) :
    collision_raw_contact_geometry(latest_contact_geometry);
  latest_contact_corridor_at_maxlevel =
    collision_contact_corridor_at_maxlevel(&initial_film_amr_geometry);
#if COLLISION_FULL_STATE_RESTART
  if (film_state.active) {
    if (!collision_refresh_film_observation()) {
      fprintf(stderr,
        "stateful film restart has no valid first-step observation at "
        "t=%.17g\n", t);
      film_state_failed = true;
      exit(5);
    }
    latest_film_timestep_limit =
      axisymmetric_film_runtime_timestep_limit(
      &film_state,
      &latest_film_observation.observation);
    if ((!isfinite(latest_film_timestep_limit) &&
         latest_film_timestep_limit != HUGE_VAL) ||
        (isfinite(latest_film_timestep_limit) &&
         latest_film_timestep_limit <= 0.)) {
      fprintf(stderr,
        "stateful film restart has invalid first-step timestep limit "
        "at t=%.17g limit=%.17g\n", t, latest_film_timestep_limit);
      film_state_failed = true;
      exit(5);
    }
  }
#endif
#endif
  latest_anchor_load = capillary_anchor_measure_axi();
  estimated_buoyancy_force =
    (rho_water*latest_anchor_load.support_volume -
     latest_anchor_load.effective_mass)*(double)GRAVITY_MAGNITUDE;
  anchor_config.feedforward_acceleration =
    -estimated_buoyancy_force/
    max(latest_anchor_load.effective_mass, 1.e-300);

  if (pid() == 0) {
    CollisionState initialized = initial;
    FILE * audit = fopen("data/initialization_audit.dat", "w");
    if (!audit) {
      perror("data/initialization_audit.dat");
      exit(2);
    }
    fprintf(audit,
      "# profile_mode requested_gap plic_gap plic_valid lower_center "
      "upper_center lower_u upper_u envelope_overlap gas_overlap "
      "nesting overfill\n"
      "%d %.17g %.17g %d %.17g %.17g %.17g %.17g %.17g %.17g "
      "%.17g %.17g\n",
      (int)UPPER_PROFILE_MODE, (double)INITIAL_GAP,
      latest_contact_geometry.gap, latest_contact_geometry.valid,
      initialized.centers[0], initialized.centers[1],
      initialized.center_velocities[0], initialized.center_velocities[1],
      initialized.envelope_overlap_volume, initialized.gas_overlap_volume,
      initialized.nesting_volume, initialized.overfill_volume);
    fclose(audit);

    FILE * fp = fopen("data/anchor_stats.dat", "w");
    if (!fp) {
      perror("data/anchor_stats.dat");
      exit(2);
    }
    fprintf(fp,
      "# t i dt cells lower_center upper_center lower_u upper_u "
      "closing_speed moment_gap plic_gap plic_valid lower_samples "
      "upper_samples anchor_active released_this_step anchor_acceleration "
      "anchor_feedforward effective_mass support_volume buoyancy_force "
      "reaction_force reaction_power constraint_work "
      "release_time release_gap saturation_count lower_gas_drift "
      "lower_oil_drift upper_gas_drift upper_oil_drift nesting_volume "
      "envelope_overlap_volume gas_overlap_volume overfill_volume "
      "maximum_speed umax_x umax_y umax_water_fraction "
      "umax_oil_fraction umax_gas_fraction "
      "lower_outer_interface_velocity_rms "
      "lower_outer_interface_velocity_max lower_gas_interface_velocity_rms "
      "lower_gas_interface_velocity_max upper_outer_interface_velocity_rms "
      "upper_outer_interface_velocity_max upper_gas_interface_velocity_rms "
      "upper_gas_interface_velocity_max projected_volume_change pressure_iterations "
      "scaled_pressure_residual viscosity_iterations viscosity_residual "
      "momentum_transport_mode invalid_cells invalid_faces transport_failed "
      "max_cell_partition_defect max_face_partition_defect "
      "max_vof_bound_defect max_partition_correction "
      "nested_plic_fallbacks nested_plic_fallbacks_max "
      "transport_cpu_seconds lower_aspect upper_aspect kinetic_energy "
      "surface_energy pressure_integral film_mean_pressure "
      "film_max_pressure film_support_volume lower_contact_delta "
      "upper_contact_delta film_cells film_distribution_valid_bins "
      "film_distribution_total_bins film_valid_area_fraction "
      "film_nonpositive_gap_bins film_min_gap film_min_cells "
      "film_p10_cells film_median_cells "
      "film_resolved_area_fraction\n");
    fclose(fp);
    fp = fopen("data/contact_amr_stats.dat", "w");
    if (!fp) {
      perror("data/contact_amr_stats.dat");
      exit(2);
    }
    fprintf(fp,
      "# t i active forced_cells lower_top upper_bottom gap finest_delta\n");
    fclose(fp);
#if ENABLE_STATEFUL_FILM
    fp = fopen("data/film_state_stats.dat", "w");
    if (!fp) {
      perror("data/film_state_stats.dat");
      exit(2);
    }
    fprintf(fp,
      "# t i dt state_status state_valid active handed_off reacquired "
      "activation_count update_count handoff_time handoff_raw_valid_area "
      "handoff_raw_crossed_bins handoff_corridor_at_maxlevel "
      "observation_status observation_valid virtual_kinematics raw_available "
      "raw_valid_area raw_crossed raw_unresolved min_state_gap min_raw_gap "
      "max_state_raw_gap_discrepancy film_dt_limit inventory_before "
      "inventory_after expected_inventory_after inventory_residual edge_flux "
      "reynolds_flux_residual reynolds_continuity_residual "
      "reynolds_pressure_residual traction_valid traction_active "
      "traction_applied "
      "requested_lower_force requested_upper_force requested_net_force "
      "actual_lower_force actual_upper_force actual_net_force "
      "lower_force_residual upper_force_residual net_force_residual "
      "force_audit_tolerance application_roundoff_faces "
      "swamped_increment_faces maximum_application_error "
      "application_force_error_l1 missing_lower_support_bins "
      "missing_upper_support_bins identity_audit_failed force_audit_failed "
      "force_audit_nonfinite actual_net_moment moment_audit_tolerance "
      "radial_first_moment_imbalance radial_first_moment_tolerance "
      "max_bin_force_imbalance step_work "
      "cumulative_work positive_film_work contact_corridor_at_maxlevel "
      "restart_restored restart_state_time restart_continuity_error\n");
    fclose(fp);
    fp = fopen("data/film_state_primitive.dat", "w");
    if (!fp) {
      perror("data/film_state_primitive.dat");
      exit(2);
    }
    fprintf(fp,
      "# t i bin active activation_count update_count reacquisition_count "
      "unresolved_count raw_crossing_count source_time handoff_inventory "
      "inventory cumulative_edge_outflow cumulative_work thickness midpoint "
      "radius width delta face_flux_inner face_flux_outer\n");
    fclose(fp);
#endif
#if UPPER_PROFILE_MODE == 1
    fp = fopen("data/upper_profile_metadata.dat", "w");
    if (!fp) {
      perror("data/upper_profile_metadata.dat");
      exit(2);
    }
    fprintf(fp,
      "profile_file %s\npoints %d\nouter_volume %.17g\n"
      "gas_volume %.17g\noil_volume %.17g\noil_fraction %.17g\n"
      "outer_centroid_local %.17g\naxial_offset %.17g\n"
      "minimum_radial_film %.17g\n",
      upper_profile_path, upper_profile.count,
      upper_profile.outer_volume, upper_profile.gas_volume,
      upper_profile.oil_volume, upper_profile.oil_fraction,
      upper_profile.outer_centroid, upper_profile_offset,
      upper_profile.minimum_radial_film);
    fclose(fp);
#endif
  }
#if COLLISION_FULL_STATE_RESTART
  append_collision_diagnostics(restart_dump_time, restart_dump_iteration);
# if ENABLE_STATEFUL_FILM
  append_film_state_primitive(restart_dump_time, restart_dump_iteration);
# endif
#endif
}

event stability (i++, last)
{
#if ENABLE_STATEFUL_FILM
  if (film_state.active && latest_film_observation.valid) {
    double film_limit =
      axisymmetric_film_runtime_timestep_limit(
      &film_state,
      &latest_film_observation.observation);
    if (isfinite(film_limit) && film_limit > 0.) {
      latest_film_timestep_limit = film_limit;
      dtmax = min(dtmax, film_limit);
    }
    else if (film_limit != HUGE_VAL)
      film_state_failed = true;
  }
#endif
}

event acceleration (i++)
{
  // The i++ event chain is called once at the scheduled terminal time. The
  // physical state already represents END_TIME, so model ledgers must not be
  // advanced into a step which the solver will never take.
  if (collision_terminal_event_time(t, (double)END_TIME))
    return 0;
  CollisionState state = measure_collision_state();
  latest_contact_geometry = measure_axisymmetric_contact_geometry(
    (double)CONTACT_SUPPORT_HALF_WIDTH);
  double release_gap = latest_contact_geometry.valid ?
                       latest_contact_geometry.gap : state.outer_gap;
  CapillaryAnchorObservation observation = {
    .time = t,
    .dt = dt,
    .outer_gap = release_gap,
    .closing_speed = state.center_velocities[0] -
                     state.center_velocities[1],
    .upper_center_velocity = state.center_velocities[1]
  };
  latest_anchor_command = capillary_anchor_update(
    &anchor_state, &anchor_config, observation);

  face vector acceleration_field = a;
  foreach_face(x)
    acceleration_field.x[] += -(double)GRAVITY_MAGNITUDE;
  latest_anchor_load = capillary_anchor_apply_axi(
    acceleration_field, latest_anchor_command.acceleration);
  capillary_anchor_record_reaction(
    &anchor_state, latest_anchor_load.reaction_force,
    state.center_velocities[1], dt);
#if ENABLE_STATEFUL_FILM
  if (!collision_refresh_film_observation())
    film_state_failed = true;
  else if (latest_film_observation.valid) {
      AxisymmetricFilmState state_before_update = film_state;
      AxisymmetricFilmState traction_state = film_state;
      AxisymmetricFilmInput traction_input = {0};
      bool just_handed_off = false;
      if (traction_state.active)
        traction_input = axisymmetric_film_state_reynolds_input(
          &traction_state, &latest_film_observation.observation);
      latest_film_state_result = axisymmetric_film_state_update(
        &film_state_config, &film_state,
        &latest_film_observation.observation, dt);
      switch (latest_film_state_result.status) {
      case AXISYMMETRIC_FILM_STATE_INACTIVE:
      case AXISYMMETRIC_FILM_STATE_REACQUIRED:
        break;
      case AXISYMMETRIC_FILM_STATE_HANDOFF:
        just_handed_off = true;
        traction_state = film_state;
        traction_input = axisymmetric_film_state_reynolds_input(
          &traction_state, &latest_film_observation.observation);
        latest_film_state_result.reynolds_solution =
          axisymmetric_reynolds_solution(&film_state_config.reynolds,
            &traction_input, traction_state.thickness);
        if (!latest_film_state_result.reynolds_solution.valid ||
            latest_film_state_result.reynolds_solution.
              maximum_flux_balance_error >
                film_state_config.flux_balance_tolerance) {
          latest_film_state_result.status =
            AXISYMMETRIC_FILM_STATE_FLUX_SOLVE_FAILED;
          latest_film_state_result.valid = false;
          collision_abort_film_transaction(&state_before_update);
          break;
        }
        // Fall through: handoff traction uses the copied, pre-evolution state.
      case AXISYMMETRIC_FILM_STATE_ACTIVE:
        if (!latest_film_state_result.valid) {
          collision_abort_film_transaction(&state_before_update);
          break;
        }
        latest_film_traction = axisymmetric_matched_film_evaluate(
          &film_state_config.reynolds, &traction_input,
          &latest_film_state_result.reynolds_solution,
          traction_state.cumulative_work, dt);
        if (!latest_film_traction.valid) {
          collision_abort_film_transaction(&state_before_update);
          break;
        }
        {
          const double work_scale = max(DBL_MIN,
            fabs(traction_state.cumulative_work) +
            fabs(latest_film_traction.step_work) +
            fabs(latest_film_traction.cumulative_work));
          const double work_tolerance = 8192.*DBL_EPSILON*work_scale;
          positive_film_work = latest_film_traction.step_work >
              work_tolerance || latest_film_traction.cumulative_work >
              work_tolerance;
          if (positive_film_work) {
            collision_abort_film_transaction(&state_before_update);
            break;
          }
        }
        latest_film_pair_ledger = axisymmetric_film_state_apply_axi(
          acceleration_field, film_owner_mask,
          lower_envelope_fraction, upper_envelope_fraction,
          &latest_film_traction, &traction_state,
          (double)CONTACT_SUPPORT_HALF_WIDTH);
        const int missing_support_bins =
          latest_film_pair_ledger.missing_lower_support_bins +
          latest_film_pair_ledger.missing_upper_support_bins;
        if (!latest_film_pair_ledger.applied ||
            missing_support_bins ||
            latest_film_pair_ledger.identity_audit_failed ||
            latest_film_pair_ledger.force_audit_failed ||
            latest_film_pair_ledger.force_audit_nonfinite) {
          collision_abort_film_transaction(&state_before_update);
          break;
        }
        film_state.cumulative_work =
          latest_film_traction.cumulative_work;
        if (just_handed_off) {
          film_handoff_time = t;
          film_handoff_raw_valid_area =
            latest_film_observation.observation.raw_valid_area_fraction;
          film_handoff_raw_crossed_bins =
            latest_film_observation.observation.raw_crossed_bins;
          film_handoff_corridor_at_maxlevel =
            latest_film_observation.observation.contact_corridor_at_maxlevel;
        }
        break;
      case AXISYMMETRIC_FILM_STATE_NONPOSITIVE_STATE:
      case AXISYMMETRIC_FILM_STATE_FLUX_SOLVE_FAILED:
      case AXISYMMETRIC_FILM_STATE_INVENTORY_MISMATCH:
      case AXISYMMETRIC_FILM_STATE_RESTART_MISMATCH:
      case AXISYMMETRIC_FILM_STATE_INVALID_CONFIG:
      case AXISYMMETRIC_FILM_STATE_INVALID_OBSERVATION:
      default:
        collision_abort_film_transaction(&state_before_update);
        break;
      }
      if (film_state.active)
        latest_film_timestep_limit =
          axisymmetric_film_runtime_timestep_limit(
            &film_state,
            &latest_film_observation.observation);
      else
        latest_film_timestep_limit = HUGE_VAL;
      if (film_state.active &&
          ((!isfinite(latest_film_timestep_limit) &&
            latest_film_timestep_limit != HUGE_VAL) ||
           (isfinite(latest_film_timestep_limit) &&
            latest_film_timestep_limit <= 0.)))
        collision_abort_film_transaction(&state_before_update);
  }
#endif
}

event projection (i++)
{
#if COLLISION_FULL_STATE_RESTART
  if (!restart_pressure_boundary_initialized) {
    p.dirty = pf.dirty = true;
    boundary((scalar *){p, pf});
    restart_pressure_boundary_initialized = true;
  }
#if RESTART_FIELD_AUDIT
  if (restart_projection_audit_count < 2) {
    long invalid_faces = 0, invalid_cells = 0;
    double maximum_face_velocity = 0., maximum_divergence = 0.;
    face vector acceleration_field = a;
    foreach_face(reduction(+:invalid_faces)
                 reduction(max:maximum_face_velocity)) {
      bool valid = isfinite(uf.x[]) && isfinite(alpha.x[]) &&
                   isfinite(acceleration_field.x[]) && isfinite(fm.x[]) &&
                   isfinite(p[]) && isfinite(p[-1]);
      if (!valid) {
        if (invalid_faces == 0)
          fprintf(stderr,
            "[COLLISION-RESTART-FACE-REJECT] pid=%d t=%.17g i=%d "
            "x=%.17g y=%.17g Delta=%.17g uf=%.17g alpha=%.17g "
            "a=%.17g fm=%.17g p_left=%.17g p_right=%.17g\n",
            pid(), t, i, x, y, Delta, uf.x[], alpha.x[],
            acceleration_field.x[], fm.x[], p[-1], p[]), fflush(stderr);
        invalid_faces++;
      }
      else
        maximum_face_velocity = max(maximum_face_velocity, fabs(uf.x[]));
    }
    foreach(reduction(+:invalid_cells) reduction(max:maximum_divergence)) {
      bool valid = isfinite(dt) && dt > 0.;
      double divergence = 0.;
      foreach_dimension() {
        valid = valid && isfinite(uf.x[]) && isfinite(uf.x[1]);
        if (valid)
          divergence += uf.x[1] - uf.x[];
      }
      if (valid) {
        divergence /= dt*Delta;
        valid = isfinite(divergence);
      }
      if (!valid)
        invalid_cells++;
      else
        maximum_divergence = max(maximum_divergence, fabs(divergence));
    }
    if (pid() == 0)
      fprintf(stderr,
        "[COLLISION-RESTART-PROJECTION-AUDIT] t=%.17g i=%d dt=%.17g "
        "invalid_faces=%ld invalid_cells=%ld ufmax=%.17g divmax=%.17g\n",
        t, i, dt, invalid_faces, invalid_cells,
        maximum_face_velocity, maximum_divergence), fflush(stderr);
    restart_projection_audit_count++;
    if (invalid_faces || invalid_cells) {
      restart_projection_input_failed = true;
      fflush(stderr);
      exit(5);
    }
  }
#endif
#endif
}

#if ENABLE_STATEFUL_FILM
static void append_film_state_reject_primitive (double output_time,
                                                int output_iteration)
{
  char path[128];
  snprintf(path, sizeof(path), "data/film_state_reject_rank%05d.dat", pid());
  FILE * fp = fopen(path, "w");
  if (!fp) {
    perror(path);
    return;
  }
  fprintf(fp,
    "# t i pid bin state_active radius width state_gap raw_gap delta "
    "lower_position upper_position lower_vn upper_vn lower_vt upper_vt "
    "raw_relative_vn raw_lower_vt raw_upper_vt pressure residual "
    "inner_flux outer_flux\n");
  AxisymmetricFilmInput failed_input = film_state.active ?
    axisymmetric_film_state_reynolds_input(
      &film_state, &latest_film_observation.observation) :
    latest_film_observation.observation.raw;
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++)
    fprintf(fp,
      "%.17g %d %d %d %d %.17g %.17g %.17g %.17g %.17g %.17g "
      "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g "
      "%.17g %.17g %.17g\n",
      output_time, output_iteration, pid(), bin, film_state.active,
      failed_input.radius[bin], failed_input.width[bin],
      failed_input.gap[bin],
      latest_film_observation.observation.raw.gap[bin],
      failed_input.delta[bin], failed_input.lower_position[bin],
      failed_input.upper_position[bin],
      latest_film_observation.observation.lower_normal_velocity[bin],
      latest_film_observation.observation.upper_normal_velocity[bin],
      latest_film_observation.observation.lower_tangential_velocity[bin],
      latest_film_observation.observation.upper_tangential_velocity[bin],
      latest_film_observation.observation.raw.relative_normal_velocity[bin],
      latest_film_observation.observation.raw.lower_tangential_velocity[bin],
      latest_film_observation.observation.raw.upper_tangential_velocity[bin],
      (double) latest_film_state_result.reynolds_solution.pressure[bin],
      latest_film_state_result.reynolds_solution.residual[bin],
      latest_film_state_result.reynolds_solution.face_flux[bin],
      latest_film_state_result.reynolds_solution.face_flux[bin + 1]);
  fclose(fp);
}
#endif

event film_state_guard (i++, last)
{
#if ENABLE_STATEFUL_FILM
  if (film_state_failed) {
    append_film_state_reject_primitive(t, i);
    fprintf(stderr,
      "[FILM-STATE-REJECT-RANK] pid=%d t=%.17g i=%d "
      "state_status=%d observation_reason=%d "
      "raw_valid_area=%.17g raw_crossed_bins=%d "
      "contact_corridor_at_maxlevel=%d\n",
      pid(), t, i, latest_film_state_result.status,
      latest_film_state_result.observation_reason,
      latest_film_observation.observation.raw_valid_area_fraction,
      latest_film_observation.observation.raw_crossed_bins,
      latest_film_observation.observation.contact_corridor_at_maxlevel);
    fflush(stderr);
    if (pid() == 0) {
      fprintf(stderr,
        "stateful film closure rejected step at t=%.17g i=%d "
        "state_status=%d observation_status=%d state_active=%d "
        "state_valid=%d observation_valid=%d raw_ready=%d "
        "raw_crossed_bins=%d nonpositive_bin=%d nonpositive_h=%.17g "
        "missing_lower_support_bins=%d missing_upper_support_bins=%d "
        "identity_audit_failed=%d force_audit_failed=%d "
        "force_audit_nonfinite=%d requested_lower_force=%.17g "
        "requested_upper_force=%.17g measured_lower_force=%.17g "
        "measured_upper_force=%.17g lower_force_residual=%.17g "
        "upper_force_residual=%.17g force_audit_tolerance=%.17g "
        "application_roundoff_faces=%ld swamped_increment_faces=%ld "
        "maximum_application_error=%.17g "
        "application_force_error_l1=%.17g "
        "positive_film_work=%d inventory_residual=%.17g "
        "flux_residual=%.17g continuity_residual=%.17g "
        "pressure_residual=%.17g dt_limit=%.17g\n",
        t, i, latest_film_state_result.status,
        latest_film_observation.status, film_state.active,
        latest_film_state_result.valid, latest_film_observation.valid,
        latest_film_samples.input_ready,
        latest_film_samples.nonpositive_gap_bins,
        latest_film_state_result.nonpositive_bin,
        latest_film_state_result.nonpositive_thickness,
        latest_film_pair_ledger.missing_lower_support_bins,
        latest_film_pair_ledger.missing_upper_support_bins,
        latest_film_pair_ledger.identity_audit_failed,
        latest_film_pair_ledger.force_audit_failed,
        latest_film_pair_ledger.force_audit_nonfinite,
        latest_film_pair_ledger.requested_lower_axial_force,
        latest_film_pair_ledger.requested_upper_axial_force,
        latest_film_pair_ledger.measured_lower_axial_force,
        latest_film_pair_ledger.measured_upper_axial_force,
        latest_film_pair_ledger.lower_force_residual,
        latest_film_pair_ledger.upper_force_residual,
        latest_film_pair_ledger.force_audit_tolerance,
        latest_film_pair_ledger.application_roundoff_faces,
        latest_film_pair_ledger.swamped_increment_faces,
        latest_film_pair_ledger.maximum_application_error,
        latest_film_pair_ledger.application_force_error_l1,
        positive_film_work,
        latest_film_state_result.inventory_closure_error,
        latest_film_state_result.reynolds_solution.maximum_flux_balance_error,
        latest_film_state_result.reynolds_solution.
          maximum_continuity_balance_error,
        latest_film_state_result.reynolds_solution.
          maximum_pressure_balance_error,
        latest_film_timestep_limit);
      fflush(stderr);
    }
    return 1;
  }
#endif
}

static void append_collision_diagnostics (double output_time,
                                          int output_iteration)
{
  CollisionState state = measure_collision_state();
  latest_contact_geometry = measure_axisymmetric_contact_geometry(
    (double)CONTACT_SUPPORT_HALF_WIDTH);
  AxisymmetricCollisionObservables observables =
    measure_axisymmetric_collision_observables(
      latest_contact_geometry, (double)CONTACT_SUPPORT_HALF_WIDTH,
      (double)SIGMA_OG, (double)SIGMA_OW);
  AxisymmetricVelocityAudit velocity_audit =
    measure_axisymmetric_velocity_audit();
  double film_cells = contact_resolution_cells(
    latest_contact_geometry.gap,
    latest_contact_geometry.lower_max_delta,
    latest_contact_geometry.upper_max_delta);
  AxisymmetricFilmDistribution film_distribution =
    measure_axisymmetric_film_distribution(
      (double)CONTACT_SUPPORT_HALF_WIDTH, 4.);
  FilmDistributionStatistics film_stats = film_distribution.statistics;
  double gas_drifts[2], oil_drifts[2];
  for (int identity = 0; identity < 2; identity++) {
    gas_drifts[identity] = fabs(
      state.gas_volumes[identity] - initial_gas_volumes[identity])/
      max(initial_gas_volumes[identity], 1.e-300);
    oil_drifts[identity] = fabs(
      state.oil_volumes[identity] - initial_oil_volumes[identity])/
      max(initial_oil_volumes[identity], 1.e-300);
  }
  if (pid() == 0) {
    FILE * fp = fopen("data/anchor_stats.dat", "a");
    fprintf(fp,
      "%.17g %d %.17g %ld "
      "%.17g %.17g %.17g %.17g %.17g "
      "%.17g %.17g %d %d %d %d %d %.17g "
      "%.17g %.17g %.17g %.17g %.17g %.17g %.17g "
      "%.17g %.17g %lu "
      "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g "
      "%.17g ",
      output_time, output_iteration, dt, grid->tn,
      state.centers[0], state.centers[1],
      state.center_velocities[0], state.center_velocities[1],
      state.center_velocities[0] - state.center_velocities[1],
      state.outer_gap, latest_contact_geometry.gap,
      latest_contact_geometry.valid,
      latest_contact_geometry.lower_samples,
      latest_contact_geometry.upper_samples,
      latest_anchor_command.active,
      latest_anchor_command.released_this_step,
      latest_anchor_command.acceleration,
      anchor_config.feedforward_acceleration,
      latest_anchor_load.effective_mass,
      latest_anchor_load.support_volume,
      estimated_buoyancy_force,
      latest_anchor_load.reaction_force,
      anchor_state.last_power, anchor_state.cumulative_work,
      anchor_state.release_time, anchor_state.release_gap,
      anchor_state.saturation_count,
      gas_drifts[0], oil_drifts[0], gas_drifts[1], oil_drifts[1],
      state.nesting_volume, state.envelope_overlap_volume,
      state.gas_overlap_volume, state.overfill_volume,
      state.maximum_speed);
    fprintf(fp, "%.17g %.17g %.17g %.17g %.17g ",
      velocity_audit.maximum_x, velocity_audit.maximum_y,
      velocity_audit.maximum_water_fraction,
      velocity_audit.maximum_oil_fraction,
      velocity_audit.maximum_gas_fraction);
    fprintf(fp, "%.17g %.17g %.17g %.17g ",
      velocity_audit.lower_outer_interface_velocity_rms,
      velocity_audit.lower_outer_interface_velocity_max,
      velocity_audit.lower_gas_interface_velocity_rms,
      velocity_audit.lower_gas_interface_velocity_max);
    fprintf(fp, "%.17g %.17g %.17g %.17g ",
      velocity_audit.upper_outer_interface_velocity_rms,
      velocity_audit.upper_outer_interface_velocity_max,
      velocity_audit.upper_gas_interface_velocity_rms,
      velocity_audit.upper_gas_interface_velocity_max);
    fprintf(fp, "%.17g %d %.17g ",
      state.projected_volume_change, mgp.i, mgp.resa*sq(dt));
    fprintf(fp,
      "%d %.17g %d %ld %ld %d %.17g %.17g %.17g %.17g %ld %ld %.17g "
      "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g "
      "%.17g %.17g %d %d %.17g %d %.17g %.17g %.17g %.17g %.17g\n",
      mgu.i, mgu.resa, (int)MOMENTUM_TRANSPORT_MODE,
      dual_momentum_invalid_cell_count,
      dual_momentum_invalid_face_count,
      dual_momentum_transport_failed,
      dual_momentum_max_cell_partition_defect,
      dual_momentum_max_face_partition_defect,
      dual_momentum_max_vof_bound_defect,
      dual_momentum_max_partition_correction,
      dual_momentum_nested_plic_fallback_count,
      dual_momentum_nested_plic_fallback_count_max,
      dual_momentum_last_transport_cpu_seconds,
      state.moment_aspects[0], state.moment_aspects[1],
      observables.kinetic_energy, observables.surface_energy,
      observables.pressure_integral, observables.film_mean_pressure,
      observables.film_max_pressure, observables.film_support_volume,
      latest_contact_geometry.lower_max_delta,
      latest_contact_geometry.upper_max_delta, film_cells,
      film_stats.valid_bins, film_stats.total_bins,
      film_stats.valid_area_fraction,
      film_distribution.nonpositive_gap_bins, film_stats.minimum_gap,
      film_stats.minimum_cells, film_stats.p10_cells,
      film_stats.median_cells, film_stats.resolved_area_fraction);
    fclose(fp);
#if ENABLE_STATEFUL_FILM
    fp = fopen("data/film_state_stats.dat", "a");
    if (!fp) {
      perror("data/film_state_stats.dat");
      exit(2);
    }
    fprintf(fp,
      "%.17g %d %.17g %d %d %d %d %d %lu %lu %.17g %.17g %d %d ",
      output_time, output_iteration, dt, latest_film_state_result.status,
      latest_film_state_result.valid, film_state.active,
      latest_film_state_result.handed_off,
      latest_film_state_result.reacquired,
      film_state.activation_count, film_state.update_count,
      film_handoff_time, film_handoff_raw_valid_area,
      film_handoff_raw_crossed_bins,
      film_handoff_corridor_at_maxlevel);
    fprintf(fp, "%d %d %d %d %.17g %d %d ",
      latest_film_observation.status, latest_film_observation.valid,
      latest_film_observation.using_virtual_kinematics,
      latest_film_observation.observation.raw_input_available,
      latest_film_observation.observation.raw_valid_area_fraction,
      latest_film_observation.observation.raw_crossed_bins,
      latest_film_state_result.raw_geometry_unresolved);
    fprintf(fp,
      "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g "
      "%.17g %.17g ",
      latest_film_state_result.minimum_state_thickness,
      latest_film_state_result.minimum_raw_gap,
      latest_film_state_result.maximum_state_raw_gap_discrepancy,
      latest_film_timestep_limit,
      latest_film_state_result.inventory_before,
      latest_film_state_result.inventory_after,
      latest_film_state_result.expected_inventory_after,
      latest_film_state_result.inventory_closure_error,
      film_state.face_flux[AXISYMMETRIC_FILM_STATE_BINS],
      latest_film_state_result.reynolds_solution.maximum_flux_balance_error,
      latest_film_state_result.reynolds_solution.
        maximum_continuity_balance_error,
      latest_film_state_result.reynolds_solution.
        maximum_pressure_balance_error);
    fprintf(fp, "%d %d %d ",
      latest_film_traction.valid, latest_film_traction.active,
      latest_film_pair_ledger.applied);
    fprintf(fp,
      "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g ",
      latest_film_pair_ledger.requested_lower_axial_force,
      latest_film_pair_ledger.requested_upper_axial_force,
      latest_film_pair_ledger.requested_net_axial_force,
      latest_film_pair_ledger.actual_lower_axial_force,
      latest_film_pair_ledger.actual_upper_axial_force,
      latest_film_pair_ledger.actual_net_axial_force,
      latest_film_pair_ledger.lower_force_residual,
      latest_film_pair_ledger.upper_force_residual,
      latest_film_pair_ledger.net_force_residual,
      latest_film_pair_ledger.force_audit_tolerance);
    fprintf(fp, "%ld %ld %.17g %.17g %d %d %d %d %d ",
      latest_film_pair_ledger.application_roundoff_faces,
      latest_film_pair_ledger.swamped_increment_faces,
      latest_film_pair_ledger.maximum_application_error,
      latest_film_pair_ledger.application_force_error_l1,
      latest_film_pair_ledger.missing_lower_support_bins,
      latest_film_pair_ledger.missing_upper_support_bins,
      latest_film_pair_ledger.identity_audit_failed,
      latest_film_pair_ledger.force_audit_failed,
      latest_film_pair_ledger.force_audit_nonfinite);
    fprintf(fp,
      "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %d %d %d ",
      latest_film_pair_ledger.actual_net_axial_moment,
      latest_film_pair_ledger.force_audit_tolerance*
        (double)CONTACT_SUPPORT_HALF_WIDTH,
      latest_film_pair_ledger.radial_first_moment_imbalance,
      latest_film_pair_ledger.radial_first_moment_tolerance,
      latest_film_pair_ledger.maximum_bin_force_imbalance,
      latest_film_traction.step_work, film_state.cumulative_work,
      positive_film_work, latest_contact_corridor_at_maxlevel,
      film_restart_restored);
    fprintf(fp, "%.17g %.17g\n",
#if COLLISION_FULL_STATE_RESTART
      (double)COLLISION_RESTART_FILM_STATE_TIME,
#else
      NAN,
#endif
      film_restart_continuity_error);
    fclose(fp);
#endif
    fp = fopen("data/contact_amr_stats.dat", "a");
    if (!fp) {
      perror("data/contact_amr_stats.dat");
      exit(2);
    }
    double finest_delta = (double)DOMAIN_LENGTH/(1 << MAXLEVEL);
    ContactAmrPolicy contact_policy = collision_contact_amr_policy();
    ContactAmrGeometry contact_geometry =
#if ENABLE_STATEFUL_FILM
      film_state.active ? axisymmetric_film_state_contact_geometry(&film_state) :
#endif
      collision_raw_contact_geometry(latest_contact_geometry);
    fprintf(fp, "%.17g %d %d %ld %.17g %.17g %.17g %.17g\n",
      output_time, output_iteration,
      contact_amr_is_active(&contact_policy, &contact_geometry),
      latest_contact_amr_forced_cells,
      contact_geometry.virtual_bounds_present ?
        contact_geometry.axial_minimum : contact_geometry.lower_top,
      contact_geometry.virtual_bounds_present ?
        contact_geometry.axial_maximum : contact_geometry.upper_bottom,
      contact_amr_geometry_minimum_gap(&contact_geometry),
      finest_delta);
    fclose(fp);
  }
}

#if COLLISION_FULL_STATE_RESTART
event diagnostics (t = RESTART_SOURCE_TIME;
                   t <= END_TIME; t += OUTPUT_INTERVAL)
#else
event diagnostics (t = 0.; t <= END_TIME; t += OUTPUT_INTERVAL)
#endif
{
  append_collision_diagnostics(t, i);
}

#if ENABLE_STATEFUL_FILM
static void append_film_state_primitive (double output_time,
                                         int output_iteration)
{
  if (pid() != 0)
    return;
  const double time_tolerance =
    128.*DBL_EPSILON*max(1., fabs(output_time));
  if (isfinite(last_film_primitive_time) &&
      fabs(output_time - last_film_primitive_time) <= time_tolerance)
    return;
  FILE * fp = fopen("data/film_state_primitive.dat", "a");
  if (!fp) {
    perror("data/film_state_primitive.dat");
    exit(2);
  }
  for (int bin = 0; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++)
    fprintf(fp,
      "%.17g %d %d %d %lu %lu %lu %lu %lu %.17g %.17g %.17g %.17g "
      "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
      output_time, output_iteration, bin, film_state.active,
      film_state.activation_count,
      film_state.update_count,
      film_state.reacquisition_consecutive_observations,
      film_state.cumulative_unresolved_observations,
      film_state.cumulative_raw_crossing_observations,
      film_state.source_time, film_state.handoff_inventory,
      film_state.inventory, film_state.cumulative_edge_outflow,
      film_state.cumulative_work, film_state.thickness[bin],
      film_state.midpoint[bin], film_state.radius[bin],
      film_state.width[bin], film_state.delta[bin],
      film_state.face_flux[bin], film_state.face_flux[bin + 1]);
  fclose(fp);
  last_film_primitive_time = output_time;
}
#endif

#if COLLISION_FULL_STATE_RESTART
event snapshots (t = RESTART_SOURCE_TIME;
                 t <= END_TIME; t += DUMP_INTERVAL)
#else
event snapshots (t = 0.; t <= END_TIME; t += DUMP_INTERVAL)
#endif
{
  char path[128];
  sprintf(path, "dumps/anchor-%09.6f", t);
  dump(file = path);
#if ENABLE_STATEFUL_FILM
  append_film_state_primitive(t, i);
#endif
}

event adapt (i++)
{
#if TREE
  adapt_wavelet((scalar *){
      lower_gas_fraction, lower_envelope_fraction,
      upper_gas_fraction, upper_envelope_fraction, u
    },
    (double[]){
      (double)FRACTION_TOLERANCE, (double)FRACTION_TOLERANCE,
      (double)FRACTION_TOLERANCE, (double)FRACTION_TOLERANCE,
      (double)VELOCITY_TOLERANCE, (double)VELOCITY_TOLERANCE
    }, MAXLEVEL, MINLEVEL);
  AxisymmetricContactGeometry contact =
    measure_axisymmetric_contact_geometry((double)CONTACT_SUPPORT_HALF_WIDTH);
  ContactAmrGeometry geometry =
#if ENABLE_STATEFUL_FILM
    film_state.active ? axisymmetric_film_state_contact_geometry(&film_state) :
#endif
    collision_raw_contact_geometry(contact);
  latest_contact_amr_forced_cells =
    collision_refine_contact_corridor(&geometry);
#if ENABLE_STATEFUL_FILM
  latest_contact_corridor_at_maxlevel =
    collision_contact_corridor_at_maxlevel(&geometry);
#endif
#endif
}

event finish (t = END_TIME)
{
  dump(file = "dumps/final");
#if ENABLE_STATEFUL_FILM
  append_film_state_primitive(t, i);
#endif
  if (pid() == 0)
    fprintf(stderr,
      "# FINISH t=%.17g i=%d anchor_status=%d release_time=%.17g\n",
      t, i, anchor_state.status, anchor_state.release_time);
}
