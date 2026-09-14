#ifndef AXISYMMETRIC_FILM_STATE_AXI_H
#define AXISYMMETRIC_FILM_STATE_AXI_H

#include <float.h>
#include <math.h>
#include <stdbool.h>

#include "axisymmetric_film_samples.h"
#include "axisymmetric_film_state.h"
#include "contact_amr_policy.h"

/* This adapter owns the Basilisk-specific support integration.  Include it
   before any non-adapter inclusion of axisymmetric_matched_film_axi.h. */
#ifndef AXISYMMETRIC_MATCHED_FILM_BASILISK_ADAPTER
# define AXISYMMETRIC_MATCHED_FILM_BASILISK_ADAPTER 1
# define AXISYMMETRIC_FILM_STATE_AXI_DEFINED_MATCHED_ADAPTER 1
#endif
#include "axisymmetric_matched_film_axi.h"
#ifdef AXISYMMETRIC_FILM_STATE_AXI_DEFINED_MATCHED_ADAPTER
# undef AXISYMMETRIC_MATCHED_FILM_BASILISK_ADAPTER
# undef AXISYMMETRIC_FILM_STATE_AXI_DEFINED_MATCHED_ADAPTER
#endif

enum {
  AXISYMMETRIC_FILM_STATE_AXI_OBSERVATION_OK = 0,
  AXISYMMETRIC_FILM_STATE_AXI_INVALID_ARGUMENT = 1,
  AXISYMMETRIC_FILM_STATE_AXI_INVALID_ACTIVE_STATE = 2,
  AXISYMMETRIC_FILM_STATE_AXI_RAW_INPUT_UNAVAILABLE = 3,
  AXISYMMETRIC_FILM_STATE_AXI_RAW_KINEMATICS_INVALID = 4,
  AXISYMMETRIC_FILM_STATE_AXI_VIRTUAL_SAMPLE_INVALID = 5
};
typedef int AxisymmetricFilmStateAxiObservationStatus;

typedef struct {
  bool valid;
  bool using_virtual_kinematics;
  AxisymmetricFilmStateAxiObservationStatus status;
  int invalid_bin;
  int invalid_side;
  int invalid_component;
  double invalid_value;
  int sampled_virtual_points;
  AxisymmetricFilmObservation observation;
} AxisymmetricFilmStateAxiObservation;

enum {
  AXISYMMETRIC_FILM_STATE_AXI_OWNER_LOWER = -1,
  AXISYMMETRIC_FILM_STATE_AXI_OWNER_NONE = 0,
  AXISYMMETRIC_FILM_STATE_AXI_OWNER_UPPER = 1
};

#define AXISYMMETRIC_FILM_STATE_AXI_IDENTITY_BOUND_TOLERANCE \
  (1024.*DBL_EPSILON)
#define AXISYMMETRIC_FILM_STATE_AXI_FORCE_AUDIT_RELATIVE_TOLERANCE \
  (8192.*DBL_EPSILON)
#define AXISYMMETRIC_FILM_STATE_AXI_RADIAL_MOMENT_RELATIVE_TOLERANCE .02

typedef struct {
  bool applied;
  bool identity_audit_failed;
  bool force_audit_failed;
  bool force_audit_nonfinite;
  double requested_lower_axial_force;
  double requested_upper_axial_force;
  double requested_net_axial_force;
  double actual_lower_axial_force;
  double actual_upper_axial_force;
  double actual_net_axial_force;
  double measured_lower_axial_force;
  double measured_upper_axial_force;
  double measured_net_axial_force;
  double lower_force_residual;
  double upper_force_residual;
  double net_force_residual;
  double force_audit_tolerance;
  double lower_radial_force_moment;
  double upper_radial_force_moment;
  double actual_net_axial_moment;
  double radial_first_moment_imbalance;
  double radial_first_moment_tolerance;
  double lower_support_volume;
  double upper_support_volume;
  double overlap_candidate_volume;
  double potential_local_cancellation_force;
  double maximum_bin_force_imbalance;
  double maximum_identity_bound_defect;
  double maximum_application_error;
  double application_force_error_l1;
  int missing_lower_support_bins;
  int missing_upper_support_bins;
  int invalid_requested_force_bins;
  long lower_owned_faces;
  long upper_owned_faces;
  long overlap_candidate_faces;
  long identity_tie_faces;
  long distance_tie_faces;
  long invalid_identity_faces;
  long invalid_property_faces;
  long identity_free_candidate_faces;
  long application_roundoff_faces;
  long swamped_increment_faces;
} AxisymmetricFilmStateAxiPairLedger;

static inline AxisymmetricFilmStateAxiObservation
axisymmetric_film_state_axi_observation_result (void)
{
  AxisymmetricFilmStateAxiObservation result = {0};
  result.invalid_bin = -1;
  result.invalid_side = -1;
  result.invalid_component = -1;
  result.invalid_value = NAN;
  return result;
}

static inline void axisymmetric_film_state_axi_copy_raw_geometry (
  AxisymmetricFilmObservation * observation,
  const AxisymmetricFilmSamples * raw_samples,
  bool contact_corridor_at_maxlevel)
{
  if (!observation)
    return;
  observation->time = t;
  observation->contact_corridor_at_maxlevel =
    contact_corridor_at_maxlevel;
  observation->raw_input_available = raw_samples &&
    raw_samples->bins == AXISYMMETRIC_FILM_STATE_BINS;
  if (!observation->raw_input_available)
    return;
  observation->raw_valid_area_fraction = raw_samples->valid_area_fraction;
  observation->raw_crossed_bins = raw_samples->nonpositive_gap_bins;
  observation->raw = raw_samples->input;
}

static inline bool axisymmetric_film_state_axi_sample_is_finite (
  double value)
{
  return value != nodata && isfinite(value);
}

static inline bool axisymmetric_film_state_axi_project_raw_kinematics (
  const AxisymmetricFilmSample * sample, double * values)
{
  if (!sample || !values || !sample->valid)
    return false;
  values[0] = sample->lower_velocity.x*sample->normal_x +
    sample->lower_velocity.r*sample->normal_r;
  values[1] = sample->upper_velocity.x*sample->normal_x +
    sample->upper_velocity.r*sample->normal_r;
  values[2] = -sample->lower_velocity.x*sample->normal_r +
    sample->lower_velocity.r*sample->normal_x;
  values[3] = -sample->upper_velocity.x*sample->normal_r +
    sample->upper_velocity.r*sample->normal_x;
  return isfinite(values[0]) && isfinite(values[1]) &&
    isfinite(values[2]) && isfinite(values[3]);
}

static inline bool axisymmetric_film_state_axi_raw_kinematics (
  const AxisymmetricFilmSamples * samples, int bin, double * values)
{
  if (!samples || !values || bin < 0 ||
      bin >= AXISYMMETRIC_FILM_STATE_BINS)
    return false;
  if (axisymmetric_film_state_axi_project_raw_kinematics(
        &samples->bin[bin], values))
    return true;

  int lower = bin - 1, upper = bin + 1;
  while (lower >= 0 && !samples->bin[lower].valid)
    lower--;
  while (upper < AXISYMMETRIC_FILM_STATE_BINS &&
         !samples->bin[upper].valid)
    upper++;
  if (lower < 0 || upper >= AXISYMMETRIC_FILM_STATE_BINS)
    return false;

  double lower_values[4], upper_values[4];
  if (!axisymmetric_film_state_axi_project_raw_kinematics(
        &samples->bin[lower], lower_values) ||
      !axisymmetric_film_state_axi_project_raw_kinematics(
        &samples->bin[upper], upper_values))
    return false;
  for (int component = 0; component < 4; component++) {
    values[component] = axisymmetric_film_samples_linear_value(
      samples->bin[lower].radius, lower_values[component],
      samples->bin[upper].radius, upper_values[component],
      samples->bin[bin].radius);
    if (!isfinite(values[component]))
      return false;
  }
  return true;
}

static inline bool axisymmetric_film_state_axi_state_is_ready (
  const AxisymmetricFilmState * state)
{
  if (!axisymmetric_film_state_active_state_is_valid(state) ||
      !axisymmetric_film_state_active_ledger_is_valid(state))
    return false;
  const double reconstructed_inventory =
    axisymmetric_film_state_inventory(state);
  const double ledger_inventory = state->handoff_inventory -
    state->cumulative_edge_outflow;
  return axisymmetric_film_state_close(reconstructed_inventory,
      state->inventory) && axisymmetric_film_state_close(state->inventory,
      ledger_inventory);
}

static inline AxisymmetricFilmStateAxiObservation
axisymmetric_film_state_observe_axi (
  const AxisymmetricFilmState * state,
  const AxisymmetricFilmSamples * raw_samples,
  bool contact_corridor_at_maxlevel)
{
  enum { bins = AXISYMMETRIC_FILM_STATE_BINS };
  AxisymmetricFilmStateAxiObservation result =
    axisymmetric_film_state_axi_observation_result();
  axisymmetric_film_state_axi_copy_raw_geometry(&result.observation,
    raw_samples, contact_corridor_at_maxlevel);

  if (!state) {
    result.status = AXISYMMETRIC_FILM_STATE_AXI_INVALID_ARGUMENT;
    return result;
  }

  if (!state->active) {
    if (!raw_samples || !raw_samples->input_ready ||
        raw_samples->bins != bins) {
      result.status = AXISYMMETRIC_FILM_STATE_AXI_RAW_INPUT_UNAVAILABLE;
      return result;
    }
    for (int bin = 0; bin < bins; bin++) {
      double kinematics[4];
      if (!axisymmetric_film_state_axi_raw_kinematics(raw_samples, bin,
            kinematics)) {
        result.status = AXISYMMETRIC_FILM_STATE_AXI_RAW_KINEMATICS_INVALID;
        result.invalid_bin = bin;
        return result;
      }
      result.observation.lower_normal_velocity[bin] = kinematics[0];
      result.observation.upper_normal_velocity[bin] = kinematics[1];
      result.observation.lower_tangential_velocity[bin] = kinematics[2];
      result.observation.upper_tangential_velocity[bin] = kinematics[3];
    }
    result.valid = true;
    result.status = AXISYMMETRIC_FILM_STATE_AXI_OBSERVATION_OK;
    return result;
  }

  if (!axisymmetric_film_state_axi_state_is_ready(state)) {
    result.status = AXISYMMETRIC_FILM_STATE_AXI_INVALID_ACTIVE_STATE;
    return result;
  }

  coord virtual_points[2*bins];
  for (int bin = 0; bin < bins; bin++) {
    const double lower_position = state->midpoint[bin] -
      .5*state->thickness[bin];
    const double upper_position = state->midpoint[bin] +
      .5*state->thickness[bin];
    virtual_points[2*bin] = (coord){lower_position, state->radius[bin]};
    virtual_points[2*bin + 1] =
      (coord){upper_position, state->radius[bin]};
  }

  double virtual_velocity[4*bins];
  for (int value = 0; value < 4*bins; value++)
    virtual_velocity[value] = nodata;
  interpolate_array((scalar *){u.x, u.y}, virtual_points, 2*bins,
                    virtual_velocity, true);
  result.using_virtual_kinematics = true;
  result.sampled_virtual_points = 2*bins;
  for (int bin = 0; bin < bins; bin++) {
    for (int side = 0; side < 2; side++)
      for (int component = 0; component < 2; component++) {
        const int index = 4*bin + 2*side + component;
        if (!axisymmetric_film_state_axi_sample_is_finite(
              virtual_velocity[index])) {
          result.status =
            AXISYMMETRIC_FILM_STATE_AXI_VIRTUAL_SAMPLE_INVALID;
          result.invalid_bin = bin;
          result.invalid_side = side;
          result.invalid_component = component;
          result.invalid_value = virtual_velocity[index];
          return result;
        }
      }
    result.observation.lower_normal_velocity[bin] =
      virtual_velocity[4*bin];
    result.observation.lower_tangential_velocity[bin] =
      virtual_velocity[4*bin + 1];
    result.observation.upper_normal_velocity[bin] =
      virtual_velocity[4*bin + 2];
    result.observation.upper_tangential_velocity[bin] =
      virtual_velocity[4*bin + 3];
  }
  result.valid = true;
  result.status = AXISYMMETRIC_FILM_STATE_AXI_OBSERVATION_OK;
  return result;
}

static inline double axisymmetric_film_state_axi_identity_fraction (
  double value)
{
  return fmin(1., fmax(0., value));
}

static inline double axisymmetric_film_state_axi_identity_bound_defect (
  double value)
{
  if (!isfinite(value))
    return HUGE_VAL;
  return max(max(0., -value), max(0., value - 1.));
}

static inline double axisymmetric_film_state_axi_support_weight (
  const AxisymmetricFilmState * state, int bin, int side,
  double axial_coordinate, double radial_coordinate, double delta)
{
  const double interface_position = state->midpoint[bin] +
    side*.5*state->thickness[bin];
  return axisymmetric_matched_film_face_kernel(
      radial_coordinate, state->radius[bin], delta)*
    axisymmetric_matched_film_face_kernel(
      axial_coordinate, interface_position, delta);
}

static inline AxisymmetricFilmStateAxiPairLedger
axisymmetric_film_state_apply_axi (
  face vector acceleration_field,
  face vector owner_mask,
  scalar lower_envelope_identity,
  scalar upper_envelope_identity,
  const AxisymmetricMatchedFilmResult * matched_result,
  const AxisymmetricFilmState * state,
  double support_half_width)
{
  enum { bins = AXISYMMETRIC_FILM_STATE_BINS };
  AxisymmetricFilmStateAxiPairLedger ledger = {0};
  foreach_face()
    owner_mask.x[] = AXISYMMETRIC_FILM_STATE_AXI_OWNER_NONE;

  if (matched_result)
    for (int bin = 0; bin < bins; bin++) {
      if (!isfinite(matched_result->ring_force[bin]) ||
          matched_result->ring_force[bin] < 0.) {
        ledger.invalid_requested_force_bins++;
        continue;
      }
      ledger.requested_lower_axial_force -= matched_result->ring_force[bin];
      ledger.requested_upper_axial_force += matched_result->ring_force[bin];
    }
  ledger.requested_net_axial_force =
    ledger.requested_lower_axial_force +
    ledger.requested_upper_axial_force;

  if (!matched_result || !matched_result->valid || !matched_result->active ||
      ledger.invalid_requested_force_bins ||
      !state ||
      !axisymmetric_film_state_axi_state_is_ready(state) ||
      !isfinite(support_half_width) || support_half_width <= 0.)
    return ledger;

  double candidate_lower_support[bins], candidate_upper_support[bins];
  double lower_support[bins], upper_support[bins];
  for (int bin = 0; bin < bins; bin++) {
    candidate_lower_support[bin] = candidate_upper_support[bin] = 0.;
    lower_support[bin] = upper_support[bin] = 0.;
  }

  double maximum_identity_bound_defect = 0.;
  long invalid_identity_faces = 0, invalid_property_faces = 0;
  long identity_free_candidate_faces = 0;
  foreach_face(x, reduction(+:candidate_lower_support[:bins])
                  reduction(+:candidate_upper_support[:bins])
                  reduction(max:maximum_identity_bound_defect)
                  reduction(+:invalid_identity_faces)
                  reduction(+:invalid_property_faces)
                  reduction(+:identity_free_candidate_faces)) {
    if (y < 0. || y >= support_half_width + Delta)
      continue;

    bool geometric_candidate = false;
    double lower_geometric_weight[bins], upper_geometric_weight[bins];
    for (int bin = 0; bin < bins; bin++) {
      lower_geometric_weight[bin] =
        axisymmetric_film_state_axi_support_weight(state, bin,
          AXISYMMETRIC_FILM_STATE_AXI_OWNER_LOWER, x, y, Delta);
      upper_geometric_weight[bin] =
        axisymmetric_film_state_axi_support_weight(state, bin,
          AXISYMMETRIC_FILM_STATE_AXI_OWNER_UPPER, x, y, Delta);
      geometric_candidate = geometric_candidate ||
        lower_geometric_weight[bin] > 0. ||
        upper_geometric_weight[bin] > 0.;
    }
    if (!geometric_candidate)
      continue;

    const double lower_identity_raw = .5*(lower_envelope_identity[] +
      lower_envelope_identity[-1]);
    const double upper_identity_raw = .5*(upper_envelope_identity[] +
      upper_envelope_identity[-1]);
    const double lower_defect =
      axisymmetric_film_state_axi_identity_bound_defect(lower_identity_raw);
    const double upper_defect =
      axisymmetric_film_state_axi_identity_bound_defect(upper_identity_raw);
    maximum_identity_bound_defect = max(maximum_identity_bound_defect,
      max(lower_defect, upper_defect));
    if (!isfinite(lower_identity_raw) || !isfinite(upper_identity_raw) ||
        lower_defect > AXISYMMETRIC_FILM_STATE_AXI_IDENTITY_BOUND_TOLERANCE ||
        upper_defect > AXISYMMETRIC_FILM_STATE_AXI_IDENTITY_BOUND_TOLERANCE) {
      invalid_identity_faces++;
      continue;
    }

    const double lower_identity =
      axisymmetric_film_state_axi_identity_fraction(lower_identity_raw);
    const double upper_identity =
      axisymmetric_film_state_axi_identity_fraction(upper_identity_raw);
    if (lower_identity == 0. && upper_identity == 0.) {
      identity_free_candidate_faces++;
      continue;
    }
    if (!isfinite(fm.x[]) || fm.x[] < 0. ||
        !isfinite(dual_compound_alpha.x[]) ||
        dual_compound_alpha.x[] < 0. ||
        (fm.x[] > 0. && dual_compound_alpha.x[] <= 0.)) {
      invalid_property_faces++;
      continue;
    }

    const double physical_dual_volume = 2.*pi*fm.x[]*sq(Delta);
    for (int bin = 0; bin < bins; bin++) {
      const double lower_weight = lower_geometric_weight[bin]*lower_identity;
      const double upper_weight = upper_geometric_weight[bin]*upper_identity;
      candidate_lower_support[bin] += physical_dual_volume*lower_weight;
      candidate_upper_support[bin] += physical_dual_volume*upper_weight;
    }
  }

  ledger.invalid_identity_faces = invalid_identity_faces;
  ledger.invalid_property_faces = invalid_property_faces;
  ledger.identity_free_candidate_faces = identity_free_candidate_faces;
  ledger.maximum_identity_bound_defect = maximum_identity_bound_defect;
  ledger.identity_audit_failed = invalid_identity_faces > 0 ||
    maximum_identity_bound_defect >
      AXISYMMETRIC_FILM_STATE_AXI_IDENTITY_BOUND_TOLERANCE;
  if (ledger.identity_audit_failed || invalid_property_faces)
    return ledger;

  face vector candidate_owner = new face vector;
  foreach_face()
    candidate_owner.x[] = AXISYMMETRIC_FILM_STATE_AXI_OWNER_NONE;

  double overlap_candidate_volume = 0.;
  double potential_local_cancellation_force = 0.;
  long lower_owned_faces = 0, upper_owned_faces = 0;
  long overlap_candidate_faces = 0, identity_tie_faces = 0;
  long distance_tie_faces = 0;
  foreach_face(x, reduction(+:lower_support[:bins])
                  reduction(+:upper_support[:bins])
                  reduction(+:overlap_candidate_volume)
                  reduction(+:potential_local_cancellation_force)
                  reduction(+:lower_owned_faces)
                  reduction(+:upper_owned_faces)
                  reduction(+:overlap_candidate_faces)
                  reduction(+:identity_tie_faces)
                  reduction(+:distance_tie_faces)) {
    if (y < 0. || y >= support_half_width + Delta)
      continue;

    const double lower_identity = axisymmetric_film_state_axi_identity_fraction(
      .5*(lower_envelope_identity[] + lower_envelope_identity[-1]));
    const double upper_identity = axisymmetric_film_state_axi_identity_fraction(
      .5*(upper_envelope_identity[] + upper_envelope_identity[-1]));
    if (lower_identity == 0. && upper_identity == 0.)
      continue;

    bool lower_candidate = false, upper_candidate = false;
    double minimum_lower_distance = HUGE_VAL;
    double minimum_upper_distance = HUGE_VAL;
    double lower_weights[bins], upper_weights[bins];
    for (int bin = 0; bin < bins; bin++) {
      lower_weights[bin] = lower_identity*
        axisymmetric_film_state_axi_support_weight(state, bin,
          AXISYMMETRIC_FILM_STATE_AXI_OWNER_LOWER, x, y, Delta);
      upper_weights[bin] = upper_identity*
        axisymmetric_film_state_axi_support_weight(state, bin,
          AXISYMMETRIC_FILM_STATE_AXI_OWNER_UPPER, x, y, Delta);
      if (lower_weights[bin] > 0.) {
        lower_candidate = true;
        minimum_lower_distance = fmin(minimum_lower_distance,
          fabs(x - (state->midpoint[bin] - .5*state->thickness[bin])));
      }
      if (upper_weights[bin] > 0.) {
        upper_candidate = true;
        minimum_upper_distance = fmin(minimum_upper_distance,
          fabs(x - (state->midpoint[bin] + .5*state->thickness[bin])));
      }
    }
    if (!lower_candidate && !upper_candidate)
      continue;

    const double physical_dual_volume = 2.*pi*fm.x[]*sq(Delta);
    int owner = AXISYMMETRIC_FILM_STATE_AXI_OWNER_NONE;
    if (lower_candidate && !upper_candidate)
      owner = AXISYMMETRIC_FILM_STATE_AXI_OWNER_LOWER;
    else if (upper_candidate && !lower_candidate)
      owner = AXISYMMETRIC_FILM_STATE_AXI_OWNER_UPPER;
    else {
      overlap_candidate_faces++;
      overlap_candidate_volume += physical_dual_volume;
      const double identity_tolerance = 128.*DBL_EPSILON*max(1.,
        max(fabs(lower_identity), fabs(upper_identity)));
      if (lower_identity > upper_identity + identity_tolerance)
        owner = AXISYMMETRIC_FILM_STATE_AXI_OWNER_LOWER;
      else if (upper_identity > lower_identity + identity_tolerance)
        owner = AXISYMMETRIC_FILM_STATE_AXI_OWNER_UPPER;
      else {
        identity_tie_faces++;
        const double distance_tolerance = 128.*DBL_EPSILON*max(Delta,
          max(minimum_lower_distance, minimum_upper_distance));
        if (minimum_lower_distance <
            minimum_upper_distance - distance_tolerance)
          owner = AXISYMMETRIC_FILM_STATE_AXI_OWNER_LOWER;
        else if (minimum_upper_distance <
                 minimum_lower_distance - distance_tolerance)
          owner = AXISYMMETRIC_FILM_STATE_AXI_OWNER_UPPER;
        else {
          distance_tie_faces++;
          owner = AXISYMMETRIC_FILM_STATE_AXI_OWNER_LOWER;
        }
      }
    }
    candidate_owner.x[] = owner;
    if (owner == AXISYMMETRIC_FILM_STATE_AXI_OWNER_LOWER)
      lower_owned_faces++;
    else if (owner == AXISYMMETRIC_FILM_STATE_AXI_OWNER_UPPER)
      upper_owned_faces++;

    double candidate_lower_density = 0., candidate_upper_density = 0.;
    for (int bin = 0; bin < bins; bin++) {
      if (owner == AXISYMMETRIC_FILM_STATE_AXI_OWNER_LOWER)
        lower_support[bin] += physical_dual_volume*lower_weights[bin];
      else if (owner == AXISYMMETRIC_FILM_STATE_AXI_OWNER_UPPER)
        upper_support[bin] += physical_dual_volume*upper_weights[bin];
      if (matched_result->ring_force[bin] == 0.)
        continue;
      if (candidate_lower_support[bin] > 0.)
        candidate_lower_density -= matched_result->ring_force[bin]/
          candidate_lower_support[bin]*lower_weights[bin];
      if (candidate_upper_support[bin] > 0.)
        candidate_upper_density += matched_result->ring_force[bin]/
          candidate_upper_support[bin]*upper_weights[bin];
    }
    if (lower_candidate && upper_candidate &&
        candidate_lower_density < 0. && candidate_upper_density > 0.)
      potential_local_cancellation_force += physical_dual_volume*
        min(-candidate_lower_density, candidate_upper_density);
  }

  ledger.lower_owned_faces = lower_owned_faces;
  ledger.upper_owned_faces = upper_owned_faces;
  ledger.overlap_candidate_faces = overlap_candidate_faces;
  ledger.identity_tie_faces = identity_tie_faces;
  ledger.distance_tie_faces = distance_tie_faces;
  ledger.overlap_candidate_volume = overlap_candidate_volume;
  ledger.potential_local_cancellation_force =
    potential_local_cancellation_force;
  for (int bin = 0; bin < bins; bin++) {
    ledger.lower_support_volume += lower_support[bin];
    ledger.upper_support_volume += upper_support[bin];
    if (matched_result->ring_force[bin] != 0. && lower_support[bin] <= 0.)
      ledger.missing_lower_support_bins++;
    if (matched_result->ring_force[bin] != 0. && upper_support[bin] <= 0.)
      ledger.missing_upper_support_bins++;
  }
  if (ledger.missing_lower_support_bins ||
      ledger.missing_upper_support_bins) {
    delete ((scalar *){candidate_owner});
    return ledger;
  }

  face vector acceleration_before = new face vector;
  face vector film_acceleration_increment = new face vector;
  foreach_face() {
    acceleration_before.x[] = acceleration_field.x[];
    film_acceleration_increment.x[] = 0.;
  }

  foreach_face()
    owner_mask.x[] = candidate_owner.x[];

  double lower_by_bin[bins], upper_by_bin[bins];
  for (int bin = 0; bin < bins; bin++)
    lower_by_bin[bin] = upper_by_bin[bin] = 0.;
  foreach_face(x, reduction(+:lower_by_bin[:bins])
                  reduction(+:upper_by_bin[:bins])) {
    const int owner = (int) owner_mask.x[];
    if (owner == AXISYMMETRIC_FILM_STATE_AXI_OWNER_NONE)
      continue;
    const double physical_dual_volume = 2.*pi*fm.x[]*sq(Delta);
    double density_force = 0.;
    for (int bin = 0; bin < bins; bin++) {
      if (matched_result->ring_force[bin] == 0.)
        continue;
      if (owner == AXISYMMETRIC_FILM_STATE_AXI_OWNER_LOWER) {
        const double identity = axisymmetric_film_state_axi_identity_fraction(
          .5*(lower_envelope_identity[] + lower_envelope_identity[-1]));
        const double weight = identity*
          axisymmetric_film_state_axi_support_weight(
            state, bin, owner, x, y, Delta);
        const double bin_density = -matched_result->ring_force[bin]/
          lower_support[bin]*weight;
        density_force += bin_density;
        lower_by_bin[bin] += bin_density*physical_dual_volume;
      }
      else {
        const double identity = axisymmetric_film_state_axi_identity_fraction(
          .5*(upper_envelope_identity[] + upper_envelope_identity[-1]));
        const double weight = identity*
          axisymmetric_film_state_axi_support_weight(
            state, bin, owner, x, y, Delta);
        const double bin_density = matched_result->ring_force[bin]/
          upper_support[bin]*weight;
        density_force += bin_density;
        upper_by_bin[bin] += bin_density*physical_dual_volume;
      }
    }
    film_acceleration_increment.x[] =
      dual_compound_alpha.x[]/(fm.x[] + SEPS)*density_force;
    acceleration_field.x[] += film_acceleration_increment.x[];
  }

  double measured_lower = 0., measured_upper = 0.;
  double measured_lower_moment = 0., measured_upper_moment = 0.;
  double maximum_application_error = 0.;
  double application_force_error_l1 = 0.;
  long nonfinite_actual_faces = 0;
  long application_roundoff_faces = 0, swamped_increment_faces = 0;
  foreach_face(x, reduction(+:measured_lower) reduction(+:measured_upper)
                  reduction(+:measured_lower_moment)
                  reduction(+:measured_upper_moment)
                  reduction(max:maximum_application_error)
                  reduction(+:application_force_error_l1)
                  reduction(+:nonfinite_actual_faces)
                  reduction(+:application_roundoff_faces)
                  reduction(+:swamped_increment_faces)) {
    const int owner = (int) owner_mask.x[];
    if (owner == AXISYMMETRIC_FILM_STATE_AXI_OWNER_NONE)
      continue;
    const double delta_acceleration = film_acceleration_increment.x[];
    const double observed_increment = acceleration_field.x[] -
      acceleration_before.x[];
    const double application_error =
      fabs(observed_increment - delta_acceleration);
    maximum_application_error = max(maximum_application_error,
      application_error);
    if (application_error > 0.)
      application_roundoff_faces++;
    if (delta_acceleration != 0. && observed_increment == 0.)
      swamped_increment_faces++;
    const double physical_dual_volume = 2.*pi*fm.x[]*sq(Delta);
    if (physical_dual_volume == 0.)
      continue;
    const double acceleration_to_force =
      (fm.x[] + SEPS)/dual_compound_alpha.x[]*physical_dual_volume;
    const double force_density = observed_increment*(fm.x[] + SEPS)/
      dual_compound_alpha.x[];
    const double actual_force = observed_increment*acceleration_to_force;
    const double force_error_bound =
      application_error*fabs(acceleration_to_force);
    if (!isfinite(delta_acceleration) || !isfinite(observed_increment) ||
        !isfinite(acceleration_to_force) || !isfinite(force_density) ||
        !isfinite(actual_force) || !isfinite(force_error_bound)) {
      nonfinite_actual_faces++;
      continue;
    }
    application_force_error_l1 += force_error_bound;
    if (owner == AXISYMMETRIC_FILM_STATE_AXI_OWNER_LOWER) {
      measured_lower += actual_force;
      measured_lower_moment += y*actual_force;
    }
    else if (owner == AXISYMMETRIC_FILM_STATE_AXI_OWNER_UPPER) {
      measured_upper += actual_force;
      measured_upper_moment += y*actual_force;
    }
    else
      nonfinite_actual_faces++;
  }

  ledger.measured_lower_axial_force = measured_lower;
  ledger.measured_upper_axial_force = measured_upper;
  ledger.measured_net_axial_force = measured_lower + measured_upper;
  ledger.lower_force_residual = measured_lower -
    ledger.requested_lower_axial_force;
  ledger.upper_force_residual = measured_upper -
    ledger.requested_upper_axial_force;
  ledger.net_force_residual = ledger.measured_net_axial_force -
    ledger.requested_net_axial_force;
  ledger.maximum_application_error = maximum_application_error;
  ledger.application_force_error_l1 = application_force_error_l1;
  ledger.application_roundoff_faces = application_roundoff_faces;
  ledger.swamped_increment_faces = swamped_increment_faces;
  const double requested_force_scale = max(DBL_MIN,
    max(fabs(ledger.requested_lower_axial_force),
        fabs(ledger.requested_upper_axial_force)));
  ledger.force_audit_tolerance =
    AXISYMMETRIC_FILM_STATE_AXI_FORCE_AUDIT_RELATIVE_TOLERANCE*
    requested_force_scale;
  ledger.force_audit_nonfinite = nonfinite_actual_faces > 0 ||
    !isfinite(ledger.measured_lower_axial_force) ||
    !isfinite(ledger.measured_upper_axial_force) ||
    !isfinite(ledger.lower_force_residual) ||
    !isfinite(ledger.upper_force_residual);
  ledger.force_audit_failed = ledger.force_audit_nonfinite ||
    fabs(ledger.lower_force_residual) > ledger.force_audit_tolerance ||
    fabs(ledger.upper_force_residual) > ledger.force_audit_tolerance;
  if (ledger.force_audit_failed) {
    foreach_face() {
      acceleration_field.x[] = acceleration_before.x[];
      owner_mask.x[] = AXISYMMETRIC_FILM_STATE_AXI_OWNER_NONE;
    }
    delete ((scalar *){acceleration_before});
    delete ((scalar *){film_acceleration_increment});
    delete ((scalar *){candidate_owner});
    return ledger;
  }

  ledger.actual_lower_axial_force = measured_lower;
  ledger.actual_upper_axial_force = measured_upper;
  ledger.actual_net_axial_force = measured_lower + measured_upper;
  ledger.lower_radial_force_moment = measured_lower_moment;
  ledger.upper_radial_force_moment = measured_upper_moment;
  // A complete axisymmetric traction ring has zero vector torque. The scalar
  // r*Fz first moment is retained separately as a support-discretization audit.
  ledger.actual_net_axial_moment = 0.;
  ledger.radial_first_moment_imbalance =
    measured_lower_moment + measured_upper_moment;
  ledger.radial_first_moment_tolerance = max(
    ledger.force_audit_tolerance*support_half_width,
    AXISYMMETRIC_FILM_STATE_AXI_RADIAL_MOMENT_RELATIVE_TOLERANCE*
      requested_force_scale*support_half_width);
  for (int bin = 0; bin < bins; bin++)
    ledger.maximum_bin_force_imbalance = max(
      ledger.maximum_bin_force_imbalance,
      fabs(lower_by_bin[bin] + upper_by_bin[bin]));
  ledger.applied = true;
  delete ((scalar *){acceleration_before});
  delete ((scalar *){film_acceleration_increment});
  delete ((scalar *){candidate_owner});
  return ledger;
}

static inline ContactAmrGeometry axisymmetric_film_state_contact_geometry (
  const AxisymmetricFilmState * state)
{
  ContactAmrGeometry geometry = {0};
  if (!axisymmetric_film_state_axi_state_is_ready(state))
    return geometry;

  geometry.valid = true;
  geometry.virtual_bounds_present = true;
  geometry.active_state_support = true;
  geometry.minimum_gap = state->thickness[0];
  geometry.axial_minimum = state->midpoint[0] - .5*state->thickness[0];
  geometry.axial_maximum = state->midpoint[0] + .5*state->thickness[0];
  for (int bin = 1; bin < AXISYMMETRIC_FILM_STATE_BINS; bin++) {
    geometry.minimum_gap = fmin(geometry.minimum_gap,
      state->thickness[bin]);
    geometry.axial_minimum = fmin(geometry.axial_minimum,
      state->midpoint[bin] - .5*state->thickness[bin]);
    geometry.axial_maximum = fmax(geometry.axial_maximum,
      state->midpoint[bin] + .5*state->thickness[bin]);
  }
  return geometry;
}

#endif
