#ifndef AXISYMMETRIC_MUSEHANE_TRACTION_AXI_H
#define AXISYMMETRIC_MUSEHANE_TRACTION_AXI_H

#include <float.h>
#include <math.h>
#include <stdbool.h>

#include "axisymmetric_musehane_support_policy.h"
#include "film_model/musehane_traction_deposition.h"

#ifndef AXISYMMETRIC_MATCHED_FILM_BINS
# define AXISYMMETRIC_MATCHED_FILM_BINS 32
#endif

typedef enum {
  AXISYMMETRIC_MUSEHANE_TRACTION_OK = 0,
  AXISYMMETRIC_MUSEHANE_TRACTION_INVALID_INPUT = 1,
  AXISYMMETRIC_MUSEHANE_TRACTION_PLAN_FAILED = 2,
  AXISYMMETRIC_MUSEHANE_TRACTION_MISSING_SUPPORT = 3,
  AXISYMMETRIC_MUSEHANE_TRACTION_INVALID_FIELD = 4,
  AXISYMMETRIC_MUSEHANE_TRACTION_DELIVERY_FAILED = 5
} AxisymmetricMusehaneTractionStatus;

typedef struct {
  AxisymmetricMusehaneTractionStatus status;
  MusehaneTractionDepositionStatus plan_status;
  bool applied;
  int missing_support_components;
  int geometric_fallback_components;
  unsigned int geometric_fallback_mask_by_ring[AXISYMMETRIC_MATCHED_FILM_BINS];
  long invalid_identity_faces;
  long invalid_property_faces;
  long swamped_increment_faces;
  double requested_lower_x_force;
  double requested_upper_x_force;
  double requested_lower_r_force;
  double requested_upper_r_force;
  /* This reaction is carried by the resolved film and is intentionally not
     deposited as an extra force by this adapter. */
  double requested_resolved_reaction_force;
  double maximum_scalar_pair_closure_residual;
  double delivered_lower_x_force;
  double delivered_upper_x_force;
  double delivered_lower_r_force;
  double delivered_upper_r_force;
  double maximum_component_force_residual;
  double maximum_bin_force_residual;
  double requested_lower_power;
  double requested_upper_power;
  double delivered_lower_power;
  double delivered_upper_power;
  double requested_power;
  double delivered_power;
  double power_residual;
  double maximum_ring_power_residual;
  double maximum_interface_velocity_mismatch;
  double maximum_application_error;
  size_t ring_count;
  double requested_lower_x_force_by_ring[AXISYMMETRIC_MATCHED_FILM_BINS];
  double requested_upper_x_force_by_ring[AXISYMMETRIC_MATCHED_FILM_BINS];
  double requested_lower_r_force_by_ring[AXISYMMETRIC_MATCHED_FILM_BINS];
  double requested_upper_r_force_by_ring[AXISYMMETRIC_MATCHED_FILM_BINS];
  double delivered_lower_x_force_by_ring[AXISYMMETRIC_MATCHED_FILM_BINS];
  double delivered_upper_x_force_by_ring[AXISYMMETRIC_MATCHED_FILM_BINS];
  double delivered_lower_r_force_by_ring[AXISYMMETRIC_MATCHED_FILM_BINS];
  double delivered_upper_r_force_by_ring[AXISYMMETRIC_MATCHED_FILM_BINS];
} AxisymmetricMusehaneTractionLedger;

typedef enum {
  AXISYMMETRIC_MUSEHANE_VELOCITY_OK = 0,
  AXISYMMETRIC_MUSEHANE_VELOCITY_INVALID_INPUT = 1,
  AXISYMMETRIC_MUSEHANE_VELOCITY_INVALID_FIELD = 2,
  AXISYMMETRIC_MUSEHANE_VELOCITY_MISSING_SUPPORT = 3
} AxisymmetricMusehaneVelocityStatus;

typedef struct {
  AxisymmetricMusehaneVelocityStatus status;
  bool observed;
  int missing_support_components;
  int geometric_fallback_components;
  int first_geometric_fallback_bin;
  unsigned int first_geometric_fallback_mask;
  unsigned int geometric_fallback_mask_by_ring[AXISYMMETRIC_MATCHED_FILM_BINS];
  int first_missing_support_bin;
  unsigned int first_missing_support_mask;
  long invalid_identity_faces;
  long invalid_velocity_faces;
  double first_missing_lower_x_support;
  double first_missing_lower_r_support;
  double first_missing_upper_x_support;
  double first_missing_upper_r_support;
  double first_missing_lower_x_geometric_support;
  double first_missing_lower_r_geometric_support;
  double first_missing_upper_x_geometric_support;
  double first_missing_upper_r_geometric_support;
  double first_missing_lower_x_max_identity;
  double first_missing_lower_r_max_identity;
  double first_missing_upper_x_max_identity;
  double first_missing_upper_r_max_identity;
  double maximum_plic_velocity_adjustment;
} AxisymmetricMusehaneVelocityLedger;

static inline double axisymmetric_musehane_traction_kernel (
  double coordinate, double center, double delta)
{
  if (!isfinite(coordinate) || !isfinite(center) ||
      !isfinite(delta) || delta <= 0.)
    return 0.;
  const double distance = fabs(coordinate - center)/delta;
  return distance < 1. ? 1. - distance : 0.;
}

static inline bool axisymmetric_musehane_identity_is_valid (double value)
{
  return isfinite(value) && value >= -128.*DBL_EPSILON &&
    value <= 1. + 128.*DBL_EPSILON;
}

static inline double axisymmetric_musehane_identity (double value)
{
  return fmin(1., fmax(0., value));
}

static inline void axisymmetric_musehane_select_side_weights (
  double lower_identity, double upper_identity,
  double lower_geometric, double upper_geometric,
  double lower_distance, double upper_distance,
  double *lower_weight, double *upper_weight)
{
  *lower_weight = lower_identity*lower_geometric;
  *upper_weight = upper_identity*upper_geometric;
  if (*lower_weight > 0. && *upper_weight > 0.) {
    const double scale = fmax(DBL_MIN,
      fmax(*lower_weight, *upper_weight));
    const double tolerance = 256.*DBL_EPSILON*scale;
    if (*lower_weight > *upper_weight + tolerance)
      *upper_weight = 0.;
    else if (*upper_weight > *lower_weight + tolerance)
      *lower_weight = 0.;
    else if (lower_distance <= upper_distance)
      *upper_weight = 0.;
    else
      *lower_weight = 0.;
  }
}

/* Observe interface tangential velocities with the exact transpose of the
   force-spreading operator used below.  This gives the discrete coupling the
   power identity F.(J u) == (J^T F).u without force rescaling. */
static inline AxisymmetricMusehaneVelocityLedger
axisymmetric_musehane_observe_support_velocity_axi (
  vector velocity,
  scalar lower_envelope_identity,
  scalar upper_envelope_identity,
  FilmRingObservation *rings,
  size_t ring_count,
  double support_half_width)
{
  enum { bins = AXISYMMETRIC_MATCHED_FILM_BINS };
  AxisymmetricMusehaneVelocityLedger ledger = {
    .status = AXISYMMETRIC_MUSEHANE_VELOCITY_INVALID_INPUT,
    .first_missing_support_bin = -1,
    .first_geometric_fallback_bin = -1
  };
  if (!rings || ring_count != bins || !isfinite(support_half_width) ||
      support_half_width <= 0.)
    return ledger;

  double support_lx[bins], support_ux[bins];
  double support_lr[bins], support_ur[bins];
  double geometric_lx[bins], geometric_ux[bins];
  double geometric_lr[bins], geometric_ur[bins];
  double identity_lx[bins], identity_ux[bins];
  double identity_lr[bins], identity_ur[bins];
  double velocity_lx[bins], velocity_ux[bins];
  double velocity_lr[bins], velocity_ur[bins];
  double geometric_velocity_lx[bins], geometric_velocity_ux[bins];
  double geometric_velocity_lr[bins], geometric_velocity_ur[bins];
  for (int bin = 0; bin < bins; bin++) {
    const FilmRingObservation ring = rings[bin];
    const double lower_norm = hypot(ring.lower_normal_x,
      ring.lower_normal_r);
    const double upper_norm = hypot(ring.upper_normal_x,
      ring.upper_normal_r);
    if (!ring.geometry_valid || !ring.lower_identity_valid ||
        !ring.upper_identity_valid || !isfinite(ring.radius) ||
        !isfinite(ring.lower_position) ||
        !isfinite(ring.upper_position) ||
        !isfinite(lower_norm) || lower_norm <= DBL_MIN ||
        !isfinite(upper_norm) || upper_norm <= DBL_MIN)
      return ledger;
    support_lx[bin] = support_ux[bin] =
      support_lr[bin] = support_ur[bin] =
      geometric_lx[bin] = geometric_ux[bin] =
      geometric_lr[bin] = geometric_ur[bin] =
      identity_lx[bin] = identity_ux[bin] =
      identity_lr[bin] = identity_ur[bin] =
      velocity_lx[bin] = velocity_ux[bin] =
      velocity_lr[bin] = velocity_ur[bin] =
      geometric_velocity_lx[bin] = geometric_velocity_ux[bin] =
      geometric_velocity_lr[bin] = geometric_velocity_ur[bin] = 0.;
  }

  long invalid_identity = 0, invalid_velocity = 0;
  foreach_face(x, reduction(+:support_lx[:bins])
                  reduction(+:support_ux[:bins])
                  reduction(+:geometric_lx[:bins])
                  reduction(+:geometric_ux[:bins])
                  reduction(max:identity_lx[:bins])
                  reduction(max:identity_ux[:bins])
                  reduction(+:velocity_lx[:bins])
                  reduction(+:velocity_ux[:bins])
                  reduction(+:geometric_velocity_lx[:bins])
                  reduction(+:geometric_velocity_ux[:bins])
                  reduction(+:invalid_identity)
                  reduction(+:invalid_velocity)) {
    if (y < 0. || y >= support_half_width + Delta)
      continue;
    const double lower_raw = .5*(lower_envelope_identity[] +
      lower_envelope_identity[-1]);
    const double upper_raw = .5*(upper_envelope_identity[] +
      upper_envelope_identity[-1]);
    const double face_velocity = .5*(velocity.x[] + velocity.x[-1]);
    if (!axisymmetric_musehane_identity_is_valid(lower_raw) ||
        !axisymmetric_musehane_identity_is_valid(upper_raw)) {
      invalid_identity++;
      continue;
    }
    if (!isfinite(fm.x[]) || fm.x[] < 0. ||
        !isfinite(face_velocity)) {
      invalid_velocity++;
      continue;
    }
    const double volume = 2.*pi*fm.x[]*sq(Delta);
    const double lower_identity = axisymmetric_musehane_identity(lower_raw);
    const double upper_identity = axisymmetric_musehane_identity(upper_raw);
    for (int bin = 0; bin < bins; bin++) {
      const FilmRingObservation ring = rings[bin];
      const double radial = axisymmetric_musehane_traction_kernel(
        y, ring.radius, Delta);
      const double lower_geometric = radial*
        axisymmetric_musehane_traction_kernel(
          x, ring.lower_position, Delta);
      const double upper_geometric = radial*
        axisymmetric_musehane_traction_kernel(
          x, ring.upper_position, Delta);
      geometric_lx[bin] += volume*lower_geometric;
      geometric_ux[bin] += volume*upper_geometric;
      geometric_velocity_lx[bin] += volume*lower_geometric*face_velocity;
      geometric_velocity_ux[bin] += volume*upper_geometric*face_velocity;
      if (lower_geometric > 0.)
        identity_lx[bin] = max(identity_lx[bin], lower_identity);
      if (upper_geometric > 0.)
        identity_ux[bin] = max(identity_ux[bin], upper_identity);
      double lw = 0., uw = 0.;
      axisymmetric_musehane_select_side_weights(
        lower_identity, upper_identity,
        lower_geometric, upper_geometric,
        fabs(x - ring.lower_position), fabs(x - ring.upper_position),
        &lw, &uw);
      support_lx[bin] += volume*lw;
      support_ux[bin] += volume*uw;
      velocity_lx[bin] += volume*lw*face_velocity;
      velocity_ux[bin] += volume*uw*face_velocity;
    }
  }
  foreach_face(y, reduction(+:support_lr[:bins])
                  reduction(+:support_ur[:bins])
                  reduction(+:geometric_lr[:bins])
                  reduction(+:geometric_ur[:bins])
                  reduction(max:identity_lr[:bins])
                  reduction(max:identity_ur[:bins])
                  reduction(+:velocity_lr[:bins])
                  reduction(+:velocity_ur[:bins])
                  reduction(+:geometric_velocity_lr[:bins])
                  reduction(+:geometric_velocity_ur[:bins])
                  reduction(+:invalid_identity)
                  reduction(+:invalid_velocity)) {
    if (y < 0. || y >= support_half_width + Delta)
      continue;
    const double lower_raw = .5*(lower_envelope_identity[] +
      lower_envelope_identity[0,-1]);
    const double upper_raw = .5*(upper_envelope_identity[] +
      upper_envelope_identity[0,-1]);
    const double face_velocity = .5*(velocity.y[] + velocity.y[0,-1]);
    if (!axisymmetric_musehane_identity_is_valid(lower_raw) ||
        !axisymmetric_musehane_identity_is_valid(upper_raw)) {
      invalid_identity++;
      continue;
    }
    if (!isfinite(fm.y[]) || fm.y[] < 0. ||
        !isfinite(face_velocity)) {
      invalid_velocity++;
      continue;
    }
    const double volume = 2.*pi*fm.y[]*sq(Delta);
    const double lower_identity = axisymmetric_musehane_identity(lower_raw);
    const double upper_identity = axisymmetric_musehane_identity(upper_raw);
    for (int bin = 0; bin < bins; bin++) {
      const FilmRingObservation ring = rings[bin];
      const double radial = axisymmetric_musehane_traction_kernel(
        y, ring.radius, Delta);
      const double lower_geometric = radial*
        axisymmetric_musehane_traction_kernel(
          x, ring.lower_position, Delta);
      const double upper_geometric = radial*
        axisymmetric_musehane_traction_kernel(
          x, ring.upper_position, Delta);
      geometric_lr[bin] += volume*lower_geometric;
      geometric_ur[bin] += volume*upper_geometric;
      geometric_velocity_lr[bin] += volume*lower_geometric*face_velocity;
      geometric_velocity_ur[bin] += volume*upper_geometric*face_velocity;
      if (lower_geometric > 0.)
        identity_lr[bin] = max(identity_lr[bin], lower_identity);
      if (upper_geometric > 0.)
        identity_ur[bin] = max(identity_ur[bin], upper_identity);
      double lw = 0., uw = 0.;
      axisymmetric_musehane_select_side_weights(
        lower_identity, upper_identity,
        lower_geometric, upper_geometric,
        fabs(x - ring.lower_position), fabs(x - ring.upper_position),
        &lw, &uw);
      support_lr[bin] += volume*lw;
      support_ur[bin] += volume*uw;
      velocity_lr[bin] += volume*lw*face_velocity;
      velocity_ur[bin] += volume*uw*face_velocity;
    }
  }
  ledger.invalid_identity_faces = invalid_identity;
  ledger.invalid_velocity_faces = invalid_velocity;
  if (invalid_identity || invalid_velocity) {
    ledger.status = AXISYMMETRIC_MUSEHANE_VELOCITY_INVALID_FIELD;
    return ledger;
  }

  double lower_velocity[bins], upper_velocity[bins];
  for (int bin = 0; bin < bins; bin++) {
    const FilmRingObservation ring = rings[bin];
    const double lower_norm = hypot(ring.lower_normal_x,
      ring.lower_normal_r);
    const double upper_norm = hypot(ring.upper_normal_x,
      ring.upper_normal_r);
    const double lower_tangent_x = -ring.lower_normal_r/lower_norm;
    const double lower_tangent_r = ring.lower_normal_x/lower_norm;
    const double upper_tangent_x = -ring.upper_normal_r/upper_norm;
    const double upper_tangent_r = ring.upper_normal_x/upper_norm;
    const AxisymmetricMusehaneSupportSelection lower_x_selection =
      axisymmetric_musehane_select_component_support(
        support_lx[bin], geometric_lx[bin], lower_tangent_x != 0.);
    const AxisymmetricMusehaneSupportSelection lower_r_selection =
      axisymmetric_musehane_select_component_support(
        support_lr[bin], geometric_lr[bin], lower_tangent_r != 0.);
    const AxisymmetricMusehaneSupportSelection upper_x_selection =
      axisymmetric_musehane_select_component_support(
        support_ux[bin], geometric_ux[bin], upper_tangent_x != 0.);
    const AxisymmetricMusehaneSupportSelection upper_r_selection =
      axisymmetric_musehane_select_component_support(
        support_ur[bin], geometric_ur[bin], upper_tangent_r != 0.);
    unsigned int fallback_mask = 0U;
    if (lower_x_selection.geometric_fallback) {
      fallback_mask |= 1U;
      support_lx[bin] = lower_x_selection.support;
      velocity_lx[bin] = geometric_velocity_lx[bin];
    }
    if (lower_r_selection.geometric_fallback) {
      fallback_mask |= 2U;
      support_lr[bin] = lower_r_selection.support;
      velocity_lr[bin] = geometric_velocity_lr[bin];
    }
    if (upper_x_selection.geometric_fallback) {
      fallback_mask |= 4U;
      support_ux[bin] = upper_x_selection.support;
      velocity_ux[bin] = geometric_velocity_ux[bin];
    }
    if (upper_r_selection.geometric_fallback) {
      fallback_mask |= 8U;
      support_ur[bin] = upper_r_selection.support;
      velocity_ur[bin] = geometric_velocity_ur[bin];
    }
    ledger.geometric_fallback_mask_by_ring[bin] = fallback_mask;
    if (fallback_mask) {
      ledger.geometric_fallback_components +=
        !!(fallback_mask & 1U) + !!(fallback_mask & 2U) +
        !!(fallback_mask & 4U) + !!(fallback_mask & 8U);
      if (ledger.first_geometric_fallback_bin < 0) {
        ledger.first_geometric_fallback_bin = bin;
        ledger.first_geometric_fallback_mask = fallback_mask;
      }
    }
    unsigned int missing_mask = 0U;
    if (!lower_x_selection.valid)
      missing_mask |= 1U;
    if (!lower_r_selection.valid)
      missing_mask |= 2U;
    if (!upper_x_selection.valid)
      missing_mask |= 4U;
    if (!upper_r_selection.valid)
      missing_mask |= 8U;
    if (missing_mask) {
      ledger.missing_support_components +=
        !!(missing_mask & 1U) + !!(missing_mask & 2U) +
        !!(missing_mask & 4U) + !!(missing_mask & 8U);
      if (ledger.first_missing_support_bin < 0) {
        ledger.first_missing_support_bin = bin;
        ledger.first_missing_support_mask = missing_mask;
        ledger.first_missing_lower_x_support = support_lx[bin];
        ledger.first_missing_lower_r_support = support_lr[bin];
        ledger.first_missing_upper_x_support = support_ux[bin];
        ledger.first_missing_upper_r_support = support_ur[bin];
        ledger.first_missing_lower_x_geometric_support = geometric_lx[bin];
        ledger.first_missing_lower_r_geometric_support = geometric_lr[bin];
        ledger.first_missing_upper_x_geometric_support = geometric_ux[bin];
        ledger.first_missing_upper_r_geometric_support = geometric_ur[bin];
        ledger.first_missing_lower_x_max_identity = identity_lx[bin];
        ledger.first_missing_lower_r_max_identity = identity_lr[bin];
        ledger.first_missing_upper_x_max_identity = identity_ux[bin];
        ledger.first_missing_upper_r_max_identity = identity_ur[bin];
      }
      continue;
    }
    lower_velocity[bin] =
      (lower_tangent_x == 0. ? 0. :
        lower_tangent_x*velocity_lx[bin]/support_lx[bin]) +
      (lower_tangent_r == 0. ? 0. :
        lower_tangent_r*velocity_lr[bin]/support_lr[bin]);
    upper_velocity[bin] =
      (upper_tangent_x == 0. ? 0. :
        upper_tangent_x*velocity_ux[bin]/support_ux[bin]) +
      (upper_tangent_r == 0. ? 0. :
        upper_tangent_r*velocity_ur[bin]/support_ur[bin]);
    if (!isfinite(lower_velocity[bin]) ||
        !isfinite(upper_velocity[bin])) {
      ledger.invalid_velocity_faces++;
      continue;
    }
    ledger.maximum_plic_velocity_adjustment = max(
      ledger.maximum_plic_velocity_adjustment,
      max(fabs(lower_velocity[bin] - ring.lower_tangential_velocity),
          fabs(upper_velocity[bin] - ring.upper_tangential_velocity)));
  }
  if (ledger.missing_support_components) {
    ledger.status = AXISYMMETRIC_MUSEHANE_VELOCITY_MISSING_SUPPORT;
    return ledger;
  }
  if (ledger.invalid_velocity_faces) {
    ledger.status = AXISYMMETRIC_MUSEHANE_VELOCITY_INVALID_FIELD;
    return ledger;
  }
  for (int bin = 0; bin < bins; bin++) {
    rings[bin].lower_tangential_velocity = lower_velocity[bin];
    rings[bin].upper_tangential_velocity = upper_velocity[bin];
  }
  ledger.status = AXISYMMETRIC_MUSEHANE_VELOCITY_OK;
  ledger.observed = true;
  return ledger;
}

static inline AxisymmetricMusehaneTractionLedger
axisymmetric_musehane_apply_traction_axi (
  face vector acceleration_field,
  face vector inverse_density,
  vector velocity,
  scalar lower_envelope_identity,
  scalar upper_envelope_identity,
  const FilmModelResult *result,
  const FilmDNSView *dns,
  double support_half_width)
{
  enum { bins = AXISYMMETRIC_MATCHED_FILM_BINS };
  AxisymmetricMusehaneTractionLedger ledger = {
    .status = AXISYMMETRIC_MUSEHANE_TRACTION_INVALID_INPUT,
    .plan_status = MUSEHANE_TRACTION_DEPOSITION_INVALID_INPUT
  };
  if (!result || !dns || result->status != FILM_MODEL_OK ||
      !result->commit_allowed || !result->tractions ||
      result->ring_count != bins || dns->ring_count != bins ||
      result->ring_count != dns->ring_count ||
      !isfinite(support_half_width) || support_half_width <= 0.)
    return ledger;

  double lower_x[bins], lower_r[bins], upper_x[bins], upper_r[bins];
  MusehaneTractionDepositionOutput plan = {
    .lower_force_x = lower_x, .lower_force_r = lower_r,
    .upper_force_x = upper_x, .upper_force_r = upper_r
  };
  const MusehaneTractionDepositionInput plan_input = {
    .ring_count = result->ring_count,
    .rings = dns->rings,
    .tractions = result->tractions
  };
  ledger.plan_status = musehane_plan_traction_deposition(&plan_input, &plan);
  if (ledger.plan_status != MUSEHANE_TRACTION_DEPOSITION_OK) {
    ledger.status = AXISYMMETRIC_MUSEHANE_TRACTION_PLAN_FAILED;
    return ledger;
  }
  ledger.requested_power = plan.requested_total_power;
  ledger.ring_count = result->ring_count;
  ledger.requested_lower_power = plan.requested_lower_power;
  ledger.requested_upper_power = plan.requested_upper_power;
  ledger.requested_resolved_reaction_force =
    plan.requested_resolved_reaction_force;
  ledger.maximum_scalar_pair_closure_residual =
    plan.maximum_scalar_pair_residual;
  for (int bin = 0; bin < bins; bin++) {
    ledger.requested_lower_x_force_by_ring[bin] = lower_x[bin];
    ledger.requested_upper_x_force_by_ring[bin] = upper_x[bin];
    ledger.requested_lower_r_force_by_ring[bin] = lower_r[bin];
    ledger.requested_upper_r_force_by_ring[bin] = upper_r[bin];
    ledger.requested_lower_x_force += lower_x[bin];
    ledger.requested_upper_x_force += upper_x[bin];
    ledger.requested_lower_r_force += lower_r[bin];
    ledger.requested_upper_r_force += upper_r[bin];
  }

  double support_lx[bins], support_ux[bins];
  double support_lr[bins], support_ur[bins];
  double geometric_lx[bins], geometric_ux[bins];
  double geometric_lr[bins], geometric_ur[bins];
  for (int bin = 0; bin < bins; bin++)
    support_lx[bin] = support_ux[bin] =
      support_lr[bin] = support_ur[bin] =
      geometric_lx[bin] = geometric_ux[bin] =
      geometric_lr[bin] = geometric_ur[bin] = 0.;

  long invalid_identity = 0, invalid_property = 0;
  foreach_face(x, reduction(+:support_lx[:bins])
                  reduction(+:support_ux[:bins])
                  reduction(+:geometric_lx[:bins])
                  reduction(+:geometric_ux[:bins])
                  reduction(+:invalid_identity)
                  reduction(+:invalid_property)) {
    if (y < 0. || y >= support_half_width + Delta)
      continue;
    const double lower_raw = .5*(lower_envelope_identity[] +
      lower_envelope_identity[-1]);
    const double upper_raw = .5*(upper_envelope_identity[] +
      upper_envelope_identity[-1]);
    if (!axisymmetric_musehane_identity_is_valid(lower_raw) ||
        !axisymmetric_musehane_identity_is_valid(upper_raw)) {
      invalid_identity++;
      continue;
    }
    if (!isfinite(fm.x[]) || fm.x[] < 0. ||
        !isfinite(inverse_density.x[]) || inverse_density.x[] < 0. ||
        (fm.x[] > 0. && inverse_density.x[] <= 0.)) {
      invalid_property++;
      continue;
    }
    const double volume = 2.*pi*fm.x[]*sq(Delta);
    const double lower_identity = axisymmetric_musehane_identity(lower_raw);
    const double upper_identity = axisymmetric_musehane_identity(upper_raw);
    for (int bin = 0; bin < bins; bin++) {
      const FilmRingObservation ring = dns->rings[bin];
      const double radial = axisymmetric_musehane_traction_kernel(
        y, ring.radius, Delta);
      const double lower_geometric = radial*
        axisymmetric_musehane_traction_kernel(
          x, ring.lower_position, Delta);
      const double upper_geometric = radial*
        axisymmetric_musehane_traction_kernel(
          x, ring.upper_position, Delta);
      double lw = 0., uw = 0.;
      axisymmetric_musehane_select_side_weights(
        lower_identity, upper_identity,
        lower_geometric, upper_geometric,
        fabs(x - ring.lower_position), fabs(x - ring.upper_position),
        &lw, &uw);
      support_lx[bin] += volume*lw;
      support_ux[bin] += volume*uw;
      geometric_lx[bin] += volume*lower_geometric;
      geometric_ux[bin] += volume*upper_geometric;
    }
  }
  foreach_face(y, reduction(+:support_lr[:bins])
                  reduction(+:support_ur[:bins])
                  reduction(+:geometric_lr[:bins])
                  reduction(+:geometric_ur[:bins])
                  reduction(+:invalid_identity)
                  reduction(+:invalid_property)) {
    if (y < 0. || y >= support_half_width + Delta)
      continue;
    const double lower_raw = .5*(lower_envelope_identity[] +
      lower_envelope_identity[0,-1]);
    const double upper_raw = .5*(upper_envelope_identity[] +
      upper_envelope_identity[0,-1]);
    if (!axisymmetric_musehane_identity_is_valid(lower_raw) ||
        !axisymmetric_musehane_identity_is_valid(upper_raw)) {
      invalid_identity++;
      continue;
    }
    if (!isfinite(fm.y[]) || fm.y[] < 0. ||
        !isfinite(inverse_density.y[]) || inverse_density.y[] < 0. ||
        (fm.y[] > 0. && inverse_density.y[] <= 0.)) {
      invalid_property++;
      continue;
    }
    const double volume = 2.*pi*fm.y[]*sq(Delta);
    const double lower_identity = axisymmetric_musehane_identity(lower_raw);
    const double upper_identity = axisymmetric_musehane_identity(upper_raw);
    for (int bin = 0; bin < bins; bin++) {
      const FilmRingObservation ring = dns->rings[bin];
      const double radial = axisymmetric_musehane_traction_kernel(
        y, ring.radius, Delta);
      const double lower_geometric = radial*
        axisymmetric_musehane_traction_kernel(
          x, ring.lower_position, Delta);
      const double upper_geometric = radial*
        axisymmetric_musehane_traction_kernel(
          x, ring.upper_position, Delta);
      double lw = 0., uw = 0.;
      axisymmetric_musehane_select_side_weights(
        lower_identity, upper_identity,
        lower_geometric, upper_geometric,
        fabs(x - ring.lower_position), fabs(x - ring.upper_position),
        &lw, &uw);
      support_lr[bin] += volume*lw;
      support_ur[bin] += volume*uw;
      geometric_lr[bin] += volume*lower_geometric;
      geometric_ur[bin] += volume*upper_geometric;
    }
  }
  ledger.invalid_identity_faces = invalid_identity;
  ledger.invalid_property_faces = invalid_property;
  if (invalid_identity || invalid_property) {
    ledger.status = AXISYMMETRIC_MUSEHANE_TRACTION_INVALID_FIELD;
    return ledger;
  }
  for (int bin = 0; bin < bins; bin++) {
    const AxisymmetricMusehaneSupportSelection selections[4] = {
      axisymmetric_musehane_select_component_support(
        support_lx[bin], geometric_lx[bin], lower_x[bin] != 0.),
      axisymmetric_musehane_select_component_support(
        support_lr[bin], geometric_lr[bin], lower_r[bin] != 0.),
      axisymmetric_musehane_select_component_support(
        support_ux[bin], geometric_ux[bin], upper_x[bin] != 0.),
      axisymmetric_musehane_select_component_support(
        support_ur[bin], geometric_ur[bin], upper_r[bin] != 0.)
    };
    double *supports[4] = {
      &support_lx[bin], &support_lr[bin], &support_ux[bin], &support_ur[bin]
    };
    unsigned int fallback_mask = 0U;
    for (int component = 0; component < 4; component++) {
      if (!selections[component].valid)
        ledger.missing_support_components++;
      else if (selections[component].geometric_fallback) {
        fallback_mask |= 1U << component;
        *supports[component] = selections[component].support;
        ledger.geometric_fallback_components++;
      }
    }
    ledger.geometric_fallback_mask_by_ring[bin] = fallback_mask;
  }
  if (ledger.missing_support_components) {
    ledger.status = AXISYMMETRIC_MUSEHANE_TRACTION_MISSING_SUPPORT;
    return ledger;
  }

  face vector before = new face vector;
  face vector lower_increment = new face vector;
  face vector upper_increment = new face vector;
  foreach_face() {
    before.x[] = acceleration_field.x[];
    lower_increment.x[] = upper_increment.x[] = 0.;
  }
  double delivered_lx[bins], delivered_ux[bins];
  double delivered_lr[bins], delivered_ur[bins];
  double delivered_lower_power_by_bin[bins];
  double delivered_upper_power_by_bin[bins];
  for (int bin = 0; bin < bins; bin++)
    delivered_lx[bin] = delivered_ux[bin] =
      delivered_lr[bin] = delivered_ur[bin] =
      delivered_lower_power_by_bin[bin] =
      delivered_upper_power_by_bin[bin] = 0.;

  foreach_face(x, reduction(+:delivered_lx[:bins])
                  reduction(+:delivered_ux[:bins])
                  reduction(+:delivered_lower_power_by_bin[:bins])
                  reduction(+:delivered_upper_power_by_bin[:bins])) {
    if (y < 0. || y >= support_half_width + Delta)
      continue;
    const double volume = 2.*pi*fm.x[]*sq(Delta);
    const double face_velocity = .5*(velocity.x[] + velocity.x[-1]);
    const double lower_identity = axisymmetric_musehane_identity(
      .5*(lower_envelope_identity[] + lower_envelope_identity[-1]));
    const double upper_identity = axisymmetric_musehane_identity(
      .5*(upper_envelope_identity[] + upper_envelope_identity[-1]));
    double lower_density = 0., upper_density = 0.;
    for (int bin = 0; bin < bins; bin++) {
      const FilmRingObservation ring = dns->rings[bin];
      const double radial = axisymmetric_musehane_traction_kernel(
        y, ring.radius, Delta);
      double lw = 0., uw = 0.;
      axisymmetric_musehane_select_side_weights(
        lower_identity, upper_identity,
        radial*axisymmetric_musehane_traction_kernel(
          x, ring.lower_position, Delta),
        radial*axisymmetric_musehane_traction_kernel(
          x, ring.upper_position, Delta),
        fabs(x - ring.lower_position), fabs(x - ring.upper_position),
        &lw, &uw);
      if (ledger.geometric_fallback_mask_by_ring[bin] & 1U)
        lw = radial*axisymmetric_musehane_traction_kernel(
          x, ring.lower_position, Delta);
      if (ledger.geometric_fallback_mask_by_ring[bin] & 4U)
        uw = radial*axisymmetric_musehane_traction_kernel(
          x, ring.upper_position, Delta);
      const double lf = lower_x[bin] == 0. ? 0. :
        lower_x[bin]/support_lx[bin]*lw;
      const double upper_force_density = upper_x[bin] == 0. ? 0. :
        upper_x[bin]/support_ux[bin]*uw;
      lower_density += lf;
      upper_density += upper_force_density;
      delivered_lx[bin] += lf*volume;
      delivered_ux[bin] += upper_force_density*volume;
      delivered_lower_power_by_bin[bin] += lf*volume*face_velocity;
      delivered_upper_power_by_bin[bin] +=
        upper_force_density*volume*face_velocity;
    }
    lower_increment.x[] = inverse_density.x[]/(fm.x[] + SEPS)*lower_density;
    upper_increment.x[] = inverse_density.x[]/(fm.x[] + SEPS)*upper_density;
  }
  foreach_face(y, reduction(+:delivered_lr[:bins])
                  reduction(+:delivered_ur[:bins])
                  reduction(+:delivered_lower_power_by_bin[:bins])
                  reduction(+:delivered_upper_power_by_bin[:bins])) {
    if (y < 0. || y >= support_half_width + Delta)
      continue;
    const double volume = 2.*pi*fm.y[]*sq(Delta);
    const double face_velocity = .5*(velocity.y[] + velocity.y[0,-1]);
    const double lower_identity = axisymmetric_musehane_identity(
      .5*(lower_envelope_identity[] + lower_envelope_identity[0,-1]));
    const double upper_identity = axisymmetric_musehane_identity(
      .5*(upper_envelope_identity[] + upper_envelope_identity[0,-1]));
    double lower_density = 0., upper_density = 0.;
    for (int bin = 0; bin < bins; bin++) {
      const FilmRingObservation ring = dns->rings[bin];
      const double radial = axisymmetric_musehane_traction_kernel(
        y, ring.radius, Delta);
      double lw = 0., uw = 0.;
      axisymmetric_musehane_select_side_weights(
        lower_identity, upper_identity,
        radial*axisymmetric_musehane_traction_kernel(
          x, ring.lower_position, Delta),
        radial*axisymmetric_musehane_traction_kernel(
          x, ring.upper_position, Delta),
        fabs(x - ring.lower_position), fabs(x - ring.upper_position),
        &lw, &uw);
      if (ledger.geometric_fallback_mask_by_ring[bin] & 2U)
        lw = radial*axisymmetric_musehane_traction_kernel(
          x, ring.lower_position, Delta);
      if (ledger.geometric_fallback_mask_by_ring[bin] & 8U)
        uw = radial*axisymmetric_musehane_traction_kernel(
          x, ring.upper_position, Delta);
      const double lf = lower_r[bin] == 0. ? 0. :
        lower_r[bin]/support_lr[bin]*lw;
      const double upper_force_density = upper_r[bin] == 0. ? 0. :
        upper_r[bin]/support_ur[bin]*uw;
      lower_density += lf;
      upper_density += upper_force_density;
      delivered_lr[bin] += lf*volume;
      delivered_ur[bin] += upper_force_density*volume;
      delivered_lower_power_by_bin[bin] += lf*volume*face_velocity;
      delivered_upper_power_by_bin[bin] +=
        upper_force_density*volume*face_velocity;
    }
    lower_increment.y[] = inverse_density.y[]/(fm.y[] + SEPS)*lower_density;
    upper_increment.y[] = inverse_density.y[]/(fm.y[] + SEPS)*upper_density;
  }

  foreach_face()
    acceleration_field.x[] += lower_increment.x[] + upper_increment.x[];
  long swamped = 0;
  double maximum_error = 0.;
  foreach_face(reduction(+:swamped) reduction(max:maximum_error)) {
    const double requested_increment = lower_increment.x[] +
      upper_increment.x[];
    const double observed_increment = acceleration_field.x[] - before.x[];
    const double error = fabs(observed_increment - requested_increment);
    maximum_error = max(maximum_error, error);
    if (requested_increment != 0. && observed_increment == 0.)
      swamped++;
    if (!isfinite(observed_increment) || !isfinite(error))
      swamped++;
  }
  ledger.swamped_increment_faces = swamped;
  ledger.maximum_application_error = maximum_error;
  for (int bin = 0; bin < bins; bin++) {
    ledger.delivered_lower_x_force_by_ring[bin] = delivered_lx[bin];
    ledger.delivered_upper_x_force_by_ring[bin] = delivered_ux[bin];
    ledger.delivered_lower_r_force_by_ring[bin] = delivered_lr[bin];
    ledger.delivered_upper_r_force_by_ring[bin] = delivered_ur[bin];
    ledger.delivered_lower_x_force += delivered_lx[bin];
    ledger.delivered_upper_x_force += delivered_ux[bin];
    ledger.delivered_lower_r_force += delivered_lr[bin];
    ledger.delivered_upper_r_force += delivered_ur[bin];
    ledger.maximum_bin_force_residual = max(
      ledger.maximum_bin_force_residual,
      max(max(fabs(delivered_lx[bin] - lower_x[bin]),
              fabs(delivered_ux[bin] - upper_x[bin])),
          max(fabs(delivered_lr[bin] - lower_r[bin]),
              fabs(delivered_ur[bin] - upper_r[bin]))));
    ledger.delivered_lower_power += delivered_lower_power_by_bin[bin];
    ledger.delivered_upper_power += delivered_upper_power_by_bin[bin];
    const double requested_lower_ring_power =
      result->tractions[bin].lower_tangential_force*
      dns->rings[bin].lower_tangential_velocity;
    const double requested_upper_ring_power =
      result->tractions[bin].upper_tangential_force*
      dns->rings[bin].upper_tangential_velocity;
    ledger.maximum_ring_power_residual = max(
      ledger.maximum_ring_power_residual,
      max(fabs(delivered_lower_power_by_bin[bin] -
               requested_lower_ring_power),
          fabs(delivered_upper_power_by_bin[bin] -
               requested_upper_ring_power)));
    if (result->tractions[bin].lower_tangential_force != 0.)
      ledger.maximum_interface_velocity_mismatch = max(
        ledger.maximum_interface_velocity_mismatch,
        fabs(delivered_lower_power_by_bin[bin]/
             result->tractions[bin].lower_tangential_force -
             dns->rings[bin].lower_tangential_velocity));
    if (result->tractions[bin].upper_tangential_force != 0.)
      ledger.maximum_interface_velocity_mismatch = max(
        ledger.maximum_interface_velocity_mismatch,
        fabs(delivered_upper_power_by_bin[bin]/
             result->tractions[bin].upper_tangential_force -
             dns->rings[bin].upper_tangential_velocity));
  }
  ledger.delivered_power = ledger.delivered_lower_power +
    ledger.delivered_upper_power;
  ledger.power_residual = ledger.delivered_power - ledger.requested_power;
  ledger.maximum_component_force_residual = max(
    max(fabs(ledger.delivered_lower_x_force -
             ledger.requested_lower_x_force),
        fabs(ledger.delivered_upper_x_force -
             ledger.requested_upper_x_force)),
    max(fabs(ledger.delivered_lower_r_force -
             ledger.requested_lower_r_force),
        fabs(ledger.delivered_upper_r_force -
             ledger.requested_upper_r_force)));
  const double force_scale = max(DBL_MIN,
    max(max(fabs(ledger.requested_lower_x_force),
            fabs(ledger.requested_upper_x_force)),
        max(fabs(ledger.requested_lower_r_force),
            fabs(ledger.requested_upper_r_force))));
  const double force_tolerance = 8192.*DBL_EPSILON*force_scale;
  const double power_scale = max(DBL_MIN,
    max(max(fabs(ledger.requested_lower_power),
            fabs(ledger.requested_upper_power)),
        max(fabs(ledger.delivered_lower_power),
            fabs(ledger.delivered_upper_power))));
  const double power_tolerance = 8192.*DBL_EPSILON*power_scale;
  if (swamped || !isfinite(ledger.delivered_power) ||
      ledger.maximum_component_force_residual > force_tolerance ||
      ledger.maximum_bin_force_residual > force_tolerance ||
      fabs(ledger.power_residual) > power_tolerance ||
      ledger.maximum_ring_power_residual > power_tolerance) {
    foreach_face()
      acceleration_field.x[] = before.x[];
    ledger.status = AXISYMMETRIC_MUSEHANE_TRACTION_DELIVERY_FAILED;
  }
  else {
    ledger.status = AXISYMMETRIC_MUSEHANE_TRACTION_OK;
    ledger.applied = true;
  }
  delete ((scalar *){before, lower_increment, upper_increment});
  return ledger;
}

#endif
