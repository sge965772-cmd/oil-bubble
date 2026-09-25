#ifndef TRACK_A_SUBGRID_OBSERVATION_POLICY_H
#define TRACK_A_SUBGRID_OBSERVATION_POLICY_H

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef enum {
  TRACK_A_OBSERVATION_INVALID = 0,
  TRACK_A_OBSERVATION_NATIVE = 1,
  TRACK_A_OBSERVATION_SUBGRID = 2
} TrackASubgridObservationStatus;

typedef struct {
  TrackASubgridObservationStatus status;
  bool native_ordering_lost;
  double native_lower_position;
  double native_upper_position;
  double native_gap;
  double midpoint;
  double lower_position;
  double upper_position;
} TrackASubgridObservationGeometry;

typedef struct {
  bool valid;
  int ring_count;
  const unsigned char *active_n;
  const unsigned char *active_np1;
  const double *thickness;
  uint64_t release_latch_mask;
} TrackASubgridObservationState;

static inline bool track_a_observation_release_is_latched (
  uint64_t latch_mask, size_t ring)
{
  return ring < 64U && (latch_mask & (UINT64_C(1) << ring)) != 0U;
}

static inline unsigned track_a_observation_release_latch_count (
  uint64_t latch_mask)
{
  unsigned count = 0U;
  while (latch_mask) {
    count += (unsigned)(latch_mask & UINT64_C(1));
    latch_mask >>= 1U;
  }
  return count;
}

/* A released ring stays latched until the native film is resolved at the
   handoff scale. The release transition wins over a same-sample native
   crossing so that a sub-grid ring cannot reactivate on the next step. */
static inline uint64_t track_a_update_observation_release_latch (
  uint64_t latch_mask, size_t ring, bool native_resolved,
  bool active_n, bool active_np1)
{
  if (ring >= 64U)
    return latch_mask;
  const uint64_t bit = UINT64_C(1) << ring;
  if (active_n && !active_np1)
    return latch_mask | bit;
  if (native_resolved || active_np1)
    return latch_mask & ~bit;
  return latch_mask;
}

/* Active rings and released-but-not-yet-native rings retain subgrid
   observation authority.  The release latch clears on native recovery. */
static inline bool track_a_subgrid_observation_is_authoritative (
  bool active_n, bool active_np1, bool release_latched, bool native_valid)
{
  return active_n || active_np1 || (release_latched && !native_valid);
}

/* Selects coupling support coordinates without changing the native PLIC
   diagnostic.  Once a ring is active, the model thickness is authoritative;
   native crossings provide only marker identity and the local midpoint. */
static inline TrackASubgridObservationGeometry
track_a_select_observation_geometry (
  int lower_samples,
  int upper_samples,
  double native_lower_position,
  double native_upper_position,
  bool native_valid,
  bool subgrid_authoritative,
  double track_a_thickness)
{
  TrackASubgridObservationGeometry result = {
    .status = TRACK_A_OBSERVATION_INVALID,
    .native_ordering_lost = true,
    .native_lower_position = native_lower_position,
    .native_upper_position = native_upper_position,
    .native_gap = native_upper_position - native_lower_position,
    .midpoint = NAN,
    .lower_position = NAN,
    .upper_position = NAN
  };
  const bool identities_available = lower_samples > 0 && upper_samples > 0 &&
    isfinite(native_lower_position) && isfinite(native_upper_position);
  if (!identities_available)
    return result;

  result.midpoint = .5*(native_lower_position + native_upper_position);
  result.native_ordering_lost = !native_valid ||
    !(result.native_gap > 0.);
  if (subgrid_authoritative) {
    if (!isfinite(track_a_thickness) || !(track_a_thickness > 0.))
      return result;
    result.status = TRACK_A_OBSERVATION_SUBGRID;
    result.lower_position = result.midpoint - .5*track_a_thickness;
    result.upper_position = result.midpoint + .5*track_a_thickness;
    return result;
  }
  if (native_valid && result.native_gap > 0.) {
    result.status = TRACK_A_OBSERVATION_NATIVE;
    result.lower_position = native_lower_position;
    result.upper_position = native_upper_position;
    return result;
  }
  return result;
}

#endif
