#ifndef AGGREGATE_FACE_PROJECTION_POLICY_H
#define AGGREGATE_FACE_PROJECTION_POLICY_H

#include <math.h>

#include "nested_face_projection.h"

typedef enum {
  AGGREGATE_FACE_PROJECTION_REJECTED = 0,
  AGGREGATE_FACE_PROJECTION_HIGH_ORDER = 1,
  AGGREGATE_FACE_PROJECTION_LOW_ORDER_FALLBACK = 2
} AggregateFaceProjectionReason;

typedef struct {
  int valid;
  int use_projected_candidate;
  int used_low_order_fallback;
  AggregateFaceProjectionReason reason;
  double transported_correction;
} AggregateFaceProjectionDecision;

/*
 * A non-admissible high-order geometric candidate is not a failed state when
 * the donor-cell partition is valid.  In that case the conservative low-order
 * flux is the admissible FCT candidate.  Non-finite geometry, VOF bound
 * violations and invalid donor states remain hard failures.
 */
static inline AggregateFaceProjectionDecision
aggregate_face_projection_decision (
  const NestedFaceProjection * projection,
  int upwind_partition_valid, int physical_partition_valid,
  double courant, double bound_tolerance, double correction_limit)
{
  AggregateFaceProjectionDecision decision = {0};
  if (!projection || !isfinite(courant) || !isfinite(bound_tolerance) ||
      !isfinite(correction_limit) || bound_tolerance < 0. ||
      correction_limit < 0. || !upwind_partition_valid ||
      !physical_partition_valid || !projection->valid ||
      !isfinite(projection->raw_bound_defect) ||
      !isfinite(projection->correction_linf) ||
      projection->raw_bound_defect > bound_tolerance)
    return decision;

  decision.valid = 1;
  decision.transported_correction =
    fabs(courant)*projection->correction_linf;
  if (projection->correction_linf <= correction_limit) {
    decision.use_projected_candidate = 1;
    decision.reason = AGGREGATE_FACE_PROJECTION_HIGH_ORDER;
  }
  else {
    decision.used_low_order_fallback = 1;
    decision.reason = AGGREGATE_FACE_PROJECTION_LOW_ORDER_FALLBACK;
  }
  return decision;
}

#endif
