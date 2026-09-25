#ifndef AXISYMMETRIC_FORCE_AUDIT_POLICY_H
#define AXISYMMETRIC_FORCE_AUDIT_POLICY_H

#include <float.h>
#include <math.h>
#include <stdbool.h>

#define AXISYMMETRIC_FORCE_AUDIT_RELATIVE_TOLERANCE (8192.*DBL_EPSILON)
#define AXISYMMETRIC_FORCE_AUDIT_DELIVERY_RELATIVE_TOLERANCE \
  (sqrt(DBL_EPSILON))

typedef struct {
  bool failed;
  bool nonfinite;
  bool swamped_increment;
  bool application_error_covered;
  bool delivered_force_retained;
  double tolerance;
} AxisymmetricForceAudit;

static inline bool axisymmetric_force_component_retained (
  double requested_force,
  double force_residual)
{
  if (requested_force == 0.)
    return force_residual == 0.;
  return fabs(force_residual) <=
    AXISYMMETRIC_FORCE_AUDIT_DELIVERY_RELATIVE_TOLERANCE*
      fabs(requested_force);
}

static inline AxisymmetricForceAudit axisymmetric_force_audit_evaluate (
  double requested_lower_force,
  double requested_upper_force,
  double lower_force_residual,
  double upper_force_residual,
  double application_force_error_l1,
  long swamped_increment_faces,
  bool upstream_nonfinite)
{
  AxisymmetricForceAudit audit = {0};
  const double requested_force_scale = fmax(DBL_MIN,
    fmax(fabs(requested_lower_force), fabs(requested_upper_force)));
  audit.nonfinite = upstream_nonfinite ||
    !isfinite(requested_lower_force) || !isfinite(requested_upper_force) ||
    !isfinite(lower_force_residual) || !isfinite(upper_force_residual) ||
    !isfinite(application_force_error_l1) || application_force_error_l1 < 0.;
  /* Bound the audited force by both reduction roundoff and the measured
     forward error from adding the acceleration increment to the face field. */
  audit.tolerance = audit.nonfinite ? NAN :
    AXISYMMETRIC_FORCE_AUDIT_RELATIVE_TOLERANCE*requested_force_scale +
    application_force_error_l1;
  audit.swamped_increment = swamped_increment_faces > 0;
  /* The per-face L1 forward-error sum is deliberately conservative and can
     be ill-conditioned relative to a near-zero startup load even when its
     signed component errors cancel. Audit the force actually delivered to
     each compound bubble directly. sqrt(epsilon) is the precision-derived
     retention floor; it rejects a lost or materially altered component while
     leaving the measured L1 sum to bound the reduction below. */
  audit.delivered_force_retained = !audit.nonfinite &&
    axisymmetric_force_component_retained(requested_lower_force,
      lower_force_residual) &&
    axisymmetric_force_component_retained(requested_upper_force,
      upper_force_residual);
  audit.application_error_covered = !audit.nonfinite &&
    fabs(lower_force_residual) <= audit.tolerance &&
    fabs(upper_force_residual) <= audit.tolerance;
  audit.failed = audit.nonfinite || !audit.delivered_force_retained ||
    !audit.application_error_covered;
  return audit;
}

#endif
