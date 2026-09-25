#ifndef AXISYMMETRIC_INITIAL_VELOCITY_CONDITIONER_AXI_H
#define AXISYMMETRIC_INITIAL_VELOCITY_CONDITIONER_AXI_H

#include "axisymmetric_initial_velocity_conditioner.h"

typedef struct {
  double lower;
  double upper;
} AxisymmetricConditioningPairVelocity;

typedef struct {
  int valid;
  double original_lower_velocity;
  double original_upper_velocity;
  double projected_lower_velocity;
  double projected_upper_velocity;
  double final_lower_velocity;
  double final_upper_velocity;
  double lower_coefficient;
  double upper_coefficient;
  double response_determinant;
  double maximum_metric_divergence;
  double maximum_velocity_gradient_scale;
  double relative_metric_divergence;
  mgstats baseline_projection;
  mgstats lower_basis_projection;
  mgstats upper_basis_projection;
} AxisymmetricInitialVelocityAudit;

static inline void axisymmetric_conditioning_fill_centered (
  face vector metric_flux, vector cell_velocity)
{
  foreach()
    foreach_dimension()
      cell_velocity.x[] =
        (metric_flux.x[] + metric_flux.x[1])/
        (fm.x[] + fm.x[1] + SEPS);
  boundary((scalar *){cell_velocity});
}

static inline AxisymmetricConditioningPairVelocity
axisymmetric_conditioning_measure_pair_velocity (
  vector cell_velocity, scalar lower_envelope, scalar upper_envelope)
{
  double lower_volume = 0., upper_volume = 0.;
  double lower_momentum = 0., upper_momentum = 0.;
  foreach(reduction(+:lower_volume) reduction(+:upper_volume)
          reduction(+:lower_momentum) reduction(+:upper_momentum)) {
    double lower_weight = axisymmetric_conditioning_clamp01(
      lower_envelope[])*dv();
    double upper_weight = axisymmetric_conditioning_clamp01(
      upper_envelope[])*dv();
    lower_volume += lower_weight;
    upper_volume += upper_weight;
    lower_momentum += lower_weight*cell_velocity.x[];
    upper_momentum += upper_weight*cell_velocity.x[];
  }
  AxisymmetricConditioningPairVelocity result = {
    .lower = lower_momentum/fmax(lower_volume, 1.e-300),
    .upper = upper_momentum/fmax(upper_volume, 1.e-300)
  };
  return result;
}

static inline AxisymmetricConditioningPairVelocity
axisymmetric_conditioning_measure_pair_flux_response (
  face vector metric_flux, scalar lower_envelope, scalar upper_envelope)
{
  double lower_volume = 0., upper_volume = 0.;
  double lower_response = 0., upper_response = 0.;
  foreach(reduction(+:lower_volume) reduction(+:upper_volume)
          reduction(+:lower_response) reduction(+:upper_response)) {
    double axial_velocity =
      (metric_flux.x[] + metric_flux.x[1])/
      (fm.x[] + fm.x[1] + SEPS);
    double lower_weight = axisymmetric_conditioning_clamp01(
      lower_envelope[])*dv();
    double upper_weight = axisymmetric_conditioning_clamp01(
      upper_envelope[])*dv();
    lower_volume += lower_weight;
    upper_volume += upper_weight;
    lower_response += lower_weight*axial_velocity;
    upper_response += upper_weight*axial_velocity;
  }
  AxisymmetricConditioningPairVelocity result = {
    .lower = lower_response/fmax(lower_volume, 1.e-300),
    .upper = upper_response/fmax(upper_volume, 1.e-300)
  };
  return result;
}

static inline void axisymmetric_conditioning_build_basis_flux (
  face vector metric_flux, double center, double core_radius,
  double outer_radius)
{
  foreach_face(x) {
    double plus = axisymmetric_conditioning_streamfunction(
      x, y + 0.5*Delta, center, 1., core_radius, outer_radius);
    double minus = axisymmetric_conditioning_streamfunction(
      x, y - 0.5*Delta, center, 1., core_radius, outer_radius);
    metric_flux.x[] = (plus - minus)/Delta;
  }
  foreach_face(y) {
    double plus = axisymmetric_conditioning_streamfunction(
      x + 0.5*Delta, y, center, 1., core_radius, outer_radius);
    double minus = axisymmetric_conditioning_streamfunction(
      x - 0.5*Delta, y, center, 1., core_radius, outer_radius);
    metric_flux.y[] = -(plus - minus)/Delta;
  }
  boundary((scalar *){metric_flux});
}

static inline double axisymmetric_conditioning_maximum_divergence (
  face vector metric_flux);

static inline double axisymmetric_conditioning_velocity_gradient_scale (
  vector cell_velocity)
{
  double maximum = 0.;
  foreach(reduction(max:maximum)) {
    double speed = 0.;
    foreach_dimension()
      speed += sq(cell_velocity.x[]);
    maximum = fmax(maximum, sqrt(speed)/Delta);
  }
  return maximum;
}

static inline mgstats axisymmetric_conditioning_project_flux (
  face vector metric_flux, face vector specific_volume,
  double projection_tolerance, int projection_iterations)
{
  mgstats no_projection = {0};
  if (axisymmetric_conditioning_maximum_divergence(metric_flux) <=
      projection_tolerance)
    return no_projection;
  scalar pressure_correction[];
  foreach()
    pressure_correction[] = 0.;
  boundary({pressure_correction});
  double saved_tolerance = TOLERANCE;
  int saved_iterations = NITERMAX;
  TOLERANCE = projection_tolerance;
  NITERMAX = projection_iterations;
  mgstats result = project(metric_flux, pressure_correction,
                           specific_volume, 1.);
  TOLERANCE = saved_tolerance;
  NITERMAX = saved_iterations;
  return result;
}

static inline double axisymmetric_conditioning_maximum_divergence (
  face vector metric_flux)
{
  double maximum = 0.;
  foreach(reduction(max:maximum)) {
    double divergence = 0.;
    foreach_dimension()
      divergence += metric_flux.x[1] - metric_flux.x[];
    divergence /= Delta*(cm[] + SEPS);
    maximum = fmax(maximum, fabs(divergence));
  }
  return maximum;
}

static inline AxisymmetricInitialVelocityAudit
axisymmetric_condition_upper_translation (
  face vector metric_flux, vector cell_velocity,
  scalar lower_envelope, scalar upper_envelope,
  face vector specific_volume,
  double lower_center, double upper_center,
  double lower_core_radius, double upper_core_radius,
  double transition_width, double projection_tolerance,
  int projection_iterations)
{
  AxisymmetricInitialVelocityAudit audit = {0};
  AxisymmetricConditioningPairVelocity original =
    axisymmetric_conditioning_measure_pair_velocity(
      cell_velocity, lower_envelope, upper_envelope);
  audit.original_lower_velocity = original.lower;
  audit.original_upper_velocity = original.upper;

  audit.baseline_projection = axisymmetric_conditioning_project_flux(
    metric_flux, specific_volume, projection_tolerance,
    projection_iterations);
  axisymmetric_conditioning_fill_centered(metric_flux, cell_velocity);
  AxisymmetricConditioningPairVelocity projected =
    axisymmetric_conditioning_measure_pair_velocity(
      cell_velocity, lower_envelope, upper_envelope);
  audit.projected_lower_velocity = projected.lower;
  audit.projected_upper_velocity = projected.upper;

  face vector lower_basis[], upper_basis[];
  axisymmetric_conditioning_build_basis_flux(
    lower_basis, lower_center, lower_core_radius,
    lower_core_radius + transition_width);
  axisymmetric_conditioning_build_basis_flux(
    upper_basis, upper_center, upper_core_radius,
    upper_core_radius + transition_width);
  audit.lower_basis_projection = axisymmetric_conditioning_project_flux(
    lower_basis, specific_volume, projection_tolerance,
    projection_iterations);
  audit.upper_basis_projection = axisymmetric_conditioning_project_flux(
    upper_basis, specific_volume, projection_tolerance,
    projection_iterations);

  AxisymmetricConditioningPairVelocity lower_response =
    axisymmetric_conditioning_measure_pair_flux_response(
      lower_basis, lower_envelope, upper_envelope);
  AxisymmetricConditioningPairVelocity upper_response =
    axisymmetric_conditioning_measure_pair_flux_response(
      upper_basis, lower_envelope, upper_envelope);
  double response[2][2] = {
    {lower_response.lower, upper_response.lower},
    {lower_response.upper, upper_response.upper}
  };
  double target_change[2] = {
    original.lower - projected.lower,
    -projected.upper
  };
  double coefficient[2] = {0., 0.};
  audit.response_determinant = response[0][0]*response[1][1] -
                               response[0][1]*response[1][0];
  if (!axisymmetric_conditioning_solve_2x2(
        response, target_change, coefficient, 1.e-10))
    return audit;
  audit.lower_coefficient = coefficient[0];
  audit.upper_coefficient = coefficient[1];

  foreach_face()
    metric_flux.x[] += coefficient[0]*lower_basis.x[] +
                       coefficient[1]*upper_basis.x[];
  boundary((scalar *){metric_flux});
  axisymmetric_conditioning_fill_centered(metric_flux, cell_velocity);
  AxisymmetricConditioningPairVelocity final =
    axisymmetric_conditioning_measure_pair_velocity(
      cell_velocity, lower_envelope, upper_envelope);
  audit.final_lower_velocity = final.lower;
  audit.final_upper_velocity = final.upper;
  audit.maximum_metric_divergence =
    axisymmetric_conditioning_maximum_divergence(metric_flux);
  audit.maximum_velocity_gradient_scale =
    axisymmetric_conditioning_velocity_gradient_scale(cell_velocity);
  audit.relative_metric_divergence =
    audit.maximum_velocity_gradient_scale > 0. ?
      audit.maximum_metric_divergence/
        audit.maximum_velocity_gradient_scale :
      (audit.maximum_metric_divergence == 0. ? 0. : HUGE_VAL);
  audit.valid = isfinite(audit.final_lower_velocity) &&
                isfinite(audit.final_upper_velocity) &&
                isfinite(audit.relative_metric_divergence) &&
                audit.baseline_projection.i < projection_iterations &&
                audit.lower_basis_projection.i < projection_iterations &&
                audit.upper_basis_projection.i < projection_iterations;
  return audit;
}

#endif
