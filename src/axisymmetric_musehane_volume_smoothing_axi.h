#ifndef AXISYMMETRIC_MUSEHANE_VOLUME_SMOOTHING_AXI_H
#define AXISYMMETRIC_MUSEHANE_VOLUME_SMOOTHING_AXI_H

#include <math.h>

#include "poisson.h"

typedef enum {
  AXISYMMETRIC_MUSEHANE_VOLUME_SMOOTHING_OK = 0,
  AXISYMMETRIC_MUSEHANE_VOLUME_SMOOTHING_INVALID_INPUT = 1,
  AXISYMMETRIC_MUSEHANE_VOLUME_SMOOTHING_NOT_CONVERGED = 2
} AxisymmetricMusehaneVolumeSmoothingStatus;

typedef struct {
  AxisymmetricMusehaneVolumeSmoothingStatus status;
  bool applied;
  int axial_iterations;
  int radial_iterations;
  double maximum_input_magnitude;
  double maximum_output_magnitude;
  double maximum_relative_residual;
} AxisymmetricMusehaneVolumeSmoothingAudit;

/* Diagnostic implementation of the volume-grid interpretation of Musehane
   et al. (2018), equation (25). In axisymmetric coordinates the Helmholtz
   operator must carry Basilisk's face and cell metrics:

     div(eta grad(g_s)) - g_s = -g,

   so alpha=eta*fm, lambda=-cm and b=-cm*g. The operation is deliberately
   separate from the 1-D surface-mesh smoother used by the film core. */
static inline AxisymmetricMusehaneVolumeSmoothingStatus
smooth_axisymmetric_musehane_pressure_driver (
  vector pressure_gradient, vector capillary_force_density,
  vector smoothed, double smoothing_length_squared,
  double relative_tolerance,
  AxisymmetricMusehaneVolumeSmoothingAudit *audit)
{
  if (!audit || !isfinite(smoothing_length_squared) ||
      smoothing_length_squared < 0. || !isfinite(relative_tolerance) ||
      relative_tolerance <= 0. ||
      pressure_gradient.x.i == capillary_force_density.x.i ||
      pressure_gradient.x.i == capillary_force_density.y.i ||
      pressure_gradient.y.i == capillary_force_density.x.i ||
      pressure_gradient.y.i == capillary_force_density.y.i ||
      smoothed.x.i == pressure_gradient.x.i ||
      smoothed.x.i == pressure_gradient.y.i ||
      smoothed.y.i == pressure_gradient.x.i ||
      smoothed.y.i == pressure_gradient.y.i ||
      smoothed.x.i == capillary_force_density.x.i ||
      smoothed.x.i == capillary_force_density.y.i ||
      smoothed.y.i == capillary_force_density.x.i ||
      smoothed.y.i == capillary_force_density.y.i) {
    if (audit)
      *audit = (AxisymmetricMusehaneVolumeSmoothingAudit) {
        .status = AXISYMMETRIC_MUSEHANE_VOLUME_SMOOTHING_INVALID_INPUT
      };
    return AXISYMMETRIC_MUSEHANE_VOLUME_SMOOTHING_INVALID_INPUT;
  }

  scalar corrected_x[], corrected_r[], rhs_x[], rhs_r[];
  scalar axial_lambda[], radial_lambda[];
  face vector smoothing_alpha[];
  double maximum_input_magnitude = 0.;
  foreach(reduction(max:maximum_input_magnitude)) {
    corrected_x[] = pressure_gradient.x[] - capillary_force_density.x[];
    corrected_r[] = pressure_gradient.y[] - capillary_force_density.y[];
    const double magnitude = hypot(corrected_x[], corrected_r[]);
    maximum_input_magnitude = max(maximum_input_magnitude, magnitude);
    smoothed.x[] = corrected_x[];
    smoothed.y[] = corrected_r[];
    axial_lambda[] = -cm[];
    /* The radial component is a component of a vector, not an axisymmetric
     * scalar. Its Laplacian contains -g_r/r^2. */
    radial_lambda[] = -cm[]*(1. + smoothing_length_squared/sq(y));
    rhs_x[] = -cm[]*corrected_x[];
    rhs_r[] = -cm[]*corrected_r[];
  }
  foreach_face()
    smoothing_alpha.x[] = fm.x[]*smoothing_length_squared;
  boundary({corrected_x, corrected_r, rhs_x, rhs_r,
    axial_lambda, radial_lambda});
  boundary((scalar *) {smoothing_alpha, smoothed});

  mgstats axial = {0}, radial = {0};
  double maximum_relative_residual = 0.;
  if (smoothing_length_squared > 0.) {
    const double scale = max(maximum_input_magnitude, 1.);
    const double absolute_tolerance = relative_tolerance*scale;
    axial = poisson(smoothed.x, rhs_x,
      alpha = smoothing_alpha, lambda = axial_lambda,
      tolerance = absolute_tolerance, nrelax = 4);
    radial = poisson(smoothed.y, rhs_r,
      alpha = smoothing_alpha, lambda = radial_lambda,
      tolerance = absolute_tolerance, nrelax = 4);
    maximum_relative_residual = max(axial.resa, radial.resa)/scale;
  }
  boundary((scalar *) {smoothed});

  double maximum_output_magnitude = 0.;
  foreach(reduction(max:maximum_output_magnitude))
    maximum_output_magnitude = max(maximum_output_magnitude,
      hypot(smoothed.x[], smoothed.y[]));

  const bool converged = isfinite(maximum_input_magnitude) &&
    isfinite(maximum_output_magnitude) &&
    isfinite(maximum_relative_residual) &&
    maximum_relative_residual <= relative_tolerance;
  *audit = (AxisymmetricMusehaneVolumeSmoothingAudit) {
    .status = converged ? AXISYMMETRIC_MUSEHANE_VOLUME_SMOOTHING_OK :
      AXISYMMETRIC_MUSEHANE_VOLUME_SMOOTHING_NOT_CONVERGED,
    .applied = smoothing_length_squared > 0.,
    .axial_iterations = axial.i,
    .radial_iterations = radial.i,
    .maximum_input_magnitude = maximum_input_magnitude,
    .maximum_output_magnitude = maximum_output_magnitude,
    .maximum_relative_residual = maximum_relative_residual
  };
  delete({corrected_x, corrected_r, rhs_x, rhs_r,
    axial_lambda, radial_lambda});
  delete((scalar *) {smoothing_alpha});
  return audit->status;
}

#endif
