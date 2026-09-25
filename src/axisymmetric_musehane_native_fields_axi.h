#ifndef AXISYMMETRIC_MUSEHANE_NATIVE_FIELDS_AXI_H
#define AXISYMMETRIC_MUSEHANE_NATIVE_FIELDS_AXI_H

#include <math.h>

typedef enum {
  AXISYMMETRIC_MUSEHANE_NATIVE_FIELD_OK = 0,
  AXISYMMETRIC_MUSEHANE_NATIVE_FIELD_INVALID_INPUT = 1,
  AXISYMMETRIC_MUSEHANE_NATIVE_FIELD_NO_CAPILLARY_INTERFACE = 2,
  AXISYMMETRIC_MUSEHANE_NATIVE_FIELD_NONFINITE = 3
} AxisymmetricMusehaneNativeFieldStatus;

typedef struct {
  int interfacial_field_count;
  long nonfinite_cells;
  double maximum_pressure_gradient;
  double maximum_capillary_force_density;
} AxisymmetricMusehaneNativeFieldAudit;

/* Build cell-centred diagnostic fields from the same face stencils used by
 * centered.h and iforce.h.  The source VOF fields are copied before clipping
 * or changing prolongation, so this observer cannot alter DNS state. */
static inline AxisymmetricMusehaneNativeFieldStatus
build_axisymmetric_musehane_native_fields (
  scalar pressure,
  vector pressure_gradient,
  vector capillary_force_density,
  AxisymmetricMusehaneNativeFieldAudit *output_audit)
{
  if (!output_audit ||
      pressure_gradient.x.i == capillary_force_density.x.i ||
      pressure_gradient.x.i == capillary_force_density.y.i ||
      pressure_gradient.y.i == capillary_force_density.x.i ||
      pressure_gradient.y.i == capillary_force_density.y.i)
    return AXISYMMETRIC_MUSEHANE_NATIVE_FIELD_INVALID_INPUT;

  face vector pressure_gradient_face[], capillary_force_face[];
  foreach_face() {
    pressure_gradient_face.x[] = fm.x[] > 0. ?
      (pressure[] - pressure[-1])/Delta : 0.;
    capillary_force_face.x[] = 0.;
  }

  int interfacial_field_count = 0;
  for (scalar source_fraction in interfaces)
    if (source_fraction.sigma) {
      interfacial_field_count++;
      scalar marker = new scalar;
      scalar potential = new scalar;
#if TREE
      marker.refine = marker.prolongation = fraction_refine;
#endif
      foreach()
        marker[] = clamp(source_fraction[], 0., 1.);
      boundary({marker});

      curvature(source_fraction, potential, source_fraction.sigma,
        add = false);
#if TREE
      marker.prolongation = pressure.prolongation;
      marker.dirty = true;
#endif
      foreach_face()
        if (marker[] != marker[-1] && fm.x[] > 0.) {
          const double potential_face =
            potential[] < nodata && potential[-1] < nodata ?
              .5*(potential[] + potential[-1]) :
            potential[] < nodata ? potential[] :
            potential[-1] < nodata ? potential[-1] : 0.;
          capillary_force_face.x[] += potential_face*
            (marker[] - marker[-1])/Delta;
        }
      delete({marker, potential});
    }

  if (interfacial_field_count == 0)
    return AXISYMMETRIC_MUSEHANE_NATIVE_FIELD_NO_CAPILLARY_INTERFACE;

  vector candidate_pressure_gradient[], candidate_capillary_force[];
  foreach() {
    foreach_dimension() {
      const double metric_sum = fm.x[] + fm.x[1];
      candidate_pressure_gradient.x[] = metric_sum > 0. ?
        (fm.x[]*pressure_gradient_face.x[] +
         fm.x[1]*pressure_gradient_face.x[1])/metric_sum : 0.;
      candidate_capillary_force.x[] = metric_sum > 0. ?
        (fm.x[]*capillary_force_face.x[] +
         fm.x[1]*capillary_force_face.x[1])/metric_sum : 0.;
    }
  }
  boundary((scalar *) {
    candidate_pressure_gradient,
    candidate_capillary_force
  });

  long nonfinite_cells = 0;
  double maximum_pressure_gradient = 0.;
  double maximum_capillary_force_density = 0.;
  foreach(reduction(+:nonfinite_cells)
          reduction(max:maximum_pressure_gradient)
          reduction(max:maximum_capillary_force_density)) {
    double pressure_norm2 = 0., capillary_norm2 = 0.;
    int finite_cell = 1;
    foreach_dimension() {
      const double pressure_component = candidate_pressure_gradient.x[];
      const double capillary_component = candidate_capillary_force.x[];
      if (!isfinite(pressure_component) || !isfinite(capillary_component))
        finite_cell = 0;
      pressure_norm2 += sq(pressure_component);
      capillary_norm2 += sq(capillary_component);
    }
    if (!isfinite(pressure_norm2) || !isfinite(capillary_norm2))
      finite_cell = 0;
    if (!finite_cell)
      nonfinite_cells++;
    else {
      maximum_pressure_gradient =
        max(maximum_pressure_gradient, sqrt(pressure_norm2));
      maximum_capillary_force_density =
        max(maximum_capillary_force_density, sqrt(capillary_norm2));
    }
  }

  if (nonfinite_cells != 0)
    return AXISYMMETRIC_MUSEHANE_NATIVE_FIELD_NONFINITE;

  foreach()
    foreach_dimension() {
      pressure_gradient.x[] = candidate_pressure_gradient.x[];
      capillary_force_density.x[] = candidate_capillary_force.x[];
    }
  boundary((scalar *) {pressure_gradient, capillary_force_density});
  *output_audit = (AxisymmetricMusehaneNativeFieldAudit) {
    .interfacial_field_count = interfacial_field_count,
    .nonfinite_cells = nonfinite_cells,
    .maximum_pressure_gradient = maximum_pressure_gradient,
    .maximum_capillary_force_density = maximum_capillary_force_density
  };
  return AXISYMMETRIC_MUSEHANE_NATIVE_FIELD_OK;
}

#endif
