#ifndef AXISYMMETRIC_COLLISION_OBSERVABLES_H
#define AXISYMMETRIC_COLLISION_OBSERVABLES_H

typedef struct {
  double kinetic_energy;
  double surface_energy;
  double pressure_integral;
  double film_mean_pressure;
  double film_max_pressure;
  double film_support_volume;
} AxisymmetricCollisionObservables;

static inline double axisymmetric_interface_area (scalar fraction_field)
{
  double area = 0.;
  face vector no_face_fractions = {{-1}};
  foreach(reduction(+:area))
    if (fraction_field[] > 1.e-6 && fraction_field[] < 1. - 1.e-6) {
      coord normal = facet_normal(point, fraction_field, no_face_fractions);
      double alpha = plane_alpha(fraction_field[], normal);
      coord segment[2];
      if (facets(normal, alpha, segment) == 2) {
        double radius0 = max(y + segment[0].y*Delta, 0.);
        double radius1 = max(y + segment[1].y*Delta, 0.);
        double length = Delta*sqrt(
          sq(segment[1].x - segment[0].x) +
          sq(segment[1].y - segment[0].y));
        area += 2.*pi*0.5*(radius0 + radius1)*length;
      }
    }
  return area;
}

static inline AxisymmetricCollisionObservables
measure_axisymmetric_collision_observables (
  AxisymmetricContactGeometry contact, double support_half_width,
  double sigma_oil_gas, double sigma_oil_water)
{
  double kinetic_energy = 0., pressure_integral = 0.;
  double film_pressure_integral = 0., film_volume = 0.;
  double film_max_pressure = -HUGE_VAL;
  foreach(reduction(+:kinetic_energy) reduction(+:pressure_integral)
          reduction(+:film_pressure_integral) reduction(+:film_volume)
          reduction(max:film_max_pressure)) {
    DualPhasePartition phase = dual_phase_partition(
      lower_gas_fraction[], lower_envelope_fraction[],
      upper_gas_fraction[], upper_envelope_fraction[]);
    double density = dual_compound_density(phase);
    kinetic_energy += 0.5*density*(sq(u.x[]) + sq(u.y[]))*dv();
    pressure_integral += p[]*dv();
    if (contact.valid && x >= contact.lower_top && x <= contact.upper_bottom &&
        y <= support_half_width && phase.water > 0.) {
      double weight = phase.water*dv();
      film_pressure_integral += p[]*weight;
      film_volume += weight;
      film_max_pressure = max(film_max_pressure, p[]);
    }
  }

  AxisymmetricCollisionObservables result = {
    .kinetic_energy = 2.*pi*kinetic_energy,
    .surface_energy =
      sigma_oil_gas*(
        axisymmetric_interface_area(lower_gas_fraction) +
        axisymmetric_interface_area(upper_gas_fraction)) +
      sigma_oil_water*(
        axisymmetric_interface_area(lower_envelope_fraction) +
        axisymmetric_interface_area(upper_envelope_fraction)),
    .pressure_integral = 2.*pi*pressure_integral,
    .film_mean_pressure = film_volume > 0. ?
      film_pressure_integral/film_volume : NAN,
    .film_max_pressure = film_volume > 0. ? film_max_pressure : NAN,
    .film_support_volume = 2.*pi*film_volume
  };
  return result;
}

#endif
