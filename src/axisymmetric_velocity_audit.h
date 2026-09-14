#ifndef AXISYMMETRIC_VELOCITY_AUDIT_H
#define AXISYMMETRIC_VELOCITY_AUDIT_H

typedef struct {
  double maximum_speed;
  double maximum_x;
  double maximum_y;
  double maximum_water_fraction;
  double maximum_oil_fraction;
  double maximum_gas_fraction;
  double lower_outer_interface_velocity_rms;
  double lower_outer_interface_velocity_max;
  double lower_gas_interface_velocity_rms;
  double lower_gas_interface_velocity_max;
  double upper_outer_interface_velocity_rms;
  double upper_outer_interface_velocity_max;
  double upper_gas_interface_velocity_rms;
  double upper_gas_interface_velocity_max;
} AxisymmetricVelocityAudit;

/* True PLIC segment area after revolution about the x axis. */
static inline double axisymmetric_plic_area_element (Point point, scalar c)
{
  if (c[] <= 1.e-6 || c[] >= 1. - 1.e-6)
    return 0.;
  coord normal = interface_normal(point, c), centroid;
  double alpha = plane_alpha(c[], normal);
  double segment_length = plane_area_center(normal, alpha, &centroid);
  double radius = max(y + centroid.y*Delta, 0.);
  return 2.*pi*radius*Delta*segment_length;
}

static inline AxisymmetricVelocityAudit
measure_axisymmetric_velocity_audit (void)
{
  double maximum_speed = 0.;
  foreach(reduction(max:maximum_speed))
    maximum_speed = max(maximum_speed, sqrt(sq(u.x[]) + sq(u.y[])));

  double peak_weight = 0., peak_x = 0., peak_y = 0.;
  double peak_water = 0., peak_oil = 0., peak_gas = 0.;
  double lower_outer_area = 0., lower_outer_v2 = 0.;
  double lower_gas_area = 0., lower_gas_v2 = 0.;
  double upper_outer_area = 0., upper_outer_v2 = 0.;
  double upper_gas_area = 0., upper_gas_v2 = 0.;
  double lower_outer_vmax = 0., lower_gas_vmax = 0.;
  double upper_outer_vmax = 0., upper_gas_vmax = 0.;
  foreach(reduction(+:peak_weight) reduction(+:peak_x)
          reduction(+:peak_y) reduction(+:peak_water)
          reduction(+:peak_oil) reduction(+:peak_gas)
          reduction(+:lower_outer_area) reduction(+:lower_outer_v2)
          reduction(+:lower_gas_area) reduction(+:lower_gas_v2)
          reduction(+:upper_outer_area) reduction(+:upper_outer_v2)
          reduction(+:upper_gas_area) reduction(+:upper_gas_v2)
          reduction(max:lower_outer_vmax) reduction(max:lower_gas_vmax)
          reduction(max:upper_outer_vmax) reduction(max:upper_gas_vmax)) {
    double speed_squared = sq(u.x[]) + sq(u.y[]);
    double speed = sqrt(speed_squared);
    if (fabs(speed - maximum_speed) <=
        1.e-12*max(1., maximum_speed)) {
      DualPhasePartition phase = dual_phase_partition(
        lower_gas_fraction[], lower_envelope_fraction[],
        upper_gas_fraction[], upper_envelope_fraction[]);
      double weight = dv();
      peak_weight += weight;
      peak_x += x*weight;
      peak_y += y*weight;
      peak_water += phase.water*weight;
      peak_oil += phase.oil*weight;
      peak_gas += phase.gas*weight;
    }

    double area = axisymmetric_plic_area_element(
      point, lower_envelope_fraction);
    lower_outer_area += area;
    lower_outer_v2 += area*speed_squared;
    if (area > 0.)
      lower_outer_vmax = max(lower_outer_vmax, speed);

    area = axisymmetric_plic_area_element(point, lower_gas_fraction);
    lower_gas_area += area;
    lower_gas_v2 += area*speed_squared;
    if (area > 0.)
      lower_gas_vmax = max(lower_gas_vmax, speed);

    area = axisymmetric_plic_area_element(
      point, upper_envelope_fraction);
    upper_outer_area += area;
    upper_outer_v2 += area*speed_squared;
    if (area > 0.)
      upper_outer_vmax = max(upper_outer_vmax, speed);

    area = axisymmetric_plic_area_element(point, upper_gas_fraction);
    upper_gas_area += area;
    upper_gas_v2 += area*speed_squared;
    if (area > 0.)
      upper_gas_vmax = max(upper_gas_vmax, speed);
  }

  AxisymmetricVelocityAudit audit = {
    .maximum_speed = maximum_speed,
    .maximum_x = peak_x/max(peak_weight, 1.e-300),
    .maximum_y = peak_y/max(peak_weight, 1.e-300),
    .maximum_water_fraction = peak_water/max(peak_weight, 1.e-300),
    .maximum_oil_fraction = peak_oil/max(peak_weight, 1.e-300),
    .maximum_gas_fraction = peak_gas/max(peak_weight, 1.e-300),
    .lower_outer_interface_velocity_rms = sqrt(
      lower_outer_v2/max(lower_outer_area, 1.e-300)),
    .lower_outer_interface_velocity_max = lower_outer_vmax,
    .lower_gas_interface_velocity_rms = sqrt(
      lower_gas_v2/max(lower_gas_area, 1.e-300)),
    .lower_gas_interface_velocity_max = lower_gas_vmax,
    .upper_outer_interface_velocity_rms = sqrt(
      upper_outer_v2/max(upper_outer_area, 1.e-300)),
    .upper_outer_interface_velocity_max = upper_outer_vmax,
    .upper_gas_interface_velocity_rms = sqrt(
      upper_gas_v2/max(upper_gas_area, 1.e-300)),
    .upper_gas_interface_velocity_max = upper_gas_vmax
  };
  return audit;
}

#endif
