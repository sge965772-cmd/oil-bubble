#ifndef PRECURSOR_STATE_H
#define PRECURSOR_STATE_H

#include "axisymmetric_projection_audit.h"

typedef struct {
  double lower_gas_volume;
  double lower_oil_volume;
  double center;
  double center_velocity;
  double moment_aspect;
  double equivalent_diameter;
  double kinetic_energy;
  double pressure_integral;
  double projected_pressure_integral;
  double nesting_volume;
  double envelope_overlap_volume;
  double gas_overlap_volume;
  double overfill_volume;
  double maximum_speed;
  double projected_volume_change;
} PrecursorState;

static inline PrecursorState measure_precursor_state (void)
{
  PrecursorState state = {0};
  double gas_volume = 0., oil_volume = 0., envelope_volume = 0.;
  double center_sum = 0., axial_second_sum = 0., radial_second_sum = 0.;
  double velocity_sum = 0., kinetic_energy = 0.;
  double pressure_integral = 0., projected_pressure_integral = 0.;
  double nesting = 0., envelope_overlap = 0., gas_overlap = 0.;
  double overfill = 0., maximum_speed = 0.;
  foreach(reduction(+:gas_volume) reduction(+:oil_volume)
          reduction(+:envelope_volume) reduction(+:center_sum)
          reduction(+:axial_second_sum) reduction(+:radial_second_sum)
          reduction(+:velocity_sum) reduction(+:kinetic_energy)
          reduction(+:pressure_integral)
          reduction(+:projected_pressure_integral)
          reduction(+:nesting) reduction(+:envelope_overlap)
          reduction(+:gas_overlap) reduction(+:overfill)
          reduction(max:maximum_speed)) {
    DualPhasePartition phase = dual_phase_partition(
      lower_gas_fraction[], lower_envelope_fraction[],
      upper_gas_fraction[], upper_envelope_fraction[]);
    double volume_element = dv();
    double lower_compound = phase.lower_gas + phase.lower_oil;
    gas_volume += phase.lower_gas*volume_element;
    oil_volume += phase.lower_oil*volume_element;
    envelope_volume += lower_compound*volume_element;
    center_sum += lower_compound*x*volume_element;
    axial_second_sum += lower_compound*sq(x)*volume_element;
    radial_second_sum += lower_compound*sq(y)*volume_element;
    velocity_sum += lower_compound*u.x[]*volume_element;
    double speed_squared = sq(u.x[]) + sq(u.y[]);
    kinetic_energy += 0.5*dual_compound_density(phase)*
                      speed_squared*volume_element;
    pressure_integral += p[]*volume_element;
    projected_pressure_integral += pf[]*volume_element;
    nesting += phase.assigned_nesting_defect*volume_element;
    envelope_overlap += phase.identity_envelope_overlap*volume_element;
    gas_overlap += phase.identity_gas_overlap*volume_element;
    overfill += phase.phase_overfill*volume_element;
    maximum_speed = max(maximum_speed, sqrt(speed_squared));
  }
  state.lower_gas_volume = 2.*pi*gas_volume;
  state.lower_oil_volume = 2.*pi*oil_volume;
  state.center = center_sum/max(envelope_volume, 1.e-300);
  state.center_velocity = velocity_sum/max(envelope_volume, 1.e-300);
  double axial_variance = axial_second_sum/max(envelope_volume, 1.e-300) -
                          sq(state.center);
  state.moment_aspect = sqrt(
    radial_second_sum/
    (2.*max(envelope_volume, 1.e-300)*max(axial_variance, 1.e-300)));
  state.equivalent_diameter = 2.*cbrt(3.*2.*pi*envelope_volume/(4.*pi));
  state.kinetic_energy = 2.*pi*kinetic_energy;
  state.pressure_integral = 2.*pi*pressure_integral;
  state.projected_pressure_integral = 2.*pi*projected_pressure_integral;
  state.nesting_volume = 2.*pi*nesting;
  state.envelope_overlap_volume = 2.*pi*envelope_overlap;
  state.gas_overlap_volume = 2.*pi*gas_overlap;
  state.overfill_volume = 2.*pi*overfill;
  state.maximum_speed = maximum_speed;
  state.projected_volume_change =
    axisymmetric_projected_cell_volume_change();
  return state;
}

#endif
