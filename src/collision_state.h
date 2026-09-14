#ifndef COLLISION_STATE_H
#define COLLISION_STATE_H

#include "axisymmetric_projection_audit.h"

typedef struct {
  double gas_volumes[2];
  double oil_volumes[2];
  double envelope_volumes[2];
  double centers[2];
  double center_velocities[2];
  double axial_radii[2];
  double moment_aspects[2];
  double outer_gap;
  double nesting_volume;
  double envelope_overlap_volume;
  double gas_overlap_volume;
  double overfill_volume;
  double maximum_speed;
  double projected_volume_change;
} CollisionState;

static inline CollisionState measure_collision_state (void)
{
  CollisionState state = {0};
  double gas0 = 0., oil0 = 0., gas1 = 0., oil1 = 0.;
  double envelope0 = 0., envelope1 = 0.;
  double center0 = 0., center1 = 0.;
  double axial20 = 0., axial21 = 0., radial20 = 0., radial21 = 0.;
  double velocity0 = 0., velocity1 = 0.;
  double nesting = 0., envelope_overlap = 0., gas_overlap = 0.;
  double overfill = 0., maximum_speed = 0.;
  foreach(reduction(+:gas0) reduction(+:oil0)
          reduction(+:gas1) reduction(+:oil1)
          reduction(+:envelope0) reduction(+:envelope1)
          reduction(+:center0) reduction(+:center1)
          reduction(+:axial20) reduction(+:axial21)
          reduction(+:radial20) reduction(+:radial21)
          reduction(+:velocity0) reduction(+:velocity1)
          reduction(+:nesting) reduction(+:envelope_overlap)
          reduction(+:gas_overlap) reduction(+:overfill)
          reduction(max:maximum_speed)) {
    DualPhasePartition phase = dual_phase_partition(
      lower_gas_fraction[], lower_envelope_fraction[],
      upper_gas_fraction[], upper_envelope_fraction[]);
    double volume_element = dv();
    double compound0 = phase.lower_gas + phase.lower_oil;
    double compound1 = phase.upper_gas + phase.upper_oil;
    gas0 += phase.lower_gas*volume_element;
    oil0 += phase.lower_oil*volume_element;
    gas1 += phase.upper_gas*volume_element;
    oil1 += phase.upper_oil*volume_element;
    envelope0 += compound0*volume_element;
    envelope1 += compound1*volume_element;
    center0 += compound0*x*volume_element;
    center1 += compound1*x*volume_element;
    axial20 += compound0*sq(x)*volume_element;
    axial21 += compound1*sq(x)*volume_element;
    radial20 += compound0*sq(y)*volume_element;
    radial21 += compound1*sq(y)*volume_element;
    velocity0 += compound0*u.x[]*volume_element;
    velocity1 += compound1*u.x[]*volume_element;
    nesting += phase.assigned_nesting_defect*volume_element;
    envelope_overlap += phase.identity_envelope_overlap*volume_element;
    gas_overlap += phase.identity_gas_overlap*volume_element;
    overfill += phase.phase_overfill*volume_element;
    maximum_speed = max(maximum_speed, sqrt(sq(u.x[]) + sq(u.y[])));
  }
  double envelopes[2] = {envelope0, envelope1};
  double centers[2] = {
    center0/max(envelope0, 1.e-300),
    center1/max(envelope1, 1.e-300)
  };
  double axial_second[2] = {axial20, axial21};
  double radial_second[2] = {radial20, radial21};
  state.gas_volumes[0] = 2.*pi*gas0;
  state.oil_volumes[0] = 2.*pi*oil0;
  state.gas_volumes[1] = 2.*pi*gas1;
  state.oil_volumes[1] = 2.*pi*oil1;
  for (int identity = 0; identity < 2; identity++) {
    state.envelope_volumes[identity] = 2.*pi*envelopes[identity];
    if (envelopes[identity] <= 1.e-250) {
      state.centers[identity] = 0.;
      state.center_velocities[identity] = 0.;
      state.axial_radii[identity] = 0.;
      state.moment_aspects[identity] = 0.;
      continue;
    }
    state.centers[identity] = centers[identity];
    state.center_velocities[identity] = identity == 0 ?
      velocity0/max(envelope0, 1.e-300) :
      velocity1/max(envelope1, 1.e-300);
    double axial_variance =
      axial_second[identity]/max(envelopes[identity], 1.e-300) -
      sq(centers[identity]);
    state.axial_radii[identity] = sqrt(5.*max(axial_variance, 1.e-300));
    state.moment_aspects[identity] = sqrt(
      radial_second[identity]/
      (2.*envelopes[identity]*
       max(axial_variance, 1.e-300)));
  }
  state.outer_gap = state.centers[1] - state.centers[0] -
                    state.axial_radii[0] - state.axial_radii[1];
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
