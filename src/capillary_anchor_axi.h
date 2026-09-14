#ifndef CAPILLARY_ANCHOR_AXI_H
#define CAPILLARY_ANCHOR_AXI_H

typedef struct {
  double effective_mass;
  double support_volume;
  double reaction_force;
} CapillaryAnchorLoad;

static inline CapillaryAnchorLoad capillary_anchor_measure_axi (void)
{
  double effective_mass = 0., support_volume = 0.;
  foreach(reduction(+:effective_mass) reduction(+:support_volume)) {
    DualPhasePartition phase = dual_phase_partition(
      lower_gas_fraction[], lower_envelope_fraction[],
      upper_gas_fraction[], upper_envelope_fraction[]);
    double mixture_density = dual_compound_density(phase);
    double support = dual_clip_unit(upper_envelope_fraction[]);
    effective_mass += mixture_density*support*dv();
    support_volume += support*dv();
  }

  CapillaryAnchorLoad load = {
    .effective_mass = 2.*pi*effective_mass,
    .support_volume = 2.*pi*support_volume,
    .reaction_force = 0.
  };
  return load;
}

static inline CapillaryAnchorLoad capillary_anchor_apply_axi (
  face vector acceleration_field,
  double commanded_acceleration)
{
  CapillaryAnchorLoad load = capillary_anchor_measure_axi();
  double reaction_force = load.effective_mass*commanded_acceleration;
  double force_density = reaction_force/
                         max(load.support_volume, 1.e-300);
  if (commanded_acceleration != 0.)
    foreach_face(x) {
      double support = dual_clip_unit(
        (upper_envelope_fraction[] + upper_envelope_fraction[-1])/2.);
      acceleration_field.x[] +=
        dual_compound_alpha.x[]/(fm.x[] + SEPS)*force_density*support;
    }

  load.reaction_force = reaction_force;
  return load;
}

#endif
