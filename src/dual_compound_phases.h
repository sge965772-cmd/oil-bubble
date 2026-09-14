#ifndef DUAL_COMPOUND_PHASES_H
#define DUAL_COMPOUND_PHASES_H

#include "vof.h"
#include "dual_phase_partition.h"

#ifndef FILTERED_PROPERTIES
# define FILTERED_PROPERTIES 0
#endif
#if FILTERED_PROPERTIES
# include "property_filter_kernel.h"
#endif

scalar lower_gas_fraction[], lower_envelope_fraction[];
scalar upper_gas_fraction[], upper_envelope_fraction[];
scalar * interfaces = {
  lower_gas_fraction,
  lower_envelope_fraction,
  upper_gas_fraction,
  upper_envelope_fraction
};

double rho_water = 998.2, mu_water = 0.89e-3;
double rho_oil = 960.0, mu_oil = 48.0e-3;
double rho_gas = 1.20, mu_gas = 1.80e-5;

face vector dual_compound_alpha[];
scalar dual_compound_rho[];

#if FILTERED_PROPERTIES
# if dimension != 2
#  error "filtered dual-compound properties currently require a 2D grid"
# endif
scalar lower_gas_property_fraction[], lower_envelope_property_fraction[];
scalar upper_gas_property_fraction[], upper_envelope_property_fraction[];
# define LOWER_GAS_PROPERTY lower_gas_property_fraction
# define LOWER_ENVELOPE_PROPERTY lower_envelope_property_fraction
# define UPPER_GAS_PROPERTY upper_gas_property_fraction
# define UPPER_ENVELOPE_PROPERTY upper_envelope_property_fraction
#else
# define LOWER_GAS_PROPERTY lower_gas_fraction
# define LOWER_ENVELOPE_PROPERTY lower_envelope_fraction
# define UPPER_GAS_PROPERTY upper_gas_fraction
# define UPPER_ENVELOPE_PROPERTY upper_envelope_fraction
#endif

static inline double dual_compound_density (DualPhasePartition p)
{
  return p.water*rho_water + p.oil*rho_oil + p.gas*rho_gas;
}

static inline double dual_compound_viscosity (DualPhasePartition p)
{
  return p.water*mu_water + p.oil*mu_oil + p.gas*mu_gas;
}

event defaults (i = 0)
{
  alpha = dual_compound_alpha;
  rho = dual_compound_rho;
  mu = new face vector;
#if TREE
  lower_gas_fraction.refine = lower_gas_fraction.prolongation = fraction_refine;
  lower_envelope_fraction.refine =
    lower_envelope_fraction.prolongation = fraction_refine;
  upper_gas_fraction.refine = upper_gas_fraction.prolongation = fraction_refine;
  upper_envelope_fraction.refine =
    upper_envelope_fraction.prolongation = fraction_refine;
# if FILTERED_PROPERTIES
  lower_gas_property_fraction.refine =
    lower_gas_property_fraction.prolongation = refine_bilinear;
  lower_envelope_property_fraction.refine =
    lower_envelope_property_fraction.prolongation = refine_bilinear;
  upper_gas_property_fraction.refine =
    upper_gas_property_fraction.prolongation = refine_bilinear;
  upper_envelope_property_fraction.refine =
    upper_envelope_property_fraction.prolongation = refine_bilinear;
# endif
#endif
}

#if FILTERED_PROPERTIES
static inline void dual_compound_update_property_filter (void)
{
  foreach() {
#define FILTER_PROPERTY(source) property_filter_2d(                    \
  source[], source[0,1], source[0,-1], source[1], source[-1],         \
  source[1,1], source[-1,1], source[1,-1], source[-1,-1])
    lower_gas_property_fraction[] = FILTER_PROPERTY(lower_gas_fraction);
    lower_envelope_property_fraction[] =
      FILTER_PROPERTY(lower_envelope_fraction);
    upper_gas_property_fraction[] = FILTER_PROPERTY(upper_gas_fraction);
    upper_envelope_property_fraction[] =
      FILTER_PROPERTY(upper_envelope_fraction);
#undef FILTER_PROPERTY
  }
  boundary((scalar *){
    lower_gas_property_fraction, lower_envelope_property_fraction,
    upper_gas_property_fraction, upper_envelope_property_fraction
  });
}
#endif

static inline void dual_compound_update_properties (void)
{
#if FILTERED_PROPERTIES
  dual_compound_update_property_filter();
#endif
  face vector muv = mu;
  foreach_face() {
    DualPhasePartition p = dual_phase_partition(
      (LOWER_GAS_PROPERTY[] + LOWER_GAS_PROPERTY[-1])/2.,
      (LOWER_ENVELOPE_PROPERTY[] + LOWER_ENVELOPE_PROPERTY[-1])/2.,
      (UPPER_GAS_PROPERTY[] + UPPER_GAS_PROPERTY[-1])/2.,
      (UPPER_ENVELOPE_PROPERTY[] + UPPER_ENVELOPE_PROPERTY[-1])/2.);
    double density = dual_compound_density(p);
    density = density > 1.e-12 ? density : 1.e-12;
    double viscosity = dual_compound_viscosity(p);
    viscosity = viscosity > 1.e-30 ? viscosity : 1.e-30;
    dual_compound_alpha.x[] = fm.x[]/density;
    muv.x[] = fm.x[]*viscosity;
  }

  foreach() {
    DualPhasePartition p = dual_phase_partition(
      LOWER_GAS_PROPERTY[], LOWER_ENVELOPE_PROPERTY[],
      UPPER_GAS_PROPERTY[], UPPER_ENVELOPE_PROPERTY[]);
    double density = dual_compound_density(p);
    dual_compound_rho[] = cm[]*(density > 1.e-12 ? density : 1.e-12);
  }
  boundary({dual_compound_rho});
}

event properties (i++)
{
  dual_compound_update_properties();
}

#endif
