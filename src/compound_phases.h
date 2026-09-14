#ifndef COMPOUND_PHASES_H
#define COMPOUND_PHASES_H

#include "vof.h"

/* A single compound bubble: gas is nested inside the compound envelope. */
scalar gas_fraction[], envelope_fraction[];
scalar * interfaces = {gas_fraction, envelope_fraction};

double rho_water = 998.2, mu_water = 0.89e-3;
double rho_oil = 960.0, mu_oil = 48.0e-3;
double rho_gas = 1.20, mu_gas = 1.80e-5;

face vector compound_alpha[];
scalar compound_rho[];

typedef struct {
  double water;
  double oil;
  double gas;
  double nesting_defect;
} CompoundFractions;

static inline CompoundFractions compound_fractions (double gas, double envelope)
{
  gas = clamp(gas, 0., 1.);
  envelope = clamp(envelope, 0., 1.);
  CompoundFractions f = {
    .water = 0.,
    .oil = 0.,
    .gas = gas,
    .nesting_defect = max(gas - envelope, 0.)
  };
  f.oil = clamp(envelope - gas, 0., 1. - gas);
  f.water = max(1. - f.gas - f.oil, 0.);
  return f;
}

static inline double compound_density (CompoundFractions f)
{
  return f.water*rho_water + f.oil*rho_oil + f.gas*rho_gas;
}

static inline double compound_viscosity (CompoundFractions f)
{
  return f.water*mu_water + f.oil*mu_oil + f.gas*mu_gas;
}

event defaults (i = 0)
{
  alpha = compound_alpha;
  rho = compound_rho;
  mu = new face vector;
#if TREE
  gas_fraction.refine = gas_fraction.prolongation = fraction_refine;
  envelope_fraction.refine = envelope_fraction.prolongation = fraction_refine;
#endif
}

event properties (i++)
{
  face vector muv = mu;
  foreach_face() {
    double gas = clamp((gas_fraction[] + gas_fraction[-1])/2., 0., 1.);
    double envelope = clamp((envelope_fraction[] + envelope_fraction[-1])/2.,
                            0., 1.);
    CompoundFractions f = compound_fractions(gas, envelope);
    double density = max(compound_density(f), 1.e-12);
    compound_alpha.x[] = fm.x[]/density;
    muv.x[] = fm.x[]*max(compound_viscosity(f), 1.e-30);
  }

  foreach() {
    CompoundFractions f = compound_fractions(gas_fraction[],
                                             envelope_fraction[]);
    compound_rho[] = cm[]*max(compound_density(f), 1.e-12);
  }
  boundary({compound_rho});
}

#endif
