#ifndef DUAL_MOMENTUM_ALGEBRA_H
#define DUAL_MOMENTUM_ALGEBRA_H

#include <math.h>

typedef struct {
  double water_base;
  double lower_envelope_correction;
  double lower_gas_correction;
  double upper_envelope_correction;
  double upper_gas_correction;
  double density;
} DualMomentumComponents;

static inline double dual_momentum_indicator_partition_defect (
  double lower_gas, double lower_envelope,
  double upper_gas, double upper_envelope)
{
  double defect = 0.;
  double values[] = {
    lower_gas, lower_envelope, upper_gas, upper_envelope
  };
  for (int index = 0; index < 4; index++) {
    defect = fmax(defect, -values[index]);
    defect = fmax(defect, values[index] - 1.);
  }
  defect = fmax(defect, lower_gas - lower_envelope);
  defect = fmax(defect, upper_gas - upper_envelope);
  defect = fmax(defect, lower_envelope + upper_envelope - 1.);
  return defect;
}

static inline int dual_momentum_indicators_are_transportable (
  double lower_gas, double lower_envelope,
  double upper_gas, double upper_envelope,
  double tolerance)
{
  if (!isfinite(lower_gas) || !isfinite(lower_envelope) ||
      !isfinite(upper_gas) || !isfinite(upper_envelope) ||
      !isfinite(tolerance) || tolerance < 0.)
    return 0;
  return dual_momentum_indicator_partition_defect(
    lower_gas, lower_envelope, upper_gas, upper_envelope) <= tolerance;
}

static inline DualMomentumComponents dual_momentum_decompose (
  double velocity,
  double lower_gas, double lower_envelope,
  double upper_gas, double upper_envelope,
  double rho_water, double rho_oil, double rho_gas)
{
  double envelope_density_increment = rho_oil - rho_water;
  double gas_density_increment = rho_gas - rho_oil;
  DualMomentumComponents components = {
    .water_base = rho_water*velocity,
    .lower_envelope_correction =
      envelope_density_increment*lower_envelope*velocity,
    .lower_gas_correction = gas_density_increment*lower_gas*velocity,
    .upper_envelope_correction =
      envelope_density_increment*upper_envelope*velocity,
    .upper_gas_correction = gas_density_increment*upper_gas*velocity,
    .density = rho_water +
      envelope_density_increment*(lower_envelope + upper_envelope) +
      gas_density_increment*(lower_gas + upper_gas)
  };
  return components;
}

static inline double dual_momentum_total (
  const DualMomentumComponents * components)
{
  return components->water_base +
         components->lower_envelope_correction +
         components->lower_gas_correction +
         components->upper_envelope_correction +
         components->upper_gas_correction;
}

static inline double dual_momentum_reconstruct_velocity (
  const DualMomentumComponents * components)
{
  return components->density > 0. ?
         dual_momentum_total(components)/components->density : NAN;
}

#endif
