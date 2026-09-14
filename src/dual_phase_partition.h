#ifndef DUAL_PHASE_PARTITION_H
#define DUAL_PHASE_PARTITION_H

typedef struct {
  double water;
  double oil;
  double gas;
  double lower_oil;
  double lower_gas;
  double upper_oil;
  double upper_gas;
  double assigned_nesting_defect;
  double identity_envelope_overlap;
  double identity_gas_overlap;
  double phase_overfill;
} DualPhasePartition;

static inline double dual_clip_unit (double value)
{
  return value < 0. ? 0. : value > 1. ? 1. : value;
}

/*
 * Map two gas-core/outer-envelope identity pairs to physical phases.
 * Invalid identity overlap is diagnosed before the physical totals are
 * bounded. Gas is retained first when an invalid overfilled state occurs;
 * callers must not interpret a nonzero defect as a valid physical solution.
 */
static inline DualPhasePartition dual_phase_partition (
  double lower_gas, double lower_envelope,
  double upper_gas, double upper_envelope)
{
  lower_gas = dual_clip_unit(lower_gas);
  lower_envelope = dual_clip_unit(lower_envelope);
  upper_gas = dual_clip_unit(upper_gas);
  upper_envelope = dual_clip_unit(upper_envelope);

  DualPhasePartition p = {0};
  p.lower_gas = lower_gas;
  p.upper_gas = upper_gas;
  p.lower_oil = dual_clip_unit(lower_envelope - lower_gas);
  p.upper_oil = dual_clip_unit(upper_envelope - upper_gas);
  p.assigned_nesting_defect =
    (lower_gas > lower_envelope ? lower_gas - lower_envelope : 0.) +
    (upper_gas > upper_envelope ? upper_gas - upper_envelope : 0.);
  p.identity_envelope_overlap = lower_envelope < upper_envelope ?
                                lower_envelope : upper_envelope;
  p.identity_gas_overlap = lower_gas < upper_gas ? lower_gas : upper_gas;

  double raw_gas = lower_gas + upper_gas;
  double raw_oil = p.lower_oil + p.upper_oil;
  double raw_total = raw_gas + raw_oil;
  p.phase_overfill = raw_total > 1. ? raw_total - 1. : 0.;
  p.gas = raw_gas < 1. ? raw_gas : 1.;
  double remaining = 1. - p.gas;
  p.oil = raw_oil < remaining ? raw_oil : remaining;
  p.water = 1. - p.gas - p.oil;
  return p;
}

#endif
