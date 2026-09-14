#ifndef POSITIVE_MATERIAL_PARTITION_H
#define POSITIVE_MATERIAL_PARTITION_H

#include <math.h>

typedef struct {
  double water;
  double lower_oil;
  double lower_gas;
  double upper_oil;
  double upper_gas;
  double bound_defect;
  double identity_defect;
  double water_defect;
  double lower_oil_defect;
  double upper_oil_defect;
  double defect;
  double correction_l1;
  int valid;
} PositiveMaterialPartition;

typedef struct {
  double water;
  double oil;
  double gas;
  double bound_defect;
  double identity_defect;
  double correction_l1;
  int valid;
} AggregatePhysicalPartition;

enum {
  POSITIVE_MATERIAL_WATER = 0,
  POSITIVE_MATERIAL_LOWER_OIL = 1,
  POSITIVE_MATERIAL_LOWER_GAS = 2,
  POSITIVE_MATERIAL_UPPER_OIL = 3,
  POSITIVE_MATERIAL_UPPER_GAS = 4,
  POSITIVE_MATERIAL_COUNT = 5
};

enum {
  AGGREGATE_MATERIAL_WATER = 0,
  AGGREGATE_MATERIAL_OIL = 1,
  AGGREGATE_MATERIAL_GAS = 2,
  AGGREGATE_MATERIAL_COUNT = 3
};

static inline double positive_material_fraction (
  const PositiveMaterialPartition * partition, int material)
{
  switch (material) {
  case POSITIVE_MATERIAL_WATER: return partition->water;
  case POSITIVE_MATERIAL_LOWER_OIL: return partition->lower_oil;
  case POSITIVE_MATERIAL_LOWER_GAS: return partition->lower_gas;
  case POSITIVE_MATERIAL_UPPER_OIL: return partition->upper_oil;
  case POSITIVE_MATERIAL_UPPER_GAS: return partition->upper_gas;
  default: return NAN;
  }
}

static inline int positive_material_dominant (
  const PositiveMaterialPartition * partition)
{
  int dominant = POSITIVE_MATERIAL_WATER;
  double maximum = partition->water;
  for (int material = 1; material < POSITIVE_MATERIAL_COUNT; material++) {
    double fraction = positive_material_fraction(partition, material);
    if (fraction > maximum) {
      maximum = fraction;
      dominant = material;
    }
  }
  return dominant;
}

static inline double positive_material_sum (
  const PositiveMaterialPartition * partition)
{
  return partition->water + partition->lower_oil + partition->lower_gas +
         partition->upper_oil + partition->upper_gas;
}

/*
 * Convert two nested gas/envelope identity pairs into five positive material
 * fractions. The uncorrected algebra sums exactly to one. Roundoff-scale
 * negative fractions are projected conservatively onto the simplex; the L1
 * correction is reported and larger defects reject the partition.
 */
static inline PositiveMaterialPartition
positive_material_partition_with_tolerances (
  double lower_gas, double lower_envelope,
  double upper_gas, double upper_envelope,
  double identity_tolerance, double bound_tolerance)
{
  PositiveMaterialPartition partition = {0};
  double inputs[] = {
    lower_gas, lower_envelope, upper_gas, upper_envelope,
    identity_tolerance, bound_tolerance
  };
  for (int index = 0; index < 6; index++)
    if (!isfinite(inputs[index])) {
      partition.bound_defect = INFINITY;
      partition.identity_defect = INFINITY;
      partition.water_defect = INFINITY;
      partition.lower_oil_defect = INFINITY;
      partition.upper_oil_defect = INFINITY;
      partition.defect = INFINITY;
      return partition;
    }
  if (identity_tolerance < 0. || bound_tolerance < 0.) {
    partition.bound_defect = INFINITY;
    partition.identity_defect = INFINITY;
    partition.water_defect = INFINITY;
    partition.lower_oil_defect = INFINITY;
    partition.upper_oil_defect = INFINITY;
    partition.defect = INFINITY;
    return partition;
  }

  double clipped[4];
  double bound_correction = 0.;
  for (int index = 0; index < 4; index++) {
    partition.bound_defect = fmax(partition.bound_defect, -inputs[index]);
    partition.bound_defect = fmax(
      partition.bound_defect, inputs[index] - 1.);
    clipped[index] = fmin(1., fmax(0., inputs[index]));
    bound_correction += fabs(clipped[index] - inputs[index]);
  }
  partition.defect = partition.bound_defect;
  if (partition.bound_defect > bound_tolerance)
    return partition;

  double values[] = {
    1. - clipped[1] - clipped[3],
    clipped[1] - clipped[0],
    clipped[0],
    clipped[3] - clipped[2],
    clipped[2]
  };
  partition.water_defect = fmax(0., -values[0]);
  partition.lower_oil_defect = fmax(0., -values[1]);
  partition.upper_oil_defect = fmax(0., -values[3]);
  for (int index = 0; index < 5; index++)
    partition.identity_defect = fmax(
      partition.identity_defect, -values[index]);
  partition.defect = fmax(
    partition.bound_defect, partition.identity_defect);
  if (partition.identity_defect > identity_tolerance)
    return partition;

  double removed_negative = 0.;
  int largest = 0;
  for (int index = 0; index < 5; index++) {
    if (values[index] < 0.) {
      removed_negative -= values[index];
      values[index] = 0.;
    }
    if (values[index] > values[largest])
      largest = index;
  }
  values[largest] -= removed_negative;
  if (values[largest] < 0.) {
    partition.identity_defect = fmax(
      partition.identity_defect, -values[largest]);
    partition.defect = fmax(
      partition.bound_defect, partition.identity_defect);
    return partition;
  }

  partition.water = values[0];
  partition.lower_oil = values[1];
  partition.lower_gas = values[2];
  partition.upper_oil = values[3];
  partition.upper_gas = values[4];
  partition.correction_l1 = bound_correction + 2.*removed_negative;
  partition.valid = 1;
  return partition;
}

static inline PositiveMaterialPartition positive_material_partition (
  double lower_gas, double lower_envelope,
  double upper_gas, double upper_envelope,
  double tolerance)
{
  return positive_material_partition_with_tolerances(
    lower_gas, lower_envelope, upper_gas, upper_envelope,
    tolerance, tolerance);
}

static inline double positive_material_density (
  const PositiveMaterialPartition * partition,
  double rho_water, double rho_oil, double rho_gas)
{
  return rho_water*partition->water +
         rho_oil*(partition->lower_oil + partition->upper_oil) +
         rho_gas*(partition->lower_gas + partition->upper_gas);
}

static inline AggregatePhysicalPartition aggregate_physical_partition (
  const PositiveMaterialPartition * identity_partition)
{
  AggregatePhysicalPartition aggregate = {0};
  if (!identity_partition)
    return aggregate;
  aggregate.water = identity_partition->water;
  aggregate.oil = identity_partition->lower_oil +
                  identity_partition->upper_oil;
  aggregate.gas = identity_partition->lower_gas +
                  identity_partition->upper_gas;
  aggregate.bound_defect = identity_partition->bound_defect;
  aggregate.identity_defect = identity_partition->identity_defect;
  aggregate.correction_l1 = identity_partition->correction_l1;
  aggregate.valid = identity_partition->valid;
  return aggregate;
}

static inline double aggregate_physical_sum (
  const AggregatePhysicalPartition * partition)
{
  return partition->water + partition->oil + partition->gas;
}

static inline double aggregate_physical_fraction (
  const AggregatePhysicalPartition * partition, int material)
{
  switch (material) {
  case AGGREGATE_MATERIAL_WATER: return partition->water;
  case AGGREGATE_MATERIAL_OIL: return partition->oil;
  case AGGREGATE_MATERIAL_GAS: return partition->gas;
  default: return NAN;
  }
}

static inline int aggregate_physical_dominant (
  const AggregatePhysicalPartition * partition)
{
  int dominant = AGGREGATE_MATERIAL_WATER;
  double maximum = partition->water;
  for (int material = 1; material < AGGREGATE_MATERIAL_COUNT; material++) {
    double fraction = aggregate_physical_fraction(partition, material);
    if (fraction > maximum) {
      maximum = fraction;
      dominant = material;
    }
  }
  return dominant;
}

static inline double aggregate_physical_density (
  const AggregatePhysicalPartition * partition,
  double rho_water, double rho_oil, double rho_gas)
{
  return rho_water*partition->water + rho_oil*partition->oil +
         rho_gas*partition->gas;
}

#endif
