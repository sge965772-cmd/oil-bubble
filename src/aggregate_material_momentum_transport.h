#ifndef AGGREGATE_MATERIAL_MOMENTUM_TRANSPORT_H
#define AGGREGATE_MATERIAL_MOMENTUM_TRANSPORT_H

#include <time.h>
#include "positive_material_partition.h"
#include "nested_face_projection.h"
#include "nested_plic_common_normal.h"
#include "aggregate_mass_weighting.h"

#if dimension != 2
# error "aggregate material momentum transport currently requires 2D axisymmetry"
#endif

#ifndef POSITIVE_MATERIAL_PARTITION_TOLERANCE
# define POSITIVE_MATERIAL_PARTITION_TOLERANCE 1.e-8
#endif
#ifndef POSITIVE_MATERIAL_VOF_BOUND_TOLERANCE
# define POSITIVE_MATERIAL_VOF_BOUND_TOLERANCE 1.e-3
#endif
#ifndef AGGREGATE_MATERIAL_FACE_PROJECTION_LIMIT
# define AGGREGATE_MATERIAL_FACE_PROJECTION_LIMIT 1.e-4
#endif
#ifndef AGGREGATE_MATERIAL_CELL_DIAGNOSTICS
# define AGGREGATE_MATERIAL_CELL_DIAGNOSTICS 0
#endif
#ifndef AGGREGATE_MATERIAL_FALLBACK_CELL_DIAGNOSTICS
# define AGGREGATE_MATERIAL_FALLBACK_CELL_DIAGNOSTICS 0
#endif
#ifndef AGGREGATE_MATERIAL_SIMPLEX_FLUX_LIMITING
# define AGGREGATE_MATERIAL_SIMPLEX_FLUX_LIMITING 0
#endif
#if AGGREGATE_MATERIAL_SIMPLEX_FLUX_LIMITING
# include "simplex_face_flux_limiter.h"
# include "aggregate_face_projection_policy.h"
#endif

long dual_momentum_invalid_cell_count = 0;
long dual_momentum_invalid_cell_count_max = 0;
long dual_momentum_invalid_face_count = 0;
long dual_momentum_invalid_face_count_max = 0;
double dual_momentum_max_face_partition_defect = 0.;
double dual_momentum_max_cell_partition_defect = 0.;
double dual_momentum_max_partition_correction = 0.;
double dual_momentum_max_vof_bound_defect = 0.;
int dual_momentum_transport_failed = 0;
double dual_momentum_last_transport_cpu_seconds = 0.;
long dual_momentum_limited_face_count = 0;
long dual_momentum_limited_face_count_max = 0;
double dual_momentum_min_flux_alpha = 1.;
double dual_momentum_min_flux_alpha_ever = 1.;
double dual_momentum_max_simplex_closure_defect = 0.;
double dual_momentum_max_flux_fraction_correction = 0.;
double dual_momentum_max_flux_courant_correction = 0.;
long dual_momentum_face_projection_fallback_count = 0;
long dual_momentum_face_projection_fallback_count_max = 0;
double dual_momentum_max_rejected_face_projection_defect = 0.;
double dual_momentum_max_rejected_face_projection_correction = 0.;
double dual_momentum_max_rejected_face_transport_correction = 0.;
long dual_momentum_nested_plic_fallback_count = 0;
long dual_momentum_nested_plic_fallback_count_max = 0;
long dual_momentum_compression_fallback_count = 0;
long dual_momentum_compression_fallback_count_max = 0;
static int aggregate_current_sweep_axis = -1;
static int aggregate_current_transport_iteration = -1;
static double aggregate_current_transport_time = -1.;

#if AGGREGATE_MATERIAL_CELL_DIAGNOSTICS
static FILE * aggregate_nested_plic_diagnostic_stream (void)
{
  static FILE * stream = NULL;
  static int attempted = 0;
  if (stream || attempted)
    return stream;
  attempted = 1;
  char path[128];
  sprintf(path, "data/nested_plic_cells.rank-%05d.csv", pid());
  stream = fopen(path, "w");
  if (!stream) {
    fprintf(stderr, "[MOMENTUM-DIAGNOSTIC-IO] pid=%d path=%s\n",
            pid(), path);
    fflush(stderr);
    return NULL;
  }
  fprintf(stream,
    "t,i,pid,sweep_axis,x,y,Delta,level,inner_fraction,outer_fraction,"
    "inner_normal_x,inner_normal_y,outer_normal_x,outer_normal_y,"
    "inner_norm_l1,outer_norm_l1,inner_raw_weight,outer_raw_weight,"
    "inner_weight,outer_weight,normal_dot,"
    "inner_alpha,outer_alpha,selected_normal_x,selected_normal_y,"
    "sliver_fallback,reason,"
    "fallback_neighbor,fallback_uses_outer,fallback_support_count,"
    "fallback_quality,"
    "xm_inner_fraction,xm_outer_fraction,xm_inner_normal_x,"
    "xm_inner_normal_y,xm_outer_normal_x,xm_outer_normal_y,"
    "xp_inner_fraction,xp_outer_fraction,xp_inner_normal_x,"
    "xp_inner_normal_y,xp_outer_normal_x,xp_outer_normal_y,"
    "ym_inner_fraction,ym_outer_fraction,ym_inner_normal_x,"
    "ym_inner_normal_y,ym_outer_normal_x,ym_outer_normal_y,"
    "yp_inner_fraction,yp_outer_fraction,yp_inner_normal_x,"
    "yp_inner_normal_y,yp_outer_normal_x,yp_outer_normal_y\n");
  fflush(stream);
  return stream;
}

static FILE * aggregate_face_reject_diagnostic_stream (void)
{
  static FILE * stream = NULL;
  static int attempted = 0;
  if (stream || attempted)
    return stream;
  attempted = 1;
  char path[160];
  sprintf(path, "data/momentum_face_rejects.rank-%05d.log", pid());
  stream = fopen(path, "w");
  if (!stream) {
    fprintf(stderr, "[MOMENTUM-FACE-DIAGNOSTIC-IO] pid=%d path=%s\n",
            pid(), path);
    fflush(stderr);
  }
  return stream;
}

static inline void aggregate_report_nested_plic_cell (
  Point point, scalar inner_fraction, scalar outer_fraction,
  vector inner_normal, vector outer_normal,
  double inner_alpha, double outer_alpha,
  NestedPlicCommonNormal common)
{
  FILE * stream = aggregate_nested_plic_diagnostic_stream();
  if (!stream)
    return;
  fprintf(stream,
    "%.17g,%d,%d,%d,%.17g,%.17g,%.17g,%d,"
    "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
    "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
    "%.17g,%.17g,%.17g,%.17g,%d,%d,%d,%d,%d,%.17g,"
    "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
    "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
    "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
    "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
    aggregate_current_transport_time,
    aggregate_current_transport_iteration, pid(),
    aggregate_current_sweep_axis, x, y, Delta, point.level,
    inner_fraction[], outer_fraction[],
    inner_normal.x[], inner_normal.y[],
    outer_normal.x[], outer_normal.y[],
    common.inner_norm_l1, common.outer_norm_l1,
    common.inner_raw_weight, common.outer_raw_weight,
    common.inner_weight, common.outer_weight, common.normal_dot,
    inner_alpha, outer_alpha, common.x, common.y,
    common.used_sliver_fallback, (int)common.reason,
    common.fallback_neighbor, common.fallback_uses_outer,
    common.fallback_support_count,
    common.fallback_quality,
    inner_fraction[-1,0], outer_fraction[-1,0],
    inner_normal.x[-1,0], inner_normal.y[-1,0],
    outer_normal.x[-1,0], outer_normal.y[-1,0],
    inner_fraction[1,0], outer_fraction[1,0],
    inner_normal.x[1,0], inner_normal.y[1,0],
    outer_normal.x[1,0], outer_normal.y[1,0],
    inner_fraction[0,-1], outer_fraction[0,-1],
    inner_normal.x[0,-1], inner_normal.y[0,-1],
    outer_normal.x[0,-1], outer_normal.y[0,-1],
    inner_fraction[0,1], outer_fraction[0,1],
    inner_normal.x[0,1], inner_normal.y[0,1],
    outer_normal.x[0,1], outer_normal.y[0,1]);
  fflush(stream);
}
#endif

static inline void aggregate_report_invalid_cell (
  const char * stage, PositiveMaterialPartition identity,
  double lower_gas, double lower_envelope,
  double upper_gas, double upper_envelope,
  double cell_x, double cell_y, double cell_delta)
{
#if AGGREGATE_MATERIAL_CELL_DIAGNOSTICS
  fprintf(stderr,
    "[MOMENTUM-CELL-REJECT] pid=%d stage=%s t=%.17g i=%d "
    "sweep_axis=%d x=%.17g y=%.17g Delta=%.17g "
    "lower_gas=%.17g lower_envelope=%.17g "
    "upper_gas=%.17g upper_envelope=%.17g "
    "identity_defect=%.17g bound_defect=%.17g correction_l1=%.17g\n",
    pid(), stage, aggregate_current_transport_time,
    aggregate_current_transport_iteration, aggregate_current_sweep_axis,
    cell_x, cell_y, cell_delta,
    lower_gas, lower_envelope, upper_gas, upper_envelope,
    identity.identity_defect, identity.bound_defect, identity.correction_l1);
  fflush(stderr);
#else
  (void)stage; (void)identity;
  (void)lower_gas; (void)lower_envelope;
  (void)upper_gas; (void)upper_envelope;
  (void)cell_x; (void)cell_y; (void)cell_delta;
#endif
}

static inline void aggregate_report_invalid_low_order_cell (
  const char * material, double initial, double low,
  double left_low_flux, double right_low_flux,
  double left_anti_flux, double right_anti_flux,
  double step, double divergence, double scale,
  double cell_x, double cell_y, double cell_delta)
{
#if AGGREGATE_MATERIAL_CELL_DIAGNOSTICS
  FILE * stream = aggregate_face_reject_diagnostic_stream();
  if (stream) {
    const double left_anti_contribution = scale*left_anti_flux;
    const double right_anti_contribution = -scale*right_anti_flux;
    fprintf(stream,
      "kind=low_order_cell pid=%d t=%.17g i=%d sweep_axis=%d "
      "x=%.17g y=%.17g Delta=%.17g material=%s "
      "initial=%.17g low=%.17g "
      "left_low_flux=%.17g right_low_flux=%.17g "
      "left_anti_flux=%.17g right_anti_flux=%.17g "
      "left_anti_contribution=%.17g right_anti_contribution=%.17g "
      "step=%.17g divergence=%.17g scale=%.17g\n",
      pid(), aggregate_current_transport_time,
      aggregate_current_transport_iteration, aggregate_current_sweep_axis,
      cell_x, cell_y, cell_delta, material,
      initial, low, left_low_flux, right_low_flux,
      left_anti_flux, right_anti_flux,
      left_anti_contribution, right_anti_contribution,
      step, divergence, scale);
    fflush(stream);
  }
#else
  (void)material; (void)initial; (void)low;
  (void)left_low_flux; (void)right_low_flux;
  (void)left_anti_flux; (void)right_anti_flux;
  (void)step; (void)divergence; (void)scale;
  (void)cell_x; (void)cell_y; (void)cell_delta;
#endif
}

static inline void aggregate_report_stage_reject (
  const char * stage, long invalid_cells, long invalid_nested_plic,
  double maximum_water_defect, double maximum_lower_oil_defect,
  double maximum_upper_oil_defect, double maximum_bound_defect)
{
#if AGGREGATE_MATERIAL_CELL_DIAGNOSTICS
  if (pid() == 0 && (invalid_cells > 0 || invalid_nested_plic > 0)) {
    fprintf(stderr,
      "[MOMENTUM-STAGE-REJECT] stage=%s t=%.17g i=%d sweep_axis=%d "
      "invalid_cells=%ld invalid_nested_plic=%ld "
      "water_defect=%.17g lower_oil_defect=%.17g "
      "upper_oil_defect=%.17g bound_defect=%.17g\n",
      stage, aggregate_current_transport_time,
      aggregate_current_transport_iteration, aggregate_current_sweep_axis,
      invalid_cells, invalid_nested_plic, maximum_water_defect,
      maximum_lower_oil_defect, maximum_upper_oil_defect,
      maximum_bound_defect);
    fflush(stderr);
  }
#else
  (void)stage; (void)invalid_cells; (void)invalid_nested_plic;
  (void)maximum_water_defect; (void)maximum_lower_oil_defect;
  (void)maximum_upper_oil_defect; (void)maximum_bound_defect;
#endif
}

attribute {
  int aggregate_material_id;
}

static inline PositiveMaterialPartition aggregate_identity_partition (
  double lower_gas, double lower_envelope,
  double upper_gas, double upper_envelope)
{
  return positive_material_partition_with_tolerances(
    lower_gas, lower_envelope, upper_gas, upper_envelope,
    (double)POSITIVE_MATERIAL_PARTITION_TOLERANCE,
    (double)POSITIVE_MATERIAL_VOF_BOUND_TOLERANCE);
}

static inline AggregatePhysicalPartition aggregate_transport_partition (
  double lower_gas, double lower_envelope,
  double upper_gas, double upper_envelope,
  PositiveMaterialPartition * identities)
{
  PositiveMaterialPartition identity = aggregate_identity_partition(
    lower_gas, lower_envelope, upper_gas, upper_envelope);
  if (identities)
    *identities = identity;
  return aggregate_physical_partition(&identity);
}

static inline double aggregate_transport_density (Point point)
{
  AggregatePhysicalPartition partition = aggregate_transport_partition(
    lower_gas_fraction[], lower_envelope_fraction[],
    upper_gas_fraction[], upper_envelope_fraction[], NULL);
  return aggregate_physical_density(
    &partition, rho_water, rho_oil, rho_gas);
}

static inline double aggregate_transport_material_density (int material)
{
  return material == AGGREGATE_MATERIAL_WATER ? rho_water :
    material == AGGREGATE_MATERIAL_OIL ? rho_oil : rho_gas;
}

#if TREE
static void aggregate_momentum_refine (Point point, scalar velocity)
{
  double parent_velocity = velocity[];
  double parent_metric = cm[];
  double parent_density = aggregate_transport_density(point);
  refine_bilinear(point, velocity);

  double child_momentum = 0., child_mass = 0.;
  foreach_child() {
    double density = aggregate_transport_density(point);
    child_momentum += cm[]*density*velocity[];
    child_mass += cm[]*density;
  }
  double target_momentum =
    (1 << dimension)*parent_metric*parent_density*parent_velocity;
  AggregateMassWeightedValue correction = aggregate_mass_weighted_value(
    target_momentum - child_momentum, child_mass);
  if (!correction.valid) {
    dual_momentum_invalid_cell_count++;
    dual_momentum_invalid_cell_count_max = max(
      dual_momentum_invalid_cell_count_max,
      dual_momentum_invalid_cell_count);
    dual_momentum_transport_failed = 1;
#if AGGREGATE_MATERIAL_CELL_DIAGNOSTICS
    fprintf(stderr,
      "[MOMENTUM-AMR-REJECT] pid=%d stage=refine "
      "x=%.17g y=%.17g Delta=%.17g target=%.17g "
      "child_momentum=%.17g child_mass=%.17g\n",
      pid(), x, y, Delta, target_momentum,
      child_momentum, child_mass);
    fflush(stderr);
#endif
  }
  foreach_child()
    velocity[] += correction.value;
}

static void aggregate_momentum_restriction (Point point, scalar velocity)
{
  double child_momentum = 0.;
  double child_mass = 0.;
  foreach_child() {
    child_mass += cm[]*aggregate_transport_density(point);
    child_momentum += cm[]*aggregate_transport_density(point)*velocity[];
  }
  AggregateMassWeightedValue restricted = aggregate_mass_weighted_value(
    child_momentum, child_mass);
  if (!restricted.valid) {
    dual_momentum_invalid_cell_count++;
    dual_momentum_invalid_cell_count_max = max(
      dual_momentum_invalid_cell_count_max,
      dual_momentum_invalid_cell_count);
    dual_momentum_transport_failed = 1;
#if AGGREGATE_MATERIAL_CELL_DIAGNOSTICS
    fprintf(stderr,
      "[MOMENTUM-AMR-REJECT] pid=%d stage=restriction "
      "x=%.17g y=%.17g Delta=%.17g "
      "child_momentum=%.17g child_mass=%.17g\n",
      pid(), x, y, Delta, child_momentum, child_mass);
    fflush(stderr);
#endif
  }
  velocity[] = restricted.value;
}

static void aggregate_move_to_front (scalar field)
{
  int position = 0;
  while (all[position].i != field.i)
    position++;
  while (position > 0 && all[position].i) {
    all[position] = all[position - 1];
    position--;
  }
  all[position] = field;
}
#endif

event defaults (i = 0)
{
  stokes = true;
#if TREE
  aggregate_move_to_front(upper_envelope_fraction);
  aggregate_move_to_front(upper_gas_fraction);
  aggregate_move_to_front(lower_envelope_fraction);
  aggregate_move_to_front(lower_gas_fraction);
  foreach_dimension() {
    u.x.refine = u.x.prolongation = aggregate_momentum_refine;
    u.x.restriction = aggregate_momentum_restriction;
    u.x.depends = list_add(u.x.depends, lower_gas_fraction);
    u.x.depends = list_add(u.x.depends, lower_envelope_fraction);
    u.x.depends = list_add(u.x.depends, upper_gas_fraction);
    u.x.depends = list_add(u.x.depends, upper_envelope_fraction);
  }
#endif
}

event stability (i++)
  dtmax = timestep(uf, dtmax);

foreach_dimension()
static double aggregate_boundary_momentum_x (
  Point neighbor, Point point, scalar momentum, bool * data)
{
  AggregatePhysicalPartition partition = aggregate_transport_partition(
    lower_gas_fraction[], lower_envelope_fraction[],
    upper_gas_fraction[], upper_envelope_fraction[], NULL);
  double fraction = aggregate_physical_fraction(
    &partition, momentum.aggregate_material_id);
  return fraction*aggregate_transport_material_density(
    momentum.aggregate_material_id)*u.x[];
}

#if TREE
foreach_dimension()
static void aggregate_prolong_momentum_x (Point point, scalar momentum)
{
  foreach_child() {
    AggregatePhysicalPartition partition = aggregate_transport_partition(
      lower_gas_fraction[], lower_envelope_fraction[],
      upper_gas_fraction[], upper_envelope_fraction[], NULL);
    double fraction = aggregate_physical_fraction(
      &partition, momentum.aggregate_material_id);
    momentum[] = fraction*aggregate_transport_material_density(
      momentum.aggregate_material_id)*u.x[];
  }
}
#endif

foreach_dimension()
static double aggregate_geometric_face_fraction_x (
  Point point, scalar indicator, vector normal, scalar intercept,
  int upwind, double direction, double courant)
{
  return (indicator[upwind] <= 0. || indicator[upwind] >= 1.) ?
    indicator[upwind] :
    rectangle_fraction(
      (coord){-direction*normal.x[upwind], normal.y[upwind]},
      intercept[upwind], (coord){-0.5, -0.5},
      (coord){direction*courant - 0.5, 0.5});
}

foreach_dimension()
static double aggregate_concentration_gradient_x (
  Point point, scalar momentum, int material)
{
  AggregatePhysicalPartition left = aggregate_transport_partition(
    lower_gas_fraction[-1], lower_envelope_fraction[-1],
    upper_gas_fraction[-1], upper_envelope_fraction[-1], NULL);
  AggregatePhysicalPartition center = aggregate_transport_partition(
    lower_gas_fraction[], lower_envelope_fraction[],
    upper_gas_fraction[], upper_envelope_fraction[], NULL);
  AggregatePhysicalPartition right = aggregate_transport_partition(
    lower_gas_fraction[1], lower_envelope_fraction[1],
    upper_gas_fraction[1], upper_envelope_fraction[1], NULL);
  if (!left.valid || !center.valid || !right.valid)
    return 0.;
  double cl = aggregate_physical_fraction(&left, material);
  double cc = aggregate_physical_fraction(&center, material);
  double cr = aggregate_physical_fraction(&right, material);
  if (cc >= 0.5 && momentum.gradient != zero) {
    if (cr >= 0.5) {
      if (cl >= 0.5) {
        if (momentum.gradient)
          return momentum.gradient(
            momentum[-1]/cl, momentum[]/cc, momentum[1]/cr)/Delta;
        return (momentum[1]/cr - momentum[-1]/cl)/(2.*Delta);
      }
      return (momentum[1]/cr - momentum[]/cc)/Delta;
    }
    if (cl >= 0.5)
      return (momentum[]/cc - momentum[-1]/cl)/Delta;
  }
  return 0.;
}

foreach_dimension()
static double aggregate_phase_momentum_face_flux_x (
  Point point, int material, scalar momentum, scalar gradient,
  int upwind, double direction, double courant,
  double face_fraction)
{
  AggregatePhysicalPartition upwind_partition =
    aggregate_transport_partition(
      lower_gas_fraction[upwind], lower_envelope_fraction[upwind],
      upper_gas_fraction[upwind], upper_envelope_fraction[upwind], NULL);
  if (!upwind_partition.valid)
    return 0.;
  double cell_fraction = aggregate_physical_fraction(
    &upwind_partition, material);
  if (cell_fraction <= 1.e-10)
    return 0.;
  double concentration = momentum[upwind]/cell_fraction +
    direction*min(1., 1. - direction*courant)*
    gradient[upwind]*Delta/2.;
  return concentration*face_fraction*uf.x[];
}

static long aggregate_reconcile_nested_plic (
  scalar inner_fraction, scalar outer_fraction,
  vector inner_normal, vector outer_normal,
  scalar inner_intercept, scalar outer_intercept)
{
  long invalid = 0;
  long fallback = 0;
  double invalid_inner_fraction = 0., invalid_outer_fraction = 0.;
  double invalid_inner_x = 0., invalid_inner_y = 0.;
  double invalid_outer_x = 0., invalid_outer_y = 0.;
  double invalid_x = 0., invalid_y = 0., invalid_delta = 0.;
  foreach(reduction(+:invalid) reduction(+:fallback)
          reduction(+:invalid_inner_fraction)
          reduction(+:invalid_outer_fraction)
          reduction(+:invalid_inner_x) reduction(+:invalid_inner_y)
          reduction(+:invalid_outer_x) reduction(+:invalid_outer_y)
          reduction(+:invalid_x) reduction(+:invalid_y)
          reduction(+:invalid_delta)) {
    bool inner_mixed = inner_fraction[] > 0. && inner_fraction[] < 1.;
    bool outer_mixed = outer_fraction[] > 0. && outer_fraction[] < 1.;
    if (inner_mixed && outer_mixed &&
        nested_plic_requires_common_normal(
          inner_fraction[], outer_fraction[]) &&
        !nested_plic_independent_cuts_are_nested(
          inner_normal.x[], inner_normal.y[], inner_intercept[],
          outer_normal.x[], outer_normal.y[], outer_intercept[])) {
      const NestedPlicNeighbor neighbors[NESTED_PLIC_NEIGHBOR_COUNT] = {
        {
          inner_fraction[-1,0], outer_fraction[-1,0],
          inner_normal.x[-1,0], inner_normal.y[-1,0],
          outer_normal.x[-1,0], outer_normal.y[-1,0]
        },
        {
          inner_fraction[1,0], outer_fraction[1,0],
          inner_normal.x[1,0], inner_normal.y[1,0],
          outer_normal.x[1,0], outer_normal.y[1,0]
        },
        {
          inner_fraction[0,-1], outer_fraction[0,-1],
          inner_normal.x[0,-1], inner_normal.y[0,-1],
          outer_normal.x[0,-1], outer_normal.y[0,-1]
        },
        {
          inner_fraction[0,1], outer_fraction[0,1],
          inner_normal.x[0,1], inner_normal.y[0,1],
          outer_normal.x[0,1], outer_normal.y[0,1]
        }
      };
      NestedPlicCommonNormal common =
        nested_plic_common_normal_with_neighbors(
        inner_fraction[], outer_fraction[],
        inner_normal.x[], inner_normal.y[],
        outer_normal.x[], outer_normal.y[], neighbors);
#if AGGREGATE_MATERIAL_CELL_DIAGNOSTICS && \
    AGGREGATE_MATERIAL_FALLBACK_CELL_DIAGNOSTICS
      bool consensus_repair =
        (common.reason == NESTED_PLIC_NEIGHBOR_FALLBACK ||
         common.reason == NESTED_PLIC_PAIRED_NEIGHBOR_CONSENSUS) &&
        common.fallback_support_count >=
          (int)NESTED_PLIC_NEIGHBOR_CONSENSUS_COUNT_MIN;
      if (common.used_sliver_fallback &&
          (AGGREGATE_MATERIAL_FALLBACK_CELL_DIAGNOSTICS >= 2 ||
           consensus_repair))
        aggregate_report_nested_plic_cell(
          point, inner_fraction, outer_fraction,
          inner_normal, outer_normal,
          inner_intercept[], outer_intercept[], common);
#endif
      if (!common.valid) {
#if AGGREGATE_MATERIAL_CELL_DIAGNOSTICS
        aggregate_report_nested_plic_cell(
          point, inner_fraction, outer_fraction,
          inner_normal, outer_normal,
          inner_intercept[], outer_intercept[], common);
#endif
        invalid_inner_fraction += inner_fraction[];
        invalid_outer_fraction += outer_fraction[];
        invalid_inner_x += inner_normal.x[];
        invalid_inner_y += inner_normal.y[];
        invalid_outer_x += outer_normal.x[];
        invalid_outer_y += outer_normal.y[];
        invalid_x += x;
        invalid_y += y;
        invalid_delta += Delta;
        invalid++;
        continue;
      }
      fallback += common.used_sliver_fallback;
      coord normal = {common.x, common.y};
      inner_normal.x[] = outer_normal.x[] = common.x;
      inner_normal.y[] = outer_normal.y[] = common.y;
      inner_intercept[] = plane_alpha(inner_fraction[], normal);
      outer_intercept[] = plane_alpha(outer_fraction[], normal);
    }
  }
  boundary((scalar *){
    inner_normal, outer_normal, inner_intercept, outer_intercept
  });
#if AGGREGATE_MATERIAL_CELL_DIAGNOSTICS
  if (pid() == 0 && invalid > 0) {
    double count = invalid;
    fprintf(stderr,
      "[MOMENTUM-NESTED-PLIC-CELL] count=%ld t=%.17g i=%d "
      "sweep_axis=%d mean_x=%.17g mean_y=%.17g mean_Delta=%.17g "
      "mean_inner_fraction=%.17g mean_outer_fraction=%.17g "
      "mean_inner_normal=(%.17g,%.17g) "
      "mean_outer_normal=(%.17g,%.17g)\n",
      invalid, aggregate_current_transport_time,
      aggregate_current_transport_iteration,
      aggregate_current_sweep_axis, invalid_x/count, invalid_y/count,
      invalid_delta/count, invalid_inner_fraction/count,
      invalid_outer_fraction/count, invalid_inner_x/count,
      invalid_inner_y/count, invalid_outer_x/count, invalid_outer_y/count);
    fflush(stderr);
  }
#endif
  dual_momentum_nested_plic_fallback_count += fallback;
  dual_momentum_nested_plic_fallback_count_max = max(
    dual_momentum_nested_plic_fallback_count_max,
    dual_momentum_nested_plic_fallback_count);
  return invalid;
}

foreach_dimension()
static void aggregate_momentum_sweep_x (
  vector water_momentum, vector oil_momentum, vector gas_momentum,
  scalar lower_gas_step, scalar lower_envelope_step,
  scalar upper_gas_step, scalar upper_envelope_step,
  vector water_compression, vector oil_compression, vector gas_compression)
{
  vector lower_gas_normal[], lower_envelope_normal[];
  vector upper_gas_normal[], upper_envelope_normal[];
  scalar lower_gas_intercept[], lower_envelope_intercept[];
  scalar upper_gas_intercept[], upper_envelope_intercept[];
  scalar lower_gas_flux[], lower_envelope_flux[];
  scalar upper_gas_flux[], upper_envelope_flux[];
  vector water_gradient[], oil_gradient[], gas_gradient[];
  scalar water_first_flux[], water_second_flux[];
  scalar oil_first_flux[], oil_second_flux[];
  scalar gas_first_flux[], gas_second_flux[];
#if AGGREGATE_MATERIAL_SIMPLEX_FLUX_LIMITING
  scalar material_water[], material_lower_oil[], material_lower_gas[];
  scalar material_upper_oil[], material_upper_gas[];
  scalar step_water[], step_lower_oil[], step_lower_gas[];
  scalar step_upper_oil[], step_upper_gas[];
  scalar low_water_flux[], low_lower_oil_flux[], low_lower_gas_flux[];
  scalar low_upper_oil_flux[], low_upper_gas_flux[];
  scalar anti_water_flux[], anti_lower_oil_flux[], anti_lower_gas_flux[];
  scalar anti_upper_oil_flux[], anti_upper_gas_flux[];
  scalar ratio_water[], ratio_lower_oil[], ratio_lower_gas[];
  scalar ratio_upper_oil[], ratio_upper_gas[];
#endif

  reconstruction(lower_gas_fraction,
                 lower_gas_normal, lower_gas_intercept);
  reconstruction(lower_envelope_fraction,
                 lower_envelope_normal, lower_envelope_intercept);
  reconstruction(upper_gas_fraction,
                 upper_gas_normal, upper_gas_intercept);
  reconstruction(upper_envelope_fraction,
                 upper_envelope_normal, upper_envelope_intercept);

  long invalid_nested_plic = 0;
  invalid_nested_plic += aggregate_reconcile_nested_plic(
    lower_gas_fraction, lower_envelope_fraction,
    lower_gas_normal, lower_envelope_normal,
    lower_gas_intercept, lower_envelope_intercept);
  invalid_nested_plic += aggregate_reconcile_nested_plic(
    upper_gas_fraction, upper_envelope_fraction,
    upper_gas_normal, upper_envelope_normal,
    upper_gas_intercept, upper_envelope_intercept);
  if (invalid_nested_plic > 0) {
    if (pid() == 0)
      fprintf(stderr,
        "[MOMENTUM-NESTED-PLIC-REJECT] t=%.17g i=%d sweep_axis=%d "
        "invalid_cells=%ld\n",
        aggregate_current_transport_time,
        aggregate_current_transport_iteration,
        aggregate_current_sweep_axis, invalid_nested_plic);
    dual_momentum_transport_failed = 1;
  }

  long invalid_cells = 0;
  double maximum_cell_defect = 0., maximum_bound_defect = 0.;
  double maximum_correction = 0.;
  double maximum_water_defect = 0., maximum_lower_oil_defect = 0.;
  double maximum_upper_oil_defect = 0.;
  foreach(reduction(+:invalid_cells) reduction(max:maximum_cell_defect)
          reduction(max:maximum_bound_defect)
          reduction(max:maximum_correction)
          reduction(max:maximum_water_defect)
          reduction(max:maximum_lower_oil_defect)
          reduction(max:maximum_upper_oil_defect)) {
    PositiveMaterialPartition identity;
    AggregatePhysicalPartition partition = aggregate_transport_partition(
      lower_gas_fraction[], lower_envelope_fraction[],
      upper_gas_fraction[], upper_envelope_fraction[], &identity);
    maximum_cell_defect = max(
      maximum_cell_defect, identity.identity_defect);
    maximum_bound_defect = max(
      maximum_bound_defect, identity.bound_defect);
    maximum_correction = max(maximum_correction, identity.correction_l1);
    maximum_water_defect = max(
      maximum_water_defect, identity.water_defect);
    maximum_lower_oil_defect = max(
      maximum_lower_oil_defect, identity.lower_oil_defect);
    maximum_upper_oil_defect = max(
      maximum_upper_oil_defect, identity.upper_oil_defect);
    if (!partition.valid) {
      aggregate_report_invalid_cell(
        "pre-sweep", identity,
        lower_gas_fraction[], lower_envelope_fraction[],
        upper_gas_fraction[], upper_envelope_fraction[], x, y, Delta);
      invalid_cells++;
    }
#if AGGREGATE_MATERIAL_SIMPLEX_FLUX_LIMITING
    material_water[] = identity.water;
    material_lower_oil[] = identity.lower_oil;
    material_lower_gas[] = identity.lower_gas;
    material_upper_oil[] = identity.upper_oil;
    material_upper_gas[] = identity.upper_gas;
    step_water[] = 1. - lower_envelope_step[] - upper_envelope_step[];
    step_lower_oil[] = lower_envelope_step[] - lower_gas_step[];
    step_lower_gas[] = lower_gas_step[];
    step_upper_oil[] = upper_envelope_step[] - upper_gas_step[];
    step_upper_gas[] = upper_gas_step[];
#endif
  }
  aggregate_report_stage_reject(
    "pre-sweep", invalid_cells, invalid_nested_plic,
    maximum_water_defect, maximum_lower_oil_defect,
    maximum_upper_oil_defect, maximum_bound_defect);
  invalid_cells += invalid_nested_plic;

  foreach() {
    water_gradient.x[] = aggregate_concentration_gradient_x(
      point, water_momentum.x, AGGREGATE_MATERIAL_WATER);
    water_gradient.y[] = aggregate_concentration_gradient_x(
      point, water_momentum.y, AGGREGATE_MATERIAL_WATER);
    oil_gradient.x[] = aggregate_concentration_gradient_x(
      point, oil_momentum.x, AGGREGATE_MATERIAL_OIL);
    oil_gradient.y[] = aggregate_concentration_gradient_x(
      point, oil_momentum.y, AGGREGATE_MATERIAL_OIL);
    gas_gradient.x[] = aggregate_concentration_gradient_x(
      point, gas_momentum.x, AGGREGATE_MATERIAL_GAS);
    gas_gradient.y[] = aggregate_concentration_gradient_x(
      point, gas_momentum.y, AGGREGATE_MATERIAL_GAS);
  }

  double maximum_face_defect = 0.;
  long invalid_faces = 0;
  long face_projection_fallbacks = 0;
  double maximum_cfl = 0.;
  double maximum_rejected_face_defect = 0.;
  double maximum_rejected_face_correction = 0.;
  double maximum_rejected_face_transport_correction = 0.;
  foreach_face(x, reduction(max:maximum_face_defect)
                   reduction(+:invalid_faces) reduction(max:maximum_cfl)
                   reduction(+:face_projection_fallbacks)
                   reduction(max:maximum_bound_defect)
                   reduction(max:maximum_correction)
                   reduction(max:maximum_rejected_face_defect)
                   reduction(max:maximum_rejected_face_correction)
                   reduction(max:maximum_rejected_face_transport_correction)) {
    double courant = uf.x[]*dt/(Delta*fm.x[] + SEPS);
    double direction = sign(courant);
    int upwind = -(direction + 1.)/2.;
    maximum_cfl = max(maximum_cfl,
                      courant*fm.x[]*direction/(cm[] + SEPS));

    double raw_lower_gas_face = aggregate_geometric_face_fraction_x(
      point, lower_gas_fraction, lower_gas_normal, lower_gas_intercept,
      upwind, direction, courant);
    double raw_lower_envelope_face = aggregate_geometric_face_fraction_x(
      point, lower_envelope_fraction, lower_envelope_normal,
      lower_envelope_intercept, upwind, direction, courant);
    double raw_upper_gas_face = aggregate_geometric_face_fraction_x(
      point, upper_gas_fraction, upper_gas_normal, upper_gas_intercept,
      upwind, direction, courant);
    double raw_upper_envelope_face = aggregate_geometric_face_fraction_x(
      point, upper_envelope_fraction, upper_envelope_normal,
      upper_envelope_intercept, upwind, direction, courant);

    NestedFaceProjection projection = nested_face_projection(
      raw_lower_gas_face, raw_lower_envelope_face,
      raw_upper_gas_face, raw_upper_envelope_face);
    bool projection_admissible = projection.valid &&
      projection.raw_bound_defect <=
        (double)POSITIVE_MATERIAL_VOF_BOUND_TOLERANCE &&
      projection.correction_linf <=
        (double)AGGREGATE_MATERIAL_FACE_PROJECTION_LIMIT;
    PositiveMaterialPartition used_identity = projection.materials;
    AggregatePhysicalPartition face =
      aggregate_physical_partition(&used_identity);
    bool hard_face_reject = !projection_admissible || !face.valid;
#if AGGREGATE_MATERIAL_SIMPLEX_FLUX_LIMITING
    PositiveMaterialPartition upwind_identity = aggregate_identity_partition(
      lower_gas_fraction[upwind], lower_envelope_fraction[upwind],
      upper_gas_fraction[upwind], upper_envelope_fraction[upwind]);
    int upwind_partition_valid = upwind_identity.valid;
    if (!upwind_partition_valid) {
      invalid_faces++;
#if AGGREGATE_MATERIAL_CELL_DIAGNOSTICS
      FILE * reject_stream = aggregate_face_reject_diagnostic_stream();
      if (reject_stream) {
        fprintf(reject_stream,
          "kind=upwind_identity pid=%d t=%.17g i=%d sweep_axis=%d "
          "x=%.17g y=%.17g Delta=%.17g courant=%.17g upwind=%d "
          "lower_gas=%.17g lower_envelope=%.17g upper_gas=%.17g "
          "upper_envelope=%.17g bound_defect=%.17g "
          "identity_defect=%.17g water_defect=%.17g "
          "lower_oil_defect=%.17g upper_oil_defect=%.17g "
          "partition_defect=%.17g correction_l1=%.17g\n",
          pid(), aggregate_current_transport_time,
          aggregate_current_transport_iteration,
          aggregate_current_sweep_axis,
          x, y, Delta, courant, upwind,
          lower_gas_fraction[upwind], lower_envelope_fraction[upwind],
          upper_gas_fraction[upwind], upper_envelope_fraction[upwind],
          upwind_identity.bound_defect, upwind_identity.identity_defect,
          upwind_identity.water_defect, upwind_identity.lower_oil_defect,
          upwind_identity.upper_oil_defect, upwind_identity.defect,
          upwind_identity.correction_l1);
        fflush(reject_stream);
      }
#endif
      fprintf(stderr,
        "[MOMENTUM-UPWIND-IDENTITY-REJECT] pid=%d t=%.17g i=%d "
        "sweep_axis=%d x=%.17g y=%.17g Delta=%.17g "
        "courant=%.17g upwind=%d lower_gas=%.17g "
        "lower_envelope=%.17g upper_gas=%.17g "
        "upper_envelope=%.17g bound_defect=%.17g "
        "identity_defect=%.17g water_defect=%.17g "
        "lower_oil_defect=%.17g upper_oil_defect=%.17g "
        "partition_defect=%.17g correction_l1=%.17g\n",
        pid(), aggregate_current_transport_time,
        aggregate_current_transport_iteration,
        aggregate_current_sweep_axis,
        x, y, Delta, courant, upwind,
        lower_gas_fraction[upwind], lower_envelope_fraction[upwind],
        upper_gas_fraction[upwind], upper_envelope_fraction[upwind],
        upwind_identity.bound_defect, upwind_identity.identity_defect,
        upwind_identity.water_defect, upwind_identity.lower_oil_defect,
        upwind_identity.upper_oil_defect, upwind_identity.defect,
        upwind_identity.correction_l1);
      fflush(stderr);
      upwind_identity = aggregate_identity_partition(0., 0., 0., 0.);
    }
    AggregateFaceProjectionDecision projection_decision =
      aggregate_face_projection_decision(
        &projection, upwind_partition_valid, face.valid, courant,
        (double)POSITIVE_MATERIAL_VOF_BOUND_TOLERANCE,
        (double)AGGREGATE_MATERIAL_FACE_PROJECTION_LIMIT);
    if (projection_decision.used_low_order_fallback) {
      face_projection_fallbacks++;
      maximum_rejected_face_defect = max(
        maximum_rejected_face_defect, projection.raw_identity_defect);
      maximum_rejected_face_correction = max(
        maximum_rejected_face_correction, projection.correction_linf);
      maximum_rejected_face_transport_correction = max(
        maximum_rejected_face_transport_correction,
        projection_decision.transported_correction);
#if AGGREGATE_MATERIAL_CELL_DIAGNOSTICS
      FILE * fallback_stream = aggregate_face_reject_diagnostic_stream();
      if (fallback_stream) {
        fprintf(fallback_stream,
          "kind=face_projection_fallback pid=%d t=%.17g i=%d "
          "sweep_axis=%d x=%.17g y=%.17g Delta=%.17g "
          "courant=%.17g upwind=%d identity_defect=%.17g "
          "correction_linf=%.17g transported_correction=%.17g\n",
          pid(), aggregate_current_transport_time,
          aggregate_current_transport_iteration,
          aggregate_current_sweep_axis,
          x, y, Delta, courant, upwind,
          projection.raw_identity_defect, projection.correction_linf,
          projection_decision.transported_correction);
        fflush(fallback_stream);
      }
#endif
    }
    PositiveMaterialPartition high_identity =
      projection_decision.use_projected_candidate ?
      projection.materials : upwind_identity;
    used_identity = high_identity;
    face = aggregate_physical_partition(&used_identity);
    hard_face_reject = !projection_decision.valid || !face.valid;
#define AGGREGATE_MATERIAL_FLUXES(name, member) do {                    \
      low_##name##_flux[] = upwind_identity.member*uf.x[];              \
      anti_##name##_flux[] =                                            \
        (high_identity.member - upwind_identity.member)*uf.x[];         \
    } while (0)
    AGGREGATE_MATERIAL_FLUXES(water, water);
    AGGREGATE_MATERIAL_FLUXES(lower_oil, lower_oil);
    AGGREGATE_MATERIAL_FLUXES(lower_gas, lower_gas);
    AGGREGATE_MATERIAL_FLUXES(upper_oil, upper_oil);
    AGGREGATE_MATERIAL_FLUXES(upper_gas, upper_gas);
#undef AGGREGATE_MATERIAL_FLUXES
#else
    /* Version 3 preserves the four independently reconstructed VOF fluxes. */
    lower_gas_flux[] = projection.vof_lower_gas*uf.x[];
    lower_envelope_flux[] = projection.vof_lower_envelope*uf.x[];
    upper_gas_flux[] = projection.vof_upper_gas*uf.x[];
    upper_envelope_flux[] = projection.vof_upper_envelope*uf.x[];
#endif

    maximum_face_defect = max(
      maximum_face_defect,
      hard_face_reject ? projection.raw_identity_defect :
      used_identity.identity_defect);
    maximum_bound_defect = max(
      maximum_bound_defect, projection.raw_bound_defect);
    maximum_correction = max(
      maximum_correction,
      hard_face_reject ? projection.materials.correction_l1 :
      used_identity.correction_l1);
    if (hard_face_reject) {
      invalid_faces++;
#if AGGREGATE_MATERIAL_CELL_DIAGNOSTICS
      FILE * reject_stream = aggregate_face_reject_diagnostic_stream();
      if (reject_stream) {
        fprintf(reject_stream,
          "kind=face_projection pid=%d t=%.17g i=%d sweep_axis=%d "
          "x=%.17g y=%.17g Delta=%.17g courant=%.17g upwind=%d "
          "lower_gas=%.17g lower_envelope=%.17g upper_gas=%.17g "
          "upper_envelope=%.17g identity_defect=%.17g "
          "bound_defect=%.17g correction_l1=%.17g "
          "correction_linf=%.17g limit=%.17g\n",
          pid(), aggregate_current_transport_time,
          aggregate_current_transport_iteration,
          aggregate_current_sweep_axis,
          x, y, Delta, courant, upwind,
          raw_lower_gas_face, raw_lower_envelope_face,
          raw_upper_gas_face, raw_upper_envelope_face,
          projection.raw_identity_defect, projection.raw_bound_defect,
          projection.materials.correction_l1,
          projection.correction_linf,
          (double)AGGREGATE_MATERIAL_FACE_PROJECTION_LIMIT);
        fflush(reject_stream);
      }
#endif
      fprintf(stderr,
        "[MOMENTUM-FACE-REJECT] pid=%d t=%.17g i=%d sweep_axis=%d "
        "x=%.17g y=%.17g Delta=%.17g courant=%.17g "
        "lower_gas=%.17g lower_envelope=%.17g "
        "upper_gas=%.17g upper_envelope=%.17g "
        "identity_defect=%.17g bound_defect=%.17g "
        "correction_l1=%.17g correction_linf=%.17g limit=%.17g\n",
        pid(), aggregate_current_transport_time,
        aggregate_current_transport_iteration,
        aggregate_current_sweep_axis,
        x, y, Delta, courant,
        raw_lower_gas_face, raw_lower_envelope_face,
        raw_upper_gas_face, raw_upper_envelope_face,
        projection.raw_identity_defect, projection.raw_bound_defect,
        projection.materials.correction_l1,
        projection.correction_linf,
        (double)AGGREGATE_MATERIAL_FACE_PROJECTION_LIMIT);
      fflush(stderr);
      PositiveMaterialPartition water = aggregate_identity_partition(
        0., 0., 0., 0.);
      face = aggregate_physical_partition(&water);
    }

#if !AGGREGATE_MATERIAL_SIMPLEX_FLUX_LIMITING
#define AGGREGATE_Q_FLUX(material, momentum, gradient, first, second) do { \
      (first)[] = aggregate_phase_momentum_face_flux_x(                 \
        point, (material), (momentum).x, (gradient).x,                  \
        upwind, direction, courant,                                    \
        aggregate_physical_fraction(&face, (material)));               \
      (second)[] = aggregate_phase_momentum_face_flux_x(                \
        point, (material), (momentum).y, (gradient).y,                  \
        upwind, direction, courant,                                    \
        aggregate_physical_fraction(&face, (material)));               \
    } while (0)
    AGGREGATE_Q_FLUX(AGGREGATE_MATERIAL_WATER,
      water_momentum, water_gradient, water_first_flux, water_second_flux);
    AGGREGATE_Q_FLUX(AGGREGATE_MATERIAL_OIL,
      oil_momentum, oil_gradient, oil_first_flux, oil_second_flux);
    AGGREGATE_Q_FLUX(AGGREGATE_MATERIAL_GAS,
      gas_momentum, gas_gradient, gas_first_flux, gas_second_flux);
#undef AGGREGATE_Q_FLUX
#endif
  }

#if AGGREGATE_MATERIAL_SIMPLEX_FLUX_LIMITING
  long invalid_low_order_cells = 0, compression_partition_fallbacks = 0;
  foreach(reduction(+:invalid_low_order_cells)
          reduction(+:compression_partition_fallbacks)) {
    double divergence = uf.x[1] - uf.x[];
    double scale = dt/(cm[]*Delta + SEPS);
    int compression_fallback_needed = 0;
#define AGGREGATE_LOW_ORDER_BUDGET(name) do {                           \
      double low = material_##name[] + scale*(                         \
        low_##name##_flux[] - low_##name##_flux[1] +                   \
        step_##name[]*divergence);                                     \
      double left_anti = scale*anti_##name##_flux[];                   \
      double right_anti = -scale*anti_##name##_flux[1];                \
      double negative_anti = min(0., left_anti) + min(0., right_anti); \
      if (!isfinite(low) || low <                                      \
          -(double)POSITIVE_MATERIAL_PARTITION_TOLERANCE) {            \
        compression_fallback_needed = 1;                               \
        ratio_##name[] = 0.;                                          \
      }                                                               \
      else                                                             \
        ratio_##name[] = negative_anti < 0. ?                          \
          clamp(max(0., low)/(-negative_anti + SEPS), 0., 1.) : 1.;   \
    } while (0)
    AGGREGATE_LOW_ORDER_BUDGET(water);
    AGGREGATE_LOW_ORDER_BUDGET(lower_oil);
    AGGREGATE_LOW_ORDER_BUDGET(lower_gas);
    AGGREGATE_LOW_ORDER_BUDGET(upper_oil);
    AGGREGATE_LOW_ORDER_BUDGET(upper_gas);
#undef AGGREGATE_LOW_ORDER_BUDGET

    if (compression_fallback_needed) {
      double initial[SIMPLEX_FLUX_MATERIALS] = {
        material_water[], material_lower_oil[], material_lower_gas[],
        material_upper_oil[], material_upper_gas[]
      };
      double left_low[SIMPLEX_FLUX_MATERIALS] = {
        low_water_flux[], low_lower_oil_flux[], low_lower_gas_flux[],
        low_upper_oil_flux[], low_upper_gas_flux[]
      };
      double right_low[SIMPLEX_FLUX_MATERIALS] = {
        low_water_flux[1], low_lower_oil_flux[1], low_lower_gas_flux[1],
        low_upper_oil_flux[1], low_upper_gas_flux[1]
      };
      double left_anti[SIMPLEX_FLUX_MATERIALS] = {
        anti_water_flux[], anti_lower_oil_flux[], anti_lower_gas_flux[],
        anti_upper_oil_flux[], anti_upper_gas_flux[]
      };
      double right_anti[SIMPLEX_FLUX_MATERIALS] = {
        anti_water_flux[1], anti_lower_oil_flux[1], anti_lower_gas_flux[1],
        anti_upper_oil_flux[1], anti_upper_gas_flux[1]
      };
      double requested[SIMPLEX_FLUX_MATERIALS] = {
        step_water[], step_lower_oil[], step_lower_gas[],
        step_upper_oil[], step_upper_gas[]
      };
      int fallback_material = 0;
      for (int material = 1; material < SIMPLEX_FLUX_MATERIALS; material++)
        if (initial[material] > initial[fallback_material])
          fallback_material = material;
      double fallback[SIMPLEX_FLUX_MATERIALS] = {0., 0., 0., 0., 0.};
      fallback[fallback_material] = 1.;
      SimplexLowOrderBudgetSet budgets = simplex_low_order_budget_set(
        initial, left_low, right_low, left_anti, right_anti,
        requested, fallback, divergence, scale,
        (double)POSITIVE_MATERIAL_PARTITION_TOLERANCE);
      if (budgets.valid && budgets.used_fallback) {
        compression_partition_fallbacks++;
        step_water[] = fallback[0];
        step_lower_oil[] = fallback[1];
        step_lower_gas[] = fallback[2];
        step_upper_oil[] = fallback[3];
        step_upper_gas[] = fallback[4];
        lower_gas_step[] = fallback[2];
        lower_envelope_step[] = fallback[1] + fallback[2];
        upper_gas_step[] = fallback[4];
        upper_envelope_step[] = fallback[3] + fallback[4];
#if AGGREGATE_MATERIAL_CELL_DIAGNOSTICS
        FILE * fallback_stream = aggregate_face_reject_diagnostic_stream();
        if (fallback_stream) {
          fprintf(fallback_stream,
            "kind=compression_partition_fallback pid=%d t=%.17g i=%d "
            "sweep_axis=%d x=%.17g y=%.17g Delta=%.17g "
            "requested=%.17g,%.17g,%.17g,%.17g,%.17g "
            "fallback_material=%d\n",
            pid(), aggregate_current_transport_time,
            aggregate_current_transport_iteration,
            aggregate_current_sweep_axis, x, y, Delta,
            requested[0], requested[1], requested[2], requested[3],
            requested[4], fallback_material);
          fflush(fallback_stream);
        }
#endif
        ratio_water[] = budgets.material[0].depletion_ratio;
        ratio_lower_oil[] = budgets.material[1].depletion_ratio;
        ratio_lower_gas[] = budgets.material[2].depletion_ratio;
        ratio_upper_oil[] = budgets.material[3].depletion_ratio;
        ratio_upper_gas[] = budgets.material[4].depletion_ratio;
      }
      else {
        invalid_low_order_cells++;
        ratio_water[] = ratio_lower_oil[] = ratio_lower_gas[] = 0.;
        ratio_upper_oil[] = ratio_upper_gas[] = 0.;
        for (int material = 0; material < SIMPLEX_FLUX_MATERIALS;
             material++)
          if (!budgets.material[material].valid)
            aggregate_report_invalid_low_order_cell(
              material == 0 ? "water" :
              material == 1 ? "lower_oil" :
              material == 2 ? "lower_gas" :
              material == 3 ? "upper_oil" : "upper_gas",
              initial[material], budgets.material[material].value,
              left_low[material], right_low[material],
              left_anti[material], right_anti[material],
              fallback[material], divergence, scale, x, y, Delta);
      }
    }
  }
  boundary((scalar *){
    ratio_water, ratio_lower_oil, ratio_lower_gas,
    ratio_upper_oil, ratio_upper_gas
  });

  long invalid_limiter_faces = 0, limited_faces = 0;
  double minimum_flux_alpha = 1., maximum_simplex_closure_defect = 0.;
  double maximum_flux_fraction_correction = 0.;
  double maximum_flux_courant_correction = 0.;
  foreach_face(x, reduction(+:invalid_limiter_faces)
                  reduction(+:limited_faces)
                  reduction(min:minimum_flux_alpha)
                  reduction(max:maximum_simplex_closure_defect)
                  reduction(max:maximum_flux_fraction_correction)
                  reduction(max:maximum_flux_courant_correction)) {
    double high[SIMPLEX_FLUX_MATERIALS], low[SIMPLEX_FLUX_MATERIALS];
    double left_ratio[SIMPLEX_FLUX_MATERIALS] = {
      ratio_water[-1], ratio_lower_oil[-1], ratio_lower_gas[-1],
      ratio_upper_oil[-1], ratio_upper_gas[-1]
    };
    double right_ratio[SIMPLEX_FLUX_MATERIALS] = {
      ratio_water[], ratio_lower_oil[], ratio_lower_gas[],
      ratio_upper_oil[], ratio_upper_gas[]
    };
    double low_flux[SIMPLEX_FLUX_MATERIALS] = {
      low_water_flux[], low_lower_oil_flux[], low_lower_gas_flux[],
      low_upper_oil_flux[], low_upper_gas_flux[]
    };
    double anti_flux[SIMPLEX_FLUX_MATERIALS] = {
      anti_water_flux[], anti_lower_oil_flux[], anti_lower_gas_flux[],
      anti_upper_oil_flux[], anti_upper_gas_flux[]
    };
    if (fabs(uf.x[]) > SEPS)
      for (int material = 0; material < SIMPLEX_FLUX_MATERIALS;
           material++) {
        low[material] = low_flux[material]/uf.x[];
        high[material] = (low_flux[material] + anti_flux[material])/uf.x[];
      }
    else {
      PositiveMaterialPartition zero_velocity = aggregate_identity_partition(
        lower_gas_fraction[], lower_envelope_fraction[],
        upper_gas_fraction[], upper_envelope_fraction[]);
      low[0] = high[0] = zero_velocity.water;
      low[1] = high[1] = zero_velocity.lower_oil;
      low[2] = high[2] = zero_velocity.lower_gas;
      low[3] = high[3] = zero_velocity.upper_oil;
      low[4] = high[4] = zero_velocity.upper_gas;
    }
    double courant = uf.x[]*dt/(Delta*fm.x[] + SEPS);
    SimplexFaceFlux limited = simplex_face_flux_limit(
      high, low, uf.x[], left_ratio, right_ratio,
      (double)POSITIVE_MATERIAL_PARTITION_TOLERANCE);
    if (!limited.valid) {
      invalid_limiter_faces++;
#if AGGREGATE_MATERIAL_CELL_DIAGNOSTICS
      FILE * reject_stream = aggregate_face_reject_diagnostic_stream();
      if (reject_stream) {
        fprintf(reject_stream,
          "kind=simplex_limiter pid=%d t=%.17g i=%d sweep_axis=%d "
          "x=%.17g y=%.17g Delta=%.17g uf=%.17g courant=%.17g "
          "reason=%d material=%d "
          "high=%.17g,%.17g,%.17g,%.17g,%.17g "
          "low=%.17g,%.17g,%.17g,%.17g,%.17g "
          "left=%.17g,%.17g,%.17g,%.17g,%.17g "
          "right=%.17g,%.17g,%.17g,%.17g,%.17g\n",
          pid(), aggregate_current_transport_time,
          aggregate_current_transport_iteration,
          aggregate_current_sweep_axis,
          x, y, Delta, uf.x[], courant,
          (int)limited.failure_reason, limited.failure_material,
          high[0], high[1], high[2], high[3], high[4],
          low[0], low[1], low[2], low[3], low[4],
          left_ratio[0], left_ratio[1], left_ratio[2],
          left_ratio[3], left_ratio[4],
          right_ratio[0], right_ratio[1], right_ratio[2],
          right_ratio[3], right_ratio[4]);
        fflush(reject_stream);
      }
#endif
      fprintf(stderr,
        "[MOMENTUM-LIMITER-REJECT] pid=%d t=%.17g i=%d sweep_axis=%d "
        "x=%.17g y=%.17g Delta=%.17g uf=%.17g courant=%.17g "
        "reason=%d material=%d "
        "high=%.17g,%.17g,%.17g,%.17g,%.17g "
        "low=%.17g,%.17g,%.17g,%.17g,%.17g "
        "left=%.17g,%.17g,%.17g,%.17g,%.17g "
        "right=%.17g,%.17g,%.17g,%.17g,%.17g\n",
        pid(), aggregate_current_transport_time,
        aggregate_current_transport_iteration,
        aggregate_current_sweep_axis,
        x, y, Delta, uf.x[], courant,
        (int)limited.failure_reason, limited.failure_material,
        high[0], high[1], high[2], high[3], high[4],
        low[0], low[1], low[2], low[3], low[4],
        left_ratio[0], left_ratio[1], left_ratio[2],
        left_ratio[3], left_ratio[4],
        right_ratio[0], right_ratio[1], right_ratio[2],
        right_ratio[3], right_ratio[4]);
      fflush(stderr);
      limited.alpha = 0.;
      limited.valid = 1;
      for (int material = 0; material < SIMPLEX_FLUX_MATERIALS;
           material++) {
        limited.fraction[material] = low[material];
        limited.correction_linf = max(limited.correction_linf,
          fabs(low[material] - high[material]));
      }
    }
    if (limited.alpha < 1. - 64.*DBL_EPSILON)
      limited_faces++;
    minimum_flux_alpha = min(minimum_flux_alpha, limited.alpha);
    maximum_simplex_closure_defect = max(
      maximum_simplex_closure_defect, limited.closure_defect);
    maximum_flux_fraction_correction = max(
      maximum_flux_fraction_correction, limited.correction_linf);

    maximum_flux_courant_correction = max(
      maximum_flux_courant_correction,
      fabs(courant)*limited.correction_linf);

    lower_gas_flux[] = limited.fraction[2]*uf.x[];
    lower_envelope_flux[] =
      (limited.fraction[1] + limited.fraction[2])*uf.x[];
    upper_gas_flux[] = limited.fraction[4]*uf.x[];
    upper_envelope_flux[] =
      (limited.fraction[3] + limited.fraction[4])*uf.x[];
    AggregatePhysicalPartition limited_face = {
      limited.fraction[0],
      limited.fraction[1] + limited.fraction[3],
      limited.fraction[2] + limited.fraction[4],
      0., 0., 0., 1
    };
    double direction = sign(courant);
    int upwind = -(direction + 1.)/2.;
#define AGGREGATE_LIMITED_Q_FLUX(material, momentum, gradient, first, second) do { \
      (first)[] = aggregate_phase_momentum_face_flux_x(                         \
        point, (material), (momentum).x, (gradient).x,                          \
        upwind, direction, courant,                                             \
        aggregate_physical_fraction(&limited_face, (material)));                \
      (second)[] = aggregate_phase_momentum_face_flux_x(                        \
        point, (material), (momentum).y, (gradient).y,                          \
        upwind, direction, courant,                                             \
        aggregate_physical_fraction(&limited_face, (material)));                \
    } while (0)
    AGGREGATE_LIMITED_Q_FLUX(AGGREGATE_MATERIAL_WATER,
      water_momentum, water_gradient, water_first_flux, water_second_flux);
    AGGREGATE_LIMITED_Q_FLUX(AGGREGATE_MATERIAL_OIL,
      oil_momentum, oil_gradient, oil_first_flux, oil_second_flux);
    AGGREGATE_LIMITED_Q_FLUX(AGGREGATE_MATERIAL_GAS,
      gas_momentum, gas_gradient, gas_first_flux, gas_second_flux);
#undef AGGREGATE_LIMITED_Q_FLUX
  }
  invalid_faces += invalid_limiter_faces;
  invalid_cells += invalid_low_order_cells;
  dual_momentum_compression_fallback_count +=
    compression_partition_fallbacks;
  dual_momentum_compression_fallback_count_max = max(
    dual_momentum_compression_fallback_count_max,
    dual_momentum_compression_fallback_count);
  dual_momentum_limited_face_count += limited_faces;
  dual_momentum_limited_face_count_max = max(
    dual_momentum_limited_face_count_max,
    dual_momentum_limited_face_count);
  dual_momentum_min_flux_alpha = min(
    dual_momentum_min_flux_alpha, minimum_flux_alpha);
  dual_momentum_min_flux_alpha_ever = min(
    dual_momentum_min_flux_alpha_ever, minimum_flux_alpha);
  dual_momentum_max_simplex_closure_defect = max(
    dual_momentum_max_simplex_closure_defect,
    maximum_simplex_closure_defect);
  dual_momentum_max_flux_fraction_correction = max(
    dual_momentum_max_flux_fraction_correction,
    maximum_flux_fraction_correction);
  dual_momentum_max_flux_courant_correction = max(
    dual_momentum_max_flux_courant_correction,
    maximum_flux_courant_correction);
#endif
  dual_momentum_face_projection_fallback_count +=
    face_projection_fallbacks;
  dual_momentum_face_projection_fallback_count_max = max(
    dual_momentum_face_projection_fallback_count_max,
    dual_momentum_face_projection_fallback_count);
  dual_momentum_max_rejected_face_projection_defect = max(
    dual_momentum_max_rejected_face_projection_defect,
    maximum_rejected_face_defect);
  dual_momentum_max_rejected_face_projection_correction = max(
    dual_momentum_max_rejected_face_projection_correction,
    maximum_rejected_face_correction);
  dual_momentum_max_rejected_face_transport_correction = max(
    dual_momentum_max_rejected_face_transport_correction,
    maximum_rejected_face_transport_correction);

  if (maximum_cfl > 0.5 + 1.e-6)
    fprintf(ferr,
      "aggregate_material_momentum_transport.h:%d: warning: "
      "VOF CFL exceeds 0.5 by %g\n",
      LINENO, maximum_cfl - 0.5), fflush(ferr);

  dual_momentum_invalid_cell_count = invalid_cells;
  dual_momentum_invalid_face_count = invalid_faces;
  dual_momentum_invalid_cell_count_max =
    max(dual_momentum_invalid_cell_count_max, invalid_cells);
  dual_momentum_invalid_face_count_max =
    max(dual_momentum_invalid_face_count_max, invalid_faces);
  dual_momentum_max_cell_partition_defect =
    max(dual_momentum_max_cell_partition_defect, maximum_cell_defect);
  dual_momentum_max_face_partition_defect =
    max(dual_momentum_max_face_partition_defect, maximum_face_defect);
  dual_momentum_max_partition_correction =
    max(dual_momentum_max_partition_correction, maximum_correction);
  dual_momentum_max_vof_bound_defect =
    max(dual_momentum_max_vof_bound_defect, maximum_bound_defect);
  if (invalid_cells > 0 || invalid_faces > 0)
    dual_momentum_transport_failed = 1;

  foreach() {
    double divergence = uf.x[1] - uf.x[];
#define AGGREGATE_Q_UPDATE(momentum, first, second, compression) do { \
      (momentum).x[] += dt*((first)[] - (first)[1] +                 \
        (compression).x[]*divergence)/(cm[]*Delta);                 \
      (momentum).y[] += dt*((second)[] - (second)[1] +              \
        (compression).y[]*divergence)/(cm[]*Delta);                 \
    } while (0)
    AGGREGATE_Q_UPDATE(water_momentum, water_first_flux,
      water_second_flux, water_compression);
    AGGREGATE_Q_UPDATE(oil_momentum, oil_first_flux,
      oil_second_flux, oil_compression);
    AGGREGATE_Q_UPDATE(gas_momentum, gas_first_flux,
      gas_second_flux, gas_compression);
#undef AGGREGATE_Q_UPDATE
    lower_gas_fraction[] += dt*(lower_gas_flux[] - lower_gas_flux[1] +
      lower_gas_step[]*divergence)/(cm[]*Delta);
    lower_envelope_fraction[] += dt*(lower_envelope_flux[] -
      lower_envelope_flux[1] + lower_envelope_step[]*divergence)/
      (cm[]*Delta);
    upper_gas_fraction[] += dt*(upper_gas_flux[] - upper_gas_flux[1] +
      upper_gas_step[]*divergence)/(cm[]*Delta);
    upper_envelope_fraction[] += dt*(upper_envelope_flux[] -
      upper_envelope_flux[1] + upper_envelope_step[]*divergence)/
      (cm[]*Delta);
  }
  boundary((scalar *){
    lower_gas_fraction, lower_envelope_fraction,
    upper_gas_fraction, upper_envelope_fraction,
    water_momentum, oil_momentum, gas_momentum
  });
}

static scalar * aggregate_saved_interfaces = NULL;

event vof (i++)
{
  clock_t transport_started = clock();
  dual_momentum_limited_face_count = 0;
  dual_momentum_min_flux_alpha = 1.;
  dual_momentum_max_flux_fraction_correction = 0.;
  dual_momentum_max_flux_courant_correction = 0.;
  dual_momentum_face_projection_fallback_count = 0;
  dual_momentum_nested_plic_fallback_count = 0;
  dual_momentum_compression_fallback_count = 0;
  aggregate_current_transport_time = t;
  aggregate_current_transport_iteration = i;
  for (scalar indicator in interfaces)
    if (indicator.tracers) {
      if (pid() == 0)
        fprintf(stderr,
          "aggregate material momentum transport does not accept "
          "external VOF tracers\n");
      dual_momentum_transport_failed = 1;
    }

  vector water_momentum[], oil_momentum[], gas_momentum[];
  vector water_compression[], oil_compression[], gas_compression[];
  scalar lower_gas_step[], lower_envelope_step[];
  scalar upper_gas_step[], upper_envelope_step[];
  foreach_dimension() {
    water_momentum.x.aggregate_material_id = AGGREGATE_MATERIAL_WATER;
    oil_momentum.x.aggregate_material_id = AGGREGATE_MATERIAL_OIL;
    gas_momentum.x.aggregate_material_id = AGGREGATE_MATERIAL_GAS;
    water_momentum.x.gradient = u.x.gradient;
    oil_momentum.x.gradient = u.x.gradient;
    gas_momentum.x.gradient = u.x.gradient;
  }
  for (scalar material_momentum in {
      water_momentum, oil_momentum, gas_momentum
    }) {
    material_momentum.depends = list_add(
      material_momentum.depends, lower_gas_fraction);
    material_momentum.depends = list_add(
      material_momentum.depends, lower_envelope_fraction);
    material_momentum.depends = list_add(
      material_momentum.depends, upper_gas_fraction);
    material_momentum.depends = list_add(
      material_momentum.depends, upper_envelope_fraction);
    foreach_dimension()
      material_momentum.v.x.i = -1;
  }
  for (int boundary_index = 0; boundary_index < nboundary;
       boundary_index++)
    foreach_dimension()
      for (scalar component in {
          water_momentum.x, oil_momentum.x, gas_momentum.x
        })
        component.boundary[boundary_index] =
          aggregate_boundary_momentum_x;
#if TREE
  foreach_dimension()
    for (scalar component in {
        water_momentum.x, oil_momentum.x, gas_momentum.x
      })
      component.prolongation = aggregate_prolong_momentum_x;
#endif

  int invalid_initial_partition = 0;
  double initial_water_defect = 0., initial_lower_oil_defect = 0.;
  double initial_upper_oil_defect = 0., initial_bound_defect = 0.;
  foreach(reduction(max:invalid_initial_partition)
          reduction(max:initial_water_defect)
          reduction(max:initial_lower_oil_defect)
          reduction(max:initial_upper_oil_defect)
          reduction(max:initial_bound_defect)) {
    PositiveMaterialPartition identity;
    AggregatePhysicalPartition partition = aggregate_transport_partition(
      lower_gas_fraction[], lower_envelope_fraction[],
      upper_gas_fraction[], upper_envelope_fraction[], &identity);
    initial_water_defect = max(
      initial_water_defect, identity.water_defect);
    initial_lower_oil_defect = max(
      initial_lower_oil_defect, identity.lower_oil_defect);
    initial_upper_oil_defect = max(
      initial_upper_oil_defect, identity.upper_oil_defect);
    initial_bound_defect = max(initial_bound_defect, identity.bound_defect);
    if (!partition.valid) {
      aggregate_report_invalid_cell(
        "initial", identity,
        lower_gas_fraction[], lower_envelope_fraction[],
        upper_gas_fraction[], upper_envelope_fraction[], x, y, Delta);
      invalid_initial_partition = 1;
      PositiveMaterialPartition water = aggregate_identity_partition(
        0., 0., 0., 0.);
      partition = aggregate_physical_partition(&water);
    }
    foreach_dimension() {
      water_momentum.x[] = partition.water*rho_water*u.x[];
      oil_momentum.x[] = partition.oil*rho_oil*u.x[];
      gas_momentum.x[] = partition.gas*rho_gas*u.x[];
    }

    int dominant = aggregate_physical_dominant(&partition);
#define AGGREGATE_COMPRESSION(material, fraction, momentum, compression) do { \
      foreach_dimension()                                                   \
        (compression).x[] = (dominant == (material) &&                     \
                             (fraction) > 1.e-10) ?                         \
          (momentum).x[]/(fraction) : 0.;                                   \
    } while (0)
    AGGREGATE_COMPRESSION(AGGREGATE_MATERIAL_WATER,
      partition.water, water_momentum, water_compression);
    AGGREGATE_COMPRESSION(AGGREGATE_MATERIAL_OIL,
      partition.oil, oil_momentum, oil_compression);
    AGGREGATE_COMPRESSION(AGGREGATE_MATERIAL_GAS,
      partition.gas, gas_momentum, gas_compression);
#undef AGGREGATE_COMPRESSION
    lower_gas_step[] = lower_gas_fraction[] > 0.5;
    lower_envelope_step[] = lower_envelope_fraction[] > 0.5;
    upper_gas_step[] = upper_gas_fraction[] > 0.5;
    upper_envelope_step[] = upper_envelope_fraction[] > 0.5;
  }
  aggregate_report_stage_reject(
    "initial", invalid_initial_partition, 0,
    initial_water_defect, initial_lower_oil_defect,
    initial_upper_oil_defect, initial_bound_defect);
  if (invalid_initial_partition)
    dual_momentum_transport_failed = 1;
  boundary((scalar *){water_momentum, oil_momentum, gas_momentum});

  void (* sweep[dimension]) (
    vector, vector, vector,
    scalar, scalar, scalar, scalar,
    vector, vector, vector);
  int direction = 0;
  foreach_dimension()
    sweep[direction++] = aggregate_momentum_sweep_x;
  for (direction = 0; direction < dimension; direction++)
  {
    aggregate_current_sweep_axis = (i + direction) % dimension;
    sweep[aggregate_current_sweep_axis](
      water_momentum, oil_momentum, gas_momentum,
      lower_gas_step, lower_envelope_step,
      upper_gas_step, upper_envelope_step,
      water_compression, oil_compression, gas_compression);
  }
  aggregate_current_sweep_axis = -1;

  long invalid_cells = 0;
  double maximum_cell_defect = 0., maximum_bound_defect = 0.;
  double maximum_correction = 0.;
  double maximum_water_defect = 0., maximum_lower_oil_defect = 0.;
  double maximum_upper_oil_defect = 0.;
  foreach(reduction(+:invalid_cells) reduction(max:maximum_cell_defect)
          reduction(max:maximum_bound_defect)
          reduction(max:maximum_correction)
          reduction(max:maximum_water_defect)
          reduction(max:maximum_lower_oil_defect)
          reduction(max:maximum_upper_oil_defect)) {
    PositiveMaterialPartition identity;
    AggregatePhysicalPartition partition = aggregate_transport_partition(
      lower_gas_fraction[], lower_envelope_fraction[],
      upper_gas_fraction[], upper_envelope_fraction[], &identity);
    maximum_cell_defect = max(
      maximum_cell_defect, identity.identity_defect);
    maximum_bound_defect = max(
      maximum_bound_defect, identity.bound_defect);
    maximum_correction = max(maximum_correction, identity.correction_l1);
    maximum_water_defect = max(
      maximum_water_defect, identity.water_defect);
    maximum_lower_oil_defect = max(
      maximum_lower_oil_defect, identity.lower_oil_defect);
    maximum_upper_oil_defect = max(
      maximum_upper_oil_defect, identity.upper_oil_defect);
    if (!partition.valid) {
      aggregate_report_invalid_cell(
        "post-sweeps", identity,
        lower_gas_fraction[], lower_envelope_fraction[],
        upper_gas_fraction[], upper_envelope_fraction[], x, y, Delta);
      invalid_cells++;
      PositiveMaterialPartition water = aggregate_identity_partition(
        0., 0., 0., 0.);
      partition = aggregate_physical_partition(&water);
    }
    double density = aggregate_physical_density(
      &partition, rho_water, rho_oil, rho_gas);
    foreach_dimension()
      u.x[] = (water_momentum.x[] + oil_momentum.x[] +
               gas_momentum.x[])/(density + SEPS);
  }
  aggregate_report_stage_reject(
    "post-sweeps", invalid_cells, 0,
    maximum_water_defect, maximum_lower_oil_defect,
    maximum_upper_oil_defect, maximum_bound_defect);
  boundary((scalar *){u});

  dual_momentum_invalid_cell_count = invalid_cells;
  dual_momentum_invalid_cell_count_max =
    max(dual_momentum_invalid_cell_count_max, invalid_cells);
  dual_momentum_max_cell_partition_defect =
    max(dual_momentum_max_cell_partition_defect, maximum_cell_defect);
  dual_momentum_max_partition_correction =
    max(dual_momentum_max_partition_correction, maximum_correction);
  dual_momentum_max_vof_bound_defect =
    max(dual_momentum_max_vof_bound_defect, maximum_bound_defect);
  if (invalid_cells > 0)
    dual_momentum_transport_failed = 1;

  aggregate_saved_interfaces = interfaces;
  interfaces = NULL;
  dual_momentum_last_transport_cpu_seconds =
    (double)(clock() - transport_started)/CLOCKS_PER_SEC;
}

event tracer_advection (i++)
{
  interfaces = aggregate_saved_interfaces;
}

event aggregate_material_momentum_guard (i++, last)
{
  if (dual_momentum_transport_failed) {
    if (pid() == 0)
      fprintf(stderr,
        "aggregate material momentum transport rejected state at "
        "t=%g i=%d: invalid_cells=%ld invalid_faces=%ld "
        "invalid_cells_max=%ld invalid_faces_max=%ld "
        "max_cell_defect=%g max_face_defect=%g bound_defect=%g "
        "correction_l1=%g\n",
        t, i, dual_momentum_invalid_cell_count,
        dual_momentum_invalid_face_count,
        dual_momentum_invalid_cell_count_max,
        dual_momentum_invalid_face_count_max,
        dual_momentum_max_cell_partition_defect,
        dual_momentum_max_face_partition_defect,
        dual_momentum_max_vof_bound_defect,
        dual_momentum_max_partition_correction);
    return 1;
  }
}

#endif
