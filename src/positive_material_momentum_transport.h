#ifndef POSITIVE_MATERIAL_MOMENTUM_TRANSPORT_H
#define POSITIVE_MATERIAL_MOMENTUM_TRANSPORT_H

#include <time.h>
#include "positive_material_partition.h"

#if dimension != 2
# error "positive material momentum transport is currently limited to 2D axisymmetry"
#endif

#ifndef POSITIVE_MATERIAL_PARTITION_TOLERANCE
# define POSITIVE_MATERIAL_PARTITION_TOLERANCE 1.e-8
#endif
#ifndef POSITIVE_MATERIAL_VOF_BOUND_TOLERANCE
# define POSITIVE_MATERIAL_VOF_BOUND_TOLERANCE 1.e-3
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
long dual_momentum_nested_plic_fallback_count = 0;
long dual_momentum_nested_plic_fallback_count_max = 0;

attribute {
  int positive_material_id;
}

static inline PositiveMaterialPartition positive_material_transport_partition (
  double lower_gas, double lower_envelope,
  double upper_gas, double upper_envelope)
{
  return positive_material_partition_with_tolerances(
    lower_gas, lower_envelope, upper_gas, upper_envelope,
    (double)POSITIVE_MATERIAL_PARTITION_TOLERANCE,
    (double)POSITIVE_MATERIAL_VOF_BOUND_TOLERANCE);
}

static inline double positive_material_transport_density (
  Point point, PositiveMaterialPartition * output)
{
  PositiveMaterialPartition partition = positive_material_transport_partition(
    lower_gas_fraction[], lower_envelope_fraction[],
    upper_gas_fraction[], upper_envelope_fraction[]);
  if (output)
    *output = partition;
  return positive_material_density(&partition, rho_water, rho_oil, rho_gas);
}

static inline double positive_material_transport_rho (int material)
{
  return material == POSITIVE_MATERIAL_WATER ? rho_water :
    (material == POSITIVE_MATERIAL_LOWER_OIL ||
     material == POSITIVE_MATERIAL_UPPER_OIL) ? rho_oil : rho_gas;
}

#if TREE
static void positive_material_momentum_refine (Point point, scalar velocity)
{
  double parent_velocity = velocity[];
  double parent_metric = cm[];
  double parent_density = positive_material_transport_density(point, NULL);
  refine_bilinear(point, velocity);

  double child_momentum = 0., child_mass = 0.;
  foreach_child() {
    double density = positive_material_transport_density(point, NULL);
    child_momentum += cm[]*density*velocity[];
    child_mass += cm[]*density;
  }
  double target_momentum =
    (1 << dimension)*parent_metric*parent_density*parent_velocity;
  double correction = (target_momentum - child_momentum)/
                      (child_mass + SEPS);
  foreach_child()
    velocity[] += correction;
}

static void positive_material_momentum_restriction (
  Point point, scalar velocity)
{
  double child_momentum = 0.;
  foreach_child()
    child_momentum += cm[]*
      positive_material_transport_density(point, NULL)*velocity[];
  double parent_mass = (1 << dimension)*(cm[] + SEPS)*
    positive_material_transport_density(point, NULL);
  velocity[] = child_momentum/(parent_mass + SEPS);
}

static void positive_material_move_to_front (scalar field)
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
  positive_material_move_to_front(upper_envelope_fraction);
  positive_material_move_to_front(upper_gas_fraction);
  positive_material_move_to_front(lower_envelope_fraction);
  positive_material_move_to_front(lower_gas_fraction);
  foreach_dimension() {
    u.x.refine = u.x.prolongation = positive_material_momentum_refine;
    u.x.restriction = positive_material_momentum_restriction;
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
static double positive_material_boundary_q_x (
  Point neighbor, Point point, scalar momentum, bool * data)
{
  PositiveMaterialPartition partition;
  positive_material_transport_density(point, &partition);
  double phase_fraction = positive_material_fraction(
    &partition, momentum.positive_material_id);
  return phase_fraction*
    positive_material_transport_rho(momentum.positive_material_id)*u.x[];
}

#if TREE
foreach_dimension()
static void positive_material_prolong_q_x (Point point, scalar momentum)
{
  foreach_child() {
    PositiveMaterialPartition partition;
    positive_material_transport_density(point, &partition);
    double phase_fraction = positive_material_fraction(
      &partition, momentum.positive_material_id);
    momentum[] = phase_fraction*
      positive_material_transport_rho(momentum.positive_material_id)*u.x[];
  }
}
#endif

foreach_dimension()
static double positive_material_geometric_face_fraction_x (
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
static double positive_material_concentration_gradient_x (
  Point point, scalar momentum, int material)
{
  PositiveMaterialPartition left = positive_material_transport_partition(
    lower_gas_fraction[-1], lower_envelope_fraction[-1],
    upper_gas_fraction[-1], upper_envelope_fraction[-1]);
  PositiveMaterialPartition center = positive_material_transport_partition(
    lower_gas_fraction[], lower_envelope_fraction[],
    upper_gas_fraction[], upper_envelope_fraction[]);
  PositiveMaterialPartition right = positive_material_transport_partition(
    lower_gas_fraction[1], lower_envelope_fraction[1],
    upper_gas_fraction[1], upper_envelope_fraction[1]);
  if (!left.valid || !center.valid || !right.valid)
    return 0.;
  double cl = positive_material_fraction(&left, material);
  double cc = positive_material_fraction(&center, material);
  double cr = positive_material_fraction(&right, material);
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
static double positive_material_tracer_flux_x (
  Point point, int material, scalar momentum, scalar gradient,
  int upwind, double direction, double courant,
  double face_fraction)
{
  PositiveMaterialPartition upwind_partition =
    positive_material_transport_partition(
    lower_gas_fraction[upwind], lower_envelope_fraction[upwind],
    upper_gas_fraction[upwind], upper_envelope_fraction[upwind]);
  if (!upwind_partition.valid)
    return 0.;
  double cell_fraction = positive_material_fraction(
    &upwind_partition, material);
  if (cell_fraction <= 1.e-10)
    return 0.;
  double concentration = momentum[upwind]/cell_fraction +
    direction*min(1., 1. - direction*courant)*
    gradient[upwind]*Delta/2.;
  return concentration*face_fraction*uf.x[];
}

foreach_dimension()
static void positive_material_momentum_sweep_x (
  vector water_momentum, vector lower_oil_momentum,
  vector lower_gas_momentum, vector upper_oil_momentum,
  vector upper_gas_momentum,
  scalar lower_gas_step, scalar lower_envelope_step,
  scalar upper_gas_step, scalar upper_envelope_step,
  vector water_compression, vector lower_oil_compression,
  vector lower_gas_compression, vector upper_oil_compression,
  vector upper_gas_compression)
{
  vector lower_gas_normal[], lower_envelope_normal[];
  vector upper_gas_normal[], upper_envelope_normal[];
  scalar lower_gas_intercept[], lower_envelope_intercept[];
  scalar upper_gas_intercept[], upper_envelope_intercept[];
  scalar lower_gas_flux[], lower_envelope_flux[];
  scalar upper_gas_flux[], upper_envelope_flux[];

  vector water_gradient[], lower_oil_gradient[], lower_gas_gradient[];
  vector upper_oil_gradient[], upper_gas_gradient[];
  scalar water_first_flux[], water_second_flux[];
  scalar lower_oil_first_flux[], lower_oil_second_flux[];
  scalar lower_gas_first_flux[], lower_gas_second_flux[];
  scalar upper_oil_first_flux[], upper_oil_second_flux[];
  scalar upper_gas_first_flux[], upper_gas_second_flux[];

  reconstruction(lower_gas_fraction,
                 lower_gas_normal, lower_gas_intercept);
  reconstruction(lower_envelope_fraction,
                 lower_envelope_normal, lower_envelope_intercept);
  reconstruction(upper_gas_fraction,
                 upper_gas_normal, upper_gas_intercept);
  reconstruction(upper_envelope_fraction,
                 upper_envelope_normal, upper_envelope_intercept);

  long invalid_cells = 0;
  double maximum_cell_defect = 0., maximum_bound_defect = 0.;
  double maximum_correction = 0.;
  foreach(reduction(+:invalid_cells) reduction(max:maximum_cell_defect)
          reduction(max:maximum_bound_defect)
          reduction(max:maximum_correction)) {
    PositiveMaterialPartition partition = positive_material_transport_partition(
      lower_gas_fraction[], lower_envelope_fraction[],
      upper_gas_fraction[], upper_envelope_fraction[]);
    maximum_cell_defect = max(
      maximum_cell_defect, partition.identity_defect);
    maximum_bound_defect = max(
      maximum_bound_defect, partition.bound_defect);
    maximum_correction = max(maximum_correction, partition.correction_l1);
    if (!partition.valid) {
      invalid_cells++;
      partition = positive_material_transport_partition(0., 0., 0., 0.);
    }
  }

  foreach() {
    water_gradient.x[] = positive_material_concentration_gradient_x(
      point, water_momentum.x, POSITIVE_MATERIAL_WATER);
    water_gradient.y[] = positive_material_concentration_gradient_x(
      point, water_momentum.y, POSITIVE_MATERIAL_WATER);
    lower_oil_gradient.x[] = positive_material_concentration_gradient_x(
      point, lower_oil_momentum.x, POSITIVE_MATERIAL_LOWER_OIL);
    lower_oil_gradient.y[] = positive_material_concentration_gradient_x(
      point, lower_oil_momentum.y, POSITIVE_MATERIAL_LOWER_OIL);
    lower_gas_gradient.x[] = positive_material_concentration_gradient_x(
      point, lower_gas_momentum.x, POSITIVE_MATERIAL_LOWER_GAS);
    lower_gas_gradient.y[] = positive_material_concentration_gradient_x(
      point, lower_gas_momentum.y, POSITIVE_MATERIAL_LOWER_GAS);
    upper_oil_gradient.x[] = positive_material_concentration_gradient_x(
      point, upper_oil_momentum.x, POSITIVE_MATERIAL_UPPER_OIL);
    upper_oil_gradient.y[] = positive_material_concentration_gradient_x(
      point, upper_oil_momentum.y, POSITIVE_MATERIAL_UPPER_OIL);
    upper_gas_gradient.x[] = positive_material_concentration_gradient_x(
      point, upper_gas_momentum.x, POSITIVE_MATERIAL_UPPER_GAS);
    upper_gas_gradient.y[] = positive_material_concentration_gradient_x(
      point, upper_gas_momentum.y, POSITIVE_MATERIAL_UPPER_GAS);
  }

  double maximum_face_defect = 0.;
  long invalid_faces = 0;
  double maximum_cfl = 0.;
  foreach_face(x, reduction(max:maximum_face_defect)
                  reduction(+:invalid_faces) reduction(max:maximum_cfl)
                  reduction(max:maximum_bound_defect)
                  reduction(max:maximum_correction)) {
    double courant = uf.x[]*dt/(Delta*fm.x[] + SEPS);
    double direction = sign(courant);
    int upwind = -(direction + 1.)/2.;
    maximum_cfl = max(maximum_cfl,
                      courant*fm.x[]*direction/(cm[] + SEPS));

    double lower_gas_face = positive_material_geometric_face_fraction_x(
      point, lower_gas_fraction, lower_gas_normal, lower_gas_intercept,
      upwind, direction, courant);
    double lower_envelope_face = positive_material_geometric_face_fraction_x(
      point, lower_envelope_fraction, lower_envelope_normal,
      lower_envelope_intercept, upwind, direction, courant);
    double upper_gas_face = positive_material_geometric_face_fraction_x(
      point, upper_gas_fraction, upper_gas_normal, upper_gas_intercept,
      upwind, direction, courant);
    double upper_envelope_face = positive_material_geometric_face_fraction_x(
      point, upper_envelope_fraction, upper_envelope_normal,
      upper_envelope_intercept, upwind, direction, courant);

    lower_gas_flux[] = lower_gas_face*uf.x[];
    lower_envelope_flux[] = lower_envelope_face*uf.x[];
    upper_gas_flux[] = upper_gas_face*uf.x[];
    upper_envelope_flux[] = upper_envelope_face*uf.x[];

    PositiveMaterialPartition face = positive_material_transport_partition(
      lower_gas_face, lower_envelope_face,
      upper_gas_face, upper_envelope_face);
    maximum_face_defect = max(
      maximum_face_defect, face.identity_defect);
    maximum_bound_defect = max(
      maximum_bound_defect, face.bound_defect);
    maximum_correction = max(maximum_correction, face.correction_l1);
    if (!face.valid) {
      invalid_faces++;
      face = positive_material_transport_partition(0., 0., 0., 0.);
    }

#define POSITIVE_Q_FLUX(material, momentum, gradient, first_flux, second_flux) do { \
      (first_flux)[] = positive_material_tracer_flux_x(                     \
        point, (material), (momentum).x, (gradient).x, upwind, direction,    \
        courant, positive_material_fraction(&face, (material)));            \
      (second_flux)[] = positive_material_tracer_flux_x(                    \
        point, (material), (momentum).y, (gradient).y, upwind, direction,    \
        courant, positive_material_fraction(&face, (material)));            \
    } while (0)
    POSITIVE_Q_FLUX(POSITIVE_MATERIAL_WATER,
      water_momentum, water_gradient, water_first_flux, water_second_flux);
    POSITIVE_Q_FLUX(POSITIVE_MATERIAL_LOWER_OIL,
      lower_oil_momentum, lower_oil_gradient,
      lower_oil_first_flux, lower_oil_second_flux);
    POSITIVE_Q_FLUX(POSITIVE_MATERIAL_LOWER_GAS,
      lower_gas_momentum, lower_gas_gradient,
      lower_gas_first_flux, lower_gas_second_flux);
    POSITIVE_Q_FLUX(POSITIVE_MATERIAL_UPPER_OIL,
      upper_oil_momentum, upper_oil_gradient,
      upper_oil_first_flux, upper_oil_second_flux);
    POSITIVE_Q_FLUX(POSITIVE_MATERIAL_UPPER_GAS,
      upper_gas_momentum, upper_gas_gradient,
      upper_gas_first_flux, upper_gas_second_flux);
#undef POSITIVE_Q_FLUX
  }

  if (maximum_cfl > 0.5 + 1.e-6)
    fprintf(ferr,
      "positive_material_momentum_transport.h:%d: warning: "
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
#define POSITIVE_Q_UPDATE(momentum, first_flux, second_flux, compression) do { \
      (momentum).x[] += dt*((first_flux)[] - (first_flux)[1] +             \
        (compression).x[]*divergence)/(cm[]*Delta);                        \
      (momentum).y[] += dt*((second_flux)[] - (second_flux)[1] +           \
        (compression).y[]*divergence)/(cm[]*Delta);                        \
    } while (0)
    POSITIVE_Q_UPDATE(water_momentum, water_first_flux, water_second_flux,
                      water_compression);
    POSITIVE_Q_UPDATE(lower_oil_momentum,
                      lower_oil_first_flux, lower_oil_second_flux,
                      lower_oil_compression);
    POSITIVE_Q_UPDATE(lower_gas_momentum,
                      lower_gas_first_flux, lower_gas_second_flux,
                      lower_gas_compression);
    POSITIVE_Q_UPDATE(upper_oil_momentum,
                      upper_oil_first_flux, upper_oil_second_flux,
                      upper_oil_compression);
    POSITIVE_Q_UPDATE(upper_gas_momentum,
                      upper_gas_first_flux, upper_gas_second_flux,
                      upper_gas_compression);
#undef POSITIVE_Q_UPDATE

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
    water_momentum, lower_oil_momentum, lower_gas_momentum,
    upper_oil_momentum, upper_gas_momentum
  });
}

static scalar * positive_material_saved_interfaces = NULL;

event vof (i++)
{
  clock_t transport_started = clock();
  for (scalar indicator in interfaces)
    if (indicator.tracers) {
      if (pid() == 0)
        fprintf(stderr,
          "positive material momentum transport does not accept "
          "external VOF tracers\n");
      dual_momentum_transport_failed = 1;
    }

  vector water_momentum[], lower_oil_momentum[], lower_gas_momentum[];
  vector upper_oil_momentum[], upper_gas_momentum[];
  vector water_compression[], lower_oil_compression[];
  vector lower_gas_compression[], upper_oil_compression[];
  vector upper_gas_compression[];
  scalar lower_gas_step[], lower_envelope_step[];
  scalar upper_gas_step[], upper_envelope_step[];
  foreach_dimension() {
    water_momentum.x.positive_material_id = POSITIVE_MATERIAL_WATER;
    lower_oil_momentum.x.positive_material_id =
      POSITIVE_MATERIAL_LOWER_OIL;
    lower_gas_momentum.x.positive_material_id =
      POSITIVE_MATERIAL_LOWER_GAS;
    upper_oil_momentum.x.positive_material_id =
      POSITIVE_MATERIAL_UPPER_OIL;
    upper_gas_momentum.x.positive_material_id =
      POSITIVE_MATERIAL_UPPER_GAS;
    water_momentum.x.gradient = u.x.gradient;
    lower_oil_momentum.x.gradient = u.x.gradient;
    lower_gas_momentum.x.gradient = u.x.gradient;
    upper_oil_momentum.x.gradient = u.x.gradient;
    upper_gas_momentum.x.gradient = u.x.gradient;
  }
  for (scalar material_momentum in {
      water_momentum, lower_oil_momentum, lower_gas_momentum,
      upper_oil_momentum, upper_gas_momentum
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
          water_momentum.x, lower_oil_momentum.x, lower_gas_momentum.x,
          upper_oil_momentum.x, upper_gas_momentum.x
        })
        component.boundary[boundary_index] = positive_material_boundary_q_x;
#if TREE
  foreach_dimension()
    for (scalar component in {
        water_momentum.x, lower_oil_momentum.x, lower_gas_momentum.x,
        upper_oil_momentum.x, upper_gas_momentum.x
      })
      component.prolongation = positive_material_prolong_q_x;
#endif

  int invalid_initial_partition = 0;
  foreach(reduction(max:invalid_initial_partition)) {
    PositiveMaterialPartition partition = positive_material_transport_partition(
      lower_gas_fraction[], lower_envelope_fraction[],
      upper_gas_fraction[], upper_envelope_fraction[]);
    if (!partition.valid) {
      invalid_initial_partition = 1;
      partition = positive_material_transport_partition(0., 0., 0., 0.);
    }
    foreach_dimension() {
      water_momentum.x[] = partition.water*rho_water*u.x[];
      lower_oil_momentum.x[] = partition.lower_oil*rho_oil*u.x[];
      lower_gas_momentum.x[] = partition.lower_gas*rho_gas*u.x[];
      upper_oil_momentum.x[] = partition.upper_oil*rho_oil*u.x[];
      upper_gas_momentum.x[] = partition.upper_gas*rho_gas*u.x[];
    }
    int dominant = positive_material_dominant(&partition);
    lower_gas_step[] = lower_gas_fraction[] > 0.5;
    lower_envelope_step[] = lower_envelope_fraction[] > 0.5;
    upper_gas_step[] = upper_gas_fraction[] > 0.5;
    upper_envelope_step[] = upper_envelope_fraction[] > 0.5;
#define POSITIVE_COMPRESSION(material, fraction, momentum, compression) do { \
      foreach_dimension()                                                   \
        (compression).x[] = (dominant == (material) &&                     \
                             (fraction) > 1.e-10) ?                         \
          (momentum).x[]/(fraction) : 0.;                                   \
    } while (0)
    POSITIVE_COMPRESSION(POSITIVE_MATERIAL_WATER,
      partition.water, water_momentum, water_compression);
    POSITIVE_COMPRESSION(POSITIVE_MATERIAL_LOWER_OIL,
      partition.lower_oil, lower_oil_momentum, lower_oil_compression);
    POSITIVE_COMPRESSION(POSITIVE_MATERIAL_LOWER_GAS,
      partition.lower_gas, lower_gas_momentum, lower_gas_compression);
    POSITIVE_COMPRESSION(POSITIVE_MATERIAL_UPPER_OIL,
      partition.upper_oil, upper_oil_momentum, upper_oil_compression);
    POSITIVE_COMPRESSION(POSITIVE_MATERIAL_UPPER_GAS,
      partition.upper_gas, upper_gas_momentum, upper_gas_compression);
#undef POSITIVE_COMPRESSION
  }
  if (invalid_initial_partition)
    dual_momentum_transport_failed = 1;
  boundary((scalar *){
    water_momentum, lower_oil_momentum, lower_gas_momentum,
    upper_oil_momentum, upper_gas_momentum
  });

  void (* sweep[dimension]) (
    vector, vector, vector, vector, vector,
    scalar, scalar, scalar, scalar,
    vector, vector, vector, vector, vector);
  int direction = 0;
  foreach_dimension()
    sweep[direction++] = positive_material_momentum_sweep_x;
  for (direction = 0; direction < dimension; direction++) {
    sweep[(i + direction) % dimension](
      water_momentum, lower_oil_momentum, lower_gas_momentum,
      upper_oil_momentum, upper_gas_momentum,
      lower_gas_step, lower_envelope_step,
      upper_gas_step, upper_envelope_step,
      water_compression, lower_oil_compression, lower_gas_compression,
      upper_oil_compression, upper_gas_compression);
  }

  long invalid_cells = 0;
  double maximum_cell_defect = 0., maximum_bound_defect = 0.;
  double maximum_correction = 0.;
  foreach(reduction(+:invalid_cells) reduction(max:maximum_cell_defect)
          reduction(max:maximum_bound_defect)
          reduction(max:maximum_correction)) {
    PositiveMaterialPartition partition = positive_material_transport_partition(
      lower_gas_fraction[], lower_envelope_fraction[],
      upper_gas_fraction[], upper_envelope_fraction[]);
    maximum_cell_defect = max(
      maximum_cell_defect, partition.identity_defect);
    maximum_bound_defect = max(
      maximum_bound_defect, partition.bound_defect);
    maximum_correction = max(maximum_correction, partition.correction_l1);
    if (!partition.valid) {
      invalid_cells++;
      partition = positive_material_transport_partition(0., 0., 0., 0.);
    }
    double density = positive_material_density(
      &partition, rho_water, rho_oil, rho_gas);
    foreach_dimension()
      u.x[] = (water_momentum.x[] + lower_oil_momentum.x[] +
               lower_gas_momentum.x[] + upper_oil_momentum.x[] +
               upper_gas_momentum.x[])/(density + SEPS);
  }
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

  positive_material_saved_interfaces = interfaces;
  interfaces = NULL;
  dual_momentum_last_transport_cpu_seconds =
    (double)(clock() - transport_started)/CLOCKS_PER_SEC;
}

event tracer_advection (i++)
{
  interfaces = positive_material_saved_interfaces;
}

event positive_material_momentum_guard (i++, last)
{
  if (dual_momentum_transport_failed) {
    if (pid() == 0)
      fprintf(stderr,
        "positive material momentum transport rejected state at "
        "t=%g i=%d: invalid_cells=%ld invalid_faces=%ld "
        "max_cell_defect=%g max_face_defect=%g bound_defect=%g "
        "correction_l1=%g\n",
        t, i, dual_momentum_invalid_cell_count,
        dual_momentum_invalid_face_count,
        dual_momentum_max_cell_partition_defect,
        dual_momentum_max_face_partition_defect,
        dual_momentum_max_vof_bound_defect,
        dual_momentum_max_partition_correction);
    return 1;
  }
}

#endif
