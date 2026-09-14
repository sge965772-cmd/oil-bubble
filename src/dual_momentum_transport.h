#ifndef DUAL_MOMENTUM_TRANSPORT_H
#define DUAL_MOMENTUM_TRANSPORT_H

/*
 * Rejected experimental adapter retained only for differential regression.
 * ADR 0004 records why it must not be included by a production solver.
 */
#ifndef ENABLE_EXPERIMENTAL_DUAL_MOMENTUM_V1
# error "dual_momentum_transport.h is a rejected experimental v1 adapter"
#endif

#include "dual_momentum_algebra.h"

#if dimension != 2
# error "dual_momentum_transport.h is currently verified only in 2D axisymmetry"
#endif

#ifndef DUAL_MOMENTUM_INDICATOR_TOLERANCE
# define DUAL_MOMENTUM_INDICATOR_TOLERANCE 1.e-8
#endif
#ifndef DUAL_MOMENTUM_FACE_PARTITION_TOLERANCE
# define DUAL_MOMENTUM_FACE_PARTITION_TOLERANCE 1.e-8
#endif

long dual_momentum_invalid_cell_count = 0;
long dual_momentum_invalid_cell_count_max = 0;
long dual_momentum_invalid_face_count = 0;
long dual_momentum_invalid_face_count_max = 0;
double dual_momentum_max_face_partition_defect = 0.;
double dual_momentum_max_cell_partition_defect = 0.;
int dual_momentum_transport_failed = 0;

static inline double dual_momentum_transport_density (Point point)
{
  DualMomentumComponents components = dual_momentum_decompose(
    0., lower_gas_fraction[], lower_envelope_fraction[],
    upper_gas_fraction[], upper_envelope_fraction[],
    rho_water, rho_oil, rho_gas);
  return components.density;
}

#if TREE
static void dual_momentum_refine (Point point, scalar velocity)
{
  double parent_velocity = velocity[];
  double parent_metric = cm[];
  double parent_density = dual_momentum_transport_density(point);
  refine_bilinear(point, velocity);

  double child_momentum = 0., child_mass = 0.;
  foreach_child() {
    double density = dual_momentum_transport_density(point);
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

static void dual_momentum_restriction (Point point, scalar velocity)
{
  double child_momentum = 0.;
  foreach_child()
    child_momentum += cm[]*dual_momentum_transport_density(point)*velocity[];
  double parent_mass =
    (1 << dimension)*(cm[] + SEPS)*dual_momentum_transport_density(point);
  velocity[] = child_momentum/(parent_mass + SEPS);
}

static void dual_momentum_move_to_front (scalar field)
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
  dual_momentum_move_to_front(upper_envelope_fraction);
  dual_momentum_move_to_front(upper_gas_fraction);
  dual_momentum_move_to_front(lower_envelope_fraction);
  dual_momentum_move_to_front(lower_gas_fraction);
  foreach_dimension() {
    u.x.refine = u.x.prolongation = dual_momentum_refine;
    u.x.restriction = dual_momentum_restriction;
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
static double dual_momentum_boundary_x (
  Point neighbor, Point point, scalar momentum, bool * data)
{
  return dual_momentum_transport_density(point)*u.x[];
}

#if TREE
foreach_dimension()
static void dual_momentum_prolong_x (Point point, scalar momentum)
{
  foreach_child()
    momentum[] = dual_momentum_transport_density(point)*u.x[];
}
#endif

foreach_dimension()
static double dual_momentum_geometric_face_fraction_x (
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
static void dual_momentum_sweep_x (vector momentum)
{
  vector lower_gas_normal[], lower_envelope_normal[];
  vector upper_gas_normal[], upper_envelope_normal[];
  scalar lower_gas_intercept[], lower_envelope_intercept[];
  scalar upper_gas_intercept[], upper_envelope_intercept[];
  scalar lower_gas_flux[], lower_envelope_flux[];
  scalar upper_gas_flux[], upper_envelope_flux[];
  scalar lower_gas_step[], lower_envelope_step[];
  scalar upper_gas_step[], upper_envelope_step[];
  scalar axial_momentum_flux[], radial_momentum_flux[];
  tensor velocity_gradient[];

  reconstruction(lower_gas_fraction,
                 lower_gas_normal, lower_gas_intercept);
  reconstruction(lower_envelope_fraction,
                 lower_envelope_normal, lower_envelope_intercept);
  reconstruction(upper_gas_fraction,
                 upper_gas_normal, upper_gas_intercept);
  reconstruction(upper_envelope_fraction,
                 upper_envelope_normal, upper_envelope_intercept);
  gradients((scalar *){u}, (vector *){velocity_gradient});

  foreach() {
    lower_gas_step[] = lower_gas_fraction[] > 0.5;
    lower_envelope_step[] = lower_envelope_fraction[] > 0.5;
    upper_gas_step[] = upper_gas_fraction[] > 0.5;
    upper_envelope_step[] = upper_envelope_fraction[] > 0.5;
  }

  double maximum_face_defect = 0.;
  long invalid_faces = 0;
  double maximum_cfl = 0.;
  foreach_face(x, reduction(max:maximum_face_defect)
                  reduction(+:invalid_faces) reduction(max:maximum_cfl)) {
    double courant = uf.x[]*dt/(Delta*fm.x[] + SEPS);
    double direction = sign(courant);
    int upwind = -(direction + 1.)/2.;
    maximum_cfl = max(maximum_cfl,
                      courant*fm.x[]*direction/(cm[] + SEPS));

    double lower_gas_face = dual_momentum_geometric_face_fraction_x(
      point, lower_gas_fraction, lower_gas_normal, lower_gas_intercept,
      upwind, direction, courant);
    double lower_envelope_face = dual_momentum_geometric_face_fraction_x(
      point, lower_envelope_fraction, lower_envelope_normal,
      lower_envelope_intercept, upwind, direction, courant);
    double upper_gas_face = dual_momentum_geometric_face_fraction_x(
      point, upper_gas_fraction, upper_gas_normal, upper_gas_intercept,
      upwind, direction, courant);
    double upper_envelope_face = dual_momentum_geometric_face_fraction_x(
      point, upper_envelope_fraction, upper_envelope_normal,
      upper_envelope_intercept, upwind, direction, courant);

    lower_gas_flux[] = lower_gas_face*uf.x[];
    lower_envelope_flux[] = lower_envelope_face*uf.x[];
    upper_gas_flux[] = upper_gas_face*uf.x[];
    upper_envelope_flux[] = upper_envelope_face*uf.x[];

    double face_defect = max(
      max(lower_gas_face - lower_envelope_face,
          upper_gas_face - upper_envelope_face),
      lower_envelope_face + upper_envelope_face - 1.);
    face_defect = max(face_defect, 0.);
    maximum_face_defect = max(maximum_face_defect, face_defect);
    if (face_defect > (double)DUAL_MOMENTUM_FACE_PARTITION_TOLERANCE)
      invalid_faces++;

    double face_density = rho_water +
      (rho_oil - rho_water)*
        (lower_envelope_face + upper_envelope_face) +
      (rho_gas - rho_oil)*(lower_gas_face + upper_gas_face);
    if (!(face_density > 0.) || !isfinite(face_density)) {
      invalid_faces++;
      face_density = max(face_density, 1.e-12);
    }
    double mass_flux = face_density*uf.x[];

    double advected_axial_velocity =
      u.x[upwind] + direction*min(1., 1. - direction*courant)*
      velocity_gradient.x.x[upwind]*Delta/2.;
    double advected_radial_velocity =
      u.y[upwind] + direction*min(1., 1. - direction*courant)*
      velocity_gradient.y.x[upwind]*Delta/2.;
#if dimension > 1
    if (fm.y[upwind] && fm.y[upwind,1]) {
      double transverse_velocity =
        (uf.y[upwind] + uf.y[upwind,1])/
        (fm.y[upwind] + fm.y[upwind,1]);
      double axial_difference = transverse_velocity < 0. ?
        u.x[upwind,1] - u.x[upwind] :
        u.x[upwind] - u.x[upwind,-1];
      double radial_difference = transverse_velocity < 0. ?
        u.y[upwind,1] - u.y[upwind] :
        u.y[upwind] - u.y[upwind,-1];
      advected_axial_velocity -=
        dt*transverse_velocity*axial_difference/(2.*Delta);
      advected_radial_velocity -=
        dt*transverse_velocity*radial_difference/(2.*Delta);
    }
#endif
    axial_momentum_flux[] = mass_flux*advected_axial_velocity;
    radial_momentum_flux[] = mass_flux*advected_radial_velocity;
  }

  if (maximum_cfl > 0.5 + 1.e-6)
    fprintf(ferr,
      "dual_momentum_transport.h:%d: warning: VOF CFL exceeds 0.5 by %g\n",
      LINENO, maximum_cfl - 0.5), fflush(ferr);
  dual_momentum_invalid_face_count = invalid_faces;
  dual_momentum_invalid_face_count_max =
    max(dual_momentum_invalid_face_count_max, invalid_faces);
  dual_momentum_max_face_partition_defect =
    max(dual_momentum_max_face_partition_defect, maximum_face_defect);
  if (invalid_faces > 0)
    dual_momentum_transport_failed = 1;

  foreach() {
    double face_divergence = uf.x[1] - uf.x[];
    double step_density = rho_water +
      (rho_oil - rho_water)*
        (lower_envelope_step[] + upper_envelope_step[]) +
      (rho_gas - rho_oil)*(lower_gas_step[] + upper_gas_step[]);
    momentum.x[] += dt*(axial_momentum_flux[] -
      axial_momentum_flux[1] + step_density*u.x[]*face_divergence)/
      (cm[]*Delta);
    momentum.y[] += dt*(radial_momentum_flux[] -
      radial_momentum_flux[1] + step_density*u.y[]*face_divergence)/
      (cm[]*Delta);

    lower_gas_fraction[] += dt*(lower_gas_flux[] - lower_gas_flux[1] +
      lower_gas_step[]*face_divergence)/(cm[]*Delta);
    lower_envelope_fraction[] += dt*(lower_envelope_flux[] -
      lower_envelope_flux[1] + lower_envelope_step[]*face_divergence)/
      (cm[]*Delta);
    upper_gas_fraction[] += dt*(upper_gas_flux[] - upper_gas_flux[1] +
      upper_gas_step[]*face_divergence)/(cm[]*Delta);
    upper_envelope_fraction[] += dt*(upper_envelope_flux[] -
      upper_envelope_flux[1] + upper_envelope_step[]*face_divergence)/
      (cm[]*Delta);
  }
  boundary((scalar *){
    lower_gas_fraction, lower_envelope_fraction,
    upper_gas_fraction, upper_envelope_fraction, momentum
  });

  foreach() {
    double density = dual_momentum_transport_density(point);
    if (!(density > 0.) || !isfinite(density))
      density = 1.e-12;
    u.x[] = momentum.x[]/density;
    u.y[] = momentum.y[]/density;
  }
  boundary((scalar *){u});
}

static scalar * dual_momentum_saved_interfaces = NULL;

event vof (i++)
{
  for (scalar indicator in interfaces)
    if (indicator.tracers) {
      if (pid() == 0)
        fprintf(stderr,
          "dual momentum transport does not accept external VOF tracers\n");
      dual_momentum_transport_failed = 1;
    }

  vector momentum[];
  for (scalar component in {momentum})
    foreach_dimension()
      component.v.x.i = -1;
  for (int boundary_index = 0; boundary_index < nboundary;
       boundary_index++)
    foreach_dimension()
      momentum.x.boundary[boundary_index] = dual_momentum_boundary_x;
#if TREE
  foreach_dimension()
    momentum.x.prolongation = dual_momentum_prolong_x;
#endif

  foreach()
    foreach_dimension()
      momentum.x[] = dual_momentum_transport_density(point)*u.x[];
  boundary((scalar *){momentum});

  void (* sweep[dimension]) (vector);
  int direction = 0;
  foreach_dimension()
    sweep[direction++] = dual_momentum_sweep_x;
  for (direction = 0; direction < dimension; direction++)
    sweep[(i + direction) % dimension](momentum);

  long invalid_cells = 0;
  double maximum_cell_defect = 0.;
  foreach(reduction(+:invalid_cells) reduction(max:maximum_cell_defect)) {
    double defect = dual_momentum_indicator_partition_defect(
      lower_gas_fraction[], lower_envelope_fraction[],
      upper_gas_fraction[], upper_envelope_fraction[]);
    maximum_cell_defect = max(maximum_cell_defect, defect);
    if (defect > (double)DUAL_MOMENTUM_INDICATOR_TOLERANCE)
      invalid_cells++;
  }
  dual_momentum_invalid_cell_count = invalid_cells;
  dual_momentum_invalid_cell_count_max =
    max(dual_momentum_invalid_cell_count_max, invalid_cells);
  dual_momentum_max_cell_partition_defect =
    max(dual_momentum_max_cell_partition_defect, maximum_cell_defect);
  if (invalid_cells > 0)
    dual_momentum_transport_failed = 1;

  dual_momentum_saved_interfaces = interfaces;
  interfaces = NULL;
}

event tracer_advection (i++)
  interfaces = dual_momentum_saved_interfaces;

event dual_momentum_transport_guard (i++, last)
{
  if (dual_momentum_transport_failed) {
    if (pid() == 0)
      fprintf(stderr,
        "dual momentum transport rejected state at t=%g i=%d: "
        "invalid_cells=%ld invalid_faces=%ld max_cell_defect=%g "
        "max_face_defect=%g\n",
        t, i, dual_momentum_invalid_cell_count,
        dual_momentum_invalid_face_count,
        dual_momentum_max_cell_partition_defect,
        dual_momentum_max_face_partition_defect);
    return 1;
  }
}

#endif
