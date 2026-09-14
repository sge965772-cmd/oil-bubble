#include "grid/quadtree.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "axisymmetric_numerical_policy.h"
#include "dual_compound_phases.h"
#ifndef MOMENTUM_TRANSPORT_VERSION
# define MOMENTUM_TRANSPORT_VERSION 3
#endif
#include "positive_material_momentum_policy.h"
#include "tension.h"
#include "separated_solver_tolerance_policy.h"
#include "precursor_state.h"

#ifndef MAXLEVEL
# define MAXLEVEL 10
#endif
#ifndef MINLEVEL
# define MINLEVEL 7
#endif
#ifndef DOMAIN_LENGTH
# define DOMAIN_LENGTH 0.080
#endif
#ifndef AXIAL_ORIGIN
# define AXIAL_ORIGIN -0.030
#endif
#ifndef INITIAL_CENTER
# define INITIAL_CENTER -0.020
#endif
#ifndef OUTER_RADIUS
# define OUTER_RADIUS 0.002
#endif
#ifndef OIL_VOLUME_FRACTION
# define OIL_VOLUME_FRACTION 0.20
#endif
#ifndef SIGMA_OW
# define SIGMA_OW 40.9e-3
#endif
#ifndef SIGMA_OG
# define SIGMA_OG 19.1e-3
#endif
#ifndef GRAVITY_MAGNITUDE
# define GRAVITY_MAGNITUDE 9.81
#endif
#ifndef GRAVITY_RAMP_TIME
# define GRAVITY_RAMP_TIME 0.
#endif
#ifndef END_TIME
# define END_TIME 0.080
#endif
#ifndef DT_MAX
# define DT_MAX 5.e-7
#endif
#ifndef OUTPUT_INTERVAL
# define OUTPUT_INTERVAL 1.e-5
#endif
#ifndef DUMP_INTERVAL
# define DUMP_INTERVAL 1.e-3
#endif
#ifndef FRACTION_TOLERANCE
# define FRACTION_TOLERANCE 1.e-6
#endif
#ifndef VELOCITY_TOLERANCE
# define VELOCITY_TOLERANCE 1.e-4
#endif
#ifndef RESTART_FIELD_AUDIT
# define RESTART_FIELD_AUDIT 0
#endif
static double inner_radius;
static double initial_gas_volume;
static double initial_oil_volume;
static const char * restart_path = NULL;
#if RESTART_FIELD_AUDIT
static int restart_field_audit_count = 0;
static int restart_projection_audit_count = 0;
#endif
static int restart_pressure_boundary_initialized = 0;
static double axial_gravity_boundary_magnitude =
  (double)GRAVITY_RAMP_TIME > 0. ? 0. : (double)GRAVITY_MAGNITUDE;

static void initialize_lower_bubble (void)
{
  fraction(lower_envelope_fraction,
           sq((double)OUTER_RADIUS) -
           sq(x - (double)INITIAL_CENTER) - sq(y));
  fraction(lower_gas_fraction,
           sq(inner_radius) - sq(x - (double)INITIAL_CENTER) - sq(y));
  foreach() {
    upper_gas_fraction[] = 0.;
    upper_envelope_fraction[] = 0.;
  }
  boundary((scalar *){
    lower_gas_fraction, lower_envelope_fraction,
    upper_gas_fraction, upper_envelope_fraction
  });
}

int main (int argc, char ** argv)
{
  if (argc > 2) {
    fprintf(stderr, "usage: %s [RESTART_DUMP]\n", argv[0]);
    return 2;
  }
  if (argc == 2)
    restart_path = argv[1];
  inner_radius = (double)OUTER_RADIUS*
    pow(1. - (double)OIL_VOLUME_FRACTION, 1./3.);
  size((double)DOMAIN_LENGTH);
  origin((double)AXIAL_ORIGIN, 0.);
  init_grid(1 << MINLEVEL);

  DT = (double)DT_MAX;
  CFL = 0.20;
  TOLERANCE = (double)PRESSURE_TOLERANCE;
  NITERMAX = 10000;

  lower_gas_fraction.sigma = (double)SIGMA_OG;
  lower_envelope_fraction.sigma = (double)SIGMA_OW;
  upper_gas_fraction.sigma = (double)SIGMA_OG;
  upper_envelope_fraction.sigma = (double)SIGMA_OW;
  p.nodump = false;
  pf.nodump = false;
  run();
  return dual_momentum_transport_failed ? 3 : 0;
}

u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);
p[left] = neumann(rho_water*axial_gravity_boundary_magnitude);
p[right] = neumann(-rho_water*axial_gravity_boundary_magnitude);
p[top] = neumann(0.);

event init (i = 0)
{
  /* centered.h marks both pressures nodump in its defaults event. */
  p.nodump = false;
  pf.nodump = false;
  if (restart_path) {
    if (!restore(file = restart_path)) {
      if (pid() == 0)
        fprintf(stderr, "precursor restart failed: %s\n", restart_path);
      exit(2);
    }
  }
  else {
    initialize_lower_bubble();
#if TREE
    for (int pass = 0; pass < MAXLEVEL - MINLEVEL + 1; pass++) {
      adapt_wavelet((scalar *){
          lower_gas_fraction, lower_envelope_fraction
        },
        (double[]){
          (double)FRACTION_TOLERANCE, (double)FRACTION_TOLERANCE
        }, MAXLEVEL, MINLEVEL);
      initialize_lower_bubble();
    }
#endif
    foreach() {
      u.x[] = u.y[] = 0.;
      p[] = -(double)rho_water*(double)GRAVITY_MAGNITUDE*x +
            2.*(double)SIGMA_OW/(double)OUTER_RADIUS*
              lower_envelope_fraction[] +
            2.*(double)SIGMA_OG/inner_radius*lower_gas_fraction[];
      pf[] = p[];
    }
  }
  if (restart_path) {
    boundary((scalar *){
      lower_gas_fraction, lower_envelope_fraction,
      upper_gas_fraction, upper_envelope_fraction, u
    });
    /* Face properties are derived fields and are not stored in dumps. */
    dual_compound_update_properties();
  }
  else
    boundary((scalar *){
      lower_gas_fraction, lower_envelope_fraction,
      upper_gas_fraction, upper_envelope_fraction, u, p, pf
    });
  PrecursorState state = measure_precursor_state();
  initial_gas_volume = state.lower_gas_volume;
  initial_oil_volume = state.lower_oil_volume;

  if (pid() == 0) {
    if (restart_path)
      fprintf(stderr,
        "# PRECURSOR_RESTART source=%s init_event_t=%.17g init_event_i=%d "
        "target_t=%.17g\n",
        restart_path, t, i, (double)END_TIME);
    FILE * fp = fopen("data/precursor_stats.dat", "w");
    if (!fp) {
      perror("data/precursor_stats.dat");
      exit(2);
    }
    fprintf(fp,
      "# t i dt cells gas_volume oil_volume gas_drift oil_drift center "
      "rise_distance center_velocity moment_aspect equivalent_diameter "
      "Re We Bo kinetic_energy pressure_integral "
      "projected_pressure_integral nesting_volume envelope_overlap_volume "
      "gas_overlap_volume overfill_volume maximum_speed "
      "projected_volume_change pressure_iterations scaled_pressure_residual "
      "viscosity_iterations viscosity_residual momentum_transport_mode "
      "invalid_cells invalid_faces transport_failed "
      "max_cell_partition_defect max_face_partition_defect "
      "max_vof_bound_defect max_partition_correction "
      "limited_faces min_flux_alpha max_simplex_closure_defect "
      "max_flux_fraction_correction max_flux_courant_correction "
      "nested_plic_fallbacks nested_plic_fallbacks_max "
      "transport_cpu_seconds\n");
    fclose(fp);
    fp = fopen("data/dump_catalog.csv", "w");
    fprintf(fp,
      "time,iteration,center,center_velocity,moment_aspect,gas_volume,"
      "oil_volume,kinetic_energy,pressure_integral,"
      "projected_pressure_integral,cells,dump\n");
    fclose(fp);
  }
}

event acceleration (i++)
{
  double ramp = (double)GRAVITY_RAMP_TIME > 0. ?
                min(t/(double)GRAVITY_RAMP_TIME, 1.) : 1.;
  axial_gravity_boundary_magnitude = (double)GRAVITY_MAGNITUDE*ramp;
  face vector acceleration_field = a;
  foreach_face(x)
    acceleration_field.x[] += -(double)GRAVITY_MAGNITUDE*ramp;
#if RESTART_FIELD_AUDIT
  if (restart_path && restart_field_audit_count < 2) {
    long nonfinite_cells = 0, nonfinite_faces = 0;
    double maximum_speed = 0., maximum_pressure = 0.;
    double minimum_alpha = HUGE, maximum_alpha = 0.;
    foreach(reduction(+:nonfinite_cells) reduction(max:maximum_speed)
            reduction(max:maximum_pressure)) {
      if (!isfinite(u.x[]) || !isfinite(u.y[]) || !isfinite(p[]) ||
          !isfinite(pf[]) || !isfinite(lower_gas_fraction[]) ||
          !isfinite(lower_envelope_fraction[]) ||
          !isfinite(upper_gas_fraction[]) ||
          !isfinite(upper_envelope_fraction[]))
        nonfinite_cells++;
      maximum_speed = max(maximum_speed, hypot(u.x[], u.y[]));
      maximum_pressure = max(maximum_pressure, fabs(p[]));
    }
    foreach_face(reduction(+:nonfinite_faces) reduction(min:minimum_alpha)
                 reduction(max:maximum_alpha)) {
      if (!isfinite(alpha.x[]) || !isfinite(acceleration_field.x[]) ||
          !isfinite(fm.x[]))
        nonfinite_faces++;
      if (fm.x[] > 0.)
        minimum_alpha = min(minimum_alpha, alpha.x[]);
      maximum_alpha = max(maximum_alpha, alpha.x[]);
    }
    if (pid() == 0)
      fprintf(stderr,
        "[RESTART-FIELD-AUDIT] t=%.17g i=%d cells_nonfinite=%ld "
        "faces_nonfinite=%ld umax=%.17g pabs_max=%.17g "
        "alpha_min=%.17g alpha_max=%.17g\n",
        t, i, nonfinite_cells, nonfinite_faces, maximum_speed,
        maximum_pressure, minimum_alpha, maximum_alpha), fflush(stderr);
    restart_field_audit_count++;
  }
#endif
}

event projection (i++)
{
  if (restart_path && !restart_pressure_boundary_initialized) {
    boundary((scalar *){p, pf});
    restart_pressure_boundary_initialized = 1;
  }
#if RESTART_FIELD_AUDIT
  if (restart_path && restart_projection_audit_count < 2) {
    long nonfinite_face_velocity = 0, nonfinite_pressure_gradient = 0;
    long nonfinite_divergence = 0;
    double maximum_face_velocity = 0., maximum_divergence = 0.;
    foreach_face(reduction(+:nonfinite_face_velocity)
                 reduction(+:nonfinite_pressure_gradient)
                 reduction(max:maximum_face_velocity)) {
      double gradient = (p[] - p[-1])/Delta;
      if (!isfinite(uf.x[]))
        nonfinite_face_velocity++;
      if (!isfinite(gradient)) {
        if (nonfinite_pressure_gradient < 2)
          fprintf(stderr,
            "[RESTART-PGRAD-REJECT] pid=%d t=%.17g i=%d "
            "x=%.17g y=%.17g Delta=%.17g p_left=%.17g "
            "p_right=%.17g alpha=%.17g fm=%.17g\n",
            pid(), t, i, x, y, Delta, p[-1], p[], alpha.x[], fm.x[]),
            fflush(stderr);
        nonfinite_pressure_gradient++;
      }
      maximum_face_velocity = max(maximum_face_velocity, fabs(uf.x[]));
    }
    foreach(reduction(+:nonfinite_divergence)
            reduction(max:maximum_divergence)) {
      double divergence = 0.;
      foreach_dimension()
        divergence += uf.x[1] - uf.x[];
      divergence /= dt*Delta;
      if (!isfinite(divergence))
        nonfinite_divergence++;
      maximum_divergence = max(maximum_divergence, fabs(divergence));
    }
    if (pid() == 0)
      fprintf(stderr,
        "[RESTART-PROJECTION-AUDIT] t=%.17g i=%d "
        "uf_nonfinite=%ld pgrad_nonfinite=%ld div_nonfinite=%ld "
        "ufmax=%.17g divmax=%.17g\n",
        t, i, nonfinite_face_velocity, nonfinite_pressure_gradient,
        nonfinite_divergence, maximum_face_velocity,
        maximum_divergence), fflush(stderr);
    restart_projection_audit_count++;
  }
#endif
}

event diagnostics (t = 0.; t <= END_TIME; t += OUTPUT_INTERVAL)
{
  PrecursorState state = measure_precursor_state();
  double gas_drift = fabs(state.lower_gas_volume - initial_gas_volume)/
                     max(initial_gas_volume, 1.e-300);
  double oil_drift = fabs(state.lower_oil_volume - initial_oil_volume)/
                     max(initial_oil_volume, 1.e-300);
  double reynolds = rho_water*fabs(state.center_velocity)*
                    state.equivalent_diameter/max(mu_water, 1.e-300);
  double weber = rho_water*sq(state.center_velocity)*
                 state.equivalent_diameter/max((double)SIGMA_OW, 1.e-300);
  double bond = (rho_water - rho_gas)*(double)GRAVITY_MAGNITUDE*
                sq(state.equivalent_diameter)/max((double)SIGMA_OW, 1.e-300);
  if (pid() == 0) {
    FILE * fp = fopen("data/precursor_stats.dat", "a");
    fprintf(fp,
      "%.17g %d %.17g %ld %.17g %.17g %.17g %.17g %.17g "
      "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g "
      "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %d %.17g ",
      t, i, dt, grid->tn, state.lower_gas_volume, state.lower_oil_volume,
      gas_drift, oil_drift, state.center,
      state.center - (double)INITIAL_CENTER, state.center_velocity,
      state.moment_aspect, state.equivalent_diameter,
      reynolds, weber, bond, state.kinetic_energy,
      state.pressure_integral, state.projected_pressure_integral,
      state.nesting_volume, state.envelope_overlap_volume,
      state.gas_overlap_volume, state.overfill_volume,
      state.maximum_speed, state.projected_volume_change,
      mgp.i, mgp.resa*sq(dt));
    fprintf(fp,
      "%d %.17g %d %ld %ld %d %.17g %.17g %.17g %.17g "
      "%ld %.17g %.17g %.17g %.17g %ld %ld %.17g\n",
      mgu.i, mgu.resa, (int)MOMENTUM_TRANSPORT_MODE,
      dual_momentum_invalid_cell_count,
      dual_momentum_invalid_face_count,
      dual_momentum_transport_failed,
      dual_momentum_max_cell_partition_defect,
      dual_momentum_max_face_partition_defect,
      dual_momentum_max_vof_bound_defect,
      dual_momentum_max_partition_correction,
      dual_momentum_limited_face_count,
      dual_momentum_min_flux_alpha,
      dual_momentum_max_simplex_closure_defect,
      dual_momentum_max_flux_fraction_correction,
      dual_momentum_max_flux_courant_correction,
      dual_momentum_nested_plic_fallback_count,
      dual_momentum_nested_plic_fallback_count_max,
      dual_momentum_last_transport_cpu_seconds);
    fclose(fp);
  }
}

event snapshots (t = 0.; t <= END_TIME; t += DUMP_INTERVAL)
{
  char path[128];
  sprintf(path, "dumps/dump-%09.6f", t);
  dump(file = path);
  PrecursorState state = measure_precursor_state();
  if (pid() == 0) {
    FILE * fp = fopen("data/dump_catalog.csv", "a");
    fprintf(fp,
      "%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%ld,%s\n",
      t, i, state.center, state.center_velocity, state.moment_aspect,
      state.lower_gas_volume, state.lower_oil_volume, state.kinetic_energy,
      state.pressure_integral, state.projected_pressure_integral,
      grid->tn, path);
    fclose(fp);
  }
}

event adapt (i++)
{
#if TREE
  adapt_wavelet((scalar *){
      lower_gas_fraction, lower_envelope_fraction, u
    },
    (double[]){
      (double)FRACTION_TOLERANCE, (double)FRACTION_TOLERANCE,
      (double)VELOCITY_TOLERANCE, (double)VELOCITY_TOLERANCE
    }, MAXLEVEL, MINLEVEL);
#endif
}

event finish (t = END_TIME)
{
  if (pid() == 0)
    fprintf(stderr, "# FINISH t=%.17g i=%d\n", t, i);
}
