#ifndef AXISYMMETRIC_MATCHED_FILM_H
#define AXISYMMETRIC_MATCHED_FILM_H

#include <float.h>
#include <math.h>
#include <stdbool.h>

#ifndef AXISYMMETRIC_MATCHED_FILM_BINS
# define AXISYMMETRIC_MATCHED_FILM_BINS 32
#endif

typedef struct {
  double radius[AXISYMMETRIC_MATCHED_FILM_BINS];
  double width[AXISYMMETRIC_MATCHED_FILM_BINS];
  double gap[AXISYMMETRIC_MATCHED_FILM_BINS];
  double delta[AXISYMMETRIC_MATCHED_FILM_BINS];
  double lower_position[AXISYMMETRIC_MATCHED_FILM_BINS];
  double upper_position[AXISYMMETRIC_MATCHED_FILM_BINS];
  double relative_normal_velocity[AXISYMMETRIC_MATCHED_FILM_BINS];
  double lower_tangential_velocity[AXISYMMETRIC_MATCHED_FILM_BINS];
  double upper_tangential_velocity[AXISYMMETRIC_MATCHED_FILM_BINS];
} AxisymmetricFilmInput;

typedef struct {
  double water_viscosity;
  double match_cells;
  double release_factor;
  int bins;
} AxisymmetricMatchedFilmConfig;

typedef struct {
  bool active;
  double cumulative_work;
} AxisymmetricMatchedFilmState;

enum {
  AXISYMMETRIC_MATCHED_FILM_NOT_EVALUATED = 0,
  AXISYMMETRIC_MATCHED_FILM_INVALID_INPUT = 1,
  AXISYMMETRIC_MATCHED_FILM_RAW_SOLVE_FAILED = 2,
  AXISYMMETRIC_MATCHED_FILM_REFERENCE_SOLVE_FAILED = 3,
  AXISYMMETRIC_MATCHED_FILM_OK = 4
};
typedef int AxisymmetricMatchedFilmStatus;

typedef struct {
  bool valid;
  bool active;
  AxisymmetricMatchedFilmStatus status;
  double raw_pressure[AXISYMMETRIC_MATCHED_FILM_BINS];
  double reference_pressure[AXISYMMETRIC_MATCHED_FILM_BINS];
  double excess_pressure[AXISYMMETRIC_MATCHED_FILM_BINS];
  double ring_force[AXISYMMETRIC_MATCHED_FILM_BINS];
  double residual[AXISYMMETRIC_MATCHED_FILM_BINS];
  double pair_net_force;
  double pair_net_moment;
  int passivity_clipped_bins;
  double passivity_clipped_force;
  double step_work;
  double cumulative_work;
} AxisymmetricMatchedFilmResult;

typedef struct {
  bool valid;
  long double pressure[AXISYMMETRIC_MATCHED_FILM_BINS];
  double residual[AXISYMMETRIC_MATCHED_FILM_BINS];
  double face_flux[AXISYMMETRIC_MATCHED_FILM_BINS + 1];
  /* Maximum absolute ring balance error, in velocity units. */
  double maximum_continuity_balance_error;
  double maximum_pressure_balance_error;
  double maximum_flux_balance_error;
} AxisymmetricReynoldsSolution;

static inline double axisymmetric_matched_film_square (double value)
{
  return value*value;
}

static inline double axisymmetric_matched_film_cube (double value)
{
  return value*value*value;
}

static inline double axisymmetric_matched_film_max (double first, double second)
{
  return first > second ? first : second;
}

static inline bool axisymmetric_matched_film_input_is_valid (
  const AxisymmetricMatchedFilmConfig * config,
  const AxisymmetricFilmInput * input, double dt)
{
  if (!config || !input || !isfinite(dt) || dt <= 0. ||
      !isfinite(config->water_viscosity) || config->water_viscosity <= 0. ||
      !isfinite(config->match_cells) || config->match_cells <= 0. ||
      !isfinite(config->release_factor) || config->release_factor < 1. ||
      config->bins != AXISYMMETRIC_MATCHED_FILM_BINS)
    return false;
  for (int bin = 0; bin < AXISYMMETRIC_MATCHED_FILM_BINS; bin++)
    if (!isfinite(input->radius[bin]) || input->radius[bin] < 0. ||
        !isfinite(input->width[bin]) || input->width[bin] <= 0. ||
        !isfinite(input->gap[bin]) || input->gap[bin] <= 0. ||
        !isfinite(input->delta[bin]) || input->delta[bin] <= 0. ||
        !isfinite(input->lower_position[bin]) ||
        !isfinite(input->upper_position[bin]) ||
        input->upper_position[bin] <= input->lower_position[bin] ||
        !isfinite(input->relative_normal_velocity[bin]) ||
        !isfinite(input->lower_tangential_velocity[bin]) ||
        !isfinite(input->upper_tangential_velocity[bin]))
      return false;
  return true;
}

static inline double axisymmetric_reynolds_pressure_residuals (
  const AxisymmetricFilmInput * input,
  const double * pressure_face_coefficient, const double * couette_face,
  const long double * pressure, double * residual)
{
  const int last = AXISYMMETRIC_MATCHED_FILM_BINS - 1;
  double maximum_error = 0.;
  for (int bin = 0; bin < AXISYMMETRIC_MATCHED_FILM_BINS; bin++) {
    const long double inner_flux = bin > 0 ?
      (long double) pressure_face_coefficient[bin]*
        (pressure[bin - 1] - pressure[bin]) +
        (long double) couette_face[bin] : 0.L;
    const int outer_face = bin + 1;
    const long double outer_flux = bin < last ?
      (long double) pressure_face_coefficient[outer_face]*
        (pressure[bin] - pressure[bin + 1]) +
        (long double) couette_face[outer_face] :
      (long double) pressure_face_coefficient[outer_face]*pressure[bin] +
        (long double) couette_face[outer_face];
    const long double ring_area = (long double) input->radius[bin]*
      (long double) input->width[bin];
    const long double pressure_residual = outer_flux - inner_flux - ring_area*
      (long double) input->relative_normal_velocity[bin];
    residual[bin] = (double) pressure_residual;
    const double balance_error = (double) (fabsl(pressure_residual)/ring_area);
    if (!isfinite(inner_flux) || !isfinite(outer_flux) ||
        !isfinite(ring_area) || ring_area <= 0.L ||
        !isfinite(pressure_residual) || !isfinite(residual[bin]) ||
        !isfinite(balance_error))
      return NAN;
    maximum_error = axisymmetric_matched_film_max(maximum_error,
      balance_error);
  }
  return maximum_error;
}

static inline AxisymmetricReynoldsSolution axisymmetric_reynolds_solution (
  const AxisymmetricMatchedFilmConfig * config,
  const AxisymmetricFilmInput * input, const double * gap)
{
  enum { bins = AXISYMMETRIC_MATCHED_FILM_BINS };
  AxisymmetricReynoldsSolution solution = {0};
  if (!config || !input || !gap ||
      !isfinite(config->water_viscosity) || config->water_viscosity <= 0. ||
      config->bins != bins)
    return solution;
  double couette_face[bins + 1] = {0.};
  double conductivity[bins] = {0.};
  double pressure_face_coefficient[bins + 1] = {0.};
  double center_couette_flux[bins];
  for (int bin = 0; bin < bins; bin++) {
    if (!isfinite(input->radius[bin]) || input->radius[bin] < 0. ||
        !isfinite(input->width[bin]) || input->width[bin] <= 0. ||
        !isfinite(gap[bin]) || gap[bin] <= 0. ||
        !isfinite(input->relative_normal_velocity[bin]) ||
        !isfinite(input->lower_tangential_velocity[bin]) ||
        !isfinite(input->upper_tangential_velocity[bin]))
      return solution;
    double inner_radius = input->radius[bin] - .5*input->width[bin];
    double expected_inner_radius = bin == 0 ? 0. :
      input->radius[bin - 1] + .5*input->width[bin - 1];
    double geometry_scale = axisymmetric_matched_film_max(
      input->width[bin], bin == 0 ? input->width[bin] :
      input->width[bin - 1]);
    if (fabs(inner_radius - expected_inner_radius) >
        128.*DBL_EPSILON*geometry_scale ||
        fabs(input->width[bin] - input->width[0]) >
        128.*DBL_EPSILON*axisymmetric_matched_film_max(
          input->width[bin], input->width[0]))
      return solution;
    center_couette_flux[bin] = .5*gap[bin]*(
      input->lower_tangential_velocity[bin] +
      input->upper_tangential_velocity[bin]);
    conductivity[bin] = axisymmetric_matched_film_cube(gap[bin])/
      (12.*config->water_viscosity);
    if (!isfinite(center_couette_flux[bin]) ||
        !isfinite(conductivity[bin]) || conductivity[bin] <= 0.)
      return solution;
  }
  for (int face = 1; face < bins; face++) {
    if (input->radius[face] <= input->radius[face - 1])
      return solution;
    double radius_face = .5*(input->radius[face - 1] +
      input->radius[face]);
    couette_face[face] = radius_face*.5*(
      center_couette_flux[face - 1] + center_couette_flux[face]);
    double face_conductivity =
      2.*conductivity[face - 1]*conductivity[face]/
      (conductivity[face - 1] + conductivity[face]);
    pressure_face_coefficient[face] = radius_face*face_conductivity/
      (input->radius[face] - input->radius[face - 1]);
    if (!isfinite(couette_face[face]) ||
        !isfinite(pressure_face_coefficient[face]) ||
        pressure_face_coefficient[face] <= 0.)
      return solution;
  }
  couette_face[bins] =
    (input->radius[bins - 1] + .5*input->width[bins - 1])*
    center_couette_flux[bins - 1];
  if (!isfinite(couette_face[bins]))
    return solution;
  const double outer_radius = input->radius[bins - 1] +
    .5*input->width[bins - 1];
  const double outer_distance = .5*input->width[bins - 1];
  const double outer_coefficient = outer_radius*conductivity[bins - 1]/
    outer_distance;
  if (!isfinite(outer_radius) || !isfinite(outer_distance) ||
      !isfinite(outer_coefficient) || outer_radius <= 0. ||
      outer_distance <= 0. || outer_coefficient <= 0.)
    return solution;
  pressure_face_coefficient[bins] = outer_coefficient;

  /* Integrate the finite-volume continuity equation first.  Reconstructing
     these fluxes from differences of large pressures loses the conservative
     identity on strongly varying sampled gaps. */
  long double face_flux_long[bins + 1];
  face_flux_long[0] = 0.L;
  solution.face_flux[0] = 0.;
  long double cumulative_flux = 0.L;
  for (int bin = 0; bin < bins; bin++) {
    const long double inner_radius = (long double) input->radius[bin] -
      .5L*(long double) input->width[bin];
    const long double outer_radius = (long double) input->radius[bin] +
      .5L*(long double) input->width[bin];
    const long double ring_area = .5L*(outer_radius*outer_radius -
      inner_radius*inner_radius);
    cumulative_flux += ring_area*
      (long double) input->relative_normal_velocity[bin];
    if (!isfinite(cumulative_flux) ||
        fabsl(cumulative_flux) > (long double) DBL_MAX)
      return solution;
    face_flux_long[bin + 1] = cumulative_flux;
    solution.face_flux[bin + 1] = (double) cumulative_flux;
  }

  /* The outer pressure is zero.  The conservative face fluxes therefore give
     a stable one-pass pressure reconstruction through the same constitutive
     law used to assemble the Reynolds operator. */
  long double pressure_long[bins];
  if (!isfinite(pressure_face_coefficient[bins]) ||
      pressure_face_coefficient[bins] <= 0.)
    return solution;
  pressure_long[bins - 1] =
    (face_flux_long[bins] - (long double) couette_face[bins])/
    (long double) pressure_face_coefficient[bins];
  for (int bin = bins - 2; bin >= 0; bin--) {
    const int face = bin + 1;
    if (!isfinite(pressure_face_coefficient[face]) ||
        pressure_face_coefficient[face] <= 0. ||
        !isfinite(pressure_long[bin + 1]))
      return solution;
    pressure_long[bin] = pressure_long[bin + 1] +
      (face_flux_long[face] - (long double) couette_face[face])/
      (long double) pressure_face_coefficient[face];
  }
  for (int bin = 0; bin < bins; bin++) {
    if (!isfinite(pressure_long[bin]) ||
        fabsl(pressure_long[bin]) > (long double) DBL_MAX)
      return solution;
    solution.pressure[bin] = pressure_long[bin];
  }
  double maximum_pressure_error = axisymmetric_reynolds_pressure_residuals(
    input, pressure_face_coefficient, couette_face, solution.pressure,
    solution.residual);
  if (!isfinite(maximum_pressure_error))
    return solution;

  for (int bin = 0; bin < bins; bin++) {
    double inner_radius = input->radius[bin] - .5*input->width[bin];
    double outer_radius = input->radius[bin] + .5*input->width[bin];
    double ring_area = .5*(axisymmetric_matched_film_square(outer_radius) -
      axisymmetric_matched_film_square(inner_radius));
    double balance =
      (solution.face_flux[bin + 1] - solution.face_flux[bin])/ring_area;
    double balance_error =
      fabs(balance - input->relative_normal_velocity[bin]);
    if (!isfinite(solution.face_flux[bin + 1]) || !isfinite(balance) ||
        !isfinite(balance_error))
      return (AxisymmetricReynoldsSolution) {0};
    solution.maximum_continuity_balance_error = axisymmetric_matched_film_max(
      solution.maximum_continuity_balance_error, balance_error);
  }
  solution.maximum_pressure_balance_error = maximum_pressure_error;
  solution.maximum_flux_balance_error = axisymmetric_matched_film_max(
    solution.maximum_continuity_balance_error,
    solution.maximum_pressure_balance_error);
  solution.valid = true;
  return solution;
}

static inline bool axisymmetric_reynolds_pressure (
  const AxisymmetricMatchedFilmConfig * config,
  const AxisymmetricFilmInput * input, const double * gap,
  double * pressure, double * residual)
{
  if (!pressure || !residual)
    return false;
  AxisymmetricReynoldsSolution solution = axisymmetric_reynolds_solution(
    config, input, gap);
  if (!solution.valid)
    return false;
  for (int bin = 0; bin < AXISYMMETRIC_MATCHED_FILM_BINS; bin++) {
    pressure[bin] = (double) solution.pressure[bin];
    residual[bin] = solution.residual[bin];
  }
  return true;
}

static inline AxisymmetricMatchedFilmResult axisymmetric_matched_film_evaluate (
  const AxisymmetricMatchedFilmConfig * config,
  const AxisymmetricFilmInput * input,
  const AxisymmetricReynoldsSolution * raw_solution,
  double cumulative_work_before, double dt)
{
  AxisymmetricMatchedFilmResult result = {0};
  if (!axisymmetric_matched_film_input_is_valid(config, input, dt) ||
      !raw_solution || !isfinite(cumulative_work_before)) {
    result.status = AXISYMMETRIC_MATCHED_FILM_INVALID_INPUT;
    return result;
  }
  if (!raw_solution->valid ||
      !isfinite(raw_solution->maximum_flux_balance_error)) {
    result.status = AXISYMMETRIC_MATCHED_FILM_RAW_SOLVE_FAILED;
    return result;
  }
  for (int face = 0; face <= AXISYMMETRIC_MATCHED_FILM_BINS; face++)
    if (!isfinite(raw_solution->face_flux[face])) {
      result.status = AXISYMMETRIC_MATCHED_FILM_RAW_SOLVE_FAILED;
      return result;
    }
  for (int bin = 0; bin < AXISYMMETRIC_MATCHED_FILM_BINS; bin++)
    if (!isfinite(raw_solution->pressure[bin]) ||
        !isfinite(raw_solution->residual[bin])) {
      result.status = AXISYMMETRIC_MATCHED_FILM_RAW_SOLVE_FAILED;
      return result;
    }

  double effective_gap[AXISYMMETRIC_MATCHED_FILM_BINS];
  for (int bin = 0; bin < AXISYMMETRIC_MATCHED_FILM_BINS; bin++) {
    effective_gap[bin] = axisymmetric_matched_film_max(input->gap[bin],
      config->match_cells*input->delta[bin]);
    result.residual[bin] = raw_solution->residual[bin];
  }
  const AxisymmetricReynoldsSolution reference_solution =
    axisymmetric_reynolds_solution(config, input, effective_gap);
  if (!reference_solution.valid) {
    result.status = AXISYMMETRIC_MATCHED_FILM_REFERENCE_SOLVE_FAILED;
    return result;
  }

  const double two_pi = 2.*acos(-1.);
  for (int bin = 0; bin < AXISYMMETRIC_MATCHED_FILM_BINS; bin++) {
    result.raw_pressure[bin] = (double) raw_solution->pressure[bin];
    result.reference_pressure[bin] =
      (double) reference_solution.pressure[bin];
    double candidate_excess_pressure = 0.;
    if (input->gap[bin] < config->match_cells*input->delta[bin])
      candidate_excess_pressure = (double) fmaxl(0.L,
        fmaxl(0.L, raw_solution->pressure[bin]) -
        fmaxl(0.L, reference_solution.pressure[bin]));
    /* The reduced closure applies normal pressure but not film shear. */
    if (candidate_excess_pressure > 0. &&
        input->relative_normal_velocity[bin] <= 0.) {
      result.passivity_clipped_bins++;
      result.passivity_clipped_force +=
        two_pi*input->radius[bin]*input->width[bin]*
        candidate_excess_pressure;
    }
    else
      result.excess_pressure[bin] = candidate_excess_pressure;
    result.ring_force[bin] = two_pi*input->radius[bin]*input->width[bin]*
      result.excess_pressure[bin];
    result.step_work -= result.ring_force[bin]*
      input->relative_normal_velocity[bin]*dt;
  }
  result.cumulative_work = cumulative_work_before + result.step_work;
  if (!isfinite(result.step_work) || !isfinite(result.cumulative_work)) {
    result.status = AXISYMMETRIC_MATCHED_FILM_INVALID_INPUT;
    return result;
  }
  result.valid = true;
  result.active = true;
  result.status = AXISYMMETRIC_MATCHED_FILM_OK;
  return result;
}

static inline AxisymmetricMatchedFilmResult axisymmetric_matched_film_update (
  const AxisymmetricMatchedFilmConfig * config,
  AxisymmetricMatchedFilmState * state,
  const AxisymmetricFilmInput * input, double dt)
{
  AxisymmetricMatchedFilmResult result = {0};
  if (!state || !axisymmetric_matched_film_input_is_valid(config, input, dt)) {
    result.status = AXISYMMETRIC_MATCHED_FILM_INVALID_INPUT;
    return result;
  }

  bool enters_match = false;
  bool leaves_match = true;
  for (int bin = 0; bin < AXISYMMETRIC_MATCHED_FILM_BINS; bin++) {
    double match_gap = config->match_cells*input->delta[bin];
    enters_match |= input->gap[bin] < match_gap;
    leaves_match &= input->gap[bin] > config->release_factor*match_gap;
  }
  if (enters_match)
    state->active = true;
  else if (leaves_match)
    state->active = false;

  result.valid = true;
  result.active = state->active;
  result.status = AXISYMMETRIC_MATCHED_FILM_OK;
  if (!state->active) {
    result.cumulative_work = state->cumulative_work;
    return result;
  }

  AxisymmetricReynoldsSolution raw_solution = axisymmetric_reynolds_solution(
    config, input, input->gap);
  result = axisymmetric_matched_film_evaluate(config, input, &raw_solution,
    state->cumulative_work, dt);
  if (result.valid)
    state->cumulative_work = result.cumulative_work;
  return result;
}

#endif
