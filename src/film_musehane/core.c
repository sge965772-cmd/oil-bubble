#include "film_musehane/core.h"
#include "film_musehane/published_equation23_operator.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static bool config_is_valid (const MusehaneConfig *config)
{
  return config && isfinite(config->water_density) &&
    config->water_density > 0. && isfinite(config->water_viscosity) &&
    config->water_viscosity > 0. &&
    isfinite(config->nonlinear_relative_tolerance) &&
    config->nonlinear_relative_tolerance > 0. &&
    isfinite(config->inventory_relative_tolerance) &&
    config->inventory_relative_tolerance > 0. &&
    config->maximum_iterations > 0 &&
    (config->discretization ==
       MUSEHANE_DISCRETIZATION_CONSERVATIVE_HARMONIC ||
     config->discretization ==
       MUSEHANE_DISCRETIZATION_PUBLISHED_EQUATION23);
}

static bool input_is_valid (
  const MusehaneConfig *config, const MusehaneStepInput *input,
  const MusehaneStepOutput *output)
{
  if (!input || !output || !input->mesh ||
      !isfinite(input->time) || input->time < 0. ||
      !isfinite(input->dt) || input->dt <= 0. ||
      input->mesh->cell_count == 0 ||
      input->mesh->cell_count > SIZE_MAX/sizeof(double) ||
      !input->mesh->cell_radius || !input->mesh->cell_width ||
      !input->mesh->cell_measure || !input->mesh->face_measure ||
      !input->active || !input->estimated_thickness ||
      !input->old_thickness ||
      !input->mean_surface_velocity_divergence ||
      !input->decontaminated_pressure_gradient_face ||
      !output->thickness || !output->pressure_flux_integral_face)
    return false;

  if (config->discretization ==
        MUSEHANE_DISCRETIZATION_PUBLISHED_EQUATION23 &&
      !input->decontaminated_pressure_gradient_cell)
    return false;

  const size_t cells = input->mesh->cell_count;
  for (size_t face = 0; face <= cells; face++)
    if (!isfinite(input->mesh->face_measure[face]) ||
        input->mesh->face_measure[face] < 0. ||
        (face > 0 && input->mesh->face_measure[face] <
          input->mesh->face_measure[face - 1]) ||
        !isfinite(input->decontaminated_pressure_gradient_face[face]))
      return false;
  for (size_t cell = 0; cell < cells; cell++)
    if (input->active[cell] > 1 ||
        !isfinite(input->mesh->cell_radius[cell]) ||
        input->mesh->cell_radius[cell] < 0. ||
        !isfinite(input->mesh->cell_width[cell]) ||
        input->mesh->cell_width[cell] <= 0. ||
        !isfinite(input->mesh->cell_measure[cell]) ||
        input->mesh->cell_measure[cell] <= 0. ||
        !isfinite(input->estimated_thickness[cell]) ||
        input->estimated_thickness[cell] <= 0. ||
        !isfinite(input->old_thickness[cell]) ||
        input->old_thickness[cell] <= 0. ||
        (config->discretization ==
           MUSEHANE_DISCRETIZATION_PUBLISHED_EQUATION23 &&
         !isfinite(input->decontaminated_pressure_gradient_cell[cell])) ||
        !isfinite(input->mean_surface_velocity_divergence[cell]))
      return false;
  return true;
}

static void face_mobility (
  const MusehaneConfig *config, const MusehaneStepInput *input,
  const double *thickness, size_t face,
  double *mobility, double *derivative_left,
  double *derivative_right)
{
  const bool has_left = face > 0 && input->active[face - 1] != 0;
  const bool has_right = face < input->mesh->cell_count &&
    input->active[face] != 0;
  const double factor = 1./(12.*config->water_viscosity);
  const double left_h = has_left ? thickness[face - 1] : 0.;
  const double right_h = has_right ? thickness[face] : 0.;
  const double left = factor*left_h*left_h*left_h;
  const double right = factor*right_h*right_h*right_h;
  *mobility = *derivative_left = *derivative_right = 0.;
  if (has_left && has_right) {
    const double sum = left + right;
    if (sum > 0.) {
      const double left_share = left/sum;
      const double right_share = right/sum;
      *mobility = 2.*left*right_share;
      *derivative_left = 6.*factor*left_h*left_h*
        right_share*right_share;
      *derivative_right = 6.*factor*right_h*right_h*
        left_share*left_share;
    }
  }
  else if (has_left) {
    *mobility = left;
    *derivative_left = 3.*factor*left_h*left_h;
  }
  else if (has_right) {
    *mobility = right;
    *derivative_right = 3.*factor*right_h*right_h;
  }
}

static void evaluate_pressure_flux (
  const MusehaneConfig *config, const MusehaneStepInput *input,
  const double *thickness, double *flux)
{
  const size_t cells = input->mesh->cell_count;
  for (size_t face = 0; face <= cells; face++) {
    double mobility, derivative_left, derivative_right;
    face_mobility(config, input, thickness, face,
      &mobility, &derivative_left, &derivative_right);
    flux[face] = input->mesh->face_measure[face]*mobility*
      input->decontaminated_pressure_gradient_face[face];
  }
}

static double active_inventory (
  const MusehaneStepInput *input, const double *thickness)
{
  double inventory = 0.;
  for (size_t cell = 0; cell < input->mesh->cell_count; cell++)
    if (input->active[cell])
      inventory += input->mesh->cell_measure[cell]*thickness[cell];
  return inventory;
}

typedef struct {
  double before, after, pressure_change, kinematic_change;
  double expected_after, closure_error, scale;
} MusehaneInventory;

static MusehaneInventory step_inventory (
  const MusehaneStepInput *input, const double *thickness,
  const double *pressure_rate)
{
  MusehaneInventory result = {
    .before = active_inventory(input, input->old_thickness),
    .after = active_inventory(input, thickness)
  };
  for (size_t cell = 0; cell < input->mesh->cell_count; cell++) {
    if (!input->active[cell])
      continue;
    result.pressure_change += input->dt*pressure_rate[cell];
    result.kinematic_change -= input->dt*
      input->mesh->cell_measure[cell]*thickness[cell]*
      input->mean_surface_velocity_divergence[cell];
  }
  result.expected_after = result.before + result.pressure_change +
    result.kinematic_change;
  result.closure_error = result.after - result.expected_after;
  result.scale = fmax(DBL_MIN,
    fmax(fabs(result.before), fabs(result.expected_after)));
  return result;
}

static bool has_active_region (const MusehaneStepInput *input)
{
  for (size_t cell = 0; cell < input->mesh->cell_count; cell++)
    if (input->active[cell])
      return true;
  return false;
}

static double equation_residual (
  const MusehaneConfig *config, const MusehaneStepInput *input,
  const double *thickness, double *flux, double *pressure_rate,
  double *pressure_lower, double *pressure_diagonal,
  double *pressure_upper, double *residual)
{
  if (config->discretization ==
      MUSEHANE_DISCRETIZATION_PUBLISHED_EQUATION23) {
    const MusehanePublishedEquation23Input operator_input = {
      .mesh = input->mesh,
      .active = input->active,
      .thickness = thickness,
      .pressure_gradient_cell =
        input->decontaminated_pressure_gradient_cell,
      .pressure_gradient_face =
        input->decontaminated_pressure_gradient_face,
      .viscosity = config->water_viscosity
    };
    MusehanePublishedEquation23Output operator_output = {
      .pressure_inventory_rate = pressure_rate,
      .pressure_flux_integral_face = flux,
      .lower_derivative = pressure_lower,
      .diagonal_derivative = pressure_diagonal,
      .upper_derivative = pressure_upper
    };
    if (musehane_published_equation23_pressure(
          &operator_input, &operator_output) != MUSEHANE_OK)
      return INFINITY;
  }
  else {
    evaluate_pressure_flux(config, input, thickness, flux);
    for (size_t cell = 0; cell < input->mesh->cell_count; cell++)
      pressure_rate[cell] = input->active[cell] ?
        flux[cell + 1] - flux[cell] : 0.;
  }
  double maximum = 0.;
  for (size_t cell = 0; cell < input->mesh->cell_count; cell++) {
    if (!input->active[cell]) {
      if (residual)
        residual[cell] = 0.;
      continue;
    }
    const double measure = input->mesh->cell_measure[cell];
    const double diagonal = measure*(1./input->dt +
      input->mean_surface_velocity_divergence[cell]);
    const double old_term = measure*input->old_thickness[cell]/input->dt;
    const double balance = diagonal*thickness[cell] - old_term -
      pressure_rate[cell];
    const double scale = fmax(DBL_MIN,
      fabs(diagonal*thickness[cell]) + fabs(old_term) +
      fabs(pressure_rate[cell]));
    if (residual)
      residual[cell] = balance;
    if (!isfinite(balance) || !isfinite(scale))
      return INFINITY;
    maximum = fmax(maximum, fabs(balance)/scale);
  }
  return maximum;
}

static void equation_residual_scale (
  const MusehaneStepInput *input, const double *thickness,
  const double *pressure_rate, double *scale)
{
  for (size_t cell = 0; cell < input->mesh->cell_count; cell++) {
    if (!input->active[cell]) {
      scale[cell] = DBL_MIN;
      continue;
    }
    const double measure = input->mesh->cell_measure[cell];
    const double diagonal = measure*(1./input->dt +
      input->mean_surface_velocity_divergence[cell]);
    const double old_term = measure*input->old_thickness[cell]/input->dt;
    scale[cell] = fmax(DBL_MIN,
      fabs(diagonal*thickness[cell]) + fabs(old_term) +
      fabs(pressure_rate[cell]));
  }
}

static double fixed_scale_residual (
  const MusehaneStepInput *input, const double *residual,
  const double *scale)
{
  double maximum = 0.;
  for (size_t cell = 0; cell < input->mesh->cell_count; cell++) {
    if (!input->active[cell])
      continue;
    const double relative = fabs(residual[cell])/scale[cell];
    if (!isfinite(relative))
      return INFINITY;
    maximum = fmax(maximum, relative);
  }
  return maximum;
}

MusehaneStatus musehane_surface_step_diagnosed (
  const MusehaneConfig *config,
  const MusehaneStepInput *input,
  MusehaneStepOutput *output,
  MusehaneCoreFailure *failure)
{
  if (failure)
    *failure = (MusehaneCoreFailure) {0};
  if (!config_is_valid(config))
    return MUSEHANE_INVALID_CONFIG;
  if (!input_is_valid(config, input, output))
    return MUSEHANE_INVALID_INPUT;
  if (!has_active_region(input))
    return MUSEHANE_NO_ACTIVE_REGION;

  const size_t cells = input->mesh->cell_count;
  double *current = malloc(cells*sizeof *current);
  double *next = malloc(cells*sizeof *next);
  double *flux = malloc((cells + 1)*sizeof *flux);
  double *trial_flux = malloc((cells + 1)*sizeof *trial_flux);
  double *pressure_rate = malloc(cells*sizeof *pressure_rate);
  double *trial_pressure_rate = malloc(cells*sizeof *trial_pressure_rate);
  double *pressure_lower = malloc(cells*sizeof *pressure_lower);
  double *pressure_diagonal = malloc(cells*sizeof *pressure_diagonal);
  double *pressure_upper = malloc(cells*sizeof *pressure_upper);
  double *lower = malloc(cells*sizeof *lower);
  double *diagonal = malloc(cells*sizeof *diagonal);
  double *upper = malloc(cells*sizeof *upper);
  double *rhs = malloc(cells*sizeof *rhs);
  double *delta = malloc(cells*sizeof *delta);
  double *trial_residual = malloc(cells*sizeof *trial_residual);
  double *line_search_scale = malloc(cells*sizeof *line_search_scale);
  if (!current || !next || !flux || !trial_flux || !pressure_rate ||
      !trial_pressure_rate || !pressure_lower || !pressure_diagonal ||
      !pressure_upper || !lower ||
      !diagonal || !upper || !rhs || !delta || !trial_residual ||
      !line_search_scale) {
    free(current);
    free(next);
    free(flux);
    free(trial_flux);
    free(pressure_rate);
    free(trial_pressure_rate);
    free(pressure_lower);
    free(pressure_diagonal);
    free(pressure_upper);
    free(lower);
    free(diagonal);
    free(upper);
    free(rhs);
    free(delta);
    free(trial_residual);
    free(line_search_scale);
    return MUSEHANE_ALLOCATION_FAILED;
  }
  memcpy(current, input->old_thickness, cells*sizeof *current);
  for (size_t cell = 0; cell < cells; cell++)
    if (!input->active[cell])
      current[cell] = input->estimated_thickness[cell];

  bool converged = false;
  unsigned iteration = 0;
  unsigned failure_backtracks = 0;
  MusehaneStatus status = MUSEHANE_NONLINEAR_NOT_CONVERGED;
  MusehaneCoreFailureKind failure_kind = MUSEHANE_CORE_FAILURE_NONE;
  const double initial_merit = equation_residual(config, input, current,
    flux, pressure_rate, pressure_lower, pressure_diagonal,
    pressure_upper, NULL);
  if (!isfinite(initial_merit)) {
    failure_kind = MUSEHANE_CORE_FAILURE_INITIAL_RESIDUAL_NONFINITE;
    goto cleanup;
  }
  for (size_t cell = 0; cell < cells; cell++) {
    if (!input->active[cell])
      continue;
    const double measure = input->mesh->cell_measure[cell];
    const double denominator = measure*(1./input->dt +
      input->mean_surface_velocity_divergence[cell]);
    if (!isfinite(denominator) || denominator <= 0.) {
      if (failure)
        *failure = (MusehaneCoreFailure) {
          .valid = true,
          .kind = MUSEHANE_CORE_FAILURE_NONPOSITIVE_DENOMINATOR,
          .cell_index = cell,
          .iteration = 1,
          .initial_merit = initial_merit,
          .final_merit = initial_merit,
          .numerator = measure*input->old_thickness[cell]/input->dt +
            pressure_rate[cell],
          .denominator = denominator,
          .old_thickness = input->old_thickness[cell],
          .current_thickness = current[cell],
          .flux_in = flux[cell],
          .flux_out = flux[cell + 1],
          .surface_velocity_divergence =
            input->mean_surface_velocity_divergence[cell]
        };
      status = MUSEHANE_NONPOSITIVE_THICKNESS;
      goto cleanup;
    }
  }
  /* Both pressure discretizations couple adjacent rings. Solve their
     backward-Euler equations together and backtrack whenever a trial
     thickness is not positive or does not reduce the discrete residual. */
  for (iteration = 1; iteration <= config->maximum_iterations; iteration++) {
    const double merit = equation_residual(config, input, current,
      flux, pressure_rate, pressure_lower, pressure_diagonal,
      pressure_upper, rhs);
    if (!isfinite(merit)) {
      failure_kind = MUSEHANE_CORE_FAILURE_ITERATION_RESIDUAL_NONFINITE;
      goto cleanup;
    }
    if (merit <= config->nonlinear_relative_tolerance) {
      const MusehaneInventory ledger =
        step_inventory(input, current, pressure_rate);
      if (isfinite(ledger.after) && isfinite(ledger.expected_after) &&
          isfinite(ledger.closure_error) &&
          fabs(ledger.closure_error) <=
            config->inventory_relative_tolerance*ledger.scale) {
        converged = true;
        break;
      }
    }
    /* The convergence scale depends on the trial thickness through both the
       accumulation and pressure terms. Freeze it within this Newton
       iteration so the line search measures descent of the residual rather
       than changes in its normalization. */
    equation_residual_scale(input, current, pressure_rate,
      line_search_scale);
    for (size_t cell = 0; cell < cells; cell++) {
      if (!input->active[cell]) {
        lower[cell] = upper[cell] = 0.;
        diagonal[cell] = 1.;
        rhs[cell] = 0.;
        continue;
      }
      const double measure = input->mesh->cell_measure[cell];
      if (config->discretization ==
          MUSEHANE_DISCRETIZATION_PUBLISHED_EQUATION23) {
        lower[cell] = -pressure_lower[cell];
        upper[cell] = -pressure_upper[cell];
        diagonal[cell] = measure*(1./input->dt +
          input->mean_surface_velocity_divergence[cell]) -
          pressure_diagonal[cell];
      }
      else {
        double inward_mobility, inward_left, inward_right;
        double outward_mobility, outward_left, outward_right;
        face_mobility(config, input, current, cell,
          &inward_mobility, &inward_left, &inward_right);
        face_mobility(config, input, current, cell + 1,
          &outward_mobility, &outward_left, &outward_right);
        const double inward_factor = input->mesh->face_measure[cell]*
          input->decontaminated_pressure_gradient_face[cell];
        const double outward_factor = input->mesh->face_measure[cell + 1]*
          input->decontaminated_pressure_gradient_face[cell + 1];
        lower[cell] = inward_factor*inward_left;
        upper[cell] = -outward_factor*outward_right;
        diagonal[cell] = measure*(1./input->dt +
          input->mean_surface_velocity_divergence[cell]) +
          inward_factor*inward_right - outward_factor*outward_left;
      }
      rhs[cell] = -rhs[cell];
    }
    for (size_t cell = 1; cell < cells; cell++) {
      if (!isfinite(diagonal[cell - 1]) ||
          fabs(diagonal[cell - 1]) <= DBL_MIN) {
        failure_kind = MUSEHANE_CORE_FAILURE_SINGULAR_NEWTON_PIVOT;
        goto cleanup;
      }
      const double factor = lower[cell]/diagonal[cell - 1];
      diagonal[cell] -= factor*upper[cell - 1];
      rhs[cell] -= factor*rhs[cell - 1];
    }
    if (!isfinite(diagonal[cells - 1]) ||
        fabs(diagonal[cells - 1]) <= DBL_MIN) {
      failure_kind = MUSEHANE_CORE_FAILURE_SINGULAR_NEWTON_PIVOT;
      goto cleanup;
    }
    delta[cells - 1] = rhs[cells - 1]/diagonal[cells - 1];
    for (size_t cell = cells - 1; cell > 0; cell--) {
      const size_t previous = cell - 1;
      if (!isfinite(diagonal[previous]) ||
          fabs(diagonal[previous]) <= DBL_MIN) {
        failure_kind = MUSEHANE_CORE_FAILURE_SINGULAR_NEWTON_PIVOT;
        goto cleanup;
      }
      delta[previous] = (rhs[previous] -
        upper[previous]*delta[cell])/diagonal[previous];
    }
    bool accepted = false;
    failure_backtracks = 0;
    for (unsigned backtrack = 0; backtrack < 60; backtrack++) {
      failure_backtracks = backtrack + 1;
      const double step = ldexp(1., -(int)backtrack);
      bool positive = true;
      for (size_t cell = 0; cell < cells; cell++) {
        next[cell] = input->active[cell] ?
          current[cell] + step*delta[cell] :
          input->estimated_thickness[cell];
        positive &= isfinite(next[cell]) && next[cell] > 0.;
      }
      if (!positive)
        continue;
      const double trial_relative_merit = equation_residual(config, input,
        next, trial_flux, trial_pressure_rate, pressure_lower,
        pressure_diagonal, pressure_upper, trial_residual);
      const double trial_line_search_merit = fixed_scale_residual(
        input, trial_residual, line_search_scale);
      if (isfinite(trial_relative_merit) &&
          isfinite(trial_line_search_merit) &&
          (trial_relative_merit <= config->nonlinear_relative_tolerance ||
           trial_line_search_merit < merit*(1. - 1.e-4*step))) {
        accepted = true;
        break;
      }
    }
    if (!accepted) {
      failure_kind = MUSEHANE_CORE_FAILURE_LINE_SEARCH_EXHAUSTED;
      goto cleanup;
    }
    failure_backtracks = 0;
    memcpy(current, next, cells*sizeof *current);
  }
  if (!converged) {
    failure_kind = MUSEHANE_CORE_FAILURE_ITERATION_LIMIT;
    goto cleanup;
  }

  const double final_merit = equation_residual(config, input, current,
    flux, pressure_rate, pressure_lower, pressure_diagonal,
    pressure_upper, NULL);
  if (!isfinite(final_merit) ||
      final_merit > config->nonlinear_relative_tolerance) {
    failure_kind = MUSEHANE_CORE_FAILURE_FINAL_RESIDUAL_REJECTED;
    goto cleanup;
  }
  const MusehaneInventory ledger =
    step_inventory(input, current, pressure_rate);
  const double inventory_before = ledger.before;
  const double inventory_after = ledger.after;
  const double pressure_inventory_change = ledger.pressure_change;
  const double kinematic_inventory_change = ledger.kinematic_change;
  double minimum_thickness = INFINITY;
  for (size_t cell = 0; cell < cells; cell++) {
    if (!input->active[cell])
      continue;
    minimum_thickness = fmin(minimum_thickness, current[cell]);
  }
  const double expected_inventory_after = ledger.expected_after;
  const double closure_error = ledger.closure_error;
  const double inventory_scale = ledger.scale;
  if (!isfinite(inventory_after) || !isfinite(expected_inventory_after) ||
      !isfinite(closure_error) || fabs(closure_error) >
        config->inventory_relative_tolerance*inventory_scale) {
    if (failure)
      *failure = (MusehaneCoreFailure) {
        .inventory_diagnostic_valid = true,
        .inventory_before = inventory_before,
        .inventory_after = inventory_after,
        .expected_inventory_after = expected_inventory_after,
        .inventory_closure_error = closure_error,
        .inventory_relative_error = fabs(closure_error)/inventory_scale,
        .inventory_relative_tolerance =
          config->inventory_relative_tolerance
      };
    status = MUSEHANE_INVENTORY_MISMATCH;
    goto cleanup;
  }

  double poiseuille_dissipation_rate = 0.;
  if (config->discretization ==
      MUSEHANE_DISCRETIZATION_PUBLISHED_EQUATION23) {
    const double mobility_factor = 1./(12.*config->water_viscosity);
    for (size_t cell = 0; cell < cells; cell++) {
      if (!input->active[cell])
        continue;
      const double gradient =
        input->decontaminated_pressure_gradient_cell[cell];
      poiseuille_dissipation_rate += input->mesh->cell_measure[cell]*
        mobility_factor*current[cell]*current[cell]*current[cell]*
        gradient*gradient;
    }
  }
  else {
    for (size_t face = 0; face <= cells; face++) {
      double mobility, derivative_left, derivative_right;
      face_mobility(config, input, current, face,
        &mobility, &derivative_left, &derivative_right);
      const double dual_width = face == 0 ?
        .5*input->mesh->cell_width[0] :
        face == cells ? .5*input->mesh->cell_width[cells - 1] :
        .5*(input->mesh->cell_width[face - 1] +
            input->mesh->cell_width[face]);
      const double gradient =
        input->decontaminated_pressure_gradient_face[face];
      poiseuille_dissipation_rate += input->mesh->face_measure[face]*
        dual_width*mobility*gradient*gradient;
    }
  }
  if (!isfinite(poiseuille_dissipation_rate) ||
      poiseuille_dissipation_rate < 0.) {
    status = MUSEHANE_INVALID_INPUT;
    goto cleanup;
  }

  double *const output_thickness = output->thickness;
  double *const output_flux = output->pressure_flux_integral_face;
  memcpy(output_thickness, current, cells*sizeof *current);
  memcpy(output_flux, flux, (cells + 1)*sizeof *flux);
  *output = (MusehaneStepOutput) {
    .thickness = output_thickness,
    .pressure_flux_integral_face = output_flux,
    .converged = true,
    .iterations = iteration,
    .minimum_thickness = minimum_thickness,
    .inventory_before = inventory_before,
    .inventory_after = inventory_after,
    .expected_inventory_after = expected_inventory_after,
    .pressure_inventory_change = pressure_inventory_change,
    .kinematic_inventory_change = kinematic_inventory_change,
    .inventory_closure_error = closure_error,
    .poiseuille_dissipation_rate = poiseuille_dissipation_rate
  };
  status = MUSEHANE_OK;

cleanup:
  if (status == MUSEHANE_NONLINEAR_NOT_CONVERGED && failure &&
      !failure->valid) {
    const double failure_merit = equation_residual(config, input, current,
      flux, pressure_rate, pressure_lower, pressure_diagonal,
      pressure_upper, NULL);
    double worst = -1.;
    for (size_t cell = 0; isfinite(failure_merit) && cell < cells; cell++) {
      if (!input->active[cell])
        continue;
      const double measure = input->mesh->cell_measure[cell];
      const double denominator = measure*(1./input->dt +
        input->mean_surface_velocity_divergence[cell]);
      const double old_term = measure*input->old_thickness[cell]/input->dt;
      const double balance = denominator*current[cell] - old_term -
        pressure_rate[cell];
      const double scale = fmax(DBL_MIN,
        fabs(denominator*current[cell]) + fabs(old_term) +
        fabs(pressure_rate[cell]));
      const double relative = fabs(balance)/scale;
      if (isfinite(relative) && relative > worst) {
        worst = relative;
        *failure = (MusehaneCoreFailure) {
          .valid = true,
          .kind = failure_kind,
          .cell_index = cell,
          .iteration = iteration > config->maximum_iterations ?
            config->maximum_iterations : iteration,
          .line_search_backtracks = failure_backtracks,
          .initial_merit = initial_merit,
          .final_merit = failure_merit,
          .numerator = old_term + pressure_rate[cell],
          .denominator = denominator,
          .old_thickness = input->old_thickness[cell],
          .current_thickness = current[cell],
          .flux_in = flux[cell],
          .flux_out = flux[cell + 1],
          .surface_velocity_divergence =
            input->mean_surface_velocity_divergence[cell],
          .equation_relative_residual = relative
        };
      }
    }
  }
  free(current);
  free(next);
  free(flux);
  free(trial_flux);
  free(pressure_rate);
  free(trial_pressure_rate);
  free(pressure_lower);
  free(pressure_diagonal);
  free(pressure_upper);
  free(lower);
  free(diagonal);
  free(upper);
  free(rhs);
  free(delta);
  free(trial_residual);
  free(line_search_scale);
  return status;
}

MusehaneStatus musehane_surface_step (
  const MusehaneConfig *config,
  const MusehaneStepInput *input,
  MusehaneStepOutput *output)
{
  return musehane_surface_step_diagnosed(config, input, output, NULL);
}
