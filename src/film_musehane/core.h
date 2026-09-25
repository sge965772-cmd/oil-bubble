#ifndef FILM_MUSEHANE_CORE_H
#define FILM_MUSEHANE_CORE_H

#include <stdbool.h>
#include <stddef.h>

typedef enum {
  MUSEHANE_OK = 0,
  MUSEHANE_INVALID_CONFIG = 1,
  MUSEHANE_INVALID_INPUT = 2,
  MUSEHANE_NONPOSITIVE_THICKNESS = 3,
  MUSEHANE_NONLINEAR_NOT_CONVERGED = 4,
  MUSEHANE_INVENTORY_MISMATCH = 5,
  MUSEHANE_ALLOCATION_FAILED = 6,
  MUSEHANE_NO_ACTIVE_REGION = 7
} MusehaneStatus;

typedef enum {
  MUSEHANE_DISCRETIZATION_CONSERVATIVE_HARMONIC = 0,
  MUSEHANE_DISCRETIZATION_PUBLISHED_EQUATION23 = 1
} MusehaneDiscretization;

typedef struct {
  /* Retained for dimensional/configuration audit. It cancels from equation
     (11) because the carrier film is treated as constant-density water. */
  double water_density;
  double water_viscosity;
  MusehaneDiscretization discretization;
  double nonlinear_relative_tolerance;
  double inventory_relative_tolerance;
  unsigned maximum_iterations;
} MusehaneConfig;

/* Measures omit the common 2*pi factor. For an axisymmetric radial mesh,
   cell_measure=(r_outer^2-r_inner^2)/2 and face_measure=r_face. */
typedef struct {
  size_t cell_count;
  const double *cell_radius;
  const double *cell_width;
  const double *cell_measure;
  const double *face_measure;
} MusehaneSurfaceMesh;

typedef struct {
  double time;
  double dt;
  const MusehaneSurfaceMesh *mesh;
  const unsigned char *active;
  /* Current geometric separation estimate h0. Inactive cells are assigned
     this value, so newly activated cells inherit the current geometry. */
  const double *estimated_thickness;
  const double *old_thickness;
  const double *mean_surface_velocity_divergence;
  const double *decontaminated_pressure_gradient_cell;
  const double *decontaminated_pressure_gradient_face;
} MusehaneStepInput;

typedef struct {
  double *thickness;
  /* A_f M grad(p') on each face, where M=h^3/(12 mu). The physical
     Poiseuille volume-flow integral has the opposite sign. */
  double *pressure_flux_integral_face;
  bool converged;
  unsigned iterations;
  double minimum_thickness;
  double inventory_before;
  double inventory_after;
  double expected_inventory_after;
  double pressure_inventory_change;
  double kinematic_inventory_change;
  double inventory_closure_error;
  /* Non-negative integral of M |grad(p')|^2 over the surface. The common
     axisymmetric 2*pi factor is omitted consistently with the mesh measures. */
  double poiseuille_dissipation_rate;
} MusehaneStepOutput;

typedef enum {
  MUSEHANE_CORE_FAILURE_NONE = 0,
  MUSEHANE_CORE_FAILURE_NONPOSITIVE_DENOMINATOR = 1,
  MUSEHANE_CORE_FAILURE_INITIAL_RESIDUAL_NONFINITE = 2,
  MUSEHANE_CORE_FAILURE_ITERATION_RESIDUAL_NONFINITE = 3,
  MUSEHANE_CORE_FAILURE_SINGULAR_NEWTON_PIVOT = 4,
  MUSEHANE_CORE_FAILURE_LINE_SEARCH_EXHAUSTED = 5,
  MUSEHANE_CORE_FAILURE_ITERATION_LIMIT = 6,
  MUSEHANE_CORE_FAILURE_FINAL_RESIDUAL_REJECTED = 7
} MusehaneCoreFailureKind;

typedef struct {
  bool valid;
  MusehaneCoreFailureKind kind;
  size_t cell_index;
  unsigned iteration;
  unsigned line_search_backtracks;
  double initial_merit;
  double final_merit;
  double numerator;
  double denominator;
  double old_thickness;
  double current_thickness;
  double flux_in;
  double flux_out;
  double surface_velocity_divergence;
  double equation_relative_residual;
  bool inventory_diagnostic_valid;
  double inventory_before;
  double inventory_after;
  double expected_inventory_after;
  double inventory_closure_error;
  double inventory_relative_error;
  double inventory_relative_tolerance;
} MusehaneCoreFailure;

/* Backward-Euler/Newton finite-volume update of Musehane et al. (2018),
   equation (11), with p' supplied by the DNS pressure-decontamination seam.
   Output buffers and metadata are committed only when the full step passes. */
MusehaneStatus musehane_surface_step (
  const MusehaneConfig *config,
  const MusehaneStepInput *input,
  MusehaneStepOutput *output);

/* On a rejected thickness update, reports the first failing ring without
   committing the output buffers. NULL disables this optional audit. */
MusehaneStatus musehane_surface_step_diagnosed (
  const MusehaneConfig *config,
  const MusehaneStepInput *input,
  MusehaneStepOutput *output,
  MusehaneCoreFailure *failure);

#endif
