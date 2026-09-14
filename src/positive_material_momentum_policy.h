#ifndef POSITIVE_MATERIAL_MOMENTUM_POLICY_H
#define POSITIVE_MATERIAL_MOMENTUM_POLICY_H

#ifndef ENABLE_POSITIVE_MATERIAL_MOMENTUM
# define ENABLE_POSITIVE_MATERIAL_MOMENTUM 1
#endif
#ifndef ENABLE_CENTERED_STOKES_CONTROL
# define ENABLE_CENTERED_STOKES_CONTROL 0
#endif
#ifndef MOMENTUM_TRANSPORT_VERSION
# if ENABLE_POSITIVE_MATERIAL_MOMENTUM
#  define MOMENTUM_TRANSPORT_VERSION 2
# else
#  define MOMENTUM_TRANSPORT_VERSION 0
# endif
#endif

#if MOMENTUM_TRANSPORT_VERSION == 4
# define AGGREGATE_MATERIAL_SIMPLEX_FLUX_LIMITING 1
# include "aggregate_material_momentum_transport.h"
# define MOMENTUM_TRANSPORT_MODE 4
#elif MOMENTUM_TRANSPORT_VERSION == 3
# include "aggregate_material_momentum_transport.h"
# define MOMENTUM_TRANSPORT_MODE 3
#elif MOMENTUM_TRANSPORT_VERSION == 2
# include "positive_material_momentum_transport.h"
# define MOMENTUM_TRANSPORT_MODE 2
#elif MOMENTUM_TRANSPORT_VERSION == 0
# define MOMENTUM_TRANSPORT_MODE 0
static const long dual_momentum_invalid_cell_count = 0;
static const long dual_momentum_invalid_face_count = 0;
static const int dual_momentum_transport_failed = 0;
static const double dual_momentum_max_cell_partition_defect = 0.;
static const double dual_momentum_max_face_partition_defect = 0.;
static const double dual_momentum_max_vof_bound_defect = 0.;
static const double dual_momentum_max_partition_correction = 0.;
static const double dual_momentum_last_transport_cpu_seconds = 0.;
static const long dual_momentum_limited_face_count = 0;
static const long dual_momentum_limited_face_count_max = 0;
static const double dual_momentum_min_flux_alpha = 1.;
static const double dual_momentum_min_flux_alpha_ever = 1.;
static const double dual_momentum_max_simplex_closure_defect = 0.;
static const double dual_momentum_max_flux_fraction_correction = 0.;
static const double dual_momentum_max_flux_courant_correction = 0.;
static const long dual_momentum_nested_plic_fallback_count = 0;
static const long dual_momentum_nested_plic_fallback_count_max = 0;

# if ENABLE_CENTERED_STOKES_CONTROL
event defaults (i = 0)
{
  stokes = true;
}

event stability (i++)
  dtmax = timestep(uf, dtmax);
# endif
#else
# error "MOMENTUM_TRANSPORT_VERSION must be 0, 2, 3 or 4"
#endif

#endif
