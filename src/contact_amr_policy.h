#ifndef CONTACT_AMR_POLICY_H
#define CONTACT_AMR_POLICY_H

#include <math.h>
#include <stdbool.h>

typedef struct {
  double finest_delta;
  double activation_cells;
  double axial_padding_cells;
  double radial_support;
  double radial_padding_cells;
} ContactAmrPolicy;

typedef struct {
  bool valid;
  double lower_top;
  double upper_bottom;
  bool virtual_bounds_present;
  bool active_state_support;
  double minimum_gap;
  double axial_minimum;
  double axial_maximum;
} ContactAmrGeometry;

static inline bool contact_amr_geometry_is_valid (
  const ContactAmrGeometry * geometry)
{
  if (!geometry || !geometry->valid)
    return false;
  if (geometry->virtual_bounds_present)
    return isfinite(geometry->minimum_gap) && geometry->minimum_gap > 0. &&
      isfinite(geometry->axial_minimum) &&
      isfinite(geometry->axial_maximum) &&
      geometry->axial_maximum >= geometry->axial_minimum;
  return isfinite(geometry->lower_top) &&
    isfinite(geometry->upper_bottom);
}

static inline double contact_amr_geometry_minimum_gap (
  const ContactAmrGeometry * geometry)
{
  if (!contact_amr_geometry_is_valid(geometry))
    return NAN;
  return geometry->virtual_bounds_present ? geometry->minimum_gap :
    geometry->upper_bottom - geometry->lower_top;
}

static inline bool contact_amr_policy_is_valid (
  const ContactAmrPolicy * policy)
{
  return policy && isfinite(policy->finest_delta) &&
    policy->finest_delta > 0. && isfinite(policy->activation_cells) &&
    policy->activation_cells > 0. &&
    isfinite(policy->axial_padding_cells) &&
    policy->axial_padding_cells >= 0. &&
    isfinite(policy->radial_support) && policy->radial_support > 0. &&
    isfinite(policy->radial_padding_cells) &&
    policy->radial_padding_cells >= 0.;
}

static inline bool contact_amr_is_active (
  const ContactAmrPolicy * policy, const ContactAmrGeometry * geometry)
{
  if (!contact_amr_policy_is_valid(policy) ||
      !contact_amr_geometry_is_valid(geometry))
    return false;
  if (geometry->virtual_bounds_present && geometry->active_state_support)
    return true;
  const double gap = contact_amr_geometry_minimum_gap(geometry);
  return gap <= policy->activation_cells*policy->finest_delta;
}

static inline bool contact_amr_contains (
  const ContactAmrPolicy * policy, const ContactAmrGeometry * geometry,
  double axial_coordinate, double radial_coordinate)
{
  if (!contact_amr_is_active(policy, geometry) ||
      !isfinite(axial_coordinate) || !isfinite(radial_coordinate) ||
      radial_coordinate < 0.)
    return false;
  double axial_padding = policy->axial_padding_cells*policy->finest_delta;
  double radial_limit = policy->radial_support +
    policy->radial_padding_cells*policy->finest_delta;
  double axial_minimum = (geometry->virtual_bounds_present ?
    geometry->axial_minimum :
    fmin(geometry->lower_top, geometry->upper_bottom)) - axial_padding;
  double axial_maximum = (geometry->virtual_bounds_present ?
    geometry->axial_maximum :
    fmax(geometry->lower_top, geometry->upper_bottom)) + axial_padding;
  return axial_coordinate >= axial_minimum &&
    axial_coordinate <= axial_maximum &&
    radial_coordinate <= radial_limit;
}

static inline bool contact_amr_cell_intersects (
  const ContactAmrPolicy * policy, const ContactAmrGeometry * geometry,
  double axial_center, double radial_center, double cell_delta)
{
  if (!contact_amr_is_active(policy, geometry) ||
      !isfinite(axial_center) || !isfinite(radial_center) ||
      !isfinite(cell_delta) || cell_delta <= 0.)
    return false;
  const double half_delta = .5*cell_delta;
  const double axial_padding =
    policy->axial_padding_cells*policy->finest_delta;
  const double radial_limit = policy->radial_support +
    policy->radial_padding_cells*policy->finest_delta;
  const double axial_minimum = (geometry->virtual_bounds_present ?
    geometry->axial_minimum :
    fmin(geometry->lower_top, geometry->upper_bottom)) - axial_padding;
  const double axial_maximum = (geometry->virtual_bounds_present ?
    geometry->axial_maximum :
    fmax(geometry->lower_top, geometry->upper_bottom)) + axial_padding;
  const double cell_axial_minimum = axial_center - half_delta;
  const double cell_axial_maximum = axial_center + half_delta;
  const double cell_radial_minimum = fmax(0., radial_center - half_delta);
  return cell_axial_maximum >= axial_minimum &&
    cell_axial_minimum <= axial_maximum &&
    cell_radial_minimum <= radial_limit;
}

#endif
