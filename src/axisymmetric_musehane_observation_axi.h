#ifndef AXISYMMETRIC_MUSEHANE_OBSERVATION_AXI_H
#define AXISYMMETRIC_MUSEHANE_OBSERVATION_AXI_H

#include <math.h>

#include "axisymmetric_film_distribution.h"
#include "film_model/contract.h"
#include "film_musehane/dns_observation.h"
#include "track_a_subgrid_observation_policy.h"

typedef enum {
  AXISYMMETRIC_MUSEHANE_OBSERVATION_OK = 0,
  AXISYMMETRIC_MUSEHANE_OBSERVATION_INVALID_GEOMETRY = 1,
  AXISYMMETRIC_MUSEHANE_OBSERVATION_INVALID_FIELD = 2,
  AXISYMMETRIC_MUSEHANE_OBSERVATION_ASSEMBLY_FAILED = 3
} AxisymmetricMusehaneObservationStatus;

typedef struct {
  AxisymmetricMusehaneObservationStatus status;
  int sampled_bins;
  int source_valid_bins;
  int interpolated_bins;
  int subgrid_support_bins;
  int native_ordering_lost_bins;
  double source_valid_area_fraction;
  MusehaneDNSObservationAudit assembly_audit;
  unsigned char native_valid[AXISYMMETRIC_MATCHED_FILM_BINS];
  FilmRingObservation rings[AXISYMMETRIC_MATCHED_FILM_BINS];
  FilmDNSView view;
} AxisymmetricMusehaneObservation;

/* Read-only Basilisk adapter. The caller constructs pressure_gradient and
   capillary_force_density from the native DNS operators at one time level. */
static inline AxisymmetricMusehaneObservationStatus
sample_axisymmetric_musehane_observation (
  double support_half_width,
  scalar water_pressure,
  vector pressure_gradient,
  vector capillary_force_density,
  double sample_time,
  double sample_dt,
  unsigned long sample_iteration,
  bool anchor_released,
  uint32_t lower_identity,
  uint32_t upper_identity,
  const TrackASubgridObservationState *subgrid_state,
  AxisymmetricMusehaneObservation *output)
{
  AxisymmetricMusehaneObservation diagnostic = {
    .status = AXISYMMETRIC_MUSEHANE_OBSERVATION_INVALID_GEOMETRY,
    .source_valid_area_fraction = NAN
  };
  if (output)
    *output = diagnostic;
  if (!output || !isfinite(support_half_width) || support_half_width <= 0. ||
      !isfinite(sample_time) || sample_time < 0. ||
      !isfinite(sample_dt) || sample_dt <= 0. ||
      lower_identity == 0U || upper_identity == 0U ||
      lower_identity == upper_identity)
    return AXISYMMETRIC_MUSEHANE_OBSERVATION_INVALID_GEOMETRY;

  AxisymmetricFilmSamples geometry = sample_axisymmetric_water_film(
    support_half_width);
  diagnostic.sampled_bins = geometry.bins;
  diagnostic.source_valid_bins = geometry.valid_bins;
  diagnostic.interpolated_bins = geometry.interpolated_bins;
  diagnostic.source_valid_area_fraction = geometry.valid_area_fraction;
  *output = diagnostic;
  if (geometry.bins != AXISYMMETRIC_MATCHED_FILM_BINS)
    return AXISYMMETRIC_MUSEHANE_OBSERVATION_INVALID_GEOMETRY;
  const int sampled_bins = geometry.bins;

  const bool subgrid_state_ready = subgrid_state && subgrid_state->valid &&
    subgrid_state->ring_count == sampled_bins && subgrid_state->active_n &&
    subgrid_state->active_np1 && subgrid_state->thickness;
  TrackASubgridObservationGeometry selected[AXISYMMETRIC_MATCHED_FILM_BINS];

  coord points[AXISYMMETRIC_MATCHED_FILM_BINS];
  coord interface_points[2*AXISYMMETRIC_MATCHED_FILM_BINS];
  int subgrid_support_bins = 0, native_ordering_lost_bins = 0;
  for (int bin = 0; bin < sampled_bins; bin++) {
    const AxisymmetricFilmSample *sample = &geometry.bin[bin];
    const bool subgrid_authoritative = subgrid_state_ready &&
      track_a_subgrid_observation_is_authoritative(
        subgrid_state->active_n[bin] != 0U,
        subgrid_state->active_np1[bin] != 0U,
        track_a_observation_release_is_latched(
          subgrid_state->release_latch_mask, (size_t)bin),
        sample->valid);
    const double model_thickness = subgrid_authoritative ?
      subgrid_state->thickness[bin] : NAN;
    selected[bin] = track_a_select_observation_geometry(
      sample->lower_sample_count, sample->upper_sample_count,
      sample->lower_point.x, sample->upper_point.x, sample->valid,
      subgrid_authoritative, model_thickness);
    if (selected[bin].status == TRACK_A_OBSERVATION_INVALID ||
        !isfinite(sample->delta) || !(sample->delta > 0.)) {
      output->status = AXISYMMETRIC_MUSEHANE_OBSERVATION_INVALID_GEOMETRY;
      return AXISYMMETRIC_MUSEHANE_OBSERVATION_INVALID_GEOMETRY;
    }
    if (selected[bin].status == TRACK_A_OBSERVATION_SUBGRID)
      subgrid_support_bins++;
    if (selected[bin].native_ordering_lost)
      native_ordering_lost_bins++;
    points[bin] = (coord) {
      .x = selected[bin].midpoint,
      .y = geometry.bin[bin].radius
    };
    interface_points[2*bin] = (coord) {
      .x = selected[bin].lower_position, .y = geometry.bin[bin].radius};
    interface_points[2*bin + 1] = (coord) {
      .x = selected[bin].upper_position, .y = geometry.bin[bin].radius};
  }
  double values[5*AXISYMMETRIC_MATCHED_FILM_BINS];
  interpolate_array((scalar *) {
      water_pressure,
      pressure_gradient.x, pressure_gradient.y,
      capillary_force_density.x, capillary_force_density.y
    }, points, sampled_bins, values, true);
  double interface_velocity[4*AXISYMMETRIC_MATCHED_FILM_BINS];
  interpolate_array((scalar *) {u.x, u.y}, interface_points,
    2*sampled_bins, interface_velocity, true);

  MusehaneDNSRingSample samples[AXISYMMETRIC_MATCHED_FILM_BINS] = {{0}};
  for (int bin = 0; bin < sampled_bins; bin++) {
    for (int field = 0; field < 5; field++)
      if (values[5*bin + field] == nodata ||
          !isfinite(values[5*bin + field])) {
        output->status = AXISYMMETRIC_MUSEHANE_OBSERVATION_INVALID_FIELD;
        return AXISYMMETRIC_MUSEHANE_OBSERVATION_INVALID_FIELD;
      }
    const AxisymmetricFilmSample sample = geometry.bin[bin];
    for (int component = 0; component < 4; component++)
      if (interface_velocity[4*bin + component] == nodata ||
          !isfinite(interface_velocity[4*bin + component])) {
        output->status = AXISYMMETRIC_MUSEHANE_OBSERVATION_INVALID_FIELD;
        return AXISYMMETRIC_MUSEHANE_OBSERVATION_INVALID_FIELD;
      }
    const double inner_radius = sample.radius - .5*sample.width;
    const double outer_radius = sample.radius + .5*sample.width;
    samples[bin] = (MusehaneDNSRingSample) {
      .radius = sample.radius,
      .width = sample.width,
      .area_weight =
        3.141592653589793238462643383279502884*
        (sq(outer_radius) - sq(inner_radius)),
      .lower_position = selected[bin].lower_position,
      .upper_position = selected[bin].upper_position,
      .local_delta = sample.delta,
      .lower_velocity = {
        .x = interface_velocity[4*bin],
        .r = interface_velocity[4*bin + 1]},
      .upper_velocity = {
        .x = interface_velocity[4*bin + 2],
        .r = interface_velocity[4*bin + 3]},
      .water_pressure = values[5*bin],
      .pressure_gradient = {
        .x = values[5*bin + 1], .r = values[5*bin + 2]},
      .capillary_force_density = {
        .x = values[5*bin + 3], .r = values[5*bin + 4]},
      .lower_to_upper_normal = {.x = 1., .r = 0.},
      .geometry_valid = true,
      .lower_identity_valid = true,
      .upper_identity_valid = true
    };
  }

  AxisymmetricMusehaneObservation candidate = {
    .status = AXISYMMETRIC_MUSEHANE_OBSERVATION_OK,
    .sampled_bins = sampled_bins,
    .source_valid_bins = geometry.valid_bins,
    .interpolated_bins = geometry.interpolated_bins,
    .subgrid_support_bins = subgrid_support_bins,
    .native_ordering_lost_bins = native_ordering_lost_bins,
    .source_valid_area_fraction = geometry.valid_area_fraction
  };
  for (int bin = 0; bin < sampled_bins; bin++)
    candidate.native_valid[bin] = geometry.bin[bin].valid ? 1U : 0U;
  MusehaneDNSRingObservation assembled[AXISYMMETRIC_MATCHED_FILM_BINS] =
    {{0}};
  if (musehane_assemble_dns_observations(
        (size_t)sampled_bins, samples, assembled,
        &candidate.assembly_audit) != MUSEHANE_DNS_OBSERVATION_OK) {
    candidate.status = AXISYMMETRIC_MUSEHANE_OBSERVATION_ASSEMBLY_FAILED;
    *output = candidate;
    return AXISYMMETRIC_MUSEHANE_OBSERVATION_ASSEMBLY_FAILED;
  }
  for (int bin = 0; bin < sampled_bins; bin++) {
    const MusehaneDNSRingObservation source = assembled[bin];
    candidate.rings[bin] = (FilmRingObservation) {
      .radius = source.radius,
      .width = source.width,
      .area_weight = source.area_weight,
      .lower_position = source.lower_position,
      .upper_position = source.upper_position,
      .gap = source.gap,
      .midpoint = source.midpoint,
      .local_delta = source.local_delta,
      .lower_normal_velocity = source.lower_normal_velocity,
      .upper_normal_velocity = source.upper_normal_velocity,
      .lower_tangential_velocity = source.lower_tangential_velocity,
      .upper_tangential_velocity = source.upper_tangential_velocity,
      .water_pressure = source.water_pressure,
      .pressure_gradient_tangent = source.pressure_gradient_tangent,
      .capillary_force_density_tangent =
        source.capillary_force_density_tangent,
      .lower_normal_x = source.lower_normal_x,
      .lower_normal_r = source.lower_normal_r,
      .upper_normal_x = source.upper_normal_x,
      .upper_normal_r = source.upper_normal_r,
      .geometry_valid = source.geometry_valid,
      .native_geometry_valid = geometry.bin[bin].valid,
      .resolved_release_latched =
        track_a_observation_release_is_latched(
          subgrid_state_ready ? subgrid_state->release_latch_mask : 0U,
          (size_t)bin),
      .lower_identity_valid = source.lower_identity_valid,
      .upper_identity_valid = source.upper_identity_valid
    };
  }
  candidate.view = (FilmDNSView) {
    .time = sample_time,
    .dt = sample_dt,
    .iteration = sample_iteration,
    .ring_count = (size_t)sampled_bins,
    .rings = NULL,
    .anchor_released = anchor_released,
    .lower_identity = lower_identity,
    .upper_identity = upper_identity
  };
  *output = candidate;
  output->view.rings = output->rings;
  return AXISYMMETRIC_MUSEHANE_OBSERVATION_OK;
}

#endif
