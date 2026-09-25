#ifndef AXISYMMETRIC_FILM_DISTRIBUTION_H
#define AXISYMMETRIC_FILM_DISTRIBUTION_H

#include <float.h>
#include "film_distribution_statistics.h"
#include "axisymmetric_film_samples.h"

#ifndef AXISYMMETRIC_FILM_RADIAL_BINS
# define AXISYMMETRIC_FILM_RADIAL_BINS 32
#endif

typedef struct {
  double support_half_width;
  int nonpositive_gap_bins;
  FilmDistributionStatistics statistics;
} AxisymmetricFilmDistribution;

static inline bool film_segment_crossing_at_radius (
  coord first,
  coord second,
  double radius,
  double * minimum_x,
  double * maximum_x)
{
  AxisymmetricFilmSegment segment = {
    .first = {.x = first.x, .r = first.y},
    .second = {.x = second.x, .r = second.y}
  };
  return axisymmetric_film_segment_crossing_at_radius(
    segment, radius, minimum_x, maximum_x);
}

static inline AxisymmetricFilmSamples
sample_axisymmetric_water_film (
  double support_half_width)
{
  enum { bins = AXISYMMETRIC_FILM_RADIAL_BINS };
  double lower[bins], upper[bins];
  double lower_delta[bins], upper_delta[bins];
  int lower_samples[bins], upper_samples[bins];
  int invalid_plic_segments = 0;
  int skipped_degenerate_plic_facets = 0;
  for (int bin = 0; bin < bins; bin++) {
    lower[bin] = -HUGE_VAL;
    upper[bin] = HUGE_VAL;
    lower_delta[bin] = upper_delta[bin] = 0.;
    lower_samples[bin] = upper_samples[bin] = 0;
  }

  face vector no_face_fractions = {{-1}};
  foreach(reduction(max:lower[:bins]) reduction(min:upper[:bins])
          reduction(max:lower_delta[:bins])
          reduction(max:upper_delta[:bins])
          reduction(+:lower_samples[:bins])
          reduction(+:upper_samples[:bins])
          reduction(+:invalid_plic_segments)
          reduction(+:skipped_degenerate_plic_facets)) {
    coord lower_segment[2], upper_segment[2];
    int lower_facets = 0, upper_facets = 0;
    if (lower_envelope_fraction[] > 1.e-6 &&
        lower_envelope_fraction[] < 1. - 1.e-6) {
      coord normal = facet_normal(
        point, lower_envelope_fraction, no_face_fractions);
      if (!isfinite(normal.x) || !isfinite(normal.y) ||
          fabs(normal.x) + fabs(normal.y) <= DBL_MIN)
        invalid_plic_segments++;
      else {
        double alpha = plane_alpha(lower_envelope_fraction[], normal);
        if (!isfinite(alpha))
          invalid_plic_segments++;
        else {
          lower_facets = facets(normal, alpha, lower_segment);
          if (lower_facets == 2) {
            for (int endpoint = 0; endpoint < 2; endpoint++) {
              lower_segment[endpoint].x = x + lower_segment[endpoint].x*Delta;
              lower_segment[endpoint].y = y + lower_segment[endpoint].y*Delta;
            }
            AxisymmetricFilmSegment segment = {
              .first = {.x = lower_segment[0].x, .r = lower_segment[0].y},
              .second = {.x = lower_segment[1].x, .r = lower_segment[1].y}
            };
            if (!axisymmetric_film_segment_is_finite(segment)) {
              lower_facets = 0;
              invalid_plic_segments++;
            }
          }
          else if (axisymmetric_film_facet_count_is_degenerate(lower_facets)) {
            lower_facets = 0;
            skipped_degenerate_plic_facets++;
          }
        }
      }
    }
    if (upper_envelope_fraction[] > 1.e-6 &&
        upper_envelope_fraction[] < 1. - 1.e-6) {
      coord normal = facet_normal(
        point, upper_envelope_fraction, no_face_fractions);
      if (!isfinite(normal.x) || !isfinite(normal.y) ||
          fabs(normal.x) + fabs(normal.y) <= DBL_MIN)
        invalid_plic_segments++;
      else {
        double alpha = plane_alpha(upper_envelope_fraction[], normal);
        if (!isfinite(alpha))
          invalid_plic_segments++;
        else {
          upper_facets = facets(normal, alpha, upper_segment);
          if (upper_facets == 2) {
            for (int endpoint = 0; endpoint < 2; endpoint++) {
              upper_segment[endpoint].x = x + upper_segment[endpoint].x*Delta;
              upper_segment[endpoint].y = y + upper_segment[endpoint].y*Delta;
            }
            AxisymmetricFilmSegment segment = {
              .first = {.x = upper_segment[0].x, .r = upper_segment[0].y},
              .second = {.x = upper_segment[1].x, .r = upper_segment[1].y}
            };
            if (!axisymmetric_film_segment_is_finite(segment)) {
              upper_facets = 0;
              invalid_plic_segments++;
            }
          }
          else if (axisymmetric_film_facet_count_is_degenerate(upper_facets)) {
            upper_facets = 0;
            skipped_degenerate_plic_facets++;
          }
        }
      }
    }

    if (lower_facets == 2 || upper_facets == 2)
      for (int bin = 0; bin < bins; bin++) {
        double radius = (bin + 0.5)*support_half_width/bins;
        double minimum_x, maximum_x;
        if (lower_facets == 2 && film_segment_crossing_at_radius(
              lower_segment[0], lower_segment[1], radius,
              &minimum_x, &maximum_x)) {
          lower[bin] = max(lower[bin], maximum_x);
          lower_delta[bin] = max(lower_delta[bin], Delta);
          lower_samples[bin]++;
        }
        if (upper_facets == 2 && film_segment_crossing_at_radius(
              upper_segment[0], upper_segment[1], radius,
              &minimum_x, &maximum_x)) {
          upper[bin] = min(upper[bin], minimum_x);
          upper_delta[bin] = max(upper_delta[bin], Delta);
          upper_samples[bin]++;
        }
      }
  }

  /* Sample the projected velocity at the same PLIC crossings.  Basilisk's
     interpolate_array performs the required MPI point-location reduction. */
  coord interface_points[2*bins];
  int sampled_bin[bins], sampled_bins = 0;
  for (int bin = 0; bin < bins; bin++) {
    if (!axisymmetric_film_crossing_pair_is_finite(
          lower_samples[bin], upper_samples[bin], lower[bin], upper[bin]))
      continue;
    double radius = (bin + .5)*support_half_width/bins;
    sampled_bin[sampled_bins] = bin;
    interface_points[2*sampled_bins] = (coord){lower[bin], radius};
    interface_points[2*sampled_bins + 1] = (coord){upper[bin], radius};
    sampled_bins++;
  }
  double compact_velocity[4*bins];
  double interface_velocity[4*bins];
  for (int component = 0; component < 4*bins; component++)
    interface_velocity[component] = nodata;
  if (sampled_bins > 0) {
    interpolate_array((scalar *){u.x, u.y}, interface_points, 2*sampled_bins,
                      compact_velocity, true);
    for (int sample = 0; sample < sampled_bins; sample++) {
      int bin = sampled_bin[sample];
      for (int component = 0; component < 4; component++)
        interface_velocity[4*bin + component] =
          compact_velocity[4*sample + component];
    }
  }

  AxisymmetricFilmSamples samples = {0};
  axisymmetric_film_samples_set_widths(&samples, support_half_width);
  samples.invalid_plic_segments = invalid_plic_segments;
  samples.skipped_degenerate_plic_facets = skipped_degenerate_plic_facets;
  double radial_step = support_half_width/bins;
  for (int bin = 0; bin < bins; bin++) {
    bool both_interfaces_sampled =
      axisymmetric_film_crossing_pair_is_finite(
        lower_samples[bin], upper_samples[bin], lower[bin], upper[bin]);
    if (both_interfaces_sampled && upper[bin] <= lower[bin]) {
      if (samples.first_nonpositive_gap_bin < 0) {
        samples.first_nonpositive_gap_bin = bin;
        samples.first_nonpositive_lower_samples = lower_samples[bin];
        samples.first_nonpositive_upper_samples = upper_samples[bin];
        samples.first_nonpositive_lower_position = lower[bin];
        samples.first_nonpositive_upper_position = upper[bin];
        samples.first_nonpositive_lower_delta = lower_delta[bin];
        samples.first_nonpositive_upper_delta = upper_delta[bin];
      }
      samples.nonpositive_gap_bins++;
    }
    bool valid = lower_samples[bin] > 0 && upper_samples[bin] > 0 &&
                 isfinite(lower[bin]) && isfinite(upper[bin]) &&
                 upper[bin] > lower[bin] &&
                 interface_velocity[4*bin] != nodata &&
                 interface_velocity[4*bin + 1] != nodata &&
                 interface_velocity[4*bin + 2] != nodata &&
                 interface_velocity[4*bin + 3] != nodata &&
                 isfinite(interface_velocity[4*bin]) &&
                 isfinite(interface_velocity[4*bin + 1]) &&
                 isfinite(interface_velocity[4*bin + 2]) &&
                 isfinite(interface_velocity[4*bin + 3]);
    AxisymmetricFilmSample * sample = &samples.bin[bin];
    sample->lower_sample_count = lower_samples[bin];
    sample->upper_sample_count = upper_samples[bin];
    sample->radius = (bin + .5)*radial_step;
    sample->width = radial_step;
    sample->valid = valid;
    if (both_interfaces_sampled) {
      sample->delta = max(lower_delta[bin], upper_delta[bin]);
      sample->normal_x = 1.;
      sample->normal_r = 0.;
      sample->lower_point = (AxisymmetricFilmPoint){
        .x = lower[bin], .r = sample->radius};
      sample->upper_point = (AxisymmetricFilmPoint){
        .x = upper[bin], .r = sample->radius};
    }
    if (valid) {
      sample->gap = upper[bin] - lower[bin];
      sample->lower_velocity = (AxisymmetricFilmVelocity){
        .x = interface_velocity[4*bin],
        .r = interface_velocity[4*bin + 1]};
      sample->upper_velocity = (AxisymmetricFilmVelocity){
        .x = interface_velocity[4*bin + 2],
        .r = interface_velocity[4*bin + 3]};
      sample->relative_normal_velocity =
        sample->lower_velocity.x - sample->upper_velocity.x;
    }
  }
  axisymmetric_film_samples_finalize(&samples);
  axisymmetric_film_samples_prepare_input(&samples, .90);
  return samples;
}

static inline AxisymmetricFilmDistribution
measure_axisymmetric_film_distribution (
  double support_half_width,
  double resolved_cell_threshold)
{
  AxisymmetricFilmSamples samples = sample_axisymmetric_water_film(
    support_half_width);
  double gaps[AXISYMMETRIC_FILM_RADIAL_BINS];
  double deltas[AXISYMMETRIC_FILM_RADIAL_BINS];
  double weights[AXISYMMETRIC_FILM_RADIAL_BINS];
  for (int bin = 0; bin < AXISYMMETRIC_FILM_RADIAL_BINS; bin++) {
    AxisymmetricFilmSample sample = samples.bin[bin];
    gaps[bin] = sample.valid ? sample.gap : NAN;
    deltas[bin] = sample.valid ? sample.delta : NAN;
    double inner_radius = bin*sample.width;
    double outer_radius = (bin + 1.)*sample.width;
    weights[bin] = sq(outer_radius) - sq(inner_radius);
  }
  AxisymmetricFilmDistribution distribution = {
    .support_half_width = support_half_width,
    .nonpositive_gap_bins = samples.nonpositive_gap_bins,
    .statistics = film_distribution_statistics(
      gaps, deltas, weights, AXISYMMETRIC_FILM_RADIAL_BINS,
      resolved_cell_threshold)
  };
  return distribution;
}

#endif
