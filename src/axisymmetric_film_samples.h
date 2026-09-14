#ifndef AXISYMMETRIC_FILM_SAMPLES_H
#define AXISYMMETRIC_FILM_SAMPLES_H

#include <float.h>
#include <math.h>
#include <stdbool.h>

#include "axisymmetric_matched_film.h"

typedef struct {
  double x;
  double r;
} AxisymmetricFilmPoint;

typedef struct {
  AxisymmetricFilmPoint first;
  AxisymmetricFilmPoint second;
} AxisymmetricFilmSegment;

typedef struct {
  double x;
  double r;
} AxisymmetricFilmVelocity;

typedef struct {
  bool valid;
  double radius;
  double width;
  double gap;
  double delta;
  double normal_x;
  double normal_r;
  AxisymmetricFilmPoint lower_point;
  AxisymmetricFilmPoint upper_point;
  AxisymmetricFilmVelocity lower_velocity;
  AxisymmetricFilmVelocity upper_velocity;
  double relative_normal_velocity;
} AxisymmetricFilmSample;

typedef struct {
  int bins;
  int valid_bins;
  int interpolated_bins;
  int nonpositive_gap_bins;
  int invalid_plic_segments;
  int skipped_degenerate_plic_facets;
  int first_nonpositive_gap_bin;
  int first_nonpositive_lower_samples;
  int first_nonpositive_upper_samples;
  double first_nonpositive_lower_position;
  double first_nonpositive_upper_position;
  double first_nonpositive_lower_delta;
  double first_nonpositive_upper_delta;
  double valid_area_fraction;
  bool input_ready;
  AxisymmetricFilmSample bin[AXISYMMETRIC_MATCHED_FILM_BINS];
  AxisymmetricFilmInput input;
} AxisymmetricFilmSamples;

static inline bool axisymmetric_film_facet_count_is_degenerate (int facets)
{
  return facets > 0 && facets != 2;
}

static inline bool axisymmetric_film_segment_is_finite (
  AxisymmetricFilmSegment segment)
{
  return isfinite(segment.first.x) && isfinite(segment.first.r) &&
    isfinite(segment.second.x) && isfinite(segment.second.r);
}

static inline bool axisymmetric_film_segment_crossing_at_radius (
  AxisymmetricFilmSegment segment, double radius,
  double * minimum_x, double * maximum_x)
{
  if (!minimum_x || !maximum_x || !isfinite(radius) ||
      !axisymmetric_film_segment_is_finite(segment))
    return false;

  const double lower_radius = fmin(segment.first.r, segment.second.r);
  const double upper_radius = fmax(segment.first.r, segment.second.r);
  const double fraction_tolerance = 32.*DBL_EPSILON;
  if (radius < lower_radius - 32.*DBL_EPSILON*
        fmax(fabs(lower_radius), fabs(upper_radius)) ||
      radius > upper_radius + 32.*DBL_EPSILON*
        fmax(fabs(lower_radius), fabs(upper_radius)))
    return false;

  const double radial_difference = segment.second.r - segment.first.r;
  if (!isfinite(radial_difference))
    return false;
  if (fabs(radial_difference) <= 32.*DBL_EPSILON*
      fmax(fabs(segment.first.r), fabs(segment.second.r))) {
    if (fabs(radius - segment.first.r) > 32.*DBL_EPSILON*
        fmax(fabs(segment.first.r), fabs(segment.second.r)))
      return false;
    *minimum_x = fmin(segment.first.x, segment.second.x);
    *maximum_x = fmax(segment.first.x, segment.second.x);
    return isfinite(*minimum_x) && isfinite(*maximum_x);
  }

  double fraction = (radius - segment.first.r)/radial_difference;
  if (!isfinite(fraction) || fraction < -fraction_tolerance ||
      fraction > 1. + fraction_tolerance)
    return false;
  fraction = fmin(1., fmax(0., fraction));
  const double axial_difference = segment.second.x - segment.first.x;
  if (!isfinite(axial_difference))
    return false;
  const double crossing = segment.first.x + fraction*axial_difference;
  if (!isfinite(crossing))
    return false;
  *minimum_x = *maximum_x = crossing;
  return true;
}

static inline bool axisymmetric_film_crossing_pair_is_finite (
  int lower_samples, int upper_samples,
  double lower_position, double upper_position)
{
  return lower_samples > 0 && upper_samples > 0 &&
    isfinite(lower_position) && isfinite(upper_position);
}

static inline bool axisymmetric_film_segment_point_at_radius (
  AxisymmetricFilmSegment segment, double radius,
  AxisymmetricFilmPoint * point)
{
  if (!point || !isfinite(radius) ||
      !axisymmetric_film_segment_is_finite(segment))
    return false;
  double lower_radius = fmin(segment.first.r, segment.second.r);
  double upper_radius = fmax(segment.first.r, segment.second.r);
  if (radius < lower_radius - 32.*DBL_EPSILON*
        fmax(fabs(segment.first.r), fabs(segment.second.r)) ||
      radius > upper_radius + 32.*DBL_EPSILON*
        fmax(fabs(segment.first.r), fabs(segment.second.r)))
    return false;
  double difference = segment.second.r - segment.first.r;
  if (fabs(difference) <= 32.*DBL_EPSILON*
      fmax(fabs(segment.first.r), fabs(segment.second.r))) {
    if (fabs(radius - segment.first.r) > 32.*DBL_EPSILON*
        fmax(fabs(segment.first.r), fabs(segment.second.r)))
      return false;
    point->x = .5*(segment.first.x + segment.second.x);
    point->r = radius;
    return true;
  }
  double fraction = (radius - segment.first.r)/difference;
  point->x = segment.first.x + fraction*(segment.second.x - segment.first.x);
  point->r = radius;
  return isfinite(point->x);
}

static inline bool axisymmetric_film_sample_from_segments (
  AxisymmetricFilmSegment lower_segment,
  AxisymmetricFilmSegment upper_segment,
  double radius, double lower_delta, double upper_delta,
  AxisymmetricFilmVelocity lower_velocity,
  AxisymmetricFilmVelocity upper_velocity,
  AxisymmetricFilmSample * sample)
{
  if (!sample || !isfinite(radius) || radius < 0. ||
      !isfinite(lower_delta) || lower_delta <= 0. ||
      !isfinite(upper_delta) || upper_delta <= 0. ||
      !isfinite(lower_velocity.x) || !isfinite(lower_velocity.r) ||
      !isfinite(upper_velocity.x) || !isfinite(upper_velocity.r))
    return false;
  AxisymmetricFilmPoint lower_point, upper_point;
  if (!axisymmetric_film_segment_point_at_radius(
        lower_segment, radius, &lower_point) ||
      !axisymmetric_film_segment_point_at_radius(
        upper_segment, radius, &upper_point))
    return false;
  double normal_x = upper_point.x - lower_point.x;
  double normal_r = upper_point.r - lower_point.r;
  double gap = hypot(normal_x, normal_r);
  if (!isfinite(gap) || gap <= 0. || normal_x <= 0.)
    return false;
  normal_x /= gap;
  normal_r /= gap;
  *sample = (AxisymmetricFilmSample){
    .valid = true,
    .radius = radius,
    .gap = gap,
    .delta = fmax(lower_delta, upper_delta),
    .normal_x = normal_x,
    .normal_r = normal_r,
    .lower_point = lower_point,
    .upper_point = upper_point,
    .lower_velocity = lower_velocity,
    .upper_velocity = upper_velocity,
    .relative_normal_velocity =
      (lower_velocity.x - upper_velocity.x)*normal_x +
      (lower_velocity.r - upper_velocity.r)*normal_r
  };
  return true;
}

static inline void axisymmetric_film_samples_set_widths (
  AxisymmetricFilmSamples * samples, double support_half_width)
{
  samples->bins = AXISYMMETRIC_MATCHED_FILM_BINS;
  samples->first_nonpositive_gap_bin = -1;
  samples->first_nonpositive_lower_position = NAN;
  samples->first_nonpositive_upper_position = NAN;
  samples->first_nonpositive_lower_delta = NAN;
  samples->first_nonpositive_upper_delta = NAN;
  double width = support_half_width/AXISYMMETRIC_MATCHED_FILM_BINS;
  for (int bin = 0; bin < AXISYMMETRIC_MATCHED_FILM_BINS; bin++) {
    samples->bin[bin].radius = (bin + .5)*width;
    samples->bin[bin].width = width;
    samples->input.radius[bin] = samples->bin[bin].radius;
    samples->input.width[bin] = width;
  }
}

static inline void axisymmetric_film_samples_finalize (
  AxisymmetricFilmSamples * samples)
{
  samples->valid_bins = 0;
  samples->interpolated_bins = 0;
  samples->input_ready = false;
  double total_area = 0., valid_area = 0.;
  for (int bin = 0; bin < AXISYMMETRIC_MATCHED_FILM_BINS; bin++) {
    AxisymmetricFilmSample * sample = &samples->bin[bin];
    double inner_radius = bin*sample->width;
    double outer_radius = (bin + 1.)*sample->width;
    double area = axisymmetric_matched_film_square(outer_radius) -
      axisymmetric_matched_film_square(inner_radius);
    total_area += area;
    if (sample->valid) {
      samples->valid_bins++;
      valid_area += area;
      samples->input.gap[bin] = sample->gap;
      samples->input.delta[bin] = sample->delta;
      samples->input.lower_position[bin] = sample->lower_point.x;
      samples->input.upper_position[bin] = sample->upper_point.x;
      samples->input.relative_normal_velocity[bin] =
        sample->relative_normal_velocity;
      samples->input.lower_tangential_velocity[bin] =
        -sample->lower_velocity.x*sample->normal_r +
        sample->lower_velocity.r*sample->normal_x;
      samples->input.upper_tangential_velocity[bin] =
        -sample->upper_velocity.x*sample->normal_r +
        sample->upper_velocity.r*sample->normal_x;
    }
    else {
      samples->input.gap[bin] = NAN;
      samples->input.delta[bin] = NAN;
      samples->input.lower_position[bin] = NAN;
      samples->input.upper_position[bin] = NAN;
      samples->input.relative_normal_velocity[bin] = NAN;
      samples->input.lower_tangential_velocity[bin] = NAN;
      samples->input.upper_tangential_velocity[bin] = NAN;
    }
  }
  samples->valid_area_fraction = total_area > 0. ? valid_area/total_area : 0.;
}

static inline double axisymmetric_film_samples_linear_value (
  double lower_radius, double lower_value,
  double upper_radius, double upper_value, double radius)
{
  double fraction = (radius - lower_radius)/(upper_radius - lower_radius);
  return lower_value + fraction*(upper_value - lower_value);
}

/* Fill only isolated, bracketed PLIC sampling holes.  Missing edge support is
   not extrapolated because it does not define a closed Reynolds domain. */
static inline bool axisymmetric_film_samples_prepare_input (
  AxisymmetricFilmSamples * samples, double minimum_valid_area_fraction)
{
  if (!samples || !isfinite(minimum_valid_area_fraction) ||
      minimum_valid_area_fraction <= 0. ||
      minimum_valid_area_fraction > 1. ||
      samples->nonpositive_gap_bins > 0 ||
      samples->invalid_plic_segments > 0 ||
      samples->valid_area_fraction < minimum_valid_area_fraction)
    return false;

  for (int bin = 0; bin < AXISYMMETRIC_MATCHED_FILM_BINS; bin++) {
    if (samples->bin[bin].valid)
      continue;
    int lower = bin - 1;
    int upper = bin + 1;
    while (lower >= 0 && !samples->bin[lower].valid)
      lower--;
    while (upper < AXISYMMETRIC_MATCHED_FILM_BINS &&
           !samples->bin[upper].valid)
      upper++;
    if (lower < 0 || upper >= AXISYMMETRIC_MATCHED_FILM_BINS)
      return false;

    const AxisymmetricFilmSample * first = &samples->bin[lower];
    const AxisymmetricFilmSample * second = &samples->bin[upper];
    double radius = samples->bin[bin].radius;
#define INTERPOLATE_INPUT(member) \
    axisymmetric_film_samples_linear_value( \
      first->radius, samples->input.member[lower], \
      second->radius, samples->input.member[upper], radius)
    samples->input.gap[bin] = INTERPOLATE_INPUT(gap);
    samples->input.delta[bin] = INTERPOLATE_INPUT(delta);
    samples->input.lower_position[bin] =
      INTERPOLATE_INPUT(lower_position);
    samples->input.upper_position[bin] =
      INTERPOLATE_INPUT(upper_position);
    samples->input.relative_normal_velocity[bin] =
      INTERPOLATE_INPUT(relative_normal_velocity);
    samples->input.lower_tangential_velocity[bin] =
      INTERPOLATE_INPUT(lower_tangential_velocity);
    samples->input.upper_tangential_velocity[bin] =
      INTERPOLATE_INPUT(upper_tangential_velocity);
#undef INTERPOLATE_INPUT
    if (!isfinite(samples->input.gap[bin]) ||
        samples->input.gap[bin] <= 0. ||
        !isfinite(samples->input.delta[bin]) ||
        samples->input.delta[bin] <= 0. ||
        !isfinite(samples->input.lower_position[bin]) ||
        !isfinite(samples->input.upper_position[bin]) ||
        samples->input.upper_position[bin] <=
          samples->input.lower_position[bin] ||
        !isfinite(samples->input.relative_normal_velocity[bin]) ||
        !isfinite(samples->input.lower_tangential_velocity[bin]) ||
        !isfinite(samples->input.upper_tangential_velocity[bin]))
      return false;
    samples->interpolated_bins++;
  }
  samples->input_ready = true;
  return true;
}

#endif
