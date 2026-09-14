#ifndef AXISYMMETRIC_CONTACT_GEOMETRY_H
#define AXISYMMETRIC_CONTACT_GEOMETRY_H

typedef struct {
  double lower_top;
  double upper_bottom;
  double gap;
  double lower_max_delta;
  double upper_max_delta;
  int lower_samples;
  int upper_samples;
  bool valid;
} AxisymmetricContactGeometry;

typedef struct {
  double minimum_x;
  double maximum_x;
  int samples;
} ContactSegmentSupportExtrema;

static inline ContactSegmentSupportExtrema contact_segment_support_extrema (
  coord first,
  coord second,
  double half_width,
  double minimum_x,
  double maximum_x,
  int samples)
{
  ContactSegmentSupportExtrema result = {
    .minimum_x = minimum_x,
    .maximum_x = maximum_x,
    .samples = samples
  };
  if (first.y <= half_width) {
    result.minimum_x = min(result.minimum_x, first.x);
    result.maximum_x = max(result.maximum_x, first.x);
    result.samples++;
  }
  if (second.y <= half_width) {
    result.minimum_x = min(result.minimum_x, second.x);
    result.maximum_x = max(result.maximum_x, second.x);
    result.samples++;
  }
  if ((first.y - half_width)*(second.y - half_width) < 0.) {
    double interpolation = (half_width - first.y)/(second.y - first.y);
    double crossing_x = first.x + interpolation*(second.x - first.x);
    result.minimum_x = min(result.minimum_x, crossing_x);
    result.maximum_x = max(result.maximum_x, crossing_x);
    result.samples++;
  }
  return result;
}

static inline AxisymmetricContactGeometry measure_axisymmetric_contact_geometry (
  double half_width)
{
  double lower_top = -HUGE_VAL*L0;
  double upper_bottom = HUGE_VAL*L0;
  int lower_samples = 0, upper_samples = 0;
  double lower_max_delta = 0., upper_max_delta = 0.;
  face vector no_face_fractions = {{-1}};

  foreach(reduction(max:lower_top) reduction(min:upper_bottom)
          reduction(max:lower_max_delta) reduction(max:upper_max_delta)
          reduction(+:lower_samples) reduction(+:upper_samples)) {
    if (lower_envelope_fraction[] > 1.e-6 &&
        lower_envelope_fraction[] < 1. - 1.e-6) {
      coord normal = facet_normal(
        point, lower_envelope_fraction, no_face_fractions);
      double alpha = plane_alpha(lower_envelope_fraction[], normal);
      coord segment[2];
      if (facets(normal, alpha, segment) == 2) {
        for (int endpoint = 0; endpoint < 2; endpoint++) {
          segment[endpoint].x = x + segment[endpoint].x*Delta;
          segment[endpoint].y = y + segment[endpoint].y*Delta;
        }
        double unused_minimum = HUGE_VAL*L0;
        int samples_before = lower_samples;
        ContactSegmentSupportExtrema support = contact_segment_support_extrema(
          segment[0], segment[1], half_width,
          unused_minimum, lower_top, lower_samples);
        lower_top = support.maximum_x;
        lower_samples = support.samples;
        if (lower_samples > samples_before)
          lower_max_delta = max(lower_max_delta, Delta);
      }
    }
    if (upper_envelope_fraction[] > 1.e-6 &&
        upper_envelope_fraction[] < 1. - 1.e-6) {
      coord normal = facet_normal(
        point, upper_envelope_fraction, no_face_fractions);
      double alpha = plane_alpha(upper_envelope_fraction[], normal);
      coord segment[2];
      if (facets(normal, alpha, segment) == 2) {
        for (int endpoint = 0; endpoint < 2; endpoint++) {
          segment[endpoint].x = x + segment[endpoint].x*Delta;
          segment[endpoint].y = y + segment[endpoint].y*Delta;
        }
        double unused_maximum = -HUGE_VAL*L0;
        int samples_before = upper_samples;
        ContactSegmentSupportExtrema support = contact_segment_support_extrema(
          segment[0], segment[1], half_width,
          upper_bottom, unused_maximum, upper_samples);
        upper_bottom = support.minimum_x;
        upper_samples = support.samples;
        if (upper_samples > samples_before)
          upper_max_delta = max(upper_max_delta, Delta);
      }
    }
  }

  AxisymmetricContactGeometry geometry = {
    .lower_top = lower_top,
    .upper_bottom = upper_bottom,
    .gap = NAN,
    .lower_max_delta = lower_max_delta,
    .upper_max_delta = upper_max_delta,
    .lower_samples = lower_samples,
    .upper_samples = upper_samples,
    .valid = lower_samples > 0 && upper_samples > 0 &&
             isfinite(lower_top) && isfinite(upper_bottom)
  };
  if (geometry.valid)
    geometry.gap = upper_bottom - lower_top;
  return geometry;
}

#endif
