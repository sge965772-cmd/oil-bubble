#ifndef FILM_DISTRIBUTION_STATISTICS_H
#define FILM_DISTRIBUTION_STATISTICS_H

#include <math.h>

#ifndef FILM_DISTRIBUTION_MAX_SAMPLES
# define FILM_DISTRIBUTION_MAX_SAMPLES 128
#endif

typedef struct {
  int total_bins;
  int valid_bins;
  double valid_area_fraction;
  double minimum_gap;
  double minimum_cells;
  double p10_cells;
  double median_cells;
  double resolved_area_fraction;
} FilmDistributionStatistics;

typedef struct {
  double cells;
  double weight;
} FilmDistributionWeightedSample;

static inline double film_distribution_weighted_quantile (
  const FilmDistributionWeightedSample * samples,
  int count,
  double total_weight,
  double quantile)
{
  double threshold, cumulative;
  int index;
  if (count <= 0 || !(total_weight > 0.) || !isfinite(total_weight))
    return NAN;
  threshold = quantile*total_weight;
  cumulative = 0.;
  for (index = 0; index < count; index++) {
    cumulative += samples[index].weight;
    if (cumulative >= threshold)
      return samples[index].cells;
  }
  return samples[count - 1].cells;
}

static inline FilmDistributionStatistics film_distribution_statistics (
  const double * gaps,
  const double * deltas,
  const double * area_weights,
  int count,
  double resolved_cell_threshold)
{
  FilmDistributionStatistics result;
  FilmDistributionWeightedSample samples[FILM_DISTRIBUTION_MAX_SAMPLES];
  FilmDistributionWeightedSample sample;
  double total_weight = 0., valid_weight = 0., resolved_weight = 0.;
  double weight, cells;
  int index, position;
  result.total_bins = count;
  result.valid_bins = 0;
  result.valid_area_fraction = 0.;
  result.minimum_gap = NAN;
  result.minimum_cells = NAN;
  result.p10_cells = NAN;
  result.median_cells = NAN;
  result.resolved_area_fraction = NAN;
  if (!gaps || !deltas || !area_weights || count <= 0 ||
      count > FILM_DISTRIBUTION_MAX_SAMPLES ||
      !(resolved_cell_threshold > 0.))
    return result;

  for (index = 0; index < count; index++) {
    weight = area_weights[index];
    if (isfinite(weight) && weight > 0.)
      total_weight += weight;
    if (!isfinite(gaps[index]) || !isfinite(deltas[index]) ||
        !(deltas[index] > 0.) || !isfinite(weight) || !(weight > 0.))
      continue;
    cells = gaps[index]/deltas[index];
    sample.cells = cells;
    sample.weight = weight;
    position = result.valid_bins;
    while (position > 0 && samples[position - 1].cells > sample.cells) {
      samples[position] = samples[position - 1];
      position--;
    }
    samples[position] = sample;
    result.valid_bins = result.valid_bins + 1;
    valid_weight += weight;
    if (cells >= resolved_cell_threshold)
      resolved_weight += weight;
    result.minimum_gap = !isfinite(result.minimum_gap) ? gaps[index] :
      fmin(result.minimum_gap, gaps[index]);
    result.minimum_cells = !isfinite(result.minimum_cells) ? cells :
      fmin(result.minimum_cells, cells);
  }

  result.valid_area_fraction = total_weight > 0. ?
    valid_weight/total_weight : 0.;
  if (result.valid_bins > 0 && valid_weight > 0.) {
    result.p10_cells = film_distribution_weighted_quantile(
      samples, result.valid_bins, valid_weight, 0.10);
    result.median_cells = film_distribution_weighted_quantile(
      samples, result.valid_bins, valid_weight, 0.50);
    result.resolved_area_fraction = resolved_weight/valid_weight;
  }
  return result;
}

#endif
