#ifndef MATCHED_FILM_SAMPLING_POLICY_H
#define MATCHED_FILM_SAMPLING_POLICY_H

#include <math.h>
#include <stdbool.h>

typedef enum {
  MATCHED_FILM_SAMPLING_INVALID = 0,
  MATCHED_FILM_SAMPLING_NOT_REQUIRED = 1,
  MATCHED_FILM_SAMPLING_REQUIRED = 2,
  MATCHED_FILM_SAMPLING_RELEASED = 3
} MatchedFilmSamplingDecision;

static inline MatchedFilmSamplingDecision matched_film_sampling_decision (
  bool film_active, bool contact_valid, double contact_gap,
  double finest_delta, double activation_cells, double match_cells,
  double release_factor)
{
  if (!isfinite(finest_delta) || finest_delta <= 0. ||
      !isfinite(activation_cells) || activation_cells <= 0. ||
      !isfinite(match_cells) || match_cells <= 0. ||
      !isfinite(release_factor) || release_factor < 1. ||
      (contact_valid && !isfinite(contact_gap)))
    return MATCHED_FILM_SAMPLING_INVALID;
  if (film_active && contact_valid &&
      contact_gap > release_factor*match_cells*finest_delta)
    return MATCHED_FILM_SAMPLING_RELEASED;
  if (film_active ||
      (contact_valid && contact_gap <= activation_cells*finest_delta))
    return MATCHED_FILM_SAMPLING_REQUIRED;
  return MATCHED_FILM_SAMPLING_NOT_REQUIRED;
}

#endif
