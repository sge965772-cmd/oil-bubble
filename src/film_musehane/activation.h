#ifndef FILM_MUSEHANE_ACTIVATION_H
#define FILM_MUSEHANE_ACTIVATION_H

#include <stddef.h>

typedef enum {
  MUSEHANE_ACTIVATION_OK = 0,
  MUSEHANE_ACTIVATION_INVALID_INPUT = 1,
  MUSEHANE_ACTIVATION_ALLOCATION_FAILED = 2
} MusehaneActivationStatus;

typedef struct {
  double activation_cells;
  double minimum_smoothed_fraction;
} MusehaneActivationConfig;

typedef struct {
  size_t cell_count;
  const double *estimated_separation;
  const double *local_delta;
  const double *lower_smoothed_fraction;
  const double *upper_smoothed_fraction;
} MusehaneActivationInput;

typedef struct {
  unsigned char *active;
  size_t active_cells;
} MusehaneActivationOutput;

/* Musehane et al. (2018), equation (22), including H(0)=1. */
MusehaneActivationStatus musehane_evaluate_activation (
  const MusehaneActivationConfig *config,
  const MusehaneActivationInput *input,
  MusehaneActivationOutput *output);

#endif
