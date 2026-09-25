#include "film_musehane/activation.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

MusehaneActivationStatus musehane_evaluate_activation (
  const MusehaneActivationConfig *config,
  const MusehaneActivationInput *input,
  MusehaneActivationOutput *output)
{
  if (!config || !input || !output || input->cell_count == 0 ||
      input->cell_count > SIZE_MAX/sizeof(unsigned char) ||
      !isfinite(config->activation_cells) || config->activation_cells <= 0. ||
      !isfinite(config->minimum_smoothed_fraction) ||
      config->minimum_smoothed_fraction < 0. ||
      config->minimum_smoothed_fraction > 1. ||
      !input->estimated_separation || !input->local_delta ||
      !input->lower_smoothed_fraction || !input->upper_smoothed_fraction ||
      !output->active)
    return MUSEHANE_ACTIVATION_INVALID_INPUT;

  unsigned char *candidate = calloc(input->cell_count, sizeof *candidate);
  if (!candidate)
    return MUSEHANE_ACTIVATION_ALLOCATION_FAILED;
  size_t active_cells = 0;
  for (size_t cell = 0; cell < input->cell_count; cell++) {
    const double h0 = input->estimated_separation[cell];
    const double delta = input->local_delta[cell];
    const double lower = input->lower_smoothed_fraction[cell];
    const double upper = input->upper_smoothed_fraction[cell];
    if (!isfinite(h0) || h0 < 0. || !isfinite(delta) || delta <= 0. ||
        !isfinite(lower) || lower < 0. || lower > 1. ||
        !isfinite(upper) || upper < 0. || upper > 1.) {
      free(candidate);
      return MUSEHANE_ACTIVATION_INVALID_INPUT;
    }
    candidate[cell] = h0 <= config->activation_cells*delta &&
      lower >= config->minimum_smoothed_fraction &&
      upper >= config->minimum_smoothed_fraction;
    active_cells += candidate[cell] != 0;
  }

  memcpy(output->active, candidate, input->cell_count*sizeof *candidate);
  output->active_cells = active_cells;
  free(candidate);
  return MUSEHANE_ACTIVATION_OK;
}
