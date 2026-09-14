#ifndef AXISYMMETRIC_MATCHED_FILM_AXI_H
#define AXISYMMETRIC_MATCHED_FILM_AXI_H

#include "axisymmetric_matched_film.h"

typedef struct {
  bool applied;
  double lower_axial_force;
  double upper_axial_force;
  double net_axial_force;
  double lower_radial_force_moment;
  double upper_radial_force_moment;
  double net_axial_moment;
  double lower_support_volume;
  double upper_support_volume;
  double maximum_bin_force_imbalance;
  int missing_support_bins;
} AxisymmetricMatchedFilmPairLedger;

static inline double axisymmetric_matched_film_face_kernel (
  double face_position, double interface_position, double delta)
{
  if (!isfinite(face_position) || !isfinite(interface_position) ||
      !isfinite(delta) || delta <= 0.)
    return 0.;
  double distance = fabs(face_position - interface_position)/delta;
  return distance < 1. ? 1. - distance : 0.;
}

static inline AxisymmetricMatchedFilmPairLedger
axisymmetric_matched_film_pair_ledger (
  const AxisymmetricMatchedFilmResult * result)
{
  AxisymmetricMatchedFilmPairLedger ledger = {0};
  if (!result)
    return ledger;
  for (int bin = 0; bin < AXISYMMETRIC_MATCHED_FILM_BINS; bin++) {
    ledger.lower_axial_force -= result->ring_force[bin];
    ledger.upper_axial_force += result->ring_force[bin];
  }
  ledger.net_axial_force = ledger.lower_axial_force + ledger.upper_axial_force;
  return ledger;
}

#ifdef AXISYMMETRIC_MATCHED_FILM_BASILISK_ADAPTER
static inline AxisymmetricMatchedFilmPairLedger
axisymmetric_matched_film_apply_axi (
  face vector acceleration_field,
  const AxisymmetricMatchedFilmResult * result,
  const AxisymmetricFilmInput * input,
  double support_half_width)
{
  AxisymmetricMatchedFilmPairLedger ledger =
    axisymmetric_matched_film_pair_ledger(result);
  if (!result || !input || !result->valid || !result->active)
    return ledger;

  enum { bins = AXISYMMETRIC_MATCHED_FILM_BINS };
  double lower_support[bins], upper_support[bins];
  for (int bin = 0; bin < bins; bin++)
    lower_support[bin] = upper_support[bin] = 0.;
  foreach_face(x, reduction(+:lower_support[:bins])
                  reduction(+:upper_support[:bins])) {
    if (y < 0. || y >= support_half_width + Delta)
      continue;
    double physical_dual_volume = 2.*pi*fm.x[]*sq(Delta);
    for (int bin = 0; bin < bins; bin++) {
      double radial_weight = axisymmetric_matched_film_face_kernel(
        y, input->radius[bin], Delta);
      double lower_weight = radial_weight*
        axisymmetric_matched_film_face_kernel(
          x, input->lower_position[bin], Delta);
      double upper_weight = radial_weight*
        axisymmetric_matched_film_face_kernel(
          x, input->upper_position[bin], Delta);
      lower_support[bin] += physical_dual_volume*lower_weight;
      upper_support[bin] += physical_dual_volume*upper_weight;
    }
  }
  ledger.lower_support_volume = 0.;
  ledger.upper_support_volume = 0.;
  for (int bin = 0; bin < bins; bin++) {
    ledger.lower_support_volume += lower_support[bin];
    ledger.upper_support_volume += upper_support[bin];
    if (result->ring_force[bin] > 0. &&
        (lower_support[bin] <= 0. || upper_support[bin] <= 0.))
      ledger.missing_support_bins++;
  }
  if (ledger.missing_support_bins)
    return ledger;

  double applied_lower = 0., applied_upper = 0.;
  double lower_radial_moment = 0., upper_radial_moment = 0.;
  double lower_by_bin[bins], upper_by_bin[bins];
  for (int bin = 0; bin < bins; bin++)
    lower_by_bin[bin] = upper_by_bin[bin] = 0.;
  foreach_face(x, reduction(+:applied_lower) reduction(+:applied_upper)
                  reduction(+:lower_radial_moment)
                  reduction(+:upper_radial_moment)
                  reduction(+:lower_by_bin[:bins])
                  reduction(+:upper_by_bin[:bins])) {
    if (y < 0. || y >= support_half_width + Delta)
      continue;
    double physical_dual_volume = 2.*pi*fm.x[]*sq(Delta);
    double lower_density_force = 0., upper_density_force = 0.;
    for (int bin = 0; bin < bins; bin++) {
      double radial_weight = axisymmetric_matched_film_face_kernel(
        y, input->radius[bin], Delta);
      double lower_weight = radial_weight*
        axisymmetric_matched_film_face_kernel(
          x, input->lower_position[bin], Delta);
      double upper_weight = radial_weight*
        axisymmetric_matched_film_face_kernel(
          x, input->upper_position[bin], Delta);
      double lower_bin_force = lower_support[bin] > 0. ?
        -result->ring_force[bin]/lower_support[bin]*lower_weight : 0.;
      double upper_bin_force = upper_support[bin] > 0. ?
        result->ring_force[bin]/upper_support[bin]*upper_weight : 0.;
      lower_density_force += lower_bin_force;
      upper_density_force += upper_bin_force;
      lower_by_bin[bin] += lower_bin_force*physical_dual_volume;
      upper_by_bin[bin] += upper_bin_force*physical_dual_volume;
    }
    acceleration_field.x[] += dual_compound_alpha.x[]/(fm.x[] + SEPS)*
      (lower_density_force + upper_density_force);
    double lower_force = lower_density_force*physical_dual_volume;
    double upper_force = upper_density_force*physical_dual_volume;
    applied_lower += lower_force;
    applied_upper += upper_force;
    lower_radial_moment += y*lower_force;
    upper_radial_moment += y*upper_force;
  }
  ledger.lower_axial_force = applied_lower;
  ledger.upper_axial_force = applied_upper;
  ledger.net_axial_force = applied_lower + applied_upper;
  ledger.lower_radial_force_moment = lower_radial_moment;
  ledger.upper_radial_force_moment = upper_radial_moment;
  ledger.net_axial_moment = lower_radial_moment + upper_radial_moment;
  for (int bin = 0; bin < bins; bin++)
    ledger.maximum_bin_force_imbalance = max(
      ledger.maximum_bin_force_imbalance,
      fabs(lower_by_bin[bin] + upper_by_bin[bin]));
  ledger.applied = true;
  return ledger;
}
#endif

#endif
