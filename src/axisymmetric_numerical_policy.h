#ifndef AXISYMMETRIC_NUMERICAL_POLICY_H
#define AXISYMMETRIC_NUMERICAL_POLICY_H

/*
 * The tolerance is a dimensional pressure-projection residual in Basilisk's
 * centered solver. 1e-11 is the least strict value that passed the corrected
 * axisymmetric cell-volume audit at L10 without the 1e-12 iteration stall.
 */
#ifndef PRESSURE_TOLERANCE
# define PRESSURE_TOLERANCE 1.e-11
#endif

#endif
