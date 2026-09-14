#ifndef SEPARATED_SOLVER_TOLERANCE_POLICY_H
#define SEPARATED_SOLVER_TOLERANCE_POLICY_H

#ifndef VISCOUS_TOLERANCE
# define VISCOUS_TOLERANCE 1.e-5
#endif

/* Basilisk's viscosity solver reads the global TOLERANCE directly, while
   project() derives its own threshold from the same value. Event inheritance
   lets each solver see its registered tolerance without forking centered.h. */
event viscous_term (i++, last)
{
  TOLERANCE = (double)VISCOUS_TOLERANCE;
}

event acceleration (i++, last)
{
  TOLERANCE = (double)PRESSURE_TOLERANCE;
}

#endif
