#define TRACK_A_MUSEHANE_TRACTION_POLICY \
  MUSEHANE_PUBLISHED_EQUATION27_TRACTION
#define TRACK_A_MUSEHANE_REPLACE_NATIVE_VISCOSITY 1
#define INITIAL_UPPER_TRANSLATION_CONDITIONING 1
static void track_a_coupled_replace_native_viscosity (void);
#define DUAL_COMPOUND_PROPERTIES_EXTENSION() \
  track_a_coupled_replace_native_viscosity()
#define COLLISION_TRACTION_EXTENSION_HEADER \
  track_a_musehane_coupled_axi.h
#include "collision_initialization_axi.c"
