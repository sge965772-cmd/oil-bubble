#ifndef TRACK_A_NO_STATEFUL_FILM_RESTART_H
#define TRACK_A_NO_STATEFUL_FILM_RESTART_H

/* Compatibility include for the collision solver's full-state DNS restart.
   Track A persists its own dynamic state in the adjacent .track-a.bin file;
   the legacy Reynolds state is deliberately disabled. */
#define COLLISION_RESTART_FILM_STATE_V4_INCLUDED 1
#define COLLISION_RESTART_FILM_STATE_SHA256 "not-applicable-track-a"
#define COLLISION_RESTART_CONTRACT_SHA256 "track-a-adjacent-sidecar-v1"
#define COLLISION_RESTART_FILM_STATE_TIME 0.

#endif
