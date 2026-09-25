#ifndef FILM_MODEL_MUSEHANE_CHECKPOINT_FILE_H
#define FILM_MODEL_MUSEHANE_CHECKPOINT_FILE_H

#include <stdint.h>

#include "film_model/musehane.h"

typedef enum {
  MUSEHANE_CHECKPOINT_FILE_OK = 0,
  MUSEHANE_CHECKPOINT_FILE_INVALID_INPUT = 1,
  MUSEHANE_CHECKPOINT_FILE_IO_FAILED = 2,
  MUSEHANE_CHECKPOINT_FILE_FORMAT_MISMATCH = 3,
  MUSEHANE_CHECKPOINT_FILE_CHECKSUM_MISMATCH = 4,
  MUSEHANE_CHECKPOINT_FILE_MODEL_FAILED = 5
} MusehaneCheckpointFileStatus;

typedef struct {
  double time;
  unsigned long iteration;
  unsigned long accepted_steps;
  size_t payload_size;
  uint64_t payload_checksum;
} MusehaneCheckpointFileAudit;

MusehaneCheckpointFileStatus musehane_checkpoint_file_write (
  const char *path, const FilmModelState *state, double time,
  unsigned long iteration, unsigned long accepted_steps,
  MusehaneCheckpointFileAudit *audit);

MusehaneCheckpointFileStatus musehane_checkpoint_file_restore (
  const char *path, const FilmModelConfig *config, double expected_time,
  unsigned long expected_iteration, FilmModelState *state,
  unsigned long *accepted_steps, MusehaneCheckpointFileAudit *audit);

#endif
