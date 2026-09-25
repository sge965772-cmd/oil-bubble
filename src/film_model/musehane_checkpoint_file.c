#include "film_model/musehane_checkpoint_file.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MUSEHANE_CHECKPOINT_FILE_MAGIC UINT64_C(0x4d55534548414e45)
#define MUSEHANE_CHECKPOINT_FILE_SCHEMA 1U

typedef struct {
  uint64_t magic;
  uint32_t schema;
  uint32_t restart_schema;
  uint32_t model_kind;
  uint32_t reserved;
  uint64_t payload_size;
  uint64_t iteration;
  uint64_t accepted_steps;
  uint64_t payload_checksum;
  double time;
} MusehaneCheckpointFileHeader;

static uint64_t checksum_bytes (const unsigned char *bytes, size_t size)
{
  uint64_t hash = UINT64_C(14695981039346656037);
  for (size_t index = 0; index < size; index++) {
    hash ^= bytes[index];
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

static int time_matches (double actual, double expected)
{
  const double scale = fmax(1., fmax(fabs(actual), fabs(expected)));
  return isfinite(actual) && isfinite(expected) &&
    fabs(actual - expected) <= 128.*DBL_EPSILON*scale;
}

MusehaneCheckpointFileStatus musehane_checkpoint_file_write (
  const char *path, const FilmModelState *state, double time,
  unsigned long iteration, unsigned long accepted_steps,
  MusehaneCheckpointFileAudit *audit)
{
  if (!path || !path[0] || !state || !isfinite(time) || time < 0. ||
      !audit)
    return MUSEHANE_CHECKPOINT_FILE_INVALID_INPUT;
  const size_t payload_size = film_model_checkpoint_size(state);
  if (payload_size == 0)
    return MUSEHANE_CHECKPOINT_FILE_MODEL_FAILED;
  unsigned char *payload = malloc(payload_size);
  if (!payload)
    return MUSEHANE_CHECKPOINT_FILE_MODEL_FAILED;
  FilmRestartRecord record = {
    .capacity = payload_size,
    .data = payload
  };
  if (film_model_checkpoint(state, &record) != FILM_MODEL_OK ||
      record.size != payload_size) {
    free(payload);
    return MUSEHANE_CHECKPOINT_FILE_MODEL_FAILED;
  }
  const uint64_t checksum = checksum_bytes(payload, payload_size);
  MusehaneCheckpointFileHeader header;
  memset(&header, 0, sizeof header);
  header.magic = MUSEHANE_CHECKPOINT_FILE_MAGIC;
  header.schema = MUSEHANE_CHECKPOINT_FILE_SCHEMA;
  header.restart_schema = record.schema_version;
  header.model_kind = record.model_kind;
  header.payload_size = payload_size;
  header.iteration = iteration;
  header.accepted_steps = accepted_steps;
  header.payload_checksum = checksum;
  header.time = time;

  const size_t path_length = strlen(path);
  if (path_length > SIZE_MAX - 5) {
    free(payload);
    return MUSEHANE_CHECKPOINT_FILE_INVALID_INPUT;
  }
  char *temporary = malloc(path_length + 5);
  if (!temporary) {
    free(payload);
    return MUSEHANE_CHECKPOINT_FILE_MODEL_FAILED;
  }
  memcpy(temporary, path, path_length);
  memcpy(temporary + path_length, ".tmp", 5);
  FILE *file = fopen(temporary, "wb");
  int written = file && fwrite(&header, sizeof header, 1, file) == 1 &&
    fwrite(payload, payload_size, 1, file) == 1 && fflush(file) == 0;
  if (file && fclose(file) != 0)
    written = 0;
  if (!written || rename(temporary, path) != 0) {
    remove(temporary);
    free(temporary);
    free(payload);
    return MUSEHANE_CHECKPOINT_FILE_IO_FAILED;
  }
  free(temporary);
  free(payload);
  *audit = (MusehaneCheckpointFileAudit) {
    .time = time,
    .iteration = iteration,
    .accepted_steps = accepted_steps,
    .payload_size = payload_size,
    .payload_checksum = checksum
  };
  return MUSEHANE_CHECKPOINT_FILE_OK;
}

MusehaneCheckpointFileStatus musehane_checkpoint_file_restore (
  const char *path, const FilmModelConfig *config, double expected_time,
  unsigned long expected_iteration, FilmModelState *state,
  unsigned long *accepted_steps, MusehaneCheckpointFileAudit *audit)
{
  if (!path || !path[0] || !config || !isfinite(expected_time) ||
      expected_time < 0. || !state || state->guard != 0U ||
      state->implementation != NULL || !accepted_steps || !audit)
    return MUSEHANE_CHECKPOINT_FILE_INVALID_INPUT;
  FILE *file = fopen(path, "rb");
  if (!file)
    return MUSEHANE_CHECKPOINT_FILE_IO_FAILED;
  MusehaneCheckpointFileHeader header;
  if (fread(&header, sizeof header, 1, file) != 1) {
    fclose(file);
    return MUSEHANE_CHECKPOINT_FILE_FORMAT_MISMATCH;
  }
  if (header.magic != MUSEHANE_CHECKPOINT_FILE_MAGIC ||
      header.schema != MUSEHANE_CHECKPOINT_FILE_SCHEMA ||
      header.restart_schema != FILM_MODEL_RESTART_SCHEMA_VERSION ||
      header.model_kind != FILM_MODEL_MUSEHANE_INVERSE ||
      header.reserved != 0U || header.payload_size == 0U ||
      header.payload_size > SIZE_MAX || header.iteration != expected_iteration ||
      header.accepted_steps > ULONG_MAX ||
      !time_matches(header.time, expected_time)) {
    fclose(file);
    return MUSEHANE_CHECKPOINT_FILE_FORMAT_MISMATCH;
  }
  const size_t payload_size = (size_t)header.payload_size;
  unsigned char *payload = malloc(payload_size);
  if (!payload) {
    fclose(file);
    return MUSEHANE_CHECKPOINT_FILE_MODEL_FAILED;
  }
  const int payload_read = fread(payload, payload_size, 1, file) == 1;
  const int no_trailing_bytes = payload_read && fgetc(file) == EOF;
  const int close_succeeded = fclose(file) == 0;
  const int complete = payload_read && no_trailing_bytes && close_succeeded;
  if (!complete) {
    free(payload);
    return MUSEHANE_CHECKPOINT_FILE_FORMAT_MISMATCH;
  }
  if (checksum_bytes(payload, payload_size) != header.payload_checksum) {
    free(payload);
    return MUSEHANE_CHECKPOINT_FILE_CHECKSUM_MISMATCH;
  }
  FilmRestartRecord record = {
    .schema_version = header.restart_schema,
    .model_kind = (FilmModelKind)header.model_kind,
    .size = payload_size,
    .capacity = payload_size,
    .data = payload
  };
  FilmModelState candidate = FILM_MODEL_STATE_INITIALIZER;
  const FilmModelStatus restored = film_model_restore(&candidate, config,
    &record);
  free(payload);
  if (restored != FILM_MODEL_OK)
    return MUSEHANE_CHECKPOINT_FILE_MODEL_FAILED;
  *state = candidate;
  *accepted_steps = (unsigned long)header.accepted_steps;
  *audit = (MusehaneCheckpointFileAudit) {
    .time = header.time,
    .iteration = (unsigned long)header.iteration,
    .accepted_steps = (unsigned long)header.accepted_steps,
    .payload_size = payload_size,
    .payload_checksum = header.payload_checksum
  };
  return MUSEHANE_CHECKPOINT_FILE_OK;
}
