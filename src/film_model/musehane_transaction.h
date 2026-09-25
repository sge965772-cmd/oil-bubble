#ifndef FILM_MODEL_MUSEHANE_TRANSACTION_H
#define FILM_MODEL_MUSEHANE_TRANSACTION_H

#include <stdbool.h>

#include "film_model/musehane.h"

typedef enum {
  MUSEHANE_MODEL_TRANSACTION_OK = 0,
  MUSEHANE_MODEL_TRANSACTION_INVALID_STATE = 1,
  MUSEHANE_MODEL_TRANSACTION_ALLOCATION_FAILED = 2,
  MUSEHANE_MODEL_TRANSACTION_CHECKPOINT_FAILED = 3,
  MUSEHANE_MODEL_TRANSACTION_RESTORE_FAILED = 4
} MusehaneModelTransactionStatus;

typedef struct {
  FilmRestartRecord checkpoint;
  bool open;
} MusehaneModelTransaction;

#define MUSEHANE_MODEL_TRANSACTION_INITIALIZER {{0}, false}

MusehaneModelTransactionStatus musehane_model_transaction_begin (
  MusehaneModelTransaction *transaction, const FilmModelState *state);

MusehaneModelTransactionStatus musehane_model_transaction_commit (
  MusehaneModelTransaction *transaction);

MusehaneModelTransactionStatus musehane_model_transaction_rollback (
  MusehaneModelTransaction *transaction, FilmModelState *state,
  const FilmModelConfig *config);

void musehane_model_transaction_destroy (
  MusehaneModelTransaction *transaction);

#endif
