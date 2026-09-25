#include "film_model/musehane_transaction.h"

#include <stdlib.h>

MusehaneModelTransactionStatus musehane_model_transaction_begin (
  MusehaneModelTransaction *transaction, const FilmModelState *state)
{
  if (!transaction || transaction->open || !state)
    return MUSEHANE_MODEL_TRANSACTION_INVALID_STATE;
  const size_t required = film_model_checkpoint_size(state);
  if (required == 0)
    return MUSEHANE_MODEL_TRANSACTION_CHECKPOINT_FAILED;
  if (transaction->checkpoint.capacity < required) {
    unsigned char *replacement = realloc(transaction->checkpoint.data,
      required);
    if (!replacement)
      return MUSEHANE_MODEL_TRANSACTION_ALLOCATION_FAILED;
    transaction->checkpoint.data = replacement;
    transaction->checkpoint.capacity = required;
  }
  if (film_model_checkpoint(state, &transaction->checkpoint) !=
      FILM_MODEL_OK)
    return MUSEHANE_MODEL_TRANSACTION_CHECKPOINT_FAILED;
  transaction->open = true;
  return MUSEHANE_MODEL_TRANSACTION_OK;
}

MusehaneModelTransactionStatus musehane_model_transaction_commit (
  MusehaneModelTransaction *transaction)
{
  if (!transaction || !transaction->open)
    return MUSEHANE_MODEL_TRANSACTION_INVALID_STATE;
  transaction->open = false;
  return MUSEHANE_MODEL_TRANSACTION_OK;
}

MusehaneModelTransactionStatus musehane_model_transaction_rollback (
  MusehaneModelTransaction *transaction, FilmModelState *state,
  const FilmModelConfig *config)
{
  if (!transaction || !transaction->open || !state || !config)
    return MUSEHANE_MODEL_TRANSACTION_INVALID_STATE;
  FilmModelState candidate = FILM_MODEL_STATE_INITIALIZER;
  if (film_model_restore(&candidate, config, &transaction->checkpoint) !=
      FILM_MODEL_OK)
    return MUSEHANE_MODEL_TRANSACTION_RESTORE_FAILED;
  film_model_destroy(state);
  *state = candidate;
  transaction->open = false;
  return MUSEHANE_MODEL_TRANSACTION_OK;
}

void musehane_model_transaction_destroy (
  MusehaneModelTransaction *transaction)
{
  if (!transaction)
    return;
  free(transaction->checkpoint.data);
  *transaction = (MusehaneModelTransaction)
    MUSEHANE_MODEL_TRANSACTION_INITIALIZER;
}
