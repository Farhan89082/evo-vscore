"""Backend interface and the shared scoring maths."""

import abc

import numpy as np


class BackendError(RuntimeError):
    """A backend could not produce a score (bad key, network, bad response...)."""


class Backend(abc.ABC):
    #: short name used on the command line and in the output
    name = "base"
    #: identifies the model + settings; results are cached per cache_id
    cache_id = "base"
    #: True for backends whose numbers are not real model predictions
    is_mock = False

    @abc.abstractmethod
    def score(self, seq):
        """Return the summed log-probability of ``seq[1:]`` given the preceding bases."""

    def close(self):
        pass


def sequence_logprob(logits, seq):
    """Sum of log P(base i+1 | bases up to i) from next-token logits.

    ``logits`` has shape (len(seq), vocab). Row i holds the scores for the token
    that follows base i, so the last row is unused. Evo 2 is character-level:
    the token id of a base is its ASCII code (A=65, C=67, G=71, T=84).
    """
    logits = np.asarray(logits, dtype=np.float64)
    if logits.ndim != 2 or logits.shape[0] != len(seq):
        raise BackendError(
            "Model returned scores of shape %s for a sequence of %d bases"
            % (tuple(logits.shape), len(seq))
        )
    if len(seq) < 2:
        return 0.0
    ids = np.frombuffer(seq[1:].encode("ascii"), dtype=np.uint8).astype(np.int64)
    if ids.max() >= logits.shape[1]:
        raise BackendError("Sequence contains a character outside the model vocabulary")
    rows = logits[:-1]
    peak = rows.max(axis=1, keepdims=True)
    log_norm = peak[:, 0] + np.log(np.exp(rows - peak).sum(axis=1))
    chosen = rows[np.arange(len(ids)), ids]
    return float((chosen - log_norm).sum())
