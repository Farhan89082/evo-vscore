"""A fake backend for tests and demos. Its numbers are NOT biology."""

import zlib

from .base import Backend


class MockBackend(Backend):
    name = "mock"
    cache_id = "mock:v1"
    is_mock = True

    def score(self, seq):
        # Deterministic pseudo log-probability from a hash of the 7 bases ending
        # at each position. A changed base changes the hash of the 7 positions
        # around it, so variants get different (but meaningless) deltas.
        total = 0.0
        for i in range(1, len(seq)):
            context = seq[max(0, i - 6) : i + 1]
            u = zlib.crc32(context.encode("ascii")) / 4294967296.0
            total -= 0.9 + 0.8 * u
        return total
