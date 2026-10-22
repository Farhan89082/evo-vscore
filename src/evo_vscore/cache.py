"""Tiny on-disk cache of sequence scores (SQLite, standard library only).

Every score is stored under (backend id, hash of the sequence), so an
interrupted run resumes without repeating calls, and a reference window shared
by several variants is only ever scored once.
"""

import hashlib
import sqlite3
from pathlib import Path


class ScoreCache:
    def __init__(self, path=None):
        """``path=None`` keeps the cache in memory only (used by --no-cache)."""
        if path is None:
            target = ":memory:"
        else:
            Path(path).expanduser().parent.mkdir(parents=True, exist_ok=True)
            target = str(Path(path).expanduser())
        self._db = sqlite3.connect(target)
        self._db.execute(
            "CREATE TABLE IF NOT EXISTS scores ("
            "backend TEXT NOT NULL, seqhash TEXT NOT NULL, ll REAL NOT NULL, "
            "PRIMARY KEY (backend, seqhash))"
        )
        self._db.commit()

    @staticmethod
    def key(seq):
        return hashlib.sha256(seq.encode("ascii")).hexdigest()

    def get(self, backend_id, seq):
        row = self._db.execute(
            "SELECT ll FROM scores WHERE backend=? AND seqhash=?",
            (backend_id, self.key(seq)),
        ).fetchone()
        return None if row is None else row[0]

    def has_key(self, backend_id, seqhash):
        row = self._db.execute(
            "SELECT 1 FROM scores WHERE backend=? AND seqhash=?", (backend_id, seqhash)
        ).fetchone()
        return row is not None

    def put(self, backend_id, seq, ll):
        self._db.execute(
            "INSERT OR REPLACE INTO scores (backend, seqhash, ll) VALUES (?, ?, ?)",
            (backend_id, self.key(seq), float(ll)),
        )
        self._db.commit()

    def close(self):
        self._db.close()
