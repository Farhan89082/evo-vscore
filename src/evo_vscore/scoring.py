"""Turn variants into reference/alternative sequences and score them.

For each variant we cut a window of reference DNA around it (default 8,192
bases), build the same window with the variant applied, and score both with a
backend. ``delta_ll = alt_ll - ref_ll``: more negative means the model finds the
variant less plausible than the reference.

Windows sit on a fixed grid (a new window every ``window/2`` bases) so variants
that are close together share one reference window and it is scored only once.
"""

from collections import Counter
from dataclasses import dataclass
from typing import Optional

from .cache import ScoreCache

_ACGT = set("ACGT")


@dataclass
class Prepared:
    status: str
    window_start: int = 0  # 0-based
    window_end: int = 0  # 0-based, exclusive
    ref_seq: str = ""
    alt_seq: str = ""


@dataclass
class ScoredVariant:
    chrom: str
    pos: int
    id: str
    ref: str
    alt: str
    variant_type: str
    status: str
    ref_ll: Optional[float] = None
    alt_ll: Optional[float] = None
    delta_ll: Optional[float] = None
    delta_mean_ll: Optional[float] = None
    window_start: Optional[int] = None  # 1-based, inclusive
    window_end: Optional[int] = None  # 1-based, inclusive
    backend: str = ""


def choose_window(pos0, ref_len, chrom_len, window):
    """Pick the grid window whose centre is closest to the variant.

    Returns (start, end), 0-based half-open, or None if the variant does not
    fit inside the window (very long deletions).
    """
    if chrom_len <= window:
        start, end = 0, chrom_len
    else:
        stride = window // 2
        centre = pos0 + ref_len / 2.0
        k = int(round((centre - window / 2.0) / stride))
        start = max(0, min(k * stride, chrom_len - window))
        end = start + window
    if pos0 < start or pos0 + ref_len > end:
        return None
    return start, end


def prepare(variant, fasta, window):
    if variant.skip_reason:
        return Prepared("skipped:" + variant.skip_reason)
    name = fasta.resolve(variant.chrom)
    if name is None:
        return Prepared("chrom_not_found")
    pos0 = variant.pos - 1
    chrom_len = fasta.length(name)
    if pos0 < 0 or pos0 + len(variant.ref) > chrom_len:
        return Prepared("pos_out_of_range")
    if fasta.fetch(name, pos0, pos0 + len(variant.ref)) != variant.ref:
        return Prepared("ref_mismatch")
    span = choose_window(pos0, len(variant.ref), chrom_len, window)
    if span is None:
        return Prepared("too_long")
    start, end = span
    ref_seq = fasta.fetch(name, start, end)
    if len(ref_seq) < 2:
        return Prepared("window_too_short")
    if "N" in variant.alt or not set(ref_seq) <= _ACGT:
        return Prepared("ambiguous_bases")
    offset = pos0 - start
    alt_seq = ref_seq[:offset] + variant.alt + ref_seq[offset + len(variant.ref):]
    return Prepared("ok", start, end, ref_seq, alt_seq)


def plan(variants, fasta, window, cache, cache_id):
    """Count what a run would need, without scoring anything (for --dry-run)."""
    statuses = Counter()
    keys = set()
    for variant in variants:
        prepared = prepare(variant, fasta, window)
        statuses[prepared.status] += 1
        if prepared.status == "ok":
            keys.add(cache.key(prepared.ref_seq))
            keys.add(cache.key(prepared.alt_seq))
    cached = sum(1 for key in keys if cache.has_key(cache_id, key))
    return {
        "total": sum(statuses.values()),
        "statuses": statuses,
        "scoreable": statuses["ok"],
        "unique_sequences": len(keys),
        "cached": cached,
        "calls_needed": len(keys) - cached,
    }


class Scorer:
    def __init__(self, fasta, backend, cache, window):
        self.fasta = fasta
        self.backend = backend
        self.cache = cache if cache is not None else ScoreCache(None)
        self.window = window
        self.cache_hits = 0
        self.backend_calls = 0

    def _ll(self, seq):
        value = self.cache.get(self.backend.cache_id, seq)
        if value is not None:
            self.cache_hits += 1
            return value
        value = self.backend.score(seq)
        self.backend_calls += 1
        self.cache.put(self.backend.cache_id, seq, value)
        return value

    def score(self, variant):
        prepared = prepare(variant, self.fasta, self.window)
        result = ScoredVariant(
            chrom=variant.chrom,
            pos=variant.pos,
            id=variant.id,
            ref=variant.ref,
            alt=variant.alt,
            variant_type=variant.variant_type,
            status=prepared.status,
            backend=self.backend.name,
        )
        if prepared.status != "ok":
            return result
        ref_ll = self._ll(prepared.ref_seq)
        alt_ll = self._ll(prepared.alt_seq)
        ref_n = max(len(prepared.ref_seq) - 1, 1)
        alt_n = max(len(prepared.alt_seq) - 1, 1)
        result.ref_ll = ref_ll
        result.alt_ll = alt_ll
        result.delta_ll = alt_ll - ref_ll
        result.delta_mean_ll = alt_ll / alt_n - ref_ll / ref_n
        result.window_start = prepared.window_start + 1
        result.window_end = prepared.window_end
        return result
