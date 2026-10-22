"""Small, dependency-free VCF reader.

Reads the first five columns (CHROM, POS, ID, REF, ALT). Multi-allelic records
are split into one :class:`Variant` per ALT allele. Records this tool cannot
score (symbolic alleles, breakends, odd characters) are still returned, with
``skip_reason`` set, so they show up in the output instead of vanishing.
"""

import gzip
from dataclasses import dataclass
from typing import Iterator, Optional

_BASES = set("ACGTN")


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int  # 1-based, as in the VCF
    id: str
    ref: str
    alt: str
    skip_reason: Optional[str] = None

    @property
    def variant_type(self):
        if self.skip_reason:
            return "OTHER"
        if len(self.ref) == 1 and len(self.alt) == 1:
            return "SNV"
        if len(self.ref) == len(self.alt):
            return "MNV"
        return "INS" if len(self.alt) > len(self.ref) else "DEL"


def _open_text(path):
    if str(path).lower().endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline=None)
    return open(path, "r", encoding="utf-8", errors="replace", newline=None)


def _skip_reason(ref, alt):
    if not ref or not set(ref) <= _BASES:
        return "unsupported_ref_allele"
    if alt in (".", "*", "") or alt.startswith("<") or "[" in alt or "]" in alt:
        return "symbolic_or_breakend_allele"
    if not set(alt) <= _BASES:
        return "unsupported_alt_allele"
    if ref == alt:
        return "ref_equals_alt"
    return None


def iter_vcf(path) -> Iterator[Variant]:
    with _open_text(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                cols = line.split()  # tolerate space-separated files
            if len(cols) < 5:
                raise ValueError("VCF line %d has fewer than 5 columns" % line_no)
            chrom, pos_text, vid, ref, alts = cols[:5]
            try:
                pos = int(pos_text)
            except ValueError:
                raise ValueError("VCF line %d has a non-numeric POS: %r" % (line_no, pos_text))
            ref = ref.upper()
            for alt in alts.split(","):
                alt = alt.upper()
                yield Variant(chrom, pos, vid, ref, alt, _skip_reason(ref, alt))
