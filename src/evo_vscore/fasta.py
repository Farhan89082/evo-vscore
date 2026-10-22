"""Minimal random-access reader for an uncompressed, indexed-style FASTA file.

Uses a samtools-style ``.fai`` index if one sits next to the FASTA; otherwise
builds the same index in memory by scanning the file once. Only sequences with
regular line lengths (the normal case) are supported.
"""

import os


class FastaError(ValueError):
    pass


_IRREGULAR = (
    "FASTA sequence '%s' has irregular line lengths; "
    "re-wrap it (for example with 'samtools faidx') and try again"
)


def _parse_fai(fai_path):
    index = {}
    with open(fai_path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 5:
                continue
            name, length, offset, linebases, linewidth = parts[:5]
            index[name] = (int(length), int(offset), int(linebases), int(linewidth))
    return index


def _build_index(path):
    index = {}
    name = None
    length = offset = 0
    linebases = linewidth = None
    last_short = False
    pos = 0

    def close_record():
        if name is not None:
            index[name] = (length, offset, linebases or 0, linewidth or 0)

    with open(path, "rb") as fh:
        for raw in fh:
            line_len = len(raw)
            if raw.startswith(b">"):
                close_record()
                header = raw[1:].split()
                name = header[0].decode("ascii", "replace") if header else ""
                length = 0
                offset = pos + line_len
                linebases = linewidth = None
                last_short = False
            else:
                if name is None:
                    if raw.strip():
                        raise FastaError("FASTA file has sequence before the first '>' header")
                    pos += line_len
                    continue
                bases = len(raw.rstrip(b"\r\n"))
                if linebases is None:
                    if bases == 0:
                        offset = pos + line_len
                        pos += line_len
                        continue
                    linebases, linewidth = bases, line_len
                elif last_short:
                    # An earlier line was shorter than the rest, so only blank
                    # trailing lines are allowed from here on.
                    if bases != 0:
                        raise FastaError(_IRREGULAR % name)
                elif bases > linebases:
                    raise FastaError(_IRREGULAR % name)
                elif bases < linebases:
                    last_short = True
                length += bases
            pos += line_len
    close_record()
    return index


class FastaReference:
    """Fetch 0-based, half-open slices of a FASTA file as upper-case strings."""

    def __init__(self, path):
        path = os.fspath(path)
        if path.lower().endswith((".gz", ".bgz")):
            raise FastaError(
                "Compressed FASTA is not supported in v0.1. Decompress it first "
                "(for example 'gunzip genome.fa.gz')."
            )
        if not os.path.exists(path):
            raise FileNotFoundError("Reference FASTA not found: %s" % path)
        self.path = path
        fai = path + ".fai"
        self._index = _parse_fai(fai) if os.path.exists(fai) else _build_index(path)
        if not self._index:
            raise FastaError("No sequences found in %s" % path)
        self._fh = open(path, "rb")

    # -- context manager / cleanup ------------------------------------------
    def close(self):
        if self._fh is not None:
            self._fh.close()
            self._fh = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    # -- lookups ------------------------------------------------------------
    def names(self):
        return list(self._index)

    def resolve(self, chrom):
        """Map a VCF contig name to the FASTA's name (tolerates a 'chr' prefix)."""
        if chrom in self._index:
            return chrom
        candidates = []
        if chrom.lower().startswith("chr"):
            candidates.append(chrom[3:])
        else:
            candidates.append("chr" + chrom)
        if chrom in ("MT", "M"):
            candidates += ["chrM", "chrMT", "MT", "M"]
        if chrom in ("chrM", "chrMT"):
            candidates += ["MT", "M"]
        for c in candidates:
            if c in self._index:
                return c
        return None

    def length(self, name):
        return self._index[name][0]

    def fetch(self, name, start, end):
        length, offset, linebases, linewidth = self._index[name]
        start = max(0, start)
        end = min(end, length)
        if start >= end or linebases == 0:
            return ""
        byte_start = offset + (start // linebases) * linewidth + (start % linebases)
        last = end - 1
        byte_end = offset + (last // linebases) * linewidth + (last % linebases) + 1
        self._fh.seek(byte_start)
        data = self._fh.read(byte_end - byte_start)
        data = data.replace(b"\n", b"").replace(b"\r", b"")
        return data.decode("ascii").upper()
