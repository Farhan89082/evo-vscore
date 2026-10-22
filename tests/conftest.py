import random
from pathlib import Path

import pytest


def make_seq(n, seed=1):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture
def genome(tmp_path):
    """A small two-contig FASTA (60-column lines) and its sequences."""
    seqs = {"chr1": make_seq(3000, 1), "chr2": make_seq(1500, 2)}
    path = tmp_path / "genome.fa"
    with open(path, "w", newline="\n") as fh:
        for name, s in seqs.items():
            fh.write(">%s description\n" % name)
            for i in range(0, len(s), 60):
                fh.write(s[i : i + 60] + "\n")
    return path, seqs


def write_vcf(path, rows):
    with open(path, "w", newline="\n") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n")
        for row in rows:
            fh.write("\t".join(str(c) for c in row) + "\n")
    return Path(path)
