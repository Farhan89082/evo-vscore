import os

from evo_vscore import fasta as fasta_module
from evo_vscore.fasta import FastaReference


def test_index_is_saved_and_reused(genome, monkeypatch):
    path, seqs = genome
    fai = str(path) + ".fai"
    assert not os.path.exists(fai)
    with FastaReference(path) as fa:
        assert fa.fetch("chr1", 100, 200) == seqs["chr1"][100:200]
    assert os.path.exists(fai)

    def boom(_):
        raise AssertionError("index should have been read from the .fai file")

    monkeypatch.setattr(fasta_module, "_build_index", boom)
    with FastaReference(path) as fa:  # second open must not rescan the FASTA
        assert fa.fetch("chr2", 10, 60) == seqs["chr2"][10:60]


def test_saved_index_matches_samtools_format(genome):
    path, _ = genome
    FastaReference(path).close()
    lines = open(str(path) + ".fai").read().splitlines()
    assert lines[0].split("\t") == ["chr1", "3000", str(len(">chr1 description\n")), "60", "61"]
    assert len(lines) == 2


def test_stale_index_is_rebuilt(genome):
    path, seqs = genome
    FastaReference(path).close()
    fai = str(path) + ".fai"
    with open(fai, "w") as fh:  # a wrong index that is OLDER than the FASTA
        fh.write("chr1\t5\t0\t60\t61\n")
    old = os.path.getmtime(path) - 100
    os.utime(fai, (old, old))
    with FastaReference(path) as fa:
        assert fa.length("chr1") == 3000
        assert fa.fetch("chr1", 0, 10) == seqs["chr1"][:10]
