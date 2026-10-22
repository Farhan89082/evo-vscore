import gzip

import pytest

from evo_vscore.fasta import FastaError, FastaReference
from evo_vscore.vcf import iter_vcf

from conftest import write_vcf


def test_fasta_fetch_matches_sequence(genome):
    path, seqs = genome
    with FastaReference(path) as fa:
        assert fa.names() == ["chr1", "chr2"]
        assert fa.length("chr1") == 3000
        for start, end in [(0, 10), (55, 130), (59, 61), (2990, 3000), (0, 3000)]:
            assert fa.fetch("chr1", start, end) == seqs["chr1"][start:end]
        assert fa.fetch("chr2", 1499, 1500) == seqs["chr2"][1499:1500]
        assert fa.fetch("chr1", 2990, 99999) == seqs["chr1"][2990:]
        assert fa.fetch("chr1", 5, 5) == ""


def test_fasta_uses_fai_when_present(genome):
    path, seqs = genome
    # a samtools-style index: name, length, offset, linebases, linewidth
    offset = len(">chr1 description\n")
    path_fai = str(path) + ".fai"
    with open(path_fai, "w") as fh:
        fh.write("chr1\t3000\t%d\t60\t61\n" % offset)
        off2 = offset + 3000 + 50 + len(">chr2 description\n")
        fh.write("chr2\t1500\t%d\t60\t61\n" % off2)
    with FastaReference(path) as fa:
        assert fa.fetch("chr1", 100, 200) == seqs["chr1"][100:200]
        assert fa.fetch("chr2", 100, 200) == seqs["chr2"][100:200]


def test_fasta_crlf_and_lowercase(tmp_path):
    path = tmp_path / "w.fa"
    path.write_bytes(b">s\r\nacgtac\r\ngtacgt\r\nac\r\n")
    with FastaReference(path) as fa:
        assert fa.length("s") == 14
        assert fa.fetch("s", 4, 10) == "ACGTAC"
        assert fa.fetch("s", 0, 14) == "ACGTACGTACGTAC"


def test_fasta_chr_prefix_resolution(genome):
    path, _ = genome
    with FastaReference(path) as fa:
        assert fa.resolve("chr1") == "chr1"
        assert fa.resolve("1") == "chr1"
        assert fa.resolve("chrZ") is None


def test_fasta_irregular_lines_rejected(tmp_path):
    path = tmp_path / "bad.fa"
    path.write_text(">s\nACGTAC\nACG\nACGTAC\n")
    with pytest.raises(FastaError):
        FastaReference(path)


def test_fasta_gz_rejected(tmp_path):
    path = tmp_path / "g.fa.gz"
    path.write_bytes(b"")
    with pytest.raises(FastaError):
        FastaReference(path)


def test_vcf_types_and_multiallelic(tmp_path):
    vcf = write_vcf(tmp_path / "a.vcf", [
        ("chr1", 10, "a", "A", "C"),
        ("chr1", 20, "b", "A", "AGT,T"),
        ("chr1", 30, "c", "ACG", "A"),
        ("chr1", 40, "d", "AC", "GT"),
    ])
    out = list(iter_vcf(vcf))
    assert [v.variant_type for v in out] == ["SNV", "INS", "SNV", "DEL", "MNV"]
    assert [v.alt for v in out] == ["C", "AGT", "T", "A", "GT"]
    assert all(v.skip_reason is None for v in out)


def test_vcf_unsupported_alleles_are_flagged_not_dropped(tmp_path):
    vcf = write_vcf(tmp_path / "b.vcf", [
        ("chr1", 10, ".", "A", "<DEL>"),
        ("chr1", 11, ".", "A", "*"),
        ("chr1", 12, ".", "A", "A[chr2:5["),
        ("chr1", 13, ".", "A", "."),
        ("chr1", 14, ".", "A", "A"),
        ("chr1", 15, ".", "A", "ACX"),
    ])
    out = list(iter_vcf(vcf))
    assert len(out) == 6
    assert all(v.skip_reason for v in out)
    assert all(v.variant_type == "OTHER" for v in out)


def test_vcf_lowercase_and_gzip(tmp_path):
    path = tmp_path / "c.vcf.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("#CHROM\tPOS\tID\tREF\tALT\nchr1\t5\t.\ta\tg\n")
    (v,) = list(iter_vcf(path))
    assert (v.ref, v.alt, v.pos) == ("A", "G", 5)


def test_vcf_bad_lines_raise(tmp_path):
    p1 = tmp_path / "d.vcf"
    p1.write_text("chr1\tfoo\t.\tA\tC\n")
    with pytest.raises(ValueError):
        list(iter_vcf(p1))
    p2 = tmp_path / "e.vcf"
    p2.write_text("chr1\t5\t.\n")
    with pytest.raises(ValueError):
        list(iter_vcf(p2))
