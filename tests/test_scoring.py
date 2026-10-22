import numpy as np
import pytest

from evo_vscore.backends import BackendError, MockBackend, sequence_logprob
from evo_vscore.cache import ScoreCache
from evo_vscore.fasta import FastaReference
from evo_vscore.scoring import Scorer, choose_window, plan, prepare
from evo_vscore.vcf import Variant


def _variant(seqs, chrom, pos, ref_len, alt):
    ref = seqs[chrom][pos - 1 : pos - 1 + ref_len]
    return Variant(chrom, pos, ".", ref, alt)


# -- window selection --------------------------------------------------------
def test_window_small_contig_is_whole_contig():
    assert choose_window(10, 1, 500, 1000) == (0, 500)


def test_window_grid_and_containment():
    window = 1000
    for pos0 in range(0, 10000, 137):
        start, end = choose_window(pos0, 1, 10000, window)
        assert end - start == window
        assert start <= pos0 < end
        # grid-aligned, except possibly the window pinned to the contig end
        assert start % (window // 2) == 0 or end == 10000


def test_window_nearby_variants_share_a_window():
    assert choose_window(5000, 1, 100000, 1000) == choose_window(5050, 1, 100000, 1000)


def test_window_too_long_variant_is_rejected():
    assert choose_window(5000, 1200, 100000, 1000) is None


# -- sequence construction ---------------------------------------------------
def test_snv_alt_sequence(genome):
    path, seqs = genome
    ref_base = seqs["chr1"][999]
    alt = "A" if ref_base != "A" else "C"
    with FastaReference(path) as fa:
        p = prepare(Variant("chr1", 1000, ".", ref_base, alt), fa, 512)
    assert p.status == "ok"
    assert len(p.ref_seq) == len(p.alt_seq) == 512
    diffs = [i for i, (a, b) in enumerate(zip(p.ref_seq, p.alt_seq)) if a != b]
    assert diffs == [999 - p.window_start]
    assert p.alt_seq[diffs[0]] == alt


def test_insertion_and_deletion_alt_sequences(genome):
    path, seqs = genome
    s = seqs["chr1"]
    with FastaReference(path) as fa:
        ins = prepare(Variant("chr1", 1500, ".", s[1499], s[1499] + "TTT"), fa, 512)
        dele = prepare(Variant("chr1", 1500, ".", s[1499:1503], s[1499]), fa, 512)
    assert len(ins.alt_seq) == len(ins.ref_seq) + 3
    assert len(dele.alt_seq) == len(dele.ref_seq) - 3
    off = 1499 - ins.window_start
    assert ins.alt_seq[off : off + 4] == s[1499] + "TTT"
    assert dele.alt_seq == dele.ref_seq[:off + 1] + dele.ref_seq[off + 4 :]


def test_variant_near_contig_ends(genome):
    path, seqs = genome
    s = seqs["chr1"]
    with FastaReference(path) as fa:
        first = prepare(Variant("chr1", 1, ".", s[0], "A" if s[0] != "A" else "C"), fa, 512)
        last = prepare(Variant("chr1", 3000, ".", s[2999], "A" if s[2999] != "A" else "C"), fa, 512)
    assert first.status == last.status == "ok"
    assert first.window_start == 0 and last.window_end == 3000


@pytest.mark.parametrize("chrom,pos,ref,alt,status", [
    ("chrX", 10, "A", "C", "chrom_not_found"),
    ("chr1", 99999, "A", "C", "pos_out_of_range"),
    ("chr1", 0, "A", "C", "pos_out_of_range"),
])
def test_prepare_bad_locations(genome, chrom, pos, ref, alt, status):
    path, _ = genome
    with FastaReference(path) as fa:
        assert prepare(Variant(chrom, pos, ".", ref, alt), fa, 512).status == status


def test_prepare_ref_mismatch_and_skips(genome):
    path, seqs = genome
    wrong = "A" if seqs["chr1"][49] != "A" else "C"
    with FastaReference(path) as fa:
        assert prepare(Variant("chr1", 50, ".", wrong, "T"), fa, 512).status == "ref_mismatch"
        v = Variant("chr1", 50, ".", seqs["chr1"][49], "<DEL>", skip_reason="symbolic_or_breakend_allele")
        assert prepare(v, fa, 512).status == "skipped:symbolic_or_breakend_allele"


def test_prepare_ambiguous_bases(tmp_path):
    fa_path = tmp_path / "n.fa"
    fa_path.write_text(">c\n" + "ACGT" * 5 + "N" + "ACGT" * 5 + "\n")
    with FastaReference(fa_path) as fa:
        assert prepare(Variant("c", 3, ".", "G", "A"), fa, 256).status == "ambiguous_bases"
        assert prepare(Variant("c", 3, ".", "G", "N"), fa, 256).status == "ambiguous_bases"


def test_too_long_deletion(genome):
    path, seqs = genome
    with FastaReference(path) as fa:
        big = seqs["chr1"][1000:1700]
        assert prepare(Variant("chr1", 1001, ".", big, big[0]), fa, 512).status == "too_long"


# -- scoring maths -----------------------------------------------------------
def test_sequence_logprob_uniform():
    seq = "ACGTACGT"
    logits = np.zeros((len(seq), 512))
    assert sequence_logprob(logits, seq) == pytest.approx((len(seq) - 1) * -np.log(512))


def test_sequence_logprob_picks_the_right_token():
    seq = "ACG"
    logits = np.full((3, 512), -50.0)
    logits[0, ord("C")] = 50.0  # after A the model is sure the next base is C
    logits[1, ord("G")] = 50.0  # after C it is sure of G
    assert sequence_logprob(logits, seq) == pytest.approx(0.0, abs=1e-9)
    logits[1, ord("G")] = -50.0
    logits[1, ord("T")] = 50.0  # now it expects T but sees G
    assert sequence_logprob(logits, seq) == pytest.approx(-100.0, abs=1e-6)


def test_sequence_logprob_shape_mismatch():
    with pytest.raises(BackendError):
        sequence_logprob(np.zeros((5, 512)), "ACGTACGT")


# -- end to end with the mock ----------------------------------------------------
def test_mock_backend_is_deterministic_and_sensitive():
    b = MockBackend()
    assert b.score("ACGTACGTAC") == b.score("ACGTACGTAC")
    assert b.score("ACGTACGTAC") != b.score("ACGTACGAAC")


class CountingBackend(MockBackend):
    cache_id = "counting:v1"

    def __init__(self):
        self.calls = 0

    def score(self, seq):
        self.calls += 1
        return super().score(seq)


def test_scorer_reuses_reference_windows_and_cache(genome):
    path, seqs = genome
    variants = [_variant(seqs, "chr1", p, 1, "A" if seqs["chr1"][p - 1] != "A" else "C") for p in (1000, 1010, 1020)]
    backend = CountingBackend()
    with FastaReference(path) as fa:
        scorer = Scorer(fa, backend, ScoreCache(None), 512)
        results = [scorer.score(v) for v in variants]
    assert all(r.status == "ok" for r in results)
    assert backend.calls == 1 + 3  # one shared reference window + three alts
    assert len({r.ref_ll for r in results}) == 1
    assert all(r.delta_ll == pytest.approx(r.alt_ll - r.ref_ll) for r in results)
    assert all(r.window_start == 1 or r.window_start > 1 for r in results)


def test_plan_counts_calls_and_cached(genome):
    path, seqs = genome
    variants = [_variant(seqs, "chr1", p, 1, "A" if seqs["chr1"][p - 1] != "A" else "C") for p in (1000, 1010)]
    cache = ScoreCache(None)
    with FastaReference(path) as fa:
        info = plan(variants, fa, 512, cache, "counting:v1")
        assert info["scoreable"] == 2 and info["unique_sequences"] == 3 and info["calls_needed"] == 3
        Scorer(fa, CountingBackend(), cache, 512).score(variants[0])
        info = plan(variants, fa, 512, cache, "counting:v1")
    assert info["cached"] == 2 and info["calls_needed"] == 1


def test_cache_persists_on_disk(tmp_path):
    path = tmp_path / "sub" / "c.sqlite"
    c = ScoreCache(path)
    c.put("b", "ACGT", -1.2345678901234567)
    c.close()
    c2 = ScoreCache(path)
    assert c2.get("b", "ACGT") == -1.2345678901234567
    assert c2.get("b", "ACGA") is None
    assert c2.get("other", "ACGT") is None
