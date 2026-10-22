import csv
import gzip

import numpy as np
import pytest

from evo_vscore import cli, validation as v
from evo_vscore.fasta import FastaReference

HEADER = ["#AlleleID", "Type", "Name", "GeneID", "GeneSymbol", "ClinicalSignificance", "ClinSigSimple",
          "Assembly", "ChromosomeAccession", "Chromosome", "Start", "Stop", "ReviewStatus",
          "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"]
ONE_STAR = "criteria provided, single submitter"


def row(vtype="Deletion", name="NM_1(G):c.100del", gene="G", sig="Pathogenic", assembly="GRCh37", chrom="1",
        review=ONE_STAR, pos="100", ref="AC", alt="A"):
    return ["1", vtype, name, "1", gene, sig, "1", assembly, "NC_1", chrom, pos, pos, review, pos, ref, alt]


def write_clinvar(path, rows):
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(HEADER)
        w.writerows(rows)


@pytest.fixture
def clinvar_file(genome, tmp_path):
    _, seqs = genome
    rows = []
    for i in range(14):  # pathogenic coding-style deletions on chr1
        pos = 200 + i * 150
        s = seqs["chr1"]
        rows.append(row("Deletion", "NM_1(G):c.%ddel" % (i + 10), "G1", "Pathogenic", chrom="1",
                        pos=str(pos), ref=s[pos - 1:pos + 3], alt=s[pos - 1]))
    for i in range(14):  # benign intronic insertions on chr2
        pos = 100 + i * 80
        b = seqs["chr2"][pos - 1]
        rows.append(row("Insertion", "NM_2(H):c.%d+5dup" % (i + 10), "H1", "Likely benign", chrom="2",
                        pos=str(pos), ref=b, alt=b + "AC"))
    s1 = seqs["chr1"]
    rows += [
        row("single nucleotide variant", pos="900", ref=s1[899], alt="A" if s1[899] != "A" else "C"),
        row(sig="Uncertain significance", pos="2900", ref=s1[2899:2901], alt=s1[2899]),
        row(sig="Benign", review="no assertion criteria provided", pos="2800", ref=s1[2799:2801], alt=s1[2799]),
        row(assembly="GRCh38", pos="2700", ref=s1[2699:2701], alt=s1[2699]),
        row(chrom="MT", pos="50", ref="AC", alt="A"),
        row(pos="2600", ref=s1[2599:2605] + "A" * 100, alt=s1[2599]),            # 100-base deletion
        row(pos="2500", ref="na", alt="na"),
        row(pos="2400", ref=s1[2399:2401], alt=s1[2399], sig="Pathogenic"),     # conflicting pair below
        row(pos="2400", ref=s1[2399:2401], alt=s1[2399], sig="Benign"),
        row(pos="200", ref=s1[199:203], alt=s1[199]),                           # exact duplicate of a kept row
    ]
    path = tmp_path / "variant_summary.txt.gz"
    write_clinvar(path, rows)
    return path


@pytest.mark.parametrize("name,expected", [
    ("NM_007294.4(BRCA1):c.5266dup (p.Gln1756fs)", "coding_approx"),
    ("NM_1(G):c.123_125del", "coding_approx"),
    ("NM_1(G):c.4986+1_4986+2del", "noncoding_approx"),
    ("NM_1(G):c.123-5del", "noncoding_approx"),
    ("NM_1(G):c.-20del", "noncoding_approx"),
    ("NM_1(G):c.*5del", "noncoding_approx"),
    ("NR_1(G):n.123del", "noncoding_approx"),
    ("", "noncoding_approx"),
])
def test_approx_region(name, expected):
    assert v.approx_region(name) == expected


def test_prepare_clinvar_filters_and_reports(clinvar_file, tmp_path):
    info = v.prepare_clinvar(clinvar_file, tmp_path / "cv", chroms=["1", "2"], max_per_class=0)
    assert info["classes"] == {"BENIGN": 14, "PATH": 14}
    d = info["dropped"]
    assert d["not_clearly_pathogenic_or_benign"] == 1
    assert d["low_review_status"] == 1
    assert d["chromosome_not_selected"] == 1  # MT
    assert d["longer_than_max_indel"] == 1
    assert d["unusable_alleles"] == 1
    assert d["conflicting_duplicate"] == 1
    assert info["regions"][("PATH", "coding_approx")] == 14
    assert info["regions"][("BENIGN", "noncoding_approx")] == 14
    labels = list(csv.DictReader(open(tmp_path / "cv" / "clinvar_labels.csv", newline="")))
    vcf = [l for l in open(tmp_path / "cv" / "clinvar.vcf") if not l.startswith("#")]
    assert len(labels) == len(vcf) == 28 and labels[0]["id"] == "c0001"
    assert {l["gene"] for l in labels} == {"G1", "H1"}


def test_prepare_clinvar_caps_and_is_reproducible(clinvar_file, tmp_path):
    a = v.prepare_clinvar(clinvar_file, tmp_path / "a", chroms=["1", "2"], max_per_class=5, seed=3)
    b = v.prepare_clinvar(clinvar_file, tmp_path / "b", chroms=["1", "2"], max_per_class=5, seed=3)
    assert a["classes"] == {"BENIGN": 5, "PATH": 5} and a["available"] == {"BENIGN": 14, "PATH": 14}
    assert (tmp_path / "a" / "clinvar.vcf").read_text() == (tmp_path / "b" / "clinvar.vcf").read_text()


def test_prepare_clinvar_min_stars(clinvar_file, tmp_path):
    # every row in the test table is 1-star, so asking for 2 stars leaves nothing
    with pytest.raises(ValueError, match="No usable"):
        v.prepare_clinvar(clinvar_file, tmp_path / "s", chroms=["1", "2"], min_stars=2)


def test_prepare_clinvar_assembly_selection(clinvar_file, tmp_path):
    info = v.prepare_clinvar(clinvar_file, tmp_path / "g38", assembly="GRCh38", chroms=["1"], max_per_class=0)
    assert info["n"] == 1 and info["classes"] == {"PATH": 1}


def test_clinvar_missing_columns_has_a_clear_message(tmp_path):
    bad = tmp_path / "bad.txt"
    bad.write_text("#A\tB\nx\ty\n")
    with pytest.raises(ValueError, match="Columns found"):
        v.read_clinvar(bad)
    with pytest.raises(FileNotFoundError):
        v.read_clinvar(tmp_path / "nope.txt")


def test_clinvar_accepts_newer_column_names(tmp_path):
    header = [("GermlineClassification" if h == "ClinicalSignificance" else
               "GermlineReviewStatus" if h == "ReviewStatus" else h) for h in HEADER]
    path = tmp_path / "new.txt"
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerow(row(pos="10", ref="AC", alt="A"))
    records, _ = v.read_clinvar(path, chroms=["1"])
    assert len(records) == 1 and records[0]["class"] == "PATH"


# -- evaluation with indel features -----------------------------------------------------
def test_evaluate_reports_types_and_per_base_score(tmp_path):
    rng = np.random.default_rng(4)
    types = ["SNV"] * 40 + ["INS"] * 40 + ["DEL"] * 40
    classes = (["PATH"] * 20 + ["BENIGN"] * 20) * 3
    with open(tmp_path / "l.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "chrom", "pos", "ref", "alt", "class", "score", "consequence"])
        for i, c in enumerate(classes):
            w.writerow(["c%03d" % i, "1", i + 1, "A", "C", c, "", ""])
    with open(tmp_path / "r.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "status", "delta_ll", "delta_mean_ll", "variant_type"])
        for i, (t, c) in enumerate(zip(types, classes)):
            d = rng.normal() - 4 * (c == "PATH")
            w.writerow(["c%03d" % i, "ok", "%.4f" % d, "%.6f" % (d / 8000), t])
    res = v.evaluate(tmp_path / "r.csv", tmp_path / "l.csv", n_boot=20, positive="PATH", drop=())
    names = [m["name"] for m in res["measures"]]
    assert names[0] == "PATH vs BENIGN (all)"
    assert any("per-base" in n for n in names)
    assert {"by variant type: SNV", "by variant type: INS", "by variant type: DEL"} <= set(names)
    assert not any("excluded" in n for n in names)
    assert res["measures"][0]["auroc"] > 0.9


def test_evaluate_wrong_positive_class_is_a_clear_error(tmp_path):
    with open(tmp_path / "l.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "chrom", "pos", "ref", "alt", "class", "score", "consequence"])
        w.writerow(["a", "1", 1, "A", "C", "FUNC", "", ""])
    with open(tmp_path / "r.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "status", "delta_ll"])
        w.writerow(["a", "ok", "-1"])
    with pytest.raises(ValueError, match="No scored variants have class"):
        v.evaluate(tmp_path / "r.csv", tmp_path / "l.csv", positive="LOF")


def test_clinvar_pipeline_with_mock_backend(clinvar_file, genome, tmp_path, capsys):
    fa, _ = genome
    assert v.main(["prepare-clinvar", "--summary-file", str(clinvar_file), "--out-dir", str(tmp_path / "cv"),
                   "--chroms", "1,2", "--max-per-class", "0"]) == 0
    out = capsys.readouterr().out
    assert "28 indels" in out and "coding_approx" in out
    assert cli.main(["--vcf", str(tmp_path / "cv" / "clinvar.vcf"), "--ref", str(fa), "--backend", "mock",
                     "--window", "512", "--no-cache", "--out", str(tmp_path / "r")]) == 0
    capsys.readouterr()
    assert v.main(["evaluate", "--results", str(tmp_path / "r.csv"), "--labels",
                   str(tmp_path / "cv" / "clinvar_labels.csv"), "--positive", "PATH", "--drop", "",
                   "--bootstrap", "10", "--summary", str(tmp_path / "s.md")]) == 0
    text = (tmp_path / "s.md").read_text()
    assert "PATH vs BENIGN (all)" in text and "per-base" in text


# -- chromosome download and parsing ----------------------------------------------------------
def test_parse_chroms():
    assert v.parse_chroms("17") == ["17"]
    assert v.parse_chroms("1-3, x") == ["1", "2", "3", "X"]
    assert len(v.parse_chroms("all")) == 23
    for bad in ("foo", "3-1", "1-30", ""):
        with pytest.raises(ValueError):
            v.parse_chroms(bad)


def test_fetch_chromosomes_joins_files_cleanly(tmp_path):
    src = tmp_path / "remote"
    src.mkdir()
    with gzip.open(src / "chr1.fa.gz", "wt") as fh:
        fh.write(">chr1\nACGTAC\nGT\n")
    with gzip.open(src / "chr2.fa.gz", "wt") as fh:
        fh.write(">chr2\nTTGA")  # no trailing newline
    with gzip.open(src / "chrX.fa.gz", "wt") as fh:
        fh.write(">chrX\nCCCC\n")
    out = tmp_path / "data" / "genome.fa"
    v.fetch_chromosomes(["1", "2", "X"], out, src.as_uri(), progress=lambda m: None)
    assert out.read_text() == ">chr1\nACGTAC\nGT\n>chr2\nTTGA\n>chrX\nCCCC\n"
    with FastaReference(out) as fa:
        assert fa.names() == ["chr1", "chr2", "chrX"]
        assert fa.fetch("chr2", 0, 4) == "TTGA" and fa.fetch("chr1", 4, 8) == "ACGT"
    assert not (tmp_path / "data" / "genome.fa.part").exists()


def test_fetch_chromosomes_failure_is_a_clean_message(tmp_path, capsys):
    code = v.main(["fetch-ref", "--out", str(tmp_path / "g.fa"), "--chroms", "5", "--base-url",
                   (tmp_path / "nothing").as_uri()])
    assert code == 2 and "error:" in capsys.readouterr().err
