import csv
import gzip
import itertools
import random
from pathlib import Path

import numpy as np
import pytest

from evo_vscore import cli, validation as v

openpyxl = pytest.importorskip("openpyxl")


# -- statistics ----------------------------------------------------------------
def brute_force_auroc(y, s):
    pos = [a for a, b in zip(s, y) if b]
    neg = [a for a, b in zip(s, y) if not b]
    wins = sum(1.0 if p > n else 0.5 if p == n else 0.0 for p, n in itertools.product(pos, neg))
    return wins / (len(pos) * len(neg))


def test_auroc_known_cases():
    assert v.auroc([1, 1, 0, 0], [0.9, 0.8, 0.2, 0.1]) == 1.0
    assert v.auroc([1, 1, 0, 0], [0.1, 0.2, 0.8, 0.9]) == 0.0
    assert v.auroc([1, 0, 1, 0], [0.5, 0.5, 0.5, 0.5]) == 0.5


def test_auroc_matches_brute_force_with_ties():
    rng = random.Random(3)
    for _ in range(20):
        y = [rng.random() < 0.4 for _ in range(60)]
        s = [round(rng.random(), 1) for _ in range(60)]  # many ties
        if any(y) and not all(y):
            assert v.auroc(y, s) == pytest.approx(brute_force_auroc(y, s))


def test_auroc_needs_both_classes():
    with pytest.raises(ValueError):
        v.auroc([1, 1, 1], [0.1, 0.2, 0.3])


def test_bootstrap_ci_brackets_estimate():
    rng = np.random.default_rng(0)
    y = rng.random(300) < 0.3
    s = rng.normal(size=300) + 1.5 * y
    low, high = v.bootstrap_ci(y, s, n_boot=200)
    assert low < v.auroc(y, s) < high


# -- spreadsheet reading -----------------------------------------------------------
HEADER = ["chromosome", "position (hg19)", "reference", "alt", "function.score.mean", "func.class", "consequence"]


def make_xlsx(path, rows, junk_rows=2):
    wb = openpyxl.Workbook()
    ws = wb.active
    for i in range(junk_rows):
        ws.append(["Supplementary Table 1 (title row %d)" % i])
    ws.append(HEADER)
    for r in rows:
        ws.append(r)
    wb.save(path)


def test_read_findlay_filters_rows(tmp_path):
    x = tmp_path / "t.xlsx"
    make_xlsx(x, [
        ["chr17", 100, "A", "G", -1.5, "LOF", "Missense"],
        ["chr17", 101, "C", "T", 0.1, "FUNC", "Synonymous"],
        ["chr17", 102, "G", "A", -0.4, "INT", "Missense"],
        ["chr17", 103, "G", "GA", 0.0, "FUNC", "Insertion"],   # not an SNV
        ["chr17", 104, "T", "C", 0.0, "unknown", "Missense"],  # unlabelled
        [None, None, None, None, None, None, None],            # incomplete
    ])
    rows, skipped = v.read_findlay_xlsx(x)
    assert [(r["pos"], r["class"]) for r in rows] == [(100, "LOF"), (101, "FUNC"), (102, "INT")]
    assert skipped == {"not_a_snv": 1, "unlabelled": 1, "incomplete": 1}


def test_read_findlay_wrong_file(tmp_path):
    wb = openpyxl.Workbook()
    wb.active.append(["a", "b"])
    wb.save(tmp_path / "w.xlsx")
    with pytest.raises(ValueError, match="header"):
        v.read_findlay_xlsx(tmp_path / "w.xlsx")
    with pytest.raises(FileNotFoundError):
        v.read_findlay_xlsx(tmp_path / "nope.xlsx")


def test_prepare_writes_vcf_and_labels_with_sampling(tmp_path):
    x = tmp_path / "t.xlsx"
    rows = [["17", 1000 + i, "ACGT"[i % 4], "ACGT"[(i + 1) % 4], 0.0, ["LOF", "FUNC", "INT"][i % 3], ""] for i in range(30)]
    make_xlsx(x, rows)
    info = v.prepare(x, tmp_path / "out", sample=12, seed=1)
    assert info["n"] == 12
    vcf = [l for l in (tmp_path / "out" / "brca1.vcf").read_text().splitlines() if not l.startswith("#")]
    labels = list(csv.DictReader(open(tmp_path / "out" / "brca1_labels.csv", newline="")))
    assert len(vcf) == len(labels) == 12
    assert [l.split("\t")[2] for l in vcf] == [r["id"] for r in labels]
    # sampling is reproducible
    again = v.prepare(x, tmp_path / "out2", sample=12, seed=1)
    assert (tmp_path / "out2" / "brca1_labels.csv").read_text() == (tmp_path / "out" / "brca1_labels.csv").read_text()
    assert again["classes"] == info["classes"]


# -- evaluation ---------------------------------------------------------------------
def write_results(path, labels, deltas, status="ok"):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "status", "delta_ll"])
        for lab, d in zip(labels, deltas):
            w.writerow([lab["id"], status, "%.4f" % d])


def write_labels(path, classes, consequences=None):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "chrom", "pos", "ref", "alt", "class", "score", "consequence"])
        for i, c in enumerate(classes):
            w.writerow(["b%04d" % i, "17", 100 + i, "A", "C", c, "", (consequences or [""] * len(classes))[i]])
    return [{"id": "b%04d" % i} for i in range(len(classes))]


def test_evaluate_detects_planted_signal_and_chance(tmp_path):
    rng = np.random.default_rng(1)
    classes = ["LOF"] * 60 + ["FUNC"] * 100 + ["INT"] * 40
    labs = write_labels(tmp_path / "l.csv", classes)
    lof = np.array([c == "LOF" for c in classes])
    good = rng.normal(size=200) - 3.0 * lof          # LOF get much more negative deltas
    write_results(tmp_path / "good.csv", labs, good)
    res = v.evaluate(tmp_path / "good.csv", tmp_path / "l.csv", n_boot=50)
    assert res["measures"][0]["auroc"] > 0.95
    assert res["scored"] == 200 and res["class_counts"] == {"FUNC": 100, "INT": 40, "LOF": 60}
    assert res["mean_delta"]["LOF"] < res["mean_delta"]["FUNC"]
    assert any("FUNC only" in m["name"] for m in res["measures"])

    write_results(tmp_path / "noise.csv", labs, rng.normal(size=200))
    res = v.evaluate(tmp_path / "noise.csv", tmp_path / "l.csv", n_boot=50)
    assert 0.35 < res["measures"][0]["auroc"] < 0.65


def test_evaluate_ignores_unscored_and_reports_statuses(tmp_path):
    labs = write_labels(tmp_path / "l.csv", ["LOF", "FUNC", "FUNC", "LOF"])
    with open(tmp_path / "r.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "status", "delta_ll"])
        w.writerows([["b0000", "ok", "-5"], ["b0001", "ok", "0.2"], ["b0002", "ref_mismatch", ""], ["b0003", "ok", "-4"]])
    res = v.evaluate(tmp_path / "r.csv", tmp_path / "l.csv", n_boot=10)
    assert res["scored"] == 3 and res["statuses"]["ref_mismatch"] == 1
    assert res["measures"][0]["auroc"] == 1.0


def test_evaluate_all_unscored_is_an_error(tmp_path):
    labs = write_labels(tmp_path / "l.csv", ["LOF", "FUNC"])
    write_results(tmp_path / "r.csv", labs, [0, 0], status="ref_mismatch")
    with pytest.raises(ValueError, match="None of the variants"):
        v.evaluate(tmp_path / "r.csv", tmp_path / "l.csv")


def test_evaluate_by_consequence_needs_enough_examples(tmp_path):
    rng = np.random.default_rng(2)
    classes = (["LOF"] * 15 + ["FUNC"] * 15) * 2 + ["LOF"] * 2 + ["FUNC"] * 2
    cons = ["Missense"] * 30 + ["Intronic"] * 30 + ["Rare"] * 4
    labs = write_labels(tmp_path / "l.csv", classes, cons)
    lof = np.array([c == "LOF" for c in classes])
    write_results(tmp_path / "r.csv", labs, rng.normal(size=len(classes)) - 2 * lof)
    res = v.evaluate(tmp_path / "r.csv", tmp_path / "l.csv", n_boot=20)
    names = [m["name"] for m in res["measures"]]
    assert "by consequence: Missense" in names and "by consequence: Intronic" in names
    assert "by consequence: Rare" not in names
    assert "AUROC" in v.format_summary(res)


# -- reference download and the command line ----------------------------------------------------
def test_fetch_reference_from_file_url(tmp_path):
    src = tmp_path / "chr17.fa.gz"
    with gzip.open(src, "wt") as fh:
        fh.write(">chr17\nACGTACGT\nACG\n")
    out = tmp_path / "data" / "chr17.fa"
    v.fetch_reference(src.as_uri(), out)
    assert out.read_text() == ">chr17\nACGTACGT\nACG\n"
    assert not (tmp_path / "data" / "chr17.fa.gz.part").exists()


def test_fetch_reference_failure_is_a_clean_message(tmp_path, capsys):
    code = v.main(["fetch-ref", "--out", str(tmp_path / "x.fa"), "--url", (tmp_path / "missing.gz").as_uri()])
    err = capsys.readouterr().err
    assert code == 2 and "error:" in err and "7-Zip" in err


def test_full_pipeline_with_mock_backend(tmp_path, capsys):
    """spreadsheet -> VCF -> evo-vscore (mock) -> evaluate, all on a synthetic chr17."""
    rng = random.Random(5)
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    fa = tmp_path / "chr17.fa"
    fa.write_text(">chr17\n" + "\n".join(seq[i:i + 50] for i in range(0, len(seq), 50)) + "\n")
    rows = []
    for i in range(60):
        pos = 500 + i * 80
        ref = seq[pos - 1]
        alt = "A" if ref != "A" else "C"
        rows.append(["chr17", pos, ref, alt, 0.0, ["LOF", "FUNC", "INT"][i % 3], "Missense"])
    x = tmp_path / "s.xlsx"
    make_xlsx(x, rows)
    assert v.main(["prepare", "--xlsx", str(x), "--out-dir", str(tmp_path / "d")]) == 0
    code = cli.main(["--vcf", str(tmp_path / "d" / "brca1.vcf"), "--ref", str(fa), "--backend", "mock",
                     "--window", "1024", "--no-cache", "--out", str(tmp_path / "r")])
    assert code == 0
    capsys.readouterr()
    summary = tmp_path / "sum.md"
    assert v.main(["evaluate", "--results", str(tmp_path / "r.csv"), "--labels",
                   str(tmp_path / "d" / "brca1_labels.csv"), "--summary", str(summary), "--bootstrap", "20"]) == 0
    text = summary.read_text()
    assert "Variants scored: 60" in text and "LOF vs FUNC+INT" in text
