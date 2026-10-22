import csv
from pathlib import Path

import pytest

from evo_vscore import cli
from evo_vscore.backends import BackendError

from conftest import write_vcf

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


def run(*args):
    return cli.main([str(a) for a in args])


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def test_example_files_end_to_end_with_mock(tmp_path):
    out = tmp_path / "res"
    code = run("--vcf", EXAMPLES / "tiny.vcf", "--ref", EXAMPLES / "tiny.fa",
               "--backend", "mock", "--out", out, "--no-cache")
    assert code == 0
    rows = read_rows(str(out) + ".csv")
    by_id = {r["id"]: r for r in rows}
    assert len(rows) == 7
    assert [by_id[i]["status"] for i in ("snv_1", "snv_2", "ins_1", "del_1", "snv_3")] == ["ok"] * 5
    assert by_id["bad_ref"]["status"] == "ref_mismatch"
    assert by_id["sv_1"]["status"].startswith("skipped:")
    assert by_id["sv_1"]["variant_type"] == "OTHER"
    assert by_id["snv_1"]["delta_ll"] != ""
    html = Path(str(out) + ".html").read_text(encoding="utf-8")
    assert "Mock backend" in html and "chr1:500" in html


def test_resume_uses_cache(tmp_path, capsys):
    cache = tmp_path / "c.sqlite"
    args = ["--vcf", EXAMPLES / "tiny.vcf", "--ref", EXAMPLES / "tiny.fa", "--backend", "mock",
            "--out", tmp_path / "r", "--cache", cache]
    assert run(*args) == 0
    capsys.readouterr()
    assert run(*args) == 0
    assert "Backend calls: 0" in capsys.readouterr().err


def test_dry_run_scores_nothing_and_needs_no_key(tmp_path, capsys, monkeypatch):
    monkeypatch.delenv("NVIDIA_API_KEY", raising=False)
    monkeypatch.delenv("NVCF_KEY", raising=False)
    out = tmp_path / "never"
    code = run("--vcf", EXAMPLES / "tiny.vcf", "--ref", EXAMPLES / "tiny.fa", "--backend", "nvidia",
               "--dry-run", "--no-cache", "--out", out)
    err = capsys.readouterr().err
    assert code == 0
    assert "Backend calls needed:     7" in err and "minutes" in err
    assert not Path(str(out) + ".csv").exists()


def test_missing_key_is_a_clean_error(tmp_path, capsys, monkeypatch):
    monkeypatch.delenv("NVIDIA_API_KEY", raising=False)
    monkeypatch.delenv("NVCF_KEY", raising=False)
    code = run("--vcf", EXAMPLES / "tiny.vcf", "--ref", EXAMPLES / "tiny.fa", "--backend", "nvidia",
               "--no-cache", "--out", tmp_path / "x")
    assert code == 2
    assert "No NVIDIA API key" in capsys.readouterr().err


def test_backend_failure_keeps_partial_results(tmp_path, monkeypatch, capsys):
    calls = {"n": 0}

    class Flaky(cli.MockBackend):
        def score(self, seq):
            calls["n"] += 1
            if calls["n"] > 3:
                raise BackendError("simulated outage")
            return super().score(seq)

    monkeypatch.setitem(cli.BACKENDS, "mock", Flaky)
    out = tmp_path / "p"
    code = run("--vcf", EXAMPLES / "tiny.vcf", "--ref", EXAMPLES / "tiny.fa", "--backend", "mock",
               "--no-cache", "--out", out)
    assert code == 2
    assert "simulated outage" in capsys.readouterr().err
    assert len(read_rows(str(out) + ".csv")) >= 1
    assert "Partial results" in Path(str(out) + ".html").read_text(encoding="utf-8")


def test_failure_before_any_score_has_no_done_or_resume_message(tmp_path, monkeypatch, capsys):
    class Broken(cli.MockBackend):
        def score(self, seq):
            raise BackendError("NVIDIA rejected the API key (HTTP 401).")

    monkeypatch.setitem(cli.BACKENDS, "mock", Broken)
    out = tmp_path / "z"
    code = run("--vcf", EXAMPLES / "tiny.vcf", "--ref", EXAMPLES / "tiny.fa", "--backend", "mock",
               "--no-cache", "--out", out)
    err = capsys.readouterr().err
    assert code == 2
    assert "Stopped: NVIDIA rejected the API key" in err
    assert "Done" not in err and "resume" not in err
    assert not Path(str(out) + ".html").exists()


def test_limit_and_bad_inputs(tmp_path, capsys):
    out = tmp_path / "l"
    assert run("--vcf", EXAMPLES / "tiny.vcf", "--ref", EXAMPLES / "tiny.fa", "--backend", "mock",
               "--no-cache", "--limit", "2", "--out", out) == 0
    assert len(read_rows(str(out) + ".csv")) == 2
    assert run("--vcf", tmp_path / "missing.vcf", "--ref", EXAMPLES / "tiny.fa", "--backend", "mock") == 2
    assert run("--vcf", EXAMPLES / "tiny.vcf", "--ref", EXAMPLES / "tiny.fa", "--backend", "mock",
               "--window", "10") == 2
    assert "error:" in capsys.readouterr().err


def test_chr_prefix_mismatch_between_vcf_and_fasta(genome, tmp_path):
    path, seqs = genome
    ref = seqs["chr1"][99]
    alt = "A" if ref != "A" else "C"
    vcf = write_vcf(tmp_path / "n.vcf", [("1", 100, "x", ref, alt)])  # VCF says "1", FASTA says "chr1"
    out = tmp_path / "o"
    assert run("--vcf", vcf, "--ref", path, "--backend", "mock", "--no-cache", "--out", out) == 0
    (row,) = read_rows(str(out) + ".csv")
    assert row["status"] == "ok" and row["chrom"] == "1"


def test_nvidia_runs_print_a_data_notice_and_never_the_key(tmp_path, monkeypatch, capsys):
    class Stub(cli.MockBackend):
        name = "nvidia"
        cache_id = "stub-nvidia:v1"
        is_mock = False

        def __init__(self, api_key, rpm=30):
            pass

    monkeypatch.setenv("NVIDIA_API_KEY", "nvapi-" + "x" * 64)
    monkeypatch.setattr(cli, "NvidiaBackend", Stub)
    code = run("--vcf", EXAMPLES / "tiny.vcf", "--ref", EXAMPLES / "tiny.fa", "--backend", "nvidia",
               "--no-cache", "--limit", "1", "--out", tmp_path / "n")
    err = capsys.readouterr().err
    assert code == 0
    assert "sent to NVIDIA" in err and "identifiable patient" in err
    assert "nvapi-" not in err


def test_mock_runs_do_not_print_the_nvidia_notice(tmp_path, capsys):
    run("--vcf", EXAMPLES / "tiny.vcf", "--ref", EXAMPLES / "tiny.fa", "--backend", "mock",
        "--no-cache", "--limit", "1", "--out", tmp_path / "m")
    assert "sent to NVIDIA" not in capsys.readouterr().err
