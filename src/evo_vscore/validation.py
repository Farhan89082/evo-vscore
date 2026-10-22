"""Validate evo-vscore on the BRCA1 saturation-genome-editing (SGE) data.

Benchmark: Findlay et al. 2018, Nature, "Accurate classification of BRCA1 variants
with saturation genome editing", Supplementary Table 1 (3,893 SNVs, hg19). Each
SNV is labelled FUNC (functional), INT (intermediate) or LOF (loss of function).
A good scorer gives LOF variants a more negative ``delta_ll``.

This module does not ship the data. You download the spreadsheet yourself from
the paper's page; see the README for the steps.

Three sub-commands (``evo-vscore-brca1 --help``):

* ``fetch-ref``  download and unpack the hg19 chr17 reference
* ``prepare``    turn the spreadsheet into ``brca1.vcf`` + ``brca1_labels.csv``
* ``evaluate``   compare an evo-vscore results CSV with the labels (AUROC)
"""

import argparse
import csv
import gzip
import os
import random
import shutil
import sys
import urllib.request
from pathlib import Path

import numpy as np

HG19_CHR17_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr17.fa.gz"

_REQUIRED = {
    "chrom": "chromosome",
    "pos": "position (hg19)",
    "ref": "reference",
    "alt": "alt",
    "class": "func.class",
}
_OPTIONAL = {"score": "function.score.mean", "consequence": "consequence"}
_CLASSES = {"FUNC", "INT", "LOF"}


# -- reference download -------------------------------------------------------
def fetch_reference(url, out_path):
    """Download a .fa.gz file and unpack it to ``out_path``."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    part = Path(str(out_path) + ".gz.part")
    with urllib.request.urlopen(url, timeout=120) as response, open(part, "wb") as fh:
        shutil.copyfileobj(response, fh)
    with gzip.open(part, "rb") as src, open(out_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    os.remove(part)
    return out_path


# -- reading the spreadsheet ----------------------------------------------------
def read_findlay_xlsx(path):
    """Return (rows, skipped_counts) from Supplementary Table 1.

    Only single-base substitutions with a FUNC/INT/LOF label are kept.
    """
    try:
        import openpyxl
    except ImportError:
        raise ValueError("Reading .xlsx needs openpyxl. Install it with:  pip install openpyxl")
    if not Path(path).exists():
        raise FileNotFoundError("Spreadsheet not found: %s" % path)

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    for sheet in workbook.worksheets:
        row_iter = sheet.iter_rows(values_only=True)
        header = None
        for number, row in enumerate(row_iter):
            cells = [str(c).strip().lower() if c is not None else "" for c in row]
            if _REQUIRED["pos"] in cells and _REQUIRED["class"] in cells:
                header = cells
                break
            if number > 25:
                break
        if header is None:
            continue
        col = {key: header.index(name) for key, name in _REQUIRED.items() if name in header}
        missing = [name for key, name in _REQUIRED.items() if key not in col]
        if missing:
            raise ValueError("Spreadsheet header lacks column(s): %s. Found: %s" % (missing, header))
        for key, name in _OPTIONAL.items():
            if name in header:
                col[key] = header.index(name)

        rows, skipped = [], {"not_a_snv": 0, "unlabelled": 0, "incomplete": 0}
        for row in row_iter:
            def cell(key):
                i = col.get(key)
                return row[i] if i is not None and i < len(row) else None

            if cell("chrom") is None or cell("pos") is None:
                skipped["incomplete"] += 1
                continue
            ref, alt = str(cell("ref")).strip().upper(), str(cell("alt")).strip().upper()
            label = str(cell("class")).strip().upper()
            if label not in _CLASSES:
                skipped["unlabelled"] += 1
                continue
            if len(ref) != 1 or len(alt) != 1 or ref == alt or not set(ref + alt) <= set("ACGT"):
                skipped["not_a_snv"] += 1
                continue
            score = cell("score")
            rows.append({
                "chrom": str(cell("chrom")).strip(),
                "pos": int(float(cell("pos"))),
                "ref": ref,
                "alt": alt,
                "class": label,
                "score": "" if score is None else score,
                "consequence": "" if cell("consequence") is None else str(cell("consequence")).strip(),
            })
        return rows, skipped
    raise ValueError(
        "Could not find the table header (a row containing 'position (hg19)' and 'func.class') "
        "in any sheet of %s. Is this Supplementary Table 1 (41586_2018_461_MOESM3_ESM.xlsx)?" % path
    )


def prepare(xlsx_path, out_dir, sample=None, seed=42):
    rows, skipped = read_findlay_xlsx(xlsx_path)
    if not rows:
        raise ValueError("No usable SNV rows were found in the spreadsheet")
    if sample is not None and sample < len(rows):
        rows = random.Random(seed).sample(rows, sample)
    rows.sort(key=lambda r: (r["chrom"], r["pos"], r["alt"]))
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for i, row in enumerate(rows, start=1):
        row["id"] = "b%04d" % i
    with open(out / "brca1.vcf", "w", newline="\n", encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.2\n##source=evo-vscore prepare (Findlay 2018 SGE SNVs, hg19)\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\n")
        for r in rows:
            fh.write("%s\t%d\t%s\t%s\t%s\n" % (r["chrom"], r["pos"], r["id"], r["ref"], r["alt"]))
    fields = ["id", "chrom", "pos", "ref", "alt", "class", "score", "consequence"]
    with open(out / "brca1_labels.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows({k: r[k] for k in fields} for r in rows)
    counts = {c: sum(1 for r in rows if r["class"] == c) for c in sorted(_CLASSES)}
    return {"n": len(rows), "classes": counts, "skipped": skipped, "dir": str(out)}


# -- statistics -------------------------------------------------------------------
def _average_ranks(x):
    _, inverse, counts = np.unique(x, return_inverse=True, return_counts=True)
    last = np.cumsum(counts)
    first = last - counts + 1
    return ((first + last) / 2.0)[inverse]


def auroc(is_positive, score):
    """Area under the ROC curve: P(random positive scores higher than random negative).

    Ties count as one half. Computed from ranks, so it needs no extra packages.
    """
    y = np.asarray(is_positive, dtype=bool)
    s = np.asarray(score, dtype=float)
    n_pos, n_neg = int(y.sum()), int((~y).sum())
    if n_pos == 0 or n_neg == 0:
        raise ValueError("AUROC needs at least one positive and one negative example")
    ranks = _average_ranks(s)
    return float((ranks[y].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def bootstrap_ci(is_positive, score, n_boot=1000, seed=0):
    y = np.asarray(is_positive, dtype=bool)
    s = np.asarray(score, dtype=float)
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        if y[idx].any() and (~y[idx]).any():
            values.append(auroc(y[idx], s[idx]))
    if not values:
        return (float("nan"), float("nan"))
    return (float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5)))


# -- evaluation ---------------------------------------------------------------------
def _read_csv(path):
    if not Path(path).exists():
        raise FileNotFoundError("File not found: %s" % path)
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _measure(label, y, s, n_boot):
    low, high = bootstrap_ci(y, s, n_boot)
    return {"name": label, "n": len(y), "n_lof": int(np.sum(y)),
            "auroc": auroc(y, s), "ci_low": low, "ci_high": high}


def evaluate(results_csv, labels_csv, n_boot=1000):
    labels = {r["id"]: r for r in _read_csv(labels_csv)}
    results = _read_csv(results_csv)
    status_counts, joined = {}, []
    for r in results:
        label = labels.get(r.get("id"))
        if label is None:
            continue
        status_counts[r["status"]] = status_counts.get(r["status"], 0) + 1
        if r["status"] == "ok" and r["delta_ll"] != "":
            joined.append((label, float(r["delta_ll"])))
    if not joined:
        raise ValueError("None of the variants in the results file were scored (statuses: %s)" % status_counts)

    classes = np.array([j[0]["class"] for j in joined])
    delta = np.array([j[1] for j in joined])
    lof = classes == "LOF"
    # a more negative delta should mean "more likely LOF", so rank by -delta
    out = {"scored": len(joined), "statuses": status_counts,
           "class_counts": {c: int((classes == c).sum()) for c in sorted(_CLASSES)},
           "mean_delta": {c: float(delta[classes == c].mean()) for c in sorted(_CLASSES) if (classes == c).any()},
           "measures": []}
    out["measures"].append(_measure("LOF vs FUNC+INT (all)", lof, -delta, n_boot))
    keep = classes != "INT"
    if lof[keep].any() and (~lof[keep]).any():
        out["measures"].append(_measure("LOF vs FUNC only (INT excluded)", lof[keep], -delta[keep], n_boot))

    groups = {}
    for (label, _), is_lof, d in zip(joined, lof, delta):
        if label.get("consequence"):
            groups.setdefault(label["consequence"], []).append((is_lof, d))
    for name, items in sorted(groups.items()):
        y = np.array([i[0] for i in items])
        if y.sum() >= 10 and (~y).sum() >= 10:
            out["measures"].append(_measure("by consequence: %s" % name, y, -np.array([i[1] for i in items]), n_boot))
    return out


def format_summary(result, backend_note=""):
    lines = ["# BRCA1 validation (Findlay 2018 SGE SNVs)", ""]
    if backend_note:
        lines += [backend_note, ""]
    lines.append("Variants scored: %d. Classes: %s." % (
        result["scored"], ", ".join("%s=%d" % kv for kv in result["class_counts"].items())))
    bad = {k: v for k, v in result["statuses"].items() if k != "ok"}
    if bad:
        lines.append("Not scored: %s." % ", ".join("%s=%d" % kv for kv in sorted(bad.items())))
    lines.append("Mean delta_ll by class: %s." % ", ".join("%s=%.3f" % kv for kv in result["mean_delta"].items()))
    lines += ["", "| Comparison | n | LOF | AUROC | 95% CI (bootstrap) |", "|---|---|---|---|---|"]
    for m in result["measures"]:
        lines.append("| %s | %d | %d | %.3f | %.3f to %.3f |" % (
            m["name"], m["n"], m["n_lof"], m["auroc"], m["ci_low"], m["ci_high"]))
    lines += ["", "AUROC 0.5 = no better than chance, 1.0 = perfect separation of LOF from the rest."]
    return "\n".join(lines) + "\n"


# -- command line ----------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(prog="evo-vscore-brca1", description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="command", required=True)

    f = sub.add_parser("fetch-ref", help="download and unpack the hg19 chr17 reference FASTA")
    f.add_argument("--out", required=True, help="where to write the unpacked FASTA, e.g. data\\chr17.fa")
    f.add_argument("--url", default=HG19_CHR17_URL)

    r = sub.add_parser("prepare", help="turn Supplementary Table 1 into a VCF plus a labels file")
    r.add_argument("--xlsx", required=True, help="41586_2018_461_MOESM3_ESM.xlsx")
    r.add_argument("--out-dir", required=True)
    r.add_argument("--sample", type=int, help="random subset of N variants (default: all)")
    r.add_argument("--seed", type=int, default=42)

    e = sub.add_parser("evaluate", help="compute AUROC from an evo-vscore results CSV")
    e.add_argument("--results", required=True)
    e.add_argument("--labels", required=True)
    e.add_argument("--summary", help="also write the summary to this .md file")
    e.add_argument("--bootstrap", type=int, default=1000)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        if args.command == "fetch-ref":
            path = fetch_reference(args.url, args.out)
            print("Reference written to %s (%.1f MB)" % (path, Path(path).stat().st_size / 1e6))
        elif args.command == "prepare":
            info = prepare(args.xlsx, args.out_dir, args.sample, args.seed)
            print("Wrote %s\\brca1.vcf and brca1_labels.csv: %d SNVs %s" % (info["dir"], info["n"], info["classes"]))
            print("Skipped rows: %s" % info["skipped"])
        else:
            result = evaluate(args.results, args.labels, args.bootstrap)
            text = format_summary(result)
            print(text)
            if args.summary:
                Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
                Path(args.summary).write_text(text, encoding="utf-8")
                print("Summary saved to %s" % args.summary)
        return 0
    except (ValueError, OSError) as error:
        print("error: %s" % error, file=sys.stderr)
        if args.command == "fetch-ref":
            print("If the download is blocked, open %s in a browser, save chr17.fa.gz, "
                  "and unpack it with 7-Zip." % HG19_CHR17_URL, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
