"""Validate evo-vscore against public benchmarks with known answers.

Two benchmarks are supported. Neither dataset ships with this tool; you download
them yourself and respect their terms.

* **BRCA1 SNVs**: Findlay et al. 2018, Nature, saturation genome editing,
  Supplementary Table 1 (3,893 SNVs, hg19), labelled FUNC / INT / LOF.
* **ClinVar indels**: pathogenic versus benign insertions and deletions from
  NCBI's ``variant_summary.txt.gz``, the same kind of comparison the Evo 2 paper
  made for non-SNVs.

A good scorer gives the damaging class (LOF, or pathogenic) a more negative
``delta_ll``. Sub-commands (``evo-vscore-brca1 --help``):

* ``fetch-ref``       download hg19 chromosome(s) from UCSC and unpack them
* ``prepare``         BRCA1 spreadsheet -> ``brca1.vcf`` + ``brca1_labels.csv``
* ``prepare-clinvar`` ClinVar table -> ``clinvar.vcf`` + ``clinvar_labels.csv``
* ``evaluate``        compare an evo-vscore results CSV with the labels (AUROC)
"""

import argparse
import csv
import gzip
import os
import random
import re
import shutil
import sys
import urllib.request
from collections import Counter
from pathlib import Path

import numpy as np

UCSC_HG19_BASE = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/"
HG19_CHR17_URL = UCSC_HG19_BASE + "chr17.fa.gz"

_REQUIRED = {
    "chrom": "chromosome",
    "pos": "position (hg19)",
    "ref": "reference",
    "alt": "alt",
    "class": "func.class",
}
_OPTIONAL = {"score": "function.score.mean", "consequence": "consequence"}
_CLASSES = {"FUNC", "INT", "LOF"}
_STANDARD_CHROMS = [str(i) for i in range(1, 23)] + ["X"]


# -- reference download -------------------------------------------------------
def fetch_reference(url, out_path):
    """Download one .fa.gz file and unpack it to ``out_path``."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    part = Path(str(out_path) + ".gz.part")
    with urllib.request.urlopen(url, timeout=120) as response, open(part, "wb") as fh:
        shutil.copyfileobj(response, fh)
    with gzip.open(part, "rb") as src, open(out_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    os.remove(part)
    return out_path


def parse_chroms(text):
    """'17' -> ['17'];  '1-3,X' -> ['1','2','3','X'];  'all' -> 1-22 and X."""
    text = text.strip()
    if text.lower() == "all":
        return list(_STANDARD_CHROMS)
    chroms = []
    for token in (t.strip() for t in text.split(",") if t.strip()):
        match = re.fullmatch(r"(\d+)-(\d+)", token)
        if match:
            low, high = int(match.group(1)), int(match.group(2))
            if low > high or high > 22:
                raise ValueError("Bad chromosome range: %r" % token)
            chroms.extend(str(i) for i in range(low, high + 1))
        elif re.fullmatch(r"\d+|[XYxy]", token):
            chroms.append(token.upper())
        else:
            raise ValueError("Bad chromosome %r (use numbers, X, Y, ranges like 1-22, or 'all')" % token)
    if not chroms:
        raise ValueError("No chromosomes given")
    return chroms


def fetch_chromosomes(chroms, out_path, base_url=UCSC_HG19_BASE, progress=print):
    """Download chr<N>.fa.gz for each chromosome and join them into one FASTA."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    base = base_url if base_url.endswith("/") else base_url + "/"
    part_gz = Path(str(out_path) + ".gz.part")
    part_out = Path(str(out_path) + ".part")
    with open(part_out, "wb") as dst:
        for chrom in chroms:
            progress("downloading chr%s ..." % chrom)
            with urllib.request.urlopen("%schr%s.fa.gz" % (base, chrom), timeout=300) as response, \
                    open(part_gz, "wb") as fh:
                shutil.copyfileobj(response, fh)
            last = b"\n"
            with gzip.open(part_gz, "rb") as src:
                while True:
                    chunk = src.read(1 << 20)
                    if not chunk:
                        break
                    dst.write(chunk)
                    last = chunk[-1:]
            if last != b"\n":  # keep the next '>' header on its own line
                dst.write(b"\n")
            os.remove(part_gz)
    os.replace(part_out, out_path)
    return out_path


# -- BRCA1: reading the spreadsheet ---------------------------------------------
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


_LABEL_FIELDS = ["id", "chrom", "pos", "ref", "alt", "class", "score", "consequence", "gene", "stars", "clinvar_type"]


def _write_benchmark_files(rows, out_dir, vcf_name, labels_name, id_prefix, source):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for i, row in enumerate(rows, start=1):
        row["id"] = "%s%04d" % (id_prefix, i)
    with open(out / vcf_name, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.2\n##source=%s\n#CHROM\tPOS\tID\tREF\tALT\n" % source)
        for r in rows:
            fh.write("%s\t%d\t%s\t%s\t%s\n" % (r["chrom"], r["pos"], r["id"], r["ref"], r["alt"]))
    with open(out / labels_name, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=_LABEL_FIELDS)
        writer.writeheader()
        writer.writerows({k: r.get(k, "") for k in _LABEL_FIELDS} for r in rows)


def prepare(xlsx_path, out_dir, sample=None, seed=42):
    rows, skipped = read_findlay_xlsx(xlsx_path)
    if not rows:
        raise ValueError("No usable SNV rows were found in the spreadsheet")
    if sample is not None and sample < len(rows):
        rows = random.Random(seed).sample(rows, sample)
    rows.sort(key=lambda r: (r["chrom"], r["pos"], r["alt"]))
    _write_benchmark_files(rows, out_dir, "brca1.vcf", "brca1_labels.csv", "b",
                           "evo-vscore prepare (Findlay 2018 SGE SNVs, hg19)")
    counts = {c: sum(1 for r in rows if r["class"] == c) for c in sorted(_CLASSES)}
    return {"n": len(rows), "classes": counts, "skipped": skipped, "dir": str(out_dir)}


# -- ClinVar indels -----------------------------------------------------------------
_STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
}
_PATHOGENIC = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}
_BENIGN = {"benign", "likely benign", "benign/likely benign"}
_INDEL_TYPES = {"deletion", "insertion", "indel", "duplication"}
_NAME_FIELDS = {
    "type": ["type"],
    "name": ["name"],
    "gene": ["genesymbol"],
    "assembly": ["assembly"],
    "chrom": ["chromosome"],
    "significance": ["clinicalsignificance", "germlineclassification"],
    "review": ["reviewstatus", "germlinereviewstatus"],
    "pos": ["positionvcf"],
    "ref": ["referenceallelevcf"],
    "alt": ["alternateallelevcf"],
}


def approx_region(name):
    """Rough coding / non-coding call from a ClinVar HGVS name (a heuristic, not an annotation).

    A plain ``c.<number>`` position counts as coding. Intronic offsets (``c.123+5``,
    ``c.123-5``), 5' and 3' UTR positions (``c.-20``, ``c.*5``) and non-coding
    transcripts count as non-coding.
    """
    match = re.search(r":c\.(\S+)", name or "")
    if not match:
        return "noncoding_approx"
    position = re.match(r"[-*+\d_?]+", match.group(1))
    if position is None:
        return "noncoding_approx"
    text = position.group(0)
    if text.startswith(("-", "*")) or "+" in text or "-" in text:
        return "noncoding_approx"
    return "coding_approx"


def _open_text(path):
    if str(path).lower().endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return open(path, "r", encoding="utf-8", errors="replace", newline="")


def read_clinvar(path, assembly="GRCh37", min_stars=1, chroms=None, max_indel=50):
    """Stream ClinVar's variant_summary table; return (records, dropped_counts)."""
    if not Path(path).exists():
        raise FileNotFoundError("ClinVar file not found: %s" % path)
    allowed = set(chroms or _STANDARD_CHROMS)
    dropped = Counter()
    records, conflicted = {}, set()
    with _open_text(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = [h.lstrip("#").strip().lower() for h in next(reader, [])]
        col = {}
        for key, names in _NAME_FIELDS.items():
            for name in names:
                if name in header:
                    col[key] = header.index(name)
                    break
        missing = [k for k in _NAME_FIELDS if k not in col]
        if missing:
            raise ValueError("ClinVar table lacks expected column(s) for %s. Columns found: %s"
                             % (missing, header))
        width = max(col.values()) + 1
        for row in reader:
            if len(row) < width or row[col["assembly"]] != assembly:
                continue
            ctype = row[col["type"]].strip().lower()
            if ctype not in _INDEL_TYPES:
                continue
            significance = row[col["significance"]].strip().lower()
            if significance in _PATHOGENIC:
                label = "PATH"
            elif significance in _BENIGN:
                label = "BENIGN"
            else:
                dropped["not_clearly_pathogenic_or_benign"] += 1
                continue
            stars = _STARS.get(row[col["review"]].strip().lower(), 0)
            if stars < min_stars:
                dropped["low_review_status"] += 1
                continue
            chrom = row[col["chrom"]].strip()
            if chrom not in allowed:
                dropped["chromosome_not_selected"] += 1
                continue
            ref, alt = row[col["ref"]].strip().upper(), row[col["alt"]].strip().upper()
            pos_text = row[col["pos"]].strip()
            if (not pos_text.isdigit() or not re.fullmatch(r"[ACGT]+", ref) or not re.fullmatch(r"[ACGT]+", alt)
                    or ref == alt or len(ref) == len(alt)):
                dropped["unusable_alleles"] += 1
                continue
            if abs(len(ref) - len(alt)) > max_indel or len(ref) > 2 * max_indel + 2:
                dropped["longer_than_max_indel"] += 1
                continue
            key = (chrom, int(pos_text), ref, alt)
            if key in conflicted:
                continue
            if key in records:
                if records[key]["class"] != label:
                    conflicted.add(key)
                    del records[key]
                    dropped["conflicting_duplicate"] += 1
                continue
            region = approx_region(row[col["name"]])
            if region == "coding_approx":  # a length change that is not a multiple of 3 shifts the reading frame
                region = "coding_frameshift_approx" if (len(alt) - len(ref)) % 3 else "coding_inframe_approx"
            records[key] = {
                "chrom": chrom, "pos": int(pos_text), "ref": ref, "alt": alt, "class": label,
                "score": "", "consequence": region, "gene": row[col["gene"]].strip(),
                "stars": stars, "clinvar_type": row[col["type"]].strip(),
            }
    return list(records.values()), dict(dropped)


def _chrom_order(chrom):
    return int(chrom) if chrom.isdigit() else 100 + ord(chrom[0])


def _sort_key(r):
    return (_chrom_order(r["chrom"]), r["pos"], r["alt"])


def prepare_clinvar(summary_path, out_dir, assembly="GRCh37", min_stars=1, chroms=None,
                    max_per_class=300, seed=42, max_indel=50, per_stratum=None):
    """Write a ClinVar indel benchmark.

    With ``per_stratum=N`` the sample is *matched*: within each variant type (coding
    frameshift, coding in-frame, non-coding) it takes the same number of pathogenic and
    benign variants, up to N each, so variant type alone cannot predict the label.
    Without it, up to ``max_per_class`` are drawn at random per class.
    """
    records, dropped = read_clinvar(summary_path, assembly, min_stars, chroms, max_indel)
    if not records:
        raise ValueError("No usable ClinVar indels were found (dropped: %s)" % dropped)
    available = Counter(r["class"] for r in records)
    strata = {}
    for r in records:
        strata.setdefault(r["consequence"], {"PATH": [], "BENIGN": []})[r["class"]].append(r)
    available_strata = {s: {c: len(v) for c, v in groups.items()} for s, groups in sorted(strata.items())}
    chosen, taken = [], {}
    if per_stratum:
        rng = random.Random(seed)
        for stratum, groups in sorted(strata.items()):
            n = min(per_stratum, len(groups["PATH"]), len(groups["BENIGN"]))
            taken[stratum] = n
            for label in ("PATH", "BENIGN"):
                pool = sorted(groups[label], key=_sort_key)
                chosen.extend(rng.sample(pool, n) if n < len(pool) else pool)
        if not chosen:
            raise ValueError("No variant type has both pathogenic and benign variants (available: %s)"
                             % available_strata)
    else:
        for label in sorted(available):
            group = sorted((r for r in records if r["class"] == label), key=_sort_key)
            if max_per_class and len(group) > max_per_class:
                group = random.Random(seed).sample(group, max_per_class)
            chosen.extend(group)
    chosen.sort(key=_sort_key)
    _write_benchmark_files(chosen, out_dir, "clinvar.vcf", "clinvar_labels.csv", "c",
                           "evo-vscore prepare-clinvar (%s indels)" % assembly)
    return {
        "n": len(chosen),
        "available": dict(available),
        "classes": dict(Counter(r["class"] for r in chosen)),
        "regions": dict(Counter((r["class"], r["consequence"]) for r in chosen)),
        "available_strata": available_strata,
        "taken_per_stratum": taken,
        "chromosomes": sorted({r["chrom"] for r in chosen}, key=_chrom_order),
        "dropped": dropped,
        "dir": str(out_dir),
    }


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
    return {"name": label, "n": len(y), "n_pos": int(np.sum(y)),
            "auroc": auroc(y, s), "ci_low": low, "ci_high": high}


def _float_or_nan(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return float("nan")


def evaluate(results_csv, labels_csv, n_boot=1000, positive="LOF", drop=("INT",)):
    """Compare scores with labels. ``positive`` is the damaging class (LOF, PATH...);
    ``drop`` lists classes to leave out of the second, stricter comparison."""
    labels = {r["id"]: r for r in _read_csv(labels_csv)}
    status_counts, joined = {}, []
    for r in _read_csv(results_csv):
        label = labels.get(r.get("id"))
        if label is None:
            continue
        status_counts[r["status"]] = status_counts.get(r["status"], 0) + 1
        if r["status"] == "ok" and r["delta_ll"] != "":
            joined.append((label, r))
    if not joined:
        raise ValueError("None of the variants in the results file were scored (statuses: %s)" % status_counts)

    classes = np.array([lab["class"] for lab, _ in joined])
    delta = np.array([float(r["delta_ll"]) for _, r in joined])
    per_base = np.array([_float_or_nan(r.get("delta_mean_ll")) for _, r in joined])
    types = np.array([r.get("variant_type", "") or "" for _, r in joined])
    is_pos = classes == positive
    if not is_pos.any():
        raise ValueError("No scored variants have class %r. Classes present: %s"
                         % (positive, sorted(set(classes))))
    names = sorted(set(classes))
    others = [c for c in names if c != positive]
    if not others:
        raise ValueError("Only class %r is present, so there is nothing to compare it with" % positive)

    out = {
        "scored": len(joined), "statuses": status_counts, "positive": positive,
        "class_counts": {c: int((classes == c).sum()) for c in names},
        "mean_delta": {c: float(delta[classes == c].mean()) for c in names},
        "measures": [],
    }
    # a more negative delta should mean "more likely damaging", so rank by -delta
    primary = "%s vs %s (all)" % (positive, "+".join(others))
    out["measures"].append(_measure(primary, is_pos, -delta, n_boot))

    distinct_types = {t for t in types if t}
    if np.isfinite(per_base).all() and (distinct_types - {"SNV"}):
        out["measures"].append(_measure(primary + " [per-base delta_mean_ll]", is_pos, -per_base, n_boot))

    dropped = [c for c in drop if c in others]
    if dropped and len(others) > len(dropped):
        keep = ~np.isin(classes, dropped)
        rest = [c for c in others if c not in dropped]
        out["measures"].append(_measure(
            "%s vs %s only (%s excluded)" % (positive, "+".join(rest), "+".join(dropped)),
            is_pos[keep], -delta[keep], n_boot))

    def add_groups(kind, values):
        for name in sorted(set(values) - {""}):
            mask = values == name
            if is_pos[mask].sum() >= 10 and (~is_pos[mask]).sum() >= 10:
                out["measures"].append(_measure("%s: %s" % (kind, name), is_pos[mask], -delta[mask], n_boot))

    if len(distinct_types) > 1:
        add_groups("by variant type", types)
    consequences = np.array([lab.get("consequence", "") or "" for lab, _ in joined])
    add_groups("by consequence", consequences)
    return out


def format_summary(result, title="Validation summary"):
    lines = ["# " + title, ""]
    lines.append("Variants scored: %d. Classes: %s. Positive (damaging) class: %s." % (
        result["scored"], ", ".join("%s=%d" % kv for kv in result["class_counts"].items()),
        result.get("positive", "LOF")))
    bad = {k: v for k, v in result["statuses"].items() if k != "ok"}
    if bad:
        lines.append("Not scored: %s." % ", ".join("%s=%d" % kv for kv in sorted(bad.items())))
    lines.append("Mean delta_ll by class: %s." % ", ".join("%s=%.3f" % kv for kv in result["mean_delta"].items()))
    lines += ["", "| Comparison | n | Damaging | AUROC | 95% CI (bootstrap) |", "|---|---|---|---|---|"]
    for m in result["measures"]:
        lines.append("| %s | %d | %d | %.3f | %.3f to %.3f |" % (
            m["name"], m["n"], m["n_pos"], m["auroc"], m["ci_low"], m["ci_high"]))
    lines += ["", "AUROC 0.5 = no better than chance, 1.0 = perfect separation of the damaging class from the rest."]
    return "\n".join(lines) + "\n"


# -- command line ----------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(prog="evo-vscore-brca1", description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="command", required=True)

    f = sub.add_parser("fetch-ref", help="download hg19 chromosome(s) from UCSC and unpack them")
    f.add_argument("--out", required=True, help="where to write the unpacked FASTA, e.g. data\\chr17.fa")
    f.add_argument("--chroms", default="17", help="e.g. 17, or 1-22,X, or all (default: %(default)s)")
    f.add_argument("--base-url", default=UCSC_HG19_BASE)
    f.add_argument("--url", help="download this single .fa.gz file instead of using --chroms")

    r = sub.add_parser("prepare", help="turn the BRCA1 Supplementary Table 1 into a VCF plus a labels file")
    r.add_argument("--xlsx", required=True, help="41586_2018_461_MOESM3_ESM.xlsx")
    r.add_argument("--out-dir", required=True)
    r.add_argument("--sample", type=int, help="random subset of N variants (default: all)")
    r.add_argument("--seed", type=int, default=42)

    c = sub.add_parser("prepare-clinvar", help="turn ClinVar's variant_summary.txt.gz into an indel benchmark")
    c.add_argument("--summary-file", required=True, help="variant_summary.txt.gz from NCBI ClinVar")
    c.add_argument("--out-dir", required=True)
    c.add_argument("--assembly", default="GRCh37", help="GRCh37 matches the hg19 reference (default)")
    c.add_argument("--chroms", default="1-22,X", help="chromosomes to keep (default: %(default)s)")
    c.add_argument("--min-stars", type=int, default=1, help="minimum ClinVar review stars (default: %(default)s)")
    c.add_argument("--max-per-class", type=int, default=300, help="cap per class when not matching, 0 = no cap (default: %(default)s)")
    c.add_argument("--per-stratum", type=int, help="MATCHED sampling: up to N pathogenic AND N benign within each "
                   "variant type (coding frameshift / coding in-frame / non-coding). Recommended.")
    c.add_argument("--max-indel", type=int, default=50, help="longest insertion/deletion in bases (default: %(default)s)")
    c.add_argument("--seed", type=int, default=42)

    e = sub.add_parser("evaluate", help="compute AUROC from an evo-vscore results CSV")
    e.add_argument("--results", required=True)
    e.add_argument("--labels", required=True)
    e.add_argument("--positive", default="LOF", help="the damaging class: LOF for BRCA1, PATH for ClinVar")
    e.add_argument("--drop", default="INT", help="classes to leave out of the stricter comparison; "
                   "use NONE for none (PowerShell drops an empty \"\" argument)")
    e.add_argument("--summary", help="also write the summary to this .md file")
    e.add_argument("--title", default="Validation summary")
    e.add_argument("--bootstrap", type=int, default=1000)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        if args.command == "fetch-ref":
            if args.url:
                path = fetch_reference(args.url, args.out)
            else:
                path = fetch_chromosomes(parse_chroms(args.chroms), args.out, args.base_url)
            print("Reference written to %s (%.1f MB)" % (path, Path(path).stat().st_size / 1e6))
        elif args.command == "prepare":
            info = prepare(args.xlsx, args.out_dir, args.sample, args.seed)
            print("Wrote %s\\brca1.vcf and brca1_labels.csv: %d SNVs %s" % (info["dir"], info["n"], info["classes"]))
            print("Skipped rows: %s" % info["skipped"])
        elif args.command == "prepare-clinvar":
            info = prepare_clinvar(args.summary_file, args.out_dir, args.assembly, args.min_stars,
                                   parse_chroms(args.chroms), args.max_per_class, args.seed, args.max_indel,
                                   args.per_stratum)
            print("Wrote %s\\clinvar.vcf and clinvar_labels.csv: %d indels %s" % (info["dir"], info["n"], info["classes"]))
            print("Available before sampling: %s" % info["available"])
            print("Available per variant type (approximate, from the HGVS name):")
            for stratum, counts in info["available_strata"].items():
                note = ("  -> using %d of each" % info["taken_per_stratum"][stratum]) if info["taken_per_stratum"] else ""
                print("  %-26s PATH=%-7d BENIGN=%-7d%s" % (stratum, counts["PATH"], counts["BENIGN"], note))
            print("Sample by class and type: %s" % {"%s/%s" % k: v for k, v in sorted(info["regions"].items())})
            print("Chromosomes: %s" % ", ".join(info["chromosomes"]))
            print("Dropped: %s" % info["dropped"])
            thin = [s for s, c in info["available_strata"].items() if min(c["PATH"], c["BENIGN"]) < 10]
            if info["taken_per_stratum"] and thin:
                print("note: these types lack 10+ of one class and cannot be compared: %s" % ", ".join(thin))
            if min(info["classes"].values()) < 30 or len(info["classes"]) < 2:
                print("warning: a class has fewer than 30 variants, so the AUROC will be very uncertain")
        else:
            drop = tuple(x.strip() for x in args.drop.split(",") if x.strip())
            result = evaluate(args.results, args.labels, args.bootstrap, args.positive, drop)
            text = format_summary(result, args.title)
            print(text)
            if args.summary:
                Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
                Path(args.summary).write_text(text, encoding="utf-8")
                print("Summary saved to %s" % args.summary)
        return 0
    except (ValueError, OSError) as error:
        print("error: %s" % error, file=sys.stderr)
        if args.command == "fetch-ref":
            print("If the download is blocked, open %s in a browser, save the chromosome file(s), "
                  "and unpack them with 7-Zip." % UCSC_HG19_BASE, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
