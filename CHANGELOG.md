# Changelog

## 0.2.2 (2026-10-22)
- Notice printed before any sequence is sent to NVIDIA's API; new "Data privacy and NVIDIA's terms"
  and "Citing and acknowledgements" sections in the README.
- Added `NOTICE`, `CITATION.cff`, `SECURITY.md`, `.gitattributes`.
- ClinVar indel validation result (matched sample, AUROC 0.826) recorded in the README.
- `evaluate --drop NONE` replaces the empty `--drop ""` example, which PowerShell drops.

## 0.2.1 (2026-10-22)
- ClinVar indel benchmark tooling (`prepare-clinvar`) with matched sampling (`--per-stratum`) by
  variant type (coding frameshift, coding in-frame, non-coding).
- `evaluate` reports AUROC by variant type and compares summed and per-base scores.
- Genome index (`.fai`) is saved and reused; `fetch-ref --chroms` downloads several chromosomes.

## 0.2.0 (2026-10-22)
- BRCA1 validation tooling (`evo-vscore-brca1 prepare / evaluate`) and the first results
  (AUROC 0.873 on all 3,893 SNVs).

## 0.1.0 (2026-10-22)
- First version: `evo-vscore` command with `mock` and `nvidia` backends, window grid with shared
  reference scores, resumable cache, CSV and HTML output.
