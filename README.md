# evo-vscore

Score DNA variants from a VCF file with a DNA language model. For each SNV or indel the tool asks
the model how plausible the DNA looks with the variant versus without it, and reports the
difference. It writes a CSV and a self-contained HTML report.

> **Status: v0.2, alpha.** Validated so far on one benchmark: BRCA1 SNVs (see
> [Validation on BRCA1](#validation-on-brca1)). Indels and other genes are **not yet validated**.
> **Research use only. Not for clinical decisions.**

The first real backend uses [Evo 2](https://github.com/ArcInstitute/evo2) (Arc Institute) through
NVIDIA's hosted API, so no GPU is needed on your computer. Evo 2: Nature (2026),
[doi:10.1038/s41586-026-10176-5](https://doi.org/10.1038/s41586-026-10176-5).

This tool only **scores** existing variants. It does not generate sequences.

## Install (Windows PowerShell, also fine on macOS/Linux)

Needs Python 3.9 or newer.

```powershell
cd D:\projects\evo-vscore
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -e ".[dev]"
pytest
```

## Try it in one minute (no key needed)

The `mock` backend produces **fake** numbers so you can see the workflow and the output files:

```powershell
evo-vscore --vcf examples\tiny.vcf --ref examples\tiny.fa --backend mock --out out\demo
```

This writes `out\demo.csv` and `out\demo.html`. The example genome is synthetic.

## Real scores with NVIDIA's hosted Evo 2 7B

1. Create a free account at <https://build.nvidia.com>, open the
   [Evo 2 7B page](https://build.nvidia.com/arc/evo2-7b-forward) and click **Get API Key**
   (the key starts with `nvapi-`).
2. Set it for the current terminal only. This keeps it out of your files and shell history:

   ```powershell
   $env:NVIDIA_API_KEY = Read-Host "Paste your key"
   ```

   (`NVCF_KEY` also works as the variable name.) **Never commit a key to GitHub.**
3. Check the size of a job before running it, then run it:

   ```powershell
   evo-vscore --vcf my_variants.vcf --ref genome.fa --backend nvidia --dry-run
   evo-vscore --vcf my_variants.vcf --ref genome.fa --backend nvidia --out out\my_run
   ```

Your own free NVIDIA account has its own rate limit (about 40 requests a minute at the time of
writing; the tool defaults to 30 to stay under it). Check NVIDIA's current terms for your use.

## Inputs

| Input | Notes |
|---|---|
| `--vcf` | Plain or `.gz`. Only CHROM, POS, ID, REF, ALT are used. Multi-allelic records are split. |
| `--ref` | Uncompressed FASTA of the **same genome build** as the VCF. A `.fai` index next to it is used if present. `chr1` and `1` naming both work. |
| `--backend` | `nvidia` for real scores, `mock` for fake test numbers. Required, so you never mistake one for the other. |

## Options

| Option | Default | What it does |
|---|---|---|
| `--out` | `results` | Output prefix: writes `PREFIX.csv` and `PREFIX.html`. |
| `--window` | 8192 | Bases of reference around each variant that are scored. |
| `--rpm` | 30 | Max API requests per minute (nvidia backend). |
| `--dry-run` | off | Count variants and API calls and estimate the time. Scores nothing, needs no key. |
| `--limit N` | all | Only the first N variants. Handy for a trial run. |
| `--cache FILE` / `--no-cache` | `~/.cache/evo-vscore/scores.sqlite` | Where finished scores are kept. |

## Output columns

| Column | Meaning |
|---|---|
| `status` | `ok`, or why the variant was not scored (see below). |
| `ref_ll`, `alt_ll` | Model log-likelihood of the reference window and of the window with the variant. |
| `delta_ll` | `alt_ll - ref_ll`. **More negative = the model finds the variant less plausible.** |
| `delta_mean_ll` | The same difference per base. For insertions and deletions the two windows have different lengths, so this is given as well. Which of the two suits indels better is a question for the validation work in v0.2. |
| `window_start`, `window_end` | 1-based, inclusive reference coordinates of the scored window. |

Variants that are not scored are still listed, with a reason: `ref_mismatch` (the VCF REF does not
match your FASTA, usually a wrong genome build), `chrom_not_found`, `pos_out_of_range`,
`ambiguous_bases` (an N in the window), `too_long` (the deletion does not fit in the window),
`skipped:symbolic_or_breakend_allele` (such as `<DEL>`).

## How it works

- One reference window is cut around each variant. Windows sit on a grid, so nearby variants share
  one reference score, which is computed once.
- Every score is cached on disk. If a run is interrupted (Ctrl-C, network loss, rate limit), run the
  **same command again** and it resumes; finished scores are not requested twice.
- Rate-limit and server errors are retried with waiting. Bad keys, missing credits and rejected
  requests stop the run with a clear message and keep the partial CSV.

## Validation on BRCA1

`evo-vscore-brca1` checks the scores against the BRCA1 saturation-genome-editing data of
[Findlay et al. 2018](https://www.nature.com/articles/s41586-018-0461-z) (3,893 SNVs, each labelled
FUNC, INT or LOF). The data is **not** included here; download Supplementary Table 1
(`41586_2018_461_MOESM3_ESM.xlsx`) from the paper's page and respect its terms.

```powershell
pip install -e ".[dev]"                                    # adds openpyxl
evo-vscore-brca1 fetch-ref --out data\chr17.fa             # hg19 chr17 from UCSC (about 80 MB)
evo-vscore-brca1 prepare --xlsx data\41586_2018_461_MOESM3_ESM.xlsx --out-dir data\brca1 --sample 400
evo-vscore --vcf data\brca1\brca1.vcf --ref data\chr17.fa --backend nvidia --dry-run
evo-vscore --vcf data\brca1\brca1.vcf --ref data\chr17.fa --backend nvidia --out out\brca1
evo-vscore-brca1 evaluate --results out\brca1.csv --labels data\brca1\brca1_labels.csv --summary out\brca1_summary.md
```

Leave out `--sample 400` to use all 3,893 SNVs (about two hours at API rates). The evaluation reports
the AUROC of LOF versus the rest, with a bootstrap confidence interval. 0.5 is chance and 1.0 is
perfect. If many variants show `ref_mismatch`, the reference build or chromosome does not match the
spreadsheet (it uses hg19).

### Results (22 Oct 2026)

All 3,893 BRCA1 SNVs, Evo 2 7B through NVIDIA's hosted API, 8,192-base windows (about 2.2 hours at
36 requests a minute). LOF is the positive class.

| Comparison | n | LOF | AUROC | 95% CI (bootstrap) |
|---|---|---|---|---|
| LOF vs FUNC+INT (all) | 3,893 | 823 | 0.873 | 0.860 to 0.886 |
| LOF vs FUNC only (INT excluded) | 3,644 | 823 | 0.886 | 0.874 to 0.898 |
| Missense only | 2,086 | 441 | 0.784 | 0.760 to 0.807 |
| Splice region only | 446 | 102 | 0.910 | 0.880 to 0.938 |
| Canonical splice only | 143 | 128 | 0.691 | 0.512 to 0.839 |

Mean `delta_ll` by class: FUNC -7.1, INT -17.7, LOF -33.6 (more negative = less plausible).
An earlier 400-variant random sample gave 0.890 (0.849 to 0.921), consistent with the full run.

For context, NVIDIA's [BioNeMo tutorial](https://docs.nvidia.com/bionemo-framework/2.6/user-guide/examples/bionemo-evo2/zeroshot_brca1/)
lists about 0.87 for Evo 2 7B on the full set. The protocol is similar but not identical, so treat that
as a sanity check, not a head-to-head comparison.

What this does **not** show: it is one gene and SNVs only (indels, the other main feature of this tool,
are not covered by this benchmark), and the labels come from a cell-line assay rather than patient
outcomes. The canonical-splice group is almost all LOF (128 of 143), so its AUROC rests on only 15
non-LOF variants and is very uncertain.

## Validation on ClinVar indels

SNVs are only half the story, so there is a second benchmark for insertions and deletions: ClinVar
variants labelled pathogenic (or likely pathogenic) versus benign (or likely benign), with at least one
review star, up to 50 bases long. This follows the kind of comparison the Evo 2 paper made for non-SNVs.
Download `variant_summary.txt.gz` from NCBI's
[ClinVar FTP area](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/) and read ClinVar's
[data-use terms](https://www.ncbi.nlm.nih.gov/clinvar/intro/) first.

```powershell
evo-vscore-brca1 prepare-clinvar --summary-file data\variant_summary.txt.gz --out-dir data\clinvar --per-stratum 100
evo-vscore-brca1 fetch-ref --out data\hg19.fa --chroms 1-22,X        # about 3 GB; slow but one-off
evo-vscore --vcf data\clinvar\clinvar.vcf --ref data\hg19.fa --backend nvidia --dry-run
evo-vscore --vcf data\clinvar\clinvar.vcf --ref data\hg19.fa --backend nvidia --rpm 36 --out out\clinvar
evo-vscore-brca1 evaluate --results out\clinvar.csv --labels data\clinvar\clinvar_labels.csv --positive PATH --drop NONE --summary out\clinvar_summary.md
```

**Why the sample is matched.** In ClinVar, pathogenic and benign indels sit in different places:
a plain random sample of 300 of each had 280 of 300 pathogenic variants in coding regions and 291 of 300
benign ones outside them. A scorer that merely tells coding DNA from non-coding DNA would look good on
that, so `--per-stratum N` takes the *same number* of pathogenic and benign variants within each variant
type (coding frameshift, coding in-frame, non-coding), up to N each. Type alone then cannot predict the
label, and the evaluation reports an AUROC within each type as well as overall. Types that lack enough
variants of one class (benign frameshifts are rare) are listed and left out.

The type is a heuristic read from the variant's HGVS name and the length change (a length that is not a
multiple of 3 counts as a frameshift), not a proper annotation. The evaluation also compares the summed
`delta_ll` with the per-base `delta_mean_ll`, since insertions and deletions change the window length.

### Results (22 Oct 2026)

600 ClinVar indels, matched: 100 pathogenic and 100 benign in each of three types. GRCh37, at least one
review star, up to 50 bases, 8,192-base windows, Evo 2 7B through NVIDIA's hosted API (47 minutes).

| Comparison | n | Pathogenic | AUROC | 95% CI (bootstrap) |
|---|---|---|---|---|
| All, summed `delta_ll` | 600 | 300 | 0.826 | 0.791 to 0.861 |
| All, per-base `delta_mean_ll` | 600 | 300 | 0.821 | 0.786 to 0.854 |
| Deletions | 397 | 219 | 0.823 | 0.781 to 0.864 |
| Insertions | 203 | 81 | 0.838 | 0.783 to 0.890 |
| Coding, frameshift | 200 | 100 | 0.667 | 0.588 to 0.744 |
| Coding, in-frame | 200 | 100 | 0.946 | 0.912 to 0.976 |
| Non-coding | 200 | 100 | 0.923 | 0.877 to 0.963 |

Mean `delta_ll`: benign -20.0, pathogenic -56.3. The summed and per-base scores perform alike (the
intervals overlap), so the summed `delta_ll` stays the default.

The model separates pathogenic from benign well for in-frame coding and non-coding indels and poorly for
frameshifts. Both classes in that group are frameshifts, so they differ in where the frameshift falls
(for example whether the transcript escapes decay) rather than in the change itself; that is a plausible
reason, not something tested here.

What this does **not** show: it is a single 600-variant sample; ClinVar labels are submitters'
interpretations rather than ground truth (the benign frameshift labels in particular may be noisy); variant
type was matched but gene and repeat context were not (many benign indels sit in repeat tracts, an easy cue);
and no comparison with the published Evo 2 non-SNV numbers is made, since the protocols differ.

## Limitations (v0.2)

- Validated only on BRCA1 SNVs so far; indels and other genes are untested. Do not draw clinical conclusions.
- Only the 7B model, through one hosted API. Its availability, limits and terms are NVIDIA's.
- Large VCFs are slow at API rates: roughly 1 to 2 calls per variant at 30 requests a minute.
- Uncompressed FASTA only. No structural variants.

## Roadmap

- **v0.2 (done):** BRCA1 SNV validation tooling and a full 3,893-variant result (above).
- **v0.2.1:** ClinVar indel validation (above): overall AUROC 0.826 on a matched sample, weakest for frameshifts.
- **v0.3:** a `local` backend for people with their own GPU, a rented GPU or a cluster, with no API limit.
- Later: genome-track outputs (bigWig/BED), a pip release, and a web interface.

## Responsible use

Use this for research. The Evo 2 authors excluded eukaryotic-virus sequences from training as a
safeguard; this tool adds no capability for designing sequences, but you remain responsible for
following the licences and terms of the model and of any API you use.

## Licence

Apache-2.0, see `LICENSE`. Evo 2 and the NVIDIA API have their own licences and terms.
