# evo-vscore

Score DNA variants from a VCF file with a DNA language model. For each SNV or indel the tool asks
the model how plausible the DNA looks with the variant versus without it, and reports the
difference. It writes a CSV and a self-contained HTML report.

> **Status: v0.1, alpha.** The scores are **not yet validated** against known pathogenic and
> benign variants (that is planned for v0.2). **Research use only. Not for clinical decisions.**

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

## Limitations (v0.1)

- Scores are unvalidated. Do not draw clinical conclusions from them.
- Only the 7B model, through one hosted API. Its availability, limits and terms are NVIDIA's.
- Large VCFs are slow at API rates: roughly 1 to 2 calls per variant at 30 requests a minute.
- Uncompressed FASTA only. No structural variants.

## Roadmap

- **v0.2:** validation on a public benchmark (for example BRCA1), plus a short results write-up.
- **v0.3:** a `local` backend for people with their own GPU, a rented GPU or a cluster, with no API limit.
- Later: genome-track outputs (bigWig/BED), a pip release, and a web interface.

## Responsible use

Use this for research. The Evo 2 authors excluded eukaryotic-virus sequences from training as a
safeguard; this tool adds no capability for designing sequences, but you remain responsible for
following the licences and terms of the model and of any API you use.

## Licence

Apache-2.0, see `LICENSE`. Evo 2 and the NVIDIA API have their own licences and terms.
