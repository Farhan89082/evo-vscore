# Security and privacy

## Reporting a problem

Please report security problems **privately**, not in a public issue. Use the
"Report a vulnerability" button on the repository's **Security** tab (GitHub private
vulnerability reporting). Include what you found, how to reproduce it, and the version
(`evo-vscore --version`). This is a small volunteer project, so replies may take a few days.

## API keys

- The tool reads your NVIDIA API key from the `NVIDIA_API_KEY` (or `NVCF_KEY`) environment
  variable. It never prints, logs or stores the key, and does not accept it as a command-line
  argument (which would end up in your shell history).
- Never commit a key to Git. The included `.gitignore` excludes common key and data files, but
  check `git status` and `git diff --cached` before every commit.
- If a key may have been exposed, revoke it on build.nvidia.com and create a new one.

## What leaves your computer

- `--backend nvidia` sends DNA sequences (a window of reference genome with the variant applied)
  to NVIDIA's servers. **Do not use it with identifiable patient or other health data.** See
  "Data privacy and NVIDIA's terms" in the README.
- `--backend mock` sends nothing anywhere.
- The local score cache (`~/.cache/evo-vscore/scores.sqlite`) stores SHA-256 hashes of sequences
  and their scores, not the sequences themselves. Result files (`.csv`, `.html`) contain your
  variants, so handle them with the same care as the VCF they came from.

## Supported versions

Only the latest release receives fixes.
