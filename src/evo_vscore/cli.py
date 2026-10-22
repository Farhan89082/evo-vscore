"""Command-line interface: ``evo-vscore --vcf V --ref FASTA --backend mock|nvidia``."""

import argparse
import sys
import time
from pathlib import Path

from . import __version__
from .backends import BACKENDS, BackendError, NvidiaBackend, MockBackend, get_api_key
from .cache import ScoreCache
from .fasta import FastaError, FastaReference
from .report import CsvWriter, write_html
from .scoring import Scorer, plan
from .vcf import iter_vcf

DEFAULT_WINDOW = 8192
MIN_WINDOW = 256
MAX_TESTED_WINDOW = 32000
DEFAULT_RPM = 30.0


def default_cache_path():
    return Path.home() / ".cache" / "evo-vscore" / "scores.sqlite"


def build_parser():
    p = argparse.ArgumentParser(
        prog="evo-vscore",
        description="Score SNVs and indels from a VCF with a DNA language model. "
        "Research use only; scores are unvalidated.",
    )
    p.add_argument("--vcf", required=True, help="input VCF (plain or .gz)")
    p.add_argument("--ref", required=True, help="reference genome FASTA (uncompressed)")
    p.add_argument("--backend", required=True, choices=sorted(BACKENDS),
                   help="'nvidia' = real Evo 2 scores via NVIDIA's hosted API; "
                        "'mock' = fake numbers for testing")
    p.add_argument("--out", default="results", help="output prefix (writes PREFIX.csv and PREFIX.html)")
    p.add_argument("--window", type=int, default=DEFAULT_WINDOW,
                   help="bases of reference context scored per variant (default %(default)s)")
    p.add_argument("--rpm", type=float, default=DEFAULT_RPM,
                   help="max API requests per minute for the nvidia backend (default %(default)s)")
    p.add_argument("--cache", help="score cache file (default: ~/.cache/evo-vscore/scores.sqlite)")
    p.add_argument("--no-cache", action="store_true", help="do not read or write the on-disk cache")
    p.add_argument("--limit", type=int, help="only process the first N variants (handy for a trial run)")
    p.add_argument("--dry-run", action="store_true",
                   help="count variants and API calls, estimate the time, score nothing")
    p.add_argument("--version", action="version", version="evo-vscore " + __version__)
    return p


def _log(message):
    print(message, file=sys.stderr, flush=True)


def _make_backend(args):
    if args.backend == "nvidia":
        return NvidiaBackend(get_api_key(), rpm=args.rpm)
    return BACKENDS[args.backend]()


def _print_plan(info, args, cache_id):
    _log("Variants read:            %d" % info["total"])
    _log("Scoreable:                %d" % info["scoreable"])
    for status, count in info["statuses"].most_common():
        if status != "ok":
            _log("  not scored: %-36s %d" % (status, count))
    _log("Unique windows/sequences: %d (%d already cached)" % (info["unique_sequences"], info["cached"]))
    _log("Backend calls needed:     %d" % info["calls_needed"])
    if args.backend == "nvidia":
        minutes = info["calls_needed"] / args.rpm
        _log("Estimated time at %g requests/min: about %.1f minutes" % (args.rpm, minutes))
    else:
        _log("The mock backend makes no API calls.")


def run(args):
    if not (MIN_WINDOW <= args.window <= 2 * MAX_TESTED_WINDOW):
        raise ValueError("--window must be between %d and %d" % (MIN_WINDOW, 2 * MAX_TESTED_WINDOW))
    if args.window > MAX_TESTED_WINDOW:
        _log("warning: windows above %d bases have not been tested with the NVIDIA API" % MAX_TESTED_WINDOW)
    if args.rpm <= 0:
        raise ValueError("--rpm must be positive")
    if not Path(args.vcf).exists():
        raise FileNotFoundError("VCF not found: %s" % args.vcf)

    variants = list(iter_vcf(args.vcf))
    if args.limit is not None:
        variants = variants[: max(args.limit, 0)]
    cache = ScoreCache(None if args.no_cache else (args.cache or default_cache_path()))
    backend_cls = BACKENDS[args.backend]

    with FastaReference(args.ref) as fasta:
        if args.dry_run:
            _print_plan(plan(variants, fasta, args.window, cache, backend_cls.cache_id), args, backend_cls.cache_id)
            return 0

        backend = _make_backend(args)
        if backend.is_mock:
            _log("NOTE: mock backend - the scores are fake and have no biological meaning.")
        scorer = Scorer(fasta, backend, cache, args.window)
        csv_path, html_path = args.out + ".csv", args.out + ".html"
        writer = CsvWriter(csv_path)
        results, total, started = [], len(variants), time.time()
        step = 1 if total <= 100 else max(1, total // 50)
        problem, code = None, 0
        try:
            for index, variant in enumerate(variants, start=1):
                result = scorer.score(variant)
                results.append(result)
                writer.write(result)
                if index % step == 0 or index == total:
                    _log("[%d/%d] %s:%d %s>%s  %s" % (index, total, variant.chrom, variant.pos,
                                                      variant.ref[:10], variant.alt[:10], result.status))
        except KeyboardInterrupt:
            problem, code = "Interrupted.", 130
        except BackendError as error:
            problem, code = "Stopped: %s" % error, 2
        finally:
            writer.close()
            backend.close()
            cache.close()
        if problem and not results:
            _log(problem)  # nothing was scored, so there are no results or resume hint to report
            return code
        write_html(results, html_path, backend.name, backend.is_mock, args.window, partial=problem is not None)

    scored = sum(1 for r in results if r.status == "ok")
    _log("%s in %.1fs: %d scored, %d skipped. Backend calls: %d, cache hits: %d."
         % ("Stopped early" if problem else "Done", time.time() - started, scored,
            len(results) - scored, scorer.backend_calls, scorer.cache_hits))
    _log("Wrote %s and %s" % (csv_path, html_path))
    if problem:
        _log(problem)
        if scorer.backend_calls:
            _log("Re-run the same command to resume; finished scores are cached.")
    return code


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        return run(args)
    except (BackendError, FastaError, ValueError, FileNotFoundError) as error:
        _log("error: %s" % error)
        return 2


if __name__ == "__main__":
    sys.exit(main())
