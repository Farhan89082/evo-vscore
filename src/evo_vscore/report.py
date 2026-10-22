"""Write results as CSV (incrementally) and as a single self-contained HTML page."""

import csv
import dataclasses
import html
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from . import __version__
from .scoring import ScoredVariant

COLUMNS = [f.name for f in dataclasses.fields(ScoredVariant)]
MAX_HTML_ROWS = 2000


def _fmt(value):
    if value is None:
        return ""
    if isinstance(value, float):
        return "%.4f" % value
    return value


class CsvWriter:
    """Writes one row at a time and flushes, so a crash keeps the finished rows."""

    def __init__(self, path):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        self._fh = open(path, "w", newline="", encoding="utf-8")
        self._writer = csv.writer(self._fh)
        self._writer.writerow(COLUMNS)

    def write(self, result):
        row = dataclasses.asdict(result)
        self._writer.writerow([_fmt(row[c]) for c in COLUMNS])
        self._fh.flush()

    def close(self):
        self._fh.close()


_CSS = """
:root{--bg:#fff;--fg:#1d2328;--muted:#5b6670;--line:#dde2e6;--card:#f6f8fa;--bar:#c2410c;--warn-bg:#fff4d6;--warn-fg:#5c4400}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#14181c;--fg:#e6eaee;--muted:#9aa5af;--line:#2b333a;--card:#1b2026;--bar:#fb923c;--warn-bg:#3a2f0b;--warn-fg:#f3dd9a}}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,Segoe UI,sans-serif}
main{max-width:1000px;margin:0 auto;padding:24px 16px 48px}
h1{font-size:22px;margin:0 0 4px}h2{font-size:17px;margin:28px 0 8px}
.muted{color:var(--muted)}
.banner{background:var(--warn-bg);color:var(--warn-fg);border-radius:8px;padding:10px 14px;margin:14px 0}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:12px 14px}
table{border-collapse:collapse;width:100%;font-size:14px}
th,td{text-align:left;padding:6px 8px;border-bottom:1px solid var(--line);vertical-align:middle}
th{color:var(--muted);font-weight:600}
td.num{text-align:right;font-variant-numeric:tabular-nums}
.bar{display:inline-block;height:8px;background:var(--bar);border-radius:4px;vertical-align:middle}
.scroll{overflow-x:auto}
code{font-family:ui-monospace,Consolas,monospace;font-size:13px}
"""


def write_html(results, path, backend_name, is_mock, window, partial=False):
    esc = html.escape
    scored = [r for r in results if r.status == "ok"]
    statuses = Counter(r.status for r in results if r.status != "ok")
    scored.sort(key=lambda r: r.delta_ll)
    shown = scored[:MAX_HTML_ROWS]
    biggest = max((-r.delta_ll for r in shown if r.delta_ll < 0), default=0.0) or 1.0

    parts = []
    parts.append("<!doctype html><html lang='en'><head><meta charset='utf-8'>")
    parts.append("<meta name='viewport' content='width=device-width,initial-scale=1'>")
    parts.append("<title>evo-vscore results</title><style>%s</style></head><body><main>" % _CSS)
    parts.append("<h1>evo-vscore results</h1>")
    parts.append(
        "<div class='muted'>evo-vscore %s &middot; backend <code>%s</code> &middot; window %d bases &middot; %s UTC</div>"
        % (esc(__version__), esc(backend_name), window, datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"))
    )
    if is_mock:
        parts.append(
            "<div class='banner'><strong>Mock backend.</strong> These numbers are fake and carry no "
            "biological meaning. Use <code>--backend nvidia</code> for real scores.</div>"
        )
    if partial:
        parts.append(
            "<div class='banner'><strong>Partial results.</strong> The run stopped early. "
            "Re-run the same command to resume; finished scores are cached.</div>"
        )
    parts.append(
        "<div class='card' style='margin-top:14px'><strong>%d</strong> variants read &middot; "
        "<strong>%d</strong> scored &middot; <strong>%d</strong> skipped</div>"
        % (len(results), len(scored), len(results) - len(scored))
    )
    parts.append(
        "<h2>How to read this</h2><div class='card'>"
        "<code>delta_ll</code> is the model's log-likelihood of the window with the variant minus the "
        "reference window. <strong>More negative = the model finds the variant less plausible</strong>, "
        "which is a hint of disruption, not a diagnosis. For insertions and deletions the two windows have "
        "different lengths, so <code>delta_mean_ll</code> (per base) is shown as well. Scores are "
        "<strong>research output, validated so far only on BRCA1 SNVs</strong>; do not use them for clinical decisions.</div>"
    )
    parts.append("<h2>Most disruptive variants</h2>")
    if not shown:
        parts.append("<p class='muted'>No variants were scored.</p>")
    else:
        parts.append("<div class='scroll'><table><thead><tr><th>Variant</th><th>ID</th><th>Type</th>"
                     "<th class='num'>delta_ll</th><th></th><th class='num'>delta_mean_ll</th></tr></thead><tbody>")
        for r in shown:
            # bars only for variants the model finds LESS plausible (negative delta)
            width = max(2, int(100 * -r.delta_ll / biggest)) if r.delta_ll < 0 else 0
            label = "%s:%d %s&gt;%s" % (esc(r.chrom), r.pos, esc(_short(r.ref)), esc(_short(r.alt)))
            parts.append(
                "<tr><td><code>%s</code></td><td>%s</td><td>%s</td><td class='num'>%.3f</td>"
                "<td><span class='bar' style='width:%dpx'></span></td><td class='num'>%.5f</td></tr>"
                % (label, esc(r.id), esc(r.variant_type), r.delta_ll, width, r.delta_mean_ll)
            )
        parts.append("</tbody></table></div>")
        if len(scored) > len(shown):
            parts.append("<p class='muted'>Showing the %d most negative of %d scored variants. "
                         "The CSV has all of them.</p>" % (len(shown), len(scored)))
    if statuses:
        parts.append("<h2>Skipped variants</h2><div class='card'><table><tbody>")
        for status, count in statuses.most_common():
            parts.append("<tr><td><code>%s</code></td><td class='num'>%d</td></tr>" % (esc(status), count))
        parts.append("</tbody></table></div>")
    parts.append("</main></body></html>")

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("".join(parts), encoding="utf-8")


def _short(allele, limit=12):
    return allele if len(allele) <= limit else allele[:limit] + "..(%d)" % len(allele)
