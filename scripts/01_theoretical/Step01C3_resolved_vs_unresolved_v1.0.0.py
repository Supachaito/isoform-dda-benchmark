#!/usr/bin/env python3
"""
Step01C3_resolved_vs_unresolved_v1.0.0.py

Characterize the strongest Step01C.2 set:

- shared_all_strict accessibility
- accessible in all 4 workflows
- theoretically primary-isoform-resolvable
- compare families resolved by >=1 workflow vs unresolved by all 4 workflows

The purpose is to test whether unresolved families simply have fewer theoretical
isoform-discriminative opportunities or weaker family-accessibility evidence.

No search engine, FASTA remapping, or classification is re-run.
"""

from __future__ import annotations
import argparse
import csv
import json
import math
import random
import statistics
from pathlib import Path
from typing import List, Dict, Optional

SCRIPT_VERSION = "1.0.0"
EXPECTED = {
    "all4_theoretically_resolvable": 59,
    "resolved_any": 7,
    "unresolved_all4": 52,
}

METRICS = [
    ("n_reference_entries", "Reference entries per family"),
    ("theoretical_primary_peptides", "Theoretical primary isoform-discriminative peptides"),
    ("min_access_evidence_peptides_across_supporting_workflows",
     "Minimum shared-all accessibility peptides across the 4 workflows"),
    ("total_access_evidence_peptide_workflow_incidences",
     "Total shared-all peptide-workflow incidences"),
]


def intv(x):
    try:
        return int(float(str(x).strip()))
    except Exception:
        return 0


def floatv(x):
    try:
        return float(str(x).strip())
    except Exception:
        return math.nan


def write_tsv(path: Path, rows: List[dict], fields: Optional[List[str]] = None):
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    if fields is None:
        fields = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def median_iqr(vals):
    vals = sorted(vals)
    if not vals:
        return math.nan, math.nan, math.nan
    med = statistics.median(vals)
    # Inclusive quartiles are stable for small n and easy to reproduce.
    if len(vals) == 1:
        return med, vals[0], vals[0]
    q = statistics.quantiles(vals, n=4, method="inclusive")
    return med, q[0], q[2]


def cliffs_delta(a, b):
    """Cliff's delta: P(a>b) - P(a<b)."""
    if not a or not b:
        return math.nan
    gt = lt = 0
    for x in a:
        for y in b:
            if x > y:
                gt += 1
            elif x < y:
                lt += 1
    return (gt - lt) / (len(a) * len(b))


def bootstrap_median_diff(a, b, nboot=10000, seed=20260812):
    """Median(resolved) - median(unresolved), percentile bootstrap CI."""
    if not a or not b:
        return math.nan, math.nan, math.nan
    rng = random.Random(seed)
    diffs = []
    for _ in range(nboot):
        aa = [rng.choice(a) for _ in range(len(a))]
        bb = [rng.choice(b) for _ in range(len(b))]
        diffs.append(statistics.median(aa) - statistics.median(bb))
    diffs.sort()
    lo = diffs[int(0.025 * (nboot - 1))]
    hi = diffs[int(0.975 * (nboot - 1))]
    obs = statistics.median(a) - statistics.median(b)
    return obs, lo, hi


def mann_whitney_if_available(a, b):
    try:
        from scipy.stats import mannwhitneyu
        res = mannwhitneyu(a, b, alternative="two-sided", method="auto")
        return float(res.statistic), float(res.pvalue), "scipy.stats.mannwhitneyu"
    except Exception:
        return math.nan, math.nan, "not_available"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--step01c2-dir", required=True)
    ap.add_argument("--outdir", default="Step01C3_resolved_vs_unresolved_results")
    args = ap.parse_args()

    indir = Path(args.step01c2_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    inp = indir / "consensus_accessibility_family_long.tsv"
    if not inp.exists():
        raise SystemExit(f"Missing input: {inp}")

    with inp.open("r", encoding="utf-8-sig", errors="replace", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    strongest = [
        r for r in rows
        if r.get("accessibility_definition") == "shared_all_strict"
        and intv(r.get("accessible_workflow_count")) == 4
        and intv(r.get("theoretically_primary_resolvable")) == 1
    ]

    resolved = [r for r in strongest if intv(r.get("resolved_by_any_workflow")) == 1]
    unresolved = [r for r in strongest if intv(r.get("unresolved_by_all_workflows")) == 1]

    validation = []
    for metric, observed, expected in [
        ("all4_theoretically_resolvable", len(strongest), EXPECTED["all4_theoretically_resolvable"]),
        ("resolved_any", len(resolved), EXPECTED["resolved_any"]),
        ("unresolved_all4", len(unresolved), EXPECTED["unresolved_all4"]),
    ]:
        validation.append({
            "metric": metric,
            "observed": observed,
            "expected": expected,
            "difference": observed - expected,
            "status": "MATCH" if observed == expected else "CHECK",
        })
    write_tsv(outdir / "Step01C3_validation.tsv", validation)

    if any(r["status"] != "MATCH" for r in validation):
        raise RuntimeError("Validation failed; do not interpret Step01C3.")

    summary = []
    for col, label in METRICS:
        a = [floatv(r.get(col)) for r in resolved]
        b = [floatv(r.get(col)) for r in unresolved]
        a = [x for x in a if not math.isnan(x)]
        b = [x for x in b if not math.isnan(x)]

        amed, aq1, aq3 = median_iqr(a)
        bmed, bq1, bq3 = median_iqr(b)
        d = cliffs_delta(a, b)
        mdiff, ci_lo, ci_hi = bootstrap_median_diff(a, b)
        u, p, test_source = mann_whitney_if_available(a, b)

        summary.append({
            "metric": col,
            "metric_label": label,
            "resolved_n": len(a),
            "resolved_median": amed,
            "resolved_q1": aq1,
            "resolved_q3": aq3,
            "unresolved_n": len(b),
            "unresolved_median": bmed,
            "unresolved_q1": bq1,
            "unresolved_q3": bq3,
            "median_difference_resolved_minus_unresolved": mdiff,
            "bootstrap95CI_low": ci_lo,
            "bootstrap95CI_high": ci_hi,
            "cliffs_delta_resolved_vs_unresolved": d,
            "mann_whitney_U": u,
            "mann_whitney_p_two_sided": p,
            "statistical_test_source": test_source,
        })

    write_tsv(outdir / "resolved_vs_unresolved_summary.tsv", summary)

    # Family-level table with a compact status.
    family_rows = []
    for r in strongest:
        x = dict(r)
        x["group"] = "resolved" if intv(r.get("resolved_by_any_workflow")) == 1 else "unresolved"
        family_rows.append(x)

    family_rows.sort(
        key=lambda r: (
            0 if r["group"] == "resolved" else 1,
            -intv(r.get("theoretical_primary_peptides")),
            -intv(r.get("min_access_evidence_peptides_across_supporting_workflows")),
            -intv(r.get("total_access_evidence_peptide_workflow_incidences")),
            str(r.get("family")),
        )
    )
    write_tsv(outdir / "all4_strict_theoretical_59_family_characteristics.tsv", family_rows)

    # Rank the 52 unresolved families: strongest evidence first.
    ranked_unresolved = sorted(
        unresolved,
        key=lambda r: (
            -intv(r.get("min_access_evidence_peptides_across_supporting_workflows")),
            -intv(r.get("total_access_evidence_peptide_workflow_incidences")),
            -intv(r.get("theoretical_primary_peptides")),
            -intv(r.get("n_reference_entries")),
            str(r.get("family")),
        ),
    )
    rank_rows = []
    for i, r in enumerate(ranked_unresolved, 1):
        rank_rows.append({
            "rank": i,
            "family": r.get("family", ""),
            "n_reference_entries": intv(r.get("n_reference_entries")),
            "theoretical_primary_peptides": intv(r.get("theoretical_primary_peptides")),
            "min_shared_all_peptides_across_4_workflows":
                intv(r.get("min_access_evidence_peptides_across_supporting_workflows")),
            "total_shared_all_peptide_workflow_incidences":
                intv(r.get("total_access_evidence_peptide_workflow_incidences")),
            "observed_primary_union": intv(r.get("observed_primary_union")),
        })
    write_tsv(outdir / "ranked_52_high_confidence_unresolved_families.tsv", rank_rows)

    # Descriptive headline.
    resolution_support = {}
    for n in range(1, 5):
        resolution_support[str(n)] = sum(
            1 for r in resolved
            if intv(r.get("primary_resolution_workflow_count")) == n
        )

    headline = {
        "n_all4_accessible_theoretically_resolvable": len(strongest),
        "n_resolved_by_any_workflow": len(resolved),
        "n_unresolved_by_all_workflows": len(unresolved),
        "unresolved_pct": 100 * len(unresolved) / len(strongest),
        "n_resolved_in_all4_workflows": resolution_support["4"],
        "resolved_in_all4_workflows_pct_of_59":
            100 * resolution_support["4"] / len(strongest),
        "resolution_support_distribution_among_7_resolved": resolution_support,
    }
    (outdir / "Step01C3_headline.json").write_text(
        json.dumps(headline, indent=2), encoding="utf-8"
    )

    # Figures.
    try:
        import matplotlib.pyplot as plt

        # One figure per metric; no subplots.
        for col, label in METRICS:
            a = [floatv(r.get(col)) for r in resolved]
            b = [floatv(r.get(col)) for r in unresolved]
            a = [x for x in a if not math.isnan(x)]
            b = [x for x in b if not math.isnan(x)]

            fig, ax = plt.subplots(figsize=(5.8, 4.8))
            ax.boxplot([a, b], labels=["Resolved\n(n=7)", "Unresolved\n(n=52)"],
                       showfliers=True)
            # jittered points, deterministic
            rng = random.Random(20260812)
            for xpos, vals in [(1, a), (2, b)]:
                xs = [xpos + rng.uniform(-0.06, 0.06) for _ in vals]
                ax.scatter(xs, vals, alpha=0.7, s=20)
            ax.set_ylabel(label)
            ax.set_title(label)
            ax.grid(axis="y", alpha=0.25)
            fig.tight_layout()
            safe = col.replace("/", "_")
            fig.savefig(outdir / f"Fig_C3_{safe}.png", dpi=300)
            fig.savefig(outdir / f"Fig_C3_{safe}.pdf")
            plt.close(fig)

    except Exception as exc:
        print(f"[warning] Figure generation skipped: {exc}", file=sys.stderr)

    print("[validation] all anchors MATCH")
    print(f"[cohort] 59 strict all-4 accessible + theoretically resolvable families")
    print(f"         resolved by any workflow: {len(resolved)}")
    print(f"         unresolved by all 4:      {len(unresolved)} ({100*len(unresolved)/len(strongest):.2f}%)")
    print(f"         resolved in all 4:        {resolution_support['4']} ({100*resolution_support['4']/len(strongest):.2f}% of 59)")
    print("\nResolved vs unresolved medians:")
    for r in summary:
        ptxt = "NA" if math.isnan(float(r["mann_whitney_p_two_sided"])) else f"{float(r['mann_whitney_p_two_sided']):.4g}"
        print(
            f"  {r['metric']}: "
            f"{r['resolved_median']} vs {r['unresolved_median']}; "
            f"Cliff's delta={float(r['cliffs_delta_resolved_vs_unresolved']):.3f}; p={ptxt}"
        )
    print(f"\nDone. Outputs: {outdir.resolve()}")


if __name__ == "__main__":
    main()
