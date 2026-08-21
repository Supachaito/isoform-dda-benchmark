#!/usr/bin/env python3
"""
Step01C2_consensus_accessibility_v1.0.0.py

Consensus accessibility across AP / FP / MM / MQ.

This analysis is downstream of the frozen Step01B v1.1.1 and Step01C v1.0.0
outputs. It does NOT alter peptide extraction, FASTA mapping, isoform-family
classification, or the theoretical tryptic catalog.

Question
--------
If a multi-entry protein family is independently confirmed as accessible by
multiple workflows, does it remain empirically isoform-unresolved?

For each Step01C accessibility definition:
  - broad_any
  - family_specific_nonprimary
  - shared_all_strict

families are stratified by accessibility support in:
  >=1, >=2, >=3, or all 4 workflows.

The main conservative endpoint is:
  family is accessible in >=N workflows
  AND theoretically primary-isoform-resolvable by trypsin
  AND NO primary isoform-discriminative peptide is observed by ANY of the four
  workflows.

Using "resolved by any workflow" makes the unresolved endpoint conservative:
a family is called unresolved only if all four pipelines fail to supply primary
isoform-discriminative evidence despite consensus family accessibility.

Additional outputs count how many workflows provide primary isoform evidence
(0,1,2,3,4), and count non-primary accessibility-evidence peptide keys for the
family-specific and shared-all definitions.

No search engine is re-run.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

SCRIPT_VERSION = "1.0.0"
WORKFLOWS = ["AP", "FP", "MM", "MQ"]
DEFINITIONS = ["broad_any", "family_specific_nonprimary", "shared_all_strict"]

# Frozen Step01C v1.0.0 pooled anchors already reproduced in the current analysis.
EXPECTED_ANCHORS = {
    "n_multi_entry_families": 10587,
    "broad_any": {
        "accessible_ge1": 4647,
        "theoretically_resolvable_ge1": 3198,
        "resolved_any_ge1": 199,
    },
    "family_specific_nonprimary": {
        "accessible_ge1": 347,
        "theoretically_resolvable_ge1": 269,
        "resolved_any_ge1": 25,
    },
    "shared_all_strict": {
        "accessible_ge1": 245,
        "theoretically_resolvable_ge1": 195,
        "resolved_any_ge1": 21,
    },
}


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def write_tsv(path: Path, rows: List[dict], fieldnames: Optional[List[str]] = None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    if fieldnames is None:
        fieldnames = list(rows[0])
    with path.open("wt", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rows)


def intv(v) -> int:
    try:
        return int(float(str(v).strip()))
    except Exception:
        return 0


def load_family_matrix(step01c_dir: Path) -> List[dict]:
    path = step01c_dir / "accessible_family_status_matrix.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Step01C family matrix not found: {path}")

    with path.open("rt", encoding="utf-8-sig", errors="replace", newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        fields = set(r.fieldnames or [])
        required = {
            "family",
            "n_reference_entries",
            "theoretical_primary_peptides",
            "theoretically_primary_resolvable",
            "observed_primary_union",
        }
        for wf in WORKFLOWS:
            required.add(f"observed_primary_{wf}")
        for d in DEFINITIONS:
            for wf in WORKFLOWS:
                required.add(f"accessible__{d}__{wf}")

        missing = required - fields
        if missing:
            raise ValueError(
                f"{path} is missing required columns: {sorted(missing)}"
            )
        return list(r)


def load_observed_classification(step01b_dir: Path) -> List[dict]:
    path = step01b_dir / "observed_full_FASTA_classification.tsv"
    if not path.exists():
        raise FileNotFoundError(
            f"Step01B full-FASTA classification not found: {path}"
        )
    with path.open("rt", encoding="utf-8-sig", errors="replace", newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        fields = set(r.fieldnames or [])
        required = {
            "peptide_key_IL_equivalent",
            "category",
            "primary_isoform_discriminative",
            "target_base_accession",
            *WORKFLOWS,
        }
        missing = required - fields
        if missing:
            raise ValueError(
                f"{path} is missing required columns: {sorted(missing)}"
            )
        return list(r)


def build_evidence_counts(obs_rows: List[dict]):
    """
    Count peptide-key evidence by family/workflow.

    broad_any cannot be reconstructed from Step01B classification alone because
    cross-family shared peptides do not carry the complete set of mapped
    families in this exported table. Binary broad accessibility therefore comes
    from the frozen Step01C family matrix.

    family_specific_nonprimary:
      any non-primary peptide with one target_base_accession.

    shared_all_strict:
      within_family_shared_all peptide.
    """
    fs_counts: Dict[str, Counter] = defaultdict(Counter)
    strict_counts: Dict[str, Counter] = defaultdict(Counter)
    primary_counts: Dict[str, Counter] = defaultdict(Counter)

    fs_keys: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    strict_keys: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    primary_keys: Dict[Tuple[str, str], Set[str]] = defaultdict(set)

    for row in obs_rows:
        key = str(row.get("peptide_key_IL_equivalent", "")).strip()
        fam = str(row.get("target_base_accession", "")).strip()
        cat = str(row.get("category", "")).strip()
        primary = str(row.get("primary_isoform_discriminative", "0")).strip() == "1"

        if not key or not fam:
            continue

        for wf in WORKFLOWS:
            if intv(row.get(wf, 0)) != 1:
                continue

            if primary:
                primary_keys[(fam, wf)].add(key)

            # This mirrors Step01C family_specific_nonprimary: target family is
            # defined and peptide is not primary. In the frozen Step01B
            # classification, a target family is assigned only for one-family
            # mappings.
            if not primary:
                fs_keys[(fam, wf)].add(key)

            if cat == "within_family_shared_all":
                strict_keys[(fam, wf)].add(key)

    for (fam, wf), keys in fs_keys.items():
        fs_counts[fam][wf] = len(keys)
    for (fam, wf), keys in strict_keys.items():
        strict_counts[fam][wf] = len(keys)
    for (fam, wf), keys in primary_keys.items():
        primary_counts[fam][wf] = len(keys)

    return fs_counts, strict_counts, primary_counts


def support_pattern(row: dict, definition: str) -> Tuple[int, str, List[str]]:
    supported = [
        wf for wf in WORKFLOWS
        if intv(row.get(f"accessible__{definition}__{wf}", 0)) == 1
    ]
    return len(supported), ",".join(supported), supported


def resolution_pattern(row: dict) -> Tuple[int, str, List[str]]:
    supported = [
        wf for wf in WORKFLOWS
        if intv(row.get(f"observed_primary_{wf}", 0)) > 0
    ]
    return len(supported), ",".join(supported), supported


def make_family_long(
    family_rows: List[dict],
    fs_counts,
    strict_counts,
    primary_counts,
) -> List[dict]:
    out = []

    for r in family_rows:
        fam = str(r["family"])
        resolution_n, resolution_list, _ = resolution_pattern(r)

        for d in DEFINITIONS:
            access_n, access_list, access_wfs = support_pattern(r, d)

            row = {
                "accessibility_definition": d,
                "family": fam,
                "n_reference_entries": intv(r["n_reference_entries"]),
                "theoretical_primary_peptides": intv(r["theoretical_primary_peptides"]),
                "theoretically_primary_resolvable": intv(r["theoretically_primary_resolvable"]),
                "accessible_workflow_count": access_n,
                "accessible_workflows": access_list,
                "primary_resolution_workflow_count": resolution_n,
                "primary_resolution_workflows": resolution_list,
                "observed_primary_union": intv(r["observed_primary_union"]),
                "resolved_by_any_workflow": int(resolution_n >= 1),
                "unresolved_by_all_workflows": int(resolution_n == 0),
            }

            for wf in WORKFLOWS:
                row[f"accessible_{wf}"] = intv(r.get(f"accessible__{d}__{wf}", 0))
                row[f"observed_primary_{wf}"] = intv(r.get(f"observed_primary_{wf}", 0))

                if d == "family_specific_nonprimary":
                    row[f"access_evidence_peptides_{wf}"] = fs_counts[fam][wf]
                elif d == "shared_all_strict":
                    row[f"access_evidence_peptides_{wf}"] = strict_counts[fam][wf]
                else:
                    # broad_any binary support is valid, but peptide-key count
                    # is not reconstructable without the original many-to-many
                    # mappings, so leave blank rather than inventing a number.
                    row[f"access_evidence_peptides_{wf}"] = ""

                row[f"primary_peptide_keys_{wf}"] = primary_counts[fam][wf]

            if d == "family_specific_nonprimary":
                vals = [fs_counts[fam][wf] for wf in access_wfs]
                row["min_access_evidence_peptides_across_supporting_workflows"] = (
                    min(vals) if vals else 0
                )
                row["total_access_evidence_peptide_workflow_incidences"] = sum(
                    fs_counts[fam][wf] for wf in WORKFLOWS
                )
            elif d == "shared_all_strict":
                vals = [strict_counts[fam][wf] for wf in access_wfs]
                row["min_access_evidence_peptides_across_supporting_workflows"] = (
                    min(vals) if vals else 0
                )
                row["total_access_evidence_peptide_workflow_incidences"] = sum(
                    strict_counts[fam][wf] for wf in WORKFLOWS
                )
            else:
                row["min_access_evidence_peptides_across_supporting_workflows"] = ""
                row["total_access_evidence_peptide_workflow_incidences"] = ""

            out.append(row)

    return out


def make_summary(family_long: List[dict]) -> List[dict]:
    rows = []

    by_def: Dict[str, List[dict]] = {
        d: [r for r in family_long if r["accessibility_definition"] == d]
        for d in DEFINITIONS
    }

    for d in DEFINITIONS:
        drows = by_def[d]
        for threshold in [1, 2, 3, 4]:
            consensus = [
                r for r in drows
                if intv(r["accessible_workflow_count"]) >= threshold
            ]
            theory = [
                r for r in consensus
                if intv(r["theoretically_primary_resolvable"]) == 1
            ]
            resolved = [
                r for r in theory
                if intv(r["resolved_by_any_workflow"]) == 1
            ]
            unresolved = [
                r for r in theory
                if intv(r["unresolved_by_all_workflows"]) == 1
            ]

            resolution_support_counts = Counter(
                intv(r["primary_resolution_workflow_count"]) for r in theory
            )

            rows.append({
                "accessibility_definition": d,
                "minimum_accessible_workflows": threshold,
                "consensus_accessible_multi_entry_families": len(consensus),
                "consensus_accessible_theoretically_resolvable_families": len(theory),
                "theoretically_resolvable_pct_of_consensus_accessible":
                    100.0 * len(theory) / len(consensus) if consensus else math.nan,
                "resolved_by_any_workflow": len(resolved),
                "unresolved_by_all_workflows": len(unresolved),
                "resolved_pct_of_consensus_accessible_theoretical":
                    100.0 * len(resolved) / len(theory) if theory else math.nan,
                "unresolved_pct_of_consensus_accessible_theoretical":
                    100.0 * len(unresolved) / len(theory) if theory else math.nan,
                "resolved_in_1_workflow": resolution_support_counts.get(1, 0),
                "resolved_in_2_workflows": resolution_support_counts.get(2, 0),
                "resolved_in_3_workflows": resolution_support_counts.get(3, 0),
                "resolved_in_4_workflows": resolution_support_counts.get(4, 0),
                "total_theoretical_primary_peptides_in_consensus_accessible_families":
                    sum(intv(r["theoretical_primary_peptides"]) for r in theory),
                "total_observed_primary_unique_family_peptide_counts":
                    sum(intv(r["observed_primary_union"]) for r in theory),
            })

    return rows


def validate_inputs(family_rows: List[dict], summary_rows: List[dict]) -> List[dict]:
    rows = []

    def add(metric, observed, expected):
        rows.append({
            "metric": metric,
            "observed": observed,
            "expected": expected,
            "difference": observed - expected,
            "status": "MATCH" if observed == expected else "CHECK",
        })

    add(
        "multi_entry_family_rows",
        len(family_rows),
        EXPECTED_ANCHORS["n_multi_entry_families"],
    )

    for d in DEFINITIONS:
        ge1 = next(
            r for r in summary_rows
            if r["accessibility_definition"] == d
            and intv(r["minimum_accessible_workflows"]) == 1
        )
        exp = EXPECTED_ANCHORS[d]
        add(
            f"{d}__accessible_ge1",
            intv(ge1["consensus_accessible_multi_entry_families"]),
            exp["accessible_ge1"],
        )
        add(
            f"{d}__theoretically_resolvable_ge1",
            intv(ge1["consensus_accessible_theoretically_resolvable_families"]),
            exp["theoretically_resolvable_ge1"],
        )
        add(
            f"{d}__resolved_any_ge1",
            intv(ge1["resolved_by_any_workflow"]),
            exp["resolved_any_ge1"],
        )

    return rows


def write_high_confidence_tables(
    outdir: Path,
    family_long: List[dict],
):
    strict = [
        r for r in family_long
        if r["accessibility_definition"] == "shared_all_strict"
        and intv(r["theoretically_primary_resolvable"]) == 1
        and intv(r["unresolved_by_all_workflows"]) == 1
    ]

    for threshold in [2, 3, 4]:
        rows = [
            r for r in strict
            if intv(r["accessible_workflow_count"]) >= threshold
        ]
        rows.sort(
            key=lambda r: (
                -intv(r["accessible_workflow_count"]),
                -intv(r["min_access_evidence_peptides_across_supporting_workflows"]),
                -intv(r["total_access_evidence_peptide_workflow_incidences"]),
                -intv(r["theoretical_primary_peptides"]),
                str(r["family"]),
            )
        )
        write_tsv(
            outdir / f"strict_accessible_ge{threshold}workflows_but_unresolved_all4.tsv",
            rows,
        )

    # The strongest set: family accessible by shared-all evidence in all 4
    # workflows, theoretically resolvable, yet no primary peptide in any workflow.
    strongest = [
        r for r in strict
        if intv(r["accessible_workflow_count"]) == 4
    ]
    strongest.sort(
        key=lambda r: (
            -intv(r["min_access_evidence_peptides_across_supporting_workflows"]),
            -intv(r["total_access_evidence_peptide_workflow_incidences"]),
            -intv(r["theoretical_primary_peptides"]),
            str(r["family"]),
        )
    )
    write_tsv(
        outdir / "strict_accessible_all4_but_isoform_unresolved.tsv",
        strongest,
    )


def make_figures(outdir: Path, summary_rows: List[dict]):
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception as exc:
        eprint(f"[warning] matplotlib unavailable; figures skipped: {exc}")
        return

    # Figure 1: unresolved percentage by consensus threshold for all definitions.
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    for d in DEFINITIONS:
        rows = sorted(
            [r for r in summary_rows if r["accessibility_definition"] == d],
            key=lambda r: intv(r["minimum_accessible_workflows"]),
        )
        x = [intv(r["minimum_accessible_workflows"]) for r in rows]
        y = [float(r["unresolved_pct_of_consensus_accessible_theoretical"]) for r in rows]
        label = {
            "broad_any": "Broad any-peptide",
            "family_specific_nonprimary": "Family-specific non-primary",
            "shared_all_strict": "Shared-all strict",
        }[d]
        ax.plot(x, y, marker="o", label=label)
    ax.set_xticks([1, 2, 3, 4])
    ax.set_xlabel("Minimum number of workflows confirming family accessibility")
    ax.set_ylabel("Empirically unresolved families (%)")
    ax.set_title("Isoform non-resolution persists with increasing accessibility consensus")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "Fig_C2_unresolved_pct_by_consensus_threshold.png", dpi=300)
    fig.savefig(outdir / "Fig_C2_unresolved_pct_by_consensus_threshold.pdf")
    plt.close(fig)

    # Figure 2: strict resolved vs unresolved counts for thresholds 1..4.
    rows = sorted(
        [r for r in summary_rows if r["accessibility_definition"] == "shared_all_strict"],
        key=lambda r: intv(r["minimum_accessible_workflows"]),
    )
    xlabels = [f"≥{intv(r['minimum_accessible_workflows'])}" for r in rows]
    resolved = [intv(r["resolved_by_any_workflow"]) for r in rows]
    unresolved = [intv(r["unresolved_by_all_workflows"]) for r in rows]

    x = list(range(len(rows)))
    width = 0.38
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.bar([i - width/2 for i in x], resolved, width=width, label="Resolved by ≥1 workflow")
    ax.bar([i + width/2 for i in x], unresolved, width=width, label="Unresolved by all 4")
    ax.set_xticks(x)
    ax.set_xticklabels(xlabels)
    ax.set_xlabel("Shared-all accessibility confirmed by workflows")
    ax.set_ylabel("Theoretically resolvable families")
    ax.set_title("Consensus-accessible families: resolved versus unresolved")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "Fig_C2_strict_consensus_resolved_unresolved.png", dpi=300)
    fig.savefig(outdir / "Fig_C2_strict_consensus_resolved_unresolved.pdf")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Step01C.2: test whether isoform non-resolution persists when family "
            "accessibility is independently supported by multiple workflows."
        )
    )
    ap.add_argument(
        "--step01c-dir",
        required=True,
        help="Step01C v1.0.0 output directory.",
    )
    ap.add_argument(
        "--step01b-dir",
        required=True,
        help="Frozen Step01B v1.1.1 output directory.",
    )
    ap.add_argument(
        "--outdir",
        default="Step01C2_consensus_accessibility_results",
    )
    args = ap.parse_args()

    step01c_dir = Path(args.step01c_dir)
    step01b_dir = Path(args.step01b_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    eprint(f"[1/5] Step01C.2 v{SCRIPT_VERSION}: loading frozen family matrix")
    family_rows = load_family_matrix(step01c_dir)
    eprint(f"      Multi-entry family rows: {len(family_rows):,}")

    eprint("[2/5] Loading frozen full-FASTA observed peptide classifications")
    obs_rows = load_observed_classification(step01b_dir)
    fs_counts, strict_counts, primary_counts = build_evidence_counts(obs_rows)

    eprint("[3/5] Building accessibility-consensus and resolution matrices")
    family_long = make_family_long(
        family_rows, fs_counts, strict_counts, primary_counts
    )
    summary_rows = make_summary(family_long)

    eprint("[4/5] Validating against frozen Step01C pooled anchors")
    validation_rows = validate_inputs(family_rows, summary_rows)
    write_tsv(outdir / "Step01C2_validation.tsv", validation_rows)

    failed = [r for r in validation_rows if r["status"] != "MATCH"]
    if failed:
        for r in failed:
            eprint(
                f"      CHECK {r['metric']}: observed={r['observed']} "
                f"expected={r['expected']}"
            )
        raise RuntimeError(
            "Step01C.2 validation failed. Inputs do not reproduce the frozen "
            "Step01C v1.0.0 anchors; no interpretation should be made."
        )
    eprint("      All pooled anchors MATCH.")

    eprint("[5/5] Writing consensus summaries, high-confidence families and figures")
    write_tsv(outdir / "consensus_accessibility_summary.tsv", summary_rows)
    write_tsv(outdir / "consensus_accessibility_family_long.tsv", family_long)
    write_high_confidence_tables(outdir, family_long)
    make_figures(outdir, summary_rows)

    # Manuscript-facing compact table for the strict definition.
    strict_rows = [
        r for r in summary_rows
        if r["accessibility_definition"] == "shared_all_strict"
    ]
    write_tsv(outdir / "strict_consensus_summary.tsv", strict_rows)

    manifest = {
        "script_version": SCRIPT_VERSION,
        "workflows": WORKFLOWS,
        "accessibility_definitions": DEFINITIONS,
        "primary_endpoint": (
            "Among families accessible in >=N workflows and theoretically "
            "primary-isoform-resolvable by trypsin, unresolved means no primary "
            "isoform-discriminative peptide was observed by ANY of AP/FP/MM/MQ."
        ),
        "why_endpoint_is_conservative": (
            "Resolution by any one workflow is sufficient to remove a family "
            "from the unresolved group."
        ),
        "inputs": {
            "step01c_dir": str(step01c_dir.resolve()),
            "step01b_dir": str(step01b_dir.resolve()),
        },
    }
    (outdir / "Step01C2_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    eprint("\nStrict shared-all consensus:")
    for r in strict_rows:
        eprint(
            f"      >= {r['minimum_accessible_workflows']} workflows: "
            f"accessible+theoretical={r['consensus_accessible_theoretically_resolvable_families']}, "
            f"resolved_any={r['resolved_by_any_workflow']}, "
            f"unresolved_all4={r['unresolved_by_all_workflows']} "
            f"({float(r['unresolved_pct_of_consensus_accessible_theoretical']):.2f}%)"
        )

    eprint(f"\nDone. Outputs: {outdir.resolve()}")


if __name__ == "__main__":
    main()
