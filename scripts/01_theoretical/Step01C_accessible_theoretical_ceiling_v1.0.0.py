#!/usr/bin/env python3
"""
Step01C_accessible_theoretical_ceiling_v1.0.0.py

Accessible theoretical ceiling for isoform-discriminative bottom-up proteomics.

This analysis CONDITIONS the theoretical tryptic isoform-discriminative space on
protein-entry families that are empirically accessible in the benchmark data.
It is deliberately downstream of the frozen Step01B v1.1.1 observed-peptide
extraction/full-FASTA classification and does not alter those definitions.

Inputs
------
1) Exact common UniProt FASTA used in the benchmark.
2) Step01B v1.1.1 output directory containing observed_peptides.tsv.
3) Step01 v1.0.3 output directory containing peptide_catalog_trypsin.tsv.gz.

Accessibility definitions
-------------------------
A) broad_any
   A multi-entry accession family is accessible if ANY observed peptide maps to
   at least one entry of that family. Cross-family shared peptides are allowed.
   This is deliberately permissive.

B) family_specific_nonprimary
   A multi-entry accession family is accessible if at least one observed peptide
   maps ONLY to that accession family (one base accession) AND the peptide is not
   primary isoform-discriminative. This establishes family accessibility without
   using the primary isoform-discriminative endpoint itself.

C) shared_all_strict
   A multi-entry accession family is accessible if at least one observed peptide
   is classified as within_family_shared_all, i.e. it maps to all reference
   entries of that family and therefore confirms family detection while carrying
   no within-family isoform discrimination. This is the most conservative and
   least circular definition.

For each definition, pooled and workflow-specific (AP/FP/MM/MQ) analyses report:
- accessible multi-entry families
- accessible families that are theoretically isoform-resolvable by trypsin
- theoretical primary isoform-discriminative peptide denominator within those
  accessible families
- observed primary isoform-discriminative peptide recovery
- families resolved by >=1 observed primary peptide
- accessible-but-isoform-unresolved families

No search engine is re-run.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

SCRIPT_VERSION = "1.0.0"
WORKFLOWS = ["AP", "FP", "MM", "MQ"]
EXPECTED_OBSERVED_PRIMARY = {"AP": 180, "FP": 101, "MM": 149, "MQ": 70}
EXPECTED_PRIMARY_UNION = 297
EXPECTED_PRIMARY_INCIDENCES = 500
EXPECTED_THEORETICAL_PRIMARY = 129631
EXPECTED_THEORETICAL_PRIMARY_FAMILIES = 7438
EXPECTED_MULTI_ENTRY_FAMILIES = 10587

ISOFORM_SUFFIX_RE = re.compile(r"-(\d+)$")


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def normalize_il(seq: str) -> str:
    return seq.upper().replace("I", "J").replace("L", "J")


def norm_col(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(x).strip().lower())


def strip_isoform_suffix(accession: str) -> str:
    return ISOFORM_SUFFIX_RE.sub("", accession)


def parse_uniprot_header(header: str) -> Tuple[str, Optional[str]]:
    token = header.split()[0]
    if "|" in token:
        parts = token.split("|")
        accession = parts[1] if len(parts) >= 2 else token
    else:
        accession = token
    m = re.search(r"(?:^|\s)GN=([^\s]+)", header)
    gene = m.group(1) if m else None
    return accession, gene


def write_tsv(path: Path, rows: List[dict], fieldnames: Optional[List[str]] = None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    if fieldnames is None:
        fieldnames = list(rows[0])
    with path.open("wt", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


@dataclass(frozen=True)
class Protein:
    accession: str
    base_accession: str
    gene: Optional[str]
    is_suffixed_isoform: bool


@dataclass(frozen=True)
class ObsClass:
    category: str
    exact_entry: bool
    primary: bool
    structural: bool
    target_family: Optional[str]
    target_gene: Optional[str]
    n_accessions: int
    n_families: int
    n_genes: int


def classify_mapping(
    idxs: Set[int],
    proteins: Sequence[Protein],
    family_totals: Dict[str, int],
    gene_totals: Dict[str, int],
) -> ObsClass:
    if not idxs:
        return ObsClass("unmapped", False, False, False, None, None, 0, 0, 0)

    entries = [proteins[i] for i in idxs]
    families = {p.base_accession for p in entries}
    real_genes = {p.gene for p in entries if p.gene}
    has_missing_gene = any(p.gene is None for p in entries)
    n_genes_effective = len(real_genes) + (1 if has_missing_gene else 0)

    if len(entries) == 1:
        p = entries[0]
        if p.is_suffixed_isoform:
            return ObsClass(
                "single_isoform_unique", True, True, True,
                p.base_accession, p.gene, 1, 1, 1 if p.gene else 0,
            )
        return ObsClass(
            "single_canonical_unique", True, False, True,
            p.base_accession, p.gene, 1, 1, 1 if p.gene else 0,
        )

    if len(families) == 1:
        fam = next(iter(families))
        gene = next(iter(real_genes)) if len(real_genes) == 1 and not has_missing_gene else None
        if len(entries) < family_totals[fam]:
            return ObsClass(
                "within_family_subset_discriminative", False, True, True,
                fam, gene, len(entries), 1, n_genes_effective,
            )
        return ObsClass(
            "within_family_shared_all", False, False, False,
            fam, gene, len(entries), 1, n_genes_effective,
        )

    if len(real_genes) == 1 and not has_missing_gene:
        gene = next(iter(real_genes))
        if len(entries) < gene_totals[gene]:
            return ObsClass(
                "same_gene_subset_discriminative", False, False, True,
                None, gene, len(entries), len(families), 1,
            )
        return ObsClass(
            "same_gene_multi_entry_shared", False, False, False,
            None, gene, len(entries), len(families), 1,
        )

    if len(real_genes) >= 2:
        return ObsClass(
            "cross_gene_shared", False, False, False,
            None, None, len(entries), len(families), n_genes_effective,
        )

    return ObsClass(
        "multi_entry_gene_unresolved", False, False, False,
        None, None, len(entries), len(families), n_genes_effective,
    )


def load_observed(step01b_dir: Path) -> Tuple[Dict[str, Set[str]], Dict[str, str]]:
    path = step01b_dir / "observed_peptides.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Step01B observed peptide table not found: {path}")

    wf_keys: Dict[str, Set[str]] = {wf: set() for wf in WORKFLOWS}
    representatives: Dict[str, str] = {}

    with path.open("rt", encoding="utf-8-sig", errors="replace", newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        fields = set(r.fieldnames or [])
        needed = {"workflow", "peptide_sequence"}
        if not needed.issubset(fields):
            raise ValueError(f"{path} must contain {sorted(needed)}; found {r.fieldnames}")
        for row in r:
            wf = str(row.get("workflow", "")).strip().upper()
            if wf not in wf_keys:
                continue
            seq = str(row.get("peptide_sequence", "")).strip().upper()
            if not seq:
                continue
            key = str(row.get("peptide_key_IL_equivalent", "")).strip().upper()
            if not key:
                key = normalize_il(seq)
            wf_keys[wf].add(key)
            representatives.setdefault(key, seq)

    return wf_keys, representatives


def map_observed_to_fasta(
    fasta: Path,
    observed_keys: Set[str],
) -> Tuple[
    List[Protein], Dict[str, Set[int]], Dict[str, ObsClass], Counter, Counter
]:
    try:
        import ahocorasick  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "Step01C requires pyahocorasick. Install with: python -m pip install pyahocorasick"
        ) from exc

    eprint(f"      Building Aho-Corasick index for {len(observed_keys):,} observed peptide keys")
    automaton = ahocorasick.Automaton()
    for key in sorted(observed_keys):
        automaton.add_word(key, key)
    automaton.make_automaton()

    proteins: List[Protein] = []
    mappings: Dict[str, Set[int]] = defaultdict(set)
    family_totals: Counter = Counter()
    gene_totals: Counter = Counter()

    header: Optional[str] = None
    chunks: List[str] = []

    def flush():
        nonlocal header, chunks
        if header is None:
            return
        seq = "".join(chunks).replace(" ", "").replace("\r", "").upper()
        if not seq:
            header = None
            chunks = []
            return

        accession, gene = parse_uniprot_header(header)
        base = strip_isoform_suffix(accession)
        idx = len(proteins)
        proteins.append(
            Protein(
                accession=accession,
                base_accession=base,
                gene=gene,
                is_suffixed_isoform=bool(ISOFORM_SUFFIX_RE.search(accession)),
            )
        )
        family_totals[base] += 1
        if gene:
            gene_totals[gene] += 1

        seq_key = normalize_il(seq)
        for _, key in automaton.iter(seq_key):
            mappings[key].add(idx)

        header = None
        chunks = []

    with fasta.open("rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
            elif header is not None:
                chunks.append(line.strip())
        flush()

    classes: Dict[str, ObsClass] = {}
    for key in observed_keys:
        classes[key] = classify_mapping(mappings.get(key, set()), proteins, family_totals, gene_totals)

    return proteins, mappings, classes, family_totals, gene_totals


def load_theoretical_primary(
    theoretical_dir: Path,
    multi_entry_families: Set[str],
) -> Tuple[Dict[str, int], Set[str], int, int]:
    path = theoretical_dir / "peptide_catalog_trypsin.tsv.gz"
    if not path.exists():
        raise FileNotFoundError(f"Theoretical trypsin catalog not found: {path}")

    theoretical_by_family: Counter = Counter()
    primary_keys_multi: Set[str] = set()
    primary_all = 0
    primary_without_multi_family = 0

    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        required = {"peptide_key", "primary_isoform_discriminative", "target_base_accession"}
        if not required.issubset(set(r.fieldnames or [])):
            raise ValueError(f"Theoretical catalog missing {sorted(required)}; found {r.fieldnames}")
        for row in r:
            if str(row.get("primary_isoform_discriminative", "0")) != "1":
                continue
            primary_all += 1
            key = row["peptide_key"]
            fam = str(row.get("target_base_accession", "")).strip()
            if fam and fam in multi_entry_families:
                theoretical_by_family[fam] += 1
                primary_keys_multi.add(key)
            else:
                primary_without_multi_family += 1

    return dict(theoretical_by_family), primary_keys_multi, primary_all, primary_without_multi_family


def family_sets_for_scope(
    scope_keys: Set[str],
    mappings: Dict[str, Set[int]],
    classes: Dict[str, ObsClass],
    proteins: Sequence[Protein],
    multi_entry_families: Set[str],
) -> Dict[str, Set[str]]:
    broad: Set[str] = set()
    family_specific_nonprimary: Set[str] = set()
    shared_all_strict: Set[str] = set()

    for key in scope_keys:
        idxs = mappings.get(key, set())
        fams = {proteins[i].base_accession for i in idxs}
        fams &= multi_entry_families
        broad |= fams

        c = classes[key]
        if c.target_family and c.target_family in multi_entry_families:
            if c.n_families == 1 and not c.primary:
                family_specific_nonprimary.add(c.target_family)
            if c.category == "within_family_shared_all":
                shared_all_strict.add(c.target_family)

    return {
        "broad_any": broad,
        "family_specific_nonprimary": family_specific_nonprimary,
        "shared_all_strict": shared_all_strict,
    }


def observed_primary_by_scope(
    scope_keys: Set[str],
    classes: Dict[str, ObsClass],
    multi_entry_families: Set[str],
) -> Tuple[Set[str], Dict[str, Set[str]]]:
    keys: Set[str] = set()
    by_family: Dict[str, Set[str]] = defaultdict(set)
    for key in scope_keys:
        c = classes[key]
        if c.primary and c.target_family and c.target_family in multi_entry_families:
            keys.add(key)
            by_family[c.target_family].add(key)
    return keys, by_family


def create_summaries(
    wf_keys: Dict[str, Set[str]],
    mappings: Dict[str, Set[int]],
    classes: Dict[str, ObsClass],
    proteins: Sequence[Protein],
    family_totals: Dict[str, int],
    theoretical_by_family: Dict[str, int],
) -> Tuple[List[dict], List[dict], Dict[Tuple[str, str], Set[str]]]:
    multi_entry_families = {f for f, n in family_totals.items() if n > 1}
    theoretical_families = set(theoretical_by_family)

    scopes: Dict[str, Set[str]] = {wf: set(wf_keys[wf]) for wf in WORKFLOWS}
    scopes["ALL"] = set().union(*(wf_keys[wf] for wf in WORKFLOWS))

    summary_rows: List[dict] = []
    family_rows: List[dict] = []
    accessible_sets: Dict[Tuple[str, str], Set[str]] = {}

    # Precompute accessibility and observed-primary family sets per scope.
    scope_access: Dict[str, Dict[str, Set[str]]] = {}
    scope_primary_keys: Dict[str, Set[str]] = {}
    scope_primary_by_family: Dict[str, Dict[str, Set[str]]] = {}

    for scope, keys in scopes.items():
        scope_access[scope] = family_sets_for_scope(
            keys, mappings, classes, proteins, multi_entry_families
        )
        pk, pf = observed_primary_by_scope(keys, classes, multi_entry_families)
        scope_primary_keys[scope] = pk
        scope_primary_by_family[scope] = pf

    for definition in ["broad_any", "family_specific_nonprimary", "shared_all_strict"]:
        for scope in ["ALL"] + WORKFLOWS:
            accessible = scope_access[scope][definition]
            accessible_sets[(definition, scope)] = accessible
            accessible_theory = accessible & theoretical_families
            theory_peptides = sum(theoretical_by_family[f] for f in accessible_theory)

            observed_family_map = scope_primary_by_family[scope]
            resolved_families = {f for f in accessible_theory if observed_family_map.get(f)}
            unresolved_families = accessible_theory - resolved_families

            observed_primary_in_accessible = set().union(
                *(observed_family_map.get(f, set()) for f in accessible_theory)
            ) if accessible_theory else set()

            summary_rows.append({
                "accessibility_definition": definition,
                "scope": scope,
                "multi_entry_family_universe": len(multi_entry_families),
                "accessible_multi_entry_families": len(accessible),
                "accessible_family_pct_of_multi_entry_universe": 100.0 * len(accessible) / len(multi_entry_families) if multi_entry_families else math.nan,
                "accessible_theoretically_resolvable_families": len(accessible_theory),
                "accessible_theoretically_resolvable_pct_of_accessible": 100.0 * len(accessible_theory) / len(accessible) if accessible else math.nan,
                "accessible_theoretical_primary_peptides": theory_peptides,
                "observed_primary_peptides_within_accessible_theoretical_families": len(observed_primary_in_accessible),
                "accessible_primary_peptide_recovery_pct": 100.0 * len(observed_primary_in_accessible) / theory_peptides if theory_peptides else math.nan,
                "resolved_accessible_theoretical_families": len(resolved_families),
                "accessible_but_isoform_unresolved_families": len(unresolved_families),
                "resolved_family_pct_of_accessible_theoretical": 100.0 * len(resolved_families) / len(accessible_theory) if accessible_theory else math.nan,
                "unresolved_family_pct_of_accessible_theoretical": 100.0 * len(unresolved_families) / len(accessible_theory) if accessible_theory else math.nan,
            })

    # One family-wide matrix, with all definitions/scopes as flags.
    for fam in sorted(multi_entry_families):
        row = {
            "family": fam,
            "n_reference_entries": family_totals[fam],
            "theoretical_primary_peptides": theoretical_by_family.get(fam, 0),
            "theoretically_primary_resolvable": int(fam in theoretical_families),
            "observed_primary_union": len(scope_primary_by_family["ALL"].get(fam, set())),
        }
        for wf in WORKFLOWS:
            row[f"observed_primary_{wf}"] = len(scope_primary_by_family[wf].get(fam, set()))
        for definition in ["broad_any", "family_specific_nonprimary", "shared_all_strict"]:
            for scope in ["ALL"] + WORKFLOWS:
                acc = fam in scope_access[scope][definition]
                row[f"accessible__{definition}__{scope}"] = int(acc)
                if fam in theoretical_families and acc:
                    resolved = bool(scope_primary_by_family[scope].get(fam))
                    row[f"status__{definition}__{scope}"] = "resolved" if resolved else "accessible_unresolved"
                elif acc:
                    row[f"status__{definition}__{scope}"] = "accessible_not_theoretically_primary_resolvable"
                else:
                    row[f"status__{definition}__{scope}"] = "not_accessible"
        family_rows.append(row)

    return summary_rows, family_rows, accessible_sets


def write_unresolved_tables(
    outdir: Path,
    family_rows: List[dict],
):
    # Primary manuscript-facing table: strict pooled definition.
    strict = [
        r for r in family_rows
        if r.get("status__shared_all_strict__ALL") == "accessible_unresolved"
    ]
    strict = sorted(strict, key=lambda r: (-int(r["theoretical_primary_peptides"]), r["family"]))
    write_tsv(outdir / "accessible_but_isoform_unresolved_families_strict_pooled.tsv", strict)

    # Long-form unresolved table for all definitions/scopes.
    rows: List[dict] = []
    for r in family_rows:
        for definition in ["broad_any", "family_specific_nonprimary", "shared_all_strict"]:
            for scope in ["ALL"] + WORKFLOWS:
                if r.get(f"status__{definition}__{scope}") == "accessible_unresolved":
                    rows.append({
                        "accessibility_definition": definition,
                        "scope": scope,
                        "family": r["family"],
                        "n_reference_entries": r["n_reference_entries"],
                        "theoretical_primary_peptides": r["theoretical_primary_peptides"],
                        "observed_primary_in_scope": r["observed_primary_union"] if scope == "ALL" else r[f"observed_primary_{scope}"],
                    })
    write_tsv(outdir / "accessible_but_isoform_unresolved_families_all_definitions.tsv", rows)


def validate(
    wf_keys: Dict[str, Set[str]],
    classes: Dict[str, ObsClass],
    family_totals: Dict[str, int],
    theoretical_by_family: Dict[str, int],
    primary_all: int,
) -> List[dict]:
    rows: List[dict] = []
    multi = {f for f, n in family_totals.items() if n > 1}

    def add(metric, observed, expected):
        rows.append({
            "metric": metric,
            "observed": observed,
            "expected": expected,
            "difference": observed - expected,
            "status": "MATCH" if observed == expected else "CHECK",
        })

    add("multi_entry_accession_families", len(multi), EXPECTED_MULTI_ENTRY_FAMILIES)
    add("theoretical_primary_peptide_keys_all", primary_all, EXPECTED_THEORETICAL_PRIMARY)
    add("theoretical_primary_multi_entry_families", len(theoretical_by_family), EXPECTED_THEORETICAL_PRIMARY_FAMILIES)

    union_primary: Set[str] = set()
    incidences = 0
    for wf in WORKFLOWS:
        obs = {
            k for k in wf_keys[wf]
            if classes[k].primary and classes[k].target_family in multi
        }
        union_primary |= obs
        incidences += len(obs)
        add(f"observed_primary_{wf}", len(obs), EXPECTED_OBSERVED_PRIMARY[wf])

    add("observed_primary_union", len(union_primary), EXPECTED_PRIMARY_UNION)
    add("observed_primary_workflow_peptide_incidences", incidences, EXPECTED_PRIMARY_INCIDENCES)
    return rows


def make_figures(outdir: Path, summary_rows: List[dict]):
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception as exc:
        eprint(f"[warning] matplotlib unavailable; figures skipped: {exc}")
        return

    # Figure 1: accessible theoretical primary peptide recovery by workflow.
    definitions = ["broad_any", "family_specific_nonprimary", "shared_all_strict"]
    for definition in definitions:
        rows = [r for r in summary_rows if r["accessibility_definition"] == definition and r["scope"] in WORKFLOWS]
        rows.sort(key=lambda r: WORKFLOWS.index(r["scope"]))
        labels = [r["scope"] for r in rows]
        vals = [r["accessible_primary_peptide_recovery_pct"] for r in rows]
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        ax.bar(labels, vals)
        ax.set_ylabel("Observed primary peptides / accessible\ntheoretical primary peptides (%)")
        ax.set_xlabel("Workflow")
        ax.set_title(f"Accessible theoretical isoform-peptide recovery: {definition}")
        ax.grid(axis="y", alpha=0.25)
        fig.tight_layout()
        fig.savefig(outdir / f"Fig_Step01C_recovery_{definition}.png", dpi=300)
        fig.savefig(outdir / f"Fig_Step01C_recovery_{definition}.pdf")
        plt.close(fig)

    # Figure 2: strict pooled accessible families resolved vs unresolved.
    strict_all = next(
        r for r in summary_rows
        if r["accessibility_definition"] == "shared_all_strict" and r["scope"] == "ALL"
    )
    labels = ["Resolved", "Accessible but\nisoform-unresolved"]
    vals = [
        strict_all["resolved_accessible_theoretical_families"],
        strict_all["accessible_but_isoform_unresolved_families"],
    ]
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    ax.bar(labels, vals)
    ax.set_ylabel("Multi-entry accession families")
    ax.set_title("Strictly accessible, theoretically resolvable families")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "Fig_Step01C_strict_resolved_vs_unresolved.png", dpi=300)
    fig.savefig(outdir / "Fig_Step01C_strict_resolved_vs_unresolved.pdf")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(
        description="Condition the theoretical tryptic isoform-discriminative space on empirically accessible accession families."
    )
    ap.add_argument("--fasta", required=True, help="Exact common UniProt FASTA used in the benchmark.")
    ap.add_argument("--step01b-dir", required=True, help="Step01B v1.1.1 output directory.")
    ap.add_argument("--theoretical-dir", required=True, help="Step01 v1.0.3 results directory.")
    ap.add_argument("--outdir", default="Step01C_accessible_results")
    args = ap.parse_args()

    fasta = Path(args.fasta)
    step01b_dir = Path(args.step01b_dir)
    theoretical_dir = Path(args.theoretical_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for p, label in [(fasta, "FASTA"), (step01b_dir, "Step01B directory"), (theoretical_dir, "Step01 directory")]:
        if not p.exists():
            raise SystemExit(f"{label} not found: {p}")

    eprint(f"[1/6] Step01C v{SCRIPT_VERSION}: loading frozen Step01B observed peptide sets")
    wf_keys, representatives = load_observed(step01b_dir)
    for wf in WORKFLOWS:
        eprint(f"      {wf}: {len(wf_keys[wf]):,} observed I/L-equivalent keys")

    eprint("[2/6] Remapping observed peptide keys to the exact common FASTA")
    union_keys = set().union(*(wf_keys[wf] for wf in WORKFLOWS))
    proteins, mappings, classes, family_totals, gene_totals = map_observed_to_fasta(fasta, union_keys)
    multi_entry_families = {f for f, n in family_totals.items() if n > 1}
    eprint(f"      Protein entries: {len(proteins):,}")
    eprint(f"      Multi-entry accession families: {len(multi_entry_families):,}")

    eprint("[3/6] Loading theoretical primary isoform-discriminative tryptic peptide space")
    theoretical_by_family, theoretical_primary_keys_multi, primary_all, primary_without_multi = load_theoretical_primary(
        theoretical_dir, multi_entry_families
    )
    eprint(f"      Theoretical primary peptide keys (all): {primary_all:,}")
    eprint(f"      Multi-entry families with >=1 theoretical primary peptide: {len(theoretical_by_family):,}")
    if primary_without_multi:
        eprint(f"      Note: {primary_without_multi:,} primary theoretical keys are not assigned to a multi-entry family and are excluded from family-conditioned denominators")

    eprint("[4/6] Computing broad and conservative accessibility-conditioned ceilings")
    summary_rows, family_rows, accessible_sets = create_summaries(
        wf_keys, mappings, classes, proteins, family_totals, theoretical_by_family
    )

    eprint("[5/6] Writing tables, validation and figures")
    write_tsv(outdir / "accessible_ceiling_summary.tsv", summary_rows)
    write_tsv(outdir / "accessible_family_status_matrix.tsv", family_rows)
    write_unresolved_tables(outdir, family_rows)

    validation_rows = validate(
        wf_keys, classes, family_totals, theoretical_by_family, primary_all
    )
    write_tsv(outdir / "Step01C_validation.tsv", validation_rows)

    manifest = {
        "script_version": SCRIPT_VERSION,
        "accessibility_definitions": {
            "broad_any": "Any observed peptide maps to at least one entry in the multi-entry accession family; cross-family shared peptides allowed.",
            "family_specific_nonprimary": "At least one observed peptide maps only to this accession family and is not primary isoform-discriminative.",
            "shared_all_strict": "At least one observed within_family_shared_all peptide maps to all entries of the family; no within-family isoform discrimination is used to establish accessibility.",
        },
        "family_level_endpoint": "accessible-but-isoform-unresolved = accessible AND theoretically primary-resolvable by trypsin AND no observed primary isoform-discriminative peptide in the same scope.",
        "important_interpretation": "The accessible theoretical denominator is conditioned on empirical family accessibility. It is not a direct measure of MS detectability and remains influenced by peptide physicochemical properties, abundance, LC-MS/MS sampling and workflow-specific identification.",
        "inputs": {
            "fasta": str(fasta.resolve()),
            "step01b_dir": str(step01b_dir.resolve()),
            "theoretical_dir": str(theoretical_dir.resolve()),
        },
    }
    (outdir / "Step01C_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    make_figures(outdir, summary_rows)

    eprint("[6/6] Self-check")
    failed = [r for r in validation_rows if r["status"] != "MATCH"]
    for r in validation_rows:
        eprint(f"      {r['metric']}: {r['observed']}/{r['expected']}  {r['status']}")

    if failed:
        eprint("\n[important] At least one frozen-analysis validation metric differs. Do not interpret Step01C until the mismatch is resolved.")
    else:
        eprint("\n[success] Frozen Step01A/Step01B anchors reproduce exactly. Step01C accessibility-conditioned results are ready for interpretation.")

    # Print manuscript-facing strict pooled result directly.
    strict = next(r for r in summary_rows if r["accessibility_definition"] == "shared_all_strict" and r["scope"] == "ALL")
    eprint("\n[strict pooled headline]")
    eprint(f"      Strictly accessible multi-entry families: {strict['accessible_multi_entry_families']:,}")
    eprint(f"      Strictly accessible + theoretically resolvable families: {strict['accessible_theoretically_resolvable_families']:,}")
    eprint(f"      Accessible theoretical primary peptides: {strict['accessible_theoretical_primary_peptides']:,}")
    eprint(f"      Observed primary peptides within those families: {strict['observed_primary_peptides_within_accessible_theoretical_families']:,}")
    eprint(f"      Accessible theoretical peptide recovery: {strict['accessible_primary_peptide_recovery_pct']:.4f}%")
    eprint(f"      Resolved families: {strict['resolved_accessible_theoretical_families']:,}")
    eprint(f"      Accessible but isoform-unresolved families: {strict['accessible_but_isoform_unresolved_families']:,} ({strict['unresolved_family_pct_of_accessible_theoretical']:.2f}%)")

    eprint(f"\nDone. Outputs: {outdir.resolve()}")


if __name__ == "__main__":
    main()
