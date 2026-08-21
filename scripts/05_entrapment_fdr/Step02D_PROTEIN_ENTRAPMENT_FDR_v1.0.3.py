#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step02D_PROTEIN_ENTRAPMENT_FDR_v1.0.0.py

Protein-level entrapment FDR evaluation for completed AP / FP / MM / MQ reruns.

Design
------
The Step02A database is a 1:1 PROTEIN-LEVEL shuffled entrapment database.
Therefore this script evaluates the completed reruns at the reported protein/
protein-group level.

Primary classification (PURE-GROUP analysis)
--------------------------------------------
Each reported protein group is classified using all group members that can be
resolved against the frozen target+entrapment FASTA:

    TARGET_ONLY
    ENTRAPMENT_ONLY
    MIXED_TARGET_ENTRAPMENT
    UNRESOLVED
    DECOY_OR_CONTAMINANT_EXCLUDED

The primary FDP calculation excludes MIXED and UNRESOLVED groups:

    lower FDP    = E / (T + E)
    combined FDP = E * (1 + 1/r) / (T + E)

where r is the effective entrapment:target protein-space ratio. For the frozen
1:1 Step02A design r_protein = entrapment_proteins / target_proteins = 1.0.

A representative-protein sensitivity analysis is also written, but the
PURE-GROUP result is the primary output.

Workflow sources
----------------
AP:
    prefer results_protein_summary.csv
    then results_proteins.csv
    then results.hdf/protein_table as fallback

FP:
    prefer combined_protein.tsv
    then protein.tsv

MM:
    prefer AllProteinGroups*.tsv/txt
    filter QValue <= 0.01 where available
    exclude D/C rows where DCT is available

MQ:
    proteinGroups.txt
    exclude Reverse == '+' and Potential contaminant == '+'

No figure is produced here. Use R after the frozen evaluation table passes QC.

Default root
------------
C:\\Users\\Supachai\\Desktop\\AphaPept_benchmark\\Benchmark_Program\\ENTRAPMENT_FDR
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from itertools import chain
from collections import Counter
from pathlib import Path

# ----------------------------------------------------------------------
# Public-release project-root resolver
# ----------------------------------------------------------------------
def _public_project_root():
    import os as _os

    env = _os.environ.get("ISOFORM_BENCHMARK_ROOT")
    if env:
        p = Path(env).expanduser().resolve()
        if not p.exists():
            raise RuntimeError(
                "ISOFORM_BENCHMARK_ROOT does not exist: " + str(p)
            )
        if p.name != "Benchmark_Program":
            raise RuntimeError(
                "ISOFORM_BENCHMARK_ROOT must point to Benchmark_Program: "
                + str(p)
            )
        return p

    script_path = Path(__file__).resolve()
    for candidate in [script_path.parent] + list(script_path.parents):
        if candidate.name == "Benchmark_Program":
            return candidate

    raise RuntimeError(
        "Could not locate Benchmark_Program. "
        "Set ISOFORM_BENCHMARK_ROOT to the Benchmark_Program folder."
    )


VERSION = "1.0.3"
WORKFLOWS = ["AP", "FP", "MM", "MQ"]

DEFAULT_ROOT = _public_project_root() / "ENTRAPMENT_FDR"

EXPECTED_STEP02A = {
    "target_proteins": 169637,
    "entrapment_proteins": 169637,
}

# =============================================================================
# General helpers
# =============================================================================

def norm_col(x):
    return re.sub(r"[^a-z0-9]+", "", str(x or "").strip().lower())

def find_col(fields, aliases):
    by_norm = {norm_col(x): x for x in fields}
    for alias in aliases:
        n = norm_col(alias)
        if n in by_norm:
            return by_norm[n]
    for alias in aliases:
        n = norm_col(alias)
        for col in fields:
            if n and n in norm_col(col):
                return col
    return None

def truthy(value):
    return str(value or "").strip().lower() in {
        "1", "true", "t", "yes", "y", "+",
        "decoy", "reverse", "contaminant"
    }

def float_or_none(value):
    try:
        if value is None or str(value).strip() == "":
            return None
        return float(value)
    except Exception:
        return None

def _raise_csv_field_limit():
    """Increase Python CSV field limit for very large protein-group fields."""
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return limit
        except OverflowError:
            limit //= 10


CSV_FIELD_LIMIT = _raise_csv_field_limit()


def read_delimited(path):
    delim = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open(
        "r",
        encoding="utf-8-sig",
        errors="replace",
        newline=""
    ) as fh:
        yield from csv.DictReader(fh, delimiter=delim)

def write_tsv(path, rows, fieldnames):
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

def locate_step02a(root):
    direct = root / "Step02A_entrapment_db_v104"
    if direct.exists():
        return direct

    hits = [
        p for p in root.rglob("*")
        if p.is_dir()
        and "step02a_entrapment_db" in p.name.lower()
    ]
    if not hits:
        raise FileNotFoundError(
            f"Could not find Step02A entrapment database under:\n{root}"
        )
    return sorted(hits, key=lambda p: len(str(p)))[0]

def locate_workflow_dir(root, workflow):
    exact = [
        p for p in root.rglob("*")
        if p.is_dir()
        and p.name.lower() == f"{workflow}_entrap".lower()
    ]
    if exact:
        return sorted(exact, key=lambda p: len(str(p)))[0]

    token = f"{workflow}_entrap".lower()
    hits = [
        p for p in root.rglob("*")
        if p.is_dir()
        and token in p.name.lower()
    ]
    if not hits:
        raise FileNotFoundError(
            f"No {workflow}_ENTRAP directory found under:\n{root}"
        )
    return sorted(hits, key=lambda p: len(str(p)))[0]

# =============================================================================
# Frozen FASTA dictionary
# =============================================================================

def parse_fasta(path):
    header = None
    seq_parts = []

    with path.open("rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)

                header = line[1:].strip()
                seq_parts = []

            elif header is not None:
                seq_parts.append(line.strip())

        if header is not None:
            yield header, "".join(seq_parts)

def is_entrapment_header(header):
    s = str(header or "").lower()
    return (
        "_p_target" in s
        or "p_target" in s
        or "entrap" in s
    )

def is_decoy_text(text):
    s = str(text or "").strip().lower()
    return (
        s.startswith("rev_")
        or "|rev_" in s
        or s.startswith("decoy_")
        or "|decoy_" in s
        or s.startswith("reverse_")
    )

def header_aliases(header):
    """
    Generate aliases for robust matching of protein IDs emitted by different
    software. The aliases are derived only from the frozen FASTA header.
    """
    aliases = set()

    raw = str(header).strip()
    if not raw:
        return aliases

    first = raw.split()[0]
    aliases.add(raw)
    aliases.add(first)

    # Pipe-delimited UniProt-like identifiers.
    if "|" in first:
        parts = [x for x in first.split("|") if x]
        aliases.update(parts)

        if len(parts) >= 2:
            aliases.add(parts[1])

    # Also add punctuation-clean token variants.
    for token in re.split(r"[;,\s|]+", first):
        token = token.strip()
        if token:
            aliases.add(token)

    return {a for a in aliases if a}

def build_fasta_alias_index(fasta):
    """
    alias -> one of:
        TARGET
        ENTRAPMENT
        COLLISION
    """
    alias_class = {}
    protein_counts = Counter()

    for header, _seq in parse_fasta(fasta):
        cls = "ENTRAPMENT" if is_entrapment_header(header) else "TARGET"
        protein_counts[cls] += 1

        for alias in header_aliases(header):
            previous = alias_class.get(alias)

            if previous is None:
                alias_class[alias] = cls
            elif previous != cls:
                alias_class[alias] = "COLLISION"

    return alias_class, protein_counts

# =============================================================================
# Protein member parsing/classification
# =============================================================================

def split_members(text):
    """
    Split a protein-group field while preserving useful UniProt tokens.

    Group formats differ among tools. We first split on common list delimiters,
    then retain each non-empty chunk. Full header chunks are allowed because
    classify_member() also detects the entrapment marker directly.
    """
    if text is None:
        return []

    s = str(text).strip()
    if not s:
        return []

    # Common group member delimiters.
    chunks = re.split(r"[;,]+", s)

    out = []
    for chunk in chunks:
        chunk = chunk.strip()
        if chunk:
            out.append(chunk)

    return out

def candidate_aliases_from_member(member):
    member = str(member).strip()
    aliases = set()

    if not member:
        return aliases

    aliases.add(member)

    first = member.split()[0]
    aliases.add(first)

    if "|" in first:
        parts = [x for x in first.split("|") if x]
        aliases.update(parts)
        if len(parts) >= 2:
            aliases.add(parts[1])

    # Protein-group outputs often contain identifiers with colon prefixes.
    if ":" in first:
        aliases.add(first.split(":")[-1])

    # Remove a native decoy prefix only for lookup diagnostics.
    for prefix in ("rev_", "REV_", "DECOY_", "decoy_"):
        if first.startswith(prefix):
            aliases.add(first[len(prefix):])

    return {x.strip() for x in aliases if x.strip()}

def classify_member(member, alias_index):
    """
    Returns:
        TARGET
        ENTRAPMENT
        DECOY
        COLLISION
        UNRESOLVED
    """
    raw = str(member or "").strip()

    if not raw:
        return "UNRESOLVED"

    if is_decoy_text(raw):
        return "DECOY"

    # Direct marker is strongest evidence.
    if is_entrapment_header(raw):
        return "ENTRAPMENT"

    classes = set()

    for alias in candidate_aliases_from_member(raw):
        cls = alias_index.get(alias)
        if cls:
            classes.add(cls)

    classes.discard("COLLISION")

    if classes == {"TARGET"}:
        return "TARGET"

    if classes == {"ENTRAPMENT"}:
        return "ENTRAPMENT"

    if "TARGET" in classes and "ENTRAPMENT" in classes:
        return "COLLISION"

    return "UNRESOLVED"

def classify_group(members, alias_index):
    """
    Pure-group primary classification.

    Decoy members are ignored only after the row-level native decoy exclusion.
    If a surviving reported group contains both real target and entrapment
    members, it is MIXED and excluded from primary T/E FDP calculation.
    """
    member_classes = [
        classify_member(x, alias_index)
        for x in members
    ]

    informative = [
        x for x in member_classes
        if x not in {"DECOY"}
    ]

    has_target = "TARGET" in informative
    has_entrap = "ENTRAPMENT" in informative
    has_unresolved = any(
        x in {"UNRESOLVED", "COLLISION"}
        for x in informative
    )

    if has_target and has_entrap:
        cls = "MIXED_TARGET_ENTRAPMENT"
    elif has_target and not has_entrap:
        cls = "TARGET_ONLY" if not has_unresolved else "TARGET_PLUS_UNRESOLVED"
    elif has_entrap and not has_target:
        cls = "ENTRAPMENT_ONLY" if not has_unresolved else "ENTRAPMENT_PLUS_UNRESOLVED"
    else:
        cls = "UNRESOLVED"

    return cls, member_classes

# =============================================================================
# AP parser
# =============================================================================

def parse_ap(workflow_dir, alpha):
    candidates = []

    # Best: protein summary after protein FDR.
    for name in [
        "results_protein_summary.csv",
        "results_proteins.csv",
    ]:
        candidates.extend(workflow_dir.rglob(name))

    candidates = sorted(
        set(candidates),
        key=lambda p: (
            0 if p.name.lower() == "results_protein_summary.csv" else 1,
            len(str(p)),
        ),
    )

    rows_out = []
    source_files = []

    if candidates:
        path = candidates[0]
        rows = read_delimited(path)
        first = next(rows, None)

        if first is None:
            raise RuntimeError(f"[AP] Empty protein result file: {path}")

        fields = list(first.keys())

        # AlphaPept results_protein_summary.csv is commonly written by
        # pandas with the protein/protein-group identifier in the unnamed
        # first index column. csv.DictReader therefore exposes that column
        # literally as "".
        #
        # Prefer the unnamed index when present; otherwise fall back to
        # conventional protein/protein-group column names.
        if "" in fields:
            group_col = ""
        else:
            group_col = find_col(
                fields,
                [
                    "protein_group",
                    "protein group",
                    "protein",
                    "protein_id",
                    "protein ids",
                    "index",
                ],
            )

        decoy_col = find_col(
            fields,
            [
                "decoy_protein",
                "decoy protein",
                "decoy",
                "reverse",
            ],
        )

        q_col = find_col(
            fields,
            [
                "q_value",
                "qvalue",
                "protein_q_value",
                "fdr",
            ],
        )

        if group_col is None:
            raise RuntimeError(
                "[AP] Could not resolve a protein/protein-group column in:\n"
                f"{path}\nColumns: {fields}\n"
                "Expected either an unnamed pandas index column ('') or a "
                "named protein/protein-group column."
            )

        # Sanity check: the selected group/index column must contain at least
        # one non-empty value before we proceed.
        preview_values = []
        preview_first = first.get(group_col)
        if preview_first is not None and str(preview_first).strip():
            preview_values.append(str(preview_first).strip())
        if not preview_values:
            raise RuntimeError(
                f"[AP] Selected protein-group column {group_col!r} is empty "
                f"in the first data row of {path}."
            )

        for i, row in enumerate(chain([first], rows), start=1):
            if decoy_col and truthy(row.get(decoy_col)):
                continue

            if q_col:
                q = float_or_none(row.get(q_col))
                if q is not None and q > alpha:
                    continue

            members = split_members(row.get(group_col))

            if not members:
                continue

            rows_out.append({
                "native_group_id": f"AP:{i}",
                "representative": members[0],
                "members": members,
                "native_q_value": row.get(q_col, "") if q_col else "",
                "native_row_text": row.get(group_col, ""),
            })

        source_files.append(path)

        return rows_out, source_files, {
            "source_type": path.name,
            "group_column": group_col,
            "decoy_column": decoy_col or "",
            "q_column": q_col or "",
        }

    # HDF fallback.
    results_hdfs = sorted(
        p for p in workflow_dir.rglob("*.hdf")
        if "result" in p.name.lower()
    )

    if not results_hdfs:
        raise FileNotFoundError(
            "[AP] Could not find results_protein_summary.csv, "
            "results_proteins.csv, or results.hdf."
        )

    try:
        import pandas as pd
    except ImportError:
        raise RuntimeError(
            "[AP] HDF fallback requires pandas. "
            "Install with: python -m pip install pandas tables"
        )

    path = results_hdfs[0]

    try:
        df = pd.read_hdf(path, "protein_table")
    except Exception as e:
        raise RuntimeError(
            f"[AP] Could not read protein_table from {path}\n{e}"
        )

    df = df.reset_index()
    fields = list(df.columns)

    group_col = find_col(
        fields,
        [
            "protein_group",
            "protein group",
            "protein",
            "index",
        ],
    )

    if group_col is None:
        raise RuntimeError(
            f"[AP] No protein-group column in HDF protein_table. "
            f"Columns: {fields}"
        )

    for i, row in df.iterrows():
        members = split_members(row.get(group_col))
        if not members:
            continue

        rows_out.append({
            "native_group_id": f"AP:{i+1}",
            "representative": members[0],
            "members": members,
            "native_q_value": "",
            "native_row_text": str(row.get(group_col, "")),
        })

    source_files.append(path)

    return rows_out, source_files, {
        "source_type": "results.hdf/protein_table",
        "group_column": group_col,
    }

# =============================================================================
# FragPipe parser
# =============================================================================

def parse_fp(workflow_dir, alpha):
    combined = sorted(workflow_dir.rglob("combined_protein.tsv"))
    individual = sorted(workflow_dir.rglob("protein.tsv"))

    paths = combined or individual

    if not paths:
        raise FileNotFoundError(
            f"[FP] No combined_protein.tsv or protein.tsv under:\n{workflow_dir}"
        )

    # Prefer one combined table. Otherwise pool protein.tsv files and deduplicate
    # later by native group text.
    if combined:
        paths = [combined[0]]

    rows_out = []
    source_files = []

    seen_group_signatures = set()

    for path in paths:
        rows = read_delimited(path)
        first = next(rows, None)

        if first is None:
            continue

        fields = list(first.keys())

        protein_col = find_col(
            fields,
            ["Protein", "Protein ID"]
        )

        indist_col = find_col(
            fields,
            [
                "Indistinguishable Proteins",
                "Indistinguishable Protein",
            ],
        )

        if protein_col is None:
            raise RuntimeError(
                f"[FP] No Protein column in {path}\nColumns: {fields}"
            )

        source_files.append(path)

        for i, row in enumerate(chain([first], rows), start=1):
            members = []

            members.extend(split_members(row.get(protein_col)))

            if indist_col:
                members.extend(split_members(row.get(indist_col)))

            # preserve order, remove exact duplicates
            members = list(dict.fromkeys(members))

            if not members:
                continue

            sig = tuple(members)

            if sig in seen_group_signatures:
                continue

            seen_group_signatures.add(sig)

            rows_out.append({
                "native_group_id": f"FP:{i}",
                "representative": members[0],
                "members": members,
                "native_q_value": "",
                "native_row_text": " ; ".join(members),
            })

    return rows_out, source_files, {
        "source_type": "combined_protein.tsv" if combined else "protein.tsv",
        "fdr_status": "native FDR-filtered protein report",
    }

# =============================================================================
# MetaMorpheus parser
# =============================================================================

def _mm_candidate_score(path, fields):
    """
    Rank MetaMorpheus protein-group-like outputs by filename and schema.
    The official/common output is AllProteinGroups.tsv, but this parser also
    supports renamed/exported variants as long as the header contains a
    protein accession/group column.
    """
    name = path.name.lower()

    score = 0

    if "allproteingroups" in name:
        score += 100
    elif "proteingroups" in name:
        score += 80
    elif "protein" in name and "group" in name:
        score += 60
    elif "protein" in name:
        score += 20

    protein_col = find_col(
        fields,
        [
            "Protein Accession",
            "Protein Accessions",
            "Protein Group",
            "Protein Groups",
            "Protein",
            "Protein IDs",
            "Protein ID",
        ],
    )

    q_col = find_col(
        fields,
        [
            "QValue",
            "Q Value",
            "q-value",
            "q_value",
            "qvalue",
        ],
    )

    dct_col = find_col(
        fields,
        [
            "Decoy/Contaminant/Target",
            "Decoy Contaminant Target",
            "DCT",
        ],
    )

    if protein_col:
        score += 50
    if q_col:
        score += 15
    if dct_col:
        score += 10

    # Prefer tabular outputs over logs/configs.
    if path.suffix.lower() in {".tsv", ".txt", ".csv"}:
        score += 5

    return score, protein_col, q_col, dct_col


def parse_mm(workflow_dir, alpha):
    """
    MetaMorpheus protein-group parser.

    Primary target:
        AllProteinGroups.tsv

    Fallback:
        Any tabular output under MM_ENTRAP whose filename/header indicate a
        protein-group table. This avoids failing solely because the output was
        renamed or exported to a different filename.

    Native filtering:
        QValue <= alpha, when present
        DCT rows beginning with D or C are excluded, when present
    """

    allowed_suffixes = {".tsv", ".txt", ".csv"}

    all_files = [
        p for p in workflow_dir.rglob("*")
        if p.is_file()
        and p.suffix.lower() in allowed_suffixes
    ]

    inspected = []
    viable = []

    for path in all_files:
        # First prioritize filenames that plausibly contain protein data.
        lname = path.name.lower()

        filename_plausible = (
            "protein" in lname
            or "group" in lname
        )

        if not filename_plausible:
            continue

        try:
            rows = read_delimited(path)
            first = next(rows, None)
        except Exception:
            continue

        if first is None:
            continue

        fields = list(first.keys())

        score, protein_col, q_col, dct_col = _mm_candidate_score(
            path,
            fields,
        )

        inspected.append({
            "path": path,
            "score": score,
            "protein_col": protein_col,
            "q_col": q_col,
            "dct_col": dct_col,
            "fields": fields,
        })

        if protein_col is not None:
            viable.append(inspected[-1])

    if not viable:
        # Helpful failure report: show filenames in MM_ENTRAP that mention
        # protein/group so we can immediately see the actual MetaMorpheus
        # output naming on this run.
        names = sorted({
            str(p.relative_to(workflow_dir))
            for p in all_files
            if (
                "protein" in p.name.lower()
                or "group" in p.name.lower()
            )
        })

        preview = "\n".join(
            f"  - {x}"
            for x in names[:80]
        )

        raise FileNotFoundError(
            "[MM] Could not find a parseable MetaMorpheus protein-group "
            "table under:\n"
            f"{workflow_dir}\n\n"
            "Files containing 'protein' or 'group':\n"
            f"{preview if preview else '  (none)'}\n\n"
            "The common MetaMorpheus output is AllProteinGroups.tsv."
        )

    # Highest schema/name score; then larger file; then shorter path.
    viable.sort(
        key=lambda x: (
            -x["score"],
            -x["path"].stat().st_size,
            len(str(x["path"])),
        )
    )

    chosen = viable[0]

    path = chosen["path"]
    protein_col = chosen["protein_col"]
    q_col = chosen["q_col"]
    dct_col = chosen["dct_col"]

    rows = read_delimited(path)
    first = next(rows, None)

    if first is None:
        raise RuntimeError(f"[MM] Selected file is empty: {path}")

    rows_out = []

    for i, row in enumerate(chain([first], rows), start=1):

        if q_col:
            q = float_or_none(row.get(q_col))

            if q is not None and q > alpha:
                continue

        if dct_col:
            dct = str(
                row.get(dct_col, "")
            ).strip().upper()

            if dct.startswith("D") or dct.startswith("C"):
                continue

        members = split_members(
            row.get(protein_col)
        )

        if not members:
            continue

        rows_out.append({
            "native_group_id": f"MM:{i}",
            "representative": members[0],
            "members": members,
            "native_q_value":
                row.get(q_col, "")
                if q_col
                else "",
            "native_row_text":
                row.get(protein_col, ""),
        })

    if not rows_out:
        raise RuntimeError(
            "[MM] A protein-group-like file was found but no target rows "
            "survived parsing/filtering.\n"
            f"Selected file: {path}\n"
            f"Protein column: {protein_col!r}\n"
            f"Q column: {q_col!r}\n"
            f"DCT column: {dct_col!r}"
        )

    meta = {
        "source_type": path.name,
        "selected_source": str(path),
        "candidate_score": chosen["score"],
        "protein_column": protein_col,
        "q_column": q_col or "",
        "dct_column": dct_col or "",
        "n_viable_candidate_tables": len(viable),
        "viable_candidates": [
            {
                "path": str(x["path"]),
                "score": x["score"],
                "protein_col": x["protein_col"],
                "q_col": x["q_col"] or "",
                "dct_col": x["dct_col"] or "",
            }
            for x in viable[:20]
        ],
    }

    return rows_out, [path], meta

# =============================================================================
# MaxQuant parser
# =============================================================================

def parse_mq(workflow_dir, alpha):
    paths = sorted(workflow_dir.rglob("proteinGroups.txt"))

    if not paths:
        raise FileNotFoundError(
            f"[MQ] No proteinGroups.txt under:\n{workflow_dir}"
        )

    paths.sort(
        key=lambda p: (
            0 if p.parent.name.lower() == "txt" else 1,
            len(str(p)),
        )
    )

    path = paths[0]

    rows = read_delimited(path)
    first = next(rows, None)

    if first is None:
        raise RuntimeError(f"[MQ] Empty file: {path}")

    fields = list(first.keys())

    group_col = find_col(
        fields,
        [
            "Protein IDs",
            "Majority protein IDs",
            "Leading razor protein",
        ],
    )

    representative_col = find_col(
        fields,
        [
            "Leading razor protein",
            "Majority protein IDs",
            "Protein IDs",
        ],
    )

    reverse_col = find_col(fields, ["Reverse"])
    contaminant_col = find_col(
        fields,
        ["Potential contaminant", "Potential Contaminant"]
    )

    if group_col is None:
        raise RuntimeError(
            f"[MQ] Could not resolve Protein IDs column.\nColumns: {fields}"
        )

    rows_out = []

    for i, row in enumerate(chain([first], rows), start=1):
        if reverse_col and truthy(row.get(reverse_col)):
            continue

        if contaminant_col and truthy(row.get(contaminant_col)):
            continue

        members = split_members(row.get(group_col))

        if not members:
            continue

        rep_members = split_members(
            row.get(representative_col)
            if representative_col
            else ""
        )

        representative = (
            rep_members[0]
            if rep_members
            else members[0]
        )

        rows_out.append({
            "native_group_id": f"MQ:{i}",
            "representative": representative,
            "members": members,
            "native_q_value": "",
            "native_row_text": row.get(group_col, ""),
        })

    return rows_out, [path], {
        "source_type": "proteinGroups.txt",
        "group_column": group_col,
        "representative_column": representative_col or "",
        "reverse_column": reverse_col or "",
        "contaminant_column": contaminant_col or "",
    }

# =============================================================================
# Evaluation
# =============================================================================

def evaluate_workflow(
    workflow,
    rows,
    alias_index,
    r_protein,
    nominal,
):
    detail_rows = []
    primary_counts = Counter()
    rep_counts = Counter()

    for row in rows:
        group_class, member_classes = classify_group(
            row["members"],
            alias_index,
        )

        primary_counts[group_class] += 1

        rep_class = classify_member(
            row["representative"],
            alias_index,
        )

        rep_counts[rep_class] += 1

        detail_rows.append({
            "workflow": workflow,
            "native_group_id": row["native_group_id"],
            "primary_group_class": group_class,
            "representative_class": rep_class,
            "representative": row["representative"],
            "n_members": len(row["members"]),
            "members": " ; ".join(row["members"]),
            "member_classes": " ; ".join(member_classes),
            "native_q_value": row.get("native_q_value", ""),
            "native_row_text": row.get("native_row_text", ""),
        })

    # Primary: only PURE target-only and entrapment-only groups.
    T = primary_counts["TARGET_ONLY"]
    E = primary_counts["ENTRAPMENT_ONLY"]

    denom = T + E

    lower = E / denom if denom else math.nan

    combined = (
        E * (1.0 + 1.0 / r_protein) / denom
        if denom and r_protein > 0
        else math.nan
    )

    # Sensitivity: representative protein only.
    T_rep = rep_counts["TARGET"]
    E_rep = rep_counts["ENTRAPMENT"]
    denom_rep = T_rep + E_rep

    lower_rep = (
        E_rep / denom_rep
        if denom_rep
        else math.nan
    )

    combined_rep = (
        E_rep * (1.0 + 1.0 / r_protein) / denom_rep
        if denom_rep and r_protein > 0
        else math.nan
    )

    def status(lower_value, combined_value):
        if not math.isfinite(lower_value) or not math.isfinite(combined_value):
            return "NOT_EVALUABLE"

        if combined_value <= nominal:
            return "EVIDENCE_CONSISTENT_WITH_CONTROL"

        if lower_value > nominal:
            return "EVIDENCE_SUGGESTING_FAILURE"

        return "INCONCLUSIVE_BOUNDS_STRADDLE_NOMINAL"

    summary = {
        "workflow": workflow,
        "reported_groups_after_native_filtering": len(rows),

        "target_only_groups": T,
        "entrapment_only_groups": E,

        "mixed_target_entrapment_groups":
            primary_counts["MIXED_TARGET_ENTRAPMENT"],

        "target_plus_unresolved_groups":
            primary_counts["TARGET_PLUS_UNRESOLVED"],

        "entrapment_plus_unresolved_groups":
            primary_counts["ENTRAPMENT_PLUS_UNRESOLVED"],

        "unresolved_groups":
            primary_counts["UNRESOLVED"],

        "primary_evaluable_T_plus_E": denom,

        "r_protein": r_protein,
        "nominal_protein_fdr": nominal,

        "primary_lower_bound_fdp": lower,
        "primary_combined_fdp": combined,

        "primary_lower_bound_pct":
            100 * lower if math.isfinite(lower) else math.nan,

        "primary_combined_fdp_pct":
            100 * combined if math.isfinite(combined) else math.nan,

        "primary_interpretation": status(lower, combined),

        "representative_target": T_rep,
        "representative_entrapment": E_rep,
        "representative_unresolved":
            rep_counts["UNRESOLVED"] + rep_counts["COLLISION"],

        "representative_lower_bound_pct":
            100 * lower_rep if math.isfinite(lower_rep) else math.nan,

        "representative_combined_fdp_pct":
            100 * combined_rep if math.isfinite(combined_rep) else math.nan,

        "representative_interpretation":
            status(lower_rep, combined_rep),
    }

    return summary, detail_rows

# =============================================================================
# Main
# =============================================================================

def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--root",
        default=str(DEFAULT_ROOT),
        help="ENTRAPMENT_FDR root folder",
    )

    ap.add_argument(
        "--alpha",
        type=float,
        default=0.01,
        help="Nominal native protein FDR (default 0.01)",
    )

    ap.add_argument(
        "--outdir",
        default=None,
        help="Optional output folder",
    )

    args = ap.parse_args()

    root = Path(args.root)
    nominal = float(args.alpha)

    if not root.exists():
        raise SystemExit(f"Root does not exist:\n{root}")

    outdir = (
        Path(args.outdir)
        if args.outdir
        else root / "Step02D_PROTEIN_ENTRAPMENT_FDR_v103"
    )

    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 100)
    print(f"STEP02D PROTEIN-LEVEL ENTRAPMENT FDR v{VERSION}")
    print("=" * 100)
    print("Root:", root)
    print("Nominal protein FDR:", f"{100*nominal:.2f}%")

    # ---------------------------------------------------------------------
    # Frozen database
    # ---------------------------------------------------------------------

    step02a = locate_step02a(root)

    fasta = (
        step02a
        / "Step02_target_plus_shuffled_entrapment_r1.fasta"
    )

    manifest_path = step02a / "Step02A_manifest.json"

    if not fasta.exists():
        raise FileNotFoundError(f"Frozen FASTA not found:\n{fasta}")

    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found:\n{manifest_path}")

    manifest = json.loads(
        manifest_path.read_text(encoding="utf-8")
    )

    pc = manifest.get("protein_counts", {})

    target_proteins = int(pc.get("target", -1))
    entrapment_proteins = int(pc.get("entrapment", -1))

    print("\n[1/5] Frozen protein-level entrapment database QC")
    print("Target proteins:", f"{target_proteins:,}")
    print("Entrapment proteins:", f"{entrapment_proteins:,}")

    if (
        target_proteins != EXPECTED_STEP02A["target_proteins"]
        or entrapment_proteins != EXPECTED_STEP02A["entrapment_proteins"]
    ):
        raise SystemExit(
            "Frozen Step02A protein counts do not match v103/v104 QC."
        )

    r_protein = entrapment_proteins / target_proteins

    print("Protein-space r =", f"{r_protein:.12f}")

    # ---------------------------------------------------------------------
    # Build FASTA aliases
    # ---------------------------------------------------------------------

    print("\n[2/5] Building frozen FASTA protein alias dictionary")

    alias_index, observed_counts = build_fasta_alias_index(fasta)

    print("Aliases:", f"{len(alias_index):,}")
    print("Observed target headers:", f"{observed_counts['TARGET']:,}")
    print("Observed entrapment headers:", f"{observed_counts['ENTRAPMENT']:,}")

    if (
        observed_counts["TARGET"] != target_proteins
        or observed_counts["ENTRAPMENT"] != entrapment_proteins
    ):
        raise SystemExit(
            "Header classification does not reproduce frozen protein counts."
        )

    # ---------------------------------------------------------------------
    # Locate workflows
    # ---------------------------------------------------------------------

    print("\n[3/5] Locating completed reruns")

    workflow_dirs = {}

    for wf in WORKFLOWS:
        workflow_dirs[wf] = locate_workflow_dir(root, wf)
        print(f"{wf}: {workflow_dirs[wf]}")

    # ---------------------------------------------------------------------
    # Parse/evaluate
    # ---------------------------------------------------------------------

    print("\n[4/5] Parsing native protein/protein-group outputs")

    parsers = {
        "AP": parse_ap,
        "FP": parse_fp,
        "MM": parse_mm,
        "MQ": parse_mq,
    }

    summary_rows = []
    detail_rows = []
    source_rows = []

    for wf in WORKFLOWS:
        print(f"\n[{wf}]")

        rows, source_files, parser_meta = parsers[wf](
            workflow_dirs[wf],
            nominal,
        )

        if source_files:
            print("source:", source_files[0])
        print("parser:", json.dumps(parser_meta, sort_keys=True))

        summary, details = evaluate_workflow(
            wf,
            rows,
            alias_index,
            r_protein,
            nominal,
        )

        summary_rows.append(summary)
        detail_rows.extend(details)

        for source in source_files:
            source_rows.append({
                "workflow": wf,
                "source_file": str(source),
                "size_mb": round(
                    source.stat().st_size / 1024 / 1024,
                    4,
                ),
            })

        source_rows.append({
            "workflow": wf,
            "source_file":
                "PARSER_METADATA="
                + json.dumps(
                    parser_meta,
                    sort_keys=True,
                ),
            "size_mb": "",
        })

        print(
            f"reported groups={summary['reported_groups_after_native_filtering']:,}"
        )

        print(
            "PURE groups: "
            f"T={summary['target_only_groups']:,} | "
            f"E={summary['entrapment_only_groups']:,} | "
            f"mixed={summary['mixed_target_entrapment_groups']:,} | "
            f"unresolved={summary['unresolved_groups']:,}"
        )

        print(
            "Protein FDP: "
            f"lower={summary['primary_lower_bound_pct']:.4f}% | "
            f"combined={summary['primary_combined_fdp_pct']:.4f}% | "
            f"{summary['primary_interpretation']}"
        )

        print(
            "Representative sensitivity: "
            f"lower={summary['representative_lower_bound_pct']:.4f}% | "
            f"combined={summary['representative_combined_fdp_pct']:.4f}% | "
            f"{summary['representative_interpretation']}"
        )

    # ---------------------------------------------------------------------
    # Save
    # ---------------------------------------------------------------------

    print("\n[5/5] Writing frozen evaluation tables")

    summary_fields = [
        "workflow",
        "reported_groups_after_native_filtering",
        "target_only_groups",
        "entrapment_only_groups",
        "mixed_target_entrapment_groups",
        "target_plus_unresolved_groups",
        "entrapment_plus_unresolved_groups",
        "unresolved_groups",
        "primary_evaluable_T_plus_E",
        "r_protein",
        "nominal_protein_fdr",
        "primary_lower_bound_fdp",
        "primary_combined_fdp",
        "primary_lower_bound_pct",
        "primary_combined_fdp_pct",
        "primary_interpretation",
        "representative_target",
        "representative_entrapment",
        "representative_unresolved",
        "representative_lower_bound_pct",
        "representative_combined_fdp_pct",
        "representative_interpretation",
    ]

    detail_fields = [
        "workflow",
        "native_group_id",
        "primary_group_class",
        "representative_class",
        "representative",
        "n_members",
        "members",
        "member_classes",
        "native_q_value",
        "native_row_text",
    ]

    source_fields = [
        "workflow",
        "source_file",
        "size_mb",
    ]

    write_tsv(
        outdir / "Step02D_protein_entrapment_summary.tsv",
        summary_rows,
        summary_fields,
    )

    write_tsv(
        outdir / "Step02D_protein_group_classification.tsv",
        detail_rows,
        detail_fields,
    )

    write_tsv(
        outdir / "Step02D_source_files.tsv",
        source_rows,
        source_fields,
    )

    report_lines = [
        f"Step02D PROTEIN-LEVEL ENTRAPMENT FDR v{VERSION}",
        f"Nominal protein FDR: {100*nominal:.2f}%",
        f"Protein-space r: {r_protein:.12f}",
        "",
        (
            "workflow\treported_groups\tT_only\tE_only\tmixed\t"
            "unresolved\tlower_FDP_pct\tcombined_FDP_pct\tinterpretation"
        ),
    ]

    for row in summary_rows:
        report_lines.append(
            f"{row['workflow']}\t"
            f"{row['reported_groups_after_native_filtering']}\t"
            f"{row['target_only_groups']}\t"
            f"{row['entrapment_only_groups']}\t"
            f"{row['mixed_target_entrapment_groups']}\t"
            f"{row['unresolved_groups']}\t"
            f"{row['primary_lower_bound_pct']:.6f}\t"
            f"{row['primary_combined_fdp_pct']:.6f}\t"
            f"{row['primary_interpretation']}"
        )

    report_path = outdir / "Step02D_PROTEIN_FDR_REPORT.txt"

    report_path.write_text(
        "\n".join(report_lines),
        encoding="utf-8",
    )

    method_record = {
        "script_version": VERSION,
        "root": str(root),
        "frozen_fasta": str(fasta),
        "target_proteins": target_proteins,
        "entrapment_proteins": entrapment_proteins,
        "r_protein": r_protein,
        "nominal_protein_fdr": nominal,
        "primary_unit": "reported protein group",
        "primary_group_rule":
            "TARGET_ONLY and ENTRAPMENT_ONLY pure groups only; "
            "MIXED and UNRESOLVED excluded from T/E FDP denominator",
        "lower_fdp_formula": "E / (T + E)",
        "combined_fdp_formula":
            "E * (1 + 1/r_protein) / (T + E)",
        "representative_analysis":
            "sensitivity analysis only",
    }

    (
        outdir / "Step02D_method_record.json"
    ).write_text(
        json.dumps(method_record, indent=2),
        encoding="utf-8",
    )

    print("\n" + "=" * 100)
    print("STEP02D COMPLETE")
    print("=" * 100)

    for row in summary_rows:
        print(
            f"{row['workflow']}: "
            f"T={row['target_only_groups']:,} | "
            f"E={row['entrapment_only_groups']:,} | "
            f"lower={row['primary_lower_bound_pct']:.4f}% | "
            f"combined={row['primary_combined_fdp_pct']:.4f}% | "
            f"{row['primary_interpretation']}"
        )

    print("\nSend back:")
    print(report_path)
    print(outdir / "Step02D_protein_entrapment_summary.tsv")
    print(outdir / "Step02D_source_files.tsv")
    print("=" * 100)

if __name__ == "__main__":
    main()
