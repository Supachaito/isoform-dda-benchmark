#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step02E_SALVAGE_MM_MQ_v1.0.0.py

Goal
----
Determine whether the EXISTING MetaMorpheus (MM) and MaxQuant (MQ) entrapment
reruns can still be used for protein-level entrapment FDR WITHOUT rerunning
the searches.

Why this is needed
------------------
The frozen entrapment FASTA has headers like:

    >tr|A0A087WVL8|A0A087WVL8_HUMAN_p_target

The middle UniProt accession is unchanged, while the third/header-name field
contains "_p_target". If a search engine reports only the accession, target
and entrapment become indistinguishable. If its output retains FASTA headers,
protein names, descriptions, or another field containing "_p_target", the
existing search can potentially be salvaged.

What this script does
---------------------
1) Scans MM_ENTRAP and MQ_ENTRAP recursively.
2) Inspects text/tabular result files (.tsv/.txt/.csv).
3) Reports every file/column containing "_p_target".
4) Specifically audits the most likely protein-level tables:
      MM: protein/group-like tables
      MQ: proteinGroups.txt
5) For each accepted protein-group row, determines whether the row contains:
      - direct entrapment marker "_p_target"
      - full FASTA/header-like information
      - only accession-like identifiers
6) Produces a verdict for each workflow:
      SALVAGEABLE_DIRECT_MARKER
      POSSIBLY_SALVAGEABLE_NEEDS_ROW_MAPPING
      NOT_SALVAGEABLE_FROM_PROTEIN_TABLE
      NO_PROTEIN_TABLE_FOUND

IMPORTANT
---------
This script does NOT recalculate FDR and does NOT modify any search output.
It only determines whether the existing MM/MQ outputs retain enough identity
information to distinguish target from entrapment proteins.

Default root:
C:\\Users\\Supachai\\Desktop\\AphaPept_benchmark\\Benchmark_Program\\ENTRAPMENT_FDR
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
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


VERSION = "1.0.0"

DEFAULT_ROOT = _public_project_root() / "ENTRAPMENT_FDR"

MARKER = "_p_target"
TEXT_SUFFIXES = {".tsv", ".txt", ".csv"}

# -----------------------------------------------------------------------------
# CSV robustness
# -----------------------------------------------------------------------------

def _raise_csv_field_limit():
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return limit
        except OverflowError:
            limit //= 10

CSV_FIELD_LIMIT = _raise_csv_field_limit()

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def norm_col(x: str) -> str:
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

def detect_delimiter(path: Path) -> str:
    if path.suffix.lower() == ".csv":
        return ","
    return "\t"

def read_dict_rows(path: Path):
    delim = detect_delimiter(path)
    with path.open(
        "r",
        encoding="utf-8-sig",
        errors="replace",
        newline=""
    ) as fh:
        yield from csv.DictReader(fh, delimiter=delim)

def locate_workflow_dir(root: Path, workflow: str) -> Path:
    exact = root / f"{workflow}_ENTRAP"
    if exact.exists():
        return exact

    hits = [
        p for p in root.rglob("*")
        if p.is_dir()
        and p.name.lower() == f"{workflow}_entrap".lower()
    ]
    if not hits:
        raise FileNotFoundError(
            f"Could not find {workflow}_ENTRAP under:\n{root}"
        )
    return sorted(hits, key=lambda p: len(str(p)))[0]

def safe_rel(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except Exception:
        return str(path)

def file_marker_count(path: Path, marker: str = MARKER, max_bytes=None):
    """
    Count lines containing marker. Streaming; no entire file in memory.
    """
    n_lines = 0
    n_marker_lines = 0
    samples = []

    try:
        with path.open(
            "r",
            encoding="utf-8-sig",
            errors="replace"
        ) as fh:
            for line in fh:
                n_lines += 1
                if marker.lower() in line.lower():
                    n_marker_lines += 1
                    if len(samples) < 3:
                        samples.append(line.rstrip("\r\n")[:1000])
    except Exception as e:
        return {
            "ok": False,
            "error": str(e),
            "n_lines": n_lines,
            "n_marker_lines": n_marker_lines,
            "samples": samples,
        }

    return {
        "ok": True,
        "error": "",
        "n_lines": n_lines,
        "n_marker_lines": n_marker_lines,
        "samples": samples,
    }

def inspect_table_marker_columns(path: Path, max_rows=200000):
    """
    Identify columns that contain _p_target and count marked rows/values.
    """
    result = {
        "path": str(path),
        "fields": [],
        "rows_scanned": 0,
        "rows_with_marker_anywhere": 0,
        "marker_columns": Counter(),
        "marker_samples": {},
        "error": "",
    }

    try:
        rows = read_dict_rows(path)
        first = next(rows, None)
        if first is None:
            return result

        fields = list(first.keys())
        result["fields"] = fields

        def process(row):
            result["rows_scanned"] += 1
            row_has = False

            for col, val in row.items():
                text = str(val or "")
                if MARKER.lower() in text.lower():
                    result["marker_columns"][col] += 1
                    row_has = True
                    result["marker_samples"].setdefault(col, [])
                    if len(result["marker_samples"][col]) < 3:
                        result["marker_samples"][col].append(text[:1000])

            if row_has:
                result["rows_with_marker_anywhere"] += 1

        process(first)

        for row in rows:
            process(row)
            if result["rows_scanned"] >= max_rows:
                break

    except Exception as e:
        result["error"] = str(e)

    result["marker_columns"] = dict(result["marker_columns"])
    return result

# -----------------------------------------------------------------------------
# Candidate protein tables
# -----------------------------------------------------------------------------

def score_mm_candidate(path: Path, fields) -> int:
    name = path.name.lower()
    score = 0

    if "allproteingroups" in name:
        score += 120
    elif "proteingroups" in name:
        score += 100
    elif "protein" in name and "group" in name:
        score += 80
    elif "protein" in name:
        score += 30

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
        ]
    )
    if protein_col:
        score += 60

    q_col = find_col(
        fields,
        ["QValue", "Q Value", "q-value", "q_value", "qvalue"]
    )
    if q_col:
        score += 15

    dct_col = find_col(
        fields,
        [
            "Decoy/Contaminant/Target",
            "Decoy Contaminant Target",
            "DCT",
        ]
    )
    if dct_col:
        score += 10

    return score

def choose_mm_protein_table(mm_dir: Path):
    candidates = []

    for path in mm_dir.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES:
            continue

        lname = path.name.lower()
        if "protein" not in lname and "group" not in lname:
            continue

        try:
            rows = read_dict_rows(path)
            first = next(rows, None)
        except Exception:
            continue

        if first is None:
            continue

        fields = list(first.keys())
        score = score_mm_candidate(path, fields)

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
            ]
        )

        if protein_col:
            candidates.append(
                {
                    "path": path,
                    "fields": fields,
                    "score": score,
                    "protein_col": protein_col,
                    "q_col": find_col(
                        fields,
                        ["QValue", "Q Value", "q-value", "q_value", "qvalue"]
                    ),
                    "dct_col": find_col(
                        fields,
                        [
                            "Decoy/Contaminant/Target",
                            "Decoy Contaminant Target",
                            "DCT",
                        ]
                    ),
                }
            )

    if not candidates:
        return None, []

    candidates.sort(
        key=lambda x: (
            -x["score"],
            -x["path"].stat().st_size,
            len(str(x["path"])),
        )
    )

    return candidates[0], candidates

def choose_mq_protein_table(mq_dir: Path):
    hits = sorted(mq_dir.rglob("proteinGroups.txt"))
    if not hits:
        return None

    hits.sort(
        key=lambda p: (
            0 if p.parent.name.lower() == "txt" else 1,
            len(str(p)),
        )
    )
    return hits[0]

# -----------------------------------------------------------------------------
# Row-level salvage audit
# -----------------------------------------------------------------------------

HEADER_LIKE_ALIASES = [
    "Fasta headers",
    "FASTA headers",
    "Fasta header",
    "FASTA header",
    "Protein names",
    "Protein name",
    "Description",
    "Descriptions",
    "Full sequence",
    "Protein Full Name",
    "Protein Name",
]

def audit_mm_table(path: Path, alpha=0.01):
    rows = read_dict_rows(path)
    first = next(rows, None)
    if first is None:
        raise RuntimeError(f"Empty MM protein table: {path}")

    fields = list(first.keys())

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
        ]
    )
    q_col = find_col(
        fields,
        ["QValue", "Q Value", "q-value", "q_value", "qvalue"]
    )
    dct_col = find_col(
        fields,
        [
            "Decoy/Contaminant/Target",
            "Decoy Contaminant Target",
            "DCT",
        ]
    )

    header_like_cols = []
    for alias in HEADER_LIKE_ALIASES:
        col = find_col(fields, [alias])
        if col and col not in header_like_cols:
            header_like_cols.append(col)

    # Also discover any column whose name itself suggests FASTA/name/description.
    for col in fields:
        n = norm_col(col)
        if any(k in n for k in ["fasta", "header", "description", "proteinname", "fullname"]):
            if col not in header_like_cols:
                header_like_cols.append(col)

    stats = Counter()
    marker_cols = Counter()
    samples = []

    def process(row, row_no):
        q = float_or_none(row.get(q_col)) if q_col else None
        if q is not None and q > alpha:
            stats["excluded_q"] += 1
            return

        if dct_col:
            dct = str(row.get(dct_col, "")).strip().upper()
            if dct.startswith("D") or dct.startswith("C"):
                stats["excluded_decoy_contaminant"] += 1
                return

        stats["accepted_rows"] += 1

        protein_text = str(row.get(protein_col, "") or "")
        row_marker_cols = []

        for col, val in row.items():
            text = str(val or "")
            if MARKER.lower() in text.lower():
                marker_cols[col] += 1
                row_marker_cols.append(col)

        if row_marker_cols:
            stats["accepted_rows_with_direct_marker"] += 1

            if len(samples) < 10:
                samples.append(
                    {
                        "row": row_no,
                        "protein_field": protein_text[:500],
                        "marker_columns": row_marker_cols,
                        "marker_values": {
                            c: str(row.get(c, ""))[:1000]
                            for c in row_marker_cols[:5]
                        }
                    }
                )

        if MARKER.lower() in protein_text.lower():
            stats["marker_in_primary_protein_column"] += 1

        header_marker = False
        for col in header_like_cols:
            if MARKER.lower() in str(row.get(col, "") or "").lower():
                header_marker = True
                break
        if header_marker:
            stats["marker_in_header_like_column"] += 1

    process(first, 1)
    for i, row in enumerate(rows, start=2):
        process(row, i)

    accepted = stats["accepted_rows"]
    marked = stats["accepted_rows_with_direct_marker"]

    if marked > 0:
        verdict = "SALVAGEABLE_DIRECT_MARKER"
    elif header_like_cols:
        verdict = "NOT_SALVAGEABLE_FROM_PROTEIN_TABLE"
    else:
        verdict = "NOT_SALVAGEABLE_FROM_PROTEIN_TABLE"

    return {
        "workflow": "MM",
        "protein_table": str(path),
        "fields": fields,
        "primary_protein_column": protein_col,
        "q_column": q_col or "",
        "dct_column": dct_col or "",
        "header_like_columns": header_like_cols,
        "accepted_rows": accepted,
        "accepted_rows_with_direct_marker": marked,
        "marker_in_primary_protein_column": stats["marker_in_primary_protein_column"],
        "marker_in_header_like_column": stats["marker_in_header_like_column"],
        "marker_columns": dict(marker_cols),
        "samples": samples,
        "verdict": verdict,
    }

def audit_mq_table(path: Path):
    rows = read_dict_rows(path)
    first = next(rows, None)
    if first is None:
        raise RuntimeError(f"Empty MQ proteinGroups.txt: {path}")

    fields = list(first.keys())

    protein_ids_col = find_col(fields, ["Protein IDs"])
    majority_col = find_col(fields, ["Majority protein IDs"])
    leading_col = find_col(fields, ["Leading razor protein"])
    reverse_col = find_col(fields, ["Reverse"])
    contaminant_col = find_col(
        fields,
        ["Potential contaminant", "Potential Contaminant"]
    )

    header_like_cols = []
    for alias in HEADER_LIKE_ALIASES:
        col = find_col(fields, [alias])
        if col and col not in header_like_cols:
            header_like_cols.append(col)

    for col in fields:
        n = norm_col(col)
        if any(k in n for k in ["fasta", "header", "description", "proteinname", "fullname"]):
            if col not in header_like_cols:
                header_like_cols.append(col)

    stats = Counter()
    marker_cols = Counter()
    samples = []

    def process(row, row_no):
        if reverse_col and truthy(row.get(reverse_col)):
            stats["excluded_reverse"] += 1
            return

        if contaminant_col and truthy(row.get(contaminant_col)):
            stats["excluded_contaminant"] += 1
            return

        stats["accepted_rows"] += 1

        row_marker_cols = []
        for col, val in row.items():
            text = str(val or "")
            if MARKER.lower() in text.lower():
                marker_cols[col] += 1
                row_marker_cols.append(col)

        if row_marker_cols:
            stats["accepted_rows_with_direct_marker"] += 1
            if len(samples) < 10:
                samples.append(
                    {
                        "row": row_no,
                        "protein_ids": str(row.get(protein_ids_col, ""))[:500] if protein_ids_col else "",
                        "majority_ids": str(row.get(majority_col, ""))[:500] if majority_col else "",
                        "leading": str(row.get(leading_col, ""))[:500] if leading_col else "",
                        "marker_columns": row_marker_cols,
                        "marker_values": {
                            c: str(row.get(c, ""))[:1000]
                            for c in row_marker_cols[:5]
                        },
                    }
                )

        for col_name, stat_name in [
            (protein_ids_col, "marker_in_protein_ids"),
            (majority_col, "marker_in_majority_ids"),
            (leading_col, "marker_in_leading_razor"),
        ]:
            if col_name and MARKER.lower() in str(row.get(col_name, "") or "").lower():
                stats[stat_name] += 1

        header_marker = False
        for col in header_like_cols:
            if MARKER.lower() in str(row.get(col, "") or "").lower():
                header_marker = True
                break
        if header_marker:
            stats["marker_in_header_like_column"] += 1

    process(first, 1)
    for i, row in enumerate(rows, start=2):
        process(row, i)

    marked = stats["accepted_rows_with_direct_marker"]

    if marked > 0:
        verdict = "SALVAGEABLE_DIRECT_MARKER"
    elif header_like_cols:
        verdict = "NOT_SALVAGEABLE_FROM_PROTEIN_TABLE"
    else:
        verdict = "NOT_SALVAGEABLE_FROM_PROTEIN_TABLE"

    return {
        "workflow": "MQ",
        "protein_table": str(path),
        "fields": fields,
        "protein_ids_column": protein_ids_col or "",
        "majority_ids_column": majority_col or "",
        "leading_razor_column": leading_col or "",
        "reverse_column": reverse_col or "",
        "contaminant_column": contaminant_col or "",
        "header_like_columns": header_like_cols,
        "accepted_rows": stats["accepted_rows"],
        "accepted_rows_with_direct_marker": marked,
        "marker_in_protein_ids": stats["marker_in_protein_ids"],
        "marker_in_majority_ids": stats["marker_in_majority_ids"],
        "marker_in_leading_razor": stats["marker_in_leading_razor"],
        "marker_in_header_like_column": stats["marker_in_header_like_column"],
        "marker_columns": dict(marker_cols),
        "samples": samples,
        "verdict": verdict,
    }

# -----------------------------------------------------------------------------
# Scan all text files for marker
# -----------------------------------------------------------------------------

def scan_workflow_files(workflow_dir: Path, workflow: str):
    rows = []

    files = [
        p for p in workflow_dir.rglob("*")
        if p.is_file()
        and p.suffix.lower() in TEXT_SUFFIXES
    ]

    for i, path in enumerate(sorted(files), start=1):
        result = file_marker_count(path)

        rows.append({
            "workflow": workflow,
            "file": str(path),
            "relative_file": safe_rel(path, workflow_dir),
            "size_mb": round(path.stat().st_size / (1024 * 1024), 4),
            "n_lines": result["n_lines"],
            "marker_lines": result["n_marker_lines"],
            "contains_p_target": 1 if result["n_marker_lines"] > 0 else 0,
            "sample_1": result["samples"][0] if len(result["samples"]) > 0 else "",
            "sample_2": result["samples"][1] if len(result["samples"]) > 1 else "",
            "sample_3": result["samples"][2] if len(result["samples"]) > 2 else "",
            "read_error": result["error"],
        })

        if i % 25 == 0:
            print(f"  {workflow}: scanned {i}/{len(files)} text files")

    return rows

# -----------------------------------------------------------------------------
# Writing
# -----------------------------------------------------------------------------

def write_tsv(path: Path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--root",
        default=str(DEFAULT_ROOT),
        help="ENTRAPMENT_FDR root"
    )

    ap.add_argument(
        "--alpha",
        type=float,
        default=0.01,
        help="MM protein QValue cutoff for accepted-row audit"
    )

    ap.add_argument(
        "--outdir",
        default=None,
        help="Optional output directory"
    )

    args = ap.parse_args()

    root = Path(args.root)
    if not root.exists():
        raise SystemExit(f"Root does not exist:\n{root}")

    outdir = (
        Path(args.outdir)
        if args.outdir
        else root / "Step02E_SALVAGE_MM_MQ_v100"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    mm_dir = locate_workflow_dir(root, "MM")
    mq_dir = locate_workflow_dir(root, "MQ")

    print("=" * 100)
    print(f"STEP02E MM/MQ SALVAGE CHECK v{VERSION}")
    print("=" * 100)
    print("Root:", root)
    print("MM:", mm_dir)
    print("MQ:", mq_dir)
    print("Marker:", MARKER)

    # -------------------------------------------------------------------------
    # 1. Whole-folder marker scan
    # -------------------------------------------------------------------------

    print("\n[1/4] Scanning MM/MQ text outputs for '_p_target'")

    file_scan = []
    file_scan.extend(scan_workflow_files(mm_dir, "MM"))
    file_scan.extend(scan_workflow_files(mq_dir, "MQ"))

    marked_files = [
        r for r in file_scan
        if r["contains_p_target"] == 1
    ]

    print(f"Files containing marker: {len(marked_files):,}")

    for r in marked_files[:20]:
        print(
            f"  {r['workflow']} | marker_lines={r['marker_lines']:,} | "
            f"{r['relative_file']}"
        )

    # -------------------------------------------------------------------------
    # 2. Select protein tables
    # -------------------------------------------------------------------------

    print("\n[2/4] Selecting native protein-level tables")

    mm_choice, mm_candidates = choose_mm_protein_table(mm_dir)
    mq_choice = choose_mq_protein_table(mq_dir)

    if mm_choice:
        print("MM selected:", mm_choice["path"])
        print("MM protein column:", mm_choice["protein_col"])
        print("MM Q column:", mm_choice["q_col"])
        print("MM DCT column:", mm_choice["dct_col"])
    else:
        print("MM: NO protein-group table found")

    if mq_choice:
        print("MQ selected:", mq_choice)
    else:
        print("MQ: NO proteinGroups.txt found")

    # -------------------------------------------------------------------------
    # 3. Row-level salvage
    # -------------------------------------------------------------------------

    print("\n[3/4] Auditing whether accepted protein groups retain entrapment identity")

    results = {}

    if mm_choice:
        results["MM"] = audit_mm_table(
            mm_choice["path"],
            alpha=args.alpha
        )
        print(
            "MM:",
            results["MM"]["verdict"],
            "| accepted=",
            f"{results['MM']['accepted_rows']:,}",
            "| rows_with_marker=",
            f"{results['MM']['accepted_rows_with_direct_marker']:,}",
        )
        print("MM marker columns:", results["MM"]["marker_columns"])
    else:
        results["MM"] = {
            "workflow": "MM",
            "verdict": "NO_PROTEIN_TABLE_FOUND",
        }

    if mq_choice:
        results["MQ"] = audit_mq_table(mq_choice)
        print(
            "MQ:",
            results["MQ"]["verdict"],
            "| accepted=",
            f"{results['MQ']['accepted_rows']:,}",
            "| rows_with_marker=",
            f"{results['MQ']['accepted_rows_with_direct_marker']:,}",
        )
        print("MQ marker columns:", results["MQ"]["marker_columns"])
    else:
        results["MQ"] = {
            "workflow": "MQ",
            "verdict": "NO_PROTEIN_TABLE_FOUND",
        }

    # -------------------------------------------------------------------------
    # 4. Write reports
    # -------------------------------------------------------------------------

    print("\n[4/4] Writing reports")

    scan_fields = [
        "workflow",
        "file",
        "relative_file",
        "size_mb",
        "n_lines",
        "marker_lines",
        "contains_p_target",
        "sample_1",
        "sample_2",
        "sample_3",
        "read_error",
    ]

    write_tsv(
        outdir / "Step02E_all_text_files_marker_scan.tsv",
        file_scan,
        scan_fields,
    )

    summary_rows = []

    for wf in ["MM", "MQ"]:
        r = results[wf]
        summary_rows.append({
            "workflow": wf,
            "verdict": r.get("verdict", ""),
            "protein_table": r.get("protein_table", ""),
            "accepted_rows": r.get("accepted_rows", ""),
            "accepted_rows_with_direct_marker":
                r.get("accepted_rows_with_direct_marker", ""),
            "header_like_columns":
                " ; ".join(r.get("header_like_columns", []))
                if isinstance(r.get("header_like_columns", []), list)
                else str(r.get("header_like_columns", "")),
            "marker_columns_json":
                json.dumps(r.get("marker_columns", {}), sort_keys=True),
        })

    write_tsv(
        outdir / "Step02E_salvage_summary.tsv",
        summary_rows,
        [
            "workflow",
            "verdict",
            "protein_table",
            "accepted_rows",
            "accepted_rows_with_direct_marker",
            "header_like_columns",
            "marker_columns_json",
        ],
    )

    (outdir / "Step02E_MM_detail.json").write_text(
        json.dumps(results["MM"], indent=2, ensure_ascii=False),
        encoding="utf-8"
    )

    (outdir / "Step02E_MQ_detail.json").write_text(
        json.dumps(results["MQ"], indent=2, ensure_ascii=False),
        encoding="utf-8"
    )

    report_lines = [
        f"Step02E MM/MQ SALVAGE CHECK v{VERSION}",
        "",
        "Purpose: determine whether existing MM/MQ protein-level outputs retain",
        "direct _p_target identity information from the original entrapment FASTA.",
        "",
    ]

    for wf in ["MM", "MQ"]:
        r = results[wf]

        report_lines.append(
            f"{wf}: {r.get('verdict', '')}"
        )

        if "accepted_rows" in r:
            report_lines.append(
                f"  accepted protein-group rows: {r.get('accepted_rows', 0)}"
            )
            report_lines.append(
                "  accepted rows containing direct _p_target marker: "
                f"{r.get('accepted_rows_with_direct_marker', 0)}"
            )
            report_lines.append(
                "  marker columns: "
                f"{json.dumps(r.get('marker_columns', {}), sort_keys=True)}"
            )
            report_lines.append(
                "  header-like columns: "
                f"{'; '.join(r.get('header_like_columns', []))}"
            )

        report_lines.append("")

    report_lines.extend([
        "Interpretation rule:",
        "  SALVAGEABLE_DIRECT_MARKER",
        "    Existing protein-level output contains direct _p_target labels;",
        "    rerunning the search may be avoidable.",
        "",
        "  NOT_SALVAGEABLE_FROM_PROTEIN_TABLE",
        "    Protein-level output does not retain direct entrapment identity.",
        "    Inspect Step02E_all_text_files_marker_scan.tsv before deciding whether",
        "    lower-level native outputs can reconstruct the mapping.",
        "",
        "  NO_PROTEIN_TABLE_FOUND",
        "    Native protein-level table could not be located.",
    ])

    report_path = outdir / "Step02E_SALVAGE_REPORT.txt"
    report_path.write_text(
        "\n".join(report_lines),
        encoding="utf-8"
    )

    print("\n" + "=" * 100)
    print("STEP02E COMPLETE")
    print("=" * 100)

    for wf in ["MM", "MQ"]:
        print(f"{wf}: {results[wf].get('verdict', '')}")

    print("\nSend back:")
    print(report_path)
    print(outdir / "Step02E_salvage_summary.tsv")
    print(outdir / "Step02E_MM_detail.json")
    print(outdir / "Step02E_MQ_detail.json")
    print("=" * 100)


if __name__ == "__main__":
    main()
