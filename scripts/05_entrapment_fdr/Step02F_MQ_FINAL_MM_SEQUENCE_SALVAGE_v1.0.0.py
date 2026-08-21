#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step02F_MQ_FINAL_MM_SEQUENCE_SALVAGE_v1.0.0.py

Purpose
-------
Finish the existing protein-level entrapment analysis WITHOUT rerunning searches:

1) MaxQuant (MQ)
   Recalculate protein-group entrapment FDP directly from proteinGroups.txt
   using the "Fasta headers" column, which retains "_p_target".

2) MetaMorpheus (MM)
   Attempt a sequence-evidence reconstruction using the "Unique Peptides" and
   "Shared Peptides" columns from AllQuantifiedProteinGroups.tsv.
   Each peptide is mapped back to the frozen target+entrapment FASTA and
   classified as TARGET_ONLY, ENTRAPMENT_ONLY, SHARED, or ABSENT.

IMPORTANT SCIENTIFIC DISTINCTION
--------------------------------
MQ result from this script is a DIRECT protein-group entrapment evaluation.

MM result is a SEQUENCE-RECONSTRUCTION DIAGNOSTIC because the native MM protein
table has lost the "_p_target" identifier. It is useful to determine whether the
existing run can be interpreted and whether rerunning MM is necessary, but it
must not silently be presented as a direct native protein-group entrapment FDR.

The script also imports AP/FP results from Step02D v103 (if present) to create a
single combined status table.

Default root
------------
C:\\Users\\Supachai\\Desktop\\AphaPept_benchmark\\Benchmark_Program\\ENTRAPMENT_FDR

Outputs
-------
Step02F_MQ_FINAL_MM_SEQUENCE_SALVAGE_v100\\
    Step02F_MQ_direct_protein_FDR.tsv
    Step02F_MQ_group_classification.tsv
    Step02F_MM_sequence_reconstruction_summary.tsv
    Step02F_MM_group_sequence_classification.tsv
    Step02F_MM_peptide_mapping.tsv
    Step02F_combined_status.tsv
    Step02F_REPORT.txt
    Step02F_method_record.json
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import Counter, deque
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
MARKER = "_p_target"

DEFAULT_ROOT = _public_project_root() / "ENTRAPMENT_FDR"

EXPECTED_TARGET_PROTEINS = 169637
EXPECTED_ENTRAPMENT_PROTEINS = 169637

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
# Generic helpers
# -----------------------------------------------------------------------------

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

def read_delimited(path):
    delim = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open(
        "r",
        encoding="utf-8-sig",
        errors="replace",
        newline=""
    ) as fh:
        yield from csv.DictReader(fh, delimiter=delim)

def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

def locate_workflow_dir(root, workflow):
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
            f"Could not find Step02A database under:\n{root}"
        )
    return sorted(hits, key=lambda p: len(str(p)))[0]

def normalize_il(seq):
    return str(seq or "").upper().replace("I", "J").replace("L", "J")

def is_entrapment_header(header):
    return MARKER in str(header or "").lower()

# -----------------------------------------------------------------------------
# FASTA
# -----------------------------------------------------------------------------

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

def count_fasta_classes(fasta):
    target = 0
    entrap = 0

    for header, _seq in parse_fasta(fasta):
        if is_entrapment_header(header):
            entrap += 1
        else:
            target += 1

    return target, entrap

# -----------------------------------------------------------------------------
# FDP
# -----------------------------------------------------------------------------

def fdp_stats(T, E, r, nominal=0.01):
    denom = T + E

    if denom <= 0:
        return {
            "lower": math.nan,
            "combined": math.nan,
            "lower_pct": math.nan,
            "combined_pct": math.nan,
            "interpretation": "NOT_EVALUABLE",
        }

    lower = E / denom
    combined = E * (1.0 + 1.0 / r) / denom

    if combined <= nominal:
        interpretation = "EVIDENCE_CONSISTENT_WITH_CONTROL"
    elif lower > nominal:
        interpretation = "EVIDENCE_SUGGESTING_FAILURE"
    else:
        interpretation = "INCONCLUSIVE_BOUNDS_STRADDLE_NOMINAL"

    return {
        "lower": lower,
        "combined": combined,
        "lower_pct": 100 * lower,
        "combined_pct": 100 * combined,
        "interpretation": interpretation,
    }

# =============================================================================
# PART A: MaxQuant direct protein-group evaluation using Fasta headers
# =============================================================================

def locate_mq_protein_groups(mq_dir):
    hits = sorted(mq_dir.rglob("proteinGroups.txt"))
    if not hits:
        raise FileNotFoundError(
            f"No proteinGroups.txt under:\n{mq_dir}"
        )

    hits.sort(
        key=lambda p: (
            0 if p.parent.name.lower() == "txt" else 1,
            len(str(p)),
        )
    )
    return hits[0]

def split_fasta_headers(text):
    """
    MaxQuant Fasta headers are normally semicolon-separated.
    Keep each header intact so _p_target can be evaluated directly.
    """
    s = str(text or "").strip()
    if not s:
        return []

    return [
        x.strip()
        for x in s.split(";")
        if x.strip()
    ]

def evaluate_mq(mq_path, r_protein, nominal):
    rows = read_delimited(mq_path)
    first = next(rows, None)

    if first is None:
        raise RuntimeError(f"Empty MQ table: {mq_path}")

    fields = list(first.keys())

    fasta_col = find_col(fields, ["Fasta headers", "FASTA headers"])
    protein_ids_col = find_col(fields, ["Protein IDs"])
    q_col = find_col(fields, ["Q-value", "Q value", "QValue"])
    reverse_col = find_col(fields, ["Reverse", "Decoy"])
    contaminant_col = find_col(
        fields,
        ["Potential contaminant", "Potential Contaminant"]
    )

    if fasta_col is None:
        raise RuntimeError(
            "[MQ] Fasta headers column not found.\n"
            f"Columns: {fields}"
        )

    counts = Counter()
    details = []

    def process(row, row_no):
        if reverse_col and truthy(row.get(reverse_col)):
            counts["excluded_decoy"] += 1
            return

        if contaminant_col and truthy(row.get(contaminant_col)):
            counts["excluded_contaminant"] += 1
            return

        q = float_or_none(row.get(q_col)) if q_col else None
        if q is not None and q > nominal:
            counts["excluded_q"] += 1
            return

        headers = split_fasta_headers(row.get(fasta_col))

        if not headers:
            cls = "UNRESOLVED"
            n_target = 0
            n_entrap = 0
        else:
            n_entrap = sum(
                1 for h in headers
                if is_entrapment_header(h)
            )
            n_target = len(headers) - n_entrap

            if n_target > 0 and n_entrap == 0:
                cls = "TARGET_ONLY"
            elif n_entrap > 0 and n_target == 0:
                cls = "ENTRAPMENT_ONLY"
            elif n_target > 0 and n_entrap > 0:
                cls = "MIXED_TARGET_ENTRAPMENT"
            else:
                cls = "UNRESOLVED"

        counts[cls] += 1
        counts["accepted_rows"] += 1

        details.append({
            "workflow": "MQ",
            "row_number": row_no,
            "protein_ids": row.get(protein_ids_col, "") if protein_ids_col else "",
            "q_value": row.get(q_col, "") if q_col else "",
            "group_class": cls,
            "n_fasta_headers": len(headers),
            "n_target_headers": n_target,
            "n_entrapment_headers": n_entrap,
            "fasta_headers": row.get(fasta_col, ""),
        })

    process(first, 1)

    for i, row in enumerate(rows, start=2):
        process(row, i)

    T = counts["TARGET_ONLY"]
    E = counts["ENTRAPMENT_ONLY"]
    stats = fdp_stats(T, E, r_protein, nominal)

    summary = {
        "workflow": "MQ",
        "analysis_type": "DIRECT_PROTEIN_GROUP_FROM_FASTA_HEADERS",
        "source_file": str(mq_path),
        "reported_groups_after_native_filtering": counts["accepted_rows"],
        "target_only_groups": T,
        "entrapment_only_groups": E,
        "mixed_target_entrapment_groups": counts["MIXED_TARGET_ENTRAPMENT"],
        "unresolved_groups": counts["UNRESOLVED"],
        "primary_evaluable_T_plus_E": T + E,
        "r_protein": r_protein,
        "nominal_protein_fdr": nominal,
        "lower_bound_fdp": stats["lower"],
        "combined_fdp": stats["combined"],
        "lower_bound_pct": stats["lower_pct"],
        "combined_fdp_pct": stats["combined_pct"],
        "interpretation": stats["interpretation"],
        "excluded_decoy": counts["excluded_decoy"],
        "excluded_contaminant": counts["excluded_contaminant"],
        "excluded_q": counts["excluded_q"],
    }

    return summary, details

# =============================================================================
# PART B: MetaMorpheus sequence-evidence reconstruction
# =============================================================================

def locate_mm_protein_table(mm_dir):
    """
    Prefer AllQuantifiedProteinGroups.tsv because Step02E already identified it.
    """
    preferred = sorted(mm_dir.rglob("AllQuantifiedProteinGroups.tsv"))
    if preferred:
        return preferred[0]

    candidates = []
    for p in mm_dir.rglob("*"):
        if not p.is_file() or p.suffix.lower() not in {".tsv", ".txt", ".csv"}:
            continue

        lname = p.name.lower()
        if "protein" not in lname or "group" not in lname:
            continue

        try:
            rows = read_delimited(p)
            first = next(rows, None)
        except Exception:
            continue

        if first is None:
            continue

        fields = list(first.keys())
        pep_unique = find_col(fields, ["Unique Peptides"])
        accession = find_col(fields, ["Protein Accession", "Protein Accessions"])

        if accession and pep_unique:
            score = 0
            if "allquantifiedproteingroups" in lname:
                score += 100
            if "allproteingroups" in lname:
                score += 80
            candidates.append((score, p))

    if not candidates:
        raise FileNotFoundError(
            f"No suitable MM protein-group table under:\n{mm_dir}"
        )

    candidates.sort(
        key=lambda x: (
            -x[0],
            -x[1].stat().st_size,
            len(str(x[1])),
        )
    )

    return candidates[0][1]

AA = set("ACDEFGHIKLMNPQRSTVWYBJOUXZ")

def clean_peptide_token(token):
    """
    Convert a peptide field token to a base sequence.

    Handles common forms such as:
      PEPTIDE
      K.PEPTIDE.R
      PEPT[Oxidation on M]IDE
      PEPTIDE[+15.99]

    Returns "" if the token does not look like an amino-acid sequence.
    """
    s = str(token or "").strip()

    if not s:
        return ""

    # Remove flanking-residue notation K.PEPTIDE.R
    m = re.fullmatch(r"[A-Za-z\-]\.([^\.]+)\.[A-Za-z\-]", s)
    if m:
        s = m.group(1)

    # Remove annotation blocks.
    s = re.sub(r"\[[^\]]*\]", "", s)
    s = re.sub(r"\([^\)]*\)", "", s)
    s = re.sub(r"\{[^\}]*\}", "", s)

    # Keep letters only.
    s = re.sub(r"[^A-Za-z]", "", s).upper()

    if len(s) < 6:
        return ""

    if any(ch not in AA for ch in s):
        return ""

    return s

def extract_peptides_from_field(text):
    """
    Robust extraction from MetaMorpheus group peptide-list fields.

    First split on common list delimiters. If a chunk still contains multiple
    whitespace-separated sequence-like tokens, test them separately.
    """
    s = str(text or "").strip()

    if not s:
        return []

    coarse = re.split(r"[;,|]+", s)

    peptides = []

    for chunk in coarse:
        chunk = chunk.strip()

        if not chunk:
            continue

        # Try entire chunk first.
        whole = clean_peptide_token(chunk)

        if whole:
            peptides.append(whole)
            continue

        # Fallback to whitespace tokens.
        for tok in chunk.split():
            pep = clean_peptide_token(tok)
            if pep:
                peptides.append(pep)

    # Preserve order; remove duplicates.
    return list(dict.fromkeys(peptides))

def load_mm_groups_and_peptides(mm_path, nominal):
    rows = read_delimited(mm_path)
    first = next(rows, None)

    if first is None:
        raise RuntimeError(f"Empty MM table: {mm_path}")

    fields = list(first.keys())

    accession_col = find_col(
        fields,
        ["Protein Accession", "Protein Accessions"]
    )
    unique_col = find_col(fields, ["Unique Peptides"])
    shared_col = find_col(fields, ["Shared Peptides"])
    q_col = find_col(
        fields,
        ["Protein QValue", "QValue", "Q Value", "q-value"]
    )
    dct_col = find_col(
        fields,
        [
            "Protein Decoy/Contaminant/Target",
            "Decoy/Contaminant/Target",
            "DCT",
        ]
    )

    if accession_col is None:
        raise RuntimeError(
            "[MM] Protein Accession column not found."
        )

    if unique_col is None and shared_col is None:
        raise RuntimeError(
            "[MM] Neither Unique Peptides nor Shared Peptides column found.\n"
            f"Columns: {fields}"
        )

    groups = []
    all_peptides = set()
    extraction_audit = []
    counters = Counter()

    def process(row, row_no):
        q = float_or_none(row.get(q_col)) if q_col else None

        if q is not None and q > nominal:
            counters["excluded_q"] += 1
            return

        if dct_col:
            dct = str(row.get(dct_col, "") or "").strip().upper()
            if dct.startswith("D") or dct.startswith("C"):
                counters["excluded_decoy_contaminant"] += 1
                return

        unique_peps = (
            extract_peptides_from_field(row.get(unique_col))
            if unique_col else []
        )
        shared_peps = (
            extract_peptides_from_field(row.get(shared_col))
            if shared_col else []
        )

        combined = list(dict.fromkeys(unique_peps + shared_peps))

        counters["accepted_rows"] += 1
        counters["extracted_peptides_total"] += len(combined)

        if not combined:
            counters["accepted_rows_no_extracted_peptide"] += 1

        all_peptides.update(combined)

        groups.append({
            "row_number": row_no,
            "protein_accession": row.get(accession_col, ""),
            "q_value": row.get(q_col, "") if q_col else "",
            "unique_raw": row.get(unique_col, "") if unique_col else "",
            "shared_raw": row.get(shared_col, "") if shared_col else "",
            "unique_peptides": unique_peps,
            "shared_peptides": shared_peps,
            "all_peptides": combined,
        })

        if len(extraction_audit) < 20:
            extraction_audit.append({
                "row_number": row_no,
                "protein_accession": row.get(accession_col, ""),
                "unique_raw": str(row.get(unique_col, "") if unique_col else "")[:1000],
                "shared_raw": str(row.get(shared_col, "") if shared_col else "")[:1000],
                "unique_extracted": " ; ".join(unique_peps[:20]),
                "shared_extracted": " ; ".join(shared_peps[:20]),
                "n_extracted": len(combined),
            })

    process(first, 1)

    for i, row in enumerate(rows, start=2):
        process(row, i)

    meta = {
        "fields": fields,
        "accession_col": accession_col,
        "unique_col": unique_col or "",
        "shared_col": shared_col or "",
        "q_col": q_col or "",
        "dct_col": dct_col or "",
        "counters": dict(counters),
        "extraction_audit": extraction_audit,
    }

    return groups, all_peptides, meta

# -----------------------------------------------------------------------------
# Pure-Python Aho-Corasick for peptide -> FASTA mapping
# -----------------------------------------------------------------------------

class ACAutomaton:
    def __init__(self):
        self.next = [{}]
        self.fail = [0]
        self.out = [[]]

    def add(self, pattern):
        state = 0

        for ch in pattern:
            nxt = self.next[state].get(ch)

            if nxt is None:
                nxt = len(self.next)
                self.next[state][ch] = nxt
                self.next.append({})
                self.fail.append(0)
                self.out.append([])

            state = nxt

        self.out[state].append(pattern)

    def build(self):
        q = deque()

        for _ch, state in self.next[0].items():
            self.fail[state] = 0
            q.append(state)

        while q:
            r = q.popleft()

            for ch, s in self.next[r].items():
                q.append(s)

                f = self.fail[r]

                while f and ch not in self.next[f]:
                    f = self.fail[f]

                self.fail[s] = self.next[f].get(ch, 0)
                self.out[s].extend(self.out[self.fail[s]])

    def find_unique(self, text):
        state = 0
        found = set()

        for ch in text:
            while state and ch not in self.next[state]:
                state = self.fail[state]

            state = self.next[state].get(ch, 0)

            if self.out[state]:
                found.update(self.out[state])

        return found

def map_peptides_to_fasta(peptides, fasta):
    """
    Map I/L-equivalent peptide keys against target and entrapment proteins.
    """
    normalized_to_original = {}

    for pep in peptides:
        key = normalize_il(pep)
        normalized_to_original.setdefault(key, set()).add(pep)

    patterns = sorted(normalized_to_original)

    ac = ACAutomaton()

    for p in patterns:
        ac.add(p)

    ac.build()

    hit = {
        p: {
            "target": False,
            "entrapment": False,
            "first_target_header": "",
            "first_entrapment_header": "",
        }
        for p in patterns
    }

    n_entries = 0

    for header, seq in parse_fasta(fasta):
        n_entries += 1

        seq_il = normalize_il(seq)
        found = ac.find_unique(seq_il)

        if found:
            entrap = is_entrapment_header(header)

            for p in found:
                rec = hit[p]

                if entrap:
                    rec["entrapment"] = True
                    if not rec["first_entrapment_header"]:
                        rec["first_entrapment_header"] = header
                else:
                    rec["target"] = True
                    if not rec["first_target_header"]:
                        rec["first_target_header"] = header

        if n_entries % 25000 == 0:
            print(f"  FASTA entries scanned: {n_entries:,}")

    mapping = {}
    rows = []

    for key, originals in normalized_to_original.items():
        rec = hit[key]

        if rec["target"] and rec["entrapment"]:
            cls = "SHARED_TARGET_ENTRAPMENT"
        elif rec["target"]:
            cls = "TARGET_ONLY"
        elif rec["entrapment"]:
            cls = "ENTRAPMENT_ONLY"
        else:
            cls = "ABSENT_FROM_FROZEN_FASTA"

        for pep in originals:
            mapping[pep] = cls

            rows.append({
                "peptide": pep,
                "peptide_IL_equivalent": key,
                "mapping_class": cls,
                "found_target": rec["target"],
                "found_entrapment": rec["entrapment"],
                "first_target_header": rec["first_target_header"],
                "first_entrapment_header": rec["first_entrapment_header"],
            })

    return mapping, rows

def evaluate_mm_reconstruction(groups, peptide_mapping, r_protein, nominal):
    counts = Counter()
    details = []

    for g in groups:
        pep_classes = Counter(
            peptide_mapping.get(p, "ABSENT_FROM_FROZEN_FASTA")
            for p in g["all_peptides"]
        )

        n_t = pep_classes["TARGET_ONLY"]
        n_e = pep_classes["ENTRAPMENT_ONLY"]
        n_s = pep_classes["SHARED_TARGET_ENTRAPMENT"]
        n_a = pep_classes["ABSENT_FROM_FROZEN_FASTA"]

        if n_t > 0 and n_e == 0:
            cls = "TARGET_ONLY_EVIDENCE"
        elif n_e > 0 and n_t == 0:
            cls = "ENTRAPMENT_ONLY_EVIDENCE"
        elif n_t > 0 and n_e > 0:
            cls = "MIXED_TARGET_ENTRAPMENT_EVIDENCE"
        else:
            cls = "AMBIGUOUS_NO_DISCRIMINATIVE_PEPTIDE"

        counts[cls] += 1

        details.append({
            "workflow": "MM",
            "row_number": g["row_number"],
            "protein_accession": g["protein_accession"],
            "q_value": g["q_value"],
            "reconstructed_group_class": cls,
            "n_all_extracted_peptides": len(g["all_peptides"]),
            "n_target_only_peptides": n_t,
            "n_entrapment_only_peptides": n_e,
            "n_shared_target_entrapment_peptides": n_s,
            "n_absent_peptides": n_a,
            "target_only_peptides": " ; ".join(
                p for p in g["all_peptides"]
                if peptide_mapping.get(p) == "TARGET_ONLY"
            ),
            "entrapment_only_peptides": " ; ".join(
                p for p in g["all_peptides"]
                if peptide_mapping.get(p) == "ENTRAPMENT_ONLY"
            ),
            "shared_peptides": " ; ".join(
                p for p in g["all_peptides"]
                if peptide_mapping.get(p) == "SHARED_TARGET_ENTRAPMENT"
            ),
            "absent_peptides": " ; ".join(
                p for p in g["all_peptides"]
                if peptide_mapping.get(p) == "ABSENT_FROM_FROZEN_FASTA"
            ),
        })

    T = counts["TARGET_ONLY_EVIDENCE"]
    E = counts["ENTRAPMENT_ONLY_EVIDENCE"]
    stats = fdp_stats(T, E, r_protein, nominal)

    total = len(groups)
    ambiguous = counts["AMBIGUOUS_NO_DISCRIMINATIVE_PEPTIDE"]
    mixed = counts["MIXED_TARGET_ENTRAPMENT_EVIDENCE"]

    summary = {
        "workflow": "MM",
        "analysis_type": "SEQUENCE_RECONSTRUCTION_DIAGNOSTIC_ONLY",
        "accepted_native_protein_groups": total,
        "target_only_evidence_groups": T,
        "entrapment_only_evidence_groups": E,
        "mixed_evidence_groups": mixed,
        "ambiguous_groups": ambiguous,
        "reconstructed_evaluable_T_plus_E": T + E,
        "pct_groups_evaluable_by_sequence_reconstruction":
            (100 * (T + E) / total) if total else math.nan,
        "pct_groups_ambiguous":
            (100 * ambiguous / total) if total else math.nan,
        "r_protein": r_protein,
        "nominal_protein_fdr": nominal,
        "diagnostic_lower_bound_fdp": stats["lower"],
        "diagnostic_combined_fdp": stats["combined"],
        "diagnostic_lower_bound_pct": stats["lower_pct"],
        "diagnostic_combined_fdp_pct": stats["combined_pct"],
        "diagnostic_interpretation": stats["interpretation"],
        "scientific_status":
            "DIAGNOSTIC_ONLY_NOT_DIRECT_NATIVE_PROTEIN_GROUP_FDR",
    }

    return summary, details

# =============================================================================
# Import AP/FP prior direct results
# =============================================================================

def locate_step02d_summary(root):
    preferred = (
        root
        / "Step02D_PROTEIN_ENTRAPMENT_FDR_v103"
        / "Step02D_protein_entrapment_summary.tsv"
    )

    if preferred.exists():
        return preferred

    hits = list(root.rglob("Step02D_protein_entrapment_summary.tsv"))

    if not hits:
        return None

    return sorted(hits, key=lambda p: len(str(p)))[-1]

def import_ap_fp(step02d_path):
    out = {}

    if step02d_path is None:
        return out

    for row in read_delimited(step02d_path):
        wf = str(row.get("workflow", "")).strip()

        if wf not in {"AP", "FP"}:
            continue

        out[wf] = {
            "workflow": wf,
            "analysis_type": "DIRECT_PROTEIN_GROUP_FROM_STEP02D",
            "target_only_groups": row.get("target_only_groups", ""),
            "entrapment_only_groups": row.get("entrapment_only_groups", ""),
            "mixed_groups": row.get("mixed_target_entrapment_groups", ""),
            "unresolved_groups": row.get("unresolved_groups", ""),
            "lower_bound_pct": row.get("primary_lower_bound_pct", ""),
            "combined_fdp_pct": row.get("primary_combined_fdp_pct", ""),
            "interpretation": row.get("primary_interpretation", ""),
            "source": str(step02d_path),
        }

    return out

# =============================================================================
# Main
# =============================================================================

def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--root",
        default=str(DEFAULT_ROOT),
        help="ENTRAPMENT_FDR root",
    )

    ap.add_argument(
        "--alpha",
        type=float,
        default=0.01,
        help="Nominal protein FDR cutoff (default 0.01)",
    )

    ap.add_argument(
        "--outdir",
        default=None,
        help="Optional output directory",
    )

    args = ap.parse_args()

    root = Path(args.root)
    nominal = float(args.alpha)

    if not root.exists():
        raise SystemExit(f"Root does not exist:\n{root}")

    outdir = (
        Path(args.outdir)
        if args.outdir
        else root / "Step02F_MQ_FINAL_MM_SEQUENCE_SALVAGE_v100"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 100)
    print(f"STEP02F MQ FINAL + MM SEQUENCE SALVAGE v{VERSION}")
    print("=" * 100)
    print("Root:", root)
    print("Nominal protein FDR:", f"{100*nominal:.2f}%")

    # -------------------------------------------------------------------------
    # Frozen database QC
    # -------------------------------------------------------------------------

    step02a = locate_step02a(root)
    fasta = step02a / "Step02_target_plus_shuffled_entrapment_r1.fasta"

    if not fasta.exists():
        raise FileNotFoundError(f"Frozen FASTA not found:\n{fasta}")

    print("\n[1/6] Frozen FASTA QC")

    target_n, entrap_n = count_fasta_classes(fasta)

    print("Target proteins:", f"{target_n:,}")
    print("Entrapment proteins:", f"{entrap_n:,}")

    if (
        target_n != EXPECTED_TARGET_PROTEINS
        or entrap_n != EXPECTED_ENTRAPMENT_PROTEINS
    ):
        raise SystemExit(
            "Frozen FASTA protein counts do not match Step02A QC."
        )

    r_protein = entrap_n / target_n
    print("Protein-space r =", f"{r_protein:.12f}")

    # -------------------------------------------------------------------------
    # MQ direct
    # -------------------------------------------------------------------------

    print("\n[2/6] MaxQuant direct protein-group entrapment evaluation")

    mq_dir = locate_workflow_dir(root, "MQ")
    mq_path = locate_mq_protein_groups(mq_dir)

    print("MQ source:", mq_path)

    mq_summary, mq_details = evaluate_mq(
        mq_path,
        r_protein,
        nominal,
    )

    print(
        "MQ PURE groups: "
        f"T={mq_summary['target_only_groups']:,} | "
        f"E={mq_summary['entrapment_only_groups']:,} | "
        f"mixed={mq_summary['mixed_target_entrapment_groups']:,} | "
        f"unresolved={mq_summary['unresolved_groups']:,}"
    )

    print(
        "MQ protein FDP: "
        f"lower={mq_summary['lower_bound_pct']:.4f}% | "
        f"combined={mq_summary['combined_fdp_pct']:.4f}% | "
        f"{mq_summary['interpretation']}"
    )

    # -------------------------------------------------------------------------
    # MM load
    # -------------------------------------------------------------------------

    print("\n[3/6] MetaMorpheus native protein table + peptide extraction")

    mm_dir = locate_workflow_dir(root, "MM")
    mm_path = locate_mm_protein_table(mm_dir)

    print("MM source:", mm_path)

    mm_groups, mm_peptides, mm_meta = load_mm_groups_and_peptides(
        mm_path,
        nominal,
    )

    print("MM accepted groups:", f"{len(mm_groups):,}")
    print("MM distinct extracted peptide sequences:", f"{len(mm_peptides):,}")
    print(
        "MM groups with no extracted peptide:",
        f"{mm_meta['counters'].get('accepted_rows_no_extracted_peptide', 0):,}"
    )

    if len(mm_peptides) == 0:
        raise RuntimeError(
            "No MM peptide sequences could be extracted. "
            "See method/audit output before proceeding."
        )

    # -------------------------------------------------------------------------
    # MM peptide -> FASTA mapping
    # -------------------------------------------------------------------------

    print("\n[4/6] Mapping MM peptide evidence to frozen target/entrapment FASTA")
    print("This scans the frozen FASTA once.")

    mm_peptide_mapping, mm_peptide_rows = map_peptides_to_fasta(
        mm_peptides,
        fasta,
    )

    map_counts = Counter(
        r["mapping_class"]
        for r in mm_peptide_rows
    )

    print(
        "Peptide mapping: "
        f"target-only={map_counts['TARGET_ONLY']:,} | "
        f"entrapment-only={map_counts['ENTRAPMENT_ONLY']:,} | "
        f"shared={map_counts['SHARED_TARGET_ENTRAPMENT']:,} | "
        f"absent={map_counts['ABSENT_FROM_FROZEN_FASTA']:,}"
    )

    # -------------------------------------------------------------------------
    # MM reconstruct group evidence
    # -------------------------------------------------------------------------

    print("\n[5/6] Reconstructing MM group identity from peptide sequence evidence")

    mm_summary, mm_details = evaluate_mm_reconstruction(
        mm_groups,
        mm_peptide_mapping,
        r_protein,
        nominal,
    )

    print(
        "MM reconstructed groups: "
        f"T-evidence={mm_summary['target_only_evidence_groups']:,} | "
        f"E-evidence={mm_summary['entrapment_only_evidence_groups']:,} | "
        f"mixed={mm_summary['mixed_evidence_groups']:,} | "
        f"ambiguous={mm_summary['ambiguous_groups']:,}"
    )

    print(
        "MM evaluable by sequence reconstruction: "
        f"{mm_summary['pct_groups_evaluable_by_sequence_reconstruction']:.2f}%"
    )

    print(
        "MM diagnostic FDP only: "
        f"lower={mm_summary['diagnostic_lower_bound_pct']:.4f}% | "
        f"combined={mm_summary['diagnostic_combined_fdp_pct']:.4f}% | "
        f"{mm_summary['diagnostic_interpretation']}"
    )

    print(
        "MM scientific status:",
        mm_summary["scientific_status"],
    )

    # -------------------------------------------------------------------------
    # Save all
    # -------------------------------------------------------------------------

    print("\n[6/6] Writing outputs")

    # MQ
    write_tsv(
        outdir / "Step02F_MQ_direct_protein_FDR.tsv",
        [mq_summary],
        list(mq_summary.keys()),
    )

    write_tsv(
        outdir / "Step02F_MQ_group_classification.tsv",
        mq_details,
        [
            "workflow",
            "row_number",
            "protein_ids",
            "q_value",
            "group_class",
            "n_fasta_headers",
            "n_target_headers",
            "n_entrapment_headers",
            "fasta_headers",
        ],
    )

    # MM
    write_tsv(
        outdir / "Step02F_MM_sequence_reconstruction_summary.tsv",
        [mm_summary],
        list(mm_summary.keys()),
    )

    write_tsv(
        outdir / "Step02F_MM_group_sequence_classification.tsv",
        mm_details,
        [
            "workflow",
            "row_number",
            "protein_accession",
            "q_value",
            "reconstructed_group_class",
            "n_all_extracted_peptides",
            "n_target_only_peptides",
            "n_entrapment_only_peptides",
            "n_shared_target_entrapment_peptides",
            "n_absent_peptides",
            "target_only_peptides",
            "entrapment_only_peptides",
            "shared_peptides",
            "absent_peptides",
        ],
    )

    write_tsv(
        outdir / "Step02F_MM_peptide_mapping.tsv",
        mm_peptide_rows,
        [
            "peptide",
            "peptide_IL_equivalent",
            "mapping_class",
            "found_target",
            "found_entrapment",
            "first_target_header",
            "first_entrapment_header",
        ],
    )

    write_tsv(
        outdir / "Step02F_MM_extraction_audit.tsv",
        mm_meta["extraction_audit"],
        [
            "row_number",
            "protein_accession",
            "unique_raw",
            "shared_raw",
            "unique_extracted",
            "shared_extracted",
            "n_extracted",
        ],
    )

    # Combined status
    prior_path = locate_step02d_summary(root)
    prior = import_ap_fp(prior_path)

    combined = []

    for wf in ["AP", "FP"]:
        if wf in prior:
            combined.append(prior[wf])

    combined.append({
        "workflow": "MQ",
        "analysis_type": "DIRECT_PROTEIN_GROUP_FROM_FASTA_HEADERS",
        "target_only_groups": mq_summary["target_only_groups"],
        "entrapment_only_groups": mq_summary["entrapment_only_groups"],
        "mixed_groups": mq_summary["mixed_target_entrapment_groups"],
        "unresolved_groups": mq_summary["unresolved_groups"],
        "lower_bound_pct": mq_summary["lower_bound_pct"],
        "combined_fdp_pct": mq_summary["combined_fdp_pct"],
        "interpretation": mq_summary["interpretation"],
        "source": str(mq_path),
    })

    combined.append({
        "workflow": "MM",
        "analysis_type": "SEQUENCE_RECONSTRUCTION_DIAGNOSTIC_ONLY",
        "target_only_groups": mm_summary["target_only_evidence_groups"],
        "entrapment_only_groups": mm_summary["entrapment_only_evidence_groups"],
        "mixed_groups": mm_summary["mixed_evidence_groups"],
        "unresolved_groups": mm_summary["ambiguous_groups"],
        "lower_bound_pct": mm_summary["diagnostic_lower_bound_pct"],
        "combined_fdp_pct": mm_summary["diagnostic_combined_fdp_pct"],
        "interpretation":
            mm_summary["diagnostic_interpretation"]
            + " [DIAGNOSTIC_ONLY]",
        "source": str(mm_path),
    })

    write_tsv(
        outdir / "Step02F_combined_status.tsv",
        combined,
        [
            "workflow",
            "analysis_type",
            "target_only_groups",
            "entrapment_only_groups",
            "mixed_groups",
            "unresolved_groups",
            "lower_bound_pct",
            "combined_fdp_pct",
            "interpretation",
            "source",
        ],
    )

    # Method record
    method_record = {
        "script_version": VERSION,
        "nominal_protein_fdr": nominal,
        "frozen_fasta": str(fasta),
        "target_proteins": target_n,
        "entrapment_proteins": entrap_n,
        "r_protein": r_protein,
        "mq_method":
            "Direct classification of each accepted MaxQuant protein group "
            "using the native 'Fasta headers' column. A group is TARGET_ONLY "
            "if all retained headers lack _p_target, ENTRAPMENT_ONLY if all "
            "retained headers contain _p_target, MIXED if both occur.",
        "mm_method":
            "Diagnostic sequence reconstruction only. Unique Peptides and "
            "Shared Peptides from the accepted MetaMorpheus protein table are "
            "mapped I/L-equivalently to the frozen target+entrapment FASTA. "
            "Group identity is reconstructed from target-only and entrapment-"
            "only peptide evidence.",
        "mm_limitation":
            "MetaMorpheus native protein identifiers lost the _p_target label; "
            "therefore the reconstructed MM FDP is not a direct native "
            "protein-group entrapment FDR and must not be reported as such "
            "without an explicit methodological qualification.",
        "lower_fdp_formula": "E / (T + E)",
        "combined_fdp_formula": "E * (1 + 1/r_protein) / (T + E)",
    }

    (
        outdir / "Step02F_method_record.json"
    ).write_text(
        json.dumps(method_record, indent=2),
        encoding="utf-8",
    )

    # Human-readable report
    report_lines = [
        f"Step02F MQ FINAL + MM SEQUENCE SALVAGE v{VERSION}",
        f"Nominal protein FDR: {100*nominal:.2f}%",
        f"Protein-space r: {r_protein:.12f}",
        "",
        "MAXQUANT — DIRECT PROTEIN-GROUP RESULT",
        f"  target-only groups: {mq_summary['target_only_groups']}",
        f"  entrapment-only groups: {mq_summary['entrapment_only_groups']}",
        f"  mixed groups: {mq_summary['mixed_target_entrapment_groups']}",
        f"  unresolved groups: {mq_summary['unresolved_groups']}",
        f"  lower FDP: {mq_summary['lower_bound_pct']:.6f}%",
        f"  combined FDP: {mq_summary['combined_fdp_pct']:.6f}%",
        f"  interpretation: {mq_summary['interpretation']}",
        "",
        "METAMORPHEUS — SEQUENCE-RECONSTRUCTION DIAGNOSTIC",
        f"  accepted native protein groups: {mm_summary['accepted_native_protein_groups']}",
        f"  target-only evidence groups: {mm_summary['target_only_evidence_groups']}",
        f"  entrapment-only evidence groups: {mm_summary['entrapment_only_evidence_groups']}",
        f"  mixed evidence groups: {mm_summary['mixed_evidence_groups']}",
        f"  ambiguous groups: {mm_summary['ambiguous_groups']}",
        (
            "  evaluable by reconstruction: "
            f"{mm_summary['pct_groups_evaluable_by_sequence_reconstruction']:.4f}%"
        ),
        f"  diagnostic lower FDP: {mm_summary['diagnostic_lower_bound_pct']:.6f}%",
        f"  diagnostic combined FDP: {mm_summary['diagnostic_combined_fdp_pct']:.6f}%",
        f"  diagnostic interpretation: {mm_summary['diagnostic_interpretation']}",
        "",
        "IMPORTANT:",
        "  MQ is direct/native protein-group entrapment evaluation.",
        "  MM is diagnostic reconstruction only because native MM output lost",
        "  the _p_target identity label.",
    ]

    report_path = outdir / "Step02F_REPORT.txt"

    report_path.write_text(
        "\n".join(report_lines),
        encoding="utf-8",
    )

    print("\n" + "=" * 100)
    print("STEP02F COMPLETE")
    print("=" * 100)
    print("MQ:", mq_summary["interpretation"])
    print(
        "MM:",
        mm_summary["diagnostic_interpretation"],
        "[DIAGNOSTIC ONLY]",
    )
    print("\nSend back these files:")
    print(report_path)
    print(outdir / "Step02F_combined_status.tsv")
    print(outdir / "Step02F_MQ_direct_protein_FDR.tsv")
    print(outdir / "Step02F_MM_sequence_reconstruction_summary.tsv")
    print(outdir / "Step02F_MM_extraction_audit.tsv")
    print("=" * 100)


if __name__ == "__main__":
    main()
