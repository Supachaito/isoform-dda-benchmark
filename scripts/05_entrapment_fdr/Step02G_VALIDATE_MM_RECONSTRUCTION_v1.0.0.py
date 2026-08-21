#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step02G_VALIDATE_MM_RECONSTRUCTION_v1.0.0.py

Purpose
-------
Final validation of the MetaMorpheus (MM) sequence-based entrapment reconstruction
produced by Step02F, before deciding whether an MM rerun is technically necessary.

This script does NOT rerun MetaMorpheus and does NOT modify any prior result.

Validation layers
-----------------
A. Internal group-consistency QC
   Checks every Step02F MM group:
     TARGET_ONLY_EVIDENCE      => target-only peptides >0, entrapment-only =0
     ENTRAPMENT_ONLY_EVIDENCE  => entrapment-only peptides >0, target-only =0
     MIXED                     => both >0
     AMBIGUOUS                 => neither >0

B. Independent FASTA remapping QC
   Re-maps ALL discriminatory peptides used by Step02F directly against the frozen
   target+entrapment FASTA, independently of the Step02F peptide-mapping table.

C. Native MM peptide-field provenance
   Re-reads AllQuantifiedProteinGroups.tsv and determines whether each
   discriminatory peptide was present in the native "Unique Peptides" field,
   "Shared Peptides" field, or both.

D. Entrapment-group support audit
   Audits all reconstructed entrapment groups individually:
     number of confirmed entrapment-only peptides
     number from native Unique Peptides
     number from native Shared Peptides
     single-peptide vs multi-peptide support

E. Automatic QC verdict
   PASS_RECONSTRUCTION_VALIDATED if:
     - zero group-consistency contradictions
     - zero independent-remapping mismatches
     - every reconstructed entrapment group has >=1 independently confirmed
       entrapment-only peptide
     - every reconstructed target group has >=1 independently confirmed
       target-only peptide
     - no discriminatory peptide is absent from the frozen FASTA

Scientific note
---------------
PASS means the post-hoc sequence reconstruction is technically internally
validated against the exact frozen FASTA. It does NOT change the fact that the
native MetaMorpheus protein table itself lost the "_p_target" label. The report
keeps this distinction explicit.

Default root
------------
C:\\Users\\Supachai\\Desktop\\AphaPept_benchmark\\Benchmark_Program\\ENTRAPMENT_FDR
"""

from __future__ import annotations

import argparse
import csv
import json
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

EXPECTED_GROUPS = 4112
EXPECTED_TARGET_GROUPS = 4086
EXPECTED_ENTRAPMENT_GROUPS = 26

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
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)

def int0(x):
    try:
        return int(float(str(x).strip()))
    except Exception:
        return 0

def split_semicolon_peptides(text):
    return [
        x.strip()
        for x in str(text or "").split(";")
        if x.strip()
    ]

def normalize_il(seq):
    return str(seq or "").upper().replace("I", "J").replace("L", "J")

def is_entrapment_header(header):
    return MARKER in str(header or "").lower()

# -----------------------------------------------------------------------------
# Locate files
# -----------------------------------------------------------------------------

def locate_step02f(root):
    direct = root / "Step02F_MQ_FINAL_MM_SEQUENCE_SALVAGE_v100"
    if direct.exists():
        return direct

    hits = [
        p for p in root.rglob("*")
        if p.is_dir()
        and p.name.lower().startswith("step02f_mq_final_mm_sequence_salvage")
    ]
    if not hits:
        raise FileNotFoundError(
            f"Could not find Step02F output folder under:\n{root}"
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

def locate_mm_dir(root):
    direct = root / "MM_ENTRAP"
    if direct.exists():
        return direct

    hits = [
        p for p in root.rglob("*")
        if p.is_dir()
        and p.name.lower() == "mm_entrap"
    ]
    if not hits:
        raise FileNotFoundError(
            f"Could not find MM_ENTRAP under:\n{root}"
        )
    return sorted(hits, key=lambda p: len(str(p)))[0]

def locate_mm_native_table(mm_dir):
    preferred = sorted(mm_dir.rglob("AllQuantifiedProteinGroups.tsv"))
    if preferred:
        return preferred[0]

    candidates = []
    for p in mm_dir.rglob("*"):
        if not p.is_file() or p.suffix.lower() not in {".tsv", ".txt", ".csv"}:
            continue
        if "protein" not in p.name.lower() or "group" not in p.name.lower():
            continue
        try:
            rows = read_delimited(p)
            first = next(rows, None)
        except Exception:
            continue
        if first is None:
            continue
        fields = list(first.keys())
        if (
            find_col(fields, ["Protein Accession", "Protein Accessions"])
            and find_col(fields, ["Unique Peptides"])
        ):
            candidates.append(p)

    if not candidates:
        raise FileNotFoundError(
            f"No MM protein-group table found under:\n{mm_dir}"
        )

    return sorted(
        candidates,
        key=lambda p: (-p.stat().st_size, len(str(p)))
    )[0]

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

# -----------------------------------------------------------------------------
# Aho-Corasick matcher
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

# -----------------------------------------------------------------------------
# Load Step02F group classifications
# -----------------------------------------------------------------------------

def load_step02f_groups(path):
    rows = []
    for row in read_delimited(path):
        rows.append(row)
    return rows

def validate_internal_group_logic(groups):
    issues = []
    class_counts = Counter()

    for row in groups:
        cls = row.get("reconstructed_group_class", "")
        class_counts[cls] += 1

        nt = int0(row.get("n_target_only_peptides"))
        ne = int0(row.get("n_entrapment_only_peptides"))
        ns = int0(row.get("n_shared_target_entrapment_peptides"))
        na = int0(row.get("n_absent_peptides"))
        nall = int0(row.get("n_all_extracted_peptides"))

        reasons = []

        if nall != nt + ne + ns + na:
            reasons.append(
                f"n_all({nall}) != T({nt})+E({ne})+shared({ns})+absent({na})"
            )

        if cls == "TARGET_ONLY_EVIDENCE":
            if not (nt > 0 and ne == 0):
                reasons.append("TARGET_ONLY_EVIDENCE logic violated")

        elif cls == "ENTRAPMENT_ONLY_EVIDENCE":
            if not (ne > 0 and nt == 0):
                reasons.append("ENTRAPMENT_ONLY_EVIDENCE logic violated")

        elif cls == "MIXED_TARGET_ENTRAPMENT_EVIDENCE":
            if not (nt > 0 and ne > 0):
                reasons.append("MIXED evidence logic violated")

        elif cls == "AMBIGUOUS_NO_DISCRIMINATIVE_PEPTIDE":
            if not (nt == 0 and ne == 0):
                reasons.append("AMBIGUOUS logic violated")

        else:
            reasons.append(f"unknown class: {cls}")

        if reasons:
            issues.append({
                "row_number": row.get("row_number", ""),
                "protein_accession": row.get("protein_accession", ""),
                "class": cls,
                "reason": " | ".join(reasons),
            })

    return issues, class_counts

# -----------------------------------------------------------------------------
# Independent remapping
# -----------------------------------------------------------------------------

def gather_discriminatory_peptides(groups):
    expected = {}

    for row in groups:
        for pep in split_semicolon_peptides(row.get("target_only_peptides")):
            prev = expected.get(pep)
            if prev and prev != "TARGET_ONLY":
                raise RuntimeError(
                    f"Peptide has conflicting Step02F classes: {pep}"
                )
            expected[pep] = "TARGET_ONLY"

        for pep in split_semicolon_peptides(row.get("entrapment_only_peptides")):
            prev = expected.get(pep)
            if prev and prev != "ENTRAPMENT_ONLY":
                raise RuntimeError(
                    f"Peptide has conflicting Step02F classes: {pep}"
                )
            expected[pep] = "ENTRAPMENT_ONLY"

    return expected

def independent_remap(discriminatory, fasta):
    norm_to_peps = {}

    for pep in discriminatory:
        key = normalize_il(pep)
        norm_to_peps.setdefault(key, set()).add(pep)

    patterns = sorted(norm_to_peps)

    ac = ACAutomaton()
    for p in patterns:
        ac.add(p)
    ac.build()

    hits = {
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

        found = ac.find_unique(normalize_il(seq))
        if not found:
            continue

        ent = is_entrapment_header(header)

        for key in found:
            rec = hits[key]

            if ent:
                rec["entrapment"] = True
                if not rec["first_entrapment_header"]:
                    rec["first_entrapment_header"] = header
            else:
                rec["target"] = True
                if not rec["first_target_header"]:
                    rec["first_target_header"] = header

        if n_entries % 25000 == 0:
            print(f"  FASTA entries scanned: {n_entries:,}")

    results = []
    mismatches = []

    for key, peps in norm_to_peps.items():
        rec = hits[key]

        if rec["target"] and rec["entrapment"]:
            observed = "SHARED_TARGET_ENTRAPMENT"
        elif rec["target"]:
            observed = "TARGET_ONLY"
        elif rec["entrapment"]:
            observed = "ENTRAPMENT_ONLY"
        else:
            observed = "ABSENT_FROM_FROZEN_FASTA"

        for pep in sorted(peps):
            expected = discriminatory[pep]
            match = expected == observed

            row = {
                "peptide": pep,
                "peptide_IL_equivalent": key,
                "expected_step02f_class": expected,
                "independent_observed_class": observed,
                "match": match,
                "found_target": rec["target"],
                "found_entrapment": rec["entrapment"],
                "first_target_header": rec["first_target_header"],
                "first_entrapment_header": rec["first_entrapment_header"],
            }

            results.append(row)

            if not match:
                mismatches.append(row)

    return results, mismatches

# -----------------------------------------------------------------------------
# Re-read native MM peptide fields
# -----------------------------------------------------------------------------

AA = set("ACDEFGHIKLMNPQRSTVWYBJOUXZ")

def clean_peptide_token(token):
    s = str(token or "").strip()

    if not s:
        return ""

    m = re.fullmatch(r"[A-Za-z\-]\.([^\.]+)\.[A-Za-z\-]", s)
    if m:
        s = m.group(1)

    s = re.sub(r"\[[^\]]*\]", "", s)
    s = re.sub(r"\([^\)]*\)", "", s)
    s = re.sub(r"\{[^\}]*\}", "", s)
    s = re.sub(r"[^A-Za-z]", "", s).upper()

    if len(s) < 6:
        return ""

    if any(ch not in AA for ch in s):
        return ""

    return s

def extract_peptides_from_field(text):
    s = str(text or "").strip()

    if not s:
        return []

    coarse = re.split(r"[;,|]+", s)
    peptides = []

    for chunk in coarse:
        chunk = chunk.strip()
        if not chunk:
            continue

        whole = clean_peptide_token(chunk)
        if whole:
            peptides.append(whole)
            continue

        for tok in chunk.split():
            pep = clean_peptide_token(tok)
            if pep:
                peptides.append(pep)

    return list(dict.fromkeys(peptides))

def load_native_mm_peptide_provenance(mm_path):
    rows = read_delimited(mm_path)
    first = next(rows, None)

    if first is None:
        raise RuntimeError(f"Empty native MM protein table: {mm_path}")

    fields = list(first.keys())

    accession_col = find_col(
        fields,
        ["Protein Accession", "Protein Accessions"]
    )
    unique_col = find_col(fields, ["Unique Peptides"])
    shared_col = find_col(fields, ["Shared Peptides"])

    if accession_col is None or unique_col is None:
        raise RuntimeError(
            "Required MM columns not found.\n"
            f"Columns: {fields}"
        )

    provenance = {}

    def process(row, row_no):
        unique_peps = set(extract_peptides_from_field(row.get(unique_col)))
        shared_peps = set(
            extract_peptides_from_field(row.get(shared_col))
            if shared_col else []
        )

        provenance[str(row_no)] = {
            "protein_accession": row.get(accession_col, ""),
            "unique": unique_peps,
            "shared": shared_peps,
        }

    process(first, 1)

    for i, row in enumerate(rows, start=2):
        process(row, i)

    return provenance

# -----------------------------------------------------------------------------
# Entrapment/target group support audit
# -----------------------------------------------------------------------------

def audit_group_support(groups, remap_lookup, native_provenance):
    rows = []
    failure_count = 0

    for row in groups:
        cls = row.get("reconstructed_group_class", "")
        row_no = str(row.get("row_number", ""))

        if cls not in {
            "TARGET_ONLY_EVIDENCE",
            "ENTRAPMENT_ONLY_EVIDENCE",
        }:
            continue

        if cls == "TARGET_ONLY_EVIDENCE":
            discriminatory_peps = split_semicolon_peptides(
                row.get("target_only_peptides")
            )
            required_class = "TARGET_ONLY"
        else:
            discriminatory_peps = split_semicolon_peptides(
                row.get("entrapment_only_peptides")
            )
            required_class = "ENTRAPMENT_ONLY"

        native = native_provenance.get(
            row_no,
            {"unique": set(), "shared": set(), "protein_accession": ""}
        )

        n_confirmed = 0
        n_native_unique = 0
        n_native_shared = 0
        n_native_both = 0
        n_native_neither = 0
        bad_peps = []

        for pep in discriminatory_peps:
            observed = remap_lookup.get(
                pep,
                "MISSING_FROM_INDEPENDENT_REMAP"
            )

            if observed == required_class:
                n_confirmed += 1
            else:
                bad_peps.append(f"{pep}:{observed}")

            in_u = pep in native["unique"]
            in_s = pep in native["shared"]

            if in_u and in_s:
                n_native_both += 1
            elif in_u:
                n_native_unique += 1
            elif in_s:
                n_native_shared += 1
            else:
                n_native_neither += 1

        technical_pass = (
            len(discriminatory_peps) > 0
            and n_confirmed == len(discriminatory_peps)
            and len(bad_peps) == 0
            and n_native_neither == 0
        )

        if not technical_pass:
            failure_count += 1

        support_tier = (
            "MULTI_PEPTIDE_SUPPORT"
            if len(discriminatory_peps) >= 2
            else "SINGLE_PEPTIDE_SUPPORT"
        )

        rows.append({
            "row_number": row_no,
            "protein_accession_step02f": row.get("protein_accession", ""),
            "protein_accession_native": native.get("protein_accession", ""),
            "reconstructed_group_class": cls,
            "required_discriminatory_class": required_class,
            "n_discriminatory_peptides": len(discriminatory_peps),
            "n_independently_confirmed": n_confirmed,
            "n_from_native_unique_field": n_native_unique,
            "n_from_native_shared_field": n_native_shared,
            "n_present_in_both_native_fields": n_native_both,
            "n_missing_from_native_peptide_fields": n_native_neither,
            "support_tier": support_tier,
            "technical_pass": technical_pass,
            "bad_peptides": " ; ".join(bad_peps),
            "discriminatory_peptides": " ; ".join(discriminatory_peps),
        })

    return rows, failure_count

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--root",
        default=str(DEFAULT_ROOT),
        help="ENTRAPMENT_FDR root",
    )

    ap.add_argument(
        "--outdir",
        default=None,
        help="Optional output directory",
    )

    args = ap.parse_args()

    root = Path(args.root)
    if not root.exists():
        raise SystemExit(f"Root does not exist:\n{root}")

    outdir = (
        Path(args.outdir)
        if args.outdir
        else root / "Step02G_VALIDATE_MM_RECONSTRUCTION_v100"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    step02f = locate_step02f(root)
    step02a = locate_step02a(root)
    mm_dir = locate_mm_dir(root)

    group_file = (
        step02f
        / "Step02F_MM_group_sequence_classification.tsv"
    )
    peptide_file = (
        step02f
        / "Step02F_MM_peptide_mapping.tsv"
    )
    fasta = (
        step02a
        / "Step02_target_plus_shuffled_entrapment_r1.fasta"
    )
    mm_native = locate_mm_native_table(mm_dir)

    for required in [group_file, peptide_file, fasta, mm_native]:
        if not required.exists():
            raise FileNotFoundError(required)

    print("=" * 100)
    print(f"STEP02G VALIDATE MM RECONSTRUCTION v{VERSION}")
    print("=" * 100)
    print("Root:", root)
    print("Step02F groups:", group_file)
    print("Step02F peptide map:", peptide_file)
    print("Frozen FASTA:", fasta)
    print("Native MM table:", mm_native)

    # -------------------------------------------------------------------------
    # 1. Internal group consistency
    # -------------------------------------------------------------------------

    print("\n[1/5] Internal Step02F group-consistency QC")

    groups = load_step02f_groups(group_file)

    issues, class_counts = validate_internal_group_logic(groups)

    print("Total groups:", f"{len(groups):,}")
    for cls, n in class_counts.items():
        print(f"  {cls}: {n:,}")
    print("Internal contradictions:", f"{len(issues):,}")

    # -------------------------------------------------------------------------
    # 2. Independent remapping
    # -------------------------------------------------------------------------

    print("\n[2/5] Independent remapping of all discriminatory peptides")

    discriminatory = gather_discriminatory_peptides(groups)

    print(
        "Distinct discriminatory peptides to re-check:",
        f"{len(discriminatory):,}"
    )

    remap_rows, remap_mismatches = independent_remap(
        discriminatory,
        fasta,
    )

    remap_lookup = {
        r["peptide"]: r["independent_observed_class"]
        for r in remap_rows
    }

    remap_counts = Counter(
        r["independent_observed_class"]
        for r in remap_rows
    )

    print(
        "Independent remap:",
        f"target-only={remap_counts['TARGET_ONLY']:,} | "
        f"entrapment-only={remap_counts['ENTRAPMENT_ONLY']:,} | "
        f"shared={remap_counts['SHARED_TARGET_ENTRAPMENT']:,} | "
        f"absent={remap_counts['ABSENT_FROM_FROZEN_FASTA']:,}"
    )
    print("Step02F-vs-independent mismatches:", f"{len(remap_mismatches):,}")

    # -------------------------------------------------------------------------
    # 3. Native MM peptide-field provenance
    # -------------------------------------------------------------------------

    print("\n[3/5] Re-reading native MM Unique/Shared Peptides fields")

    native_provenance = load_native_mm_peptide_provenance(mm_native)

    print("Native MM rows indexed:", f"{len(native_provenance):,}")

    # -------------------------------------------------------------------------
    # 4. Group support audit
    # -------------------------------------------------------------------------

    print("\n[4/5] Group-by-group support audit")

    support_rows, support_failures = audit_group_support(
        groups,
        remap_lookup,
        native_provenance,
    )

    entrap_rows = [
        r for r in support_rows
        if r["reconstructed_group_class"]
        == "ENTRAPMENT_ONLY_EVIDENCE"
    ]
    target_rows = [
        r for r in support_rows
        if r["reconstructed_group_class"]
        == "TARGET_ONLY_EVIDENCE"
    ]

    entrap_single = sum(
        1 for r in entrap_rows
        if r["n_discriminatory_peptides"] == 1
    )
    entrap_multi = sum(
        1 for r in entrap_rows
        if r["n_discriminatory_peptides"] >= 2
    )

    target_single = sum(
        1 for r in target_rows
        if r["n_discriminatory_peptides"] == 1
    )
    target_multi = sum(
        1 for r in target_rows
        if r["n_discriminatory_peptides"] >= 2
    )

    print(
        "Entrapment groups:",
        f"{len(entrap_rows):,} | "
        f"single-peptide={entrap_single:,} | "
        f"multi-peptide={entrap_multi:,}"
    )
    print(
        "Target groups:",
        f"{len(target_rows):,} | "
        f"single-peptide={target_single:,} | "
        f"multi-peptide={target_multi:,}"
    )
    print("Group technical failures:", f"{support_failures:,}")

    # -------------------------------------------------------------------------
    # 5. Final verdict + outputs
    # -------------------------------------------------------------------------

    print("\n[5/5] Final validation verdict")

    hard_fail_reasons = []

    if len(groups) != EXPECTED_GROUPS:
        hard_fail_reasons.append(
            f"Expected {EXPECTED_GROUPS} groups, observed {len(groups)}"
        )

    if class_counts["TARGET_ONLY_EVIDENCE"] != EXPECTED_TARGET_GROUPS:
        hard_fail_reasons.append(
            "Target-only group count changed: "
            f"{class_counts['TARGET_ONLY_EVIDENCE']}"
        )

    if class_counts["ENTRAPMENT_ONLY_EVIDENCE"] != EXPECTED_ENTRAPMENT_GROUPS:
        hard_fail_reasons.append(
            "Entrapment-only group count changed: "
            f"{class_counts['ENTRAPMENT_ONLY_EVIDENCE']}"
        )

    if issues:
        hard_fail_reasons.append(
            f"{len(issues)} internal group-classification contradictions"
        )

    if remap_mismatches:
        hard_fail_reasons.append(
            f"{len(remap_mismatches)} independent peptide-remap mismatches"
        )

    if remap_counts["ABSENT_FROM_FROZEN_FASTA"] > 0:
        hard_fail_reasons.append(
            f"{remap_counts['ABSENT_FROM_FROZEN_FASTA']} discriminatory peptides "
            "absent from frozen FASTA"
        )

    if remap_counts["SHARED_TARGET_ENTRAPMENT"] > 0:
        hard_fail_reasons.append(
            f"{remap_counts['SHARED_TARGET_ENTRAPMENT']} supposedly discriminative "
            "peptides independently map to both target and entrapment"
        )

    if support_failures > 0:
        hard_fail_reasons.append(
            f"{support_failures} group-level technical support failures"
        )

    if any(r["n_independently_confirmed"] < 1 for r in entrap_rows):
        hard_fail_reasons.append(
            "At least one entrapment group lacks an independently confirmed "
            "entrapment-only peptide"
        )

    if any(r["n_independently_confirmed"] < 1 for r in target_rows):
        hard_fail_reasons.append(
            "At least one target group lacks an independently confirmed "
            "target-only peptide"
        )

    verdict = (
        "PASS_RECONSTRUCTION_VALIDATED"
        if not hard_fail_reasons
        else "FAIL_RECONSTRUCTION_VALIDATION"
    )

    print("VERDICT:", verdict)

    if hard_fail_reasons:
        for x in hard_fail_reasons:
            print("  FAIL:", x)

    # Write outputs
    write_tsv(
        outdir / "Step02G_internal_consistency_issues.tsv",
        issues,
        ["row_number", "protein_accession", "class", "reason"],
    )

    write_tsv(
        outdir / "Step02G_independent_peptide_remap.tsv",
        remap_rows,
        [
            "peptide",
            "peptide_IL_equivalent",
            "expected_step02f_class",
            "independent_observed_class",
            "match",
            "found_target",
            "found_entrapment",
            "first_target_header",
            "first_entrapment_header",
        ],
    )

    write_tsv(
        outdir / "Step02G_independent_remap_mismatches.tsv",
        remap_mismatches,
        [
            "peptide",
            "peptide_IL_equivalent",
            "expected_step02f_class",
            "independent_observed_class",
            "match",
            "found_target",
            "found_entrapment",
            "first_target_header",
            "first_entrapment_header",
        ],
    )

    write_tsv(
        outdir / "Step02G_group_support_audit.tsv",
        support_rows,
        [
            "row_number",
            "protein_accession_step02f",
            "protein_accession_native",
            "reconstructed_group_class",
            "required_discriminatory_class",
            "n_discriminatory_peptides",
            "n_independently_confirmed",
            "n_from_native_unique_field",
            "n_from_native_shared_field",
            "n_present_in_both_native_fields",
            "n_missing_from_native_peptide_fields",
            "support_tier",
            "technical_pass",
            "bad_peptides",
            "discriminatory_peptides",
        ],
    )

    write_tsv(
        outdir / "Step02G_entrapment_groups_audit.tsv",
        entrap_rows,
        [
            "row_number",
            "protein_accession_step02f",
            "protein_accession_native",
            "reconstructed_group_class",
            "required_discriminatory_class",
            "n_discriminatory_peptides",
            "n_independently_confirmed",
            "n_from_native_unique_field",
            "n_from_native_shared_field",
            "n_present_in_both_native_fields",
            "n_missing_from_native_peptide_fields",
            "support_tier",
            "technical_pass",
            "bad_peptides",
            "discriminatory_peptides",
        ],
    )

    summary = {
        "script_version": VERSION,
        "verdict": verdict,
        "total_groups": len(groups),
        "target_only_groups": class_counts["TARGET_ONLY_EVIDENCE"],
        "entrapment_only_groups": class_counts["ENTRAPMENT_ONLY_EVIDENCE"],
        "mixed_groups": class_counts["MIXED_TARGET_ENTRAPMENT_EVIDENCE"],
        "ambiguous_groups": class_counts["AMBIGUOUS_NO_DISCRIMINATIVE_PEPTIDE"],
        "internal_contradictions": len(issues),
        "distinct_discriminatory_peptides": len(discriminatory),
        "independent_remap_mismatches": len(remap_mismatches),
        "independent_target_only_peptides": remap_counts["TARGET_ONLY"],
        "independent_entrapment_only_peptides": remap_counts["ENTRAPMENT_ONLY"],
        "independent_shared_peptides": remap_counts["SHARED_TARGET_ENTRAPMENT"],
        "independent_absent_peptides": remap_counts["ABSENT_FROM_FROZEN_FASTA"],
        "group_technical_failures": support_failures,
        "entrapment_groups_total": len(entrap_rows),
        "entrapment_groups_single_peptide_support": entrap_single,
        "entrapment_groups_multi_peptide_support": entrap_multi,
        "target_groups_single_peptide_support": target_single,
        "target_groups_multi_peptide_support": target_multi,
        "hard_fail_reasons": " | ".join(hard_fail_reasons),
        "scientific_status":
            "POST_HOC_SEQUENCE_RECONSTRUCTION_VALIDATION; "
            "NATIVE_MM_PROTEIN_TABLE_DID_NOT_RETAIN_ENTRAPMENT_LABEL",
    }

    write_tsv(
        outdir / "Step02G_validation_summary.tsv",
        [summary],
        list(summary.keys()),
    )

    (outdir / "Step02G_method_record.json").write_text(
        json.dumps(
            {
                "script_version": VERSION,
                "step02f_group_file": str(group_file),
                "step02f_peptide_mapping_file": str(peptide_file),
                "frozen_fasta": str(fasta),
                "native_mm_table": str(mm_native),
                "validation_layers": [
                    "internal group-consistency QC",
                    "independent discriminatory-peptide remapping against frozen FASTA",
                    "native MM Unique/Shared peptide-field provenance",
                    "group-by-group technical support audit",
                ],
                "pass_rule":
                    "Zero internal contradictions, zero independent-remap mismatches, "
                    "zero absent/shared discriminatory peptides, zero group technical "
                    "failures, and >=1 independently confirmed discriminatory peptide "
                    "for every target-only and entrapment-only reconstructed group.",
                "scientific_limitation":
                    "A PASS validates the post-hoc sequence reconstruction technically. "
                    "It does not recreate the lost _p_target identifier in the native "
                    "MetaMorpheus protein-group output.",
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    report_lines = [
        f"Step02G VALIDATE MM RECONSTRUCTION v{VERSION}",
        "",
        f"VERDICT: {verdict}",
        "",
        "GROUP COUNTS",
        f"  total groups: {len(groups)}",
        f"  target-only groups: {class_counts['TARGET_ONLY_EVIDENCE']}",
        f"  entrapment-only groups: {class_counts['ENTRAPMENT_ONLY_EVIDENCE']}",
        f"  mixed groups: {class_counts['MIXED_TARGET_ENTRAPMENT_EVIDENCE']}",
        f"  ambiguous groups: {class_counts['AMBIGUOUS_NO_DISCRIMINATIVE_PEPTIDE']}",
        "",
        "INTERNAL CONSISTENCY",
        f"  contradictions: {len(issues)}",
        "",
        "INDEPENDENT FASTA REMAP",
        f"  distinct discriminatory peptides: {len(discriminatory)}",
        f"  target-only: {remap_counts['TARGET_ONLY']}",
        f"  entrapment-only: {remap_counts['ENTRAPMENT_ONLY']}",
        f"  shared target/entrapment: {remap_counts['SHARED_TARGET_ENTRAPMENT']}",
        f"  absent from frozen FASTA: {remap_counts['ABSENT_FROM_FROZEN_FASTA']}",
        f"  Step02F-vs-independent mismatches: {len(remap_mismatches)}",
        "",
        "ENTRAPMENT GROUP SUPPORT",
        f"  total entrapment groups: {len(entrap_rows)}",
        f"  single-peptide support: {entrap_single}",
        f"  multi-peptide support: {entrap_multi}",
        "",
        "TARGET GROUP SUPPORT",
        f"  total target groups: {len(target_rows)}",
        f"  single-peptide support: {target_single}",
        f"  multi-peptide support: {target_multi}",
        "",
        f"GROUP TECHNICAL FAILURES: {support_failures}",
        "",
        "SCIENTIFIC STATUS",
        "  This validates the post-hoc sequence reconstruction against the exact",
        "  frozen entrapment FASTA. The native MetaMorpheus protein table itself",
        "  still did not retain the _p_target label.",
    ]

    if hard_fail_reasons:
        report_lines.extend([
            "",
            "FAIL REASONS",
            *[f"  - {x}" for x in hard_fail_reasons],
        ])

    report_path = outdir / "Step02G_VALIDATION_REPORT.txt"
    report_path.write_text(
        "\n".join(report_lines),
        encoding="utf-8",
    )

    print("\n" + "=" * 100)
    print("STEP02G COMPLETE")
    print("=" * 100)
    print("Send back:")
    print(report_path)
    print(outdir / "Step02G_validation_summary.tsv")
    print(outdir / "Step02G_entrapment_groups_audit.tsv")
    print("=" * 100)


if __name__ == "__main__":
    main()
