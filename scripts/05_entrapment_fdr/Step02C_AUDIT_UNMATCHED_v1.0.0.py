#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step02C_AUDIT_UNMATCHED_v1.0.0.py

Purpose
-------
Audit the "unmatched" peptide keys reported by Step02B before making any FDR claim.

For every unmatched I/L-equivalent peptide key:
1) Search the frozen target+entrapment FASTA directly.
2) Determine whether the sequence is present anywhere in the FASTA.
3) If present in FASTA but absent from the frozen tryptic peptide-space lookup,
   classify it as a digestion/search-space mismatch candidate.
4) If absent from the frozen target+entrapment FASTA,
   classify it as external-to-database (commonly contaminant/external sequence,
   but this script does NOT automatically call it a contaminant).

This is QC/data processing only. No figure is generated.

Expected input root
-------------------
C:\\Users\\Supachai\\Desktop\\AphaPept_benchmark\\Benchmark_Program\\ENTRAPMENT_FDR

Expected Step02B input
----------------------
Step02B_HOMEWORK_CHECK_v100\\Step02B_discovery_classes.tsv

Outputs
-------
Step02C_UNMATCHED_AUDIT_v100\\
    Step02C_unmatched_detail.tsv
    Step02C_unmatched_summary.tsv
    Step02C_AUDIT_REPORT.txt
"""

from __future__ import annotations

import argparse
import csv
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

DEFAULT_ROOT = _public_project_root() / "ENTRAPMENT_FDR"

def normalize_il(seq: str) -> str:
    return seq.upper().replace("I", "J").replace("L", "J")

def is_entrapment_header(header: str) -> bool:
    s = header.lower()
    return ("_p_target" in s) or ("p_target" in s) or ("entrap" in s)

def parse_fasta(path: Path):
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

# -------------------------------------------------------------------------
# Small pure-Python Aho-Corasick automaton
# -------------------------------------------------------------------------

class ACAutomaton:
    def __init__(self):
        self.next = [{}]
        self.fail = [0]
        self.out = [[]]

    def add(self, pattern: str):
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

        for ch, state in self.next[0].items():
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

    def find_unique(self, text: str):
        state = 0
        found = set()

        for ch in text:
            while state and ch not in self.next[state]:
                state = self.fail[state]

            state = self.next[state].get(ch, 0)

            if self.out[state]:
                found.update(self.out[state])

        return found

def locate_step02b_discovery(root: Path) -> Path:
    direct = (
        root
        / "Step02B_HOMEWORK_CHECK_v100"
        / "Step02B_discovery_classes.tsv"
    )
    if direct.exists():
        return direct

    hits = list(root.rglob("Step02B_discovery_classes.tsv"))
    if not hits:
        raise FileNotFoundError(
            "Step02B_discovery_classes.tsv was not found under:\n"
            f"{root}"
        )

    return sorted(hits, key=lambda p: len(str(p)))[0]

def locate_fasta(root: Path) -> Path:
    direct = (
        root
        / "Step02A_entrapment_db_v104"
        / "Step02_target_plus_shuffled_entrapment_r1.fasta"
    )
    if direct.exists():
        return direct

    hits = list(
        root.rglob("Step02_target_plus_shuffled_entrapment_r1.fasta")
    )
    if not hits:
        raise FileNotFoundError(
            "Frozen target+entrapment FASTA was not found under:\n"
            f"{root}"
        )

    return sorted(hits, key=lambda p: len(str(p)))[0]

def load_unmatched(path: Path):
    by_workflow = {}
    all_patterns = set()

    with path.open(
        "r",
        encoding="utf-8-sig",
        errors="replace",
        newline=""
    ) as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        required = {
            "workflow",
            "class",
            "peptide_key_IL_equivalent",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Missing columns in {path.name}: {sorted(missing)}"
            )

        for row in reader:
            if row["class"] != "unmatched":
                continue

            wf = row["workflow"].strip()
            pep = row["peptide_key_IL_equivalent"].strip().upper()

            if not pep:
                continue

            by_workflow.setdefault(wf, set()).add(pep)
            all_patterns.add(pep)

    return by_workflow, all_patterns

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--root",
        default=str(DEFAULT_ROOT),
        help="ENTRAPMENT_FDR root folder",
    )
    ap.add_argument(
        "--outdir",
        default=None,
        help="Optional output folder",
    )
    args = ap.parse_args()

    root = Path(args.root)

    if not root.exists():
        raise SystemExit(f"Root does not exist:\n{root}")

    discovery_file = locate_step02b_discovery(root)
    fasta = locate_fasta(root)

    outdir = (
        Path(args.outdir)
        if args.outdir
        else root / "Step02C_UNMATCHED_AUDIT_v100"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 100)
    print(f"STEP02C UNMATCHED AUDIT v{VERSION}")
    print("=" * 100)
    print("Root:", root)
    print("Step02B classes:", discovery_file)
    print("Frozen FASTA:", fasta)

    by_wf, patterns = load_unmatched(discovery_file)

    print("\n[1/3] Unmatched keys loaded")
    for wf in ["AP", "FP", "MM", "MQ"]:
        print(f"{wf}: {len(by_wf.get(wf, set())):,}")
    print(f"Union unmatched keys: {len(patterns):,}")

    if not patterns:
        print("\nNo unmatched peptide keys. Nothing to audit.")
        return

    print("\n[2/3] Building peptide matcher")
    ac = ACAutomaton()
    for pep in sorted(patterns):
        ac.add(pep)
    ac.build()

    matches = {
        pep: {
            "found_target": False,
            "found_entrapment": False,
            "first_target_header": "",
            "first_entrapment_header": "",
        }
        for pep in patterns
    }

    print("[3/3] Scanning frozen target+entrapment FASTA once")

    n_proteins = 0
    for header, seq in parse_fasta(fasta):
        n_proteins += 1

        # FASTA sequence is normalized I/L before matching because Step02B
        # discovery keys were normalized the same way.
        seq_il = normalize_il(seq)
        found = ac.find_unique(seq_il)

        if not found:
            continue

        ent = is_entrapment_header(header)

        for pep in found:
            rec = matches[pep]

            if ent:
                rec["found_entrapment"] = True
                if not rec["first_entrapment_header"]:
                    rec["first_entrapment_header"] = header
            else:
                rec["found_target"] = True
                if not rec["first_target_header"]:
                    rec["first_target_header"] = header

        if n_proteins % 25000 == 0:
            print(f"  scanned {n_proteins:,} FASTA entries")

    detail_rows = []
    summary = []

    for wf in ["AP", "FP", "MM", "MQ"]:
        keys = sorted(by_wf.get(wf, set()))
        counts = Counter()

        for pep in keys:
            rec = matches[pep]

            if rec["found_target"] and rec["found_entrapment"]:
                location = "PRESENT_IN_TARGET_AND_ENTRAPMENT_FASTA"
                interpretation = (
                    "Sequence occurs in both FASTA classes but was outside "
                    "the frozen tryptic lookup; inspect digestion/search rules."
                )
            elif rec["found_target"]:
                location = "PRESENT_IN_TARGET_FASTA"
                interpretation = (
                    "Sequence exists in real-target proteins but was outside "
                    "the frozen tryptic lookup; likely digestion/search-space "
                    "rule mismatch (e.g. cleavage, terminal processing, length)."
                )
            elif rec["found_entrapment"]:
                location = "PRESENT_IN_ENTRAPMENT_FASTA"
                interpretation = (
                    "Sequence exists in entrapment proteins but was outside "
                    "the frozen tryptic lookup; likely digestion/search-space "
                    "rule mismatch."
                )
            else:
                location = "ABSENT_FROM_FROZEN_FASTA"
                interpretation = (
                    "Sequence is external to the frozen target+entrapment "
                    "database. Common causes include software-added "
                    "contaminants or another external sequence source. "
                    "Do not automatically label it contaminant without "
                    "checking the native output metadata."
                )

            counts[location] += 1

            detail_rows.append({
                "workflow": wf,
                "peptide_key_IL_equivalent": pep,
                "length": len(pep),
                "location_class": location,
                "found_in_target_fasta": rec["found_target"],
                "found_in_entrapment_fasta": rec["found_entrapment"],
                "first_target_header": rec["first_target_header"],
                "first_entrapment_header": rec["first_entrapment_header"],
                "interpretation": interpretation,
            })

        total = len(keys)

        summary.append({
            "workflow": wf,
            "total_unmatched": total,
            "present_target_only_or_target_plus_entrapment": (
                counts["PRESENT_IN_TARGET_FASTA"]
                + counts["PRESENT_IN_TARGET_AND_ENTRAPMENT_FASTA"]
            ),
            "present_entrapment_only": counts["PRESENT_IN_ENTRAPMENT_FASTA"],
            "absent_from_frozen_fasta": counts["ABSENT_FROM_FROZEN_FASTA"],
            "pct_absent_from_frozen_fasta": (
                100.0 * counts["ABSENT_FROM_FROZEN_FASTA"] / total
                if total else 0.0
            ),
        })

    detail_path = outdir / "Step02C_unmatched_detail.tsv"
    with detail_path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as fh:
        fields = [
            "workflow",
            "peptide_key_IL_equivalent",
            "length",
            "location_class",
            "found_in_target_fasta",
            "found_in_entrapment_fasta",
            "first_target_header",
            "first_entrapment_header",
            "interpretation",
        ]
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(detail_rows)

    summary_path = outdir / "Step02C_unmatched_summary.tsv"
    with summary_path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as fh:
        fields = [
            "workflow",
            "total_unmatched",
            "present_target_only_or_target_plus_entrapment",
            "present_entrapment_only",
            "absent_from_frozen_fasta",
            "pct_absent_from_frozen_fasta",
        ]
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary)

    report_lines = [
        f"Step02C UNMATCHED AUDIT v{VERSION}",
        f"Frozen FASTA entries scanned: {n_proteins}",
        "",
        "workflow\ttotal_unmatched\tpresent_target\tpresent_entrapment_only\tabsent_from_frozen_fasta\tpct_absent",
    ]

    for row in summary:
        report_lines.append(
            f"{row['workflow']}\t"
            f"{row['total_unmatched']}\t"
            f"{row['present_target_only_or_target_plus_entrapment']}\t"
            f"{row['present_entrapment_only']}\t"
            f"{row['absent_from_frozen_fasta']}\t"
            f"{row['pct_absent_from_frozen_fasta']:.2f}"
        )

    report_path = outdir / "Step02C_AUDIT_REPORT.txt"
    report_path.write_text(
        "\n".join(report_lines),
        encoding="utf-8",
    )

    print("\n" + "=" * 100)
    print("STEP02C COMPLETE")
    print("=" * 100)

    for row in summary:
        print(
            f"{row['workflow']}: unmatched={row['total_unmatched']:,} | "
            f"present in target={row['present_target_only_or_target_plus_entrapment']:,} | "
            f"entrapment-only={row['present_entrapment_only']:,} | "
            f"absent from FASTA={row['absent_from_frozen_fasta']:,} "
            f"({row['pct_absent_from_frozen_fasta']:.1f}%)"
        )

    print("\nSend back these two small files:")
    print(report_path)
    print(summary_path)
    print("=" * 100)

if __name__ == "__main__":
    main()
