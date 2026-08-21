#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

import argparse
import csv
import json
import sqlite3
from collections import Counter
from pathlib import Path

STANDARD = set("ACDEFGHIKLMNPQRSTVWY")

EXPECTED = {
    "target_proteins": 169637,
    "entrapment_proteins": 169637,
    "target_only_peptide_keys": 3158816,
    "entrapment_only_peptide_keys": 3169203,
    "shared_target_entrapment_peptide_keys": 5507,
}

def parse_fasta(path: Path):
    header = None
    seq = []
    with path.open("rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq).upper()
                header = line[1:].strip()
                seq = []
            elif header is not None:
                seq.append(line.strip())
        if header is not None:
            yield header, "".join(seq).upper()

def is_entrapment_header(header: str) -> bool:
    s = header.lower()
    return ("_p_target" in s) or ("p_target" in s) or ("entrap" in s)

def normalize_il(seq: str) -> str:
    return seq.replace("I", "J").replace("L", "J")

def standard_segments(seq: str):
    start = 0
    for i, aa in enumerate(seq):
        if aa not in STANDARD:
            if i > start:
                yield seq[start:i]
            start = i + 1
    if start < len(seq):
        yield seq[start:]

def tryptic_boundaries(seq: str):
    cuts = [0]
    for i, aa in enumerate(seq):
        if aa in "KR" and (i + 1 == len(seq) or seq[i + 1] != "P"):
            cuts.append(i + 1)
    if cuts[-1] != len(seq):
        cuts.append(len(seq))
    return cuts

def digest(seq: str, missed: int = 2, min_len: int = 7, max_len: int = 50):
    for segment in standard_segments(seq):
        cuts = tryptic_boundaries(segment)
        n = len(cuts) - 1
        for i in range(n):
            for mc in range(missed + 1):
                j = i + mc + 1
                if j > n:
                    break
                pep = segment[cuts[i]:cuts[j]]
                if min_len <= len(pep) <= max_len:
                    yield normalize_il(pep)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    fasta = Path(args.fasta)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    db = outdir / "Step02_entrapment_peptide_space.sqlite"
    if db.exists():
        db.unlink()

    con = sqlite3.connect(str(db))
    cur = con.cursor()
    cur.execute("PRAGMA journal_mode=WAL")
    cur.execute("PRAGMA synchronous=NORMAL")
    cur.execute(
        "CREATE TABLE peptide_space ("
        "peptide_key TEXT PRIMARY KEY, "
        "is_target INTEGER NOT NULL DEFAULT 0, "
        "is_entrapment INTEGER NOT NULL DEFAULT 0)"
    )

    proteins = Counter()
    nonstandard = Counter()
    batch = []

    def flush():
        nonlocal batch
        if not batch:
            return
        cur.executemany(
            "INSERT INTO peptide_space(peptide_key,is_target,is_entrapment) "
            "VALUES(?,?,?) "
            "ON CONFLICT(peptide_key) DO UPDATE SET "
            "is_target=MAX(is_target,excluded.is_target), "
            "is_entrapment=MAX(is_entrapment,excluded.is_entrapment)",
            batch,
        )
        con.commit()
        batch = []

    for header, seq in parse_fasta(fasta):
        ent = is_entrapment_header(header)
        proteins["entrapment" if ent else "target"] += 1

        for aa in set(seq) - STANDARD:
            nonstandard[aa] += 1

        t, e = ((0, 1) if ent else (1, 0))
        local = set()

        for key in digest(seq):
            if key in local:
                continue
            local.add(key)
            batch.append((key, t, e))
            if len(batch) >= 100000:
                flush()

    flush()

    counts = {}
    queries = {
        "target_only": "is_target=1 AND is_entrapment=0",
        "entrapment_only": "is_target=0 AND is_entrapment=1",
        "shared": "is_target=1 AND is_entrapment=1",
        "union": "1=1",
    }

    for label, where in queries.items():
        cur.execute(f"SELECT COUNT(*) FROM peptide_space WHERE {where}")
        counts[label] = int(cur.fetchone()[0])

    con.close()

    r = (
        counts["entrapment_only"] / counts["target_only"]
        if counts["target_only"] else float("nan")
    )

    summary = {
        "target_proteins": int(proteins["target"]),
        "entrapment_proteins": int(proteins["entrapment"]),
        "target_only_peptide_keys": counts["target_only"],
        "entrapment_only_peptide_keys": counts["entrapment_only"],
        "shared_target_entrapment_peptide_keys": counts["shared"],
        "union_peptide_keys": counts["union"],
        "effective_r_entrapment_over_target": r,
    }

    qc = []
    for key, expected in EXPECTED.items():
        observed = summary[key]
        qc.append({
            "metric": key,
            "expected_v103": expected,
            "observed_v104": observed,
            "match": observed == expected,
        })

    qc_ok = all(x["match"] for x in qc)

    with (outdir / "Step02A_space_summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        for key, value in summary.items():
            w.writerow([key, value])

    with (outdir / "Step02A_reproduction_QC.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["metric", "expected_v103", "observed_v104", "match"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(qc)

    manifest = {
        "database_fasta": str(fasta.resolve()),
        "lookup_digest": {
            "enzyme": "trypsin; K/R except before P",
            "missed_cleavages": 2,
            "min_length": 7,
            "max_length": 50,
            "IL_equivalent": True,
        },
        "protein_counts": dict(proteins),
        "nonstandard_residue_entry_occurrences": dict(nonstandard),
        "peptide_space_counts": counts,
        "effective_r_entrapment_over_target": r,
        "matches_previous_v103_QC": qc_ok,
    }

    (outdir / "Step02A_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    print("\nSTEP02A LOOKUP COMPLETE")
    for key, value in summary.items():
        print(f"{key}: {value}")

    print("\nREPRODUCTION QC AGAINST PREVIOUS v103")
    for row in qc:
        flag = "PASS" if row["match"] else "FAIL"
        print(
            f"{flag:4s} {row['metric']}: "
            f"expected={row['expected_v103']} observed={row['observed_v104']}"
        )

    if not qc_ok:
        raise SystemExit(
            "\nQC mismatch versus the previously successful v103 database. "
            "Do NOT rerun AP/FP/MM/MQ until this is resolved."
        )

    print("\nALL FROZEN v103 QC COUNTS REPRODUCED.")
    print("Database is ready for AP / FP / MM / MQ entrapment reruns.")

if __name__ == "__main__":
    main()
