#!/usr/bin/env python3
# Python 3.9 compatible

import argparse
import csv
from pathlib import Path

PRIVATE = "C:" + chr(92) + "Users" + chr(92) + "Supachai"
TEXT_EXT = {".py", ".ps1", ".r", ".md", ".txt", ".yaml", ".yml", ".json", ".tsv", ".csv"}
VALID = {"COPIED", "COPIED_IDENTICAL_DUPLICATES"}

def read_text(path):
    with open(str(path), "r", encoding="utf-8-sig", errors="replace") as fh:
        return fh.read()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release-root", required=True, type=Path)
    ap.add_argument("--allow-placeholders", action="store_true")
    args = ap.parse_args()

    root = args.release_root.resolve()
    problems = []
    warnings = []

    collection = root / "COLLECTION_REPORT.tsv"
    if not collection.exists():
        problems.append("COLLECTION_REPORT.tsv missing.")
    else:
        with open(str(collection), "r", encoding="utf-8-sig", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        for row in rows:
            if (
                (row.get("requirement") or "").strip() == "required"
                and (row.get("status") or "").strip() not in VALID
            ):
                problems.append(
                    "Unresolved required script: %s [%s]"
                    % (row.get("source_pattern"), row.get("status"))
                )

    syntax_failures = []
    for path in sorted((root / "scripts").rglob("*.py")):
        try:
            compile(read_text(path), str(path), "exec")
        except Exception as exc:
            syntax_failures.append((str(path.relative_to(root)), repr(exc)))

    if syntax_failures:
        problems.append(
            "%d Python source file(s) have syntax errors."
            % len(syntax_failures)
        )

    private_hits = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(root)
        if rel.parts and rel.parts[0] == "provenance":
            continue
        if path.name in {
            "COLLECTION_REPORT.tsv",
            "PATH_AUDIT.tsv",
            "PORTABILITY_REPORT.tsv",
            "PORTABILITY_REPORT_v1.0.3.tsv",
            "PORTABILITY_REPORT_v1.0.4.tsv",
            "REPAIR_REPORT_v1.1.0.tsv",
            "SYNTAX_AUDIT_v1.1.0.tsv",
            "PRIVATE_PATH_AUDIT_v1.1.0.tsv",
        }:
            continue
        if path.suffix.lower() not in TEXT_EXT:
            continue
        for i, line in enumerate(read_text(path).splitlines(), 1):
            if PRIVATE in line:
                private_hits.append((str(rel), i, line.strip()))

    if private_hits:
        problems.append(
            "%d private user-path match(es) remain."
            % len(private_hits)
        )

    cav = root / "CODE_AVAILABILITY.md"
    placeholders = []
    if cav.exists():
        t = read_text(cav)
        for label, token in [
            ("GitHub", "<GITHUB_REPOSITORY_URL>"),
            ("Zenodo", "<ZENODO_DOI_OR_ARCHIVE_URL>"),
            ("data accession", "<ACCESSION>"),
        ]:
            if token in t:
                placeholders.append(label)

    if placeholders:
        msg = "Placeholders remain: " + ", ".join(placeholders)
        if args.allow_placeholders:
            warnings.append(msg)
        else:
            problems.append(msg)

    print("=" * 78)
    print("PUBLIC RELEASE VALIDATION v1.1.0")
    print("=" * 78)
    print("Release root:", root)
    print("Python syntax failures:", len(syntax_failures))
    print("Private path matches:", len(private_hits))

    if syntax_failures:
        print("\nSYNTAX FAILURES:")
        for name, err in syntax_failures:
            print(" -", name)
            print("   ", err)

    if private_hits:
        print("\nPRIVATE PATHS:")
        for name, line_no, line in private_hits[:30]:
            print(" - %s:%d: %s" % (name, line_no, line))

    if warnings:
        print("\nWARNINGS:")
        for x in warnings:
            print(" -", x)

    if problems:
        print("\nNOT READY:")
        for x in problems:
            print(" -", x)
        raise SystemExit(2)

    if args.allow_placeholders:
        print("\nPASS (PRE-PUBLIC): code, syntax, and path checks pass.")
    else:
        print("\nPASS: code and release metadata checks pass.")

if __name__ == "__main__":
    main()
