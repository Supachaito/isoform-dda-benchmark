#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ISOFORM_KMIN_INFERENCE_v1.1.0.py

Relax isoform inference resolution without relaxing peptide confidence.
Uses manuscript-defined primary isoform-discriminative evidence:
  1) single_isoform_unique
  2) within_family_subset_discriminative

For each workflow x run x UniProt isoform family, observed peptide-to-compatible
entry sets are classified as EXACT_RESOLVED, SUBSET_RESOLVED, or
MULTIPLE_ISOFORMS_REQUIRED using an exact minimum hitting-set (Kmin) analysis.

Cell-line-level inference requires the same inference in >=2/3 replicates.
1/3 observations remain SPARSE. No imputation. MBR OFF only.
"""

from __future__ import annotations

import os
import re
import json
from pathlib import Path
from functools import lru_cache
from collections import Counter
from typing import Iterable

import numpy as np
import pandas as pd
import h5py

VERSION = "1.1.0"

# ----------------------------------------------------------------------
# ROBUST PROJECT PATH DISCOVERY
#
# Recommended location for this script:
# <PROJECT_ROOT>\ENTRAPMENT_FDR
#
# The script does NOT require the current PowerShell folder.
# It derives Benchmark_Program from the location of this .py file.
# ----------------------------------------------------------------------

SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT_PATH.parent

# Walk upward from the script location until Benchmark_Program is found.
ROOT = None
for candidate in [SCRIPT_DIR, *SCRIPT_DIR.parents]:
    if candidate.name == "Benchmark_Program":
        ROOT = candidate
        break

if ROOT is None:
    raise RuntimeError(
        "Could not locate the Benchmark_Program folder from this script location.\n"
        f"Script:\n{SCRIPT_PATH}\n\n"
        "Place this script anywhere inside:\n"
        r"<PROJECT_ROOT>"
    )

# Frozen normalized peptide mapping.
MASTER = ROOT / "MASTER_TABLES_FINAL"
STEP1D = MASTER / "STEP1D_FINAL_NORMALIZED"

# Manuscript supplementary-figure workspace.
SUPP_ROOT = (
    ROOT
    / "MANUSCRIPT_REVISION_20260813"
    / "MAIN_FIGURES_REBUILD_20260815_V01"
    / "SUPPLEMENTARY_FIGURES"
)

if not SUPP_ROOT.exists():
    raise RuntimeError(
        "SUPPLEMENTARY_FIGURES folder does not exist:\n"
        + str(SUPP_ROOT)
    )

# Existing isoform overview workspace.
ISOFORM_OVERVIEW = SUPP_ROOT / "SUPP_FIG_ISOFORM_OVERVIEW_V01"
ISOFORM_OVERVIEW_DATA = ISOFORM_OVERVIEW / "DATA"

# New Kmin analysis output.
OUT = SUPP_ROOT / "ISOFORM_KMIN_V01"
OUT.mkdir(parents=True, exist_ok=True)

print("Script:")
print(SCRIPT_PATH)
print("\nDerived project root:")
print(ROOT)
print("\nOutput folder:")
print(OUT)

PROGRAMS = ["AP", "FP", "MM", "MQ"]
PROGRAM_FOLDERS = {
    "AP": ROOT / "AP_MBR_OFF",
    "FP": ROOT / "FP_MBR_OFF_LFQ",
    "MM": ROOT / "MM_MBR_OFF",
    "MQ": ROOT / "MQ_MBR_OFF",
}
CELL_LINES = ["C33A", "SiHa", "HeLa"]
REPLICATES = [1, 2, 3]
SAMPLES = [f"{cell}_{rep}" for cell in CELL_LINES for rep in REPLICATES]
DISK_CELL = {"C33A": "C33A", "SiHa": "SIHA", "HeLa": "HELA"}

EXPECTED_FROZEN_COUNTS = {"AP": 23767, "FP": 22674, "MM": 24543, "MQ": 18962}
EXPECTED_SINGLE_UNIQUE = {"AP": 69, "FP": 13, "MM": 42, "MQ": 14}
EXPECTED_PRIMARY = {"AP": 180, "FP": 101, "MM": 149, "MQ": 70}
PRIMARY_CATEGORIES = {"single_isoform_unique", "within_family_subset_discriminative"}


def die(message: str) -> None:
    raise RuntimeError(message)


def clean_peptide(x: object) -> str:
    if x is None or pd.isna(x):
        return ""
    s = str(x).upper().strip()
    if "|" in s:
        return ""
    m = re.match(r"^[A-Z\-]\.([^\.]+)\.[A-Z\-]$", s)
    if m:
        s = m.group(1)
    return "".join(re.findall(r"[A-Z]", s))


def il_key(x: object) -> str:
    return clean_peptide(x).replace("I", "L")


def normalize_program(x: object) -> str:
    s = str(x).upper().strip()
    lookup = {
        "ALPHAPEPT": "AP", "AP": "AP",
        "FRAGPIPE": "FP", "MSFRAGGER": "FP", "FRAGPIPE/MSFRAGGER": "FP", "FP": "FP",
        "METAMORPHEUS": "MM", "MM": "MM",
        "MAXQUANT": "MQ", "MQ": "MQ",
    }
    return lookup.get(s, s)


def find_col(columns: Iterable[object], candidates: Iterable[str], contains: bool = False) -> str | None:
    cols = list(columns)
    exact = {str(c).strip().lower(): str(c) for c in cols}
    for cand in candidates:
        key = cand.strip().lower()
        if key in exact:
            return exact[key]
    if contains:
        for c in cols:
            lc = str(c).strip().lower()
            if any(cand.strip().lower() in lc for cand in candidates):
                return str(c)
    return None


def sample_from_text(x: object) -> str | None:
    if x is None:
        return None
    m = re.search(r"(?i)(C33A|HELA|SIHA)[\s_\-\.]*([123])(?:\D|$)", str(x))
    if not m:
        return None
    cell = {"C33A": "C33A", "HELA": "HeLa", "SIHA": "SiHa"}[m.group(1).upper()]
    return f"{cell}_{m.group(2)}"


def disk_sample(sample: str) -> str:
    cell, rep = sample.rsplit("_", 1)
    return f"{DISK_CELL[cell]}_{rep}"


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path, low_memory=False)
    if suffix in {".tsv", ".txt", ".psmtsv"}:
        return pd.read_csv(path, sep="\t", low_memory=False)
    die(f"Unsupported table extension:\n{path}")


def decode_hdf_strings(values: np.ndarray) -> list[str]:
    out = []
    for v in values:
        if isinstance(v, (bytes, np.bytes_)):
            out.append(v.decode("utf-8", errors="replace"))
        else:
            out.append(str(v))
    return out


UNIPROT_CORE = r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
PLAIN_ACC_RE = re.compile(rf"(?<![A-Z0-9])({UNIPROT_CORE})(?:-(\d+))?(?![A-Z0-9])")
WRAPPED_ACC_RE = re.compile(rf"(?:sp|tr)\|({UNIPROT_CORE}(?:-\d+)?)\|", re.I)


def extract_accessions(x: object) -> list[str]:
    if x is None or pd.isna(x):
        return []
    s = str(x).upper()
    found = []
    for m in WRAPPED_ACC_RE.finditer(s):
        found.append(m.group(1).upper())
    for m in PLAIN_ACC_RE.finditer(s):
        base = m.group(1).upper()
        iso = m.group(2)
        found.append(f"{base}-{iso}" if iso is not None else base)
    return sorted(set(found))


def base_accession(accession: str) -> str:
    return re.sub(r"-\d+$", "", str(accession))



# ----------------------------------------------------------------------
# Authoritative FASTA-derived isoform-family universe
# ----------------------------------------------------------------------

FASTA_BASENAME = "uniprotkb_proteome_UP000005640_2026_08_04.fasta"


def locate_frozen_fasta(root: Path) -> Path:
    """
    Locate the exact frozen FASTA used for this manuscript.
    Multiple copies are permitted only if byte-identical.
    """
    candidates = sorted(
        {
            p.resolve()
            for p in root.rglob(FASTA_BASENAME)
            if p.is_file()
        }
    )

    if not candidates:
        die(
            "Frozen FASTA not found under project root:\n"
            + str(root)
            + "\nExpected filename:\n"
            + FASTA_BASENAME
        )

    if len(candidates) == 1:
        return candidates[0]

    import hashlib

    def sha256_file(path: Path) -> str:
        h = hashlib.sha256()
        with path.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest()

    hashes = {
        p: sha256_file(p)
        for p in candidates
    }

    if len(set(hashes.values())) != 1:
        detail = "\n".join(
            f"{p}  sha256={h}"
            for p, h in hashes.items()
        )
        die(
            "Multiple non-identical copies of the frozen FASTA were found. "
            "Stopped rather than guessing:\n"
            + detail
        )

    return sorted(
        candidates,
        key=lambda p: (
            len(str(p)),
            str(p).lower(),
        )
    )[0]


def fasta_accession_from_header(header: str) -> str | None:
    """Parse accession from a standard UniProt FASTA header."""
    h = header.strip()

    if not h.startswith(">"):
        return None

    body = h[1:]

    m = re.match(
        r"^(?:sp|tr)\|([^|]+)\|",
        body,
        flags=re.I
    )

    if m:
        return m.group(1).strip()

    token = body.split(None, 1)[0].strip()
    return token if token else None


def build_fasta_family_universe(fasta_path: Path) -> pd.DataFrame:
    """
    Build complete accession-family membership directly from the frozen FASTA.
    This prevents family size from being underestimated when some annotated
    isoforms have no discriminative peptide in the observed mapping.
    """
    rows = []

    with fasta_path.open(
        "r",
        encoding="utf-8",
        errors="replace"
    ) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue

            acc = fasta_accession_from_header(line)

            if not acc:
                continue

            rows.append(
                {
                    "BaseAccession": base_accession(acc),
                    "EntryAccession": acc,
                }
            )

    entries = pd.DataFrame(rows).drop_duplicates()

    if entries.empty:
        die(
            "No FASTA entries could be parsed from:\n"
            + str(fasta_path)
        )

    family = (
        entries
        .groupby(
            "BaseAccession",
            as_index=False
        )
        .agg(
            FastaFamilyEntries=(
                "EntryAccession",
                lambda x: ";".join(
                    sorted(set(map(str, x)))
                )
            ),
            NFastaFamilyEntries=(
                "EntryAccession",
                "nunique"
            ),
        )
    )

    return family


PROGRAM_COL_CANDS = ["Program", "Workflow", "Software", "Method"]
PEPTIDE_COL_CANDS = ["Peptide", "Sequence", "PeptideSequence", "BaseSequence", "Base Sequence", "Peptide_Key", "PeptideKey", "ILKey", "I_L_Key"]
CATEGORY_COL_CANDS = ["Category", "PeptideCategory", "MappingCategory", "GeneAwareCategory", "Peptide_Class", "EvidenceClass", "Class"]
ACCESSION_COL_CANDS = ["CompatibleAccessions", "Compatible Accessions", "CompatibleProteinAccessions", "Compatible Protein Accessions", "ProteinAccessions", "Protein Accessions", "MappedProteins", "Mapped Proteins", "Accessions", "AccessionSet", "Accession Set", "ProteinEntries", "Protein Entries"]


def inspect_mapping_candidate(path: Path) -> dict | None:
    try:
        df = read_table(path)
    except Exception:
        return None
    pcol = find_col(df.columns, PROGRAM_COL_CANDS)
    pepcol = find_col(df.columns, PEPTIDE_COL_CANDS)
    catcol = find_col(df.columns, CATEGORY_COL_CANDS, contains=True)
    acccol = find_col(df.columns, ACCESSION_COL_CANDS, contains=True)
    if not all([pcol, pepcol, catcol, acccol]):
        return None
    z = pd.DataFrame({
        "Program": df[pcol].map(normalize_program),
        "PeptideRaw": df[pepcol],
        "Category": df[catcol].astype(str).str.strip(),
        "AccessionRaw": df[acccol],
    })
    z["PeptideKey"] = z["PeptideRaw"].map(il_key)
    z = z[z["Program"].isin(PROGRAMS) & z["PeptideKey"].ne("")].copy()
    if z.empty:
        return None
    conflicts = z.groupby(["Program", "PeptideKey"])["Category"].nunique(dropna=False)
    counts = z.drop_duplicates(["Program", "PeptideKey"]).groupby("Program")["PeptideKey"].nunique().to_dict()
    return {"path": path, "data": z, "conflicts": int((conflicts > 1).sum()), "counts": counts}


def load_frozen_mapping():
    if not STEP1D.exists():
        die("STEP1D folder does not exist:\n" + str(STEP1D))
    candidates = []
    for p in STEP1D.rglob("*"):
        if p.is_file() and p.suffix.lower() in {".csv", ".tsv", ".txt"}:
            item = inspect_mapping_candidate(p)
            if item is not None:
                candidates.append(item)
    if not candidates:
        die("No valid mapping table found under:\n" + str(STEP1D))
    inventory, exact = [], []
    for c in candidates:
        inventory.append({"File": str(c["path"]), "CategoryConflicts": c["conflicts"], **{p: c["counts"].get(p, 0) for p in PROGRAMS}})
        if c["conflicts"] == 0 and all(c["counts"].get(p, -1) == EXPECTED_FROZEN_COUNTS[p] for p in PROGRAMS):
            exact.append(c)
    pd.DataFrame(inventory).to_csv(OUT / "00_mapping_inventory.csv", index=False)
    if not exact:
        die("No mapping table reproduced all frozen peptide counts. Stopped rather than guessing.")
    if len(exact) > 1:
        signatures = []
        for c in exact:
            sig = c["data"][["Program", "PeptideKey", "Category"]].drop_duplicates().sort_values(["Program", "PeptideKey", "Category"]).reset_index(drop=True)
            signatures.append(sig)
        if not all(signatures[0].equals(x) for x in signatures[1:]):
            die("Multiple frozen mapping tables match counts but differ in normalized content. Stopped rather than choosing.")
        exact.sort(key=lambda c: (len(str(c["path"])), str(c["path"]).lower()))
    chosen = exact[0]
    return chosen["data"].copy(), chosen["path"]


def build_compatibility(
    raw: pd.DataFrame,
    fasta_family: pd.DataFrame
):
    """
    Reconstruct peptide -> compatible-entry sets from the frozen mapping,
    while deriving full accession-family membership from the frozen FASTA.

    v1.1.0 change:
    family size is NOT reconstructed from observed peptide mappings.
    """
    records = []

    for (program, peptide), g in raw.groupby(
        ["Program", "PeptideKey"],
        sort=False
    ):
        categories = sorted(
            set(
                g["Category"]
                .astype(str)
                .str.strip()
            )
        )

        if len(categories) != 1:
            die(
                "Category conflict after frozen-table selection:\n"
                f"{program} {peptide} {categories}"
            )

        category = categories[0]
        accs = set()

        for value in g["AccessionRaw"]:
            accs.update(
                extract_accessions(value)
            )

        if not accs:
            continue

        bases = {
            base_accession(x)
            for x in accs
        }

        # Cross-family peptides are excluded from this family-level analysis.
        if len(bases) != 1:
            continue

        base = next(iter(bases))

        records.append(
            {
                "Program": program,
                "PeptideKey": peptide,
                "Category": category,
                "BaseAccession": base,
                "CompatibleEntries": ";".join(
                    sorted(accs)
                ),
                "NCompatibleEntries": len(accs),
            }
        )

    all_family = pd.DataFrame(records)

    if all_family.empty:
        die(
            "No within-family peptide compatibility records reconstructed."
        )

    # Join complete family membership from the exact frozen FASTA.
    all_family = all_family.merge(
        fasta_family,
        on="BaseAccession",
        how="left",
        validate="many_to_one"
    )

    missing = all_family["NFastaFamilyEntries"].isna()

    if missing.any():
        bad = (
            all_family.loc[
                missing,
                "BaseAccession"
            ]
            .drop_duplicates()
            .tolist()
        )

        die(
            "Mapping contains accession families absent from the frozen FASTA: "
            + ", ".join(map(str, bad[:20]))
        )

    all_family["NFastaFamilyEntries"] = (
        all_family["NFastaFamilyEntries"]
        .astype(int)
    )

    # Validate that each compatible entry is really a member of that FASTA family.
    family_lookup = {
        row.BaseAccession:
            set(
                str(row.FastaFamilyEntries)
                .split(";")
            )
        for row in fasta_family.itertuples()
    }

    invalid_rows = []

    for row in all_family.itertuples():
        mapped = set(
            str(row.CompatibleEntries)
            .split(";")
        )

        allowed = family_lookup.get(
            row.BaseAccession,
            set()
        )

        if not mapped.issubset(allowed):
            invalid_rows.append(
                {
                    "Program": row.Program,
                    "PeptideKey": row.PeptideKey,
                    "BaseAccession": row.BaseAccession,
                    "MappedEntries": ";".join(
                        sorted(mapped)
                    ),
                    "FastaFamilyEntries": ";".join(
                        sorted(allowed)
                    ),
                }
            )

    if invalid_rows:
        pd.DataFrame(
            invalid_rows
        ).to_csv(
            OUT / "ERROR_mapping_entries_not_in_fasta.csv",
            index=False
        )

        die(
            "Some peptide-compatible entries were absent from the frozen FASTA. "
            "See ERROR_mapping_entries_not_in_fasta.csv"
        )

    primary_all = all_family[
        all_family["Category"]
        .isin(
            PRIMARY_CATEGORIES
        )
    ].copy()

    # Frozen anchors must still reproduce exactly before any family filter.
    got_primary = (
        primary_all
        .groupby("Program")["PeptideKey"]
        .nunique()
        .to_dict()
    )

    mismatch = {
        p: (
            got_primary.get(p, 0),
            EXPECTED_PRIMARY[p]
        )
        for p in PROGRAMS
        if got_primary.get(p, 0) != EXPECTED_PRIMARY[p]
    }

    if mismatch:
        die(
            "Primary isoform-discriminative counts do not reproduce "
            f"frozen anchors: {mismatch}"
        )

    got_unique = (
        primary_all[
            primary_all["Category"]
            .eq(
                "single_isoform_unique"
            )
        ]
        .groupby("Program")["PeptideKey"]
        .nunique()
        .to_dict()
    )

    mismatch_unique = {
        p: (
            got_unique.get(p, 0),
            EXPECTED_SINGLE_UNIQUE[p]
        )
        for p in PROGRAMS
        if got_unique.get(p, 0) != EXPECTED_SINGLE_UNIQUE[p]
    }

    if mismatch_unique:
        die(
            "single_isoform_unique counts do not reproduce "
            f"frozen anchors: {mismatch_unique}"
        )

    # Keep only annotated multi-entry families, now defined from FASTA.
    primary = primary_all[
        primary_all["NFastaFamilyEntries"] >= 2
    ].copy()

    # Downstream code uses NFamilyEntries; now this is authoritative.
    primary["NFamilyEntries"] = (
        primary["NFastaFamilyEntries"]
    )

    primary.to_csv(
        OUT / "01_primary_compatibility.csv",
        index=False
    )

    fasta_family.to_csv(
        OUT / "01b_fasta_family_universe.csv",
        index=False
    )

    retention = (
        primary
        .groupby("Program")["PeptideKey"]
        .nunique()
        .reindex(
            PROGRAMS,
            fill_value=0
        )
        .rename(
            "PrimaryPeptidesInMultiEntryFastaFamily"
        )
        .reset_index()
    )

    retention[
        "FrozenPrimaryPeptides"
    ] = retention["Program"].map(
        EXPECTED_PRIMARY
    )

    retention[
        "RetentionFraction"
    ] = (
        retention[
            "PrimaryPeptidesInMultiEntryFastaFamily"
        ]
        /
        retention[
            "FrozenPrimaryPeptides"
        ]
    )

    retention.to_csv(
        OUT / "01c_primary_retention_after_fasta_family_filter.csv",
        index=False
    )

    return primary, fasta_family


def collapse_quant_rows(q: pd.DataFrame) -> pd.DataFrame:
    q = q[q["PeptideKey"].ne("") & q["Sample"].isin(SAMPLES) & np.isfinite(q["Intensity"]) & q["Intensity"].gt(0)].copy()
    return q.groupby(["PeptideKey", "Sample"], as_index=False).agg(Intensity=("Intensity", "median"))


def extract_ap() -> pd.DataFrame:
    rows = []
    for sample in SAMPLES:
        path = PROGRAM_FOLDERS["AP"] / f"{disk_sample(sample)}.ms_data.hdf"
        if not path.exists():
            die("AP file missing:\n" + str(path))
        with h5py.File(path, "r") as h:
            required = ["peptide_fdr/sequence_naked", "peptide_fdr/ms1_int_sum_area"]
            missing = [x for x in required if x not in h]
            if missing:
                die(f"AP HDF missing {missing}:\n{path}")
            seq = decode_hdf_strings(h["peptide_fdr/sequence_naked"][:])
            intensity = np.asarray(h["peptide_fdr/ms1_int_sum_area"][:], dtype=float)
            n = len(seq)
            keep = np.ones(n, dtype=bool)
            if "peptide_fdr/q_value" in h:
                qv = np.asarray(h["peptide_fdr/q_value"][:], dtype=float)
                if len(qv) == n:
                    keep &= np.isfinite(qv) & (qv <= 0.01)
            if "peptide_fdr/target" in h:
                target = np.asarray(h["peptide_fdr/target"][:])
                if len(target) == n:
                    keep &= target.astype(bool)
            if "peptide_fdr/decoy" in h:
                decoy = np.asarray(h["peptide_fdr/decoy"][:])
                if len(decoy) == n:
                    keep &= ~decoy.astype(bool)
            z = pd.DataFrame({"PeptideKey": [il_key(x) for x in seq], "Sample": sample, "Intensity": intensity})
            z = z[keep & z["PeptideKey"].ne("") & np.isfinite(z["Intensity"]) & z["Intensity"].gt(0)].copy()
            rows.append(z)
    return collapse_quant_rows(pd.concat(rows, ignore_index=True))


def fp_peptide_path(sample: str) -> Path:
    exact = PROGRAM_FOLDERS["FP"] / disk_sample(sample) / "peptide.tsv"
    if exact.exists():
        return exact
    token = disk_sample(sample).upper()
    candidates = sorted(set(p for p in PROGRAM_FOLDERS["FP"].rglob("peptide.tsv") if token in str(p.parent).upper()))
    if len(candidates) != 1:
        die(f"FP expected exactly one peptide.tsv for {sample}; found {len(candidates)}:\n" + "\n".join(map(str, candidates)))
    return candidates[0]


def extract_fp() -> pd.DataFrame:
    rows = []
    for sample in SAMPLES:
        path = fp_peptide_path(sample)
        df = pd.read_csv(path, sep="\t", low_memory=False)
        seq_col = find_col(df.columns, ["Peptide", "Sequence", "Base Sequence"])
        int_col = find_col(df.columns, ["Intensity"])
        if seq_col is None or int_col is None:
            die("FP required columns missing:\n" + str(path))
        z = pd.DataFrame({"PeptideKey": df[seq_col].map(il_key), "Sample": sample, "Intensity": pd.to_numeric(df[int_col], errors="coerce")})
        z = z[z["PeptideKey"].ne("") & np.isfinite(z["Intensity"]) & z["Intensity"].gt(0)].copy()
        rows.append(z)
    return collapse_quant_rows(pd.concat(rows, ignore_index=True))


def unique_file(root: Path, name: str) -> Path:
    candidates = sorted(set(root.rglob(name)))
    if len(candidates) != 1:
        die(f"Expected exactly one {name} under:\n{root}\nFound {len(candidates)}:\n" + "\n".join(map(str, candidates)))
    return candidates[0]


def extract_wide_table(path: Path, sequence_candidates: list[str], filter_maxquant_flags: bool = False) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if filter_maxquant_flags:
        reverse_col = find_col(df.columns, ["Reverse"])
        contam_col = find_col(df.columns, ["Potential contaminant", "Potential Contaminant"])
        if reverse_col is not None:
            df = df[~df[reverse_col].astype(str).str.strip().eq("+")].copy()
        if contam_col is not None:
            df = df[~df[contam_col].astype(str).str.strip().eq("+")].copy()
    seq_col = find_col(df.columns, sequence_candidates)
    if seq_col is None:
        die("Sequence column not found:\n" + str(path))
    sample_columns = {}
    for c in df.columns:
        if "intensity" not in str(c).lower():
            continue
        sample = sample_from_text(c)
        if sample not in SAMPLES:
            continue
        score = 2 if str(c).lower().startswith("intensity") else 1
        if sample not in sample_columns or score > sample_columns[sample][0]:
            sample_columns[sample] = (score, str(c))
    missing = [s for s in SAMPLES if s not in sample_columns]
    if missing:
        die(f"Missing intensity columns for {missing}:\n{path}\nResolved columns:\n{sample_columns}")
    peptide_keys = df[seq_col].map(il_key)
    rows = []
    for sample in SAMPLES:
        col = sample_columns[sample][1]
        vals = pd.to_numeric(df[col], errors="coerce")
        z = pd.DataFrame({"PeptideKey": peptide_keys, "Sample": sample, "Intensity": vals})
        z = z[z["PeptideKey"].ne("") & np.isfinite(z["Intensity"]) & z["Intensity"].gt(0)].copy()
        rows.append(z)
    return collapse_quant_rows(pd.concat(rows, ignore_index=True))


def extract_mm() -> pd.DataFrame:
    path = unique_file(PROGRAM_FOLDERS["MM"], "AllQuantifiedPeptides.tsv")
    return extract_wide_table(path, ["Base Sequence", "BaseSequence", "Sequence", "Peptide"], False)


def extract_mq() -> pd.DataFrame:
    exact = PROGRAM_FOLDERS["MQ"] / "txt" / "peptides.txt"
    path = exact if exact.exists() else unique_file(PROGRAM_FOLDERS["MQ"], "peptides.txt")
    return extract_wide_table(path, ["Sequence", "Peptide", "Base Sequence"], True)


def extract_program(program: str) -> pd.DataFrame:
    if program == "AP": return extract_ap()
    if program == "FP": return extract_fp()
    if program == "MM": return extract_mm()
    if program == "MQ": return extract_mq()
    die(f"Unknown program: {program}")


def simplify_constraints(constraints: list[set[str]]) -> tuple[frozenset[str], ...]:
    uniq = {frozenset(x) for x in constraints if len(x) > 0}
    keep = []
    ordered = sorted(uniq, key=lambda x: (len(x), tuple(sorted(x))))
    for s in ordered:
        if not any(k.issubset(s) for k in keep):
            keep.append(s)
    return tuple(keep)


@lru_cache(maxsize=None)
def solve_hitting_set_cached(state: tuple[frozenset[str], ...]) -> frozenset[str]:
    if len(state) == 0:
        return frozenset()
    constraint = min(state, key=lambda x: (len(x), tuple(sorted(x))))
    best = None
    for candidate in sorted(constraint):
        remaining = simplify_constraints([set(s) for s in state if candidate not in s])
        solution = set(solve_hitting_set_cached(remaining)) | {candidate}
        if best is None or len(solution) < len(best):
            best = solution
        if best is not None and len(best) == 1:
            break
    return frozenset() if best is None else frozenset(best)


def solve_inference(constraints: list[set[str]]) -> dict:
    constraints = [set(x) for x in constraints if len(x) > 0]
    if not constraints:
        return {"Kmin": 0, "InferenceClass": "NONE", "IntersectionEntries": "", "NIntersectionEntries": 0, "MinSolution": "", "NConstraints": 0}
    simplified = simplify_constraints(constraints)
    intersection = set(simplified[0])
    for s in simplified[1:]:
        intersection &= set(s)
    if intersection:
        inference = "EXACT_RESOLVED" if len(intersection) == 1 else "SUBSET_RESOLVED"
        return {"Kmin": 1, "InferenceClass": inference, "IntersectionEntries": ";".join(sorted(intersection)), "NIntersectionEntries": len(intersection), "MinSolution": sorted(intersection)[0], "NConstraints": len(simplified)}
    solution = solve_hitting_set_cached(simplified)
    kmin = len(solution)
    if kmin < 2:
        die(f"Internal Kmin inconsistency: empty intersection but Kmin={kmin}")
    return {"Kmin": kmin, "InferenceClass": "MULTIPLE_ISOFORMS_REQUIRED", "IntersectionEntries": "", "NIntersectionEntries": 0, "MinSolution": ";".join(sorted(solution)), "NConstraints": len(simplified)}


def build_run_inference(primary: pd.DataFrame, quant_by_program: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for program in PROGRAMS:
        print(f"\n[{program}] building run-level isoform inference")
        compat = primary[primary["Program"].eq(program)].copy()
        if compat.empty:
            continue
        q = quant_by_program[program].copy()
        q = q[q["PeptideKey"].isin(set(compat["PeptideKey"]))].copy()
        observed = q.merge(compat[["PeptideKey", "Category", "BaseAccession", "CompatibleEntries", "NCompatibleEntries", "NFamilyEntries"]], on="PeptideKey", how="inner", validate="many_to_one")
        for (sample, base), g in observed.groupby(["Sample", "BaseAccession"], sort=False):
            constraints = [set(str(value).split(";")) for value in g["CompatibleEntries"]]
            result = solve_inference(constraints)
            n_unique = int(g.loc[g["Category"].eq("single_isoform_unique"), "PeptideKey"].nunique())
            n_subset = int(g.loc[g["Category"].eq("within_family_subset_discriminative"), "PeptideKey"].nunique())
            all_entries = sorted(set().union(*constraints))
            cell, rep = sample.rsplit("_", 1)
            exact_by_combination = int(result["InferenceClass"] == "EXACT_RESOLVED" and n_unique == 0)
            latent_multi = int(result["InferenceClass"] == "MULTIPLE_ISOFORMS_REQUIRED")
            rows.append({
                "Program": program, "Sample": sample, "CellLine": cell, "Replicate": int(rep), "BaseAccession": base,
                "NFamilyEntries": int(g["NFamilyEntries"].iloc[0]), "NObservedPrimaryPeptides": int(g["PeptideKey"].nunique()),
                "NSingleUniqueObserved": n_unique, "NSubsetObserved": n_subset, "ObservedCandidateEntries": ";".join(all_entries),
                **result, "ExactResolvedByCombination": exact_by_combination, "LatentMultiIsoformEvidence": latent_multi,
            })
    run = pd.DataFrame(rows)
    if run.empty:
        die("No run-level primary isoform evidence found.")
    return run


def complete_run_grid(run: pd.DataFrame, primary: pd.DataFrame) -> pd.DataFrame:
    family_program = primary[["Program", "BaseAccession", "NFamilyEntries"]].drop_duplicates()
    records = []
    for _, r in family_program.iterrows():
        for sample in SAMPLES:
            cell, rep = sample.rsplit("_", 1)
            records.append({"Program": r["Program"], "Sample": sample, "CellLine": cell, "Replicate": int(rep), "BaseAccession": r["BaseAccession"], "NFamilyEntries": int(r["NFamilyEntries"])})
    grid = pd.DataFrame(records).merge(run, on=["Program", "Sample", "CellLine", "Replicate", "BaseAccession", "NFamilyEntries"], how="left", validate="one_to_one")
    for col in ["NObservedPrimaryPeptides", "NSingleUniqueObserved", "NSubsetObserved", "Kmin", "NIntersectionEntries", "NConstraints", "ExactResolvedByCombination", "LatentMultiIsoformEvidence"]:
        if col in grid.columns:
            grid[col] = pd.to_numeric(grid[col], errors="coerce").fillna(0).astype(int)
    for col in ["ObservedCandidateEntries", "IntersectionEntries", "MinSolution"]:
        grid[col] = grid[col].fillna("").astype(str)
    grid["InferenceClass"] = grid["InferenceClass"].fillna("NONE")
    return grid


def mode_with_count(values: list[str]):
    values = [x for x in values if isinstance(x, str) and x != ""]
    if not values:
        return "", 0
    counts = Counter(values)
    value, n = sorted(counts.items(), key=lambda x: (-x[1], x[0]))[0]
    return value, n


def build_cellline_inference(grid: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (program, base, cell), g in grid.groupby(["Program", "BaseAccession", "CellLine"], sort=False):
        g = g.sort_values("Replicate")
        classes = g["InferenceClass"].tolist()
        n_any = int(sum(x != "NONE" for x in classes))
        n_multi = int(sum(x == "MULTIPLE_ISOFORMS_REQUIRED" for x in classes))
        n_exact = int(sum(x == "EXACT_RESOLVED" for x in classes))
        n_subset = int(sum(x == "SUBSET_RESOLVED" for x in classes))
        exact_candidate, exact_candidate_n = mode_with_count(g.loc[g["InferenceClass"].eq("EXACT_RESOLVED"), "IntersectionEntries"].astype(str).tolist())
        subset_candidate, subset_candidate_n = mode_with_count(g.loc[g["InferenceClass"].eq("SUBSET_RESOLVED"), "IntersectionEntries"].astype(str).tolist())
        if n_multi >= 2:
            cell_class = "MULTIPLE_ISOFORMS_REQUIRED"
        elif exact_candidate_n >= 2:
            cell_class = "EXACT_RESOLVED"
        elif subset_candidate_n >= 2:
            cell_class = "SUBSET_RESOLVED"
        elif n_any >= 2:
            cell_class = "VARIABLE_DISCRIMINATIVE"
        elif n_any == 1:
            cell_class = "SPARSE"
        else:
            cell_class = "NONE"
        observed_k = g.loc[g["Kmin"] > 0, "Kmin"].to_numpy(dtype=float)
        rows.append({
            "Program": program, "BaseAccession": base, "CellLine": cell, "NFamilyEntries": int(g["NFamilyEntries"].iloc[0]),
            "CellClass": cell_class, "NRepAny": n_any, "NRepExact": n_exact, "NRepSubset": n_subset, "NRepMulti": n_multi,
            "RepConsistentExactEntry": exact_candidate, "ExactEntryRepCount": exact_candidate_n,
            "RepConsistentSubset": subset_candidate, "SubsetRepCount": subset_candidate_n,
            "MedianKminObserved": float(np.median(observed_k)) if len(observed_k) else 0.0,
            "MaxKminObserved": int(np.max(observed_k)) if len(observed_k) else 0,
            "NRepExactByCombination": int(g["ExactResolvedByCombination"].sum()), "NRepLatentMulti": int(g["LatentMultiIsoformEvidence"].sum()),
            "ReplicatePattern": " | ".join(f"R{int(r.Replicate)}:{r.InferenceClass}(K={int(r.Kmin)})" for r in g.itertuples()),
        })
    return pd.DataFrame(rows)


def build_family_summary(cell: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for base, g in cell.groupby("BaseAccession"):
        counts = g["CellClass"].value_counts().to_dict()
        informative = g[~g["CellClass"].isin(["NONE", "SPARSE"])]
        rows.append({
            "BaseAccession": base,
            "NFamilyEntries": int(g["NFamilyEntries"].max()),
            "NWorkflowCellCombinations": len(g),
            "NAnyEvidence": int((g["CellClass"] != "NONE").sum()),
            "NRepSupportedInformative": int(len(informative)),
            "NExact": int(counts.get("EXACT_RESOLVED", 0)),
            "NSubset": int(counts.get("SUBSET_RESOLVED", 0)),
            "NMulti": int(counts.get("MULTIPLE_ISOFORMS_REQUIRED", 0)),
            "NVariable": int(counts.get("VARIABLE_DISCRIMINATIVE", 0)),
            "NSparse": int(counts.get("SPARSE", 0)),
            "NNone": int(counts.get("NONE", 0)),
            "NExactCombinationReplicates": int(g["NRepExactByCombination"].sum()),
            "NLatentMultiReplicates": int(g["NRepLatentMulti"].sum()),
        })
    out = pd.DataFrame(rows)
    out["PriorityScore"] = out["NMulti"] * 100 + out["NExact"] * 30 + out["NSubset"] * 20 + out["NVariable"] * 10 + out["NExactCombinationReplicates"] * 3 + out["NLatentMultiReplicates"] * 3 + out["NSparse"]
    return out.sort_values(["PriorityScore", "NRepSupportedInformative", "BaseAccession"], ascending=[False, False, True]).reset_index(drop=True)


def build_cellline_contrasts(
    cell: pd.DataFrame
) -> pd.DataFrame:
    """
    Compare replicate-supported cell-line inference states within workflow.

    v1.1.0 fix:
    exact/subset accession labels are appended ONLY when the candidate itself
    has >=2/3 replicate support. One-replicate candidates no longer contaminate
    a replicate-supported CellClass.
    """
    tmp = cell.copy()

    def supported_description(row) -> str:
        cls = str(
            row["CellClass"]
        )

        if cls == "EXACT_RESOLVED":
            if (
                int(row["ExactEntryRepCount"]) >= 2
                and str(row["RepConsistentExactEntry"]) != ""
            ):
                return (
                    "EXACT_RESOLVED["
                    + str(row["RepConsistentExactEntry"])
                    + "]"
                )

            return "EXACT_RESOLVED"

        if cls == "SUBSET_RESOLVED":
            if (
                int(row["SubsetRepCount"]) >= 2
                and str(row["RepConsistentSubset"]) != ""
            ):
                return (
                    "SUBSET_RESOLVED["
                    + str(row["RepConsistentSubset"])
                    + "]"
                )

            return "SUBSET_RESOLVED"

        return cls

    tmp["Description"] = tmp.apply(
        supported_description,
        axis=1
    )

    wide = (
        tmp[
            [
                "Program",
                "BaseAccession",
                "CellLine",
                "Description",
            ]
        ]
        .pivot_table(
            index=[
                "Program",
                "BaseAccession",
            ],
            columns="CellLine",
            values="Description",
            aggfunc="first",
            fill_value="NONE"
        )
        .reset_index()
    )

    for cell_name in CELL_LINES:
        if cell_name not in wide.columns:
            wide[cell_name] = "NONE"

    def contrast_flag(row):
        vals = [
            str(row[x])
            for x in CELL_LINES
        ]

        informative = [
            x
            for x in vals
            if x not in {
                "NONE",
                "SPARSE",
            }
        ]

        return int(
            len(informative) >= 2
            and len(set(informative)) >= 2
        )

    wide[
        "DifferentRepSupportedCellLinePattern"
    ] = wide.apply(
        contrast_flag,
        axis=1
    )

    return wide.sort_values(
        [
            "DifferentRepSupportedCellLinePattern",
            "Program",
            "BaseAccession",
        ],
        ascending=[
            False,
            True,
            True,
        ]
    )


def build_figure_matrix(cell: pd.DataFrame, family_summary: pd.DataFrame) -> pd.DataFrame:
    x = cell.merge(family_summary[["BaseAccession", "PriorityScore"]], on="BaseAccession", how="left")
    x["ProgramCell"] = x["Program"] + "_" + x["CellLine"]
    x["Display"] = x["CellClass"]
    x["DisplayCode"] = x["CellClass"].map({"NONE": 0, "SPARSE": 1, "VARIABLE_DISCRIMINATIVE": 2, "SUBSET_RESOLVED": 3, "EXACT_RESOLVED": 4, "MULTIPLE_ISOFORMS_REQUIRED": 5})
    return x.sort_values(["PriorityScore", "BaseAccession", "Program", "CellLine"], ascending=[False, True, True, True])


def print_summary(primary: pd.DataFrame, run: pd.DataFrame, cell: pd.DataFrame, family: pd.DataFrame):
    print("\n" + "=" * 96)
    print("PRIMARY ISOFORM-DISCRIMINATIVE PEPTIDES")
    print("=" * 96)
    print(primary.groupby(["Program", "Category"])["PeptideKey"].nunique().to_string())
    print("\n" + "=" * 96)
    print("RUN-LEVEL INFERENCE")
    print("=" * 96)
    print(run.groupby(["Program", "InferenceClass"]).size().rename("N").to_string())
    print("\n" + "=" * 96)
    print("REPLICATE-SUPPORTED CELL-LINE INFERENCE")
    print("=" * 96)
    c = cell.groupby(["Program", "CellClass"]).size().rename("N").reset_index()
    print(c.to_string(index=False))
    print("\n" + "=" * 96)
    print("KEY COUNTS")
    print("=" * 96)
    print("Families with replicate-supported MULTIPLE_ISOFORMS_REQUIRED:", int(family["NMulti"].gt(0).sum()))
    print("Families with replicate-supported EXACT_RESOLVED:", int(family["NExact"].gt(0).sum()))
    print("Families with replicate-supported SUBSET_RESOLVED:", int(family["NSubset"].gt(0).sum()))
    print("Families with exact inference rescued by peptide combination:", int((family["NExactCombinationReplicates"] > 0).sum()))
    print("\nTop informative families:")
    cols = ["BaseAccession", "NFamilyEntries", "NExact", "NSubset", "NMulti", "NVariable", "NSparse", "NExactCombinationReplicates", "NLatentMultiReplicates"]
    print(family[cols].head(30).to_string(index=False))


def main():
    print("=" * 96)
    print("ISOFORM KMIN INFERENCE v1.1.0")
    print("Relax inference resolution — NOT peptide confidence")
    print("=" * 96)

    raw_mapping, mapping_source = load_frozen_mapping()
    print("\nFrozen mapping:")
    print(mapping_source)

    fasta_path = locate_frozen_fasta(ROOT)

    print("\nFrozen FASTA:")
    print(fasta_path)

    fasta_family = build_fasta_family_universe(
        fasta_path
    )

    primary, _ = build_compatibility(
        raw_mapping,
        fasta_family
    )

    print(
        "\nPrimary within-family discriminative peptides successfully "
        "reconstructed using the full FASTA-derived family universe."
    )

    quant_by_program = {}
    for program in PROGRAMS:
        print(f"\n[{program}] extracting MBR-OFF peptide quantitative observations...")
        q = extract_program(program)
        print(f"  positive peptide-sample observations: {len(q):,}")
        quant_by_program[program] = q

    run_observed = build_run_inference(primary, quant_by_program)
    run_observed.to_csv(OUT / "02_run_inference.csv", index=False)

    grid = complete_run_grid(run_observed, primary)
    grid.to_csv(OUT / "02b_run_grid_complete.csv", index=False)

    cell = build_cellline_inference(grid)
    cell.to_csv(OUT / "03_cellline_inference.csv", index=False)

    cell_counts = cell.groupby(["Program", "CellLine", "CellClass"]).size().rename("N").reset_index()
    cell_counts.to_csv(OUT / "04_cellline_counts.csv", index=False)

    family = build_family_summary(cell)
    family.to_csv(OUT / "05_family_summary.csv", index=False)

    contrasts = build_cellline_contrasts(cell)
    contrasts.to_csv(OUT / "06_cellline_contrasts.csv", index=False)

    figure_matrix = build_figure_matrix(cell, family)
    figure_matrix.to_csv(OUT / "07_figure_matrix.csv", index=False)

    rescued_exact = run_observed[run_observed["ExactResolvedByCombination"].eq(1)].copy()
    rescued_exact.to_csv(OUT / "08_rescued_exact.csv", index=False)

    multi_required = run_observed[run_observed["InferenceClass"].eq("MULTIPLE_ISOFORMS_REQUIRED")].copy()
    multi_required = multi_required.sort_values(["Kmin", "NObservedPrimaryPeptides", "Program", "BaseAccession", "Sample"], ascending=[False, False, True, True, True])
    multi_required.to_csv(OUT / "09_multi_required.csv", index=False)

    robust_multi = cell[cell["CellClass"].eq("MULTIPLE_ISOFORMS_REQUIRED")].copy()
    robust_multi.to_csv(OUT / "10_robust_multi.csv", index=False)

    informative = cell[cell["CellClass"].isin(["EXACT_RESOLVED", "SUBSET_RESOLVED", "MULTIPLE_ISOFORMS_REQUIRED", "VARIABLE_DISCRIMINATIVE"])].copy()
    informative.to_csv(OUT / "11_rep_supported_informative.csv", index=False)

    method_record = {
        "version": VERSION,
        "mapping_source": str(mapping_source),
        "fasta_source": str(fasta_path),
        "family_universe": (
            "complete accession-family membership derived directly from "
            "the frozen FASTA rather than reconstructed from observed peptide mappings"
        ),
        "primary_categories": sorted(PRIMARY_CATEGORIES),
        "peptide_representation": "stripped sequence with I/L equivalence",
        "search_context": "MBR OFF",
        "quantitative_observation": "positive peptide intensity",
        "imputation": False,
        "run_level_exact": "intersection of all observed peptide-compatible family-entry sets contains exactly one entry",
        "run_level_subset": "intersection contains >1 compatible family entries",
        "run_level_multiple_required": "intersection is empty; exact minimum hitting-set size Kmin >=2",
        "cell_line_replicate_support": "same inference supported in >=2/3 replicates",
        "sparse_definition": "evidence in exactly 1/3 replicates",
        "interpretation_warning": "Kmin>=2 indicates peptide evidence cannot be explained by one annotated UniProt family entry; this is not direct proof of intact full-length proteoforms.",
    }
    (OUT / "METHOD_RECORD.json").write_text(json.dumps(method_record, indent=2), encoding="utf-8")

    print_summary(primary, run_observed, cell, family)
    print("\n" + "=" * 96)
    print("DONE")
    print("=" * 96)
    print("Output folder:")
    print(OUT)
    print("\nMost important files:")
    for name in ["03_cellline_inference.csv", "05_family_summary.csv", "06_cellline_contrasts.csv", "08_rescued_exact.csv", "09_multi_required.csv", "10_robust_multi.csv"]:
        print(" ", name)
    print("=" * 96)


if __name__ == "__main__":
    main()
