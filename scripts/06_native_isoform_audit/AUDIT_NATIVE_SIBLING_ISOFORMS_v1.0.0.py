#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
AUDIT_NATIVE_SIBLING_ISOFORMS_v1.0.0.py

Goal
----
Audit the NATIVE outputs of AP / FP / MM / MQ directly, rather than relying
only on the normalized mapping table.

Question
--------
For the SAME UniProt base accession, which explicit numeric-suffix isoforms
(e.g. Q15149-8, Q15149-9) are actually retained in each program's own output?

Programs / project folders
--------------------------
AP : Benchmark_Program/AP_MBR_OFF
FP : Benchmark_Program/FP_MBR_OFF_LFQ
MM : Benchmark_Program/MM_MBR_OFF
MQ : Benchmark_Program/MQ_MBR_OFF

Preferred native sources
------------------------
AP : *.ms_data.hdf (string datasets are inspected; sample comes from filename)
FP : peptide.tsv in per-sample folders
MM : AllQuantifiedPeptides.tsv
MQ : peptides.txt

The script also inspects likely protein-group output tables for explicit suffix
IDs, but keeps source provenance so protein-group and peptide-level sightings
are not conflated.

Outputs
-------
.../SUPPLEMENTARY_FIGURES/NATIVE_SIBLING_ISOFORM_AUDIT_V01/

00_file_inventory.csv
01_native_isoform_occurrences.csv
02_native_program_isoforms.csv
03_sibling_families.csv
04_sibling_program_matrix.csv
05_sibling_sample_matrix.csv

Interpretation
--------------
- NativePresent means the explicit UniProt suffix string is retained in that
  program's own native output.
- Sample-aware presence is reported only when the source file/row/intensity
  columns allow sample assignment.
- This audit does NOT claim full-length proteoform proof.
"""

from __future__ import annotations

import re
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

from collections import defaultdict

import numpy as np
import pandas as pd
import h5py


# ======================================================================
# 1. Project paths
# ======================================================================

SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT_PATH.parent

ROOT = None
for candidate in [SCRIPT_DIR, *SCRIPT_DIR.parents]:
    if candidate.name == "Benchmark_Program":
        ROOT = candidate
        break

if ROOT is None:
    ROOT = _public_project_root()

if not ROOT.exists():
    raise RuntimeError(f"Benchmark_Program not found:\n{ROOT}")

PROGRAM_DIRS = {
    "AP": ROOT / "AP_MBR_OFF",
    "FP": ROOT / "FP_MBR_OFF_LFQ",
    "MM": ROOT / "MM_MBR_OFF",
    "MQ": ROOT / "MQ_MBR_OFF",
}

SUPP_ROOT = (
    ROOT
    / "MANUSCRIPT_REVISION_20260813"
    / "MAIN_FIGURES_REBUILD_20260815_V01"
    / "SUPPLEMENTARY_FIGURES"
)

OUT = SUPP_ROOT / "NATIVE_SIBLING_ISOFORM_AUDIT_V01"
OUT.mkdir(parents=True, exist_ok=True)

PROGRAMS = ["AP", "FP", "MM", "MQ"]
CELL_LINES = ["C33A", "SiHa", "HeLa"]


# ======================================================================
# 2. UniProt accession parsing
# ======================================================================

UNIPROT_CORE = (
    r"(?:"
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]"
    r"|"
    r"[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"
    r")"
)

SUFFIX_RE = re.compile(
    rf"(?<![A-Z0-9])({UNIPROT_CORE}-\d+)(?![A-Z0-9])",
    re.I,
)

WRAPPED_SUFFIX_RE = re.compile(
    rf"(?:sp|tr)\|({UNIPROT_CORE}-\d+)\|",
    re.I,
)


def extract_suffix_accessions(value: object) -> list[str]:
    if value is None:
        return []

    if isinstance(value, float) and np.isnan(value):
        return []

    s = str(value).upper()

    found = set()

    for m in WRAPPED_SUFFIX_RE.finditer(s):
        found.add(m.group(1).upper())

    for m in SUFFIX_RE.finditer(s):
        found.add(m.group(1).upper())

    return sorted(found)


def base_accession(acc: str) -> str:
    return re.sub(r"-\d+$", "", str(acc))


# ======================================================================
# 3. Sample parsing
# ======================================================================

def parse_sample(text: object) -> str | None:
    if text is None:
        return None

    s = str(text)

    # C33A_1 / C33A-1 / C33A.1
    m = re.search(
        r"(?i)(C33A|SIHA|HELA)[_\-\.\s]*([123])(?:\D|$)",
        s,
    )

    if not m:
        return None

    cell = {
        "C33A": "C33A",
        "SIHA": "SiHa",
        "HELA": "HeLa",
    }[m.group(1).upper()]

    return f"{cell}_{m.group(2)}"


def sample_cell(sample: str | None) -> str:
    if not sample:
        return ""

    return sample.rsplit("_", 1)[0]


# ======================================================================
# 4. General utilities
# ======================================================================

def find_col(
    columns,
    candidates: list[str],
    contains: bool = False,
) -> str | None:
    cols = [str(c) for c in columns]
    exact = {c.strip().lower(): c for c in cols}

    for candidate in candidates:
        key = candidate.strip().lower()
        if key in exact:
            return exact[key]

    if contains:
        for c in cols:
            lc = c.strip().lower()
            if any(x.lower() in lc for x in candidates):
                return c

    return None


def safe_read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(
        path,
        sep=sep,
        low_memory=False,
    )


# ======================================================================
# 5. Record collector
# ======================================================================

records: list[dict] = []
inventory: list[dict] = []


def add_record(
    program: str,
    isoform: str,
    source_file: Path,
    source_type: str,
    source_column: str = "",
    sample: str | None = None,
    detail: str = "",
):
    records.append(
        {
            "Program": program,
            "BaseAccession": base_accession(isoform),
            "ExactIsoform": isoform,
            "Sample": sample or "",
            "CellLine": sample_cell(sample),
            "SourceType": source_type,
            "SourceFile": str(source_file),
            "SourceColumn": source_column,
            "Detail": detail,
        }
    )


# ======================================================================
# 6. AlphaPept native HDF audit
# ======================================================================

def audit_ap():
    root = PROGRAM_DIRS["AP"]

    if not root.exists():
        raise RuntimeError(f"AP folder not found:\n{root}")

    hdfs = sorted(root.rglob("*.ms_data.hdf"))

    if not hdfs:
        raise RuntimeError(f"No AlphaPept .ms_data.hdf found under:\n{root}")

    for path in hdfs:
        sample = parse_sample(path.name) or parse_sample(str(path.parent))
        datasets_scanned = 0
        suffix_hits = 0

        try:
            with h5py.File(path, "r") as h:

                def visitor(name, obj):
                    nonlocal datasets_scanned, suffix_hits

                    if not isinstance(obj, h5py.Dataset):
                        return

                    # String/object/fixed-string datasets only.
                    kind = obj.dtype.kind
                    if kind not in {"O", "S", "U"}:
                        return

                    datasets_scanned += 1

                    try:
                        arr = obj[()]
                    except Exception:
                        return

                    vals = np.asarray(arr).ravel()

                    # Avoid pathological gigantic string datasets.
                    if vals.size > 5_000_000:
                        return

                    for value in vals:
                        if isinstance(value, (bytes, np.bytes_)):
                            s = value.decode("utf-8", errors="replace")
                        else:
                            s = str(value)

                        accessions = extract_suffix_accessions(s)

                        for acc in accessions:
                            suffix_hits += 1
                            add_record(
                                "AP",
                                acc,
                                path,
                                "native_hdf_string",
                                source_column=name,
                                sample=sample,
                            )

                h.visititems(visitor)

        except Exception as e:
            inventory.append(
                {
                    "Program": "AP",
                    "File": str(path),
                    "PreferredSource": True,
                    "Status": f"ERROR: {e}",
                    "StringDatasetsScanned": 0,
                    "SuffixHits": 0,
                }
            )
            continue

        inventory.append(
            {
                "Program": "AP",
                "File": str(path),
                "PreferredSource": True,
                "Status": "OK",
                "StringDatasetsScanned": datasets_scanned,
                "SuffixHits": suffix_hits,
            }
        )


# ======================================================================
# 7. Generic text-table audit
# ======================================================================

PROTEINISH_NAMES = [
    "protein",
    "proteins",
    "accession",
    "accessions",
    "protein group",
    "protein groups",
    "protein accession",
    "protein accessions",
    "protein ids",
    "protein id",
    "leading razor protein",
    "leading proteins",
    "fasta headers",
]

SAMPLE_COL_NAMES = [
    "raw file",
    "rawfile",
    "experiment",
    "file name",
    "filename",
    "spectra file",
    "run",
    "sample",
]


def likely_accession_columns(df: pd.DataFrame) -> list[str]:
    cols = []

    # Name-based selection.
    for c in df.columns:
        lc = str(c).lower()
        if any(token in lc for token in PROTEINISH_NAMES):
            cols.append(str(c))

    # Content-based rescue from first 200 rows.
    preview = df.head(200)

    for c in df.columns:
        c = str(c)

        if c in cols:
            continue

        if preview[c].dtype.kind not in {"O", "U", "S"}:
            continue

        joined = " ".join(
            preview[c]
            .dropna()
            .astype(str)
            .head(50)
            .tolist()
        )

        if extract_suffix_accessions(joined):
            cols.append(c)

    return sorted(set(cols))


def identify_sample_column(df: pd.DataFrame) -> str | None:
    return find_col(
        df.columns,
        SAMPLE_COL_NAMES,
        contains=False,
    )


def intensity_sample_columns(df: pd.DataFrame) -> dict[str, str]:
    """
    Return mapping {column_name: sample} for wide positive-intensity columns.
    """
    out = {}

    for c in df.columns:
        lc = str(c).lower()

        if not any(
            token in lc
            for token in [
                "intensity",
                "abundance",
                "quantity",
            ]
        ):
            continue

        sample = parse_sample(c)

        if sample:
            out[str(c)] = sample

    return out


def audit_text_file(
    program: str,
    path: Path,
    source_type: str,
    preferred: bool,
):
    try:
        df = safe_read_table(path)
    except Exception as e:
        inventory.append(
            {
                "Program": program,
                "File": str(path),
                "PreferredSource": preferred,
                "Status": f"ERROR: {e}",
                "Rows": 0,
                "AccessionColumns": "",
                "SuffixHits": 0,
            }
        )
        return

    acc_cols = likely_accession_columns(df)

    if not acc_cols:
        inventory.append(
            {
                "Program": program,
                "File": str(path),
                "PreferredSource": preferred,
                "Status": "NO_ACCESSION_COLUMN",
                "Rows": len(df),
                "AccessionColumns": "",
                "SuffixHits": 0,
            }
        )
        return

    sample_col = identify_sample_column(df)
    wide_intensity = intensity_sample_columns(df)
    path_sample = parse_sample(str(path))

    suffix_hits = 0

    # Work row-wise only over accession columns.
    for idx, row in df.iterrows():
        accs = set()

        for c in acc_cols:
            accs.update(
                extract_suffix_accessions(
                    row[c]
                )
            )

        if not accs:
            continue

        row_sample = None

        if sample_col is not None:
            row_sample = parse_sample(
                row[sample_col]
            )

        if row_sample is None:
            row_sample = path_sample

        if row_sample is not None:
            # Sample-aware long/per-sample table.
            for acc in sorted(accs):
                suffix_hits += 1
                add_record(
                    program,
                    acc,
                    path,
                    source_type,
                    source_column=";".join(acc_cols),
                    sample=row_sample,
                    detail=(
                        f"row={idx}"
                    ),
                )

        elif wide_intensity:
            # Wide quantified table: attach accession to samples with positive signal.
            positive_samples = []

            for c, sample in wide_intensity.items():
                val = pd.to_numeric(
                    pd.Series([row[c]]),
                    errors="coerce",
                ).iloc[0]

                if pd.notna(val) and float(val) > 0:
                    positive_samples.append(sample)

            if positive_samples:
                for acc in sorted(accs):
                    for sample in sorted(set(positive_samples)):
                        suffix_hits += 1
                        add_record(
                            program,
                            acc,
                            path,
                            source_type,
                            source_column=";".join(acc_cols),
                            sample=sample,
                            detail=f"row={idx};wide_positive",
                        )
            else:
                # Native accession retained, but no sample-aware quantitative value.
                for acc in sorted(accs):
                    suffix_hits += 1
                    add_record(
                        program,
                        acc,
                        path,
                        source_type,
                        source_column=";".join(acc_cols),
                        sample=None,
                        detail=f"row={idx};no_positive_wide_value",
                    )

        else:
            # Identification-level only.
            for acc in sorted(accs):
                suffix_hits += 1
                add_record(
                    program,
                    acc,
                    path,
                    source_type,
                    source_column=";".join(acc_cols),
                    sample=None,
                    detail=f"row={idx}",
                )

    inventory.append(
        {
            "Program": program,
            "File": str(path),
            "PreferredSource": preferred,
            "Status": "OK",
            "Rows": len(df),
            "AccessionColumns": ";".join(acc_cols),
            "SuffixHits": suffix_hits,
        }
    )


# ======================================================================
# 8. FragPipe
# ======================================================================

def audit_fp():
    root = PROGRAM_DIRS["FP"]

    if not root.exists():
        raise RuntimeError(f"FP folder not found:\n{root}")

    # Preferred peptide-level native files.
    peptide_files = sorted(root.rglob("peptide.tsv"))

    for path in peptide_files:
        audit_text_file(
            "FP",
            path,
            "native_peptide",
            preferred=True,
        )

    # Protein-level audit for suffix retention.
    names = {
        "protein.tsv",
        "combined_protein.tsv",
        "protein_group.tsv",
        "protein_groups.tsv",
    }

    extra = sorted(
        {
            p
            for p in root.rglob("*.tsv")
            if p.name.lower() in names
        }
    )

    for path in extra:
        audit_text_file(
            "FP",
            path,
            "native_protein_group",
            preferred=False,
        )


# ======================================================================
# 9. MetaMorpheus
# ======================================================================

def audit_mm():
    root = PROGRAM_DIRS["MM"]

    if not root.exists():
        raise RuntimeError(f"MM folder not found:\n{root}")

    preferred_names = {
        "allquantifiedpeptides.tsv",
    }

    protein_names = {
        "allquantifiedproteingroups.tsv",
        "allproteingroups.tsv",
        "proteingroups.tsv",
    }

    for path in sorted(root.rglob("*.tsv")):
        name = path.name.lower()

        if name in preferred_names:
            audit_text_file(
                "MM",
                path,
                "native_peptide",
                preferred=True,
            )

        elif name in protein_names:
            audit_text_file(
                "MM",
                path,
                "native_protein_group",
                preferred=False,
            )


# ======================================================================
# 10. MaxQuant
# ======================================================================

def audit_mq():
    root = PROGRAM_DIRS["MQ"]

    if not root.exists():
        raise RuntimeError(f"MQ folder not found:\n{root}")

    preferred = sorted(
        set(root.rglob("peptides.txt"))
    )

    for path in preferred:
        audit_text_file(
            "MQ",
            path,
            "native_peptide",
            preferred=True,
        )

    for target_name in [
        "proteinGroups.txt",
        "evidence.txt",
    ]:
        for path in sorted(
            set(root.rglob(target_name))
        ):
            audit_text_file(
                "MQ",
                path,
                "native_protein_group"
                if target_name.lower().startswith("protein")
                else "native_evidence",
                preferred=False,
            )


# ======================================================================
# 11. Run native audit
# ======================================================================

print("=" * 100)
print("NATIVE SIBLING ISOFORM AUDIT")
print("=" * 100)

audit_ap()
audit_fp()
audit_mm()
audit_mq()

inventory_df = pd.DataFrame(inventory)
inventory_df.to_csv(
    OUT / "00_file_inventory.csv",
    index=False,
)

occ = pd.DataFrame(records)

if occ.empty:
    raise RuntimeError(
        "No explicit numeric-suffix UniProt accessions were retained in the "
        "audited native outputs."
    )

occ = occ.drop_duplicates().sort_values(
    [
        "Program",
        "BaseAccession",
        "ExactIsoform",
        "Sample",
        "SourceType",
        "SourceFile",
    ]
)

occ.to_csv(
    OUT / "01_native_isoform_occurrences.csv",
    index=False,
)


# ======================================================================
# 12. Program-level native suffix presence
# ======================================================================

program_isoforms = (
    occ
    .groupby(
        [
            "Program",
            "BaseAccession",
            "ExactIsoform",
        ],
        as_index=False,
    )
    .agg(
        NNativeSourceFiles=(
            "SourceFile",
            "nunique",
        ),
        NNativeSourceTypes=(
            "SourceType",
            "nunique",
        ),
        NSamplesAssigned=(
            "Sample",
            lambda x: len(
                {
                    s
                    for s in x
                    if str(s) != ""
                }
            ),
        ),
        Samples=(
            "Sample",
            lambda x: ";".join(
                sorted(
                    {
                        str(s)
                        for s in x
                        if str(s) != ""
                    }
                )
            ),
        ),
        SourceTypes=(
            "SourceType",
            lambda x: ";".join(
                sorted(set(map(str, x)))
            ),
        ),
    )
)

program_isoforms["NativePresent"] = True

program_isoforms.to_csv(
    OUT / "02_native_program_isoforms.csv",
    index=False,
)


# ======================================================================
# 13. Same-base sibling families from NATIVE outputs
# ======================================================================

family = (
    program_isoforms
    .groupby(
        "BaseAccession",
        as_index=False,
    )
    .agg(
        NExactSuffixIsoforms=(
            "ExactIsoform",
            "nunique",
        ),
        ExactSuffixIsoforms=(
            "ExactIsoform",
            lambda x: ";".join(
                sorted(set(map(str, x)))
            ),
        ),
        NPrograms=(
            "Program",
            "nunique",
        ),
        Programs=(
            "Program",
            lambda x: ";".join(
                sorted(set(map(str, x)))
            ),
        ),
    )
)

siblings = family[
    family["NExactSuffixIsoforms"] >= 2
].copy()

siblings = siblings.sort_values(
    [
        "NExactSuffixIsoforms",
        "NPrograms",
        "BaseAccession",
    ],
    ascending=[
        False,
        False,
        True,
    ],
)

siblings.to_csv(
    OUT / "03_sibling_families.csv",
    index=False,
)

sibling_bases = set(
    siblings["BaseAccession"]
)


# ======================================================================
# 14. Program matrix for same-base sibling suffixes
# ======================================================================

sib_isoforms = (
    program_isoforms[
        program_isoforms["BaseAccession"].isin(
            sibling_bases
        )
    ][
        [
            "BaseAccession",
            "ExactIsoform",
        ]
    ]
    .drop_duplicates()
)

full_program = (
    sib_isoforms
    .assign(_key=1)
    .merge(
        pd.DataFrame(
            {
                "Program": PROGRAMS,
                "_key": 1,
            }
        ),
        on="_key",
    )
    .drop(columns="_key")
)

matrix_long = full_program.merge(
    program_isoforms[
        [
            "Program",
            "BaseAccession",
            "ExactIsoform",
            "NativePresent",
            "NSamplesAssigned",
            "Samples",
            "SourceTypes",
        ]
    ],
    on=[
        "Program",
        "BaseAccession",
        "ExactIsoform",
    ],
    how="left",
)

matrix_long["NativePresent"] = (
    matrix_long["NativePresent"]
    .fillna(False)
    .astype(bool)
)

matrix_long["NSamplesAssigned"] = (
    matrix_long["NSamplesAssigned"]
    .fillna(0)
    .astype(int)
)

matrix_long["Samples"] = (
    matrix_long["Samples"]
    .fillna("")
)

matrix_long["SourceTypes"] = (
    matrix_long["SourceTypes"]
    .fillna("")
)

matrix_long.to_csv(
    OUT / "04_sibling_program_matrix.csv",
    index=False,
)


# ======================================================================
# 15. Sample matrix
# ======================================================================

sample_presence = (
    occ[
        occ["BaseAccession"].isin(
            sibling_bases
        )
        & occ["Sample"].ne("")
    ]
    .groupby(
        [
            "BaseAccession",
            "ExactIsoform",
            "Program",
            "Sample",
            "CellLine",
        ],
        as_index=False,
    )
    .agg(
        NNativeSources=(
            "SourceFile",
            "nunique",
        ),
        SourceTypes=(
            "SourceType",
            lambda x: ";".join(
                sorted(set(map(str, x)))
            ),
        ),
    )
)

sample_presence.to_csv(
    OUT / "05_sibling_sample_matrix.csv",
    index=False,
)


# ======================================================================
# 16. Console summary
# ======================================================================

print("\nNative files audited:")
print(
    inventory_df
    .groupby(
        [
            "Program",
            "Status",
        ]
    )
    .size()
    .rename("NFiles")
    .to_string()
)

print("\nExplicit suffix isoforms retained natively:")
print(
    program_isoforms
    .groupby("Program")["ExactIsoform"]
    .nunique()
    .reindex(PROGRAMS, fill_value=0)
    .to_string()
)

print(
    "\nSame-base families with >=2 different explicit suffix isoforms:"
)

if siblings.empty:
    print("NONE")
else:
    print(
        siblings.to_string(
            index=False
        )
    )

print(
    "\nCross-program sibling suffix matrix:"
)

if not siblings.empty:
    console = (
        matrix_long
        .assign(
            Mark=np.where(
                matrix_long["NativePresent"],
                "YES",
                "-",
            )
        )
        .pivot_table(
            index=[
                "BaseAccession",
                "ExactIsoform",
            ],
            columns="Program",
            values="Mark",
            aggfunc="first",
            fill_value="-",
        )
        .reset_index()
    )

    for p in PROGRAMS:
        if p not in console.columns:
            console[p] = "-"

    print(
        console[
            [
                "BaseAccession",
                "ExactIsoform",
                *PROGRAMS,
            ]
        ].to_string(
            index=False
        )
    )

print("\nOutput folder:")
print(OUT)
print("=" * 100)
