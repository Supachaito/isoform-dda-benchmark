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

from collections import Counter, defaultdict
import pandas as pd
import numpy as np
import h5py
import re

ROOT = _public_project_root()

MASTER = ROOT / "MASTER_TABLES_FINAL"

STEP1B = (
    MASTER /
    "STEP1B_MAPPING_QC"
)

AP_FOLDER = (
    ROOT /
    "AP_MBR_OFF"
)

MM_FOLDER = (
    ROOT /
    "MM_MBR_OFF"
)

OUT = (
    MASTER /
    "STEP1C_ANOMALY_RESOLUTION"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# HELPERS
# ============================================================

def dec(x):
    if isinstance(x, bytes):
        return x.decode(
            "utf-8",
            errors="replace"
        )
    return str(x)


def il_norm(x):
    return (
        str(x)
        .upper()
        .replace("I", "J")
        .replace("L", "J")
    )


def find_col(cols, names):

    def n(x):
        return re.sub(
            r"[^a-z0-9]+",
            "",
            str(x).lower()
        )

    lookup = {
        n(c): c
        for c in cols
    }

    for x in names:
        if n(x) in lookup:
            return lookup[n(x)]

    return None


def read_tsv(path):

    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        low_memory=False,
        on_bad_lines="skip"
    )


# ============================================================
# 1. AP UNMAPPED — USE db_idx DIRECTLY
# ============================================================

print()
print("=" * 100)
print("1. ALPHAPEPT UNMAPPED: db_idx DIAGNOSTIC")
print("=" * 100)


unmapped_file = (
    STEP1B /
    "03_Unmapped_Peptide_Audit.csv"
)

unmapped = pd.read_csv(
    unmapped_file,
    dtype=str
)


ap_unmapped = set(
    unmapped.loc[
        unmapped["Program"] == "AP",
        "Peptide"
    ]
    .dropna()
    .astype(str)
    .str.upper()
)


print(
    "AP unmapped peptides:",
    len(ap_unmapped)
)


# ------------------------------------------------------------
# Read AlphaPept database peptide sequences
# ------------------------------------------------------------

dbfile = (
    AP_FOLDER /
    "database.hdf"
)


with h5py.File(
    dbfile,
    "r"
) as db:

    db_peptides = [
        dec(x).upper()
        for x in db[
            "peptides/sequences"
        ][:]
    ]


print(
    "AlphaPept database peptide entries:",
    len(db_peptides)
)


# ------------------------------------------------------------
# Examine peptide_fdr row → db_idx → database peptide
# ------------------------------------------------------------

ap_rows = []


for file in sorted(
    AP_FOLDER.glob(
        "*.ms_data.hdf"
    )
):

    sample = file.name.replace(
        ".ms_data.hdf",
        ""
    )

    with h5py.File(
        file,
        "r"
    ) as h5:

        naked = [
            dec(x).upper()
            for x in h5[
                "peptide_fdr/sequence_naked"
            ][:]
        ]

        modified = (
            [
                dec(x)
                for x in h5[
                    "peptide_fdr/sequence"
                ][:]
            ]
            if (
                "peptide_fdr/sequence"
                in h5
            )
            else [""] * len(naked)
        )

        idxs = h5[
            "peptide_fdr/db_idx"
        ][:]

        qvals = (
            h5[
                "peptide_fdr/q_value"
            ][:]
            if (
                "peptide_fdr/q_value"
                in h5
            )
            else [""] * len(naked)
        )


    for i, pep in enumerate(
        naked
    ):

        if pep not in ap_unmapped:
            continue

        db_idx = int(
            idxs[i]
        )

        db_pep = (
            db_peptides[db_idx]
            if (
                0 <= db_idx
                < len(db_peptides)
            )
            else ""
        )


        if pep == db_pep:

            relation = (
                "EXACT_SAME"
            )

        elif (
            il_norm(pep)
            ==
            il_norm(db_pep)
        ):

            relation = (
                "IL_EQUIVALENT"
            )

        elif (
            pep in db_pep
        ):

            relation = (
                "OUTPUT_IS_SUBSTRING_OF_DB"
            )

        elif (
            db_pep in pep
        ):

            relation = (
                "DB_IS_SUBSTRING_OF_OUTPUT"
            )

        elif (
            len(pep)
            ==
            len(db_pep)
        ):

            relation = (
                "SAME_LENGTH_DIFFERENT_SEQUENCE"
            )

        else:

            relation = (
                "DIFFERENT_SEQUENCE_AND_LENGTH"
            )


        # positions that differ
        mismatch_positions = []

        if len(pep) == len(
            db_pep
        ):

            mismatch_positions = [
                f"{j+1}:{a}>{b}"

                for j, (a, b)
                in enumerate(
                    zip(
                        pep,
                        db_pep
                    )
                )

                if a != b
            ]


        ap_rows.append({
            "Sample":
                sample,

            "PeptideFDR_SequenceNaked":
                pep,

            "PeptideFDR_Modified":
                modified[i],

            "DB_Index":
                db_idx,

            "DatabasePeptideSequence":
                db_pep,

            "OutputLength":
                len(pep),

            "DatabasePeptideLength":
                len(db_pep),

            "Relation":
                relation,

            "MismatchPositions":
                ";".join(
                    mismatch_positions
                ),

            "QValue":
                qvals[i]
        })


ap_detail = pd.DataFrame(
    ap_rows
)


if not ap_detail.empty:

    ap_unique = (
        ap_detail
        .drop_duplicates(
            subset=[
                "PeptideFDR_SequenceNaked",
                "DatabasePeptideSequence"
            ]
        )
    )

    ap_summary = (
        ap_unique[
            "Relation"
        ]
        .value_counts()
        .rename_axis(
            "Relation"
        )
        .reset_index(
            name="DistinctPairs"
        )
    )

else:

    ap_unique = pd.DataFrame()
    ap_summary = pd.DataFrame()


print()
print(
    "AP anomalous evidence rows:",
    len(ap_detail)
)

print(
    "Distinct sequence/db-sequence pairs:",
    len(ap_unique)
)

print()

if not ap_summary.empty:

    print(
        ap_summary.to_string(
            index=False
        )
    )


# ============================================================
# 2. METAMORPHEUS INVALID BASE SEQUENCES
# ============================================================

print()
print("=" * 100)
print("2. METAMORPHEUS INVALID BASE-SEQUENCE DIAGNOSTIC")
print("=" * 100)


invalid_file = (
    STEP1B /
    "04_Invalid_Peptide_Audit.csv"
)


invalid = pd.read_csv(
    invalid_file,
    dtype=str
)


mm_invalid = set(
    invalid[
        "Peptide"
    ]
    .dropna()
    .astype(str)
)


print(
    "MM invalid peptide strings:",
    len(mm_invalid)
)


psm_files = list(
    MM_FOLDER.rglob(
        "AllPSMs.psmtsv"
    )
)


if not psm_files:

    raise RuntimeError(
        "MetaMorpheus AllPSMs.psmtsv not found"
    )


mm = read_tsv(
    psm_files[0]
)


base_col = find_col(
    mm.columns,
    ["Base Sequence"]
)

full_col = find_col(
    mm.columns,
    ["Full Sequence"]
)

ambiguity_col = find_col(
    mm.columns,
    ["Ambiguity Level"]
)

accession_col = find_col(
    mm.columns,
    ["Accession"]
)

target_col = find_col(
    mm.columns,
    [
        "Decoy/Contaminant/Target"
    ]
)

q_col = find_col(
    mm.columns,
    ["QValue"]
)

file_col = find_col(
    mm.columns,
    ["File Name"]
)

scan_col = find_col(
    mm.columns,
    ["Scan Number"]
)


qnum = pd.to_numeric(
    mm[q_col],
    errors="coerce"
)


keep = (
    qnum.le(0.01)
    &
    mm[target_col]
    .fillna("")
    .astype(str)
    .str.strip()
    .eq("T")
)


mm_filtered = mm.loc[
    keep
].copy()


mm_bad = mm_filtered[
    mm_filtered[
        base_col
    ]
    .fillna("")
    .astype(str)
    .isin(
        mm_invalid
    )
].copy()


def special_pattern(x):

    x = str(x)

    chars = sorted(
        {
            c
            for c in x
            if not (
                "A" <= c <= "Z"
            )
        }
    )

    return "".join(
        chars
    )


mm_rows = []


for _, r in mm_bad.iterrows():

    base = str(
        r[base_col]
    )

    mm_rows.append({
        "BaseSequence":
            base,

        "FullSequence":
            (
                r[full_col]
                if full_col
                else ""
            ),

        "AmbiguityLevel":
            (
                r[ambiguity_col]
                if ambiguity_col
                else ""
            ),

        "Accession":
            (
                r[accession_col]
                if accession_col
                else ""
            ),

        "QValue":
            (
                r[q_col]
                if q_col
                else ""
            ),

        "File":
            (
                r[file_col]
                if file_col
                else ""
            ),

        "Scan":
            (
                r[scan_col]
                if scan_col
                else ""
            ),

        "SpecialCharacters":
            special_pattern(
                base
            ),

        "ContainsPipe":
            "|" in base,

        "ContainsSemicolon":
            ";" in base,

        "ContainsBracket":
            any(
                x in base
                for x in "[](){}"
            ),

        "ContainsSlash":
            "/" in base
    })


mm_detail = pd.DataFrame(
    mm_rows
)


if not mm_detail.empty:

    mm_unique = (
        mm_detail
        .drop_duplicates(
            subset=[
                "BaseSequence"
            ]
        )
    )


    mm_char_summary = (
        mm_unique[
            "SpecialCharacters"
        ]
        .value_counts(
            dropna=False
        )
        .rename_axis(
            "SpecialCharacters"
        )
        .reset_index(
            name="DistinctBaseSequences"
        )
    )


    if "AmbiguityLevel" in mm_unique:

        mm_ambiguity_summary = (
            mm_unique[
                "AmbiguityLevel"
            ]
            .value_counts(
                dropna=False
            )
            .rename_axis(
                "AmbiguityLevel"
            )
            .reset_index(
                name="DistinctBaseSequences"
            )
        )

    else:

        mm_ambiguity_summary = (
            pd.DataFrame()
        )

else:

    mm_unique = pd.DataFrame()
    mm_char_summary = pd.DataFrame()
    mm_ambiguity_summary = pd.DataFrame()


print()
print(
    "Filtered target 1%-FDR rows:",
    len(
        mm_filtered
    )
)

print(
    "Rows containing invalid Base Sequence:",
    len(
        mm_detail
    )
)

print(
    "Distinct invalid Base Sequences:",
    len(
        mm_unique
    )
)


print()
print("SPECIAL CHARACTER PATTERNS")

if not mm_char_summary.empty:

    print(
        mm_char_summary.to_string(
            index=False
        )
    )


print()
print("AMBIGUITY LEVEL")

if not mm_ambiguity_summary.empty:

    print(
        mm_ambiguity_summary.to_string(
            index=False
        )
    )


# ============================================================
# 3. EXAMPLE OUTPUT
# ============================================================

print()
print("=" * 100)
print("EXAMPLES")
print("=" * 100)


print()
print("AP examples:")

if not ap_unique.empty:

    print(
        ap_unique[
            [
                "PeptideFDR_SequenceNaked",
                "DatabasePeptideSequence",
                "Relation",
                "MismatchPositions"
            ]
        ]
        .head(20)
        .to_string(
            index=False
        )
    )


print()
print("MetaMorpheus examples:")

if not mm_unique.empty:

    print(
        mm_unique[
            [
                "BaseSequence",
                "FullSequence",
                "AmbiguityLevel",
                "SpecialCharacters"
            ]
        ]
        .head(20)
        .to_string(
            index=False
        )
    )


# ============================================================
# EXPORT
# ============================================================

outputs = {
    "01_AP_Unmapped_dbidx_Detail.csv":
        ap_detail,

    "02_AP_Unmapped_dbidx_UniquePairs.csv":
        ap_unique,

    "03_AP_Unmapped_Relation_Summary.csv":
        ap_summary,

    "04_MM_InvalidPSM_Detail.csv":
        mm_detail,

    "05_MM_Invalid_UniqueSequences.csv":
        mm_unique,

    "06_MM_Invalid_CharacterSummary.csv":
        mm_char_summary,

    "07_MM_Invalid_AmbiguitySummary.csv":
        mm_ambiguity_summary
}


for name, df in outputs.items():

    df.to_csv(
        OUT / name,
        index=False,
        encoding="utf-8-sig"
    )


print()
print("=" * 100)
print("OUTPUT")
print("=" * 100)

print(OUT)

print()
print("STEP 1C COMPLETE")

