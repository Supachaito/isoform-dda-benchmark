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

from itertools import combinations
import pandas as pd
import numpy as np
import h5py
import re
import sys
import traceback

# ============================================================
# STEP 4A V5
# FINAL MBR-OFF QUANTITATIVE BENCHMARK
#
# AP = AlphaPept
# FP = FragPipe
# MM = MetaMorpheus
# MQ = MaxQuant
#
# Primary quantitative unit:
# peptide sequence
#
# Primary isoform-specific quantification:
# single_isoform_unique peptides ONLY
#
# subset-discriminative peptides:
# retained as isoform-subset quantitative entities
# ============================================================


# ============================================================
# 1. PATHS
# ============================================================

ROOT = _public_project_root()

MASTER = ROOT / "MASTER_TABLES_FINAL"

STEP1D = (
    MASTER /
    "STEP1D_FINAL_NORMALIZED"
)

MAPPING_FILE = (
    STEP1D /
    "03_FINAL_PeptideMappingDetail.csv"
)

OUT = (
    MASTER /
    "STEP4A_QUANT_BENCHMARK_FINAL"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)

DIRS = {

    "AP":
        ROOT /
        "AP_MBR_OFF",

    "FP":
        ROOT /
        "FP_MBR_OFF_LFQ",

    "MM":
        ROOT /
        "MM_MBR_OFF",

    "MQ":
        ROOT /
        "MQ_MBR_OFF"
}


# ============================================================
# 2. CONSTANTS
# ============================================================

PROGRAMS = [
    "AP",
    "FP",
    "MM",
    "MQ"
]

SAMPLES = [
    "C33A_1",
    "C33A_2",
    "C33A_3",
    "HELA_1",
    "HELA_2",
    "HELA_3",
    "SIHA_1",
    "SIHA_2",
    "SIHA_3"
]

CELL_LINE = {

    "C33A_1": "C33A",
    "C33A_2": "C33A",
    "C33A_3": "C33A",

    "HELA_1": "HeLa",
    "HELA_2": "HeLa",
    "HELA_3": "HeLa",

    "SIHA_1": "SiHa",
    "SIHA_2": "SiHa",
    "SIHA_3": "SiHa"
}


UNIQUE_CLASS = (
    "single_isoform_unique"
)

SUBSET_CLASS = (
    "within_family_subset_discriminative"
)

ISO_CLASSES = {
    UNIQUE_CLASS,
    SUBSET_CLASS
}


# ============================================================
# 3. HELPERS
# ============================================================

def dec(x):

    if isinstance(
        x,
        (
            bytes,
            np.bytes_
        )
    ):

        return x.decode(
            "utf-8",
            errors="replace"
        )

    return str(x)


def norm_token(x):

    return re.sub(
        r"[^A-Z0-9]",
        "",
        str(x).upper()
    )


def sample_from_text(text):

    t = norm_token(
        text
    )

    for sample in SAMPLES:

        if norm_token(
            sample
        ) in t:

            return sample

    return None


def clean_intensity(series):

    x = pd.to_numeric(
        series,
        errors="coerce"
    )

    return x.where(
        x > 0,
        np.nan
    )


def split_accessions(value):

    if pd.isna(value):
        return []

    return [
        x.strip()
        for x in str(value).split(";")
        if x.strip()
    ]


def is_isoform(accession):

    accession = str(
        accession
    ).strip()

    if "-" not in accession:
        return False

    base, suffix = accession.rsplit(
        "-",
        1
    )

    return (
        bool(base)
        and
        suffix.isdigit()
        and
        int(suffix) >= 1
    )


def target_bool(values):

    out = []

    for x in values:

        if isinstance(
            x,
            (
                bool,
                np.bool_
            )
        ):

            out.append(
                bool(x)
            )

        elif isinstance(
            x,
            (
                int,
                float,
                np.integer,
                np.floating
            )
        ):

            out.append(
                float(x) != 0
            )

        else:

            s = (
                dec(x)
                .strip()
                .lower()
            )

            out.append(
                s in {
                    "true",
                    "t",
                    "1",
                    "yes",
                    "target"
                }
            )

    return np.asarray(
        out,
        dtype=bool
    )


def find_first(
    folder,
    patterns
):

    hits = []

    for pattern in patterns:

        hits.extend(
            folder.rglob(
                pattern
            )
        )

    hits = sorted(
        {
            x
            for x in hits
            if x.is_file()
        },
        key=lambda x: (
            len(x.parts),
            str(x)
        )
    )

    if not hits:
        return None

    return hits[0]


def choose_column(
    df,
    candidates
):

    for c in candidates:

        if c in df.columns:
            return c

    return None


# ============================================================
# 4. AP HDF HELPERS
# ============================================================

def hdf_columns(
    h5,
    object_name
):

    if object_name not in h5:
        return []

    obj = h5[
        object_name
    ]

    if isinstance(
        obj,
        h5py.Group
    ):

        return list(
            obj.keys()
        )

    if isinstance(
        obj,
        h5py.Dataset
    ):

        if obj.dtype.names:

            return list(
                obj.dtype.names
            )

    return []


def hdf_read(
    h5,
    object_name,
    column
):

    obj = h5[
        object_name
    ]

    if isinstance(
        obj,
        h5py.Group
    ):

        if column not in obj:

            raise KeyError(
                f"{object_name}/{column}"
            )

        return np.asarray(
            obj[column][:]
        )


    if isinstance(
        obj,
        h5py.Dataset
    ):

        if (
            obj.dtype.names
            and
            column in obj.dtype.names
        ):

            return np.asarray(
                obj[:][column]
            )


    raise KeyError(
        f"{object_name}/{column}"
    )


def choose_ap_intensity_field(
    h5
):

    columns = hdf_columns(
        h5,
        "feature_table"
    )

    priorities = [
        "ms1_int_sum_area",
        "ms1_int_sum_apex",
        "ms1_int_max_area",
        "ms1_int_max_apex",
        "ms1_int_sum",
        "ms1_int_apex",
        "intensity"
    ]

    for field in priorities:

        if field in columns:

            return field

    raise RuntimeError(
        "AP: no supported MS1 intensity field. "
        +
        str(columns)
    )


# ============================================================
# 5. LOAD FINAL COMMON-REFERENCE MAPPING
# ============================================================

print()
print("=" * 130)
print("STEP 4A V5 — FINAL QUANTITATIVE BENCHMARK")
print("PRIMARY BRANCH: MBR OFF")
print("=" * 130)


if not MAPPING_FILE.exists():

    raise FileNotFoundError(
        MAPPING_FILE
    )


mapping = pd.read_csv(
    MAPPING_FILE,
    dtype=str,
    low_memory=False
)


required_mapping = {
    "Program",
    "Peptide",
    "GeneAwareClass",
    "MappedAccessions"
}


missing = (
    required_mapping
    -
    set(mapping.columns)
)


if missing:

    raise RuntimeError(
        "Mapping file missing columns: "
        +
        str(missing)
    )


mapping = mapping[
    [
        "Program",
        "Peptide",
        "GeneAwareClass",
        "MappedAccessions"
    ]
].copy()


mapping["Program"] = (
    mapping["Program"]
    .fillna("")
    .str.upper()
    .str.strip()
)


mapping["Peptide"] = (
    mapping["Peptide"]
    .fillna("")
    .str.upper()
    .str.strip()
)


# ============================================================
# 6. VERIFY PEPTIDE MAPPING CONSISTENCY
# ============================================================

mapping_qc = (
    mapping
    .groupby(
        "Peptide"
    )
    .agg(
        NClass=(
            "GeneAwareClass",
            "nunique"
        ),
        NMapping=(
            "MappedAccessions",
            "nunique"
        )
    )
    .reset_index()
)


mapping_conflict = mapping_qc[
    (mapping_qc["NClass"] > 1)
    |
    (mapping_qc["NMapping"] > 1)
]


mapping_conflict.to_csv(
    OUT /
    "00_MappingConflict_QC.csv",
    index=False,
    encoding="utf-8-sig"
)


if len(mapping_conflict) > 0:

    raise RuntimeError(
        "Conflicting peptide mappings detected: "
        +
        str(len(mapping_conflict))
    )


annotation = (
    mapping[
        [
            "Peptide",
            "GeneAwareClass",
            "MappedAccessions"
        ]
    ]
    .drop_duplicates(
        subset=[
            "Peptide"
        ]
    )
)


final_sets = {

    p:
        set(
            mapping.loc[
                mapping["Program"] == p,
                "Peptide"
            ]
        )

    for p in PROGRAMS
}


print()
print("Frozen final common-reference peptide sets:")

for p in PROGRAMS:

    print(
        f"  {p}: "
        f"{len(final_sets[p]):,}"
    )


# ============================================================
# 7. EXTRACT AP
# ============================================================

def extract_ap():

    frames = []
    audit_rows = []

    files = sorted(
        DIRS["AP"].rglob(
            "*.ms_data.hdf"
        )
    )


    if not files:

        raise RuntimeError(
            "AP: no ms_data.hdf files found"
        )


    for path in files:

        sample = sample_from_text(
            str(path)
        )

        if sample is None:
            continue


        with h5py.File(
            path,
            "r"
        ) as h5:

            pfdr_cols = hdf_columns(
                h5,
                "peptide_fdr"
            )


            required = [
                "sequence_naked",
                "q_value",
                "target",
                "feature_idx"
            ]


            missing_fields = [
                x
                for x in required
                if x not in pfdr_cols
            ]


            if missing_fields:

                raise RuntimeError(
                    f"{path.name}: "
                    f"missing peptide_fdr fields "
                    f"{missing_fields}"
                )


            peptides = np.asarray(
                [
                    dec(x)
                    .strip()
                    .upper()

                    for x in hdf_read(
                        h5,
                        "peptide_fdr",
                        "sequence_naked"
                    )
                ],
                dtype=object
            )


            qvalue = np.asarray(
                hdf_read(
                    h5,
                    "peptide_fdr",
                    "q_value"
                ),
                dtype=float
            )


            target = target_bool(
                hdf_read(
                    h5,
                    "peptide_fdr",
                    "target"
                )
            )


            feature_raw = np.asarray(
                hdf_read(
                    h5,
                    "peptide_fdr",
                    "feature_idx"
                )
            )


            intensity_field = (
                choose_ap_intensity_field(
                    h5
                )
            )


            feature_intensity = np.asarray(
                hdf_read(
                    h5,
                    "feature_table",
                    intensity_field
                ),
                dtype=float
            )


        feature_numeric = pd.to_numeric(
            pd.Series(
                feature_raw
            ),
            errors="coerce"
        ).to_numpy()


        valid_feature = np.isfinite(
            feature_numeric
        )


        feature_idx = np.full(
            len(feature_numeric),
            -1,
            dtype=int
        )


        feature_idx[
            valid_feature
        ] = (
            feature_numeric[
                valid_feature
            ]
            .astype(int)
        )


        valid_feature &= (
            feature_idx >= 0
        )

        valid_feature &= (
            feature_idx
            <
            len(feature_intensity)
        )


        keep = (
            target
            &
            np.isfinite(
                qvalue
            )
            &
            (
                qvalue <= 0.01
            )
            &
            valid_feature
        )


        x = pd.DataFrame({

            "Peptide":
                peptides[
                    keep
                ],

            "FeatureIdx":
                feature_idx[
                    keep
                ]
        })


        x[
            "RawIntensity"
        ] = (
            feature_intensity[
                x[
                    "FeatureIdx"
                ].to_numpy()
            ]
        )


        x[
            "RawIntensity"
        ] = clean_intensity(
            x[
                "RawIntensity"
            ]
        )


        x = x.dropna(
            subset=[
                "RawIntensity"
            ]
        )


        x = x[
            x[
                "Peptide"
            ].isin(
                final_sets[
                    "AP"
                ]
            )
        ]


        # Do not count same MS1 feature twice
        x = x.drop_duplicates(
            subset=[
                "Peptide",
                "FeatureIdx"
            ]
        )


        x = (
            x
            .groupby(
                "Peptide",
                as_index=False
            )[
                "RawIntensity"
            ]
            .sum()
        )


        x["Program"] = "AP"
        x["Sample"] = sample

        x["QuantSource"] = (
            "feature_table/"
            +
            intensity_field
        )


        frames.append(
            x[
                [
                    "Program",
                    "Sample",
                    "Peptide",
                    "RawIntensity",
                    "QuantSource"
                ]
            ]
        )


        audit_rows.append({

            "Program":
                "AP",

            "Sample":
                sample,

            "SourceFile":
                str(path),

            "QuantSource":
                (
                    "feature_table/"
                    +
                    intensity_field
                ),

            "PositiveQuantPeptides":
                x[
                    "Peptide"
                ].nunique()
        })


    if not frames:

        raise RuntimeError(
            "AP: no quantitative data produced"
        )


    return (
        pd.concat(
            frames,
            ignore_index=True
        ),
        pd.DataFrame(
            audit_rows
        )
    )


# ============================================================
# 8. EXTRACT FP
#
# Primary:
# FP_MBR_OFF_LFQ/<sample>/peptide.tsv::Intensity
#
# Fallback:
# ion.tsv top-3 ion intensity sum
# ============================================================

def extract_fp():

    frames = []
    audit_rows = []


    for sample in SAMPLES:

        folder = (
            DIRS["FP"]
            /
            sample
        )


        peptide_file = (
            folder /
            "peptide.tsv"
        )

        ion_file = (
            folder /
            "ion.tsv"
        )


        if not folder.exists():

            raise RuntimeError(
                f"FP: folder missing: {folder}"
            )


        use_native = False
        native_positive = 0
        final_overlap = 0


        # ----------------------------------------------------
        # Native peptide.tsv
        # ----------------------------------------------------

        if peptide_file.exists():

            df = pd.read_csv(
                peptide_file,
                sep="\t",
                low_memory=False
            )


            seq_col = choose_column(
                df,
                [
                    "Peptide",
                    "Peptide Sequence",
                    "Sequence"
                ]
            )


            if (
                seq_col is not None
                and
                "Intensity"
                in df.columns
            ):

                df["PeptideClean"] = (
                    df[seq_col]
                    .astype(str)
                    .str.upper()
                    .str.strip()
                )


                df[
                    "IntensityClean"
                ] = clean_intensity(
                    df[
                        "Intensity"
                    ]
                )


                native_positive = int(
                    df[
                        "IntensityClean"
                    ]
                    .notna()
                    .sum()
                )


                valid = df[
                    df[
                        "IntensityClean"
                    ]
                    .notna()
                ]


                final_overlap = int(
                    valid[
                        "PeptideClean"
                    ]
                    .isin(
                        final_sets[
                            "FP"
                        ]
                    )
                    .sum()
                )


                if final_overlap > 0:

                    q = pd.DataFrame({

                        "Peptide":
                            valid[
                                "PeptideClean"
                            ],

                        "RawIntensity":
                            valid[
                                "IntensityClean"
                            ]
                    })


                    q = (
                        q
                        .groupby(
                            "Peptide",
                            as_index=False
                        )[
                            "RawIntensity"
                        ]
                        .sum()
                    )


                    q = q[
                        q[
                            "Peptide"
                        ].isin(
                            final_sets[
                                "FP"
                            ]
                        )
                    ].copy()


                    if len(q) > 0:

                        use_native = True

                        source = (
                            "peptide.tsv::Intensity"
                        )


        # ----------------------------------------------------
        # Fallback ion.tsv
        # ----------------------------------------------------

        if not use_native:

            if not ion_file.exists():

                raise RuntimeError(
                    f"FP {sample}: "
                    f"neither usable peptide.tsv "
                    f"nor ion.tsv"
                )


            ion = pd.read_csv(
                ion_file,
                sep="\t",
                low_memory=False
            )


            seq_col = choose_column(
                ion,
                [
                    "Peptide Sequence",
                    "Peptide",
                    "Sequence"
                ]
            )


            if seq_col is None:

                raise RuntimeError(
                    f"FP {sample}: "
                    f"ion.tsv sequence column missing"
                )


            if "Intensity" not in ion.columns:

                raise RuntimeError(
                    f"FP {sample}: "
                    f"ion.tsv Intensity missing"
                )


            ion["Peptide"] = (
                ion[seq_col]
                .astype(str)
                .str.upper()
                .str.strip()
            )


            ion[
                "IonIntensity"
            ] = clean_intensity(
                ion[
                    "Intensity"
                ]
            )


            ion = ion.dropna(
                subset=[
                    "IonIntensity"
                ]
            )


            if len(ion) == 0:

                raise RuntimeError(
                    f"FP {sample}: "
                    f"zero positive ion intensity"
                )


            # top-three ions per stripped peptide
            ion_top3 = (
                ion
                .sort_values(
                    [
                        "Peptide",
                        "IonIntensity"
                    ],
                    ascending=[
                        True,
                        False
                    ]
                )
                .groupby(
                    "Peptide",
                    group_keys=False
                )
                .head(3)
            )


            q = (
                ion_top3
                .groupby(
                    "Peptide",
                    as_index=False
                )[
                    "IonIntensity"
                ]
                .sum()
                .rename(
                    columns={
                        "IonIntensity":
                            "RawIntensity"
                    }
                )
            )


            q = q[
                q[
                    "Peptide"
                ].isin(
                    final_sets[
                        "FP"
                    ]
                )
            ].copy()


            if len(q) == 0:

                raise RuntimeError(
                    f"FP {sample}: "
                    f"ion quant has zero overlap "
                    f"with final peptide set"
                )


            source = (
                "ion.tsv::TOP3_ion_sum"
            )


        q["Program"] = "FP"
        q["Sample"] = sample
        q["QuantSource"] = source


        frames.append(
            q[
                [
                    "Program",
                    "Sample",
                    "Peptide",
                    "RawIntensity",
                    "QuantSource"
                ]
            ]
        )


        audit_rows.append({

            "Program":
                "FP",

            "Sample":
                sample,

            "SourceFile":
                (
                    str(peptide_file)
                    if use_native
                    else
                    str(ion_file)
                ),

            "QuantSource":
                source,

            "NativePositiveRows":
                native_positive,

            "NativeFinalOverlapRows":
                final_overlap,

            "PositiveQuantPeptides":
                q[
                    "Peptide"
                ].nunique()
        })


    return (
        pd.concat(
            frames,
            ignore_index=True
        ),
        pd.DataFrame(
            audit_rows
        )
    )


# ============================================================
# 9. EXTRACT MM
# ============================================================

def extract_mm():

    path = find_first(
        DIRS["MM"],
        [
            "AllQuantifiedPeptides.tsv",
            "*QuantifiedPeptides*.tsv"
        ]
    )


    if path is None:

        raise RuntimeError(
            "MM: quantified peptide file not found"
        )


    df = pd.read_csv(
        path,
        sep="\t",
        low_memory=False
    )


    seq_col = choose_column(
        df,
        [
            "Sequence",
            "Base Sequence",
            "BaseSequence",
            "Peptide",
            "Peptide Sequence"
        ]
    )


    if seq_col is None:

        raise RuntimeError(
            "MM: sequence column not found"
        )


    df["Peptide"] = (
        df[seq_col]
        .astype(str)
        .str.upper()
        .str.strip()
    )


    # exclude MM ambiguous pipe sequences
    df = df[
        ~df[
            "Peptide"
        ].str.contains(
            r"\|",
            regex=True,
            na=False
        )
    ]


    df = df[
        df[
            "Peptide"
        ].isin(
            final_sets[
                "MM"
            ]
        )
    ].copy()


    frames = []
    audit_rows = []


    for sample in SAMPLES:

        preferred = (
            "Intensity_"
            +
            sample
        )


        if preferred in df.columns:

            col = preferred

        else:

            candidates = [

                c

                for c in df.columns

                if (
                    norm_token(sample)
                    in
                    norm_token(c)
                )
                and
                (
                    "INTENSITY"
                    in
                    norm_token(c)
                )
            ]


            if not candidates:

                raise RuntimeError(
                    f"MM: quant column not found "
                    f"for {sample}"
                )


            col = candidates[0]


        q = pd.DataFrame({

            "Peptide":
                df[
                    "Peptide"
                ],

            "RawIntensity":
                clean_intensity(
                    df[
                        col
                    ]
                )
        })


        q = q.dropna(
            subset=[
                "RawIntensity"
            ]
        )


        q = (
            q
            .groupby(
                "Peptide",
                as_index=False
            )[
                "RawIntensity"
            ]
            .sum()
        )


        q["Program"] = "MM"
        q["Sample"] = sample
        q["QuantSource"] = col


        frames.append(
            q[
                [
                    "Program",
                    "Sample",
                    "Peptide",
                    "RawIntensity",
                    "QuantSource"
                ]
            ]
        )


        audit_rows.append({

            "Program":
                "MM",

            "Sample":
                sample,

            "SourceFile":
                str(path),

            "QuantSource":
                col,

            "PositiveQuantPeptides":
                q[
                    "Peptide"
                ].nunique()
        })


    return (
        pd.concat(
            frames,
            ignore_index=True
        ),
        pd.DataFrame(
            audit_rows
        )
    )


# ============================================================
# 10. EXTRACT MQ
# ============================================================

def extract_mq():

    path = find_first(
        DIRS["MQ"],
        [
            "peptides.txt"
        ]
    )


    if path is None:

        raise RuntimeError(
            "MQ: peptides.txt not found"
        )


    df = pd.read_csv(
        path,
        sep="\t",
        low_memory=False
    )


    if "Sequence" not in df.columns:

        raise RuntimeError(
            "MQ: Sequence column missing"
        )


    if "Reverse" in df.columns:

        df = df[
            df[
                "Reverse"
            ].fillna("")
            !=
            "+"
        ]


    if (
        "Potential contaminant"
        in df.columns
    ):

        df = df[
            df[
                "Potential contaminant"
            ].fillna("")
            !=
            "+"
        ]


    df["Peptide"] = (
        df[
            "Sequence"
        ]
        .astype(str)
        .str.upper()
        .str.strip()
    )


    df = df[
        df[
            "Peptide"
        ].isin(
            final_sets[
                "MQ"
            ]
        )
    ].copy()


    frames = []
    audit_rows = []


    for sample in SAMPLES:

        preferred = (
            "Intensity "
            +
            sample
        )


        if preferred in df.columns:

            col = preferred

        else:

            candidates = [

                c

                for c in df.columns

                if (
                    norm_token(sample)
                    in
                    norm_token(c)
                )
                and
                (
                    "INTENSITY"
                    in
                    norm_token(c)
                )
                and
                (
                    "LFQ"
                    not in
                    norm_token(c)
                )
            ]


            if not candidates:

                raise RuntimeError(
                    f"MQ: intensity column not found "
                    f"for {sample}"
                )


            col = candidates[0]


        q = pd.DataFrame({

            "Peptide":
                df[
                    "Peptide"
                ],

            "RawIntensity":
                clean_intensity(
                    df[
                        col
                    ]
                )
        })


        q = q.dropna(
            subset=[
                "RawIntensity"
            ]
        )


        q = (
            q
            .groupby(
                "Peptide",
                as_index=False
            )[
                "RawIntensity"
            ]
            .sum()
        )


        q["Program"] = "MQ"
        q["Sample"] = sample
        q["QuantSource"] = col


        frames.append(
            q[
                [
                    "Program",
                    "Sample",
                    "Peptide",
                    "RawIntensity",
                    "QuantSource"
                ]
            ]
        )


        audit_rows.append({

            "Program":
                "MQ",

            "Sample":
                sample,

            "SourceFile":
                str(path),

            "QuantSource":
                col,

            "PositiveQuantPeptides":
                q[
                    "Peptide"
                ].nunique()
        })


    return (
        pd.concat(
            frames,
            ignore_index=True
        ),
        pd.DataFrame(
            audit_rows
        )
    )


# ============================================================
# 11. RUN EXTRACTION
# ============================================================

extractors = {

    "AP":
        extract_ap,

    "FP":
        extract_fp,

    "MM":
        extract_mm,

    "MQ":
        extract_mq
}


quant_frames = []
audit_frames = []
error_rows = []


for program in PROGRAMS:

    print()
    print(
        f"Extracting {program}..."
    )


    try:

        qdf, adf = (
            extractors[
                program
            ]()
        )


        quant_frames.append(
            qdf
        )

        audit_frames.append(
            adf
        )


        print(
            f"  observations = "
            f"{len(qdf):,}"
        )

        print(
            f"  distinct peptides = "
            f"{qdf['Peptide'].nunique():,}"
        )


    except Exception as e:

        error_rows.append({

            "Program":
                program,

            "Error":
                repr(e),

            "Traceback":
                traceback.format_exc()
        })


        print(
            "  ERROR:",
            repr(e)
        )


error_df = pd.DataFrame(
    error_rows
)


error_df.to_csv(
    OUT /
    "01_ExtractionErrors.csv",
    index=False,
    encoding="utf-8-sig"
)


if len(error_rows) > 0:

    print()
    print("=" * 130)
    print("EXTRACTION FAILED")
    print("=" * 130)

    print(
        error_df[
            [
                "Program",
                "Error"
            ]
        ].to_string(
            index=False
        )
    )

    sys.exit(1)


# ============================================================
# 12. COMBINE ALL FOUR PROGRAMS
# ============================================================

quant = pd.concat(
    quant_frames,
    ignore_index=True
)


audit = pd.concat(
    audit_frames,
    ignore_index=True
)


# ============================================================
# 13. STRICT PROGRAM × SAMPLE QC
# ============================================================

sample_qc = (
    quant
    .groupby(
        [
            "Program",
            "Sample"
        ],
        as_index=False
    )
    .agg(
        QuantifiedPeptides=(
            "Peptide",
            "nunique"
        ),

        QuantObservations=(
            "Peptide",
            "size"
        )
    )
)


expected_pairs = {
    (p, s)
    for p in PROGRAMS
    for s in SAMPLES
}


observed_positive = {

    (
        row.Program,
        row.Sample
    )

    for row
    in sample_qc.itertuples()

    if row.QuantifiedPeptides > 0
}


missing_pairs = (
    expected_pairs
    -
    observed_positive
)


sample_qc.to_csv(
    OUT /
    "02_ProgramSample_ExtractionQC.csv",
    index=False,
    encoding="utf-8-sig"
)


audit.to_csv(
    OUT /
    "03_QuantSource_ByProgramSample.csv",
    index=False,
    encoding="utf-8-sig"
)


if missing_pairs:

    print()
    print(
        "Missing or zero-quant program/sample pairs:"
    )

    print(
        sorted(
            missing_pairs
        )
    )

    sys.exit(1)


# ============================================================
# 14. ANNOTATE COMMON-REFERENCE MAPPING
# ============================================================

quant = quant.merge(
    annotation,
    on="Peptide",
    how="left"
)


quant["CellLine"] = (
    quant[
        "Sample"
    ].map(
        CELL_LINE
    )
)


unmapped_quant = quant[
    quant[
        "GeneAwareClass"
    ].isna()
]


unmapped_quant.to_csv(
    OUT /
    "04_UnmappedQuantPeptides_QC.csv",
    index=False,
    encoding="utf-8-sig"
)


if len(unmapped_quant) > 0:

    raise RuntimeError(
        "Quantified peptides missing final mapping: "
        +
        str(len(unmapped_quant))
    )


# ============================================================
# 15. LOG2 TRANSFORM
# ============================================================

quant["Log2Intensity"] = np.log2(
    quant[
        "RawIntensity"
    ]
)


# ============================================================
# 16. WITHIN-PROGRAM SAMPLE MEDIAN NORMALIZATION
# ============================================================

sample_medians = (
    quant
    .groupby(
        [
            "Program",
            "Sample"
        ],
        as_index=False
    )[
        "Log2Intensity"
    ]
    .median()
    .rename(
        columns={
            "Log2Intensity":
                "SampleMedianLog2"
        }
    )
)


quant = quant.merge(
    sample_medians,
    on=[
        "Program",
        "Sample"
    ],
    how="left"
)


quant["NormLog2"] = (
    quant[
        "Log2Intensity"
    ]
    -
    quant[
        "SampleMedianLog2"
    ]
)


quant["NormLinear"] = np.power(
    2.0,
    quant[
        "NormLog2"
    ]
)


# ============================================================
# 17. PEPTIDE-CENTERED ABUNDANCE
# ============================================================

peptide_center = (
    quant
    .groupby(
        [
            "Program",
            "Peptide"
        ],
        as_index=False
    )[
        "NormLog2"
    ]
    .median()
    .rename(
        columns={
            "NormLog2":
                "PeptideMedianWithinProgram"
        }
    )
)


quant = quant.merge(
    peptide_center,
    on=[
        "Program",
        "Peptide"
    ],
    how="left"
)


quant[
    "PeptideCenteredLog2"
] = (
    quant[
        "NormLog2"
    ]
    -
    quant[
        "PeptideMedianWithinProgram"
    ]
)


# ============================================================
# 18. QUANTIFICATION COVERAGE
# ============================================================

coverage_rows = []


for program in PROGRAMS:

    px = quant[
        quant[
            "Program"
        ]
        ==
        program
    ]


    layers = {

        "AllPeptides":
            final_sets[
                program
            ],

        "IsoformDiscriminative":
            set(
                mapping.loc[
                    (
                        mapping[
                            "Program"
                        ]
                        ==
                        program
                    )
                    &
                    (
                        mapping[
                            "GeneAwareClass"
                        ].isin(
                            ISO_CLASSES
                        )
                    ),
                    "Peptide"
                ]
            ),

        "SingleIsoformUnique":
            set(
                mapping.loc[
                    (
                        mapping[
                            "Program"
                        ]
                        ==
                        program
                    )
                    &
                    (
                        mapping[
                            "GeneAwareClass"
                        ]
                        ==
                        UNIQUE_CLASS
                    ),
                    "Peptide"
                ]
            )
    }


    for layer, eligible in layers.items():

        counts = (
            px[
                px[
                    "Peptide"
                ].isin(
                    eligible
                )
            ]
            .groupby(
                "Peptide"
            )[
                "Sample"
            ]
            .nunique()
        )


        coverage_rows.append({

            "Program":
                program,

            "Layer":
                layer,

            "EligiblePeptides":
                len(
                    eligible
                ),

            "QuantifiedAtLeast1":
                int(
                    (
                        counts >= 1
                    ).sum()
                ),

            "QuantifiedAtLeast3of9":
                int(
                    (
                        counts >= 3
                    ).sum()
                ),

            "QuantifiedAtLeast6of9":
                int(
                    (
                        counts >= 6
                    ).sum()
                ),

            "Quantified9of9":
                int(
                    (
                        counts == 9
                    ).sum()
                )
        })


coverage = pd.DataFrame(
    coverage_rows
)


# ============================================================
# 19. REPLICATE CV
#
# CV calculated in normalized linear space
# require >=2 quantified replicates
# ============================================================

cv_rows = []


for (
    program,
    cellline,
    peptide
), x in quant.groupby(
    [
        "Program",
        "CellLine",
        "Peptide"
    ]
):

    values = (
        x[
            "NormLinear"
        ]
        .dropna()
        .to_numpy()
    )


    if len(values) < 2:
        continue


    mean_value = np.mean(
        values
    )


    if (
        not np.isfinite(
            mean_value
        )
        or
        mean_value <= 0
    ):

        continue


    cv = (
        np.std(
            values,
            ddof=1
        )
        /
        mean_value
        *
        100.0
    )


    cls = (
        x[
            "GeneAwareClass"
        ]
        .iloc[0]
    )


    cv_rows.append({

        "Program":
            program,

        "CellLine":
            cellline,

        "Peptide":
            peptide,

        "GeneAwareClass":
            cls,

        "Layer":
            (
                "IsoformDiscriminative"

                if cls
                in ISO_CLASSES

                else
                "OtherPeptide"
            ),

        "NReplicatesQuantified":
            len(values),

        "CV_Percent":
            cv
    })


cv_df = pd.DataFrame(
    cv_rows
)


# ============================================================
# 20. CROSS-SOFTWARE CORRELATION PER SAMPLE
# ============================================================

correlation_rows = []


for sample in SAMPLES:

    sx = quant[
        quant[
            "Sample"
        ]
        ==
        sample
    ]


    for layer in [
        "AllPeptides",
        "IsoformDiscriminative"
    ]:

        if layer == "AllPeptides":

            lx = sx

        else:

            lx = sx[
                sx[
                    "GeneAwareClass"
                ].isin(
                    ISO_CLASSES
                )
            ]


        for p1, p2 in combinations(
            PROGRAMS,
            2
        ):

            a = (
                lx[
                    lx[
                        "Program"
                    ]
                    ==
                    p1
                ][
                    [
                        "Peptide",
                        "NormLog2"
                    ]
                ]
                .rename(
                    columns={
                        "NormLog2":
                            "A"
                    }
                )
            )


            b = (
                lx[
                    lx[
                        "Program"
                    ]
                    ==
                    p2
                ][
                    [
                        "Peptide",
                        "NormLog2"
                    ]
                ]
                .rename(
                    columns={
                        "NormLog2":
                            "B"
                    }
                )
            )


            m = (
                a
                .merge(
                    b,
                    on="Peptide",
                    how="inner"
                )
                .dropna()
            )


            if len(m) >= 3:

                pearson = (
                    m[
                        [
                            "A",
                            "B"
                        ]
                    ]
                    .corr(
                        method="pearson"
                    )
                    .iloc[
                        0,
                        1
                    ]
                )


                spearman = (
                    m[
                        [
                            "A",
                            "B"
                        ]
                    ]
                    .corr(
                        method="spearman"
                    )
                    .iloc[
                        0,
                        1
                    ]
                )

            else:

                pearson = np.nan
                spearman = np.nan


            correlation_rows.append({

                "Sample":
                    sample,

                "CellLine":
                    CELL_LINE[
                        sample
                    ],

                "Layer":
                    layer,

                "Program1":
                    p1,

                "Program2":
                    p2,

                "CommonPeptides":
                    len(m),

                "PearsonR":
                    pearson,

                "SpearmanR":
                    spearman
            })


correlation = pd.DataFrame(
    correlation_rows
)


# ============================================================
# 21. CELL-LINE MEAN ABUNDANCE
# ============================================================

cell_means = (
    quant
    .groupby(
        [
            "Program",
            "Peptide",
            "GeneAwareClass",
            "CellLine"
        ],
        as_index=False
    )
    .agg(
        MeanNormLog2=(
            "NormLog2",
            "mean"
        ),

        NReplicates=(
            "NormLog2",
            "count"
        )
    )
)


mean_wide = (
    cell_means
    .pivot_table(
        index=[
            "Program",
            "Peptide",
            "GeneAwareClass"
        ],
        columns="CellLine",
        values="MeanNormLog2"
    )
    .reset_index()
)


n_wide = (
    cell_means
    .pivot_table(
        index=[
            "Program",
            "Peptide",
            "GeneAwareClass"
        ],
        columns="CellLine",
        values="NReplicates"
    )
    .reset_index()
)


for cl in [
    "C33A",
    "HeLa",
    "SiHa"
]:

    if cl not in mean_wide.columns:

        mean_wide[
            cl
        ] = np.nan


    if cl not in n_wide.columns:

        n_wide[
            cl
        ] = 0


mean_wide = mean_wide.rename(
    columns={

        "C33A":
            "Mean_C33A",

        "HeLa":
            "Mean_HeLa",

        "SiHa":
            "Mean_SiHa"
    }
)


n_wide = n_wide.rename(
    columns={

        "C33A":
            "N_C33A",

        "HeLa":
            "N_HeLa",

        "SiHa":
            "N_SiHa"
    }
)


fold_change = mean_wide.merge(
    n_wide[
        [
            "Program",
            "Peptide",
            "GeneAwareClass",
            "N_C33A",
            "N_HeLa",
            "N_SiHa"
        ]
    ],
    on=[
        "Program",
        "Peptide",
        "GeneAwareClass"
    ],
    how="left"
)


def calc_fc(
    row,
    mean_a,
    mean_b,
    n_a,
    n_b
):

    if (
        row[
            n_a
        ]
        >=
        2
        and
        row[
            n_b
        ]
        >=
        2
        and
        pd.notna(
            row[
                mean_a
            ]
        )
        and
        pd.notna(
            row[
                mean_b
            ]
        )
    ):

        return (
            row[
                mean_a
            ]
            -
            row[
                mean_b
            ]
        )

    return np.nan


fold_change[
    "log2FC_C33A_vs_HeLa"
] = fold_change.apply(
    lambda r:
        calc_fc(
            r,
            "Mean_C33A",
            "Mean_HeLa",
            "N_C33A",
            "N_HeLa"
        ),
    axis=1
)


fold_change[
    "log2FC_SiHa_vs_HeLa"
] = fold_change.apply(
    lambda r:
        calc_fc(
            r,
            "Mean_SiHa",
            "Mean_HeLa",
            "N_SiHa",
            "N_HeLa"
        ),
    axis=1
)


fold_change[
    "log2FC_SiHa_vs_C33A"
] = fold_change.apply(
    lambda r:
        calc_fc(
            r,
            "Mean_SiHa",
            "Mean_C33A",
            "N_SiHa",
            "N_C33A"
        ),
    axis=1
)


# ============================================================
# 22. FOLD-CHANGE CONCORDANCE
# ============================================================

fc_rows = []


contrast_columns = [
    "log2FC_C33A_vs_HeLa",
    "log2FC_SiHa_vs_HeLa",
    "log2FC_SiHa_vs_C33A"
]


for contrast in contrast_columns:

    for layer in [
        "AllPeptides",
        "IsoformDiscriminative"
    ]:

        if layer == "AllPeptides":

            fx = fold_change

        else:

            fx = fold_change[
                fold_change[
                    "GeneAwareClass"
                ].isin(
                    ISO_CLASSES
                )
            ]


        for p1, p2 in combinations(
            PROGRAMS,
            2
        ):

            a = (
                fx[
                    fx[
                        "Program"
                    ]
                    ==
                    p1
                ][
                    [
                        "Peptide",
                        contrast
                    ]
                ]
                .rename(
                    columns={
                        contrast:
                            "A"
                    }
                )
            )


            b = (
                fx[
                    fx[
                        "Program"
                    ]
                    ==
                    p2
                ][
                    [
                        "Peptide",
                        contrast
                    ]
                ]
                .rename(
                    columns={
                        contrast:
                            "B"
                    }
                )
            )


            m = (
                a
                .merge(
                    b,
                    on="Peptide",
                    how="inner"
                )
                .dropna()
            )


            if len(m) >= 3:

                pearson = (
                    m[
                        [
                            "A",
                            "B"
                        ]
                    ]
                    .corr(
                        method="pearson"
                    )
                    .iloc[
                        0,
                        1
                    ]
                )


                spearman = (
                    m[
                        [
                            "A",
                            "B"
                        ]
                    ]
                    .corr(
                        method="spearman"
                    )
                    .iloc[
                        0,
                        1
                    ]
                )


                same_direction = (
                    (
                        np.sign(
                            m["A"]
                        )
                        ==
                        np.sign(
                            m["B"]
                        )
                    )
                    .mean()
                    *
                    100
                )

            else:

                pearson = np.nan
                spearman = np.nan
                same_direction = np.nan


            fc_rows.append({

                "Contrast":
                    contrast,

                "Layer":
                    layer,

                "Program1":
                    p1,

                "Program2":
                    p2,

                "CommonPeptides":
                    len(m),

                "PearsonR":
                    pearson,

                "SpearmanR":
                    spearman,

                "SameDirectionPercent":
                    same_direction
            })


fc_concordance = pd.DataFrame(
    fc_rows
)


# ============================================================
# 23. ISOFORM-SPECIFIC QUANTIFICATION
#
# ONLY single_isoform_unique peptides
# ============================================================

def get_single_isoform(
    value
):

    accessions = sorted(
        {
            x

            for x
            in split_accessions(
                value
            )

            if is_isoform(
                x
            )
        }
    )


    if len(
        accessions
    ) == 1:

        return accessions[0]

    return None


unique_quant = quant[
    quant[
        "GeneAwareClass"
    ]
    ==
    UNIQUE_CLASS
].copy()


unique_quant[
    "IsoformAccession"
] = unique_quant[
    "MappedAccessions"
].apply(
    get_single_isoform
)


unique_quant = unique_quant[
    unique_quant[
        "IsoformAccession"
    ].notna()
].copy()


isoform_quant = (
    unique_quant
    .groupby(
        [
            "Program",
            "Sample",
            "CellLine",
            "IsoformAccession"
        ],
        as_index=False
    )
    .agg(
        IsoformLog2Abundance=(
            "NormLog2",
            "median"
        ),

        UniquePeptideCount=(
            "Peptide",
            "nunique"
        ),

        SupportingPeptides=(
            "Peptide",
            lambda x:
                ";".join(
                    sorted(
                        set(x)
                    )
                )
        )
    )
)


isoform_quant[
    "EvidenceLevel"
] = np.where(
    isoform_quant[
        "UniquePeptideCount"
    ]
    >=
    2,

    ">=2 single-isoform-unique peptides",

    "1 single-isoform-unique peptide"
)


# ============================================================
# 24. ISOFORM-SUBSET QUANTIFICATION
# ============================================================

subset_quant_raw = quant[
    quant[
        "GeneAwareClass"
    ]
    ==
    SUBSET_CLASS
].copy()


subset_quant_raw[
    "IsoformSubsetID"
] = subset_quant_raw[
    "MappedAccessions"
].apply(
    lambda x:
        "|".join(
            sorted(
                set(
                    split_accessions(
                        x
                    )
                )
            )
        )
)


subset_quant = (
    subset_quant_raw
    .groupby(
        [
            "Program",
            "Sample",
            "CellLine",
            "IsoformSubsetID"
        ],
        as_index=False
    )
    .agg(
        SubsetLog2Abundance=(
            "NormLog2",
            "median"
        ),

        DiscriminativePeptideCount=(
            "Peptide",
            "nunique"
        ),

        SupportingPeptides=(
            "Peptide",
            lambda x:
                ";".join(
                    sorted(
                        set(x)
                    )
                )
        )
    )
)


# ============================================================
# 25. ISOFORM SUMMARY
# ============================================================

isoform_summary = (
    isoform_quant
    .groupby(
        "Program",
        as_index=False
    )
    .agg(
        IsoformsQuantified=(
            "IsoformAccession",
            "nunique"
        ),

        SampleIsoformObservations=(
            "IsoformAccession",
            "size"
        ),

        IsoformsWithAtLeast2UniquePeptides=(
            "UniquePeptideCount",
            lambda x:
                int(
                    (
                        x >= 2
                    ).sum()
                )
        )
    )
)


# ============================================================
# 26. CORRELATION SUMMARY
# ============================================================

correlation_summary = (
    correlation
    .groupby(
        [
            "Layer",
            "Program1",
            "Program2"
        ],
        as_index=False
    )
    .agg(
        MedianCommonPeptides=(
            "CommonPeptides",
            "median"
        ),

        MedianPearsonR=(
            "PearsonR",
            "median"
        ),

        MedianSpearmanR=(
            "SpearmanR",
            "median"
        )
    )
)


# ============================================================
# 27. CV SUMMARY
# ============================================================

cv_summary = (
    cv_df
    .groupby(
        [
            "Program",
            "Layer"
        ],
        as_index=False
    )
    .agg(
        N=(
            "CV_Percent",
            "size"
        ),

        MedianCV=(
            "CV_Percent",
            "median"
        ),

        MeanCV=(
            "CV_Percent",
            "mean"
        )
    )
)


# ============================================================
# 28. EXPORT ALL FINAL TABLES
# ============================================================

quant.to_csv(
    OUT /
    "05_PeptideIntensity_Normalized_Long.csv",
    index=False,
    encoding="utf-8-sig"
)


coverage.to_csv(
    OUT /
    "06_QuantificationCoverage.csv",
    index=False,
    encoding="utf-8-sig"
)


cv_df.to_csv(
    OUT /
    "07_ReplicateCV_PeptideLevel.csv",
    index=False,
    encoding="utf-8-sig"
)


cv_summary.to_csv(
    OUT /
    "08_ReplicateCV_Summary.csv",
    index=False,
    encoding="utf-8-sig"
)


correlation.to_csv(
    OUT /
    "09_CrossSoftware_PeptideCorrelation.csv",
    index=False,
    encoding="utf-8-sig"
)


correlation_summary.to_csv(
    OUT /
    "10_CrossSoftware_CorrelationSummary.csv",
    index=False,
    encoding="utf-8-sig"
)


fold_change.to_csv(
    OUT /
    "11_Peptide_CellLineFoldChanges.csv",
    index=False,
    encoding="utf-8-sig"
)


fc_concordance.to_csv(
    OUT /
    "12_CrossSoftware_FoldChangeConcordance.csv",
    index=False,
    encoding="utf-8-sig"
)


isoform_quant.to_csv(
    OUT /
    "13_IsoformSpecificQuant_SingleUniqueOnly.csv",
    index=False,
    encoding="utf-8-sig"
)


subset_quant.to_csv(
    OUT /
    "14_IsoformSubsetQuant.csv",
    index=False,
    encoding="utf-8-sig"
)


isoform_summary.to_csv(
    OUT /
    "15_IsoformQuantification_Summary.csv",
    index=False,
    encoding="utf-8-sig"
)


sample_medians.to_csv(
    OUT /
    "16_Normalization_SampleMedians.csv",
    index=False,
    encoding="utf-8-sig"
)


# ============================================================
# 29. FINAL STRICT QC
# ============================================================

duplicate_qc = (
    quant
    .groupby(
        [
            "Program",
            "Sample",
            "Peptide"
        ]
    )
    .size()
    .reset_index(
        name="N"
    )
)


duplicate_qc = duplicate_qc[
    duplicate_qc[
        "N"
    ]
    >
    1
]


duplicate_qc.to_csv(
    OUT /
    "17_DuplicateProgramSamplePeptide_QC.csv",
    index=False,
    encoding="utf-8-sig"
)


qc_rows = [

    {
        "Check":
            "FourProgramsExtracted",

        "Value":
            quant[
                "Program"
            ].nunique(),

        "Expected":
            4,

        "PASS":
            (
                quant[
                    "Program"
                ].nunique()
                ==
                4
            )
    },

    {
        "Check":
            "ProgramSamplePairs",

        "Value":
            len(
                observed_positive
            ),

        "Expected":
            36,

        "PASS":
            (
                len(
                    observed_positive
                )
                ==
                36
            )
    },

    {
        "Check":
            "UnmappedQuantPeptides",

        "Value":
            len(
                unmapped_quant
            ),

        "Expected":
            0,

        "PASS":
            (
                len(
                    unmapped_quant
                )
                ==
                0
            )
    },

    {
        "Check":
            "DuplicateProgramSamplePeptide",

        "Value":
            len(
                duplicate_qc
            ),

        "Expected":
            0,

        "PASS":
            (
                len(
                    duplicate_qc
                )
                ==
                0
            )
    },

    {
        "Check":
            "MappingConflicts",

        "Value":
            len(
                mapping_conflict
            ),

        "Expected":
            0,

        "PASS":
            (
                len(
                    mapping_conflict
                )
                ==
                0
            )
    }
]


qc = pd.DataFrame(
    qc_rows
)


qc.to_csv(
    OUT /
    "18_FINAL_QC.csv",
    index=False,
    encoding="utf-8-sig"
)


final_status = (
    "PASS"
    if qc[
        "PASS"
    ].all()
    else
    "FAIL"
)


# ============================================================
# 30. WRITE STATUS FILE
# ============================================================

with open(
    OUT /
    "19_FINAL_STATUS.txt",
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
        "STEP 4A V5 FINAL QUANTITATIVE BENCHMARK\n"
    )

    fh.write(
        "Primary branch: MBR OFF\n"
    )

    fh.write(
        f"FINAL STATUS: {final_status}\n"
    )

    fh.write(
        "FP source: FP_MBR_OFF_LFQ\n"
    )

    fh.write(
        "Isoform-specific accession quantification "
        "uses single_isoform_unique peptides only.\n"
    )

    fh.write(
        "Subset-discriminative peptides are quantified "
        "as isoform-subset entities.\n"
    )


# ============================================================
# 31. PRINT RESULTS
# ============================================================

print()
print("=" * 130)
print("PROGRAM × SAMPLE QUANTIFICATION")
print("=" * 130)

print(
    sample_qc.to_string(
        index=False
    )
)


print()
print("=" * 130)
print("QUANTIFICATION COVERAGE")
print("=" * 130)

print(
    coverage.to_string(
        index=False
    )
)


print()
print("=" * 130)
print("REPLICATE CV SUMMARY")
print("=" * 130)

print(
    cv_summary.to_string(
        index=False
    )
)


print()
print("=" * 130)
print("CROSS-SOFTWARE CORRELATION SUMMARY")
print("=" * 130)

print(
    correlation_summary.to_string(
        index=False
    )
)


print()
print("=" * 130)
print("FOLD-CHANGE CONCORDANCE")
print("=" * 130)

print(
    fc_concordance.to_string(
        index=False
    )
)


print()
print("=" * 130)
print("ISOFORM-SPECIFIC QUANTIFICATION")
print("=" * 130)

print(
    isoform_summary.to_string(
        index=False
    )
)


print()
print("=" * 130)
print("FINAL QC")
print("=" * 130)

print(
    qc.to_string(
        index=False
    )
)


print()
print("=" * 130)
print(
    "FINAL STATUS:",
    final_status
)
print("=" * 130)


print()
print("OUTPUT:")
print(OUT)

print()
print("STEP 4A V5 COMPLETE")

