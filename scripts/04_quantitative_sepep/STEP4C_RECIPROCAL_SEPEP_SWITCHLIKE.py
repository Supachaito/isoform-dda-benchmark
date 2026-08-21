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

import pandas as pd
import numpy as np


# ============================================================
# STEP 4C
# Reciprocal within-gene SEPEP regulation
#
# IMPORTANT:
# This identifies "switch-like" or reciprocal isoform-associated
# quantitative events.
#
# It does NOT by itself prove biological isoform switching.
# ============================================================


# ============================================================
# 1. PATHS
# ============================================================

ROOT = _public_project_root()

MASTER = ROOT / "MASTER_TABLES_FINAL"

STEP4B = (
    MASTER /
    "STEP4B_SEPEPQUANT_INSPIRED"
)

INPUT = (
    STEP4B /
    "09_Gene_vs_SEPEP_FoldChange.csv"
)

OUT = (
    MASTER /
    "STEP4C_RECIPROCAL_SWITCHLIKE"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


PROGRAMS = [
    "AP",
    "FP",
    "MM",
    "MQ"
]


CONTRASTS = [
    "C33A_vs_HeLa",
    "SiHa_vs_HeLa",
    "SiHa_vs_C33A"
]


# ============================================================
# 2. LOAD
# ============================================================

print()
print("=" * 120)
print("STEP 4C — RECIPROCAL WITHIN-GENE SEPEP REGULATION")
print("=" * 120)


if not INPUT.exists():

    raise FileNotFoundError(
        INPUT
    )


df = pd.read_csv(
    INPUT,
    low_memory=False
)


required = {
    "Program",
    "Gene",
    "SEPEPID",
    "SEPEPClass",
    "Contrast",
    "SEPEP_log2FC",
    "Gene_log2FC"
}


missing = required - set(
    df.columns
)


if missing:

    raise RuntimeError(
        f"Missing columns: {sorted(missing)}"
    )


# ============================================================
# 3. NORMALIZE
# ============================================================

df["Program"] = (
    df["Program"]
    .astype(str)
    .str.upper()
    .str.strip()
)


df["Gene"] = (
    df["Gene"]
    .astype(str)
    .str.strip()
)


df["SEPEPID"] = (
    df["SEPEPID"]
    .astype(str)
    .str.strip()
)


df["SEPEPClass"] = (
    df["SEPEPClass"]
    .astype(str)
    .str.upper()
    .str.strip()
)


df["Contrast"] = (
    df["Contrast"]
    .astype(str)
    .str.strip()
)


df["SEPEP_log2FC"] = pd.to_numeric(
    df["SEPEP_log2FC"],
    errors="coerce"
)


df["Gene_log2FC"] = pd.to_numeric(
    df["Gene_log2FC"],
    errors="coerce"
)


df = df[
    df["Program"].isin(
        PROGRAMS
    )
].copy()


df = df[
    df["Contrast"].isin(
        CONTRASTS
    )
].copy()


df = df[
    df["SEPEPClass"].isin(
        [
            "C2",
            "C3"
        ]
    )
].copy()


df = df.dropna(
    subset=[
        "Gene",
        "SEPEPID",
        "SEPEP_log2FC"
    ]
)


# ============================================================
# 4. OPTIONAL REPLICATE FILTER
#
# Step4B already required >=2 replicates for FC.
# Re-enforce if columns exist.
# ============================================================

if (
    "N_A" in df.columns
    and
    "N_B" in df.columns
):

    df["N_A"] = pd.to_numeric(
        df["N_A"],
        errors="coerce"
    )

    df["N_B"] = pd.to_numeric(
        df["N_B"],
        errors="coerce"
    )

    df = df[
        (df["N_A"] >= 2)
        &
        (df["N_B"] >= 2)
    ].copy()


# ============================================================
# 5. CHECK DUPLICATE SEPEP FC
# ============================================================

dup = (
    df
    .groupby(
        [
            "Program",
            "Contrast",
            "Gene",
            "SEPEPID"
        ]
    )
    .size()
    .reset_index(
        name="N"
    )
)


dup_bad = dup[
    dup["N"] > 1
]


if len(
    dup_bad
) > 0:

    dup_bad.to_csv(
        OUT /
        "00_Duplicate_SEPEP_FC_QC.csv",
        index=False
    )

    raise RuntimeError(
        f"Duplicate SEPEP fold-change rows: "
        f"{len(dup_bad)}"
    )


# ============================================================
# 6. FIND WITHIN-GENE MAXIMUM AND MINIMUM SEPEP
#
# For every:
# Program × Contrast × Gene
#
# Require >=2 distinct C2/C3 structural units.
# ============================================================

rows = []


for (
    program,
    contrast,
    gene
), x in df.groupby(
    [
        "Program",
        "Contrast",
        "Gene"
    ],
    dropna=False
):

    x = (
        x
        .drop_duplicates(
            subset=[
                "SEPEPID"
            ]
        )
        .copy()
    )


    n_sepeps = x[
        "SEPEPID"
    ].nunique()


    if n_sepeps < 2:

        continue


    up_idx = x[
        "SEPEP_log2FC"
    ].idxmax()


    down_idx = x[
        "SEPEP_log2FC"
    ].idxmin()


    up = x.loc[
        up_idx
    ]


    down = x.loc[
        down_idx
    ]


    max_fc = float(
        up[
            "SEPEP_log2FC"
        ]
    )


    min_fc = float(
        down[
            "SEPEP_log2FC"
        ]
    )


    gene_fc = float(
        np.nanmedian(
            x[
                "Gene_log2FC"
            ].values
        )
    )


    up_class = str(
        up[
            "SEPEPClass"
        ]
    )


    down_class = str(
        down[
            "SEPEPClass"
        ]
    )


    # --------------------------------------------------------
    # Pair class
    # --------------------------------------------------------

    if (
        up_class == "C2"
        and
        down_class == "C2"
    ):

        pair_class = "C2-C2"


    elif (
        "C2" in [
            up_class,
            down_class
        ]
    ):

        pair_class = "C2-C3"


    else:

        pair_class = "C3-C3"


    # --------------------------------------------------------
    # Three candidate thresholds
    #
    # Any reciprocal:
    #   >0 and <0
    #
    # Moderate:
    #   >= +0.5 and <= -0.5
    #
    # Strong primary:
    #   >= +1 and <= -1
    # --------------------------------------------------------

    reciprocal_any = (
        max_fc > 0
        and
        min_fc < 0
    )


    reciprocal_05 = (
        max_fc >= 0.5
        and
        min_fc <= -0.5
    )


    reciprocal_1 = (
        max_fc >= 1.0
        and
        min_fc <= -1.0
    )


    rows.append({

        "Program":
            program,

        "Contrast":
            contrast,

        "Gene":
            gene,

        "EligibleSEPEPs":
            n_sepeps,

        "MaxSEPEP":
            up[
                "SEPEPID"
            ],

        "MaxClass":
            up_class,

        "MaxSEPEP_log2FC":
            max_fc,

        "MinSEPEP":
            down[
                "SEPEPID"
            ],

        "MinClass":
            down_class,

        "MinSEPEP_log2FC":
            min_fc,

        "PairClass":
            pair_class,

        "Gene_log2FC":
            gene_fc,

        "WithinGeneSpan":
            max_fc
            -
            min_fc,

        "ReciprocalAny":
            reciprocal_any,

        "Reciprocal05":
            reciprocal_05,

        "StrongReciprocal":
            reciprocal_1
    })


eligible = pd.DataFrame(
    rows
)


if len(
    eligible
) == 0:

    raise RuntimeError(
        "No genes with >=2 C2/C3 SEPEPs."
    )


# ============================================================
# 7. PRIMARY STRONG SWITCH-LIKE CANDIDATES
# ============================================================

strong = eligible[
    eligible[
        "StrongReciprocal"
    ]
].copy()


strong[
    "EventID"
] = (
    strong[
        "Gene"
    ]
    +
    " | "
    +
    strong[
        "Contrast"
    ]
)


strong[
    "PairKey"
] = (
    strong[
        "MaxSEPEP"
    ]
    +
    "__UP__"
    +
    strong[
        "MinSEPEP"
    ]
    +
    "__DOWN"
)


# ============================================================
# 8. SUMMARY BY SOFTWARE
# ============================================================

summary_rows = []


for program in PROGRAMS:

    e = eligible[
        eligible[
            "Program"
        ]
        ==
        program
    ]


    s = strong[
        strong[
            "Program"
        ]
        ==
        program
    ]


    row = {

        "Program":
            program,

        "EligibleGeneContrasts":
            len(e),

        "AnyReciprocal":
            int(
                e[
                    "ReciprocalAny"
                ].sum()
            ),

        "Reciprocal_AbsFC0.5":
            int(
                e[
                    "Reciprocal05"
                ].sum()
            ),

        "StrongReciprocal_AbsFC1":
            len(s),

        "StrongFractionPercent":
            (
                100
                *
                len(s)
                /
                len(e)
                if len(e) > 0
                else np.nan
            ),

        "Strong_C2C2":
            int(
                (
                    s[
                        "PairClass"
                    ]
                    ==
                    "C2-C2"
                ).sum()
            ),

        "Strong_C2C3":
            int(
                (
                    s[
                        "PairClass"
                    ]
                    ==
                    "C2-C3"
                ).sum()
            ),

        "Strong_C3C3":
            int(
                (
                    s[
                        "PairClass"
                    ]
                    ==
                    "C3-C3"
                ).sum()
            ),

        "StrongWithAnyC2":
            int(
                s[
                    "PairClass"
                ]
                .isin(
                    [
                        "C2-C2",
                        "C2-C3"
                    ]
                )
                .sum()
            )
    }


    summary_rows.append(
        row
    )


software_summary = pd.DataFrame(
    summary_rows
)


# ============================================================
# 9. CONSENSUS AT GENE × CONTRAST LEVEL
#
# Same biological event may use different exact SEPEPs
# in different software.
#
# This is the primary cross-software consensus definition.
# ============================================================

consensus = (
    strong
    .groupby(
        [
            "Gene",
            "Contrast"
        ],
        as_index=False
    )
    .agg(

        SoftwareN=(
            "Program",
            "nunique"
        ),

        Programs=(
            "Program",
            lambda x:
                "|".join(
                    sorted(
                        set(x)
                    )
                )
        ),

        PairClasses=(
            "PairClass",
            lambda x:
                "|".join(
                    sorted(
                        set(x)
                    )
                )
        ),

        MedianMaxSEPEP_log2FC=(
            "MaxSEPEP_log2FC",
            "median"
        ),

        MedianMinSEPEP_log2FC=(
            "MinSEPEP_log2FC",
            "median"
        ),

        MedianGene_log2FC=(
            "Gene_log2FC",
            "median"
        ),

        MedianWithinGeneSpan=(
            "WithinGeneSpan",
            "median"
        ),

        C2InvolvedSoftwareN=(
            "PairClass",
            lambda x:
                int(
                    pd.Series(x)
                    .isin(
                        [
                            "C2-C2",
                            "C2-C3"
                        ]
                    )
                    .sum()
                )
        )
    )
)


consensus = consensus.sort_values(
    [
        "SoftwareN",
        "C2InvolvedSoftwareN",
        "MedianWithinGeneSpan"
    ],
    ascending=[
        False,
        False,
        False
    ]
)


consensus[
    "EventID"
] = (
    consensus[
        "Gene"
    ]
    +
    " | "
    +
    consensus[
        "Contrast"
    ]
)


# ============================================================
# 10. CONSENSUS SUPPORT DISTRIBUTION
# ============================================================

support_summary = (
    consensus
    .groupby(
        "SoftwareN",
        as_index=False
    )
    .agg(
        Events=(
            "EventID",
            "nunique"
        )
    )
)


support_summary = (
    pd.DataFrame({
        "SoftwareN": [
            1,
            2,
            3,
            4
        ]
    })
    .merge(
        support_summary,
        on="SoftwareN",
        how="left"
    )
)


support_summary[
    "Events"
] = (
    support_summary[
        "Events"
    ]
    .fillna(0)
    .astype(int)
)


# ============================================================
# 11. STRICT PAIR CONSENSUS
#
# Same exact up-SEPEP + down-SEPEP pair
# across software.
#
# More stringent than gene-level event consensus.
# ============================================================

strict_pair = (
    strong
    .groupby(
        [
            "Gene",
            "Contrast",
            "PairKey",
            "MaxSEPEP",
            "MinSEPEP"
        ],
        as_index=False
    )
    .agg(

        SoftwareN=(
            "Program",
            "nunique"
        ),

        Programs=(
            "Program",
            lambda x:
                "|".join(
                    sorted(
                        set(x)
                    )
                )
        ),

        PairClasses=(
            "PairClass",
            lambda x:
                "|".join(
                    sorted(
                        set(x)
                    )
                )
        ),

        MedianMaxFC=(
            "MaxSEPEP_log2FC",
            "median"
        ),

        MedianMinFC=(
            "MinSEPEP_log2FC",
            "median"
        ),

        MedianSpan=(
            "WithinGeneSpan",
            "median"
        )
    )
)


strict_pair = strict_pair.sort_values(
    [
        "SoftwareN",
        "MedianSpan"
    ],
    ascending=[
        False,
        False
    ]
)


# ============================================================
# 12. SOFTWARE SUPPORT MATRIX
# ============================================================

matrix_binary = (
    strong[
        [
            "EventID",
            "Program"
        ]
    ]
    .drop_duplicates()
    .assign(
        Detected=1
    )
    .pivot(
        index="EventID",
        columns="Program",
        values="Detected"
    )
    .fillna(0)
    .astype(int)
    .reset_index()
)


for p in PROGRAMS:

    if p not in matrix_binary.columns:

        matrix_binary[
            p
        ] = 0


matrix_binary = matrix_binary[
    [
        "EventID",
        "AP",
        "FP",
        "MM",
        "MQ"
    ]
]


# Pair-class matrix

matrix_class = (
    strong[
        [
            "EventID",
            "Program",
            "PairClass"
        ]
    ]
    .drop_duplicates(
        subset=[
            "EventID",
            "Program"
        ]
    )
    .pivot(
        index="EventID",
        columns="Program",
        values="PairClass"
    )
    .reset_index()
)


for p in PROGRAMS:

    if p not in matrix_class.columns:

        matrix_class[
            p
        ] = "Not detected"


matrix_class = matrix_class[
    [
        "EventID",
        "AP",
        "FP",
        "MM",
        "MQ"
    ]
]


matrix_class = matrix_class.fillna(
    "Not detected"
)


# ============================================================
# 13. TOP CONSENSUS EVENTS
#
# Prioritize:
# 4/4 > 3/4 > 2/4 > 1/4
# C2 involvement
# span
# ============================================================

top_consensus = (
    consensus[
        consensus[
            "SoftwareN"
        ]
        >=
        2
    ]
    .copy()
)


top_consensus = top_consensus.sort_values(
    [
        "SoftwareN",
        "C2InvolvedSoftwareN",
        "MedianWithinGeneSpan"
    ],
    ascending=[
        False,
        False,
        False
    ]
)


# ============================================================
# 14. EXPORT
# ============================================================

eligible.to_csv(
    OUT /
    "01_Eligible_MultiSEPEP_GeneContrasts.csv",
    index=False,
    encoding="utf-8-sig"
)


strong.to_csv(
    OUT /
    "02_Strong_Reciprocal_SwitchLike_BySoftware.csv",
    index=False,
    encoding="utf-8-sig"
)


software_summary.to_csv(
    OUT /
    "03_SwitchLike_Summary_BySoftware.csv",
    index=False,
    encoding="utf-8-sig"
)


consensus.to_csv(
    OUT /
    "04_StrongReciprocal_Consensus_GeneContrast.csv",
    index=False,
    encoding="utf-8-sig"
)


support_summary.to_csv(
    OUT /
    "05_Consensus_Support_Distribution.csv",
    index=False,
    encoding="utf-8-sig"
)


strict_pair.to_csv(
    OUT /
    "06_Strict_SEPEP_Pair_Consensus.csv",
    index=False,
    encoding="utf-8-sig"
)


matrix_binary.to_csv(
    OUT /
    "07_Software_Support_Matrix.csv",
    index=False,
    encoding="utf-8-sig"
)


matrix_class.to_csv(
    OUT /
    "08_Software_PairClass_Matrix.csv",
    index=False,
    encoding="utf-8-sig"
)


top_consensus.to_csv(
    OUT /
    "09_Top_Consensus_SwitchLike_Candidates.csv",
    index=False,
    encoding="utf-8-sig"
)


# ============================================================
# 15. QC
# ============================================================

qc = pd.DataFrame({

    "Check": [
        "ProgramsPresent",
        "ContrastsPresent",
        "DuplicateProgramGeneContrast",
        "EligibleGeneContrasts",
        "StrongReciprocalEvents",
        "Consensus2plusSoftware",
        "Consensus3plusSoftware",
        "Consensus4Software"
    ],

    "Value": [
        eligible[
            "Program"
        ].nunique(),

        eligible[
            "Contrast"
        ].nunique(),

        int(
            eligible.duplicated(
                subset=[
                    "Program",
                    "Gene",
                    "Contrast"
                ]
            ).sum()
        ),

        len(
            eligible
        ),

        len(
            strong
        ),

        int(
            (
                consensus[
                    "SoftwareN"
                ]
                >=
                2
            ).sum()
        ),

        int(
            (
                consensus[
                    "SoftwareN"
                ]
                >=
                3
            ).sum()
        ),

        int(
            (
                consensus[
                    "SoftwareN"
                ]
                ==
                4
            ).sum()
        )
    ]
})


qc.to_csv(
    OUT /
    "10_STEP4C_QC.csv",
    index=False,
    encoding="utf-8-sig"
)


status = (
    "PASS"
    if (
        eligible[
            "Program"
        ].nunique()
        ==
        4
        and
        eligible[
            "Contrast"
        ].nunique()
        ==
        3
        and
        eligible.duplicated(
            subset=[
                "Program",
                "Gene",
                "Contrast"
            ]
        ).sum()
        ==
        0
    )
    else
    "FAIL"
)


with open(
    OUT /
    "11_STEP4C_STATUS.txt",
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
        f"STEP4C_STATUS={status}\n"
    )

    fh.write(
        "Primary strong switch-like definition: "
        ">=2 C2/C3 SEPEPs within a gene; "
        "max log2FC >= +1 and min log2FC <= -1.\n"
    )

    fh.write(
        "Interpret as reciprocal SEPEP regulation / "
        "switch-like candidate, not confirmed biological isoform switching.\n"
    )


# ============================================================
# 16. PRINT RESULTS
# ============================================================

print()
print("=" * 120)
print("SWITCH-LIKE SUMMARY BY SOFTWARE")
print("=" * 120)

print(
    software_summary.to_string(
        index=False
    )
)


print()
print("=" * 120)
print("CROSS-SOFTWARE CONSENSUS")
print("=" * 120)

print(
    support_summary.to_string(
        index=False
    )
)


print()
print("=" * 120)
print("TOP CONSENSUS SWITCH-LIKE CANDIDATES")
print("=" * 120)

if len(
    top_consensus
) > 0:

    print(
        top_consensus.head(
            20
        ).to_string(
            index=False
        )
    )

else:

    print(
        "No >=2-software strong reciprocal candidates."
    )


print()
print("=" * 120)
print("STRICT EXACT SEPEP-PAIR CONSENSUS")
print("=" * 120)

print(
    strict_pair.head(
        20
    ).to_string(
        index=False
    )
)


print()
print("=" * 120)
print("FINAL STATUS:", status)
print("=" * 120)

print()
print("OUTPUT:")
print(OUT)

