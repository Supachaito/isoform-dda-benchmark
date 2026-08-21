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

from itertools import product
import pandas as pd
import numpy as np


# =====================================================================
# PATHS
# =====================================================================

ROOT = _public_project_root()

MASTER = ROOT / "MASTER_TABLES_FINAL"

STEP2A = (
    MASTER /
    "STEP2A_CELL_LINE_ANALYSIS"
)

SUPPORT_FILE = (
    STEP2A /
    "03_Peptide_ReplicateSupport.csv"
)

ROBUST_SUMMARY_FILE = (
    STEP2A /
    "02_CellLine_Summary_Robust2of3.csv"
)

SPECIFIC_FILE = (
    STEP2A /
    "08_Strict_CellLineSpecific_Peptides.csv"
)

SPECIFIC_ISOFORM_FILE = (
    STEP2A /
    "10_Strict_CellLineSpecific_IsoformEvidence.csv"
)

OUT = (
    MASTER /
    "STEP2B_CROSS_SOFTWARE_CONCORDANCE"
)

FIGDIR = OUT / "FIGURES"

OUT.mkdir(
    parents=True,
    exist_ok=True
)

FIGDIR.mkdir(
    parents=True,
    exist_ok=True
)


# =====================================================================
# CONSTANTS
# =====================================================================

PROGRAMS = [
    "AP",
    "FP",
    "MM",
    "MQ"
]

CELL_LINES = [
    "C33A",
    "HeLa",
    "SiHa"
]

ISO_CLASSES = {
    "single_isoform_unique",
    "within_family_subset_discriminative"
}

UNIQUE_CLASS = (
    "single_isoform_unique"
)


# =====================================================================
# HELPERS
# =====================================================================

def is_isoform(accession):

    if pd.isna(accession):
        return False

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


def split_accessions(value):

    if pd.isna(value):
        return []

    return [
        x.strip()

        for x in str(
            value
        ).split(";")

        if x.strip()
    ]


def combo_name(program_set):

    return "+".join(
        [
            p
            for p in PROGRAMS
            if p in program_set
        ]
    )


def make_membership(
    df,
    item_col,
    layer_name
):

    rows = []

    for cell_line in CELL_LINES:

        x = df[
            df["CellLine"]
            ==
            cell_line
        ]


        for item in sorted(
            set(
                x[
                    item_col
                ].dropna()
            )
        ):

            programs = set(
                x.loc[
                    x[
                        item_col
                    ]
                    ==
                    item,
                    "Program"
                ]
            )


            rows.append({
                "Layer":
                    layer_name,

                "CellLine":
                    cell_line,

                item_col:
                    item,

                "AP":
                    "AP" in programs,

                "FP":
                    "FP" in programs,

                "MM":
                    "MM" in programs,

                "MQ":
                    "MQ" in programs,

                "ProgramCount":
                    len(
                        programs
                    ),

                "Combination":
                    combo_name(
                        programs
                    )
            })


    return pd.DataFrame(
        rows
    )


def intersection_summary(
    membership,
    item_col
):

    if membership.empty:

        return pd.DataFrame(
            columns=[
                "Layer",
                "CellLine",
                "Combination",
                "ProgramCount",
                "ItemCount"
            ]
        )


    return (
        membership
        .groupby(
            [
                "Layer",
                "CellLine",
                "Combination",
                "ProgramCount"
            ],
            dropna=False
        )
        .size()
        .reset_index(
            name="ItemCount"
        )
        .sort_values(
            [
                "Layer",
                "CellLine",
                "ProgramCount",
                "ItemCount"
            ],
            ascending=[
                True,
                True,
                False,
                False
            ]
        )
    )


def consensus_summary(
    membership
):

    if membership.empty:

        return pd.DataFrame()


    rows = []


    for (
        layer,
        cell_line
    ), x in membership.groupby(
        [
            "Layer",
            "CellLine"
        ]
    ):

        total = len(
            x
        )


        for n in [
            1,
            2,
            3,
            4
        ]:

            count = int(
                (
                    x[
                        "ProgramCount"
                    ]
                    ==
                    n
                ).sum()
            )


            rows.append({
                "Layer":
                    layer,

                "CellLine":
                    cell_line,

                "SoftwareCount":
                    n,

                "ItemCount":
                    count,

                "PercentOfUnion":
                    (
                        100
                        *
                        count
                        /
                        total

                        if total
                        else 0
                    )
            })


    return pd.DataFrame(
        rows
    )


def pairwise_concordance(
    df,
    item_col,
    layer_name
):

    rows = []


    for cell_line in CELL_LINES:

        x = df[
            df[
                "CellLine"
            ]
            ==
            cell_line
        ]


        sets = {}

        for program in PROGRAMS:

            sets[
                program
            ] = set(
                x.loc[
                    x[
                        "Program"
                    ]
                    ==
                    program,
                    item_col
                ]
                .dropna()
            )


        for p1 in PROGRAMS:

            for p2 in PROGRAMS:

                s1 = sets[
                    p1
                ]

                s2 = sets[
                    p2
                ]


                intersection = len(
                    s1
                    &
                    s2
                )


                union = len(
                    s1
                    |
                    s2
                )


                minimum = min(
                    len(
                        s1
                    ),
                    len(
                        s2
                    )
                )


                jaccard = (
                    intersection
                    /
                    union

                    if union
                    else np.nan
                )


                overlap_coefficient = (
                    intersection
                    /
                    minimum

                    if minimum
                    else np.nan
                )


                rows.append({
                    "Layer":
                        layer_name,

                    "CellLine":
                        cell_line,

                    "Program1":
                        p1,

                    "Program2":
                        p2,

                    "N_Program1":
                        len(
                            s1
                        ),

                    "N_Program2":
                        len(
                            s2
                        ),

                    "Intersection":
                        intersection,

                    "Union":
                        union,

                    "Jaccard":
                        jaccard,

                    "OverlapCoefficient":
                        overlap_coefficient
                })


    return pd.DataFrame(
        rows
    )


# =====================================================================
# LOAD STEP 2A DATA
# =====================================================================

print()
print("=" * 120)
print("STEP 2B — CROSS-SOFTWARE CONCORDANCE")
print("ROBUST EVIDENCE = DETECTED IN >=2/3 REPLICATES")
print("=" * 120)

print()
print("Input:")
print(SUPPORT_FILE)


support = pd.read_csv(
    SUPPORT_FILE,
    dtype=str,
    low_memory=False
)


required = {
    "Program",
    "CellLine",
    "Peptide",
    "ReplicateSupport",
    "GeneAwareClass",
    "MappedAccessions"
}


missing = (
    required
    -
    set(
        support.columns
    )
)


if missing:

    raise RuntimeError(
        f"Missing columns: {missing}"
    )


for col in [
    "Program",
    "CellLine",
    "Peptide",
    "GeneAwareClass"
]:

    support[
        col
    ] = (
        support[
            col
        ]
        .fillna("")
        .astype(str)
        .str.strip()
    )


support[
    "Program"
] = support[
    "Program"
].str.upper()


support[
    "Peptide"
] = support[
    "Peptide"
].str.upper()


support[
    "ReplicateSupport"
] = pd.to_numeric(
    support[
        "ReplicateSupport"
    ],
    errors="coerce"
)


# =====================================================================
# BASIC QC
# =====================================================================

duplicates = (
    support
    .duplicated(
        subset=[
            "Program",
            "CellLine",
            "Peptide"
        ]
    )
    .sum()
)


bad_rep_support = support[
    ~support[
        "ReplicateSupport"
    ].isin(
        [
            1,
            2,
            3
        ]
    )
]


observed_pairs = (
    support[
        [
            "Program",
            "CellLine"
        ]
    ]
    .drop_duplicates()
)


# =====================================================================
# ROBUST SETS >=2/3
# =====================================================================

robust = support[
    support[
        "ReplicateSupport"
    ]
    >=
    2
].copy()


robust_all = (
    robust[
        [
            "Program",
            "CellLine",
            "Peptide",
            "ReplicateSupport",
            "GeneAwareClass",
            "MappedAccessions"
        ]
    ]
    .drop_duplicates()
)


robust_iso_peptides = robust_all[
    robust_all[
        "GeneAwareClass"
    ].isin(
        ISO_CLASSES
    )
].copy()


robust_unique_peptides = robust_all[
    robust_all[
        "GeneAwareClass"
    ]
    ==
    UNIQUE_CLASS
].copy()


print()
print(
    "Robust peptide records:",
    len(
        robust_all
    )
)

print(
    "Robust isoform-discriminative records:",
    len(
        robust_iso_peptides
    )
)


# =====================================================================
# EXPAND ROBUST ISOFORM ACCESSIONS
# =====================================================================

isoform_rows = []


for _, r in (
    robust_iso_peptides
    .iterrows()
):

    for accession in split_accessions(
        r[
            "MappedAccessions"
        ]
    ):

        if not is_isoform(
            accession
        ):
            continue


        isoform_rows.append({

            "Program":
                r[
                    "Program"
                ],

            "CellLine":
                r[
                    "CellLine"
                ],

            "IsoformAccession":
                accession,

            "Peptide":
                r[
                    "Peptide"
                ],

            "EvidenceClass":
                r[
                    "GeneAwareClass"
                ],

            "ReplicateSupport":
                r[
                    "ReplicateSupport"
                ]
        })


robust_isoform_detail = pd.DataFrame(
    isoform_rows
)


if robust_isoform_detail.empty:

    raise RuntimeError(
        "No robust isoform evidence generated."
    )


robust_implicated_isoforms = (
    robust_isoform_detail[
        [
            "Program",
            "CellLine",
            "IsoformAccession"
        ]
    ]
    .drop_duplicates()
)


# =====================================================================
# UNIQUELY RESOLVED ISOFORMS
# at least one robust single-isoform-unique peptide
# =====================================================================

robust_unique_isoforms = (
    robust_isoform_detail[
        robust_isoform_detail[
            "EvidenceClass"
        ]
        ==
        UNIQUE_CLASS
    ][
        [
            "Program",
            "CellLine",
            "IsoformAccession"
        ]
    ]
    .drop_duplicates()
)


# =====================================================================
# MEMBERSHIP MATRICES
# =====================================================================

membership_all = make_membership(
    robust_all,
    "Peptide",
    "Robust_AllPeptides"
)


membership_iso_peptide = make_membership(
    robust_iso_peptides,
    "Peptide",
    "Robust_IsoformDiscriminativePeptides"
)


membership_unique_isoform = make_membership(
    robust_unique_isoforms,
    "IsoformAccession",
    "Robust_UniquelyResolvedIsoforms"
)


membership_implicated_isoform = make_membership(
    robust_implicated_isoforms,
    "IsoformAccession",
    "Robust_DiscriminativelyImplicatedIsoforms"
)


# =====================================================================
# INTERSECTION SUMMARIES
# =====================================================================

intersection_all = intersection_summary(
    membership_all,
    "Peptide"
)


intersection_iso_peptide = intersection_summary(
    membership_iso_peptide,
    "Peptide"
)


intersection_unique_isoform = intersection_summary(
    membership_unique_isoform,
    "IsoformAccession"
)


intersection_implicated_isoform = intersection_summary(
    membership_implicated_isoform,
    "IsoformAccession"
)


intersection_combined = pd.concat(
    [
        intersection_all,
        intersection_iso_peptide,
        intersection_unique_isoform,
        intersection_implicated_isoform
    ],
    ignore_index=True
)


# =====================================================================
# CONSENSUS 1/4, 2/4, 3/4, 4/4
# =====================================================================

consensus_combined = pd.concat(
    [
        consensus_summary(
            membership_all
        ),

        consensus_summary(
            membership_iso_peptide
        ),

        consensus_summary(
            membership_unique_isoform
        ),

        consensus_summary(
            membership_implicated_isoform
        )
    ],
    ignore_index=True
)


# =====================================================================
# PAIRWISE JACCARD / OVERLAP COEFFICIENT
# =====================================================================

pairwise_all = pairwise_concordance(
    robust_all,
    "Peptide",
    "Robust_AllPeptides"
)


pairwise_iso_peptide = pairwise_concordance(
    robust_iso_peptides,
    "Peptide",
    "Robust_IsoformDiscriminativePeptides"
)


pairwise_unique_isoform = pairwise_concordance(
    robust_unique_isoforms,
    "IsoformAccession",
    "Robust_UniquelyResolvedIsoforms"
)


pairwise_implicated_isoform = pairwise_concordance(
    robust_implicated_isoforms,
    "IsoformAccession",
    "Robust_DiscriminativelyImplicatedIsoforms"
)


pairwise_combined = pd.concat(
    [
        pairwise_all,
        pairwise_iso_peptide,
        pairwise_unique_isoform,
        pairwise_implicated_isoform
    ],
    ignore_index=True
)


# =====================================================================
# ROBUST COUNT SUMMARY PER SOFTWARE × CELL LINE
# =====================================================================

count_rows = []


for program in PROGRAMS:

    for cell_line in CELL_LINES:

        all_x = robust_all[
            (
                robust_all[
                    "Program"
                ]
                ==
                program
            )
            &
            (
                robust_all[
                    "CellLine"
                ]
                ==
                cell_line
            )
        ]


        iso_x = robust_iso_peptides[
            (
                robust_iso_peptides[
                    "Program"
                ]
                ==
                program
            )
            &
            (
                robust_iso_peptides[
                    "CellLine"
                ]
                ==
                cell_line
            )
        ]


        unique_iso_x = robust_unique_isoforms[
            (
                robust_unique_isoforms[
                    "Program"
                ]
                ==
                program
            )
            &
            (
                robust_unique_isoforms[
                    "CellLine"
                ]
                ==
                cell_line
            )
        ]


        implicated_iso_x = robust_implicated_isoforms[
            (
                robust_implicated_isoforms[
                    "Program"
                ]
                ==
                program
            )
            &
            (
                robust_implicated_isoforms[
                    "CellLine"
                ]
                ==
                cell_line
            )
        ]


        count_rows.append({

            "Program":
                program,

            "CellLine":
                cell_line,

            "RobustAllPeptides":
                all_x[
                    "Peptide"
                ].nunique(),

            "RobustIsoformDiscriminativePeptides":
                iso_x[
                    "Peptide"
                ].nunique(),

            "RobustUniquelyResolvedIsoforms":
                unique_iso_x[
                    "IsoformAccession"
                ].nunique(),

            "RobustDiscriminativelyImplicatedIsoforms":
                implicated_iso_x[
                    "IsoformAccession"
                ].nunique()
        })


robust_counts = pd.DataFrame(
    count_rows
)


# =====================================================================
# CHECK AGAINST STEP 2A ROBUST SUMMARY
# =====================================================================

step2a_summary = pd.read_csv(
    ROBUST_SUMMARY_FILE
)


check = robust_counts.merge(
    step2a_summary[
        [
            "Program",
            "CellLine",
            "DistinctPeptides",
            "TotalIsoformDiscriminative"
        ]
    ],
    on=[
        "Program",
        "CellLine"
    ],
    how="left"
)


check[
    "AllPeptideDifference"
] = (
    check[
        "RobustAllPeptides"
    ]
    -
    check[
        "DistinctPeptides"
    ]
)


check[
    "IsoformPeptideDifference"
] = (
    check[
        "RobustIsoformDiscriminativePeptides"
    ]
    -
    check[
        "TotalIsoformDiscriminative"
    ]
)


checksum_failures = int(
    (
        (
            check[
                "AllPeptideDifference"
            ]
            !=
            0
        )
        |
        (
            check[
                "IsoformPeptideDifference"
            ]
            !=
            0
        )
    ).sum()
)


# =====================================================================
# STRICT CELL-LINE-SPECIFIC SOFTWARE CONSENSUS
# =====================================================================

specific = pd.read_csv(
    SPECIFIC_FILE,
    dtype=str,
    low_memory=False
)


for col in [
    "Program",
    "CellLine",
    "Peptide",
    "GeneAwareClass"
]:

    specific[
        col
    ] = (
        specific[
            col
        ]
        .fillna("")
        .astype(str)
        .str.strip()
    )


specific[
    "Program"
] = specific[
    "Program"
].str.upper()


specific[
    "Peptide"
] = specific[
    "Peptide"
].str.upper()


specific_all_membership = make_membership(
    specific,
    "Peptide",
    "Strict_CellLineSpecific_AllPeptides"
)


specific_iso = specific[
    specific[
        "GeneAwareClass"
    ].isin(
        ISO_CLASSES
    )
].copy()


specific_iso_membership = make_membership(
    specific_iso,
    "Peptide",
    "Strict_CellLineSpecific_IsoformPeptides"
)


specific_consensus = pd.concat(
    [
        consensus_summary(
            specific_all_membership
        ),

        consensus_summary(
            specific_iso_membership
        )
    ],
    ignore_index=True
)


specific_intersections = pd.concat(
    [
        intersection_summary(
            specific_all_membership,
            "Peptide"
        ),

        intersection_summary(
            specific_iso_membership,
            "Peptide"
        )
    ],
    ignore_index=True
)


# =====================================================================
# STRICT CELL-LINE-SPECIFIC ISOFORM ACCESSION CONSENSUS
# =====================================================================

if SPECIFIC_ISOFORM_FILE.exists():

    specific_isoform = pd.read_csv(
        SPECIFIC_ISOFORM_FILE,
        dtype=str,
        low_memory=False
    )


    if not specific_isoform.empty:

        for col in [
            "Program",
            "CellLine",
            "IsoformAccession"
        ]:

            specific_isoform[
                col
            ] = (
                specific_isoform[
                    col
                ]
                .fillna("")
                .astype(str)
                .str.strip()
            )


        specific_isoform[
            "Program"
        ] = specific_isoform[
            "Program"
        ].str.upper()


        specific_isoform_membership = make_membership(
            specific_isoform,
            "IsoformAccession",
            "Strict_CellLineSpecific_ImplicatedIsoforms"
        )


        specific_isoform_consensus = (
            consensus_summary(
                specific_isoform_membership
            )
        )

    else:

        specific_isoform_membership = (
            pd.DataFrame()
        )

        specific_isoform_consensus = (
            pd.DataFrame()
        )

else:

    specific_isoform_membership = (
        pd.DataFrame()
    )

    specific_isoform_consensus = (
        pd.DataFrame()
    )


# =====================================================================
# PRIMARY CROSS-SOFTWARE SUMMARY
# =====================================================================

summary_rows = []


for cell_line in CELL_LINES:

    for layer, membership in [

        (
            "Robust_AllPeptides",
            membership_all
        ),

        (
            "Robust_IsoformDiscriminativePeptides",
            membership_iso_peptide
        ),

        (
            "Robust_UniquelyResolvedIsoforms",
            membership_unique_isoform
        ),

        (
            "Robust_DiscriminativelyImplicatedIsoforms",
            membership_implicated_isoform
        )
    ]:

        x = membership[
            membership[
                "CellLine"
            ]
            ==
            cell_line
        ]


        union_n = len(
            x
        )


        all_four = int(
            (
                x[
                    "ProgramCount"
                ]
                ==
                4
            ).sum()
        )


        at_least_three = int(
            (
                x[
                    "ProgramCount"
                ]
                >=
                3
            ).sum()
        )


        at_least_two = int(
            (
                x[
                    "ProgramCount"
                ]
                >=
                2
            ).sum()
        )


        single_software = int(
            (
                x[
                    "ProgramCount"
                ]
                ==
                1
            ).sum()
        )


        summary_rows.append({

            "CellLine":
                cell_line,

            "Layer":
                layer,

            "UnionItems":
                union_n,

            "SharedAll4":
                all_four,

            "SharedAll4Percent":
                (
                    100
                    *
                    all_four
                    /
                    union_n

                    if union_n
                    else 0
                ),

            "SharedAtLeast3":
                at_least_three,

            "SharedAtLeast3Percent":
                (
                    100
                    *
                    at_least_three
                    /
                    union_n

                    if union_n
                    else 0
                ),

            "SharedAtLeast2":
                at_least_two,

            "SharedAtLeast2Percent":
                (
                    100
                    *
                    at_least_two
                    /
                    union_n

                    if union_n
                    else 0
                ),

            "SingleSoftwareOnly":
                single_software,

            "SingleSoftwareOnlyPercent":
                (
                    100
                    *
                    single_software
                    /
                    union_n

                    if union_n
                    else 0
                )
        })


primary_summary = pd.DataFrame(
    summary_rows
)


# =====================================================================
# QC
# =====================================================================

status = "PASS"


if duplicates != 0:
    status = "REVIEW"


if len(
    bad_rep_support
) != 0:
    status = "REVIEW"


if len(
    observed_pairs
) != 12:
    status = "REVIEW"


if checksum_failures != 0:
    status = "REVIEW"


qc = pd.DataFrame(
    [
        {
            "Metric":
                "InputSupportRows",

            "Value":
                len(
                    support
                )
        },

        {
            "Metric":
                "DuplicateProgramCellLinePeptides",

            "Value":
                duplicates
        },

        {
            "Metric":
                "InvalidReplicateSupportRows",

            "Value":
                len(
                    bad_rep_support
                )
        },

        {
            "Metric":
                "ObservedProgramCellLinePairs",

            "Value":
                len(
                    observed_pairs
                )
        },

        {
            "Metric":
                "Step2A_ChecksumFailures",

            "Value":
                checksum_failures
        },

        {
            "Metric":
                "RobustCriterion",

            "Value":
                "ReplicateSupport >= 2 of 3"
        },

        {
            "Metric":
                "FinalStatus",

            "Value":
                status
        }
    ]
)


# =====================================================================
# EXPORT TABLES
# =====================================================================

outputs = {

    "01_Robust_Counts_ByProgramCellLine.csv":
        robust_counts,

    "02_Membership_Robust_AllPeptides.csv":
        membership_all,

    "03_Membership_Robust_IsoformDiscriminativePeptides.csv":
        membership_iso_peptide,

    "04_Membership_Robust_UniquelyResolvedIsoforms.csv":
        membership_unique_isoform,

    "05_Membership_Robust_AllImplicatedIsoforms.csv":
        membership_implicated_isoform,

    "06_IntersectionSummary_AllLayers.csv":
        intersection_combined,

    "07_Consensus_1to4Software_AllLayers.csv":
        consensus_combined,

    "08_Pairwise_Jaccard_AllLayers.csv":
        pairwise_combined,

    "09_Robust_IsoformEvidence_Detail.csv":
        robust_isoform_detail,

    "10_Primary_CrossSoftware_Summary.csv":
        primary_summary,

    "11_StrictSpecific_Membership_AllPeptides.csv":
        specific_all_membership,

    "12_StrictSpecific_Membership_IsoformPeptides.csv":
        specific_iso_membership,

    "13_StrictSpecific_Consensus.csv":
        specific_consensus,

    "14_StrictSpecific_IntersectionSummary.csv":
        specific_intersections,

    "15_StrictSpecific_IsoformMembership.csv":
        specific_isoform_membership,

    "16_StrictSpecific_IsoformConsensus.csv":
        specific_isoform_consensus,

    "17_STEP2A_Checksum.csv":
        check,

    "18_STEP2B_QC.csv":
        qc
}


for name, df in outputs.items():

    df.to_csv(
        OUT / name,
        index=False,
        encoding="utf-8-sig"
    )


# =====================================================================
# PRINT RESULTS
# =====================================================================

print()
print("=" * 130)
print("ROBUST COUNTS — >=2/3 REPLICATES")
print("=" * 130)

print(
    robust_counts.to_string(
        index=False
    )
)


print()
print("=" * 130)
print("PRIMARY CROSS-SOFTWARE CONCORDANCE")
print("=" * 130)

print(
    primary_summary.to_string(
        index=False
    )
)


print()
print("=" * 130)
print("ROBUST ISOFORM-DISCRIMINATIVE PEPTIDES — SOFTWARE CONSENSUS")
print("=" * 130)

print(
    consensus_combined[
        consensus_combined[
            "Layer"
        ]
        ==
        "Robust_IsoformDiscriminativePeptides"
    ].to_string(
        index=False
    )
)


print()
print("=" * 130)
print("ROBUST IMPLICATED ISOFORMS — SOFTWARE CONSENSUS")
print("=" * 130)

print(
    consensus_combined[
        consensus_combined[
            "Layer"
        ]
        ==
        "Robust_DiscriminativelyImplicatedIsoforms"
    ].to_string(
        index=False
    )
)


print()
print("=" * 130)
print("STRICT CELL-LINE-SPECIFIC ISOFORM PEPTIDES — CONSENSUS")
print("=" * 130)

print(
    specific_consensus[
        specific_consensus[
            "Layer"
        ]
        ==
        "Strict_CellLineSpecific_IsoformPeptides"
    ].to_string(
        index=False
    )
)


# =====================================================================
# VISUALIZATION
# =====================================================================

plot_status = "PASS"


try:

    import matplotlib.pyplot as plt


    # -----------------------------------------------------------------
    # FIGURE A
    # Consensus distribution of robust isoform-discriminative peptides
    # -----------------------------------------------------------------

    plot_df = consensus_combined[
        consensus_combined[
            "Layer"
        ]
        ==
        "Robust_IsoformDiscriminativePeptides"
    ]


    fig, ax = plt.subplots(
        figsize=(
            8,
            5
        )
    )


    x = np.arange(
        len(
            CELL_LINES
        )
    )

    width = 0.18


    for i, nsoftware in enumerate(
        [
            1,
            2,
            3,
            4
        ]
    ):

        values = []

        for cell_line in CELL_LINES:

            row = plot_df[
                (
                    plot_df[
                        "CellLine"
                    ]
                    ==
                    cell_line
                )
                &
                (
                    plot_df[
                        "SoftwareCount"
                    ]
                    ==
                    nsoftware
                )
            ]


            values.append(
                int(
                    row.iloc[0][
                        "ItemCount"
                    ]
                )
                if len(
                    row
                )
                else 0
            )


        ax.bar(
            x
            +
            (
                i
                -
                1.5
            )
            *
            width,
            values,
            width,
            label=f"{nsoftware}/4 software"
        )


    ax.set_xticks(
        x
    )

    ax.set_xticklabels(
        CELL_LINES
    )

    ax.set_ylabel(
        "Robust isoform-discriminative peptides"
    )

    ax.set_xlabel(
        "Cell line"
    )

    ax.set_title(
        "Cross-software consensus of reproducible isoform evidence"
    )

    ax.legend(
        frameon=False
    )

    fig.tight_layout()

    fig.savefig(
        FIGDIR /
        "Figure3A_RobustIsoformPeptide_SoftwareConsensus.png",
        dpi=300
    )

    fig.savefig(
        FIGDIR /
        "Figure3A_RobustIsoformPeptide_SoftwareConsensus.pdf"
    )

    plt.close(
        fig
    )


    # -----------------------------------------------------------------
    # FIGURE B
    # Consensus distribution of implicated isoforms
    # -----------------------------------------------------------------

    plot_df = consensus_combined[
        consensus_combined[
            "Layer"
        ]
        ==
        "Robust_DiscriminativelyImplicatedIsoforms"
    ]


    fig, ax = plt.subplots(
        figsize=(
            8,
            5
        )
    )


    x = np.arange(
        len(
            CELL_LINES
        )
    )

    width = 0.18


    for i, nsoftware in enumerate(
        [
            1,
            2,
            3,
            4
        ]
    ):

        values = []

        for cell_line in CELL_LINES:

            row = plot_df[
                (
                    plot_df[
                        "CellLine"
                    ]
                    ==
                    cell_line
                )
                &
                (
                    plot_df[
                        "SoftwareCount"
                    ]
                    ==
                    nsoftware
                )
            ]


            values.append(
                int(
                    row.iloc[0][
                        "ItemCount"
                    ]
                )
                if len(
                    row
                )
                else 0
            )


        ax.bar(
            x
            +
            (
                i
                -
                1.5
            )
            *
            width,
            values,
            width,
            label=f"{nsoftware}/4 software"
        )


    ax.set_xticks(
        x
    )

    ax.set_xticklabels(
        CELL_LINES
    )

    ax.set_ylabel(
        "Discriminatively implicated isoforms"
    )

    ax.set_title(
        "Cross-software consensus of reproducibly supported isoforms"
    )

    ax.legend(
        frameon=False
    )

    fig.tight_layout()

    fig.savefig(
        FIGDIR /
        "Figure3B_RobustIsoformAccession_SoftwareConsensus.png",
        dpi=300
    )

    fig.savefig(
        FIGDIR /
        "Figure3B_RobustIsoformAccession_SoftwareConsensus.pdf"
    )

    plt.close(
        fig
    )


    # -----------------------------------------------------------------
    # JACCARD HEATMAPS — ROBUST ISOFORM PEPTIDES
    # -----------------------------------------------------------------

    for cell_line in CELL_LINES:

        xdf = pairwise_iso_peptide[
            pairwise_iso_peptide[
                "CellLine"
            ]
            ==
            cell_line
        ]


        matrix = np.zeros(
            (
                4,
                4
            )
        )


        for i, p1 in enumerate(
            PROGRAMS
        ):

            for j, p2 in enumerate(
                PROGRAMS
            ):

                row = xdf[
                    (
                        xdf[
                            "Program1"
                        ]
                        ==
                        p1
                    )
                    &
                    (
                        xdf[
                            "Program2"
                        ]
                        ==
                        p2
                    )
                ]


                matrix[
                    i,
                    j
                ] = float(
                    row.iloc[0][
                        "Jaccard"
                    ]
                )


        fig, ax = plt.subplots(
            figsize=(
                5.5,
                5
            )
        )


        image = ax.imshow(
            matrix,
            vmin=0,
            vmax=1
        )


        ax.set_xticks(
            np.arange(
                4
            )
        )

        ax.set_xticklabels(
            PROGRAMS
        )

        ax.set_yticks(
            np.arange(
                4
            )
        )

        ax.set_yticklabels(
            PROGRAMS
        )


        for i in range(
            4
        ):

            for j in range(
                4
            ):

                ax.text(
                    j,
                    i,
                    f"{matrix[i,j]:.2f}",
                    ha="center",
                    va="center"
                )


        ax.set_title(
            f"{cell_line}: robust isoform-peptide Jaccard"
        )

        fig.colorbar(
            image,
            ax=ax,
            label="Jaccard index"
        )

        fig.tight_layout()

        fig.savefig(
            FIGDIR /
            f"Figure3C_Jaccard_RobustIsoformPeptides_{cell_line}.png",
            dpi=300
        )

        fig.savefig(
            FIGDIR /
            f"Figure3C_Jaccard_RobustIsoformPeptides_{cell_line}.pdf"
        )

        plt.close(
            fig
        )


    # -----------------------------------------------------------------
    # UPSET-LIKE FUNCTION
    # -----------------------------------------------------------------

    def make_upset_like(
        intersection_df,
        cell_line,
        layer,
        filename,
        title
    ):

        xdf = intersection_df[
            (
                intersection_df[
                    "CellLine"
                ]
                ==
                cell_line
            )
            &
            (
                intersection_df[
                    "Layer"
                ]
                ==
                layer
            )
        ].copy()


        if xdf.empty:
            return


        xdf = xdf.sort_values(
            [
                "ItemCount",
                "ProgramCount"
            ],
            ascending=[
                False,
                False
            ]
        )


        # all 15 combinations maximum for four programs
        xdf = xdf.head(
            15
        ).reset_index(
            drop=True
        )


        fig = plt.figure(
            figsize=(
                10,
                7
            )
        )


        gs = fig.add_gridspec(
            2,
            1,
            height_ratios=[
                3,
                1.6
            ],
            hspace=0.05
        )


        ax_bar = fig.add_subplot(
            gs[0]
        )

        ax_mat = fig.add_subplot(
            gs[1]
        )


        xpos = np.arange(
            len(
                xdf
            )
        )


        ax_bar.bar(
            xpos,
            xdf[
                "ItemCount"
            ].values
        )


        ax_bar.set_ylabel(
            "Intersection size"
        )

        ax_bar.set_xticks(
            []
        )

        ax_bar.set_title(
            title
        )


        for xidx, value in enumerate(
            xdf[
                "ItemCount"
            ].values
        ):

            ax_bar.text(
                xidx,
                value,
                str(
                    int(
                        value
                    )
                ),
                ha="center",
                va="bottom",
                fontsize=8
            )


        ax_mat.set_xlim(
            -0.5,
            len(
                xdf
            )
            -
            0.5
        )

        ax_mat.set_ylim(
            len(
                PROGRAMS
            )
            -
            0.5,
            -0.5
        )


        ax_mat.set_yticks(
            np.arange(
                len(
                    PROGRAMS
                )
            )
        )

        ax_mat.set_yticklabels(
            PROGRAMS
        )

        ax_mat.set_xticks(
            xpos
        )

        ax_mat.set_xticklabels(
            [
                ""
                for _ in xpos
            ]
        )


        for xidx, combination in enumerate(
            xdf[
                "Combination"
            ]
        ):

            members = set(
                combination.split(
                    "+"
                )
            )


            active_y = []


            for yidx, program in enumerate(
                PROGRAMS
            ):

                if program in members:

                    ax_mat.scatter(
                        xidx,
                        yidx,
                        s=55
                    )

                    active_y.append(
                        yidx
                    )

                else:

                    ax_mat.scatter(
                        xidx,
                        yidx,
                        s=20,
                        alpha=0.25
                    )


            if len(
                active_y
            ) >= 2:

                ax_mat.plot(
                    [
                        xidx,
                        xidx
                    ],
                    [
                        min(
                            active_y
                        ),
                        max(
                            active_y
                        )
                    ],
                    linewidth=1.5
                )


        ax_mat.spines[
            "top"
        ].set_visible(
            False
        )

        ax_mat.spines[
            "right"
        ].set_visible(
            False
        )

        ax_mat.spines[
            "bottom"
        ].set_visible(
            False
        )


        fig.savefig(
            FIGDIR /
            f"{filename}.png",
            dpi=300,
            bbox_inches="tight"
        )

        fig.savefig(
            FIGDIR /
            f"{filename}.pdf",
            bbox_inches="tight"
        )

        plt.close(
            fig
        )


    # -----------------------------------------------------------------
    # UPSET-LIKE — ROBUST ISOFORM PEPTIDES
    # -----------------------------------------------------------------

    for cell_line in CELL_LINES:

        make_upset_like(
            intersection_combined,
            cell_line,
            "Robust_IsoformDiscriminativePeptides",
            f"Figure3D_UpSet_RobustIsoformPeptides_{cell_line}",
            f"{cell_line}: reproducible isoform-discriminative peptide overlap"
        )


    # -----------------------------------------------------------------
    # UPSET-LIKE — ROBUST IMPLICATED ISOFORMS
    # -----------------------------------------------------------------

    for cell_line in CELL_LINES:

        make_upset_like(
            intersection_combined,
            cell_line,
            "Robust_DiscriminativelyImplicatedIsoforms",
            f"Figure3E_UpSet_RobustIsoforms_{cell_line}",
            f"{cell_line}: reproducibly implicated isoform overlap"
        )


    # -----------------------------------------------------------------
    # STRICT CELL-LINE-SPECIFIC CONSENSUS
    # -----------------------------------------------------------------

    xdf = specific_consensus[
        specific_consensus[
            "Layer"
        ]
        ==
        "Strict_CellLineSpecific_IsoformPeptides"
    ]


    fig, ax = plt.subplots(
        figsize=(
            8,
            5
        )
    )


    x = np.arange(
        len(
            CELL_LINES
        )
    )

    width = 0.18


    for i, nsoftware in enumerate(
        [
            1,
            2,
            3,
            4
        ]
    ):

        values = []

        for cell_line in CELL_LINES:

            row = xdf[
                (
                    xdf[
                        "CellLine"
                    ]
                    ==
                    cell_line
                )
                &
                (
                    xdf[
                        "SoftwareCount"
                    ]
                    ==
                    nsoftware
                )
            ]


            values.append(
                int(
                    row.iloc[0][
                        "ItemCount"
                    ]
                )
                if len(
                    row
                )
                else 0
            )


        ax.bar(
            x
            +
            (
                i
                -
                1.5
            )
            *
            width,
            values,
            width,
            label=f"{nsoftware}/4 software"
        )


    ax.set_xticks(
        x
    )

    ax.set_xticklabels(
        CELL_LINES
    )

    ax.set_ylabel(
        "Strict cell-line-specific\nisoform-discriminative peptides"
    )

    ax.set_title(
        "Software consensus for cell-line-associated isoform evidence"
    )

    ax.legend(
        frameon=False
    )

    fig.tight_layout()

    fig.savefig(
        FIGDIR /
        "Figure4_StrictSpecific_IsoformPeptideConsensus.png",
        dpi=300
    )

    fig.savefig(
        FIGDIR /
        "Figure4_StrictSpecific_IsoformPeptideConsensus.pdf"
    )

    plt.close(
        fig
    )


except Exception as e:

    plot_status = (
        "PLOTS FAILED: "
        +
        repr(
            e
        )
    )


# =====================================================================
# STATUS FILE
# =====================================================================

with open(
    OUT /
    "19_STEP2B_STATUS.txt",
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
        "STEP 2B — CROSS-SOFTWARE CONCORDANCE\n\n"
    )

    fh.write(
        f"FINAL STATUS: {status}\n"
    )

    fh.write(
        f"PLOT STATUS: {plot_status}\n\n"
    )

    fh.write(
        "Primary evidence criterion: replicate support >=2/3 within each cell line.\n"
    )

    fh.write(
        "Isoform-discriminative peptide classes: single_isoform_unique + "
        "within_family_subset_discriminative.\n"
    )

    fh.write(
        "Uniquely resolved isoforms require at least one robust "
        "single_isoform_unique peptide.\n"
    )

    fh.write(
        "Discriminatively implicated isoforms include isoforms supported by "
        "single-isoform-unique or subset-discriminative peptides.\n"
    )


# =====================================================================
# FINAL
# =====================================================================

print()
print("=" * 130)
print("STEP 2B QC")
print("=" * 130)

print(
    qc.to_string(
        index=False
    )
)


print()
print("=" * 130)
print("FINAL STATUS:", status)
print("PLOT STATUS :", plot_status)
print("=" * 130)

print()
print("OUTPUT:")
print(OUT)

print()
print("FIGURES:")
print(FIGDIR)

print()
print("STEP 2B COMPLETE")

