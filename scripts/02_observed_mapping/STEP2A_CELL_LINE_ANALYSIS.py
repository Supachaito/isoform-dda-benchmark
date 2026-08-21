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
import pandas as pd
import numpy as np
import re


# =====================================================================
# PATHS
# =====================================================================

ROOT = _public_project_root()

MASTER = ROOT / "MASTER_TABLES_FINAL"

STEP1D = (
    MASTER /
    "STEP1D_FINAL_NORMALIZED"
)

INPUT = (
    STEP1D /
    "03_FINAL_PeptideMappingDetail.csv"
)

OUT = (
    MASTER /
    "STEP2A_CELL_LINE_ANALYSIS"
)

FIGDIR = (
    OUT /
    "FIGURES"
)

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


SAMPLE_META = {

    "C33A_1": (
        "C33A",
        1
    ),

    "C33A_2": (
        "C33A",
        2
    ),

    "C33A_3": (
        "C33A",
        3
    ),

    "HELA_1": (
        "HeLa",
        1
    ),

    "HELA_2": (
        "HeLa",
        2
    ),

    "HELA_3": (
        "HeLa",
        3
    ),

    "SIHA_1": (
        "SiHa",
        1
    ),

    "SIHA_2": (
        "SiHa",
        2
    ),

    "SIHA_3": (
        "SiHa",
        3
    )
}


ISO_CLASSES = {
    "single_isoform_unique",
    "within_family_subset_discriminative"
}


ISO_RE = re.compile(
    r"^(.+?)-([1-9][0-9]*)$"
)


# =====================================================================
# HELPERS
# =====================================================================

def is_isoform(acc):

    if pd.isna(acc):
        return False

    return (
        ISO_RE.fullmatch(
            str(acc).strip()
        )
        is not None
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


def suffix_accession_union(series):

    result = set()

    for value in series:

        for acc in split_accessions(
            value
        ):

            if is_isoform(
                acc
            ):

                result.add(
                    acc
                )

    return result


def summarize_set(
    support,
    minimum_support
):

    rows = []


    for program in PROGRAMS:

        for cell_line in CELL_LINES:

            x = support[
                (
                    support[
                        "Program"
                    ]
                    ==
                    program
                )
                &
                (
                    support[
                        "CellLine"
                    ]
                    ==
                    cell_line
                )
                &
                (
                    support[
                        "ReplicateSupport"
                    ]
                    >=
                    minimum_support
                )
            ].copy()


            counts = (
                x[
                    "GeneAwareClass"
                ]
                .value_counts()
                .to_dict()
            )


            single = x[
                x[
                    "GeneAwareClass"
                ]
                ==
                "single_isoform_unique"
            ]


            subset = x[
                x[
                    "GeneAwareClass"
                ]
                ==
                "within_family_subset_discriminative"
            ]


            single_isoforms = (
                suffix_accession_union(
                    single[
                        "MappedAccessions"
                    ]
                )
            )


            subset_isoforms = (
                suffix_accession_union(
                    subset[
                        "MappedAccessions"
                    ]
                )
            )


            supported_isoforms = (
                single_isoforms
                |
                subset_isoforms
            )


            total = x[
                "Peptide"
            ].nunique()


            single_n = counts.get(
                "single_isoform_unique",
                0
            )


            subset_n = counts.get(
                "within_family_subset_discriminative",
                0
            )


            iso_disc = (
                single_n
                +
                subset_n
            )


            rows.append({

                "Program":
                    program,

                "CellLine":
                    cell_line,

                "MinimumReplicateSupport":
                    minimum_support,

                "DistinctPeptides":
                    total,

                "SingleIsoformUnique":
                    single_n,

                "SubsetDiscriminative":
                    subset_n,

                "TotalIsoformDiscriminative":
                    iso_disc,

                "IsoformDiscriminativePercent":
                    (
                        100
                        *
                        iso_disc
                        /
                        total

                        if total
                        else 0
                    ),

                "WithinFamilySharedAll":
                    counts.get(
                        "within_family_shared_all",
                        0
                    ),

                "SameGeneMultiEntryShared":
                    counts.get(
                        "same_gene_multi_entry_shared",
                        0
                    ),

                "CrossGeneShared":
                    counts.get(
                        "cross_gene_shared",
                        0
                    ),

                "PartiallyGeneResolved":
                    counts.get(
                        "multi_entry_partially_gene_resolved",
                        0
                    ),

                "GeneUnresolved":
                    counts.get(
                        "multi_entry_gene_unresolved",
                        0
                    ),

                "SingleCanonicalUnique":
                    counts.get(
                        "single_canonical_unique",
                        0
                    ),

                "DistinctSingleIsoformAccessions":
                    len(
                        single_isoforms
                    ),

                "DistinctSubsetSupportedIsoforms":
                    len(
                        subset_isoforms
                    ),

                "DistinctDiscriminativelySupportedIsoforms":
                    len(
                        supported_isoforms
                    )
            })


    return pd.DataFrame(
        rows
    )


def complete_support_distribution(
    df
):

    rows = []


    for program in PROGRAMS:

        for cell_line in CELL_LINES:

            x = df[
                (
                    df[
                        "Program"
                    ]
                    ==
                    program
                )
                &
                (
                    df[
                        "CellLine"
                    ]
                    ==
                    cell_line
                )
            ]


            total = len(
                x
            )


            for n in [
                1,
                2,
                3
            ]:

                count = int(
                    (
                        x[
                            "ReplicateSupport"
                        ]
                        ==
                        n
                    ).sum()
                )


                rows.append({

                    "Program":
                        program,

                    "CellLine":
                        cell_line,

                    "ReplicateSupport":
                        n,

                    "PeptideCount":
                        count,

                    "Percent":
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


# =====================================================================
# LOAD FINAL STEP 1D DATA
# =====================================================================

print()
print("=" * 115)
print("STEP 2A — CELL-LINE AND REPLICATE-RESOLVED BENCHMARK")
print("=" * 115)

print()
print("Input:")
print(INPUT)


if not INPUT.exists():

    raise FileNotFoundError(
        INPUT
    )


detail = pd.read_csv(
    INPUT,
    dtype=str,
    low_memory=False
)


required = {
    "Program",
    "Sample",
    "Peptide",
    "GeneAwareClass",
    "MappedAccessions"
}


missing = (
    required
    -
    set(
        detail.columns
    )
)


if missing:

    raise RuntimeError(
        f"Missing required columns: {missing}"
    )


detail = detail[
    [
        "Program",
        "Sample",
        "Peptide",
        "GeneAwareClass",
        "MappedAccessions"
    ]
].copy()


for col in [
    "Program",
    "Sample",
    "Peptide",
    "GeneAwareClass"
]:

    detail[
        col
    ] = (
        detail[
            col
        ]
        .fillna("")
        .astype(str)
        .str.strip()
    )


detail[
    "Program"
] = detail[
    "Program"
].str.upper()


detail[
    "Peptide"
] = detail[
    "Peptide"
].str.upper()


detail[
    "SampleKey"
] = (
    detail[
        "Sample"
    ]
    .str.upper()
)


# =====================================================================
# SAMPLE METADATA
# =====================================================================

unknown_samples = sorted(
    set(
        detail[
            "SampleKey"
        ]
    )
    -
    set(
        SAMPLE_META
    )
)


if unknown_samples:

    raise RuntimeError(
        "Unknown sample names: "
        +
        str(
            unknown_samples
        )
    )


detail[
    "CellLine"
] = detail[
    "SampleKey"
].map(
    lambda x:
    SAMPLE_META[x][0]
)


detail[
    "Replicate"
] = detail[
    "SampleKey"
].map(
    lambda x:
    SAMPLE_META[x][1]
)


detail = detail[
    detail[
        "Program"
    ].isin(
        PROGRAMS
    )
].copy()


detail = (
    detail
    .drop_duplicates(
        subset=[
            "Program",
            "SampleKey",
            "Peptide"
        ]
    )
    .reset_index(
        drop=True
    )
)


print()
print(
    "Final peptide observations:",
    len(detail)
)


# =====================================================================
# CHECK PEPTIDE MAPPING CONSISTENCY
# =====================================================================

mapping_conflicts = (
    detail
    .groupby(
        "Peptide"
    )
    .agg(
        NClasses=(
            "GeneAwareClass",
            "nunique"
        ),
        NMappings=(
            "MappedAccessions",
            "nunique"
        )
    )
    .reset_index()
)


mapping_conflicts = mapping_conflicts[
    (
        mapping_conflicts[
            "NClasses"
        ]
        >
        1
    )
    |
    (
        mapping_conflicts[
            "NMappings"
        ]
        >
        1
    )
]


# =====================================================================
# REPLICATE COVERAGE QC
# =====================================================================

replicate_qc = (
    detail[
        [
            "Program",
            "CellLine",
            "Replicate"
        ]
    ]
    .drop_duplicates()
    .groupby(
        [
            "Program",
            "CellLine"
        ]
    )
    .size()
    .reset_index(
        name="ObservedReplicates"
    )
)


# =====================================================================
# PEPTIDE MAPPING TABLE
# =====================================================================

peptide_mapping = (
    detail[
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


# =====================================================================
# REPLICATE SUPPORT PER PROGRAM × CELL LINE × PEPTIDE
# =====================================================================

presence = (
    detail[
        [
            "Program",
            "CellLine",
            "Replicate",
            "Peptide"
        ]
    ]
    .drop_duplicates()
)


support = (
    presence
    .groupby(
        [
            "Program",
            "CellLine",
            "Peptide"
        ]
    )[
        "Replicate"
    ]
    .nunique()
    .reset_index(
        name="ReplicateSupport"
    )
)


support = support.merge(
    peptide_mapping,
    on="Peptide",
    how="left"
)


# =====================================================================
# SUMMARY — ANY REPLICATE
# =====================================================================

summary_any = summarize_set(
    support,
    minimum_support=1
)


# =====================================================================
# SUMMARY — ROBUST >=2/3
# =====================================================================

summary_robust = summarize_set(
    support,
    minimum_support=2
)


# =====================================================================
# REPLICATE SUPPORT DISTRIBUTION — ALL PEPTIDES
# =====================================================================

rep_support_all = (
    complete_support_distribution(
        support
    )
)


# =====================================================================
# REPLICATE SUPPORT DISTRIBUTION — ISOFORM-DISCRIMINATIVE
# =====================================================================

iso_support = support[
    support[
        "GeneAwareClass"
    ].isin(
        ISO_CLASSES
    )
].copy()


rep_support_iso = (
    complete_support_distribution(
        iso_support
    )
)


# =====================================================================
# CELL-LINE PATTERNS PER SOFTWARE
# ROBUST PRESENCE = >=2/3
# =====================================================================

pivot = (
    support
    .pivot_table(
        index=[
            "Program",
            "Peptide"
        ],
        columns="CellLine",
        values="ReplicateSupport",
        aggfunc="max",
        fill_value=0
    )
    .reset_index()
)


for cell_line in CELL_LINES:

    if cell_line not in pivot.columns:

        pivot[
            cell_line
        ] = 0


pivot = pivot.merge(
    peptide_mapping,
    on="Peptide",
    how="left"
)


for cell_line in CELL_LINES:

    pivot[
        f"{cell_line}_Robust"
    ] = (
        pivot[
            cell_line
        ]
        >=
        2
    )


def make_pattern(row):

    present = [

        cell_line

        for cell_line
        in CELL_LINES

        if row[
            f"{cell_line}_Robust"
        ]
    ]


    if len(
        present
    ) == 0:

        return (
            "No_2of3_CellLine"
        )


    return (
        "+".join(
            present
        )
    )


pivot[
    "RobustCellLinePattern"
] = pivot.apply(
    make_pattern,
    axis=1
)


pattern_all = (
    pivot
    .groupby(
        [
            "Program",
            "RobustCellLinePattern"
        ]
    )
    .size()
    .reset_index(
        name="PeptideCount"
    )
)


pattern_iso = (
    pivot[
        pivot[
            "GeneAwareClass"
        ].isin(
            ISO_CLASSES
        )
    ]
    .groupby(
        [
            "Program",
            "RobustCellLinePattern"
        ]
    )
    .size()
    .reset_index(
        name="PeptideCount"
    )
)


# =====================================================================
# STRICT CELL-LINE-SPECIFIC PEPTIDES
#
# Definition:
# target cell line >= 2/3 replicates
# AND both other cell lines == 0/3
# =====================================================================

specific_rows = []


for _, r in pivot.iterrows():

    supports = {

        cell_line:
            int(
                r[
                    cell_line
                ]
            )

        for cell_line
        in CELL_LINES
    }


    for target in CELL_LINES:

        others = [
            x
            for x in CELL_LINES
            if x != target
        ]


        if (
            supports[
                target
            ]
            >=
            2
            and
            all(
                supports[
                    other
                ]
                ==
                0

                for other
                in others
            )
        ):

            specific_rows.append({

                "Program":
                    r[
                        "Program"
                    ],

                "CellLine":
                    target,

                "Peptide":
                    r[
                        "Peptide"
                    ],

                "GeneAwareClass":
                    r[
                        "GeneAwareClass"
                    ],

                "MappedAccessions":
                    r[
                        "MappedAccessions"
                    ],

                "TargetReplicateSupport":
                    supports[
                        target
                    ],

                "C33A_Support":
                    supports[
                        "C33A"
                    ],

                "HeLa_Support":
                    supports[
                        "HeLa"
                    ],

                "SiHa_Support":
                    supports[
                        "SiHa"
                    ]
            })


specific = pd.DataFrame(
    specific_rows
)


if specific.empty:

    specific = pd.DataFrame(
        columns=[
            "Program",
            "CellLine",
            "Peptide",
            "GeneAwareClass",
            "MappedAccessions",
            "TargetReplicateSupport",
            "C33A_Support",
            "HeLa_Support",
            "SiHa_Support"
        ]
    )


# =====================================================================
# STRICT CELL-LINE-SPECIFIC SUMMARY
# =====================================================================

specific_summary_rows = []


for program in PROGRAMS:

    for cell_line in CELL_LINES:

        x = specific[
            (
                specific[
                    "Program"
                ]
                ==
                program
            )
            &
            (
                specific[
                    "CellLine"
                ]
                ==
                cell_line
            )
        ]


        iso_x = x[
            x[
                "GeneAwareClass"
            ].isin(
                ISO_CLASSES
            )
        ]


        specific_summary_rows.append({

            "Program":
                program,

            "CellLine":
                cell_line,

            "StrictSpecific_AllPeptides":
                x[
                    "Peptide"
                ].nunique(),

            "StrictSpecific_IsoformDiscriminativePeptides":
                iso_x[
                    "Peptide"
                ].nunique(),

            "StrictSpecific_SingleIsoformUnique":
                iso_x.loc[
                    iso_x[
                        "GeneAwareClass"
                    ]
                    ==
                    "single_isoform_unique",
                    "Peptide"
                ].nunique(),

            "StrictSpecific_SubsetDiscriminative":
                iso_x.loc[
                    iso_x[
                        "GeneAwareClass"
                    ]
                    ==
                    "within_family_subset_discriminative",
                    "Peptide"
                ].nunique()
        })


specific_summary = pd.DataFrame(
    specific_summary_rows
)


# =====================================================================
# CELL-LINE-SPECIFIC ISOFORM EVIDENCE
# =====================================================================

iso_specific_rows = []


for _, r in specific.iterrows():

    if (
        r[
            "GeneAwareClass"
        ]
        not in
        ISO_CLASSES
    ):

        continue


    for acc in split_accessions(
        r[
            "MappedAccessions"
        ]
    ):

        if not is_isoform(
            acc
        ):

            continue


        iso_specific_rows.append({

            "Program":
                r[
                    "Program"
                ],

            "CellLine":
                r[
                    "CellLine"
                ],

            "IsoformAccession":
                acc,

            "Peptide":
                r[
                    "Peptide"
                ],

            "EvidenceClass":
                r[
                    "GeneAwareClass"
                ],

            "TargetReplicateSupport":
                r[
                    "TargetReplicateSupport"
                ]
        })


specific_isoforms = pd.DataFrame(
    iso_specific_rows
)


if specific_isoforms.empty:

    specific_isoforms = pd.DataFrame(
        columns=[
            "Program",
            "CellLine",
            "IsoformAccession",
            "Peptide",
            "EvidenceClass",
            "TargetReplicateSupport"
        ]
    )


# =====================================================================
# CONSENSUS OF STRICT CELL-LINE-SPECIFIC PEPTIDES ACROSS SOFTWARE
# =====================================================================

if not specific.empty:

    consensus = (
        specific
        .groupby(
            [
                "CellLine",
                "Peptide"
            ]
        )
        .agg(
            SoftwareCount=(
                "Program",
                "nunique"
            ),
            Programs=(
                "Program",
                lambda x:
                ";".join(
                    sorted(
                        set(x)
                    )
                )
            ),
            GeneAwareClass=(
                "GeneAwareClass",
                "first"
            ),
            MappedAccessions=(
                "MappedAccessions",
                "first"
            )
        )
        .reset_index()
    )

else:

    consensus = pd.DataFrame(
        columns=[
            "CellLine",
            "Peptide",
            "SoftwareCount",
            "Programs",
            "GeneAwareClass",
            "MappedAccessions"
        ]
    )


consensus_summary = (
    consensus
    .groupby(
        [
            "CellLine",
            "SoftwareCount"
        ]
    )
    .size()
    .reset_index(
        name="PeptideCount"
    )
)


# =====================================================================
# ROBUST ISOFORM EVIDENCE BY CELL LINE
# =====================================================================

robust_iso_rows = []


for program in PROGRAMS:

    for cell_line in CELL_LINES:

        x = support[
            (
                support[
                    "Program"
                ]
                ==
                program
            )
            &
            (
                support[
                    "CellLine"
                ]
                ==
                cell_line
            )
            &
            (
                support[
                    "ReplicateSupport"
                ]
                >=
                2
            )
            &
            (
                support[
                    "GeneAwareClass"
                ].isin(
                    ISO_CLASSES
                )
            )
        ]


        for _, r in x.iterrows():

            for acc in split_accessions(
                r[
                    "MappedAccessions"
                ]
            ):

                if not is_isoform(
                    acc
                ):

                    continue


                robust_iso_rows.append({

                    "Program":
                        program,

                    "CellLine":
                        cell_line,

                    "IsoformAccession":
                        acc,

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


robust_iso_evidence = pd.DataFrame(
    robust_iso_rows
)


# =====================================================================
# SOFTWARE × CELL-LINE MATRICES
# =====================================================================

plot_total = (
    summary_any[
        [
            "Program",
            "CellLine",
            "DistinctPeptides"
        ]
    ]
    .copy()
)


plot_iso = (
    summary_any[
        [
            "Program",
            "CellLine",
            "SingleIsoformUnique",
            "SubsetDiscriminative",
            "TotalIsoformDiscriminative"
        ]
    ]
    .copy()
)


plot_fraction = (
    summary_robust[
        [
            "Program",
            "CellLine",
            "IsoformDiscriminativePercent"
        ]
    ]
    .copy()
)


# =====================================================================
# QC
# =====================================================================

observed_pairs = len(
    replicate_qc
)


all_three_reps = bool(
    (
        replicate_qc[
            "ObservedReplicates"
        ]
        ==
        3
    ).all()
)


status = (
    "PASS"

    if (
        len(
            unknown_samples
        )
        ==
        0
        and
        len(
            mapping_conflicts
        )
        ==
        0
        and
        observed_pairs
        ==
        12
        and
        all_three_reps
    )

    else
    "REVIEW"
)


qc = pd.DataFrame(
    [
        {
            "Metric":
                "InputRows",

            "Value":
                len(
                    detail
                )
        },
        {
            "Metric":
                "Programs",

            "Value":
                detail[
                    "Program"
                ].nunique()
        },
        {
            "Metric":
                "CellLines",

            "Value":
                detail[
                    "CellLine"
                ].nunique()
        },
        {
            "Metric":
                "ExpectedProgramCellLinePairs",

            "Value":
                12
        },
        {
            "Metric":
                "ObservedProgramCellLinePairs",

            "Value":
                observed_pairs
        },
        {
            "Metric":
                "AllProgramCellLinesHave3Replicates",

            "Value":
                all_three_reps
        },
        {
            "Metric":
                "PeptideMappingConflicts",

            "Value":
                len(
                    mapping_conflicts
                )
        },
        {
            "Metric":
                "StrictSpecificCriterion",

            "Value":
                ">=2/3 target cell-line replicates AND 0/3 in each other cell line"
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

    "01_CellLine_Summary_AnyReplicate.csv":
        summary_any,

    "02_CellLine_Summary_Robust2of3.csv":
        summary_robust,

    "03_Peptide_ReplicateSupport.csv":
        support,

    "04_ReplicateSupport_AllPeptides.csv":
        rep_support_all,

    "05_ReplicateSupport_IsoformDiscriminative.csv":
        rep_support_iso,

    "06_Robust_CellLinePatterns_AllPeptides.csv":
        pattern_all,

    "07_Robust_CellLinePatterns_IsoformDiscriminative.csv":
        pattern_iso,

    "08_Strict_CellLineSpecific_Peptides.csv":
        specific,

    "09_Strict_CellLineSpecific_Summary.csv":
        specific_summary,

    "10_Strict_CellLineSpecific_IsoformEvidence.csv":
        specific_isoforms,

    "11_Strict_CellLineSpecific_ConsensusAcrossSoftware.csv":
        consensus,

    "12_Strict_CellLineSpecific_ConsensusSummary.csv":
        consensus_summary,

    "13_Robust_IsoformEvidence_ByCellLine.csv":
        robust_iso_evidence,

    "14_PlotData_TotalPeptides_ByCellLine.csv":
        plot_total,

    "15_PlotData_IsoformEvidence_ByCellLine.csv":
        plot_iso,

    "16_PlotData_DiscriminativeFraction_Heatmap.csv":
        plot_fraction,

    "17_PlotData_ReplicateSupport_All.csv":
        rep_support_all,

    "18_PlotData_ReplicateSupport_Isoform.csv":
        rep_support_iso,

    "19_Replicate_QC.csv":
        replicate_qc,

    "20_PeptideMapping_Conflict_QC.csv":
        mapping_conflicts,

    "21_STEP2A_QC.csv":
        qc
}


for name, df in outputs.items():

    df.to_csv(
        OUT /
        name,
        index=False,
        encoding="utf-8-sig"
    )


# =====================================================================
# PRINT PRIMARY RESULTS
# =====================================================================

print()
print("=" * 125)
print("CELL-LINE SUMMARY — ANY OF 3 REPLICATES")
print("=" * 125)


show = [
    "Program",
    "CellLine",
    "DistinctPeptides",
    "SingleIsoformUnique",
    "SubsetDiscriminative",
    "TotalIsoformDiscriminative",
    "IsoformDiscriminativePercent",
    "DistinctDiscriminativelySupportedIsoforms"
]


print(
    summary_any[
        show
    ].to_string(
        index=False
    )
)


print()
print("=" * 125)
print("CELL-LINE SUMMARY — ROBUST >=2/3 REPLICATES")
print("=" * 125)

print(
    summary_robust[
        show
    ].to_string(
        index=False
    )
)


print()
print("=" * 125)
print("STRICT CELL-LINE-SPECIFIC")
print("Definition: >=2/3 in target cell line AND 0/3 in both other cell lines")
print("=" * 125)

print(
    specific_summary.to_string(
        index=False
    )
)


print()
print("=" * 125)
print("REPLICATE SUPPORT — ISOFORM-DISCRIMINATIVE PEPTIDES")
print("=" * 125)

print(
    rep_support_iso.to_string(
        index=False
    )
)


# =====================================================================
# CREATE PRELIMINARY FIGURES
# =====================================================================

try:

    import matplotlib.pyplot as plt


    # --------------------------------------------------------------
    # FIGURE 2A
    # TOTAL DISTINCT PEPTIDES BY CELL LINE
    # --------------------------------------------------------------

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


    for i, program in enumerate(
        PROGRAMS
    ):

        vals = []

        for cell_line in CELL_LINES:

            row = summary_any[
                (
                    summary_any[
                        "Program"
                    ]
                    ==
                    program
                )
                &
                (
                    summary_any[
                        "CellLine"
                    ]
                    ==
                    cell_line
                )
            ]

            vals.append(
                int(
                    row.iloc[0][
                        "DistinctPeptides"
                    ]
                )
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
            vals,
            width,
            label=program
        )


    ax.set_xticks(
        x
    )

    ax.set_xticklabels(
        CELL_LINES
    )

    ax.set_ylabel(
        "Distinct peptides"
    )

    ax.set_xlabel(
        "Cell line"
    )

    ax.set_title(
        "Peptide identification depth by cell line"
    )

    ax.legend(
        frameon=False
    )

    fig.tight_layout()

    fig.savefig(
        FIGDIR /
        "Figure2A_TotalPeptides_ByCellLine.png",
        dpi=300
    )

    fig.savefig(
        FIGDIR /
        "Figure2A_TotalPeptides_ByCellLine.pdf"
    )

    plt.close(
        fig
    )


    # --------------------------------------------------------------
    # FIGURE 2B
    # ISOFORM-DISCRIMINATIVE PEPTIDES BY CELL LINE
    # --------------------------------------------------------------

    combinations = [

        (
            cell_line,
            program
        )

        for cell_line
        in CELL_LINES

        for program
        in PROGRAMS
    ]


    labels = [

        f"{cell}\n{program}"

        for cell, program
        in combinations
    ]


    single_vals = []

    subset_vals = []


    for cell_line, program in combinations:

        row = summary_any[
            (
                summary_any[
                    "Program"
                ]
                ==
                program
            )
            &
            (
                summary_any[
                    "CellLine"
                ]
                ==
                cell_line
            )
        ].iloc[0]


        single_vals.append(
            int(
                row[
                    "SingleIsoformUnique"
                ]
            )
        )

        subset_vals.append(
            int(
                row[
                    "SubsetDiscriminative"
                ]
            )
        )


    fig, ax = plt.subplots(
        figsize=(
            11,
            5
        )
    )


    xpos = np.arange(
        len(
            combinations
        )
    )


    ax.bar(
        xpos,
        single_vals,
        label="Single-isoform unique"
    )


    ax.bar(
        xpos,
        subset_vals,
        bottom=single_vals,
        label="Subset-discriminative"
    )


    ax.set_xticks(
        xpos
    )

    ax.set_xticklabels(
        labels,
        rotation=45,
        ha="right"
    )

    ax.set_ylabel(
        "Isoform-discriminative peptides"
    )

    ax.set_title(
        "Isoform-discriminative peptide evidence by cell line"
    )

    ax.legend(
        frameon=False
    )

    fig.tight_layout()

    fig.savefig(
        FIGDIR /
        "Figure2B_IsoformDiscriminative_ByCellLine.png",
        dpi=300
    )

    fig.savefig(
        FIGDIR /
        "Figure2B_IsoformDiscriminative_ByCellLine.pdf"
    )

    plt.close(
        fig
    )


    # --------------------------------------------------------------
    # FIGURE 2C
    # ROBUST DISCRIMINATIVE FRACTION HEATMAP
    # --------------------------------------------------------------

    matrix = np.zeros(
        (
            len(
                PROGRAMS
            ),
            len(
                CELL_LINES
            )
        )
    )


    for i, program in enumerate(
        PROGRAMS
    ):

        for j, cell_line in enumerate(
            CELL_LINES
        ):

            row = summary_robust[
                (
                    summary_robust[
                        "Program"
                    ]
                    ==
                    program
                )
                &
                (
                    summary_robust[
                        "CellLine"
                    ]
                    ==
                    cell_line
                )
            ].iloc[0]


            matrix[
                i,
                j
            ] = float(
                row[
                    "IsoformDiscriminativePercent"
                ]
            )


    fig, ax = plt.subplots(
        figsize=(
            6,
            5
        )
    )


    image = ax.imshow(
        matrix,
        aspect="auto"
    )


    ax.set_xticks(
        np.arange(
            len(
                CELL_LINES
            )
        )
    )

    ax.set_xticklabels(
        CELL_LINES
    )

    ax.set_yticks(
        np.arange(
            len(
                PROGRAMS
            )
        )
    )

    ax.set_yticklabels(
        PROGRAMS
    )


    for i in range(
        len(
            PROGRAMS
        )
    ):

        for j in range(
            len(
                CELL_LINES
            )
        ):

            ax.text(
                j,
                i,
                f"{matrix[i,j]:.2f}%",
                ha="center",
                va="center"
            )


    ax.set_title(
        "Robust isoform-discriminative fraction (≥2/3 replicates)"
    )

    fig.colorbar(
        image,
        ax=ax,
        label="% of robust peptides"
    )

    fig.tight_layout()

    fig.savefig(
        FIGDIR /
        "Figure2C_Robust_DiscriminativeFraction_Heatmap.png",
        dpi=300
    )

    fig.savefig(
        FIGDIR /
        "Figure2C_Robust_DiscriminativeFraction_Heatmap.pdf"
    )

    plt.close(
        fig
    )


    # --------------------------------------------------------------
    # FIGURE 2D
    # REPLICATE REPRODUCIBILITY — ALL PEPTIDES
    # --------------------------------------------------------------

    combos = [

        (
            cell_line,
            program
        )

        for cell_line
        in CELL_LINES

        for program
        in PROGRAMS
    ]


    fig, ax = plt.subplots(
        figsize=(
            11,
            5
        )
    )


    xpos = np.arange(
        len(
            combos
        )
    )


    bottom = np.zeros(
        len(
            combos
        )
    )


    for support_n in [
        1,
        2,
        3
    ]:

        values = []

        for cell_line, program in combos:

            row = rep_support_all[
                (
                    rep_support_all[
                        "Program"
                    ]
                    ==
                    program
                )
                &
                (
                    rep_support_all[
                        "CellLine"
                    ]
                    ==
                    cell_line
                )
                &
                (
                    rep_support_all[
                        "ReplicateSupport"
                    ]
                    ==
                    support_n
                )
            ]

            values.append(
                float(
                    row.iloc[0][
                        "Percent"
                    ]
                )
            )


        ax.bar(
            xpos,
            values,
            bottom=bottom,
            label=f"{support_n}/3 replicates"
        )


        bottom += np.array(
            values
        )


    ax.set_xticks(
        xpos
    )

    ax.set_xticklabels(
        [
            f"{c}\n{p}"
            for c, p in combos
        ],
        rotation=45,
        ha="right"
    )

    ax.set_ylabel(
        "Peptides (%)"
    )

    ax.set_ylim(
        0,
        100
    )

    ax.set_title(
        "Replicate reproducibility of peptide identifications"
    )

    ax.legend(
        frameon=False
    )

    fig.tight_layout()

    fig.savefig(
        FIGDIR /
        "Figure2D_ReplicateReproducibility_AllPeptides.png",
        dpi=300
    )

    fig.savefig(
        FIGDIR /
        "Figure2D_ReplicateReproducibility_AllPeptides.pdf"
    )

    plt.close(
        fig
    )


    # --------------------------------------------------------------
    # FIGURE 2E
    # REPLICATE REPRODUCIBILITY — ISOFORM-DISCRIMINATIVE
    # --------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(
            11,
            5
        )
    )


    xpos = np.arange(
        len(
            combos
        )
    )


    bottom = np.zeros(
        len(
            combos
        )
    )


    for support_n in [
        1,
        2,
        3
    ]:

        values = []

        for cell_line, program in combos:

            row = rep_support_iso[
                (
                    rep_support_iso[
                        "Program"
                    ]
                    ==
                    program
                )
                &
                (
                    rep_support_iso[
                        "CellLine"
                    ]
                    ==
                    cell_line
                )
                &
                (
                    rep_support_iso[
                        "ReplicateSupport"
                    ]
                    ==
                    support_n
                )
            ]

            values.append(
                float(
                    row.iloc[0][
                        "Percent"
                    ]
                )
            )


        ax.bar(
            xpos,
            values,
            bottom=bottom,
            label=f"{support_n}/3 replicates"
        )


        bottom += np.array(
            values
        )


    ax.set_xticks(
        xpos
    )

    ax.set_xticklabels(
        [
            f"{c}\n{p}"
            for c, p in combos
        ],
        rotation=45,
        ha="right"
    )

    ax.set_ylabel(
        "Isoform-discriminative peptides (%)"
    )

    ax.set_ylim(
        0,
        100
    )

    ax.set_title(
        "Replicate reproducibility of isoform-discriminative evidence"
    )

    ax.legend(
        frameon=False
    )

    fig.tight_layout()

    fig.savefig(
        FIGDIR /
        "Figure2E_ReplicateReproducibility_IsoformEvidence.png",
        dpi=300
    )

    fig.savefig(
        FIGDIR /
        "Figure2E_ReplicateReproducibility_IsoformEvidence.pdf"
    )

    plt.close(
        fig
    )


    # --------------------------------------------------------------
    # FIGURE 2F
    # STRICT CELL-LINE-SPECIFIC ISOFORM-DISCRIMINATIVE PEPTIDES
    # --------------------------------------------------------------

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


    for i, program in enumerate(
        PROGRAMS
    ):

        vals = []

        for cell_line in CELL_LINES:

            row = specific_summary[
                (
                    specific_summary[
                        "Program"
                    ]
                    ==
                    program
                )
                &
                (
                    specific_summary[
                        "CellLine"
                    ]
                    ==
                    cell_line
                )
            ].iloc[0]


            vals.append(
                int(
                    row[
                        "StrictSpecific_IsoformDiscriminativePeptides"
                    ]
                )
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
            vals,
            width,
            label=program
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
        "Cell-line-specific isoform evidence"
    )

    ax.legend(
        frameon=False
    )

    fig.tight_layout()

    fig.savefig(
        FIGDIR /
        "Figure2F_Strict_CellLineSpecific_IsoformEvidence.png",
        dpi=300
    )

    fig.savefig(
        FIGDIR /
        "Figure2F_Strict_CellLineSpecific_IsoformEvidence.pdf"
    )

    plt.close(
        fig
    )


    plot_status = "PASS"


except Exception as e:

    plot_status = (
        "PLOTS FAILED: "
        +
        repr(
            e
        )
    )


# =====================================================================
# WRITE STATUS
# =====================================================================

with open(
    OUT /
    "22_STEP2A_STATUS.txt",
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
        "STEP 2A CELL-LINE ANALYSIS\n\n"
    )

    fh.write(
        f"FINAL STATUS: {status}\n"
    )

    fh.write(
        f"PLOT STATUS: {plot_status}\n\n"
    )

    fh.write(
        "Strict cell-line-specific definition:\n"
    )

    fh.write(
        "Detected in >=2/3 replicates of the target cell line "
        "AND detected in 0/3 replicates of each of the other two cell lines.\n"
    )


# =====================================================================
# FINAL PRINT
# =====================================================================

print()
print("=" * 125)
print("STEP 2A QC")
print("=" * 125)

print(
    qc.to_string(
        index=False
    )
)


print()
print("=" * 125)
print("FINAL STATUS:", status)
print("PLOT STATUS :", plot_status)
print("=" * 125)

print()
print("OUTPUT DIRECTORY:")
print(OUT)

print()
print("FIGURES:")
print(FIGDIR)

print()
print("STEP 2A COMPLETE")

