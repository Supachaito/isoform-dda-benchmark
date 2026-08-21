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
# STEP 4D
# Structural overlap analysis of reciprocal SEPEP events
#
# PURPOSE
# -------
# For each strong reciprocal event:
#
#   SEPEP-UP   log2FC >= +1
#   SEPEP-DOWN log2FC <= -1
#
# compare their exact UniProt accession memberships.
#
# Relationship:
#
#   DISJOINT
#       no shared accessions
#
#   PARTIAL_OVERLAP
#       some shared accessions, neither set contains the other
#
#   NESTED
#       one accession set is a strict subset of the other
#
#   IDENTICAL
#       should not normally occur because distinct SEPEPs
#       are defined by distinct accession sets
#
# IMPORTANT:
# This remains "switch-like structural regulation".
# It is NOT automatically confirmed biological isoform switching.
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

STEP4C = (
    MASTER /
    "STEP4C_RECIPROCAL_SWITCHLIKE"
)

META_FILE = (
    STEP4B /
    "01_SEPEP_StructuralGroup_Metadata.csv"
)

STRONG_FILE = (
    STEP4C /
    "02_Strong_Reciprocal_SwitchLike_BySoftware.csv"
)

STRICT_PAIR_FILE = (
    STEP4C /
    "06_Strict_SEPEP_Pair_Consensus.csv"
)

OUT = (
    MASTER /
    "STEP4D_STRUCTURAL_OVERLAP"
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


# ============================================================
# 2. HELPERS
# ============================================================

def parse_accession_set(value):

    if pd.isna(value):
        return set()

    value = str(value).strip()

    if value == "":
        return set()

    return {
        x.strip()
        for x in value.split("|")
        if x.strip()
    }


def accession_string(values):

    return ";".join(
        sorted(values)
    )


def classify_relation(
    up_set,
    down_set
):

    if len(up_set) == 0 or len(down_set) == 0:

        return {
            "StructuralRelation":
                "Unresolved",

            "IntersectionN":
                np.nan,

            "UnionN":
                np.nan,

            "Jaccard":
                np.nan,

            "OverlapCoefficient":
                np.nan,

            "UpUniqueN":
                np.nan,

            "DownUniqueN":
                np.nan,

            "NestedDirection":
                ""
        }


    intersection = (
        up_set
        &
        down_set
    )

    union = (
        up_set
        |
        down_set
    )


    intersection_n = len(
        intersection
    )

    union_n = len(
        union
    )


    jaccard = (
        intersection_n
        /
        union_n
        if union_n > 0
        else np.nan
    )


    min_n = min(
        len(up_set),
        len(down_set)
    )


    overlap_coefficient = (
        intersection_n
        /
        min_n
        if min_n > 0
        else np.nan
    )


    up_unique_n = len(
        up_set
        -
        down_set
    )


    down_unique_n = len(
        down_set
        -
        up_set
    )


    # --------------------------------------------------------
    # IDENTICAL
    # --------------------------------------------------------

    if up_set == down_set:

        relation = "Identical"

        nested_direction = ""


    # --------------------------------------------------------
    # DISJOINT
    # strongest structural separation
    # --------------------------------------------------------

    elif intersection_n == 0:

        relation = "Disjoint"

        nested_direction = ""


    # --------------------------------------------------------
    # NESTED
    # --------------------------------------------------------

    elif up_set < down_set:

        relation = "Nested"

        nested_direction = "UpSubsetOfDown"


    elif down_set < up_set:

        relation = "Nested"

        nested_direction = "DownSubsetOfUp"


    # --------------------------------------------------------
    # PARTIAL OVERLAP
    # --------------------------------------------------------

    else:

        relation = "Partial overlap"

        nested_direction = ""


    return {

        "StructuralRelation":
            relation,

        "IntersectionN":
            intersection_n,

        "UnionN":
            union_n,

        "Jaccard":
            jaccard,

        "OverlapCoefficient":
            overlap_coefficient,

        "UpUniqueN":
            up_unique_n,

        "DownUniqueN":
            down_unique_n,

        "NestedDirection":
            nested_direction
    }


def resolution_level(
    pair_class
):

    if pair_class == "C2-C2":

        return (
            "Exact-accession vs exact-accession"
        )


    if pair_class == "C2-C3":

        return (
            "Exact-accession vs isoform-subset"
        )


    if pair_class == "C3-C3":

        return (
            "Isoform-subset vs isoform-subset"
        )


    return "Other"


# ============================================================
# 3. LOAD INPUTS
# ============================================================

print()
print("=" * 120)
print("STEP 4D — STRUCTURAL OVERLAP OF RECIPROCAL SEPEP EVENTS")
print("=" * 120)


if not META_FILE.exists():

    raise FileNotFoundError(
        META_FILE
    )


if not STRONG_FILE.exists():

    raise FileNotFoundError(
        STRONG_FILE
    )


meta = pd.read_csv(
    META_FILE,
    low_memory=False
)


strong = pd.read_csv(
    STRONG_FILE,
    low_memory=False
)


# ============================================================
# 4. VALIDATE METADATA
# ============================================================

required_meta = {
    "SEPEPID",
    "StructuralKey",
    "Gene",
    "SEPEPClass"
}


missing_meta = (
    required_meta
    -
    set(meta.columns)
)


if missing_meta:

    raise RuntimeError(
        f"Missing metadata columns: "
        f"{sorted(missing_meta)}"
    )


required_strong = {
    "Program",
    "Gene",
    "Contrast",
    "MaxSEPEP",
    "MaxClass",
    "MaxSEPEP_log2FC",
    "MinSEPEP",
    "MinClass",
    "MinSEPEP_log2FC",
    "PairClass",
    "Gene_log2FC",
    "WithinGeneSpan"
}


missing_strong = (
    required_strong
    -
    set(strong.columns)
)


if missing_strong:

    raise RuntimeError(
        f"Missing strong-event columns: "
        f"{sorted(missing_strong)}"
    )


# ============================================================
# 5. NORMALIZE
# ============================================================

meta["SEPEPID"] = (
    meta["SEPEPID"]
    .astype(str)
    .str.strip()
)


meta["Gene"] = (
    meta["Gene"]
    .astype(str)
    .str.strip()
)


strong["Program"] = (
    strong["Program"]
    .astype(str)
    .str.upper()
    .str.strip()
)


strong["Gene"] = (
    strong["Gene"]
    .astype(str)
    .str.strip()
)


strong["MaxSEPEP"] = (
    strong["MaxSEPEP"]
    .astype(str)
    .str.strip()
)


strong["MinSEPEP"] = (
    strong["MinSEPEP"]
    .astype(str)
    .str.strip()
)


# ============================================================
# 6. CHECK UNIQUE SEPEP METADATA
# ============================================================

meta_dup = (
    meta
    .groupby(
        "SEPEPID"
    )
    .size()
    .reset_index(
        name="N"
    )
)


meta_dup_bad = meta_dup[
    meta_dup["N"] > 1
]


if len(
    meta_dup_bad
) > 0:

    meta_dup_bad.to_csv(
        OUT /
        "00_Duplicate_SEPEP_Metadata_QC.csv",
        index=False
    )

    raise RuntimeError(
        f"Duplicate SEPEP metadata IDs: "
        f"{len(meta_dup_bad)}"
    )


meta_small = (
    meta[
        [
            "SEPEPID",
            "StructuralKey",
            "Gene",
            "SEPEPClass"
        ]
    ]
    .drop_duplicates(
        subset=[
            "SEPEPID"
        ]
    )
)


# ============================================================
# 7. MERGE UP-SEPEP METADATA
# ============================================================

up_meta = meta_small.rename(
    columns={

        "SEPEPID":
            "MaxSEPEP",

        "StructuralKey":
            "UpStructuralKey",

        "Gene":
            "UpMetaGene",

        "SEPEPClass":
            "UpMetaClass"
    }
)


x = strong.merge(
    up_meta,
    on="MaxSEPEP",
    how="left"
)


# ============================================================
# 8. MERGE DOWN-SEPEP METADATA
# ============================================================

down_meta = meta_small.rename(
    columns={

        "SEPEPID":
            "MinSEPEP",

        "StructuralKey":
            "DownStructuralKey",

        "Gene":
            "DownMetaGene",

        "SEPEPClass":
            "DownMetaClass"
    }
)


x = x.merge(
    down_meta,
    on="MinSEPEP",
    how="left"
)


# ============================================================
# 9. MAPPING QC
# ============================================================

x[
    "MissingUpMetadata"
] = x[
    "UpStructuralKey"
].isna()


x[
    "MissingDownMetadata"
] = x[
    "DownStructuralKey"
].isna()


missing_mapping = x[
    x[
        "MissingUpMetadata"
    ]
    |
    x[
        "MissingDownMetadata"
    ]
]


if len(
    missing_mapping
) > 0:

    missing_mapping.to_csv(
        OUT /
        "01_Missing_SEPEP_Metadata_QC.csv",
        index=False
    )


# ============================================================
# 10. GENE-CONSISTENCY QC
# ============================================================

x[
    "UpGeneMismatch"
] = (
    x[
        "UpMetaGene"
    ].notna()
    &
    (
        x[
            "UpMetaGene"
        ]
        !=
        x[
            "Gene"
        ]
    )
)


x[
    "DownGeneMismatch"
] = (
    x[
        "DownMetaGene"
    ].notna()
    &
    (
        x[
            "DownMetaGene"
        ]
        !=
        x[
            "Gene"
        ]
    )
)


gene_mismatch = x[
    x[
        "UpGeneMismatch"
    ]
    |
    x[
        "DownGeneMismatch"
    ]
]


if len(
    gene_mismatch
) > 0:

    gene_mismatch.to_csv(
        OUT /
        "02_GeneMismatch_QC.csv",
        index=False
    )


# ============================================================
# 11. CALCULATE ACCESSION-SET RELATIONSHIP
# ============================================================

records = []


for row in x.itertuples(
    index=False
):

    up_set = parse_accession_set(
        row.UpStructuralKey
    )


    down_set = parse_accession_set(
        row.DownStructuralKey
    )


    relation = classify_relation(
        up_set,
        down_set
    )


    record = row._asdict()


    record.update({

        "UpAccessions":
            accession_string(
                up_set
            ),

        "DownAccessions":
            accession_string(
                down_set
            ),

        "UpAccessionN":
            len(up_set),

        "DownAccessionN":
            len(down_set),

        "SharedAccessions":
            accession_string(
                up_set
                &
                down_set
            ),

        "UpOnlyAccessions":
            accession_string(
                up_set
                -
                down_set
            ),

        "DownOnlyAccessions":
            accession_string(
                down_set
                -
                up_set
            ),

        "ResolutionLevel":
            resolution_level(
                row.PairClass
            )
    })


    record.update(
        relation
    )


    records.append(
        record
    )


overlap = pd.DataFrame(
    records
)


# ============================================================
# 12. STRUCTURAL INTERPRETATION
# ============================================================

def interpretation(row):

    rel = row[
        "StructuralRelation"
    ]

    pair = row[
        "PairClass"
    ]


    if rel == "Disjoint":

        if pair == "C2-C2":

            return (
                "Strongest: reciprocal exact-accession separation"
            )


        if pair == "C2-C3":

            return (
                "Strong: exact accession versus non-overlapping "
                "isoform subset"
            )


        if pair == "C3-C3":

            return (
                "Strong subset-level evidence: reciprocal "
                "non-overlapping isoform subsets"
            )


    if rel == "Partial overlap":

        return (
            "Intermediate: reciprocal structural units share "
            "some isoform accessions"
        )


    if rel == "Nested":

        return (
            "Weak for switching: one reciprocal structural "
            "unit is contained within the other"
        )


    if rel == "Identical":

        return (
            "Invalid structural contrast: identical accession sets"
        )


    return "Unresolved"


overlap[
    "StructuralInterpretation"
] = overlap.apply(
    interpretation,
    axis=1
)


# ============================================================
# 13. RELATION SUMMARY BY SOFTWARE
# ============================================================

relation_counts = (
    overlap
    .groupby(
        [
            "Program",
            "StructuralRelation"
        ],
        as_index=False
    )
    .agg(
        Events=(
            "Gene",
            "size"
        )
    )
)


totals = (
    overlap
    .groupby(
        "Program",
        as_index=False
    )
    .agg(
        TotalStrongEvents=(
            "Gene",
            "size"
        )
    )
)


relation_summary = relation_counts.merge(
    totals,
    on="Program",
    how="left"
)


relation_summary[
    "Percent"
] = (
    100
    *
    relation_summary[
        "Events"
    ]
    /
    relation_summary[
        "TotalStrongEvents"
    ]
)


# ============================================================
# 14. SOFTWARE-LEVEL STRUCTURAL QUALITY SUMMARY
# ============================================================

software_quality = []


for program in PROGRAMS:

    p = overlap[
        overlap[
            "Program"
        ]
        ==
        program
    ]


    if len(p) == 0:

        continue


    software_quality.append({

        "Program":
            program,

        "StrongEvents":
            len(p),

        "DisjointEvents":
            int(
                (
                    p[
                        "StructuralRelation"
                    ]
                    ==
                    "Disjoint"
                ).sum()
            ),

        "DisjointPercent":
            100
            *
            (
                p[
                    "StructuralRelation"
                ]
                ==
                "Disjoint"
            ).mean(),

        "PartialOverlapEvents":
            int(
                (
                    p[
                        "StructuralRelation"
                    ]
                    ==
                    "Partial overlap"
                ).sum()
            ),

        "NestedEvents":
            int(
                (
                    p[
                        "StructuralRelation"
                    ]
                    ==
                    "Nested"
                ).sum()
            ),

        "MedianJaccard":
            p[
                "Jaccard"
            ].median(),

        "MedianWithinGeneSpan":
            p[
                "WithinGeneSpan"
            ].median(),

        "DisjointWithAnyC2":
            int(
                (
                    (
                        p[
                            "StructuralRelation"
                        ]
                        ==
                        "Disjoint"
                    )
                    &
                    p[
                        "PairClass"
                    ].isin(
                        [
                            "C2-C2",
                            "C2-C3"
                        ]
                    )
                ).sum()
            )
    })


software_quality = pd.DataFrame(
    software_quality
)


# ============================================================
# 15. GENE × CONTRAST CONSENSUS
# ============================================================

consensus_rows = []


for (
    gene,
    contrast
), g in overlap.groupby(
    [
        "Gene",
        "Contrast"
    ]
):

    software_n = g[
        "Program"
    ].nunique()


    disjoint_n = int(
        (
            g[
                "StructuralRelation"
            ]
            ==
            "Disjoint"
        ).sum()
    )


    partial_n = int(
        (
            g[
                "StructuralRelation"
            ]
            ==
            "Partial overlap"
        ).sum()
    )


    nested_n = int(
        (
            g[
                "StructuralRelation"
            ]
            ==
            "Nested"
        ).sum()
    )


    unresolved_n = int(
        (
            g[
                "StructuralRelation"
            ]
            ==
            "Unresolved"
        ).sum()
    )


    relation_set = sorted(
        set(
            g[
                "StructuralRelation"
            ].dropna()
        )
    )


    if len(
        relation_set
    ) == 1:

        relation_consensus = (
            relation_set[0]
        )

    else:

        relation_consensus = (
            "Mixed"
        )


    consensus_rows.append({

        "Gene":
            gene,

        "Contrast":
            contrast,

        "SoftwareN":
            software_n,

        "Programs":
            "|".join(
                sorted(
                    set(
                        g[
                            "Program"
                        ]
                    )
                )
            ),

        "RelationConsensus":
            relation_consensus,

        "RelationClasses":
            "|".join(
                relation_set
            ),

        "DisjointSoftwareN":
            disjoint_n,

        "PartialOverlapSoftwareN":
            partial_n,

        "NestedSoftwareN":
            nested_n,

        "UnresolvedSoftwareN":
            unresolved_n,

        "AllSupportingSoftwareDisjoint":
            (
                disjoint_n
                ==
                software_n
            ),

        "MedianJaccard":
            g[
                "Jaccard"
            ].median(),

        "MaxJaccard":
            g[
                "Jaccard"
            ].max(),

        "MedianMaxSEPEP_log2FC":
            g[
                "MaxSEPEP_log2FC"
            ].median(),

        "MedianMinSEPEP_log2FC":
            g[
                "MinSEPEP_log2FC"
            ].median(),

        "MedianGene_log2FC":
            g[
                "Gene_log2FC"
            ].median(),

        "MedianWithinGeneSpan":
            g[
                "WithinGeneSpan"
            ].median(),

        "C2InvolvedSoftwareN":
            int(
                g[
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
    })


consensus = pd.DataFrame(
    consensus_rows
)


consensus = consensus.sort_values(
    [
        "SoftwareN",
        "AllSupportingSoftwareDisjoint",
        "DisjointSoftwareN",
        "C2InvolvedSoftwareN",
        "MedianWithinGeneSpan"
    ],
    ascending=[
        False,
        False,
        False,
        False,
        False
    ]
)


# ============================================================
# 16. TOP STRUCTURALLY SEPARATED CONSENSUS EVENTS
#
# strongest class:
#   >=2 software
#   ALL supporting workflows classify UP/DOWN sets as disjoint
# ============================================================

top_disjoint = consensus[
    (
        consensus[
            "SoftwareN"
        ]
        >=
        2
    )
    &
    (
        consensus[
            "AllSupportingSoftwareDisjoint"
        ]
    )
].copy()


top_disjoint = top_disjoint.sort_values(
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
# 17. ACCESSION-MEMBERSHIP LONG TABLE
#
# useful for supplement and manual candidate inspection
# ============================================================

membership_rows = []


for row in overlap.itertuples(
    index=False
):

    up_set = parse_accession_set(
        row.UpStructuralKey
    )


    down_set = parse_accession_set(
        row.DownStructuralKey
    )


    all_acc = sorted(
        up_set
        |
        down_set
    )


    for accession in all_acc:

        membership_rows.append({

            "Program":
                row.Program,

            "Gene":
                row.Gene,

            "Contrast":
                row.Contrast,

            "PairClass":
                row.PairClass,

            "StructuralRelation":
                row.StructuralRelation,

            "Accession":
                accession,

            "InUpSEPEP":
                accession in up_set,

            "InDownSEPEP":
                accession in down_set,

            "MembershipClass":
                (
                    "Shared"
                    if (
                        accession in up_set
                        and
                        accession in down_set
                    )
                    else
                    "Up only"
                    if accession in up_set
                    else
                    "Down only"
                )
        })


membership_long = pd.DataFrame(
    membership_rows
)


# ============================================================
# 18. STRICT EXACT-PAIR CONSENSUS ANNOTATION
# ============================================================

strict_annotated = pd.DataFrame()


if STRICT_PAIR_FILE.exists():

    strict = pd.read_csv(
        STRICT_PAIR_FILE,
        low_memory=False
    )


    if (
        "MaxSEPEP" in strict.columns
        and
        "MinSEPEP" in strict.columns
    ):

        strict = strict.merge(
            up_meta,
            on="MaxSEPEP",
            how="left"
        )


        strict = strict.merge(
            down_meta,
            on="MinSEPEP",
            how="left"
        )


        strict_rows = []


        for row in strict.itertuples(
            index=False
        ):

            up_set = parse_accession_set(
                row.UpStructuralKey
            )


            down_set = parse_accession_set(
                row.DownStructuralKey
            )


            relation = classify_relation(
                up_set,
                down_set
            )


            rec = row._asdict()


            rec.update({

                "UpAccessions":
                    accession_string(
                        up_set
                    ),

                "DownAccessions":
                    accession_string(
                        down_set
                    )
            })


            rec.update(
                relation
            )


            strict_rows.append(
                rec
            )


        strict_annotated = pd.DataFrame(
            strict_rows
        )


# ============================================================
# 19. EXPORT
# ============================================================

overlap.to_csv(
    OUT /
    "01_StrongReciprocal_StructuralOverlap_BySoftware.csv",
    index=False,
    encoding="utf-8-sig"
)


relation_summary.to_csv(
    OUT /
    "02_RelationSummary_BySoftware.csv",
    index=False,
    encoding="utf-8-sig"
)


software_quality.to_csv(
    OUT /
    "03_StructuralQuality_Summary_BySoftware.csv",
    index=False,
    encoding="utf-8-sig"
)


consensus.to_csv(
    OUT /
    "04_Consensus_StructuralOverlap_GeneContrast.csv",
    index=False,
    encoding="utf-8-sig"
)


top_disjoint.to_csv(
    OUT /
    "05_Top_Disjoint_Consensus_Candidates.csv",
    index=False,
    encoding="utf-8-sig"
)


membership_long.to_csv(
    OUT /
    "06_AccessionMembership_Long.csv",
    index=False,
    encoding="utf-8-sig"
)


if len(
    strict_annotated
) > 0:

    strict_annotated.to_csv(
        OUT /
        "07_StrictPair_StructuralOverlap.csv",
        index=False,
        encoding="utf-8-sig"
    )


# ============================================================
# 20. QC
# ============================================================

qc = pd.DataFrame({

    "Check": [
        "StrongInputRows",
        "OutputOverlapRows",
        "MissingUpMetadata",
        "MissingDownMetadata",
        "GeneMismatchRows",
        "IdenticalAccessionSets",
        "ProgramsPresent",
        "Consensus2plusDisjoint",
        "Consensus3plusDisjoint",
        "Consensus4Disjoint"
    ],

    "Value": [

        len(
            strong
        ),

        len(
            overlap
        ),

        int(
            overlap[
                "MissingUpMetadata"
            ].sum()
        ),

        int(
            overlap[
                "MissingDownMetadata"
            ].sum()
        ),

        len(
            gene_mismatch
        ),

        int(
            (
                overlap[
                    "StructuralRelation"
                ]
                ==
                "Identical"
            ).sum()
        ),

        overlap[
            "Program"
        ].nunique(),

        int(
            (
                (
                    consensus[
                        "SoftwareN"
                    ]
                    >=
                    2
                )
                &
                (
                    consensus[
                        "AllSupportingSoftwareDisjoint"
                ]
                )
            ).sum()
        ),

        int(
            (
                (
                    consensus[
                        "SoftwareN"
                    ]
                    >=
                    3
                )
                &
                (
                    consensus[
                        "AllSupportingSoftwareDisjoint"
                ]
                )
            ).sum()
        ),

        int(
            (
                (
                    consensus[
                        "SoftwareN"
                    ]
                    ==
                    4
                )
                &
                (
                    consensus[
                        "AllSupportingSoftwareDisjoint"
                ]
                )
            ).sum()
        )
    ]
})


qc.to_csv(
    OUT /
    "08_STEP4D_QC.csv",
    index=False,
    encoding="utf-8-sig"
)


status = (
    "PASS"
    if (
        len(
            strong
        )
        ==
        len(
            overlap
        )
        and
        overlap[
            "MissingUpMetadata"
        ].sum()
        ==
        0
        and
        overlap[
            "MissingDownMetadata"
        ].sum()
        ==
        0
        and
        len(
            gene_mismatch
        )
        ==
        0
        and
        overlap[
            "Program"
        ].nunique()
        ==
        4
    )
    else
    "FAIL"
)


with open(
    OUT /
    "09_STEP4D_STATUS.txt",
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
        f"STEP4D_STATUS={status}\n"
    )

    fh.write(
        "DISJOINT = reciprocal SEPEP accession sets share zero "
        "protein accessions; strongest structural separation.\n"
    )

    fh.write(
        "PARTIAL OVERLAP = accession sets overlap but neither "
        "contains the other.\n"
    )

    fh.write(
        "NESTED = one accession set is a strict subset of the other.\n"
    )

    fh.write(
        "Do not equate C3-C3 disjoint regulation with confirmed "
        "exact isoform switching.\n"
    )


# ============================================================
# 21. PRINT RESULTS
# ============================================================

print()
print("=" * 120)
print("STRUCTURAL RELATIONSHIP BY SOFTWARE")
print("=" * 120)

print(
    relation_summary.to_string(
        index=False
    )
)


print()
print("=" * 120)
print("STRUCTURAL QUALITY SUMMARY")
print("=" * 120)

print(
    software_quality.to_string(
        index=False
    )
)


print()
print("=" * 120)
print("TOP CONSENSUS EVENTS WITH DISJOINT ACCESSION SETS")
print("=" * 120)

if len(
    top_disjoint
) > 0:

    print(
        top_disjoint.head(
            30
        ).to_string(
            index=False
        )
    )

else:

    print(
        "No >=2-software fully disjoint consensus events."
    )


if len(
    strict_annotated
) > 0:

    print()
    print("=" * 120)
    print("TOP STRICT EXACT-PAIR STRUCTURAL RELATIONSHIPS")
    print("=" * 120)

    print(
        strict_annotated.head(
            30
        ).to_string(
            index=False
        )
    )


print()
print("=" * 120)
print("QC")
print("=" * 120)

print(
    qc.to_string(
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

