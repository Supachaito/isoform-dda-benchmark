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

STEP2A = (
    MASTER /
    "STEP2A_CELL_LINE_ANALYSIS"
)

STEP2B = (
    MASTER /
    "STEP2B_CROSS_SOFTWARE_CONCORDANCE"
)

SUPPORT_FILE = (
    STEP2A /
    "03_Peptide_ReplicateSupport.csv"
)

STEP2B_MEMBERSHIP = (
    STEP2B /
    "03_Membership_Robust_IsoformDiscriminativePeptides.csv"
)

FASTA = Path(
    r"<BENCHMARK_ROOT>"
    r"\Output_AP_cano_only\DB"
    r"\uniprotkb_proteome_UP000005640_2026_08_04.fasta"
)

OUT = (
    MASTER /
    "STEP2C_CONSENSUS_CANDIDATES"
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

SUBSET_CLASS = (
    "within_family_subset_discriminative"
)


# =====================================================================
# HELPERS
# =====================================================================

def il_norm(seq):

    return (
        str(seq)
        .upper()
        .replace("I", "J")
        .replace("L", "J")
    )


def is_isoform(acc):

    if pd.isna(acc):
        return False

    acc = str(
        acc
    ).strip()

    if "-" not in acc:
        return False

    base, suffix = acc.rsplit(
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


def parse_accession(header):

    h = header.strip()

    if h.startswith(">"):
        h = h[1:]

    first = h.split()[0]

    parts = first.split("|")

    if (
        len(parts) >= 2
        and
        parts[0].lower()
        in {
            "sp",
            "tr"
        }
    ):
        return parts[1].strip()

    return first.strip()


def read_fasta(path):

    rows = []

    header = None
    seq_parts = []


    def flush():

        nonlocal header
        nonlocal seq_parts

        if header is None:
            return

        accession = parse_accession(
            header
        )

        sequence = (
            "".join(
                seq_parts
            )
            .replace(
                " ",
                ""
            )
            .upper()
        )

        gene_match = re.search(
            r"(?:^|\s)GN=([^\s]+)",
            header
        )

        gene = (
            gene_match.group(1)
            if gene_match
            else ""
        )

        rows.append({
            "Accession":
                accession,

            "Gene":
                gene,

            "Sequence":
                sequence,

            "ProteinLength":
                len(
                    sequence
                )
        })


    with open(
        path,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                flush()

                header = line
                seq_parts = []

            else:

                seq_parts.append(
                    line
                )

        flush()


    return pd.DataFrame(
        rows
    )


def tier_order(value):

    lookup = {
        "Tier1":
            1,

        "Tier2":
            2,

        "Tier3":
            3,

        "Exploratory":
            4,

        "Not_CellLineSpecific":
            5
    }

    return lookup.get(
        value,
        99
    )


def strongest_tier(values):

    vals = [
        x
        for x in values
        if pd.notna(x)
    ]

    if not vals:
        return ""

    return sorted(
        vals,
        key=tier_order
    )[0]


# =====================================================================
# LOAD FASTA
# =====================================================================

print()
print("=" * 120)
print("STEP 2C — CONSENSUS ISOFORM CANDIDATE RANKING")
print("=" * 120)

print()
print("Loading FASTA annotation...")


fasta = read_fasta(
    FASTA
)


acc_to_gene = dict(
    zip(
        fasta[
            "Accession"
        ],
        fasta[
            "Gene"
        ]
    )
)


acc_to_seq = dict(
    zip(
        fasta[
            "Accession"
        ],
        fasta[
            "Sequence"
        ]
    )
)


acc_to_length = dict(
    zip(
        fasta[
            "Accession"
        ],
        fasta[
            "ProteinLength"
        ]
    )
)


print(
    "FASTA entries:",
    len(
        fasta
    )
)

print(
    "Gene annotated:",
    (
        fasta[
            "Gene"
        ]
        .astype(str)
        .str.len()
        .gt(0)
        .sum()
    )
)


# =====================================================================
# LOAD SUPPORT DATA
# =====================================================================

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
] = (
    support[
        "Program"
    ]
    .str.upper()
)


support[
    "Peptide"
] = (
    support[
        "Peptide"
    ]
    .str.upper()
)


support[
    "ReplicateSupport"
] = pd.to_numeric(
    support[
        "ReplicateSupport"
    ],
    errors="coerce"
).fillna(
    0
).astype(
    int
)


# =====================================================================
# PEPTIDE → FIXED MAPPING
# =====================================================================

mapping = (
    support[
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


mapping_conflicts = (
    support
    .groupby(
        "Peptide"
    )
    .agg(
        ClassCount=(
            "GeneAwareClass",
            "nunique"
        ),
        MappingCount=(
            "MappedAccessions",
            "nunique"
        )
    )
    .reset_index()
)


mapping_conflicts = mapping_conflicts[
    (
        mapping_conflicts[
            "ClassCount"
        ]
        >
        1
    )
    |
    (
        mapping_conflicts[
            "MappingCount"
        ]
        >
        1
    )
]


# =====================================================================
# COMPLETE SUPPORT CUBE
#
# peptide × program × cell line
# missing = 0
# =====================================================================

all_peptides = sorted(
    set(
        support[
            "Peptide"
        ]
    )
)


support_lookup = {

    (
        r[
            "Peptide"
        ],
        r[
            "Program"
        ],
        r[
            "CellLine"
        ]
    ):
        int(
            r[
                "ReplicateSupport"
            ]
        )

    for _, r
    in support.iterrows()
}


def get_support(
    peptide,
    program,
    cell_line
):

    return support_lookup.get(
        (
            peptide,
            program,
            cell_line
        ),
        0
    )


# =====================================================================
# ROBUST ISOFORM PEPTIDE UNION PER CELL LINE
# =====================================================================

robust_iso = support[
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
].copy()


candidate_pairs = (
    robust_iso[
        [
            "CellLine",
            "Peptide"
        ]
    ]
    .drop_duplicates()
)


# =====================================================================
# BUILD PEPTIDE CANDIDATE TABLE
# =====================================================================

candidate_rows = []


for _, candidate in (
    candidate_pairs
    .iterrows()
):

    cell_line = candidate[
        "CellLine"
    ]

    peptide = candidate[
        "Peptide"
    ]


    map_row = mapping[
        mapping[
            "Peptide"
        ]
        ==
        peptide
    ].iloc[0]


    evidence_class = map_row[
        "GeneAwareClass"
    ]


    mapped_accessions = split_accessions(
        map_row[
            "MappedAccessions"
        ]
    )


    isoform_accessions = [
        x
        for x in mapped_accessions
        if is_isoform(
            x
        )
    ]


    genes = sorted(
        {
            acc_to_gene.get(
                x,
                ""
            )
            for x in mapped_accessions
            if acc_to_gene.get(
                x,
                ""
            )
        }
    )


    target_support = {

        program:
            get_support(
                peptide,
                program,
                cell_line
            )

        for program
        in PROGRAMS
    }


    robust_programs = [
        program

        for program
        in PROGRAMS

        if target_support[
            program
        ]
        >=
        2
    ]


    full3_programs = [
        program

        for program
        in PROGRAMS

        if target_support[
            program
        ]
        ==
        3
    ]


    other_lines = [
        x
        for x in CELL_LINES
        if x != cell_line
    ]


    other_robust_hits = []

    other_any_hits = []


    for other_cell in other_lines:

        for program in PROGRAMS:

            s = get_support(
                peptide,
                program,
                other_cell
            )


            if s >= 2:

                other_robust_hits.append(
                    (
                        other_cell,
                        program,
                        s
                    )
                )


            if s >= 1:

                other_any_hits.append(
                    (
                        other_cell,
                        program,
                        s
                    )
                )


    # ---------------------------------------------------------------
    # Old strict definition at software level:
    # target >=2/3 AND exact 0/3 in each other cell line
    # ---------------------------------------------------------------

    strict_specific_programs = []


    for program in PROGRAMS:

        if target_support[
            program
        ] < 2:
            continue


        if all(
            get_support(
                peptide,
                program,
                other_cell
            )
            ==
            0

            for other_cell
            in other_lines
        ):

            strict_specific_programs.append(
                program
            )


    robust_software_count = len(
        robust_programs
    )


    full3_count = len(
        full3_programs
    )


    strict_count = len(
        strict_specific_programs
    )


    no_robust_outside = (
        len(
            other_robust_hits
        )
        ==
        0
    )


    # ---------------------------------------------------------------
    # CONSERVATIVE CELL-LINE PRIORITY TIERS
    # ---------------------------------------------------------------

    if not no_robust_outside:

        tier = (
            "Not_CellLineSpecific"
        )


    elif (
        evidence_class
        ==
        UNIQUE_CLASS
        and
        robust_software_count
        >=
        3
    ):

        tier = "Tier1"


    elif (
        (
            evidence_class
            ==
            SUBSET_CLASS
            and
            robust_software_count
            >=
            3
        )
        or
        (
            evidence_class
            ==
            UNIQUE_CLASS
            and
            robust_software_count
            >=
            2
        )
    ):

        tier = "Tier2"


    elif robust_software_count >= 2:

        tier = "Tier3"


    else:

        tier = (
            "Exploratory"
        )


    candidate_rows.append({

        "CellLine":
            cell_line,

        "Peptide":
            peptide,

        "EvidenceClass":
            evidence_class,

        "GeneSymbols":
            ";".join(
                genes
            ),

        "MappedAccessions":
            ";".join(
                mapped_accessions
            ),

        "IsoformAccessions":
            ";".join(
                isoform_accessions
            ),

        "MappedAccessionCount":
            len(
                mapped_accessions
            ),

        "IsoformAccessionCount":
            len(
                isoform_accessions
            ),

        "AP_Support":
            target_support[
                "AP"
            ],

        "FP_Support":
            target_support[
                "FP"
            ],

        "MM_Support":
            target_support[
                "MM"
            ],

        "MQ_Support":
            target_support[
                "MQ"
            ],

        "RobustSoftwareCount":
            robust_software_count,

        "RobustPrograms":
            ";".join(
                robust_programs
            ),

        "ThreeOfThreeSoftwareCount":
            full3_count,

        "ThreeOfThreePrograms":
            ";".join(
                full3_programs
            ),

        "StrictSpecificSoftwareCount":
            strict_count,

        "StrictSpecificPrograms":
            ";".join(
                strict_specific_programs
            ),

        "OtherCellRobustHitCount":
            len(
                other_robust_hits
            ),

        "OtherCellAnyHitCount":
            len(
                other_any_hits
            ),

        "NoRobustEvidenceInOtherCellLines":
            no_robust_outside,

        "PriorityTier":
            tier
    })


candidates = pd.DataFrame(
    candidate_rows
)


# =====================================================================
# SORT / RANK
# =====================================================================

candidates[
    "TierOrder"
] = candidates[
    "PriorityTier"
].map(
    tier_order
)


candidates[
    "EvidenceOrder"
] = np.where(
    candidates[
        "EvidenceClass"
    ]
    ==
    UNIQUE_CLASS,
    0,
    1
)


candidates = candidates.sort_values(
    [
        "CellLine",
        "TierOrder",
        "EvidenceOrder",
        "RobustSoftwareCount",
        "ThreeOfThreeSoftwareCount",
        "StrictSpecificSoftwareCount",
        "OtherCellAnyHitCount",
        "Peptide"
    ],
    ascending=[
        True,
        True,
        True,
        False,
        False,
        False,
        True,
        True
    ]
).reset_index(
    drop=True
)


candidates[
    "PriorityRankWithinCellLine"
] = (
    candidates
    .groupby(
        "CellLine"
    )
    .cumcount()
    +
    1
)


# =====================================================================
# CELL-LINE-SPECIFIC RANKED TABLE
# =====================================================================

specific_candidates = candidates[
    candidates[
        "PriorityTier"
    ]
    !=
    "Not_CellLineSpecific"
].copy()


specific_candidates[
    "SpecificRank"
] = (
    specific_candidates
    .groupby(
        "CellLine"
    )
    .cumcount()
    +
    1
)


# =====================================================================
# TIER SUMMARY
# =====================================================================

tier_rows = []


for cell_line in CELL_LINES:

    x = candidates[
        candidates[
            "CellLine"
        ]
        ==
        cell_line
    ]


    for tier in [
        "Tier1",
        "Tier2",
        "Tier3",
        "Exploratory",
        "Not_CellLineSpecific"
    ]:

        y = x[
            x[
                "PriorityTier"
            ]
            ==
            tier
        ]


        tier_rows.append({

            "CellLine":
                cell_line,

            "PriorityTier":
                tier,

            "PeptideCount":
                y[
                    "Peptide"
                ].nunique(),

            "SingleIsoformUnique":
                y.loc[
                    y[
                        "EvidenceClass"
                    ]
                    ==
                    UNIQUE_CLASS,
                    "Peptide"
                ].nunique(),

            "SubsetDiscriminative":
                y.loc[
                    y[
                        "EvidenceClass"
                    ]
                    ==
                    SUBSET_CLASS,
                    "Peptide"
                ].nunique()
        })


tier_summary = pd.DataFrame(
    tier_rows
)


# =====================================================================
# EXPAND PEPTIDE → ISOFORM ACCESSION
#
# IMPORTANT:
# subset-discriminative peptides implicate a SET of possible isoforms;
# they do NOT uniquely resolve each accession.
# =====================================================================

iso_rows = []


for _, r in candidates.iterrows():

    accessions = split_accessions(
        r[
            "IsoformAccessions"
        ]
    )


    for accession in accessions:

        gene = acc_to_gene.get(
            accession,
            ""
        )


        iso_rows.append({

            "CellLine":
                r[
                    "CellLine"
                ],

            "IsoformAccession":
                accession,

            "Gene":
                gene,

            "Peptide":
                r[
                    "Peptide"
                ],

            "EvidenceClass":
                r[
                    "EvidenceClass"
                ],

            "PriorityTier":
                r[
                    "PriorityTier"
                ],

            "RobustSoftwareCount":
                r[
                    "RobustSoftwareCount"
                ],

            "RobustPrograms":
                r[
                    "RobustPrograms"
                ],

            "ThreeOfThreeSoftwareCount":
                r[
                    "ThreeOfThreeSoftwareCount"
                ],

            "NoRobustEvidenceInOtherCellLines":
                r[
                    "NoRobustEvidenceInOtherCellLines"
                ]
        })


iso_detail = pd.DataFrame(
    iso_rows
)


# =====================================================================
# ISOFORM-LEVEL SUMMARY
# =====================================================================

iso_summary_rows = []


if not iso_detail.empty:

    for (
        cell_line,
        accession
    ), x in iso_detail.groupby(
        [
            "CellLine",
            "IsoformAccession"
        ]
    ):

        unique_peptides = set(
            x.loc[
                x[
                    "EvidenceClass"
                ]
                ==
                UNIQUE_CLASS,
                "Peptide"
            ]
        )


        subset_peptides = set(
            x.loc[
                x[
                    "EvidenceClass"
                ]
                ==
                SUBSET_CLASS,
                "Peptide"
            ]
        )


        all_peps = set(
            x[
                "Peptide"
            ]
        )


        programs = set()


        for value in x[
            "RobustPrograms"
        ]:

            for p in str(
                value
            ).split(";"):

                p = p.strip()

                if p:
                    programs.add(
                        p
                    )


        best_tier = strongest_tier(
            x[
                "PriorityTier"
            ]
        )


        resolution = (
            "UniquelyResolved"

            if len(
                unique_peptides
            )
            >
            0

            else
            "SubsetImplicated"
        )


        iso_summary_rows.append({

            "CellLine":
                cell_line,

            "Gene":
                acc_to_gene.get(
                    accession,
                    ""
                ),

            "IsoformAccession":
                accession,

            "ResolutionStatus":
                resolution,

            "BestPriorityTier":
                best_tier,

            "TotalDiscriminativePeptides":
                len(
                    all_peps
                ),

            "SingleIsoformUniquePeptides":
                len(
                    unique_peptides
                ),

            "SubsetDiscriminativePeptides":
                len(
                    subset_peptides
                ),

            "SupportingSoftwareCount":
                len(
                    programs
                ),

            "SupportingPrograms":
                ";".join(
                    sorted(
                        programs
                    )
                ),

            "MaximumPeptideSoftwareConsensus":
                int(
                    x[
                        "RobustSoftwareCount"
                    ].max()
                ),

            "MaximumThreeOfThreeConsensus":
                int(
                    x[
                        "ThreeOfThreeSoftwareCount"
                    ].max()
                ),

            "AllCandidatePeptides":
                ";".join(
                    sorted(
                        all_peps
                    )
                )
        })


iso_summary = pd.DataFrame(
    iso_summary_rows
)


if not iso_summary.empty:

    iso_summary[
        "TierOrder"
    ] = iso_summary[
        "BestPriorityTier"
    ].map(
        tier_order
    )


    iso_summary[
        "ResolutionOrder"
    ] = np.where(
        iso_summary[
            "ResolutionStatus"
        ]
        ==
        "UniquelyResolved",
        0,
        1
    )


    iso_summary = iso_summary.sort_values(
        [
            "CellLine",
            "TierOrder",
            "ResolutionOrder",
            "SupportingSoftwareCount",
            "TotalDiscriminativePeptides",
            "IsoformAccession"
        ],
        ascending=[
            True,
            True,
            True,
            False,
            False,
            True
        ]
    ).reset_index(
        drop=True
    )


    iso_summary[
        "IsoformRankWithinCellLine"
    ] = (
        iso_summary
        .groupby(
            "CellLine"
        )
        .cumcount()
        +
        1
    )


# =====================================================================
# PEPTIDE POSITION MAPPING
# FOR TIER 1–3 CANDIDATES
# =====================================================================

position_rows = []


position_candidates = candidates[
    candidates[
        "PriorityTier"
    ].isin(
        [
            "Tier1",
            "Tier2",
            "Tier3"
        ]
    )
]


for _, r in (
    position_candidates
    .iterrows()
):

    peptide = r[
        "Peptide"
    ]

    normalized_peptide = il_norm(
        peptide
    )


    for accession in split_accessions(
        r[
            "MappedAccessions"
        ]
    ):

        protein_sequence = acc_to_seq.get(
            accession,
            ""
        )


        if not protein_sequence:
            continue


        normalized_protein = il_norm(
            protein_sequence
        )


        starts = []

        pos = normalized_protein.find(
            normalized_peptide
        )


        while pos != -1:

            starts.append(
                pos
                +
                1
            )

            pos = normalized_protein.find(
                normalized_peptide,
                pos
                +
                1
            )


        for start in starts:

            end = (
                start
                +
                len(
                    peptide
                )
                -
                1
            )


            position_rows.append({

                "CellLine":
                    r[
                        "CellLine"
                    ],

                "PriorityTier":
                    r[
                        "PriorityTier"
                    ],

                "Peptide":
                    peptide,

                "EvidenceClass":
                    r[
                        "EvidenceClass"
                    ],

                "Gene":
                    acc_to_gene.get(
                        accession,
                        ""
                    ),

                "Accession":
                    accession,

                "IsIsoform":
                    is_isoform(
                        accession
                    ),

                "PeptideStart":
                    start,

                "PeptideEnd":
                    end,

                "ProteinLength":
                    acc_to_length.get(
                        accession,
                        ""
                    )
            })


position_df = pd.DataFrame(
    position_rows
)


# =====================================================================
# SOFTWARE SUPPORT MATRIX
# =====================================================================

matrix_cols = [
    "CellLine",
    "PriorityRankWithinCellLine",
    "PriorityTier",
    "GeneSymbols",
    "IsoformAccessions",
    "Peptide",
    "EvidenceClass",
    "AP_Support",
    "FP_Support",
    "MM_Support",
    "MQ_Support",
    "RobustSoftwareCount",
    "ThreeOfThreeSoftwareCount",
    "StrictSpecificSoftwareCount",
    "NoRobustEvidenceInOtherCellLines"
]


support_matrix = candidates[
    matrix_cols
].copy()


# =====================================================================
# TOP CANDIDATES
# =====================================================================

top_rows = []


for cell_line in CELL_LINES:

    x = specific_candidates[
        specific_candidates[
            "CellLine"
        ]
        ==
        cell_line
    ].head(
        15
    )


    top_rows.append(
        x
    )


top_candidates = pd.concat(
    top_rows,
    ignore_index=True
)


# =====================================================================
# STEP 2B CHECKSUM
# =====================================================================

checksum_rows = []


if STEP2B_MEMBERSHIP.exists():

    membership = pd.read_csv(
        STEP2B_MEMBERSHIP,
        dtype=str
    )


    for cell_line in CELL_LINES:

        expected = membership[
            membership[
                "CellLine"
            ]
            ==
            cell_line
        ][
            "Peptide"
        ].nunique()


        observed = candidates[
            candidates[
                "CellLine"
            ]
            ==
            cell_line
        ][
            "Peptide"
        ].nunique()


        checksum_rows.append({

            "CellLine":
                cell_line,

            "Step2B_RobustIsoformPeptideUnion":
                expected,

            "Step2C_RobustIsoformPeptideUnion":
                observed,

            "Difference":
                observed
                -
                expected
        })


checksum = pd.DataFrame(
    checksum_rows
)


checksum_failures = (
    int(
        (
            checksum[
                "Difference"
            ]
            !=
            0
        ).sum()
    )
    if not checksum.empty
    else 0
)


# =====================================================================
# QC
# =====================================================================

status = "PASS"


if len(
    mapping_conflicts
) != 0:

    status = "REVIEW"


if checksum_failures != 0:

    status = "REVIEW"


qc = pd.DataFrame(
    [
        {
            "Metric":
                "RobustIsoformCandidateRows",

            "Value":
                len(
                    candidates
                )
        },

        {
            "Metric":
                "MappingConflicts",

            "Value":
                len(
                    mapping_conflicts
                )
        },

        {
            "Metric":
                "Step2B_ChecksumFailures",

            "Value":
                checksum_failures
        },

        {
            "Metric":
                "Tier1Definition",

            "Value":
                "single-isoform unique; >=3 robust software; no >=2/3 evidence in other cell lines"
        },

        {
            "Metric":
                "Tier2Definition",

            "Value":
                "subset-discriminative >=3 software OR single-isoform unique >=2 software; no robust evidence outside target cell line"
        },

        {
            "Metric":
                "Tier3Definition",

            "Value":
                "isoform-discriminative >=2 software; no robust evidence outside target cell line"
        },

        {
            "Metric":
                "ExploratoryDefinition",

            "Value":
                "robust in one software only; no robust evidence outside target cell line"
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
# EXPORT
# =====================================================================

outputs = {

    "01_Ranked_All_RobustIsoformPeptideCandidates.csv":
        candidates,

    "02_Ranked_CellLineSpecific_PeptideCandidates.csv":
        specific_candidates,

    "03_CandidateTier_Summary.csv":
        tier_summary,

    "04_IsoformCandidate_EvidenceDetail.csv":
        iso_detail,

    "05_Ranked_IsoformCandidates.csv":
        iso_summary,

    "06_SoftwareSupport_Matrix.csv":
        support_matrix,

    "07_Top15_Candidates_PerCellLine.csv":
        top_candidates,

    "08_PeptidePositions_ForTier1to3.csv":
        position_df,

    "09_STEP2B_Checksum.csv":
        checksum,

    "10_STEP2C_QC.csv":
        qc,

    "11_MappingConflict_QC.csv":
        mapping_conflicts
}


for name, df in outputs.items():

    df.to_csv(
        OUT /
        name,
        index=False,
        encoding="utf-8-sig"
    )


# =====================================================================
# PRINT RESULTS
# =====================================================================

print()
print("=" * 125)
print("CANDIDATE TIER SUMMARY")
print("=" * 125)

print(
    tier_summary.to_string(
        index=False
    )
)


print()
print("=" * 125)
print("TOP CELL-LINE-ASSOCIATED PEPTIDE CANDIDATES")
print("=" * 125)


show_cols = [
    "CellLine",
    "SpecificRank",
    "PriorityTier",
    "GeneSymbols",
    "IsoformAccessions",
    "Peptide",
    "EvidenceClass",
    "AP_Support",
    "FP_Support",
    "MM_Support",
    "MQ_Support",
    "RobustSoftwareCount",
    "ThreeOfThreeSoftwareCount",
    "StrictSpecificSoftwareCount"
]


print(
    specific_candidates[
        show_cols
    ]
    .groupby(
        "CellLine",
        group_keys=False
    )
    .head(
        10
    )
    .to_string(
        index=False
    )
)


print()
print("=" * 125)
print("TOP ISOFORM-LEVEL CANDIDATES")
print("=" * 125)


if not iso_summary.empty:

    iso_show = [
        "CellLine",
        "IsoformRankWithinCellLine",
        "BestPriorityTier",
        "Gene",
        "IsoformAccession",
        "ResolutionStatus",
        "TotalDiscriminativePeptides",
        "SingleIsoformUniquePeptides",
        "SubsetDiscriminativePeptides",
        "SupportingSoftwareCount",
        "SupportingPrograms"
    ]


    print(
        iso_summary[
            iso_show
        ]
        .groupby(
            "CellLine",
            group_keys=False
        )
        .head(
            10
        )
        .to_string(
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
    # FIGURE 5A — TIER COUNTS
    # -----------------------------------------------------------------

    tiers_to_plot = [
        "Tier1",
        "Tier2",
        "Tier3",
        "Exploratory"
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

    width = 0.2


    for i, tier in enumerate(
        tiers_to_plot
    ):

        values = []

        for cell_line in CELL_LINES:

            row = tier_summary[
                (
                    tier_summary[
                        "CellLine"
                    ]
                    ==
                    cell_line
                )
                &
                (
                    tier_summary[
                        "PriorityTier"
                    ]
                    ==
                    tier
                )
            ]


            values.append(
                int(
                    row.iloc[0][
                        "PeptideCount"
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
            label=tier
        )


    ax.set_xticks(
        x
    )

    ax.set_xticklabels(
        CELL_LINES
    )

    ax.set_ylabel(
        "Cell-line-associated isoform-discriminative peptides"
    )

    ax.set_title(
        "Prioritization of reproducible isoform candidates"
    )

    ax.legend(
        frameon=False
    )

    fig.tight_layout()

    fig.savefig(
        FIGDIR /
        "Figure5A_CandidateTierCounts.png",
        dpi=300
    )

    fig.savefig(
        FIGDIR /
        "Figure5A_CandidateTierCounts.pdf"
    )

    plt.close(
        fig
    )


    # -----------------------------------------------------------------
    # FIGURE 5B — SOFTWARE × REPLICATE SUPPORT MATRIX
    # TOP SPECIFIC CANDIDATES
    # -----------------------------------------------------------------

    selected = []


    for cell_line in CELL_LINES:

        xdf = specific_candidates[
            specific_candidates[
                "CellLine"
            ]
            ==
            cell_line
        ].head(
            10
        )

        selected.append(
            xdf
        )


    selected = pd.concat(
        selected,
        ignore_index=True
    )


    if not selected.empty:

        matrix = selected[
            [
                "AP_Support",
                "FP_Support",
                "MM_Support",
                "MQ_Support"
            ]
        ].astype(
            float
        ).values


        labels = []


        for _, r in selected.iterrows():

            gene = (
                r[
                    "GeneSymbols"
                ]
                if r[
                    "GeneSymbols"
                ]
                else "NA"
            )


            acc = (
                r[
                    "IsoformAccessions"
                ]
                if r[
                    "IsoformAccessions"
                ]
                else "NA"
            )


            labels.append(
                f"{r['CellLine']} | {gene} | {acc} | {r['Peptide']}"
            )


        fig_height = max(
            7,
            0.35
            *
            len(
                selected
            )
        )


        fig, ax = plt.subplots(
            figsize=(
                8,
                fig_height
            )
        )


        image = ax.imshow(
            matrix,
            aspect="auto",
            vmin=0,
            vmax=3
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
                len(
                    labels
                )
            )
        )

        ax.set_yticklabels(
            labels,
            fontsize=7
        )


        for i in range(
            matrix.shape[0]
        ):

            for j in range(
                matrix.shape[1]
            ):

                ax.text(
                    j,
                    i,
                    str(
                        int(
                            matrix[
                                i,
                                j
                            ]
                        )
                    ),
                    ha="center",
                    va="center",
                    fontsize=8
                )


        ax.set_title(
            "Replicate support for prioritized isoform-discriminative peptides"
        )


        fig.colorbar(
            image,
            ax=ax,
            label="Replicates detected (0–3)"
        )


        fig.tight_layout()


        fig.savefig(
            FIGDIR /
            "Figure5B_TopCandidate_SoftwareReplicateMatrix.png",
            dpi=300,
            bbox_inches="tight"
        )


        fig.savefig(
            FIGDIR /
            "Figure5B_TopCandidate_SoftwareReplicateMatrix.pdf",
            bbox_inches="tight"
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
# STATUS
# =====================================================================

with open(
    OUT /
    "12_STEP2C_STATUS.txt",
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
        "STEP 2C — CONSENSUS ISOFORM CANDIDATE RANKING\n\n"
    )

    fh.write(
        f"FINAL STATUS: {status}\n"
    )

    fh.write(
        f"PLOT STATUS: {plot_status}\n\n"
    )

    fh.write(
        "Tier 1: single-isoform unique peptide, robust in >=3/4 software, "
        "with no robust >=2/3 evidence in another cell line.\n"
    )

    fh.write(
        "Tier 2: subset-discriminative peptide robust in >=3/4 software "
        "or single-isoform unique peptide robust in >=2/4 software, "
        "with no robust evidence in another cell line.\n"
    )

    fh.write(
        "Tier 3: isoform-discriminative peptide robust in >=2/4 software "
        "with no robust evidence in another cell line.\n"
    )

    fh.write(
        "Exploratory: robust in one software only with no robust evidence "
        "in another cell line.\n\n"
    )

    fh.write(
        "Subset-discriminative evidence implicates candidate isoforms "
        "but does not uniquely resolve an individual isoform accession.\n"
    )


print()
print("=" * 125)
print("STEP 2C QC")
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
print("OUTPUT:")
print(OUT)

print()
print("FIGURES:")
print(FIGDIR)

print()
print("STEP 2C COMPLETE")

