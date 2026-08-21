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
import re


# ============================================================
# STEP 4B
# SEPepQuant-inspired structural-equivalence quantification
#
# IMPORTANT:
# - This is NOT claimed to be the official SEPepQuant implementation.
# - Structural equivalence follows the same core definition:
#   peptides mapping to exactly the same protein-accession set
#   form one structural quantification unit.
# - Starting peptides have already passed Step1D target/FDR filtering.
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

STEP4A = (
    MASTER /
    "STEP4A_QUANT_BENCHMARK_FINAL"
)

MAPPING_FILE = (
    STEP1D /
    "03_FINAL_PeptideMappingDetail.csv"
)

QUANT_FILE = (
    STEP4A /
    "05_PeptideIntensity_Normalized_Long.csv"
)

FASTA_FILE = Path(
    r"<BENCHMARK_ROOT>"
    r"\Output_AP_cano_only\DB"
    r"\uniprotkb_proteome_UP000005640_2026_08_04.fasta"
)

OUT = (
    MASTER /
    "STEP4B_SEPEPQUANT_INSPIRED"
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


# ============================================================
# 2. HELPERS
# ============================================================

def split_accessions(value):

    if pd.isna(value):
        return []

    return sorted(
        {
            x.strip()
            for x in str(value).split(";")
            if x.strip()
        }
    )


def pearson_safe(a, b):

    if len(a) < 3:
        return np.nan

    if np.nanstd(a) == 0:
        return np.nan

    if np.nanstd(b) == 0:
        return np.nan

    return pd.Series(a).corr(
        pd.Series(b),
        method="pearson"
    )


def spearman_safe(a, b):

    if len(a) < 3:
        return np.nan

    return pd.Series(a).corr(
        pd.Series(b),
        method="spearman"
    )


# ============================================================
# 3. PARSE FASTA
# accession -> gene
# gene -> all protein entries in benchmark FASTA
# ============================================================

print()
print("=" * 120)
print("STEP 4B — SEPepQuant-inspired structural-equivalence quantification")
print("=" * 120)


if not FASTA_FILE.exists():
    raise FileNotFoundError(
        FASTA_FILE
    )


acc_to_gene = {}
gene_to_acc = {}


with open(
    FASTA_FILE,
    "r",
    encoding="utf-8",
    errors="replace"
) as fh:

    for line in fh:

        if not line.startswith(">"):
            continue

        header = line[1:].strip()

        first = header.split()[0]

        if "|" in first:

            parts = first.split("|")

            if len(parts) >= 2:
                accession = parts[1]

            else:
                accession = first

        else:
            accession = first


        m = re.search(
            r"\bGN=([^\s]+)",
            header
        )


        if m:

            gene = m.group(1).strip()

            acc_to_gene[
                accession
            ] = gene

            gene_to_acc.setdefault(
                gene,
                set()
            ).add(
                accession
            )


print(
    f"FASTA accessions with GN: "
    f"{len(acc_to_gene):,}"
)

print(
    f"FASTA genes: "
    f"{len(gene_to_acc):,}"
)


# ============================================================
# 4. LOAD FINAL PEPTIDE MAPPING
# ============================================================

mapping = pd.read_csv(
    MAPPING_FILE,
    dtype=str,
    low_memory=False
)


required = {
    "Program",
    "Peptide",
    "GeneAwareClass",
    "MappedAccessions"
}


missing = (
    required
    -
    set(mapping.columns)
)


if missing:

    raise RuntimeError(
        f"Missing mapping columns: {missing}"
    )


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
# 5. VERIFY THAT PEPTIDE -> ACCESSION SET IS STABLE
# ============================================================

mapping_check = (
    mapping
    .groupby(
        "Peptide"
    )
    .agg(
        NMapping=(
            "MappedAccessions",
            "nunique"
        )
    )
    .reset_index()
)


conflicts = mapping_check[
    mapping_check[
        "NMapping"
    ] > 1
]


if len(conflicts) > 0:

    conflicts.to_csv(
        OUT /
        "00_PeptideMappingConflict.csv",
        index=False
    )

    raise RuntimeError(
        f"Peptide mapping conflicts: "
        f"{len(conflicts)}"
    )


peptide_map = (
    mapping[
        [
            "Peptide",
            "MappedAccessions"
        ]
    ]
    .drop_duplicates(
        subset=[
            "Peptide"
        ]
    )
    .copy()
)


# ============================================================
# 6. STRUCTURAL-EQUIVALENCE DEFINITION
#
# Same exact accession set = same SEPEP-like unit
# ============================================================

def annotate_structural_group(
    value
):

    accessions = split_accessions(
        value
    )


    if len(accessions) == 0:

        return pd.Series({
            "StructuralKey": "",
            "MappedGenes": "",
            "NAccessions": 0,
            "NGenes": 0,
            "Gene": "",
            "GeneProteinEntries": np.nan,
            "SEPEPClass": "CU"
        })


    genes = sorted(
        {
            acc_to_gene[a]
            for a in accessions
            if a in acc_to_gene
        }
    )


    unresolved = [
        a
        for a in accessions
        if a not in acc_to_gene
    ]


    key = "|".join(
        accessions
    )


    if unresolved:

        return pd.Series({
            "StructuralKey":
                key,

            "MappedGenes":
                "|".join(
                    genes
                ),

            "NAccessions":
                len(accessions),

            "NGenes":
                len(genes),

            "Gene":
                genes[0]
                if len(genes) == 1
                else "",

            "GeneProteinEntries":
                np.nan,

            "SEPEPClass":
                "CU"
        })


    # --------------------------------------------------------
    # C5 = peptide group maps across genes
    # --------------------------------------------------------

    if len(genes) > 1:

        return pd.Series({
            "StructuralKey":
                key,

            "MappedGenes":
                "|".join(
                    genes
                ),

            "NAccessions":
                len(accessions),

            "NGenes":
                len(genes),

            "Gene":
                "",

            "GeneProteinEntries":
                np.nan,

            "SEPEPClass":
                "C5"
        })


    gene = genes[0]

    all_gene_acc = gene_to_acc[
        gene
    ]

    n_gene_entries = len(
        all_gene_acc
    )


    # --------------------------------------------------------
    # C1 = gene represented by one protein entry
    # --------------------------------------------------------

    if n_gene_entries == 1:

        cls = "C1"


    # --------------------------------------------------------
    # C2 = multi-isoform gene,
    #      group resolves one exact protein entry
    # --------------------------------------------------------

    elif len(accessions) == 1:

        cls = "C2"


    # --------------------------------------------------------
    # C4 = maps to every protein entry of the gene
    # --------------------------------------------------------

    elif set(accessions) == all_gene_acc:

        cls = "C4"


    # --------------------------------------------------------
    # C3 = maps to subset of protein entries
    # --------------------------------------------------------

    else:

        cls = "C3"


    return pd.Series({
        "StructuralKey":
            key,

        "MappedGenes":
            gene,

        "NAccessions":
            len(accessions),

        "NGenes":
            1,

        "Gene":
            gene,

        "GeneProteinEntries":
            n_gene_entries,

        "SEPEPClass":
            cls
    })


annotations = peptide_map[
    "MappedAccessions"
].apply(
    annotate_structural_group
)


peptide_map = pd.concat(
    [
        peptide_map.reset_index(
            drop=True
        ),
        annotations.reset_index(
            drop=True
        )
    ],
    axis=1
)


# ============================================================
# 7. CREATE DETERMINISTIC SEPEP IDS
# ============================================================

group_meta = (
    peptide_map[
        [
            "StructuralKey",
            "MappedGenes",
            "Gene",
            "NAccessions",
            "NGenes",
            "GeneProteinEntries",
            "SEPEPClass"
        ]
    ]
    .drop_duplicates(
        subset=[
            "StructuralKey"
        ]
    )
    .copy()
)


group_meta[
    "SEPEPID"
] = ""


# C1-C4: order within gene

single_gene = group_meta[
    group_meta[
        "SEPEPClass"
    ].isin(
        [
            "C1",
            "C2",
            "C3",
            "C4"
        ]
    )
].copy()


single_gene = single_gene.sort_values(
    [
        "Gene",
        "SEPEPClass",
        "StructuralKey"
    ]
)


single_gene[
    "WithinGeneOrder"
] = (
    single_gene
    .groupby(
        "Gene"
    )
    .cumcount()
    +
    1
)


single_gene[
    "SEPEPID"
] = (
    single_gene[
        "Gene"
    ]
    +
    "_SEPEP."
    +
    single_gene[
        "WithinGeneOrder"
    ].astype(str)
    +
    "_"
    +
    single_gene[
        "SEPEPClass"
    ]
)


group_meta.loc[
    single_gene.index,
    "SEPEPID"
] = single_gene[
    "SEPEPID"
]


# C5

multi = group_meta[
    group_meta[
        "SEPEPClass"
    ]
    ==
    "C5"
].sort_values(
    "StructuralKey"
)


for i, idx in enumerate(
    multi.index,
    start=1
):

    group_meta.loc[
        idx,
        "SEPEPID"
    ] = (
        f"Multiple_SEPEP.{i}_C5"
    )


# unresolved

unresolved = group_meta[
    group_meta[
        "SEPEPClass"
    ]
    ==
    "CU"
].sort_values(
    "StructuralKey"
)


for i, idx in enumerate(
    unresolved.index,
    start=1
):

    group_meta.loc[
        idx,
        "SEPEPID"
    ] = (
        f"Unresolved_SEPEP.{i}_CU"
    )


peptide_map = peptide_map.merge(
    group_meta[
        [
            "StructuralKey",
            "SEPEPID"
        ]
    ],
    on="StructuralKey",
    how="left"
)


# ============================================================
# 8. NUMBER OF PEPTIDES PER STRUCTURAL GROUP
# ============================================================

peptide_counts = (
    peptide_map
    .groupby(
        "SEPEPID",
        as_index=False
    )
    .agg(
        StructuralPeptides=(
            "Peptide",
            "nunique"
        )
    )
)


group_meta = group_meta.merge(
    peptide_counts,
    on="SEPEPID",
    how="left"
)


# ============================================================
# 9. LOAD STEP4A QUANTITATIVE DATA
# ============================================================

quant = pd.read_csv(
    QUANT_FILE,
    low_memory=False
)


required_quant = {
    "Program",
    "Sample",
    "Peptide",
    "NormLog2"
}


missing_quant = (
    required_quant
    -
    set(quant.columns)
)


if missing_quant:

    raise RuntimeError(
        "Quant table missing: "
        +
        str(missing_quant)
    )


quant[
    "Program"
] = (
    quant[
        "Program"
    ]
    .astype(str)
    .str.upper()
    .str.strip()
)


quant[
    "Peptide"
] = (
    quant[
        "Peptide"
    ]
    .astype(str)
    .str.upper()
    .str.strip()
)


quant[
    "NormLog2"
] = pd.to_numeric(
    quant[
        "NormLog2"
    ],
    errors="coerce"
)


quant = quant.dropna(
    subset=[
        "NormLog2"
    ]
)


# ============================================================
# 10. ADD STRUCTURAL GROUP
# ============================================================

quant = quant.merge(
    peptide_map[
        [
            "Peptide",
            "SEPEPID",
            "StructuralKey",
            "MappedGenes",
            "Gene",
            "SEPEPClass"
        ]
    ],
    on="Peptide",
    how="left"
)


missing_sepep = quant[
    quant[
        "SEPEPID"
    ].isna()
]


if len(missing_sepep) > 0:

    missing_sepep.to_csv(
        OUT /
        "01_Missing_SEPEP_mapping.csv",
        index=False
    )

    raise RuntimeError(
        f"Quant observations without SEPEP: "
        f"{len(missing_sepep)}"
    )


quant[
    "CellLine"
] = quant[
    "Sample"
].map(
    CELL_LINE
)


# ============================================================
# 11. SEPEP ABUNDANCE
#
# median of normalized peptide abundance
# ============================================================

sepep_quant = (
    quant
    .groupby(
        [
            "Program",
            "Sample",
            "CellLine",
            "SEPEPID",
            "StructuralKey",
            "MappedGenes",
            "Gene",
            "SEPEPClass"
        ],
        dropna=False,
        as_index=False
    )
    .agg(
        SEPEP_Log2Abundance=(
            "NormLog2",
            "median"
        ),

        PeptidesQuantified=(
            "Peptide",
            "nunique"
        )
    )
)


# ============================================================
# 12. GENE-LEVEL ABUNDANCE
#
# Conservative gene reference:
# use only single-gene peptides (C1-C4)
# exclude C5 multi-gene and CU unresolved
# ============================================================

gene_source = quant[
    quant[
        "SEPEPClass"
    ].isin(
        [
            "C1",
            "C2",
            "C3",
            "C4"
        ]
    )
].copy()


gene_quant = (
    gene_source
    .groupby(
        [
            "Program",
            "Sample",
            "CellLine",
            "Gene"
        ],
        as_index=False
    )
    .agg(
        Gene_Log2Abundance=(
            "NormLog2",
            "median"
        ),

        GenePeptidesQuantified=(
            "Peptide",
            "nunique"
        ),

        GeneSEPEPsQuantified=(
            "SEPEPID",
            "nunique"
        )
    )
)


# ============================================================
# 13. MULTIPLE QUANTIFICATION UNITS PER GENE
# ============================================================

multiunit = (
    sepep_quant[
        sepep_quant[
            "SEPEPClass"
        ].isin(
            [
                "C1",
                "C2",
                "C3",
                "C4"
            ]
        )
    ]
    .groupby(
        [
            "Program",
            "Gene"
        ],
        as_index=False
    )
    .agg(
        QuantifiedSEPEPs=(
            "SEPEPID",
            "nunique"
        )
    )
)


multiunit_summary = (
    multiunit
    .groupby(
        "Program",
        as_index=False
    )
    .agg(
        GenesQuantified=(
            "Gene",
            "nunique"
        ),

        GenesWith2plusSEPEPs=(
            "QuantifiedSEPEPs",
            lambda x:
                int(
                    (
                        x >= 2
                    ).sum()
                )
        ),

        GenesWith3plusSEPEPs=(
            "QuantifiedSEPEPs",
            lambda x:
                int(
                    (
                        x >= 3
                    ).sum()
                )
        )
    )
)


# ============================================================
# 14. CROSS-SOFTWARE CORRELATION
# peptide -> SEPEP -> gene
# ============================================================

def pairwise_sample_correlations(
    data,
    entity,
    value,
    level
):

    rows = []


    for sample in SAMPLES:

        sx = data[
            data[
                "Sample"
            ]
            ==
            sample
        ]


        for p1, p2 in combinations(
            PROGRAMS,
            2
        ):

            a = (
                sx[
                    sx[
                        "Program"
                    ]
                    ==
                    p1
                ][
                    [
                        entity,
                        value
                    ]
                ]
                .rename(
                    columns={
                        value:
                            "A"
                    }
                )
            )


            b = (
                sx[
                    sx[
                        "Program"
                    ]
                    ==
                    p2
                ][
                    [
                        entity,
                        value
                    ]
                ]
                .rename(
                    columns={
                        value:
                            "B"
                    }
                )
            )


            m = (
                a
                .merge(
                    b,
                    on=entity,
                    how="inner"
                )
                .dropna()
            )


            rows.append({

                "Level":
                    level,

                "Sample":
                    sample,

                "CellLine":
                    CELL_LINE[
                        sample
                    ],

                "Program1":
                    p1,

                "Program2":
                    p2,

                "CommonEntities":
                    len(m),

                "PearsonR":
                    pearson_safe(
                        m["A"],
                        m["B"]
                    ),

                "SpearmanR":
                    spearman_safe(
                        m["A"],
                        m["B"]
                    )
            })


    return pd.DataFrame(
        rows
    )


pep_cor = pairwise_sample_correlations(
    quant,
    "Peptide",
    "NormLog2",
    "Peptide"
)


sepep_cor = pairwise_sample_correlations(
    sepep_quant,
    "SEPEPID",
    "SEPEP_Log2Abundance",
    "SEPEP"
)


gene_cor = pairwise_sample_correlations(
    gene_quant,
    "Gene",
    "Gene_Log2Abundance",
    "Gene"
)


level_cor = pd.concat(
    [
        pep_cor,
        sepep_cor,
        gene_cor
    ],
    ignore_index=True
)


level_cor_summary = (
    level_cor
    .groupby(
        [
            "Level",
            "Program1",
            "Program2"
        ],
        as_index=False
    )
    .agg(
        MedianCommonEntities=(
            "CommonEntities",
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
# 15. CELL-LINE FOLD CHANGE FUNCTION
# ============================================================

CONTRASTS = [
    (
        "C33A",
        "HeLa",
        "C33A_vs_HeLa"
    ),
    (
        "SiHa",
        "HeLa",
        "SiHa_vs_HeLa"
    ),
    (
        "SiHa",
        "C33A",
        "SiHa_vs_C33A"
    )
]


def make_fold_changes(
    data,
    entity_columns,
    value_column
):

    means = (
        data
        .groupby(
            [
                "Program"
            ]
            +
            entity_columns
            +
            [
                "CellLine"
            ],
            as_index=False
        )
        .agg(
            MeanValue=(
                value_column,
                "mean"
            ),

            NReplicates=(
                value_column,
                "count"
            )
        )
    )


    rows = []


    group_columns = (
        [
            "Program"
        ]
        +
        entity_columns
    )


    for keys, x in means.groupby(
        group_columns,
        dropna=False
    ):

        if not isinstance(
            keys,
            tuple
        ):

            keys = (
                keys,
            )


        meta = dict(
            zip(
                group_columns,
                keys
            )
        )


        lookup = {
            row.CellLine: (
                row.MeanValue,
                row.NReplicates
            )
            for row in x.itertuples()
        }


        for a, b, label in CONTRASTS:

            if (
                a not in lookup
                or
                b not in lookup
            ):

                continue


            mean_a, n_a = lookup[a]
            mean_b, n_b = lookup[b]


            if (
                n_a < 2
                or
                n_b < 2
            ):

                continue


            row = meta.copy()

            row.update({

                "Contrast":
                    label,

                "MeanA":
                    mean_a,

                "MeanB":
                    mean_b,

                "N_A":
                    n_a,

                "N_B":
                    n_b,

                "log2FC":
                    mean_a
                    -
                    mean_b
            })


            rows.append(
                row
            )


    return pd.DataFrame(
        rows
    )


# ============================================================
# 16. SEPEP AND GENE FOLD CHANGES
# ============================================================

single_gene_sepep = sepep_quant[
    sepep_quant[
        "SEPEPClass"
    ].isin(
        [
            "C1",
            "C2",
            "C3",
            "C4"
        ]
    )
].copy()


sepep_fc = make_fold_changes(
    single_gene_sepep,
    [
        "SEPEPID",
        "Gene",
        "SEPEPClass"
    ],
    "SEPEP_Log2Abundance"
)


gene_fc = make_fold_changes(
    gene_quant,
    [
        "Gene"
    ],
    "Gene_Log2Abundance"
)


gene_fc = gene_fc.rename(
    columns={
        "log2FC":
            "Gene_log2FC"
    }
)


sepep_fc = sepep_fc.rename(
    columns={
        "log2FC":
            "SEPEP_log2FC"
    }
)


# ============================================================
# 17. GENE vs SEPEP
# ============================================================

gene_vs_sepep = sepep_fc.merge(
    gene_fc[
        [
            "Program",
            "Gene",
            "Contrast",
            "Gene_log2FC"
        ]
    ],
    on=[
        "Program",
        "Gene",
        "Contrast"
    ],
    how="inner"
)


gene_vs_sepep[
    "Delta_SEPEP_minus_Gene"
] = (
    gene_vs_sepep[
        "SEPEP_log2FC"
    ]
    -
    gene_vs_sepep[
        "Gene_log2FC"
    ]
)


gene_vs_sepep[
    "AbsDelta"
] = np.abs(
    gene_vs_sepep[
        "Delta_SEPEP_minus_Gene"
    ]
)


gene_vs_sepep[
    "OppositeDirection"
] = (
    np.sign(
        gene_vs_sepep[
            "SEPEP_log2FC"
        ]
    )
    !=
    np.sign(
        gene_vs_sepep[
            "Gene_log2FC"
        ]
    )
)


# ============================================================
# 18. CONSENSUS DIVERGENT SEPEPS
#
# same structural group supported by >=2 software
# ============================================================

consensus = (
    gene_vs_sepep
    .groupby(
        [
            "SEPEPID",
            "Gene",
            "SEPEPClass",
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

        MedianGene_log2FC=(
            "Gene_log2FC",
            "median"
        ),

        MedianSEPEP_log2FC=(
            "SEPEP_log2FC",
            "median"
        ),

        MedianDelta=(
            "Delta_SEPEP_minus_Gene",
            "median"
        ),

        OppositeDirectionSoftwareN=(
            "OppositeDirection",
            "sum"
        )
    )
)


consensus[
    "AbsMedianDelta"
] = np.abs(
    consensus[
        "MedianDelta"
    ]
)


consensus[
    "DirectionAgreement"
] = (
    consensus[
        "OppositeDirectionSoftwareN"
    ]
    /
    consensus[
        "SoftwareN"
    ]
)


consensus_ranked = (
    consensus[
        consensus[
            "SoftwareN"
        ]
        >=
        2
    ]
    .sort_values(
        [
            "AbsMedianDelta",
            "SoftwareN"
        ],
        ascending=[
            False,
            False
        ]
    )
)


# ============================================================
# 19. CLASS SUMMARY BY PROGRAM
# ============================================================

class_summary = (
    sepep_quant
    .groupby(
        [
            "Program",
            "SEPEPClass"
        ],
        as_index=False
    )
    .agg(
        QuantifiedSEPEPs=(
            "SEPEPID",
            "nunique"
        )
    )
)


# ============================================================
# 20. EXPORT
# ============================================================

group_meta.to_csv(
    OUT /
    "01_SEPEP_StructuralGroup_Metadata.csv",
    index=False,
    encoding="utf-8-sig"
)


peptide_map.to_csv(
    OUT /
    "02_Peptide_to_SEPEP_Mapping.csv",
    index=False,
    encoding="utf-8-sig"
)


sepep_quant.to_csv(
    OUT /
    "03_SEPEP_Quantification_Long.csv",
    index=False,
    encoding="utf-8-sig"
)


gene_quant.to_csv(
    OUT /
    "04_Gene_Quantification_Long.csv",
    index=False,
    encoding="utf-8-sig"
)


multiunit.to_csv(
    OUT /
    "05_MultiSEPEP_Genes.csv",
    index=False,
    encoding="utf-8-sig"
)


multiunit_summary.to_csv(
    OUT /
    "06_MultiSEPEP_Gene_Summary.csv",
    index=False,
    encoding="utf-8-sig"
)


level_cor.to_csv(
    OUT /
    "07_CrossSoftware_LevelCorrelation.csv",
    index=False,
    encoding="utf-8-sig"
)


level_cor_summary.to_csv(
    OUT /
    "08_CrossSoftware_LevelCorrelation_Summary.csv",
    index=False,
    encoding="utf-8-sig"
)


gene_vs_sepep.to_csv(
    OUT /
    "09_Gene_vs_SEPEP_FoldChange.csv",
    index=False,
    encoding="utf-8-sig"
)


consensus_ranked.to_csv(
    OUT /
    "10_Consensus_Divergent_SEPEPs.csv",
    index=False,
    encoding="utf-8-sig"
)


class_summary.to_csv(
    OUT /
    "11_SEPEP_Class_Summary.csv",
    index=False,
    encoding="utf-8-sig"
)


# ============================================================
# 21. QC
# ============================================================

program_sample_pairs = (
    sepep_quant[
        [
            "Program",
            "Sample"
        ]
    ]
    .drop_duplicates()
)


qc = pd.DataFrame({

    "Check": [
        "Programs",
        "ProgramSamplePairs",
        "QuantWithoutSEPEP",
        "PeptideMappingConflicts",
        "UnresolvedStructuralGroups"
    ],

    "Value": [
        sepep_quant[
            "Program"
        ].nunique(),

        len(
            program_sample_pairs
        ),

        len(
            missing_sepep
        ),

        len(
            conflicts
        ),

        int(
            (
                group_meta[
                    "SEPEPClass"
                ]
                ==
                "CU"
            ).sum()
        )
    ],

    "Expected": [
        4,
        36,
        0,
        0,
        "report-only"
    ]
})


qc.to_csv(
    OUT /
    "12_STEP4B_QC.csv",
    index=False,
    encoding="utf-8-sig"
)


status = (
    "PASS"
    if (
        sepep_quant[
            "Program"
        ].nunique()
        ==
        4
        and
        len(
            program_sample_pairs
        )
        ==
        36
        and
        len(
            missing_sepep
        )
        ==
        0
        and
        len(
            conflicts
        )
        ==
        0
    )
    else
    "FAIL"
)


# ============================================================
# 22. PRINT RESULTS
# ============================================================

print()
print("=" * 120)
print("SEPEP CLASS SUMMARY")
print("=" * 120)

print(
    class_summary.to_string(
        index=False
    )
)


print()
print("=" * 120)
print("GENES WITH MULTIPLE QUANTIFICATION UNITS")
print("=" * 120)

print(
    multiunit_summary.to_string(
        index=False
    )
)


print()
print("=" * 120)
print("PEPTIDE -> SEPEP -> GENE CROSS-SOFTWARE CONCORDANCE")
print("=" * 120)

print(
    level_cor_summary.to_string(
        index=False
    )
)


print()
print("=" * 120)
print("TOP CONSENSUS SEPEP-vs-GENE DIVERGENCE")
print("=" * 120)

print(
    consensus_ranked.head(
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

