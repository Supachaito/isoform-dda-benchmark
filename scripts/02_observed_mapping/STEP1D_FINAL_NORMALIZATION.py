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
import h5py
import re


# ============================================================
# PATHS
# ============================================================

ROOT = _public_project_root()

MASTER = ROOT / "MASTER_TABLES_FINAL"

AP_DIR = ROOT / "AP_MBR_OFF"

PRIMARY_EVIDENCE = (
    MASTER /
    "01_PRIMARY_ID_OFF_Evidence.csv"
)

PRIMARY_UNIQUE = (
    MASTER /
    "02_PRIMARY_ID_OFF_UniquePeptides.csv"
)

STEP1B = (
    MASTER /
    "STEP1B_MAPPING_QC"
)

GENE_MAP_FILE = (
    STEP1B /
    "06_GeneAware_PeptideMapping.csv"
)

OUT = (
    MASTER /
    "STEP1D_FINAL_NORMALIZED"
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
    "C33A_1", "C33A_2", "C33A_3",
    "HELA_1", "HELA_2", "HELA_3",
    "SIHA_1", "SIHA_2", "SIHA_3"
]

ISO_RE = re.compile(
    r"^(.+?)-([1-9][0-9]*)$"
)

AA_RE = re.compile(
    r"^[A-Z]+$"
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


def is_isoform(acc):

    return (
        ISO_RE.fullmatch(
            str(acc)
        )
        is not None
    )


def sample_from_filename(path):

    name = path.name.replace(
        ".ms_data.hdf",
        ""
    )

    for sample in SAMPLES:

        if sample.upper() in name.upper():
            return sample

    return name


def target_to_bool(values):

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

            continue

        if isinstance(
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

            continue

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

    return np.array(
        out,
        dtype=bool
    )


def accession_union(series):

    result = set()

    for value in series:

        if pd.isna(value):
            continue

        for acc in str(
            value
        ).split(";"):

            acc = acc.strip()

            if acc:
                result.add(
                    acc
                )

    return result


# ============================================================
# LOAD EXISTING NORMALIZED EVIDENCE
# ============================================================

print()
print("=" * 110)
print("STEP 1D — FINAL PRIMARY NORMALIZATION")
print("=" * 110)

print()
print("Loading existing master evidence...")


evidence = pd.read_csv(
    PRIMARY_EVIDENCE,
    dtype=str,
    low_memory=False
)


for col in [
    "Program",
    "Sample",
    "Peptide"
]:

    if col not in evidence.columns:

        raise RuntimeError(
            f"Missing required evidence column: {col}"
        )


evidence["Program"] = (
    evidence["Program"]
    .fillna("")
    .str.strip()
)

evidence["Sample"] = (
    evidence["Sample"]
    .fillna("")
    .str.strip()
)

evidence["Peptide"] = (
    evidence["Peptide"]
    .fillna("")
    .astype(str)
    .str.strip()
    .str.upper()
)


# ============================================================
# 1. RE-EXTRACT ALPHAPEPT CORRECTLY
#    TARGET + QVALUE <= 0.01
# ============================================================

print()
print("=" * 110)
print("1. ALPHAPEPT FINAL FILTER: TARGET + q <= 0.01")
print("=" * 110)


ap_files = sorted(
    AP_DIR.rglob(
        "*.ms_data.hdf"
    )
)


if len(ap_files) == 0:

    raise RuntimeError(
        "No AlphaPept *.ms_data.hdf files found"
    )


ap_rows = []
ap_qc_rows = []


for file in ap_files:

    sample = sample_from_filename(
        file
    )


    with h5py.File(
        file,
        "r"
    ) as h5:

        required = [
            "peptide_fdr/sequence_naked",
            "peptide_fdr/q_value",
            "peptide_fdr/target"
        ]


        missing = [
            x
            for x in required
            if x not in h5
        ]


        if missing:

            raise RuntimeError(
                f"{file.name}: missing {missing}"
            )


        peptides = np.array(
            [
                dec(x)
                .strip()
                .upper()

                for x in h5[
                    "peptide_fdr/sequence_naked"
                ][:]
            ],
            dtype=object
        )


        qvals = np.asarray(
            h5[
                "peptide_fdr/q_value"
            ][:]
        ).astype(float)


        targets = target_to_bool(
            h5[
                "peptide_fdr/target"
            ][:]
        )


    if not (
        len(peptides)
        ==
        len(qvals)
        ==
        len(targets)
    ):

        raise RuntimeError(
            f"{file.name}: AlphaPept array lengths differ"
        )


    raw_n = len(
        peptides
    )

    target_n = int(
        targets.sum()
    )


    keep = (
        targets
        &
        np.isfinite(
            qvals
        )
        &
        (
            qvals <= 0.01
        )
    )


    kept_n = int(
        keep.sum()
    )


    kept_peptides = peptides[
        keep
    ]


    valid_aa = np.array(
        [
            bool(
                AA_RE.fullmatch(
                    x
                )
            )
            for x in kept_peptides
        ]
    )


    invalid_after_filter = int(
        (
            ~valid_aa
        ).sum()
    )


    for pep, q in zip(
        peptides[keep],
        qvals[keep]
    ):

        ap_rows.append({
            "Program":
                "AP",

            "Sample":
                sample,

            "Peptide":
                pep,

            "QValue":
                float(q)
        })


    ap_qc_rows.append({
        "Sample":
            sample,

        "RawPeptideFDRRows":
            raw_n,

        "TargetRows":
            target_n,

        "Target_Q01_Rows":
            kept_n,

        "Distinct_Target_Q01_Peptides":
            len(
                set(
                    kept_peptides
                )
            ),

        "InvalidSequenceRows_AfterFilter":
            invalid_after_filter
    })


ap_evidence = pd.DataFrame(
    ap_rows
)


ap_qc = pd.DataFrame(
    ap_qc_rows
)


print(
    ap_qc.to_string(
        index=False
    )
)


print()
print(
    "AP raw peptide_fdr rows       :",
    ap_qc[
        "RawPeptideFDRRows"
    ].sum()
)

print(
    "AP final target q<=0.01 rows  :",
    ap_qc[
        "Target_Q01_Rows"
    ].sum()
)

print(
    "AP final distinct peptides    :",
    ap_evidence[
        "Peptide"
    ].nunique()
)


# ============================================================
# 2. FP — KEEP EXISTING FILTERED PRIMARY SET
# ============================================================

fp_evidence = (
    evidence[
        evidence["Program"]
        ==
        "FP"
    ][
        [
            "Program",
            "Sample",
            "Peptide"
        ]
    ]
    .copy()
)


# ============================================================
# 3. MM — KEEP TARGET 1%-FDR SET,
#    REMOVE PIPE-DELIMITED AMBIGUOUS BASE SEQUENCES
# ============================================================

print()
print("=" * 110)
print("2. METAMORPHEUS FINAL UNAMBIGUOUS SEQUENCE FILTER")
print("=" * 110)


mm_raw = (
    evidence[
        evidence["Program"]
        ==
        "MM"
    ][
        [
            "Program",
            "Sample",
            "Peptide"
        ]
    ]
    .copy()
)


mm_pipe_mask = (
    mm_raw[
        "Peptide"
    ]
    .str.contains(
        r"\|",
        regex=True,
        na=False
    )
)


mm_invalid_mask = (
    ~mm_raw[
        "Peptide"
    ].str.fullmatch(
        r"[A-Z]+",
        na=False
    )
)


mm_ambiguous_rows = mm_raw[
    mm_pipe_mask
].copy()


mm_final = mm_raw[
    ~mm_pipe_mask
    &
    ~mm_invalid_mask
].copy()


mm_qc = pd.DataFrame(
    [
        {
            "Metric":
                "Input_Target_Q01_PSMRows",

            "Value":
                len(mm_raw)
        },
        {
            "Metric":
                "PipeDelimited_Ambiguous_PSMRows",

            "Value":
                int(
                    mm_pipe_mask.sum()
                )
        },
        {
            "Metric":
                "Distinct_PipeDelimited_Strings",

            "Value":
                mm_ambiguous_rows[
                    "Peptide"
                ].nunique()
        },
        {
            "Metric":
                "Final_Unambiguous_PSMRows",

            "Value":
                len(mm_final)
        },
        {
            "Metric":
                "Final_Unambiguous_DistinctPeptides",

            "Value":
                mm_final[
                    "Peptide"
                ].nunique()
        }
    ]
)


print(
    mm_qc.to_string(
        index=False
    )
)


# ============================================================
# 4. MQ — KEEP EXISTING TARGET/NON-CONTAMINANT SET
# ============================================================

mq_evidence = (
    evidence[
        evidence["Program"]
        ==
        "MQ"
    ][
        [
            "Program",
            "Sample",
            "Peptide"
        ]
    ]
    .copy()
)


# ============================================================
# SOFTWARE-ACCEPTED OBSERVATIONS
# ============================================================

ap_obs = (
    ap_evidence[
        [
            "Program",
            "Sample",
            "Peptide"
        ]
    ]
    .drop_duplicates()
)


fp_obs = (
    fp_evidence
    .drop_duplicates()
)


mm_obs = (
    mm_final
    .drop_duplicates()
)


mq_obs = (
    mq_evidence
    .drop_duplicates()
)


accepted = pd.concat(
    [
        ap_obs,
        fp_obs,
        mm_obs,
        mq_obs
    ],
    ignore_index=True
)


accepted = (
    accepted
    .drop_duplicates(
        subset=[
            "Program",
            "Sample",
            "Peptide"
        ]
    )
    .reset_index(
        drop=True
    )
)


# ============================================================
# LOAD EXISTING I/L-EQUIVALENT GENE-AWARE MAPPING
# ============================================================

print()
print("=" * 110)
print("3. COMMON-REFERENCE MAPPING FILTER")
print("=" * 110)


gene_map = pd.read_csv(
    GENE_MAP_FILE,
    dtype=str,
    low_memory=False
)


gene_map[
    "Peptide"
] = (
    gene_map[
        "Peptide"
    ]
    .fillna("")
    .astype(str)
    .str.strip()
    .str.upper()
)


accepted_detail = accepted.merge(
    gene_map,
    on="Peptide",
    how="left"
)


accepted_detail[
    "ExclusionReason"
] = ""


missing_mapping = (
    accepted_detail[
        "GeneAwareClass"
    ]
    .isna()
)


invalid_mapping = (
    accepted_detail[
        "GeneAwareClass"
    ]
    .eq(
        "invalid_sequence"
    )
)


unmapped_mapping = (
    accepted_detail[
        "GeneAwareClass"
    ]
    .eq(
        "unmapped"
    )
)


accepted_detail.loc[
    missing_mapping,
    "ExclusionReason"
] = (
    "not_present_in_STEP1B_mapping_table"
)


accepted_detail.loc[
    invalid_mapping,
    "ExclusionReason"
] = (
    "invalid_sequence"
)


accepted_detail.loc[
    unmapped_mapping,
    "ExclusionReason"
] = (
    "not_mapped_to_common_target_FASTA"
)


excluded = accepted_detail[
    accepted_detail[
        "ExclusionReason"
    ]
    !=
    ""
].copy()


final_detail = accepted_detail[
    accepted_detail[
        "ExclusionReason"
    ]
    ==
    ""
].copy()


print()
print("SOFTWARE-ACCEPTED vs COMMON-REFERENCE-MAPPED")

comparison_qc_rows = []


for program in PROGRAMS:

    acc = (
        accepted.loc[
            accepted[
                "Program"
            ]
            ==
            program,
            "Peptide"
        ]
        .nunique()
    )


    mapped = (
        final_detail.loc[
            final_detail[
                "Program"
            ]
            ==
            program,
            "Peptide"
        ]
        .nunique()
    )


    excl = (
        excluded.loc[
            excluded[
                "Program"
            ]
            ==
            program,
            "Peptide"
        ]
        .nunique()
    )


    comparison_qc_rows.append({
        "Program":
            program,

        "SoftwareAcceptedDistinctPeptides":
            acc,

        "CommonReferenceMappedDistinctPeptides":
            mapped,

        "ExcludedDistinctPeptides":
            excl,

        "ExcludedPercent":
            (
                100 * excl / acc
                if acc
                else 0
            )
    })


comparison_qc = pd.DataFrame(
    comparison_qc_rows
)


print(
    comparison_qc.to_string(
        index=False
    )
)


# ============================================================
# FINAL SUMMARY BY PROGRAM
# ============================================================

def summarize_subset(
    df,
    program,
    sample="ALL_9_RUNS"
):

    x = df[
        df["Program"]
        ==
        program
    ]


    if sample != "ALL_9_RUNS":

        x = x[
            x["Sample"]
            ==
            sample
        ]


    x = (
        x
        .drop_duplicates(
            subset=[
                "Peptide"
            ]
        )
    )


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


    single_isoforms = set()

    for value in single[
        "MappedAccessions"
    ].fillna(""):

        for acc in str(
            value
        ).split(";"):

            acc = acc.strip()

            if (
                acc
                and
                is_isoform(
                    acc
                )
            ):
                single_isoforms.add(
                    acc
                )


    subset_isoforms = set()

    for value in subset[
        "MappedAccessions"
    ].fillna(""):

        for acc in str(
            value
        ).split(";"):

            acc = acc.strip()

            if (
                acc
                and
                is_isoform(
                    acc
                )
            ):
                subset_isoforms.add(
                    acc
                )


    supported = (
        single_isoforms
        |
        subset_isoforms
    )


    total = x[
        "Peptide"
    ].nunique()


    iso_disc = (
        counts.get(
            "single_isoform_unique",
            0
        )
        +
        counts.get(
            "within_family_subset_discriminative",
            0
        )
    )


    return {
        "Program":
            program,

        "Sample":
            sample,

        "DistinctPeptides":
            total,

        "SingleIsoformUnique":
            counts.get(
                "single_isoform_unique",
                0
            ),

        "SubsetDiscriminative":
            counts.get(
                "within_family_subset_discriminative",
                0
            ),

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
                supported
            )
    }


summary_all = pd.DataFrame(
    [
        summarize_subset(
            final_detail,
            p
        )
        for p in PROGRAMS
    ]
)


summary_sample_rows = []

for program in PROGRAMS:

    for sample in SAMPLES:

        summary_sample_rows.append(
            summarize_subset(
                final_detail,
                program,
                sample
            )
        )


summary_sample = pd.DataFrame(
    summary_sample_rows
)


# ============================================================
# SUPPORTED ISOFORM ACCESSION TABLE
# ============================================================

supported_rows = []


for program in PROGRAMS:

    x = (
        final_detail[
            final_detail[
                "Program"
            ]
            ==
            program
        ]
        .drop_duplicates(
            subset=[
                "Peptide"
            ]
        )
    )


    for _, r in x.iterrows():

        cls = r[
            "GeneAwareClass"
        ]


        if cls not in {
            "single_isoform_unique",
            "within_family_subset_discriminative"
        }:

            continue


        for acc in str(
            r[
                "MappedAccessions"
            ]
        ).split(";"):

            acc = acc.strip()

            if not (
                acc
                and
                is_isoform(
                    acc
                )
            ):

                continue


            supported_rows.append({
                "Program":
                    program,

                "IsoformAccession":
                    acc,

                "Peptide":
                    r[
                        "Peptide"
                    ],

                "EvidenceClass":
                    cls
            })


supported_detail = (
    pd.DataFrame(
        supported_rows
    )
    .drop_duplicates()
)


# ============================================================
# PEPTIDE MEMBERSHIP — ALL
# ============================================================

final_unique = (
    final_detail[
        [
            "Program",
            "Peptide",
            "GeneAwareClass"
        ]
    ]
    .drop_duplicates()
)


all_peptide_union = sorted(
    set(
        final_unique[
            "Peptide"
        ]
    )
)


membership_rows = []


for peptide in all_peptide_union:

    programs = set(
        final_unique.loc[
            final_unique[
                "Peptide"
            ]
            ==
            peptide,
            "Program"
        ]
    )


    membership_rows.append({
        "Peptide":
            peptide,

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
            )
    })


membership_all = pd.DataFrame(
    membership_rows
)


# ============================================================
# PEPTIDE MEMBERSHIP — ISOFORM-DISCRIMINATIVE ONLY
# ============================================================

iso_classes = {
    "single_isoform_unique",
    "within_family_subset_discriminative"
}


iso_peps = final_unique[
    final_unique[
        "GeneAwareClass"
    ].isin(
        iso_classes
    )
]


iso_union = sorted(
    set(
        iso_peps[
            "Peptide"
        ]
    )
)


iso_membership_rows = []


for peptide in iso_union:

    programs = set(
        iso_peps.loc[
            iso_peps[
                "Peptide"
            ]
            ==
            peptide,
            "Program"
        ]
    )


    iso_membership_rows.append({
        "Peptide":
            peptide,

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
            )
    })


membership_iso = pd.DataFrame(
    iso_membership_rows
)


# ============================================================
# ISOFORM ACCESSION MEMBERSHIP
# ============================================================

isoform_membership_rows = []


if not supported_detail.empty:

    for acc in sorted(
        set(
            supported_detail[
                "IsoformAccession"
            ]
        )
    ):

        programs = set(
            supported_detail.loc[
                supported_detail[
                    "IsoformAccession"
                ]
                ==
                acc,
                "Program"
            ]
        )


        isoform_membership_rows.append({
            "IsoformAccession":
                acc,

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
            )
        })


isoform_membership = pd.DataFrame(
    isoform_membership_rows
)


# ============================================================
# PLOT-READY COMPOSITION
# ============================================================

plot_categories = [
    "single_isoform_unique",
    "within_family_subset_discriminative",
    "within_family_shared_all",
    "same_gene_multi_entry_shared",
    "cross_gene_shared",
    "multi_entry_partially_gene_resolved",
    "multi_entry_gene_unresolved",
    "single_canonical_unique"
]


plot_rows = []


for program in PROGRAMS:

    x = (
        final_detail[
            final_detail[
                "Program"
            ]
            ==
            program
        ]
        .drop_duplicates(
            subset=[
                "Peptide"
            ]
        )
    )


    total = len(
        x
    )


    for cls in plot_categories:

        n = int(
            (
                x[
                    "GeneAwareClass"
                ]
                ==
                cls
            ).sum()
        )


        plot_rows.append({
            "Program":
                program,

            "Category":
                cls,

            "PeptideCount":
                n,

            "Percent":
                (
                    100
                    *
                    n
                    /
                    total
                    if total
                    else 0
                )
        })


plot_composition = pd.DataFrame(
    plot_rows
)


# ============================================================
# EXCLUSION SUMMARY
# ============================================================

exclusion_summary = (
    excluded[
        [
            "Program",
            "Peptide",
            "GeneAwareClass",
            "ExclusionReason"
        ]
    ]
    .drop_duplicates()
)


# ============================================================
# FINAL STATUS
# ============================================================

missing_after_merge = int(
    final_detail[
        "GeneAwareClass"
    ].isna().sum()
)


invalid_in_final = int(
    (
        final_detail[
            "GeneAwareClass"
        ]
        ==
        "invalid_sequence"
    ).sum()
)


unmapped_in_final = int(
    (
        final_detail[
            "GeneAwareClass"
        ]
        ==
        "unmapped"
    ).sum()
)


if (
    missing_after_merge == 0
    and
    invalid_in_final == 0
    and
    unmapped_in_final == 0
):

    status = "PASS"

else:

    status = "REVIEW"


# ============================================================
# EXPORT
# ============================================================

outputs = {
    "01_FINAL_SoftwareAccepted_PeptideObservations.csv":
        accepted,

    "02_FINAL_CommonReference_PeptideObservations.csv":
        final_detail[
            [
                "Program",
                "Sample",
                "Peptide"
            ]
        ].drop_duplicates(),

    "03_FINAL_PeptideMappingDetail.csv":
        final_detail,

    "04_FINAL_Summary_ByProgram.csv":
        summary_all,

    "05_FINAL_Summary_ByProgramSample.csv":
        summary_sample,

    "06_FINAL_SupportedIsoformEvidence.csv":
        supported_detail,

    "07_FINAL_PeptideMembership_All.csv":
        membership_all,

    "08_FINAL_PeptideMembership_IsoformDiscriminative.csv":
        membership_iso,

    "09_FINAL_IsoformMembership.csv":
        isoform_membership,

    "10_PlotData_GeneAwareComposition.csv":
        plot_composition,

    "11_FINAL_ComparisonQC.csv":
        comparison_qc,

    "12_FINAL_ExcludedPeptides.csv":
        exclusion_summary,

    "13_AP_FinalFilter_QC.csv":
        ap_qc,

    "14_MM_Ambiguity_QC.csv":
        mm_qc
}


for name, df in outputs.items():

    df.to_csv(
        OUT / name,
        index=False,
        encoding="utf-8-sig"
    )


status_text = f"""
STEP 1D FINAL NORMALIZATION

FINAL STATUS: {status}

Rules:
- AP: target == True AND q_value <= 0.01
- FP: existing workflow-filtered MBR-OFF PSM set
- MM: existing target QValue<=0.01 set; pipe-delimited ambiguous Base Sequences excluded from primary unambiguous peptide comparison
- MQ: existing target/non-contaminant MBR-OFF msms set
- Common comparison: peptide must map to common target FASTA using STEP1B I/L-equivalent gene-aware mapping

Missing mappings retained in final: {missing_after_merge}
Invalid sequences retained in final: {invalid_in_final}
Unmapped sequences retained in final: {unmapped_in_final}
"""


with open(
    OUT /
    "15_FINAL_FREEZE_STATUS.txt",
    "w",
    encoding="utf-8"
) as fh:

    fh.write(
        status_text
    )


# ============================================================
# PRINT
# ============================================================

print()
print("=" * 125)
print("FINAL SOFTWARE-ACCEPTED / COMMON-REFERENCE QC")
print("=" * 125)

print(
    comparison_qc.to_string(
        index=False
    )
)


print()
print("=" * 125)
print("FINAL PRIMARY ISOFORM BENCHMARK")
print("=" * 125)


show = [
    "Program",
    "DistinctPeptides",
    "SingleIsoformUnique",
    "SubsetDiscriminative",
    "TotalIsoformDiscriminative",
    "IsoformDiscriminativePercent",
    "WithinFamilySharedAll",
    "SameGeneMultiEntryShared",
    "CrossGeneShared",
    "DistinctSingleIsoformAccessions",
    "DistinctSubsetSupportedIsoforms",
    "DistinctDiscriminativelySupportedIsoforms"
]


print(
    summary_all[
        show
    ].to_string(
        index=False
    )
)


print()
print("=" * 125)
print("AP FINAL FILTER QC")
print("=" * 125)

print(
    ap_qc.to_string(
        index=False
    )
)


print()
print("=" * 125)
print("MM AMBIGUITY QC")
print("=" * 125)

print(
    mm_qc.to_string(
        index=False
    )
)


print()
print("=" * 125)
print("FINAL STATUS:", status)
print("=" * 125)

print()
print("OUTPUT:")
print(OUT)

print()
print("STEP 1D COMPLETE")

