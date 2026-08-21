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

from collections import defaultdict, deque
import pandas as pd
import numpy as np
import h5py
import re
import time


# =====================================================================
# PATHS
# =====================================================================

ROOT = _public_project_root()

MASTER = ROOT / "MASTER_TABLES_FINAL"

STEP1 = (
    MASTER /
    "STEP1_COMMON_FASTA_MAPPING"
)

FASTA = Path(
    r"<BENCHMARK_ROOT>"
    r"\Output_AP_cano_only\DB"
    r"\uniprotkb_proteome_UP000005640_2026_08_04.fasta"
)

AP_DB = (
    ROOT /
    "AP_MBR_OFF" /
    "database.hdf"
)

PRIMARY_UNIQUE = (
    MASTER /
    "02_PRIMARY_ID_OFF_UniquePeptides.csv"
)

PRIMARY_EVIDENCE = (
    MASTER /
    "01_PRIMARY_ID_OFF_Evidence.csv"
)

STEP1_MAPPING = (
    STEP1 /
    "02_CommonFASTA_PeptideMapping.csv"
)

OUT = (
    MASTER /
    "STEP1B_MAPPING_QC"
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


# =====================================================================
# HELPERS
# =====================================================================

ISO_RE = re.compile(
    r"^(.+?)-([1-9][0-9]*)$"
)

AZ_RE = re.compile(
    r"^[A-Z]+$"
)

STANDARD20 = set(
    "ACDEFGHIKLMNPQRSTVWY"
)


def dec(x):

    if isinstance(x, bytes):
        return x.decode(
            "utf-8",
            errors="replace"
        )

    return str(x)


def parse_accession(header):

    h = str(header).strip()

    if h.startswith(">"):
        h = h[1:]

    first = h.split()[0]

    parts = first.split("|")

    if (
        len(parts) >= 2
        and parts[0].lower()
        in {"sp", "tr"}
    ):
        return parts[1].strip()

    return first.strip()


def clean_accession(x):

    x = dec(x).strip()

    if "|" in x:

        parts = x.split("|")

        if (
            len(parts) >= 2
            and parts[0].lower()
            in {"sp", "tr"}
        ):
            return parts[1].strip()

    return x


def base_accession(acc):

    m = ISO_RE.fullmatch(
        acc
    )

    if m:
        return m.group(1)

    return acc


def is_isoform(acc):

    return (
        ISO_RE.fullmatch(acc)
        is not None
    )


def il_normalize(seq):

    return (
        seq
        .replace("I", "J")
        .replace("L", "J")
    )


def nonstandard_letters(seq):

    return "".join(
        sorted(
            set(seq)
            -
            STANDARD20
        )
    )


def nonletter_chars(seq):

    return "".join(
        sorted(
            {
                x
                for x in seq
                if not (
                    "A" <= x <= "Z"
                )
            }
        )
    )


# =====================================================================
# READ FASTA WITH GN= INFORMATION
# =====================================================================

def read_fasta(path):

    rows = []

    header = None
    parts = []

    def flush():

        nonlocal header
        nonlocal parts

        if header is None:
            return

        seq = (
            "".join(parts)
            .replace(" ", "")
            .upper()
        )

        acc = parse_accession(
            header
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

        if acc and seq:

            rows.append({
                "Accession":
                    acc,

                "Sequence":
                    seq,

                "Gene":
                    gene,

                "Header":
                    header
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
                parts = []

            else:

                parts.append(
                    line
                )

        flush()

    return pd.DataFrame(
        rows
    )


# =====================================================================
# AHO-CORASICK
# =====================================================================

def build_automaton(patterns):

    goto = [{}]
    fail = [0]
    output = [[]]

    for idx, pattern in enumerate(
        patterns
    ):

        state = 0

        for aa in pattern:

            if aa not in goto[state]:

                goto[state][aa] = len(
                    goto
                )

                goto.append({})
                fail.append(0)
                output.append([])

            state = goto[state][aa]

        output[state].append(
            idx
        )

    q = deque()

    for state in goto[0].values():

        fail[state] = 0
        q.append(
            state
        )

    while q:

        r = q.popleft()

        for aa, s in goto[r].items():

            q.append(s)

            f = fail[r]

            while (
                f != 0
                and aa not in goto[f]
            ):
                f = fail[f]

            fail[s] = (
                goto[f].get(
                    aa,
                    0
                )
            )

            output[s].extend(
                output[
                    fail[s]
                ]
            )

    return (
        goto,
        fail,
        output
    )


def map_patterns(
    patterns,
    sequence_to_accessions,
    label
):

    print()
    print("=" * 105)
    print(label)
    print("=" * 105)

    t0 = time.time()

    goto, fail, output = (
        build_automaton(
            patterns
        )
    )

    mapped = [
        set()
        for _ in patterns
    ]

    total = len(
        sequence_to_accessions
    )

    print(
        "Patterns          :",
        len(patterns)
    )

    print(
        "Protein sequences :",
        total
    )

    print(
        "Automaton states  :",
        len(goto)
    )

    for n, (
        sequence,
        accessions
    ) in enumerate(
        sequence_to_accessions.items(),
        1
    ):

        state = 0
        found = set()

        for aa in sequence:

            while (
                state != 0
                and aa not in goto[state]
            ):
                state = fail[state]

            state = (
                goto[state].get(
                    aa,
                    0
                )
            )

            if output[state]:

                found.update(
                    output[state]
                )

        for idx in found:

            mapped[idx].update(
                accessions
            )

        if n % 20000 == 0:

            print(
                f"  scanned {n:,}/{total:,}"
            )

    print(
        "Finished in",
        round(
            time.time() - t0,
            1
        ),
        "sec"
    )

    return mapped


# =====================================================================
# START
# =====================================================================

start_all = time.time()

print()
print("=" * 105)
print("STEP 1B — MAPPING QC + GENE-AWARE CLASSIFICATION")
print("=" * 105)


# =====================================================================
# FASTA
# =====================================================================

fasta_df = read_fasta(
    FASTA
)

accession_to_seq = dict(
    zip(
        fasta_df["Accession"],
        fasta_df["Sequence"]
    )
)

accession_to_gene = dict(
    zip(
        fasta_df["Accession"],
        fasta_df["Gene"]
    )
)


family_members = defaultdict(
    set
)

for acc in accession_to_seq:

    family_members[
        base_accession(acc)
    ].add(
        acc
    )


print()
print("FASTA entries :", len(fasta_df))
print(
    "GN= available:",
    (
        fasta_df["Gene"]
        .astype(str)
        .str.len()
        .gt(0)
        .sum()
    )
)


# =====================================================================
# 1. AP DATABASE vs FASTA CONTENT QC
# =====================================================================

print()
print("=" * 105)
print("1. FASTA vs ALPHAPEPT DATABASE CONTENT QC")
print("=" * 105)

with h5py.File(
    AP_DB,
    "r"
) as h5:

    ap_ids = [
        clean_accession(x)
        for x in h5["proteins/id"][:]
    ]

    ap_seqs = [
        dec(x).strip().upper()
        for x in h5["proteins/sequence"][:]
    ]


ap_dict = dict(
    zip(
        ap_ids,
        ap_seqs
    )
)


fasta_ids = set(
    accession_to_seq
)

ap_id_set = set(
    ap_dict
)


common_ids = (
    fasta_ids
    &
    ap_id_set
)

only_fasta = (
    fasta_ids
    -
    ap_id_set
)

only_ap = (
    ap_id_set
    -
    fasta_ids
)


sequence_mismatches = []

for acc in common_ids:

    if (
        accession_to_seq[acc]
        !=
        ap_dict[acc]
    ):

        sequence_mismatches.append({
            "Accession":
                acc,

            "FASTA_Length":
                len(
                    accession_to_seq[acc]
                ),

            "APDB_Length":
                len(
                    ap_dict[acc]
                ),

            "FASTA_Sequence":
                accession_to_seq[acc],

            "APDB_Sequence":
                ap_dict[acc]
        })


apdb_qc = pd.DataFrame(
    [
        {
            "Metric":
                "FASTA_Accessions",

            "Value":
                len(fasta_ids)
        },
        {
            "Metric":
                "APDB_Accessions",

            "Value":
                len(ap_id_set)
        },
        {
            "Metric":
                "Common_Accessions",

            "Value":
                len(common_ids)
        },
        {
            "Metric":
                "Only_FASTA",

            "Value":
                len(only_fasta)
        },
        {
            "Metric":
                "Only_APDB",

            "Value":
                len(only_ap)
        },
        {
            "Metric":
                "Sequence_Mismatches",

            "Value":
                len(sequence_mismatches)
        }
    ]
)


print(
    apdb_qc.to_string(
        index=False
    )
)


# =====================================================================
# PRIMARY PEPTIDES
# =====================================================================

primary = pd.read_csv(
    PRIMARY_UNIQUE,
    dtype=str
)


primary = (
    primary[
        [
            "Program",
            "Sample",
            "Peptide"
        ]
    ]
    .copy()
)


primary["Peptide"] = (
    primary["Peptide"]
    .fillna("")
    .astype(str)
    .str.strip()
    .str.upper()
)


primary = primary[
    primary["Peptide"] != ""
]


primary = primary.drop_duplicates(
    subset=[
        "Program",
        "Sample",
        "Peptide"
    ]
)


all_peptides = sorted(
    set(
        primary["Peptide"]
    )
)


# =====================================================================
# READ STEP 1 RESULT
# =====================================================================

step1_map = pd.read_csv(
    STEP1_MAPPING,
    dtype=str
)


step1_map["Peptide"] = (
    step1_map["Peptide"]
    .fillna("")
    .str.upper()
)


# =====================================================================
# 2. INVALID / NON-STANDARD PEPTIDE AUDIT
# =====================================================================

print()
print("=" * 105)
print("2. INVALID / NON-STANDARD PEPTIDE AUDIT")
print("=" * 105)


syntax_rows = []


for peptide in all_peptides:

    programs = sorted(
        set(
            primary.loc[
                primary["Peptide"]
                == peptide,
                "Program"
            ]
        )
    )

    syntax_rows.append({
        "Peptide":
            peptide,

        "Length":
            len(peptide),

        "AZ_Only":
            bool(
                AZ_RE.fullmatch(
                    peptide
                )
            ),

        "Standard20Only":
            set(peptide).issubset(
                STANDARD20
            ),

        "NonLetterCharacters":
            nonletter_chars(
                peptide
            ),

        "NonStandardLetters":
            nonstandard_letters(
                peptide
            ),

        "Programs":
            ";".join(
                programs
            )
    })


syntax_df = pd.DataFrame(
    syntax_rows
)


invalid_df = syntax_df[
    ~syntax_df["AZ_Only"]
].copy()


nonstandard_df = syntax_df[
    (
        syntax_df["AZ_Only"]
    )
    &
    (
        ~syntax_df["Standard20Only"]
    )
].copy()


print(
    "Non-AZ peptide strings :",
    len(invalid_df)
)

print(
    "AZ but non-standard AA :",
    len(nonstandard_df)
)


# =====================================================================
# EVIDENCE CONTEXT
# =====================================================================

evidence = pd.read_csv(
    PRIMARY_EVIDENCE,
    dtype=str,
    low_memory=False
)


evidence["Peptide"] = (
    evidence["Peptide"]
    .fillna("")
    .astype(str)
    .str.strip()
    .str.upper()
)


def summarize_evidence_context(
    peptide
):

    x = evidence[
        evidence["Peptide"]
        == peptide
    ]

    if x.empty:

        return {
            "EvidenceRows": 0,
            "Samples": "",
            "ModifiedForms": "",
            "SourceFiles": ""
        }

    return {
        "EvidenceRows":
            len(x),

        "Samples":
            ";".join(
                sorted(
                    set(
                        x["Sample"]
                        .dropna()
                        .astype(str)
                    )
                )
            ),

        "ModifiedForms":
            ";".join(
                sorted(
                    set(
                        x[
                            "ModifiedPeptide"
                        ]
                        .dropna()
                        .astype(str)
                    )
                )
            )[:20],

        "SourceFiles":
            ";".join(
                sorted(
                    set(
                        x[
                            "SourceFile"
                        ]
                        .dropna()
                        .astype(str)
                    )
                )
            )[:20]
    }


# =====================================================================
# 3. RE-MAP ALL PEPTIDES TO FASTA USING I/L EQUIVALENCE
# =====================================================================

valid_peptides = [
    x
    for x in all_peptides
    if AZ_RE.fullmatch(x)
]


normalized_to_original = defaultdict(
    set
)


for peptide in valid_peptides:

    normalized_to_original[
        il_normalize(
            peptide
        )
    ].add(
        peptide
    )


patterns = sorted(
    normalized_to_original
)


fasta_il_seq_to_acc = defaultdict(
    set
)


for acc, seq in accession_to_seq.items():

    fasta_il_seq_to_acc[
        il_normalize(
            seq
        )
    ].add(
        acc
    )


mapped_raw = map_patterns(
    patterns,
    fasta_il_seq_to_acc,
    "3. I/L-EQUIVALENT COMMON FASTA REMAPPING"
)


peptide_to_accessions = {}


for i, pattern in enumerate(
    patterns
):

    accs = mapped_raw[i]

    for original in (
        normalized_to_original[
            pattern
        ]
    ):

        peptide_to_accessions[
            original
        ] = set(
            accs
        )


# =====================================================================
# GENE-AWARE CLASSIFICATION
# =====================================================================

def gene_aware_classification(
    peptide,
    accessions
):

    if not AZ_RE.fullmatch(
        peptide
    ):

        return {
            "GeneAwareClass":
                "invalid_sequence",

            "BaseEntryCount":
                0,

            "KnownGeneCount":
                0,

            "KnownGenes":
                "",

            "UnknownGeneAccessions":
                0
        }


    if not accessions:

        return {
            "GeneAwareClass":
                "unmapped",

            "BaseEntryCount":
                0,

            "KnownGeneCount":
                0,

            "KnownGenes":
                "",

            "UnknownGeneAccessions":
                0
        }


    bases = {
        base_accession(x)
        for x in accessions
    }


    known_genes = {
        accession_to_gene.get(
            x,
            ""
        )
        for x in accessions
        if accession_to_gene.get(
            x,
            ""
        )
    }


    unknown_gene_n = sum(
        1
        for x in accessions
        if not accession_to_gene.get(
            x,
            ""
        )
    )


    # ---------------------------------------------------------
    # one exact accession
    # ---------------------------------------------------------

    if len(accessions) == 1:

        only = next(
            iter(accessions)
        )

        if is_isoform(
            only
        ):

            cls = (
                "single_isoform_unique"
            )

        else:

            cls = (
                "single_canonical_unique"
            )


    # ---------------------------------------------------------
    # one UniProt accession family
    # ---------------------------------------------------------

    elif len(bases) == 1:

        base = next(
            iter(bases)
        )

        if (
            accessions
            ==
            family_members[
                base
            ]
        ):

            cls = (
                "within_family_shared_all"
            )

        else:

            cls = (
                "within_family_subset_discriminative"
            )


    # ---------------------------------------------------------
    # multiple base accessions
    # ---------------------------------------------------------

    else:

        if len(
            known_genes
        ) >= 2:

            cls = (
                "cross_gene_shared"
            )

        elif (
            len(known_genes) == 1
            and unknown_gene_n == 0
        ):

            cls = (
                "same_gene_multi_entry_shared"
            )

        elif (
            len(known_genes) == 0
        ):

            cls = (
                "multi_entry_gene_unresolved"
            )

        else:

            cls = (
                "multi_entry_partially_gene_resolved"
            )


    return {
        "GeneAwareClass":
            cls,

        "BaseEntryCount":
            len(bases),

        "KnownGeneCount":
            len(
                known_genes
            ),

        "KnownGenes":
            ";".join(
                sorted(
                    known_genes
                )
            ),

        "UnknownGeneAccessions":
            unknown_gene_n
    }


gene_rows = []


for peptide in all_peptides:

    accs = peptide_to_accessions.get(
        peptide,
        set()
    )


    result = gene_aware_classification(
        peptide,
        accs
    )


    gene_rows.append({
        "Peptide":
            peptide,

        "PeptideLength":
            len(peptide),

        "MappedAccessionCount":
            len(accs),

        "MappedAccessions":
            ";".join(
                sorted(
                    accs
                )
            ),

        **result
    })


gene_map = pd.DataFrame(
    gene_rows
)


# =====================================================================
# 4. UNMAPPED AUDIT
# =====================================================================

print()
print("=" * 105)
print("4. UNMAPPED PEPTIDE AUDIT")
print("=" * 105)


detail = primary.merge(
    gene_map,
    on="Peptide",
    how="left"
)


unmapped_obs = detail[
    detail["GeneAwareClass"]
    ==
    "unmapped"
].copy()


unmapped_peptides = (
    unmapped_obs[
        [
            "Program",
            "Peptide"
        ]
    ]
    .drop_duplicates()
)


# =====================================================================
# CHECK AP UNMAPPED AGAINST AP DATABASE ITSELF
# =====================================================================

ap_unmapped = sorted(
    set(
        unmapped_peptides.loc[
            unmapped_peptides[
                "Program"
            ] == "AP",
            "Peptide"
        ]
    )
)


ap_exact_found = {}
ap_il_found = {}


if ap_unmapped:

    # ----------------------------------------
    # exact
    # ----------------------------------------

    ap_seq_to_acc = defaultdict(
        set
    )

    for acc, seq in ap_dict.items():

        ap_seq_to_acc[
            seq
        ].add(
            acc
        )


    exact_results = map_patterns(
        ap_unmapped,
        ap_seq_to_acc,
        "4A. AP UNMAPPED → AP DATABASE EXACT"
    )


    for i, pep in enumerate(
        ap_unmapped
    ):

        ap_exact_found[
            pep
        ] = len(
            exact_results[i]
        ) > 0


    # ----------------------------------------
    # I/L equivalent
    # ----------------------------------------

    norm_to_original = defaultdict(
        set
    )

    for pep in ap_unmapped:

        norm_to_original[
            il_normalize(
                pep
            )
        ].add(
            pep
        )


    ap_patterns = sorted(
        norm_to_original
    )


    ap_il_seq_to_acc = defaultdict(
        set
    )


    for acc, seq in ap_dict.items():

        ap_il_seq_to_acc[
            il_normalize(
                seq
            )
        ].add(
            acc
        )


    il_results = map_patterns(
        ap_patterns,
        ap_il_seq_to_acc,
        "4B. AP UNMAPPED → AP DATABASE I/L"
    )


    for i, pattern in enumerate(
        ap_patterns
    ):

        found = (
            len(
                il_results[i]
            )
            > 0
        )

        for original in (
            norm_to_original[
                pattern
            ]
        ):

            ap_il_found[
                original
            ] = found


# =====================================================================
# BUILD UNMAPPED AUDIT TABLE
# =====================================================================

unmapped_audit_rows = []


for _, r in (
    unmapped_peptides
    .iterrows()
):

    program = r[
        "Program"
    ]

    peptide = r[
        "Peptide"
    ]


    context = (
        summarize_evidence_context(
            peptide
        )
    )


    step1_row = step1_map[
        step1_map["Peptide"]
        ==
        peptide
    ]


    exact_category = (
        step1_row.iloc[0][
            "ExactCategory"
        ]
        if not step1_row.empty
        else ""
    )


    il_category = (
        step1_row.iloc[0][
            "ILCategory"
        ]
        if not step1_row.empty
        else ""
    )


    unmapped_audit_rows.append({
        "Program":
            program,

        "Peptide":
            peptide,

        "Length":
            len(peptide),

        "Standard20Only":
            set(peptide).issubset(
                STANDARD20
            ),

        "NonStandardLetters":
            nonstandard_letters(
                peptide
            ),

        "STEP1_ExactCategory":
            exact_category,

        "STEP1_ILCategory":
            il_category,

        "FoundInAPDatabase_Exact":
            (
                ap_exact_found.get(
                    peptide,
                    ""
                )
                if program == "AP"
                else ""
            ),

        "FoundInAPDatabase_IL":
            (
                ap_il_found.get(
                    peptide,
                    ""
                )
                if program == "AP"
                else ""
            ),

        **context
    })


unmapped_audit = pd.DataFrame(
    unmapped_audit_rows
)


# =====================================================================
# GENE-AWARE SUMMARY BY PROGRAM
# =====================================================================

summary_rows = []


for program in PROGRAMS:

    x = (
        detail[
            detail["Program"]
            ==
            program
        ]
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


    summary_rows.append({
        "Program":
            program,

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

        "Unmapped":
            counts.get(
                "unmapped",
                0
            ),

        "Invalid":
            counts.get(
                "invalid_sequence",
                0
            ),

        "UnmappedPercent":
            (
                100
                *
                counts.get(
                    "unmapped",
                    0
                )
                /
                total
                if total
                else 0
            ),

        "IsoformDiscriminativePercent":
            (
                100
                *
                iso_disc
                /
                total
                if total
                else 0
            )
    })


summary_df = pd.DataFrame(
    summary_rows
)


# =====================================================================
# CROSS-ENTRY BREAKDOWN
# =====================================================================

cross_classes = {
    "same_gene_multi_entry_shared",
    "cross_gene_shared",
    "multi_entry_partially_gene_resolved",
    "multi_entry_gene_unresolved"
}


cross_detail = detail[
    detail[
        "GeneAwareClass"
    ].isin(
        cross_classes
    )
].copy()


cross_summary = (
    cross_detail
    .drop_duplicates(
        subset=[
            "Program",
            "Peptide"
        ]
    )
    .groupby(
        [
            "Program",
            "GeneAwareClass"
        ]
    )
    .size()
    .reset_index(
        name="PeptideCount"
    )
)


# =====================================================================
# PLOT-READY: COMPOSITION
# =====================================================================

plot_classes = [
    "single_isoform_unique",
    "within_family_subset_discriminative",
    "within_family_shared_all",
    "same_gene_multi_entry_shared",
    "cross_gene_shared",
    "multi_entry_partially_gene_resolved",
    "multi_entry_gene_unresolved",
    "single_canonical_unique",
    "unmapped",
    "invalid_sequence"
]


plot_rows = []


for program in PROGRAMS:

    x = (
        detail[
            detail["Program"]
            ==
            program
        ]
        .drop_duplicates(
            subset=[
                "Peptide"
            ]
        )
    )


    total = len(x)


    for category in plot_classes:

        n = int(
            (
                x[
                    "GeneAwareClass"
                ]
                ==
                category
            ).sum()
        )


        plot_rows.append({
            "Program":
                program,

            "Category":
                category,

            "PeptideCount":
                n,

            "Percent":
                (
                    100 * n / total
                    if total
                    else 0
                )
        })


plot_composition = pd.DataFrame(
    plot_rows
)


# =====================================================================
# PLOT-READY: UNMAPPED
# =====================================================================

plot_unmapped = (
    summary_df[
        [
            "Program",
            "DistinctPeptides",
            "Unmapped",
            "UnmappedPercent"
        ]
    ]
    .copy()
)


# =====================================================================
# INVALID BY PROGRAM
# =====================================================================

invalid_program_rows = []


invalid_peptide_set = set(
    invalid_df[
        "Peptide"
    ]
)


for program in PROGRAMS:

    x = (
        primary[
            (
                primary[
                    "Program"
                ]
                ==
                program
            )
            &
            (
                primary[
                    "Peptide"
                ].isin(
                    invalid_peptide_set
                )
            )
        ]
    )


    invalid_program_rows.append({
        "Program":
            program,

        "DistinctInvalidPeptides":
            x[
                "Peptide"
            ].nunique()
    })


invalid_program_df = pd.DataFrame(
    invalid_program_rows
)


# =====================================================================
# LENGTH DISTRIBUTION — PLOT READY
# =====================================================================

length_plot = (
    detail[
        [
            "Program",
            "Peptide",
            "PeptideLength",
            "GeneAwareClass"
        ]
    ]
    .drop_duplicates(
        subset=[
            "Program",
            "Peptide"
        ]
    )
    .copy()
)


length_plot[
    "MappedStatus"
] = np.where(
    length_plot[
        "GeneAwareClass"
    ].eq(
        "unmapped"
    ),
    "Unmapped",
    "Mapped"
)


# =====================================================================
# MASTER QC SUMMARY
# =====================================================================

ap_unmapped_found_exact = sum(
    1
    for x in ap_unmapped
    if ap_exact_found.get(
        x,
        False
    )
)


ap_unmapped_found_il = sum(
    1
    for x in ap_unmapped
    if ap_il_found.get(
        x,
        False
    )
)


qc_rows = [
    {
        "Metric":
            "FASTA_Accessions",

        "Value":
            len(fasta_ids)
    },
    {
        "Metric":
            "APDB_Accessions",

        "Value":
            len(ap_id_set)
    },
    {
        "Metric":
            "FASTA_APDB_OnlyFASTA",

        "Value":
            len(only_fasta)
    },
    {
        "Metric":
            "FASTA_APDB_OnlyAPDB",

        "Value":
            len(only_ap)
    },
    {
        "Metric":
            "FASTA_APDB_SequenceMismatches",

        "Value":
            len(sequence_mismatches)
    },
    {
        "Metric":
            "UnionPrimaryPeptides",

        "Value":
            len(all_peptides)
    },
    {
        "Metric":
            "NonAZ_InvalidPeptides",

        "Value":
            len(invalid_df)
    },
    {
        "Metric":
            "AZ_NonStandardAA_Peptides",

        "Value":
            len(nonstandard_df)
    },
    {
        "Metric":
            "AP_Unmapped_IL",

        "Value":
            len(ap_unmapped)
    },
    {
        "Metric":
            "AP_Unmapped_FoundBackInAPDB_Exact",

        "Value":
            ap_unmapped_found_exact
    },
    {
        "Metric":
            "AP_Unmapped_FoundBackInAPDB_IL",

        "Value":
            ap_unmapped_found_il
    }
]


qc_df = pd.DataFrame(
    qc_rows
)


# =====================================================================
# WRITE
# =====================================================================

outputs = {
    "01_FASTA_vs_APDB_Content_QC.csv":
        apdb_qc,

    "02_FASTA_vs_APDB_SequenceMismatches.csv":
        pd.DataFrame(
            sequence_mismatches
        ),

    "03_Unmapped_Peptide_Audit.csv":
        unmapped_audit,

    "04_Invalid_Peptide_Audit.csv":
        invalid_df,

    "05_NonStandardAA_Peptide_Audit.csv":
        nonstandard_df,

    "06_GeneAware_PeptideMapping.csv":
        gene_map,

    "07_GeneAware_PeptideDetail_ByProgram.csv":
        detail,

    "08_GeneAware_Summary_ByProgram.csv":
        summary_df,

    "09_CrossEntry_GeneBreakdown.csv":
        cross_summary,

    "10_PlotData_GeneAware_Composition.csv":
        plot_composition,

    "11_PlotData_Unmapped_Rates.csv":
        plot_unmapped,

    "12_PlotData_PeptideLength_Mapped_vs_Unmapped.csv":
        length_plot,

    "13_Invalid_Counts_ByProgram.csv":
        invalid_program_df,

    "14_STEP1B_QC_Summary.csv":
        qc_df
}


for name, df in outputs.items():

    df.to_csv(
        OUT / name,
        index=False,
        encoding="utf-8-sig"
    )


# =====================================================================
# PRINT FINAL RESULTS
# =====================================================================

print()
print("=" * 120)
print("FASTA vs AP DATABASE")
print("=" * 120)

print(
    apdb_qc.to_string(
        index=False
    )
)


print()
print("=" * 120)
print("INVALID / NON-STANDARD")
print("=" * 120)

print(
    invalid_program_df.to_string(
        index=False
    )
)

print()

print(
    "Union non-AZ invalid:",
    len(invalid_df)
)

print(
    "Union AZ/non-standard amino acid:",
    len(nonstandard_df)
)


print()
print("=" * 120)
print("GENE-AWARE PRIMARY SUMMARY")
print("=" * 120)

show = [
    "Program",
    "DistinctPeptides",
    "SingleIsoformUnique",
    "SubsetDiscriminative",
    "TotalIsoformDiscriminative",
    "WithinFamilySharedAll",
    "SameGeneMultiEntryShared",
    "CrossGeneShared",
    "PartiallyGeneResolved",
    "GeneUnresolved",
    "Unmapped",
    "UnmappedPercent"
]


print(
    summary_df[
        show
    ].to_string(
        index=False
    )
)


print()
print("=" * 120)
print("AP UNMAPPED DIAGNOSTIC")
print("=" * 120)

print(
    "AP unmapped against common FASTA       :",
    len(ap_unmapped)
)

print(
    "Found back in AP database — exact      :",
    ap_unmapped_found_exact
)

print(
    "Found back in AP database — I/L equiv. :",
    ap_unmapped_found_il
)


print()
print("=" * 120)
print("CROSS-ENTRY → GENE-AWARE BREAKDOWN")
print("=" * 120)

print(
    cross_summary.to_string(
        index=False
    )
)


print()
print("OUTPUT DIRECTORY:")
print(OUT)

print()
print("Files created:")

for name in outputs:
    print(
        " ",
        OUT / name
    )


print()
print(
    "Total runtime:",
    round(
        time.time()
        -
        start_all,
        1
    ),
    "sec"
)

print()
print("STEP 1B COMPLETE")

