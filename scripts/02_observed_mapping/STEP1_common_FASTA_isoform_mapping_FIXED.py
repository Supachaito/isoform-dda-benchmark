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
import re
import time


# ============================================================
# PATHS
# ============================================================

ROOT = _public_project_root()

MASTER = ROOT / "MASTER_TABLES_FINAL"

PRIMARY_FILE = (
    MASTER /
    "02_PRIMARY_ID_OFF_UniquePeptides.csv"
)

FASTA_DIR = _public_project_root().parent / "Output_AP_cano_only" / "DB"

FASTA_STEM = "uniprotkb_proteome_UP000005640_2026_08_04"

OUT = (
    MASTER /
    "STEP1_COMMON_FASTA_MAPPING"
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


# ============================================================
# FIND FASTA
# ============================================================

candidates = []

for p in FASTA_DIR.glob(
    FASTA_STEM + "*"
):
    if p.is_file():

        # accept normal FASTA extensions,
        # but also accept file with no extension
        if (
            p.suffix.lower()
            in {
                ".fasta",
                ".fa",
                ".fas",
                ".faa"
            }
            or p.suffix == ""
        ):
            candidates.append(p)


if len(candidates) == 0:

    raise FileNotFoundError(
        "Could not find FASTA beginning with:\n"
        + str(
            FASTA_DIR /
            FASTA_STEM
        )
    )


if len(candidates) > 1:

    print()
    print("Multiple candidate FASTA files found:")

    for x in candidates:
        print(" ", x)

    raise RuntimeError(
        "More than one FASTA candidate found. "
        "Please keep only the intended FASTA or specify the full filename."
    )


FASTA = candidates[0]


# ============================================================
# HELPERS
# ============================================================

ISO_RE = re.compile(
    r"^(.+?)-([1-9][0-9]*)$"
)

VALID_PEPTIDE_RE = re.compile(
    r"^[A-Z]+$"
)


def parse_accession(header):

    h = header.strip()

    if h.startswith(">"):
        h = h[1:]

    first = h.split()[0]

    # UniProt format:
    # sp|P12345|NAME
    # sp|P12345-2|NAME

    parts = first.split("|")

    if (
        len(parts) >= 2
        and
        parts[0].lower()
        in {"sp", "tr"}
    ):
        return parts[1].strip()

    return first.strip()


def is_isoform(acc):

    return (
        ISO_RE.fullmatch(acc)
        is not None
    )


def base_accession(acc):

    m = ISO_RE.fullmatch(acc)

    if m:
        return m.group(1)

    return acc


def is_decoy_or_contaminant(acc):

    x = acc.upper()

    return x.startswith(
        (
            "REV__",
            "REV_",
            "REVERSE_",
            "DECOY_",
            "DECOY-",
            "CON__",
            "CONTAM_",
            "CONTAMINANT_"
        )
    )


def il_normalize(seq):

    # Conventional MS/MS does not reliably distinguish I/L.
    # Collapse both to J for conservative mapping.

    return (
        seq
        .replace("I", "J")
        .replace("L", "J")
    )


def read_fasta(path):

    output = []

    header = None
    seq_parts = []

    def flush():

        nonlocal header
        nonlocal seq_parts

        if header is None:
            return

        acc = parse_accession(
            header
        )

        seq = (
            "".join(seq_parts)
            .replace(" ", "")
            .upper()
        )

        if acc and seq:

            output.append(
                (
                    acc,
                    seq
                )
            )

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

    return output


def compact_list(values, limit=100):

    vals = sorted(values)

    if len(vals) <= limit:
        return ";".join(vals)

    return (
        ";".join(
            vals[:limit]
        )
        +
        f";...TRUNCATED_TOTAL={len(vals)}"
    )


# ============================================================
# START
# ============================================================

start_all = time.time()

print()
print("=" * 110)
print("STEP 1 — COMMON FASTA PEPTIDE → ISOFORM MAPPING")
print("=" * 110)

print()
print("FASTA:")
print(FASTA)

print()
print("Primary peptide table:")
print(PRIMARY_FILE)


# ============================================================
# READ FASTA
# ============================================================

raw_reference = read_fasta(
    FASTA
)

raw_count = len(
    raw_reference
)


reference = []

decoy_count = 0


for acc, seq in raw_reference:

    if is_decoy_or_contaminant(
        acc
    ):

        decoy_count += 1

    else:

        reference.append(
            (
                acc,
                seq
            )
        )


# accession -> sequence(s)

accession_sequences = defaultdict(
    set
)


for acc, seq in reference:

    accession_sequences[
        acc
    ].add(
        seq
    )


all_accessions = set(
    accession_sequences
)


isoform_accessions = {
    x
    for x in all_accessions
    if is_isoform(x)
}


family_members = defaultdict(
    set
)


for acc in all_accessions:

    family_members[
        base_accession(acc)
    ].add(
        acc
    )


families_with_isoforms = {

    base

    for base, members
    in family_members.items()

    if any(
        is_isoform(x)
        for x in members
    )
}


conflicting_accessions = {

    acc: seqs

    for acc, seqs
    in accession_sequences.items()

    if len(seqs) > 1
}


print()
print("=" * 110)
print("FASTA QC")
print("=" * 110)

print(
    "Raw FASTA entries              :",
    raw_count
)

print(
    "Decoy/contaminant excluded     :",
    decoy_count
)

print(
    "Target entries                 :",
    len(reference)
)

print(
    "Distinct target accessions     :",
    len(all_accessions)
)

print(
    "Suffix-bearing isoform entries :",
    len(isoform_accessions)
)

print(
    "Protein families               :",
    len(family_members)
)

print(
    "Families containing isoforms   :",
    len(families_with_isoforms)
)

print(
    "Conflicting duplicate accessions:",
    len(conflicting_accessions)
)


# expected benchmark FASTA sanity check

if (
    raw_count != 169637
    or
    len(isoform_accessions) != 22131
):

    print()
    print(
        "WARNING:"
    )

    print(
        "FASTA counts differ from the previously audited "
        "169,637 total entries / 22,131 isoform headers."
    )


# ============================================================
# READ PRIMARY MBR-OFF PEPTIDES
# ============================================================

primary = pd.read_csv(
    PRIMARY_FILE,
    dtype=str
)


required = {
    "Program",
    "Sample",
    "Peptide"
}


missing = (
    required -
    set(primary.columns)
)


if missing:

    raise RuntimeError(
        f"Missing required columns: {missing}"
    )


primary = primary[
    [
        "Program",
        "Sample",
        "Peptide"
    ]
].copy()


primary["Program"] = (
    primary["Program"]
    .fillna("")
    .str.strip()
)


primary["Sample"] = (
    primary["Sample"]
    .fillna("")
    .str.strip()
)


primary["Peptide"] = (
    primary["Peptide"]
    .fillna("")
    .str.strip()
    .str.upper()
)


primary = primary[
    primary["Peptide"] != ""
]


primary = (
    primary
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


all_peptides = sorted(
    set(
        primary["Peptide"]
    )
)


valid_peptides = [
    x
    for x in all_peptides
    if VALID_PEPTIDE_RE.fullmatch(x)
]


invalid_peptides = sorted(
    set(all_peptides)
    -
    set(valid_peptides)
)


print()
print("=" * 110)
print("PRIMARY MBR-OFF PEPTIDES")
print("=" * 110)

print(
    "Program/sample observations :",
    len(primary)
)

print(
    "Union distinct peptides     :",
    len(all_peptides)
)

print(
    "Valid AA peptides           :",
    len(valid_peptides)
)

print(
    "Invalid peptide strings     :",
    len(invalid_peptides)
)


for program in PROGRAMS:

    n = primary.loc[
        primary["Program"] == program,
        "Peptide"
    ].nunique()

    print(
        f"{program:>3} distinct peptides        : {n:,}"
    )


# ============================================================
# AHO-CORASICK
# ============================================================

def build_automaton(
    patterns
):

    goto = [
        {}
    ]

    fail = [
        0
    ]

    output = [
        []
    ]


    # build trie

    for idx, pattern in enumerate(
        patterns
    ):

        state = 0

        for aa in pattern:

            if aa not in goto[state]:

                goto[state][aa] = len(
                    goto
                )

                goto.append(
                    {}
                )

                fail.append(
                    0
                )

                output.append(
                    []
                )

            state = goto[
                state
            ][aa]

        output[
            state
        ].append(
            idx
        )


    # failure links

    q = deque()


    for state in goto[0].values():

        fail[state] = 0
        q.append(state)


    while q:

        r = q.popleft()

        for aa, s in goto[r].items():

            q.append(
                s
            )

            f = fail[r]

            while (
                f != 0
                and
                aa not in goto[f]
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


def map_peptides(
    patterns,
    sequence_to_accessions,
    label
):

    print()
    print("=" * 110)
    print(label)
    print("=" * 110)

    t0 = time.time()

    goto, fail, output = (
        build_automaton(
            patterns
        )
    )


    print(
        "Peptide patterns  :",
        len(patterns)
    )

    print(
        "Automaton states  :",
        len(goto)
    )

    print(
        "Protein sequences :",
        len(sequence_to_accessions)
    )


    mapped = [
        set()
        for _ in patterns
    ]


    total = len(
        sequence_to_accessions
    )


    for n, (
        protein_sequence,
        accessions
    ) in enumerate(
        sequence_to_accessions.items(),
        1
    ):

        state = 0

        # avoid adding same peptide many times
        # when it appears repeatedly in one protein
        found = set()


        for aa in protein_sequence:

            while (
                state != 0
                and
                aa not in goto[state]
            ):

                state = fail[
                    state
                ]

            state = (
                goto[state]
                .get(
                    aa,
                    0
                )
            )


            if output[state]:

                found.update(
                    output[state]
                )


        for pep_idx in found:

            mapped[
                pep_idx
            ].update(
                accessions
            )


        if n % 10000 == 0:

            print(
                f"  scanned {n:,} / {total:,}"
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


# ============================================================
# EXACT REFERENCE
# ============================================================

exact_seq_to_accessions = defaultdict(
    set
)


for acc, seqs in accession_sequences.items():

    for seq in seqs:

        exact_seq_to_accessions[
            seq
        ].add(
            acc
        )


exact_mapped_raw = map_peptides(
    valid_peptides,
    exact_seq_to_accessions,
    "EXACT AMINO-ACID MAPPING"
)


exact_map = {

    peptide:
    exact_mapped_raw[i]

    for i, peptide
    in enumerate(
        valid_peptides
    )
}


# ============================================================
# I/L-EQUIVALENT REFERENCE
# ============================================================

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


il_patterns = sorted(
    normalized_to_original
)


il_seq_to_accessions = defaultdict(
    set
)


for acc, seqs in accession_sequences.items():

    for seq in seqs:

        il_seq_to_accessions[
            il_normalize(
                seq
            )
        ].add(
            acc
        )


il_mapped_raw = map_peptides(
    il_patterns,
    il_seq_to_accessions,
    "I/L-EQUIVALENT MAPPING — PRIMARY"
)


il_map = {}


for i, norm_peptide in enumerate(
    il_patterns
):

    accs = il_mapped_raw[
        i
    ]

    for original in normalized_to_original[
        norm_peptide
    ]:

        il_map[
            original
        ] = set(
            accs
        )


# ============================================================
# CLASSIFICATION
# ============================================================

def classify(
    accessions
):

    if not accessions:

        return (
            "unmapped",
            set(),
            set()
        )


    bases = {
        base_accession(x)
        for x in accessions
    }


    suffixes = {
        x
        for x in accessions
        if is_isoform(x)
    }


    # unique to one exact accession

    if len(accessions) == 1:

        acc = next(
            iter(accessions)
        )

        if is_isoform(acc):

            return (
                "single_isoform_unique",
                bases,
                suffixes
            )

        return (
            "single_canonical_unique",
            bases,
            suffixes
        )


    # maps across separate UniProt entries

    if len(bases) > 1:

        return (
            "cross_entry_shared",
            bases,
            suffixes
        )


    # all hits belong to one canonical/isoform family

    base = next(
        iter(bases)
    )


    complete_family = (
        family_members[
            base
        ]
    )


    if accessions == complete_family:

        return (
            "within_family_shared_all",
            bases,
            suffixes
        )


    return (
        "within_family_subset_discriminative",
        bases,
        suffixes
    )


# ============================================================
# PEPTIDE MAPPING TABLE
# ============================================================

mapping_rows = []


for peptide in all_peptides:

    if peptide in invalid_peptides:

        mapping_rows.append({
            "Peptide":
                peptide,

            "ExactCategory":
                "invalid_sequence",

            "ExactMappedAccessionCount":
                0,

            "ExactMappedSuffixAccessionCount":
                0,

            "ExactMappedSuffixAccessions":
                "",

            "ILCategory":
                "invalid_sequence",

            "ILMappedAccessionCount":
                0,

            "ILMappedBaseEntryCount":
                0,

            "ILMappedSuffixAccessionCount":
                0,

            "ILMappedAccessions":
                "",

            "ILMappedSuffixAccessions":
                ""
        })

        continue


    exact_acc = exact_map.get(
        peptide,
        set()
    )


    il_acc = il_map.get(
        peptide,
        set()
    )


    exact_category, exact_bases, exact_suffixes = (
        classify(
            exact_acc
        )
    )


    il_category, il_bases, il_suffixes = (
        classify(
            il_acc
        )
    )


    mapping_rows.append({
        "Peptide":
            peptide,

        "ExactCategory":
            exact_category,

        "ExactMappedAccessionCount":
            len(exact_acc),

        "ExactMappedSuffixAccessionCount":
            len(exact_suffixes),

        "ExactMappedSuffixAccessions":
            ";".join(
                sorted(
                    exact_suffixes
                )
            ),

        "ILCategory":
            il_category,

        "ILMappedAccessionCount":
            len(il_acc),

        "ILMappedBaseEntryCount":
            len(il_bases),

        "ILMappedSuffixAccessionCount":
            len(il_suffixes),

        "ILMappedAccessions":
            compact_list(
                il_acc
            ),

        "ILMappedSuffixAccessions":
            ";".join(
                sorted(
                    il_suffixes
                )
            )
    })


mapping_df = pd.DataFrame(
    mapping_rows
)


# ============================================================
# JOIN PROGRAM + SAMPLE
# ============================================================

detail = primary.merge(
    mapping_df,
    on="Peptide",
    how="left"
)


# ============================================================
# SUMMARY
# ============================================================

def accession_union(
    series
):

    result = set()

    for value in series:

        if pd.isna(value):
            continue

        for x in str(
            value
        ).split(";"):

            x = x.strip()

            if (
                x
                and
                not x.startswith(
                    "...TRUNCATED"
                )
            ):

                result.add(
                    x
                )

    return result


def summarize(
    category_col,
    suffix_col,
    mode
):

    rows = []


    for program in PROGRAMS:

        pdata = detail[
            detail["Program"] == program
        ]


        for sample in (
            SAMPLES
            +
            ["ALL_9_RUNS"]
        ):

            if sample == "ALL_9_RUNS":

                x = (
                    pdata
                    .drop_duplicates(
                        subset=[
                            "Peptide"
                        ]
                    )
                )

            else:

                x = pdata[
                    pdata["Sample"] == sample
                ]


            counts = (
                x[
                    category_col
                ]
                .value_counts()
                .to_dict()
            )


            single_isoforms = accession_union(

                x.loc[
                    x[category_col]
                    ==
                    "single_isoform_unique",

                    suffix_col
                ]
            )


            subset_isoforms = accession_union(

                x.loc[
                    x[category_col]
                    ==
                    "within_family_subset_discriminative",

                    suffix_col
                ]
            )


            discriminative_isoforms = (
                single_isoforms
                |
                subset_isoforms
            )


            single_n = counts.get(
                "single_isoform_unique",
                0
            )

            subset_n = counts.get(
                "within_family_subset_discriminative",
                0
            )


            rows.append({
                "MappingMode":
                    mode,

                "Program":
                    program,

                "Sample":
                    sample,

                "DistinctPeptides":
                    x["Peptide"].nunique(),

                "SingleIsoformUnique":
                    single_n,

                "SubsetDiscriminative":
                    subset_n,

                "TotalIsoformDiscriminativePeptides":
                    single_n + subset_n,

                "WithinFamilySharedAll":
                    counts.get(
                        "within_family_shared_all",
                        0
                    ),

                "SingleCanonicalUnique":
                    counts.get(
                        "single_canonical_unique",
                        0
                    ),

                "CrossEntryShared":
                    counts.get(
                        "cross_entry_shared",
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
                        discriminative_isoforms
                    )
            })


    return pd.DataFrame(
        rows
    )


summary_il = summarize(
    "ILCategory",
    "ILMappedSuffixAccessions",
    "IL_EQUIVALENT_PRIMARY"
)


summary_exact = summarize(
    "ExactCategory",
    "ExactMappedSuffixAccessions",
    "EXACT_SEQUENCE_SENSITIVITY"
)


# ============================================================
# EXACT vs I/L EFFECT
# ============================================================

reclass = (
    mapping_df[
        [
            "Peptide",
            "ExactCategory",
            "ILCategory"
        ]
    ]
    .copy()
)


reclass[
    "Changed"
] = (
    reclass[
        "ExactCategory"
    ]
    !=
    reclass[
        "ILCategory"
    ]
)


reclass_summary = (
    reclass
    .groupby(
        [
            "ExactCategory",
            "ILCategory"
        ],
        dropna=False
    )
    .size()
    .reset_index(
        name="PeptideCount"
    )
    .sort_values(
        "PeptideCount",
        ascending=False
    )
)


# ============================================================
# QC
# ============================================================

exact_unmapped = int(
    (
        mapping_df["ExactCategory"]
        ==
        "unmapped"
    ).sum()
)


il_unmapped = int(
    (
        mapping_df["ILCategory"]
        ==
        "unmapped"
    ).sum()
)


changed = int(
    reclass["Changed"].sum()
)


qc_df = pd.DataFrame(
    [
        {
            "Metric":
                "FASTA_Path",

            "Value":
                str(FASTA)
        },
        {
            "Metric":
                "FASTA_RawEntries",

            "Value":
                raw_count
        },
        {
            "Metric":
                "FASTA_TargetEntries",

            "Value":
                len(reference)
        },
        {
            "Metric":
                "FASTA_DistinctAccessions",

            "Value":
                len(all_accessions)
        },
        {
            "Metric":
                "FASTA_IsoformAccessions",

            "Value":
                len(isoform_accessions)
        },
        {
            "Metric":
                "FamiliesWithIsoforms",

            "Value":
                len(families_with_isoforms)
        },
        {
            "Metric":
                "UnionPrimaryPeptides",

            "Value":
                len(all_peptides)
        },
        {
            "Metric":
                "Exact_Unmapped",

            "Value":
                exact_unmapped
        },
        {
            "Metric":
                "IL_Unmapped",

            "Value":
                il_unmapped
        },
        {
            "Metric":
                "Exact_vs_IL_CategoryChanges",

            "Value":
                changed
        }
    ]
)


# ============================================================
# EXPORT
# ============================================================

outputs = {

    "01_FASTA_Mapping_QC.csv":
        qc_df,

    "02_CommonFASTA_PeptideMapping.csv":
        mapping_df,

    "03_PrimaryPeptide_IsoformDetail.csv":
        detail,

    "04_IsoformEvidenceSummary_IL_PRIMARY.csv":
        summary_il,

    "05_IsoformEvidenceSummary_EXACT_Sensitivity.csv":
        summary_exact,

    "06_Exact_vs_IL_Reclassification.csv":
        reclass_summary
}


for name, df in outputs.items():

    df.to_csv(
        OUT / name,
        index=False,
        encoding="utf-8-sig"
    )


# ============================================================
# PRINT RESULT
# ============================================================

show_cols = [
    "Program",
    "DistinctPeptides",
    "SingleIsoformUnique",
    "SubsetDiscriminative",
    "TotalIsoformDiscriminativePeptides",
    "WithinFamilySharedAll",
    "CrossEntryShared",
    "DistinctSingleIsoformAccessions",
    "DistinctSubsetSupportedIsoforms",
    "DistinctDiscriminativelySupportedIsoforms",
    "Unmapped"
]


print()
print("=" * 125)
print("PRIMARY RESULT — I/L-EQUIVALENT COMMON FASTA")
print("=" * 125)


print(
    summary_il.loc[
        summary_il["Sample"]
        ==
        "ALL_9_RUNS",
        show_cols
    ].to_string(
        index=False
    )
)


print()
print("=" * 125)
print("EXACT-SEQUENCE SENSITIVITY")
print("=" * 125)


print(
    summary_exact.loc[
        summary_exact["Sample"]
        ==
        "ALL_9_RUNS",
        show_cols
    ].to_string(
        index=False
    )
)


print()
print("=" * 125)
print("QC")
print("=" * 125)

print(
    qc_df.to_string(
        index=False
    )
)


print()
print(
    "Peptides changing category after I/L equivalence:",
    changed
)


print(
    "Total runtime:",
    round(
        time.time() - start_all,
        1
    ),
    "sec"
)


print()
print("OUTPUT DIRECTORY:")
print(OUT)

print()
print("STEP 1 COMPLETE")

