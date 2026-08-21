#!/usr/bin/env python3
"""
isoform_resolvability.py

Theoretical isoform-resolvability ceiling and complementary-protease rescue
for isoform-aware bottom-up proteomics.

Primary design choices:
1) UniProt isoform suffixes are preserved (e.g. P12345 and P12345-2 are distinct).
2) Peptide sequence keys can use I/L equivalence (I,L -> J), matching the manuscript.
3) The "accession family" is defined as the UniProt accession after removing a terminal -<integer> isoform suffix.
4) By default, compatibility is enzyme-generated compatibility. With --sequence-remap,
   candidate discriminative peptides are conservatively remapped by raw sequence containment
   across the entire normalized proteome using pyahocorasick.
5) Search constraints that are not yet verified in the manuscript are command-line parameters.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

SCRIPT_VERSION = "1.0.3"

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
NONSTANDARD_SPLIT_RE = re.compile(r"[^ACDEFGHIKLMNPQRSTVWY]+")
ISOFORM_SUFFIX_RE = re.compile(r"-(\d+)$")

# ----------------------------
# FASTA parsing / metadata
# ----------------------------

@dataclass(frozen=True)
class ProteinEntry:
    accession: str
    base_accession: str
    gene: Optional[str]
    sequence: str
    sequence_key: str
    is_suffixed_isoform: bool


def normalize_il(seq: str, il_equivalent: bool = True) -> str:
    seq = seq.upper()
    if il_equivalent:
        return seq.replace("I", "J").replace("L", "J")
    return seq


def strip_isoform_suffix(accession: str) -> str:
    return ISOFORM_SUFFIX_RE.sub("", accession)


def parse_uniprot_header(header: str) -> Tuple[str, Optional[str]]:
    """
    Supports standard UniProt FASTA headers such as:
      >sp|P12345-2|NAME_HUMAN ... GN=GENE ...
      >tr|A0A...|... GN=GENE ...
    Falls back to the first whitespace-delimited token when pipes are absent.
    """
    token = header.split()[0]
    if "|" in token:
        parts = token.split("|")
        accession = parts[1] if len(parts) >= 2 else token
    else:
        accession = token

    m = re.search(r"(?:^|\s)GN=([^\s]+)", header)
    gene = m.group(1) if m else None
    return accession, gene


def read_fasta(path: Path, il_equivalent: bool) -> List[ProteinEntry]:
    proteins: List[ProteinEntry] = []
    header: Optional[str] = None
    chunks: List[str] = []

    def flush():
        nonlocal header, chunks
        if header is None:
            return
        seq = "".join(chunks).replace(" ", "").replace("\r", "").upper()
        if not seq:
            header = None
            chunks = []
            return
        # UniProt entries may contain X (unknown), U (selenocysteine), or other
        # non-standard residue symbols. Keep the entry in the reference space, but
        # digest standard-amino-acid segments only. Non-standard residues therefore
        # act as hard sequence barriers: no theoretical peptide may contain or span
        # them. This is conservative and avoids inventing a definite peptide sequence.
        accession, gene = parse_uniprot_header(header)
        proteins.append(
            ProteinEntry(
                accession=accession,
                base_accession=strip_isoform_suffix(accession),
                gene=gene,
                sequence=seq,
                sequence_key=normalize_il(seq, il_equivalent),
                is_suffixed_isoform=bool(ISOFORM_SUFFIX_RE.search(accession)),
            )
        )
        header = None
        chunks = []

    with path.open("rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
            elif header is not None:
                chunks.append(line.strip())
        flush()

    # Accession uniqueness is important because exact entry resolution assumes distinct entries.
    accessions = [p.accession for p in proteins]
    dup = [a for a, n in Counter(accessions).items() if n > 1]
    if dup:
        raise ValueError(
            f"Duplicate accession(s) found in FASTA, e.g. {dup[:5]}. "
            "Deduplicate or rename before running isoform-resolvability analysis."
        )
    return proteins


def summarize_nonstandard_residues(proteins: Sequence[ProteinEntry]) -> Tuple[int, Dict[str, int]]:
    affected_entries = 0
    counts: Counter = Counter()
    for p in proteins:
        bad = [aa for aa in p.sequence if aa not in STANDARD_AA]
        if bad:
            affected_entries += 1
            counts.update(bad)
    return affected_entries, dict(sorted(counts.items()))


def sha256sum(path: Path, blocksize: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            b = fh.read(blocksize)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def reference_space_universe(proteins: Sequence[ProteinEntry]) -> dict:
    """Summarize accession-family and annotated-gene denominators."""
    family_to_accessions: Dict[str, Set[str]] = defaultdict(set)
    gene_to_accessions: Dict[str, Set[str]] = defaultdict(set)
    gene_to_families: Dict[str, Set[str]] = defaultdict(set)

    for p in proteins:
        family_to_accessions[p.base_accession].add(p.accession)
        if p.gene:
            gene_to_accessions[p.gene].add(p.accession)
            gene_to_families[p.gene].add(p.base_accession)

    multi_entry_families = {
        fam for fam, accs in family_to_accessions.items() if len(accs) > 1
    }
    multi_entry_genes = {
        gene for gene, accs in gene_to_accessions.items() if len(accs) > 1
    }
    multi_family_genes = {
        gene for gene, fams in gene_to_families.items() if len(fams) > 1
    }

    return {
        "n_total_accession_families": len(family_to_accessions),
        "n_multi_entry_accession_families": len(multi_entry_families),
        "n_annotated_genes": len(gene_to_accessions),
        "n_multi_entry_annotated_genes": len(multi_entry_genes),
        "n_multi_family_annotated_genes": len(multi_family_genes),
        "multi_entry_families": multi_entry_families,
        "multi_entry_genes": multi_entry_genes,
    }


def write_reference_space_and_rescue_rates(
    proteins: Sequence[ProteinEntry],
    results: Dict[str, "ProteaseResult"],
    outdir: Path,
    baseline: str = "trypsin",
) -> None:
    """Write denominators, resolution rates, rescue rates, and protease-combination coverage."""
    u = reference_space_universe(proteins)

    ref_path = outdir / "reference_space_summary.tsv"
    with ref_path.open("wt", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        for key in [
            "n_total_accession_families",
            "n_multi_entry_accession_families",
            "n_annotated_genes",
            "n_multi_entry_annotated_genes",
            "n_multi_family_annotated_genes",
        ]:
            w.writerow([key, u[key]])

    fam_den = u["n_multi_entry_accession_families"]
    gene_den = u["n_multi_entry_annotated_genes"]

    rate_path = outdir / "protease_resolution_rates.tsv"
    with rate_path.open("wt", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "protease",
            "primary_discriminative_families",
            "primary_discriminative_family_pct_of_multi_entry_universe",
            "exact_entry_families",
            "exact_entry_family_pct_of_multi_entry_universe",
            "primary_discriminative_genes",
            "primary_discriminative_gene_pct_of_multi_entry_annotated_universe",
            "exact_entry_genes",
            "exact_entry_gene_pct_of_multi_entry_annotated_universe",
        ])
        for name, r in results.items():
            w.writerow([
                name,
                len(r.primary_discriminative_families),
                100.0 * len(r.primary_discriminative_families) / fam_den if fam_den else math.nan,
                len(r.exact_entry_families),
                100.0 * len(r.exact_entry_families) / fam_den if fam_den else math.nan,
                len(r.primary_discriminative_genes),
                100.0 * len(r.primary_discriminative_genes) / gene_den if gene_den else math.nan,
                len(r.exact_entry_genes),
                100.0 * len(r.exact_entry_genes) / gene_den if gene_den else math.nan,
            ])

    if baseline in results:
        b = results[baseline]
        unresolved_primary_families = u["multi_entry_families"] - b.primary_discriminative_families
        unresolved_exact_families = u["multi_entry_families"] - b.exact_entry_families
        unresolved_primary_genes = u["multi_entry_genes"] - b.primary_discriminative_genes
        unresolved_exact_genes = u["multi_entry_genes"] - b.exact_entry_genes

        rescue_rate_path = outdir / "protease_rescue_rates_vs_trypsin.tsv"
        with rescue_rate_path.open("wt", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow([
                "protease",
                "trypsin_unresolved_primary_families",
                "rescued_primary_families",
                "rescued_primary_family_pct_of_trypsin_unresolved",
                "trypsin_unresolved_exact_families",
                "rescued_exact_families",
                "rescued_exact_family_pct_of_trypsin_unresolved",
                "trypsin_unresolved_primary_genes",
                "rescued_primary_genes",
                "rescued_primary_gene_pct_of_trypsin_unresolved",
                "trypsin_unresolved_exact_genes",
                "rescued_exact_genes",
                "rescued_exact_gene_pct_of_trypsin_unresolved",
            ])
            for name, r in results.items():
                if name == baseline:
                    continue
                rf = len(r.primary_discriminative_families & unresolved_primary_families)
                ref = len(r.exact_entry_families & unresolved_exact_families)
                rg = len(r.primary_discriminative_genes & unresolved_primary_genes)
                reg = len(r.exact_entry_genes & unresolved_exact_genes)
                w.writerow([
                    name,
                    len(unresolved_primary_families),
                    rf,
                    100.0 * rf / len(unresolved_primary_families) if unresolved_primary_families else math.nan,
                    len(unresolved_exact_families),
                    ref,
                    100.0 * ref / len(unresolved_exact_families) if unresolved_exact_families else math.nan,
                    len(unresolved_primary_genes),
                    rg,
                    100.0 * rg / len(unresolved_primary_genes) if unresolved_primary_genes else math.nan,
                    len(unresolved_exact_genes),
                    reg,
                    100.0 * reg / len(unresolved_exact_genes) if unresolved_exact_genes else math.nan,
                ])

    combo_path = outdir / "protease_combination_coverage.tsv"
    names = list(results)
    combos = []
    if baseline in results:
        for name in names:
            if name != baseline:
                combos.append((f"{baseline}+{name}", [baseline, name]))
    if len(names) > 1:
        combos.append(("all_proteases", names))

    with combo_path.open("wt", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "combination",
            "n_proteases",
            "primary_discriminative_families_union",
            "primary_discriminative_family_pct_of_multi_entry_universe",
            "exact_entry_families_union",
            "exact_entry_family_pct_of_multi_entry_universe",
            "primary_discriminative_genes_union",
            "primary_discriminative_gene_pct_of_multi_entry_annotated_universe",
            "exact_entry_genes_union",
            "exact_entry_gene_pct_of_multi_entry_annotated_universe",
        ])
        for label, members in combos:
            pf = set().union(*(results[n].primary_discriminative_families for n in members))
            ef = set().union(*(results[n].exact_entry_families for n in members))
            pg = set().union(*(results[n].primary_discriminative_genes for n in members))
            eg = set().union(*(results[n].exact_entry_genes for n in members))
            w.writerow([
                label,
                len(members),
                len(pf),
                100.0 * len(pf) / fam_den if fam_den else math.nan,
                len(ef),
                100.0 * len(ef) / fam_den if fam_den else math.nan,
                len(pg),
                100.0 * len(pg) / gene_den if gene_den else math.nan,
                len(eg),
                100.0 * len(eg) / gene_den if gene_den else math.nan,
            ])


# ----------------------------
# Protease definitions
# ----------------------------

def _not_before_proline(seq: str, i: int) -> bool:
    return i + 1 >= len(seq) or seq[i + 1] != "P"


def cleave_trypsin(seq: str, i: int) -> bool:
    return seq[i] in {"K", "R"} and _not_before_proline(seq, i)


def cleave_lysc(seq: str, i: int) -> bool:
    # Lys-C: C-terminal cleavage at K. Kept explicit rather than silently applying trypsin's P rule.
    return seq[i] == "K"


def cleave_argc(seq: str, i: int) -> bool:
    return seq[i] == "R"


def cleave_gluc_e(seq: str, i: int) -> bool:
    # Glu-C E-specific mode. Buffer-dependent D cleavage is a separate explicit option.
    return seq[i] == "E"


def cleave_gluc_de(seq: str, i: int) -> bool:
    return seq[i] in {"D", "E"}


def cleave_chymo_fyw(seq: str, i: int) -> bool:
    # High-specificity conservative chymotrypsin model: F/Y/W, not before P.
    return seq[i] in {"F", "Y", "W"} and _not_before_proline(seq, i)


PROTEASES: Dict[str, Callable[[str, int], bool]] = {
    "trypsin": cleave_trypsin,
    "lysc": cleave_lysc,
    "argc": cleave_argc,
    "gluc_e": cleave_gluc_e,
    "gluc_de": cleave_gluc_de,
    "chymo_fyw": cleave_chymo_fyw,
}


def cleavage_boundaries(seq: str, cleavage_rule: Callable[[str, int], bool]) -> List[int]:
    boundaries = [0]
    for i in range(len(seq)):
        if cleavage_rule(seq, i):
            boundaries.append(i + 1)
    if boundaries[-1] != len(seq):
        boundaries.append(len(seq))
    # Avoid duplicate terminal boundary if cleavage occurs at the final residue.
    out = [boundaries[0]]
    for b in boundaries[1:]:
        if b != out[-1]:
            out.append(b)
    return out


def digest_sequence(
    seq: str,
    cleavage_rule: Callable[[str, int], bool],
    max_missed_cleavages: int,
    min_length: int,
    max_length: int,
) -> Iterator[str]:
    # Treat non-standard residues as hard barriers. Each standard-amino-acid segment
    # is digested independently, so no peptide can contain or span X/U/O/etc.
    for segment in NONSTANDARD_SPLIT_RE.split(seq):
        if not segment:
            continue
        b = cleavage_boundaries(segment, cleavage_rule)
        for start_idx in range(len(b) - 1):
            for missed in range(max_missed_cleavages + 1):
                end_idx = start_idx + missed + 1
                if end_idx >= len(b):
                    break
                peptide = segment[b[start_idx] : b[end_idx]]
                n = len(peptide)
                if n < min_length:
                    continue
                if n > max_length:
                    # Increasing missed cleavages only increases peptide length.
                    break
                yield peptide


# ----------------------------
# Structural classification
# ----------------------------

@dataclass
class Classification:
    category: str
    resolution_level: str
    n_accessions: int
    n_base_accessions: int
    n_genes: int
    target_base: Optional[str]
    target_gene: Optional[str]
    exact_entry: bool
    primary_isoform_discriminative: bool
    structural_discriminative: bool


def classify_mapping(
    entry_indices: Set[int],
    proteins: Sequence[ProteinEntry],
    base_totals: Dict[str, int],
    gene_totals: Dict[str, int],
) -> Classification:
    if not entry_indices:
        return Classification(
            "unmapped", "unmapped", 0, 0, 0, None, None, False, False, False
        )

    entries = [proteins[i] for i in entry_indices]
    bases = {p.base_accession for p in entries}

    # Missing GN is treated as unresolved rather than pooling all missing-GN entries.
    real_genes = {p.gene for p in entries if p.gene}
    has_missing_gene = any(p.gene is None for p in entries)
    n_genes_effective = len(real_genes) + (1 if has_missing_gene else 0)

    if len(entries) == 1:
        p = entries[0]
        if p.is_suffixed_isoform:
            return Classification(
                "single_isoform_unique",
                "exact_suffixed_isoform_entry",
                1, 1, 1 if p.gene else 0,
                p.base_accession, p.gene,
                True, True, True,
            )
        return Classification(
            "single_canonical_unique",
            "exact_unsuffixed_entry",
            1, 1, 1 if p.gene else 0,
            p.base_accession, p.gene,
            True, False, True,
        )

    if len(bases) == 1:
        base = next(iter(bases))
        gene = next(iter(real_genes)) if len(real_genes) == 1 and not has_missing_gene else None
        if len(entries) < base_totals[base]:
            return Classification(
                "within_family_subset_discriminative",
                "restricted_accession_family_subset",
                len(entries), 1, n_genes_effective,
                base, gene,
                False, True, True,
            )
        return Classification(
            "within_family_shared_all",
            "shared_across_accession_family",
            len(entries), 1, n_genes_effective,
            base, gene,
            False, False, False,
        )

    if len(real_genes) == 1 and not has_missing_gene:
        gene = next(iter(real_genes))
        if len(entries) < gene_totals[gene]:
            return Classification(
                "same_gene_subset_discriminative",
                "restricted_same_gene_multi_accession_subset",
                len(entries), len(bases), 1,
                None, gene,
                False, False, True,
            )
        return Classification(
            "same_gene_multi_entry_shared",
            "shared_across_all_considered_gene_entries",
            len(entries), len(bases), 1,
            None, gene,
            False, False, False,
        )

    if len(real_genes) >= 2:
        return Classification(
            "cross_gene_shared",
            "cross_gene_ambiguous",
            len(entries), len(bases), n_genes_effective,
            None, None,
            False, False, False,
        )

    return Classification(
        "multi_entry_gene_unresolved",
        "gene_annotation_unresolved",
        len(entries), len(bases), n_genes_effective,
        None, None,
        False, False, False,
    )


# ----------------------------
# Optional conservative sequence remapping
# ----------------------------

def sequence_remap_candidates(
    candidate_peptides: Sequence[str],
    proteins: Sequence[ProteinEntry],
    batch_size: int,
) -> Dict[str, Set[int]]:
    """
    Remap candidate peptide keys by sequence containment across the full I/L-normalized proteome.
    This only needs to be done for initially discriminative peptides because adding compatible
    entries cannot turn an already-shared peptide into a more discriminative peptide.
    """
    try:
        import ahocorasick  # type: ignore
    except ImportError as e:
        raise RuntimeError(
            "--sequence-remap requires the free 'pyahocorasick' package.\n"
            "Install it with: pip install pyahocorasick"
        ) from e

    remapped: Dict[str, Set[int]] = {}
    candidate_peptides = list(candidate_peptides)

    for start in range(0, len(candidate_peptides), batch_size):
        batch = candidate_peptides[start : start + batch_size]
        A = ahocorasick.Automaton()
        for pep in batch:
            A.add_word(pep, pep)
        A.make_automaton()

        local: Dict[str, Set[int]] = {pep: set() for pep in batch}
        for idx, p in enumerate(proteins):
            for _, pep in A.iter(p.sequence_key):
                local[pep].add(idx)

        remapped.update(local)
        del A
        del local

    return remapped


# ----------------------------
# One-protease analysis
# ----------------------------

@dataclass
class ProteaseResult:
    protease: str
    n_theoretical_peptides: int
    class_counts: Counter
    n_exact_entry_peptides: int
    n_primary_discriminative_peptides: int
    n_structural_discriminative_peptides: int
    exact_entry_peptides: Set[str]
    primary_discriminative_peptides: Set[str]
    structural_discriminative_peptides: Set[str]
    exact_entry_families: Set[str]
    primary_discriminative_families: Set[str]
    structural_discriminative_families: Set[str]
    exact_entry_genes: Set[str]
    primary_discriminative_genes: Set[str]
    structural_discriminative_genes: Set[str]
    peptide_to_class: Dict[str, Classification]


def analyze_protease(
    name: str,
    proteins: Sequence[ProteinEntry],
    max_missed_cleavages: int,
    min_length: int,
    max_length: int,
    il_equivalent: bool,
    outdir: Path,
    sequence_remap: bool,
    remap_batch_size: int,
    keep_peptide_details: bool,
) -> ProteaseResult:
    rule = PROTEASES[name]

    # peptide key -> compatible entries under the enzyme-generated theoretical space
    pep_to_entries: Dict[str, Set[int]] = defaultdict(set)

    for idx, p in enumerate(proteins):
        # Deduplicate repeated sequence occurrences inside one protein entry.
        per_protein: Set[str] = set()
        for pep in digest_sequence(
            p.sequence, rule, max_missed_cleavages, min_length, max_length
        ):
            per_protein.add(normalize_il(pep, il_equivalent))
        for key in per_protein:
            pep_to_entries[key].add(idx)

    base_totals = Counter(p.base_accession for p in proteins)
    gene_totals = Counter(p.gene for p in proteins if p.gene)

    candidates: List[str] = []

    for pep, idxs in pep_to_entries.items():
        c = classify_mapping(idxs, proteins, base_totals, gene_totals)
        # Sequence-remap only isoform-relevant discriminative candidates. Remapping can only
        # add compatible entries, so already-shared peptides cannot become more discriminative.
        isoform_relevant_exact = False
        if c.exact_entry and len(idxs) == 1:
            p0 = proteins[next(iter(idxs))]
            isoform_relevant_exact = (
                base_totals[p0.base_accession] > 1
                or (p0.gene is not None and gene_totals[p0.gene] > 1)
            )
        if c.primary_isoform_discriminative or c.category == "same_gene_subset_discriminative" or isoform_relevant_exact:
            candidates.append(pep)

    overrides: Dict[str, Set[int]] = {}
    if sequence_remap and candidates:
        overrides = sequence_remap_candidates(candidates, proteins, remap_batch_size)

    class_counts: Counter = Counter()
    exact_entry_peptides: Set[str] = set()
    primary_peptides: Set[str] = set()
    structural_peptides: Set[str] = set()

    exact_families: Set[str] = set()
    primary_families: Set[str] = set()
    structural_families: Set[str] = set()

    exact_genes: Set[str] = set()
    primary_genes: Set[str] = set()
    structural_genes: Set[str] = set()

    final_classes: Dict[str, Classification] = {} if keep_peptide_details else {}

    catalog_path = outdir / f"peptide_catalog_{name}.tsv.gz"
    with gzip.open(catalog_path, "wt", newline="", encoding="utf-8") as gz:
        writer = csv.writer(gz, delimiter="\t")
        writer.writerow([
            "peptide_key",
            "category",
            "resolution_level",
            "n_compatible_accessions",
            "n_base_accessions",
            "n_genes",
            "target_base_accession",
            "target_gene",
            "exact_entry",
            "primary_isoform_discriminative",
            "structural_discriminative",
            "mapping_mode",
            "compatible_accessions_if_discriminative",
        ])

        for pep in sorted(pep_to_entries):
            idxs = overrides.get(pep, pep_to_entries[pep])
            c = classify_mapping(idxs, proteins, base_totals, gene_totals)
            if keep_peptide_details:
                final_classes[pep] = c
            class_counts[c.category] += 1

            if c.exact_entry:
                if keep_peptide_details:
                    exact_entry_peptides.add(pep)
            if c.primary_isoform_discriminative:
                primary_peptides.add(pep)
            if c.structural_discriminative and keep_peptide_details:
                structural_peptides.add(pep)

            # Family/gene resolvability is meaningful only when >1 reference entry exists.
            if c.target_base and base_totals[c.target_base] > 1:
                if c.exact_entry:
                    exact_families.add(c.target_base)
                if c.primary_isoform_discriminative:
                    primary_families.add(c.target_base)
                if c.structural_discriminative:
                    structural_families.add(c.target_base)

            if c.target_gene and gene_totals[c.target_gene] > 1:
                if c.exact_entry:
                    exact_genes.add(c.target_gene)
                if c.primary_isoform_discriminative:
                    primary_genes.add(c.target_gene)
                if c.structural_discriminative:
                    structural_genes.add(c.target_gene)

            writer.writerow([
                pep,
                c.category,
                c.resolution_level,
                c.n_accessions,
                c.n_base_accessions,
                c.n_genes,
                c.target_base or "",
                c.target_gene or "",
                int(c.exact_entry),
                int(c.primary_isoform_discriminative),
                int(c.structural_discriminative),
                "sequence_containment" if pep in overrides else "enzyme_generated",
                ";".join(sorted(proteins[i].accession for i in idxs)) if c.structural_discriminative else "",
            ])

    n_exact_count = sum(v for k, v in class_counts.items() if k in {"single_isoform_unique", "single_canonical_unique"})
    n_primary_count = class_counts.get("single_isoform_unique", 0) + class_counts.get("within_family_subset_discriminative", 0)
    n_structural_count = (
        n_exact_count
        + class_counts.get("within_family_subset_discriminative", 0)
        + class_counts.get("same_gene_subset_discriminative", 0)
    )

    return ProteaseResult(
        protease=name,
        n_theoretical_peptides=len(pep_to_entries),
        class_counts=class_counts,
        n_exact_entry_peptides=n_exact_count,
        n_primary_discriminative_peptides=n_primary_count,
        n_structural_discriminative_peptides=n_structural_count,
        exact_entry_peptides=exact_entry_peptides,
        primary_discriminative_peptides=primary_peptides,
        structural_discriminative_peptides=structural_peptides,
        exact_entry_families=exact_families,
        primary_discriminative_families=primary_families,
        structural_discriminative_families=structural_families,
        exact_entry_genes=exact_genes,
        primary_discriminative_genes=primary_genes,
        structural_discriminative_genes=structural_genes,
        peptide_to_class=final_classes,
    )


# ----------------------------
# Observed peptide integration
# ----------------------------

def read_observed(path: Path, il_equivalent: bool) -> List[dict]:
    rows = []
    with path.open("rt", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"workflow", "peptide_sequence"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(
                f"Observed TSV must contain columns {sorted(required)}; found {reader.fieldnames}"
            )
        for row in reader:
            seq = (row.get("peptide_sequence") or "").strip().upper()
            if not seq:
                continue
            row["peptide_key"] = normalize_il(seq, il_equivalent)
            rows.append(row)
    return rows


def write_observed_recovery(
    observed_rows: List[dict],
    trypsin_result: ProteaseResult,
    outdir: Path,
) -> None:
    workflow_to_peps: Dict[str, Set[str]] = defaultdict(set)
    for row in observed_rows:
        workflow_to_peps[row["workflow"]].add(row["peptide_key"])

    theo_all = set(trypsin_result.peptide_to_class)
    theo_exact = trypsin_result.exact_entry_peptides
    theo_primary = trypsin_result.primary_discriminative_peptides
    theo_struct = trypsin_result.structural_discriminative_peptides

    out = outdir / "observed_recovery_by_workflow.tsv"
    with out.open("wt", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "workflow",
            "observed_unique_peptides",
            "observed_in_theoretical_tryptic_space",
            "observed_exact_entry",
            "observed_primary_isoform_discriminative",
            "observed_structural_discriminative",
            "theoretical_exact_entry_total",
            "theoretical_primary_isoform_discriminative_total",
            "theoretical_structural_discriminative_total",
            "recovery_exact_entry_pct",
            "recovery_primary_isoform_discriminative_pct",
            "recovery_structural_discriminative_pct",
        ])
        for wf in sorted(workflow_to_peps):
            obs = workflow_to_peps[wf]
            n_exact = len(obs & theo_exact)
            n_primary = len(obs & theo_primary)
            n_struct = len(obs & theo_struct)
            w.writerow([
                wf,
                len(obs),
                len(obs & theo_all),
                n_exact,
                n_primary,
                n_struct,
                len(theo_exact),
                len(theo_primary),
                len(theo_struct),
                100.0 * n_exact / len(theo_exact) if theo_exact else math.nan,
                100.0 * n_primary / len(theo_primary) if theo_primary else math.nan,
                100.0 * n_struct / len(theo_struct) if theo_struct else math.nan,
            ])

    # Cross-workflow support for theoretical primary-discriminative peptides.
    support = Counter()
    for wf, peps in workflow_to_peps.items():
        for pep in peps & theo_primary:
            support[pep] += 1

    out2 = outdir / "observed_primary_discriminative_support.tsv"
    with out2.open("wt", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["peptide_key", "n_supporting_workflows", "category"])
        for pep in sorted(theo_primary):
            w.writerow([
                pep,
                support.get(pep, 0),
                trypsin_result.peptide_to_class[pep].category,
            ])


# ----------------------------
# Rescue summaries / plots
# ----------------------------

def write_summary(results: Dict[str, ProteaseResult], outdir: Path) -> None:
    summary_path = outdir / "protease_summary.tsv"
    all_categories = sorted(
        set().union(*(set(r.class_counts) for r in results.values()))
    )
    with summary_path.open("wt", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        header = [
            "protease",
            "n_theoretical_unique_peptide_keys",
            "n_exact_entry_peptides",
            "n_primary_isoform_discriminative_peptides",
            "n_structural_discriminative_peptides",
            "n_exact_entry_families",
            "n_primary_discriminative_families",
            "n_structural_discriminative_families",
            "n_exact_entry_genes",
            "n_primary_discriminative_genes",
            "n_structural_discriminative_genes",
        ] + [f"class__{c}" for c in all_categories]
        w.writerow(header)
        for name, r in results.items():
            w.writerow([
                name,
                r.n_theoretical_peptides,
                r.n_exact_entry_peptides,
                r.n_primary_discriminative_peptides,
                r.n_structural_discriminative_peptides,
                len(r.exact_entry_families),
                len(r.primary_discriminative_families),
                len(r.structural_discriminative_families),
                len(r.exact_entry_genes),
                len(r.primary_discriminative_genes),
                len(r.structural_discriminative_genes),
            ] + [r.class_counts.get(c, 0) for c in all_categories])


def write_rescue(results: Dict[str, ProteaseResult], baseline: str, outdir: Path) -> None:
    if baseline not in results:
        return
    b = results[baseline]
    path = outdir / "protease_rescue_vs_trypsin.tsv"
    with path.open("wt", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "protease",
            "new_primary_discriminative_peptides_vs_trypsin",
            "rescued_primary_discriminative_families",
            "rescued_exact_entry_families",
            "rescued_structural_discriminative_families",
            "rescued_primary_discriminative_genes",
            "rescued_exact_entry_genes",
            "rescued_structural_discriminative_genes",
            "union_primary_discriminative_families",
            "union_exact_entry_families",
            "union_primary_discriminative_genes",
            "union_exact_entry_genes",
        ])
        for name, r in results.items():
            if name == baseline:
                continue
            w.writerow([
                name,
                len(r.primary_discriminative_peptides - b.primary_discriminative_peptides),
                len(r.primary_discriminative_families - b.primary_discriminative_families),
                len(r.exact_entry_families - b.exact_entry_families),
                len(r.structural_discriminative_families - b.structural_discriminative_families),
                len(r.primary_discriminative_genes - b.primary_discriminative_genes),
                len(r.exact_entry_genes - b.exact_entry_genes),
                len(r.structural_discriminative_genes - b.structural_discriminative_genes),
                len(r.primary_discriminative_families | b.primary_discriminative_families),
                len(r.exact_entry_families | b.exact_entry_families),
                len(r.primary_discriminative_genes | b.primary_discriminative_genes),
                len(r.exact_entry_genes | b.exact_entry_genes),
            ])

    # Per-family rescue matrix.
    all_families = set().union(*(r.structural_discriminative_families for r in results.values()))
    matrix = outdir / "protease_rescue_family_matrix.tsv"
    with matrix.open("wt", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        names = list(results)
        w.writerow(
            ["base_accession"]
            + [f"{n}__exact_entry" for n in names]
            + [f"{n}__primary_discriminative" for n in names]
            + [f"{n}__structural_discriminative" for n in names]
        )
        for fam in sorted(all_families):
            w.writerow(
                [fam]
                + [int(fam in results[n].exact_entry_families) for n in names]
                + [int(fam in results[n].primary_discriminative_families) for n in names]
                + [int(fam in results[n].structural_discriminative_families) for n in names]
            )


def make_plots(results: Dict[str, ProteaseResult], outdir: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; skipping plots.", file=sys.stderr)
        return

    names = list(results)
    primary = [results[n].n_primary_discriminative_peptides for n in names]
    exact = [results[n].n_exact_entry_peptides for n in names]
    structural = [results[n].n_structural_discriminative_peptides for n in names]

    x = list(range(len(names)))
    width = 0.26
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar([v - width for v in x], exact, width=width, label="Exact-entry")
    ax.bar(x, primary, width=width, label="Primary isoform-discriminative")
    ax.bar([v + width for v in x], structural, width=width, label="Structural-discriminative")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=25, ha="right")
    ax.set_ylabel("Theoretical unique peptide keys")
    ax.set_title("Theoretical isoform-resolvability by protease")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "Fig_theoretical_resolvability_by_protease.png", dpi=300)
    fig.savefig(outdir / "Fig_theoretical_resolvability_by_protease.pdf")
    plt.close(fig)

    if "trypsin" in results and len(results) > 1:
        b = results["trypsin"]
        alts = [n for n in names if n != "trypsin"]
        rescued_family = [
            len(results[n].primary_discriminative_families - b.primary_discriminative_families)
            for n in alts
        ]
        rescued_gene = [
            len(results[n].primary_discriminative_genes - b.primary_discriminative_genes)
            for n in alts
        ]
        x = list(range(len(alts)))
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.bar([v - 0.2 for v in x], rescued_family, width=0.4, label="Accession families")
        ax.bar([v + 0.2 for v in x], rescued_gene, width=0.4, label="Genes")
        ax.set_xticks(x)
        ax.set_xticklabels(alts, rotation=25, ha="right")
        ax.set_ylabel("Resolved by alternative protease, not trypsin")
        ax.set_title("Complementary-protease rescue of trypsin-unresolved space")
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / "Fig_alternative_protease_rescue.png", dpi=300)
        fig.savefig(outdir / "Fig_alternative_protease_rescue.pdf")
        plt.close(fig)


# ----------------------------
# CLI
# ----------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Compute theoretical isoform-resolvability and alternative-protease rescue."
    )
    p.add_argument("--fasta", required=True, type=Path, help="Isoform-inclusive UniProt FASTA.")
    p.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    p.add_argument(
        "--proteases",
        nargs="+",
        default=["trypsin", "lysc", "gluc_e", "chymo_fyw"],
        choices=sorted(PROTEASES),
        help="Proteases to evaluate."
    )
    p.add_argument(
        "--missed-cleavages",
        required=True,
        type=int,
        help="Maximum missed cleavages. Use the verified search setting from the final logs."
    )
    p.add_argument(
        "--min-length",
        required=True,
        type=int,
        help="Minimum peptide length. Use the verified search setting."
    )
    p.add_argument(
        "--max-length",
        required=True,
        type=int,
        help="Maximum peptide length. Use the verified search setting."
    )
    p.add_argument(
        "--no-il-equivalence",
        action="store_true",
        help="Do not collapse I/L. The manuscript primary mapping uses I/L equivalence, so normally leave this off."
    )
    p.add_argument(
        "--sequence-remap",
        action="store_true",
        help=(
            "Conservatively remap initially discriminative peptide keys by sequence containment "
            "across the entire proteome. Requires pyahocorasick."
        )
    )
    p.add_argument(
        "--remap-batch-size",
        type=int,
        default=100000,
        help="Pattern batch size for --sequence-remap."
    )
    p.add_argument(
        "--observed",
        type=Path,
        default=None,
        help=(
            "Optional tab-delimited observed peptide table with columns workflow and peptide_sequence. "
            "Used to calculate empirical recovery of the theoretical tryptic space."
        )
    )
    return p


def main() -> None:
    args = build_parser().parse_args()
    if args.missed_cleavages < 0:
        raise SystemExit("--missed-cleavages must be >= 0")
    if args.min_length < 1 or args.max_length < args.min_length:
        raise SystemExit("Require 1 <= --min-length <= --max-length")
    if "trypsin" not in args.proteases and args.observed:
        raise SystemExit("--observed recovery requires trypsin in --proteases")

    args.outdir.mkdir(parents=True, exist_ok=True)
    il_equivalent = not args.no_il_equivalence

    print(f"[1/4] Reading FASTA: {args.fasta}", file=sys.stderr)
    proteins = read_fasta(args.fasta, il_equivalent)
    print(f"      Protein entries: {len(proteins):,}", file=sys.stderr)
    n_nonstandard_entries, nonstandard_counts = summarize_nonstandard_residues(proteins)
    if n_nonstandard_entries:
        formatted = ", ".join(f"{k}:{v:,}" for k, v in nonstandard_counts.items())
        print(
            f"      Non-standard residues retained as hard barriers in "
            f"{n_nonstandard_entries:,} entries ({formatted})",
            file=sys.stderr,
        )

    manifest = {
        "script_version": SCRIPT_VERSION,
        "fasta": str(args.fasta.resolve()),
        "fasta_sha256": sha256sum(args.fasta),
        "n_protein_entries": len(proteins),
        "n_suffixed_isoform_entries": sum(p.is_suffixed_isoform for p in proteins),
        "n_entries_with_gene_name": sum(p.gene is not None for p in proteins),
        "n_entries_with_nonstandard_residues": n_nonstandard_entries,
        "nonstandard_residue_counts": nonstandard_counts,
        "nonstandard_residue_handling": (
            "Residues outside the 20 standard amino acids are retained in FASTA entries "
            "but treated as hard barriers during in-silico digestion; no theoretical peptide "
            "may contain or span them."
        ),
        "il_equivalent": il_equivalent,
        "proteases": args.proteases,
        "max_missed_cleavages": args.missed_cleavages,
        "min_peptide_length": args.min_length,
        "max_peptide_length": args.max_length,
        "sequence_remap": args.sequence_remap,
        "sequence_remap_batch_size": args.remap_batch_size if args.sequence_remap else None,
        "classification_note": (
            "Primary isoform-discriminative = single_isoform_unique + "
            "within_family_subset_discriminative, matching the manuscript's primary definition. "
            "single_canonical_unique is retained separately as exact-entry evidence."
        ),
        "protease_rule_note": {
            "trypsin": "cleave C-terminal K/R, not before P",
            "lysc": "cleave C-terminal K",
            "argc": "cleave C-terminal R",
            "gluc_e": "cleave C-terminal E",
            "gluc_de": "cleave C-terminal D/E",
            "chymo_fyw": "conservative high-specificity model: cleave C-terminal F/Y/W, not before P",
        }
    }
    with (args.outdir / "analysis_manifest.json").open("wt", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)

    results: Dict[str, ProteaseResult] = {}
    for i, name in enumerate(args.proteases, start=1):
        print(
            f"[2/4] Protease {i}/{len(args.proteases)}: {name}",
            file=sys.stderr
        )
        results[name] = analyze_protease(
            name=name,
            proteins=proteins,
            max_missed_cleavages=args.missed_cleavages,
            min_length=args.min_length,
            max_length=args.max_length,
            il_equivalent=il_equivalent,
            outdir=args.outdir,
            sequence_remap=args.sequence_remap,
            remap_batch_size=args.remap_batch_size,
            keep_peptide_details=(name == "trypsin" and args.observed is not None),
        )

    print("[3/4] Writing summary/rescue tables and figures", file=sys.stderr)
    write_summary(results, args.outdir)
    write_rescue(results, "trypsin", args.outdir)
    write_reference_space_and_rescue_rates(proteins, results, args.outdir, baseline="trypsin")
    make_plots(results, args.outdir)

    if args.observed:
        print("[4/4] Integrating observed peptide sets", file=sys.stderr)
        obs = read_observed(args.observed, il_equivalent)
        write_observed_recovery(obs, results["trypsin"], args.outdir)
    else:
        print("[4/4] No --observed table supplied; theoretical analysis complete.", file=sys.stderr)

    print(f"Done. Outputs: {args.outdir.resolve()}", file=sys.stderr)


if __name__ == "__main__":
    main()
