#!/usr/bin/env python3
"""
Step01B_extract_and_compare_observed_v1.0.0.py

Purpose
-------
1) Extract the frozen MBR-OFF peptide sequence sets from:
   - AlphaPept (peptide_fdr in *.ms_data.hdf or exported peptide_fdr table)
   - FragPipe/MSFragger (psm.tsv)
   - MetaMorpheus (AllPSMs*.tsv / *.psmtsv)
   - MaxQuant (msms.txt)
2) Apply the manuscript's workflow-specific target/FDR rules where they are
   explicitly available from the output.
3) Normalize peptide sequences into the same I/L-equivalent key used by the
   theoretical Step01 analysis (I and L -> J).
4) Remap observed peptides directly against the exact common UniProt FASTA, independent of enzyme assumptions, and reproduce the manuscript accession-set classes.
5) Separately compare observed peptides with peptide_catalog_trypsin.tsv.gz from Step01 v1.0.3.
6) Write manuscript-reproduction and strict-theoretical recovery summaries.

This script does NOT re-run any search engine.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

SCRIPT_VERSION = "1.1.1"

EXPECTED = {
    # Expected common-reference peptide counts from the current manuscript.
    # The category expectations are the manuscript's primary
    # single_isoform_unique + within_family_subset_discriminative breakdown.
    "AP": {
        "common_reference": 23767,
        "single_isoform_unique": 69,
        "within_family_subset_discriminative": 111,
        "primary_total": 180,
    },
    "FP": {
        "common_reference": 22674,
        "single_isoform_unique": 13,
        "within_family_subset_discriminative": 88,
        "primary_total": 101,
    },
    "MM": {
        "common_reference": 24543,
        "single_isoform_unique": 42,
        "within_family_subset_discriminative": 107,
        "primary_total": 149,
    },
    "MQ": {
        "common_reference": 18962,
        "single_isoform_unique": 14,
        "within_family_subset_discriminative": 56,
        "primary_total": 70,
    },
}

TRUE_VALUES = {"1", "true", "t", "yes", "y", "+", "decoy", "contaminant"}
FALSE_VALUES = {"0", "false", "f", "no", "n", "-", "", "target"}


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def norm_col(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(x).strip().lower())


def normalize_il(seq: str) -> str:
    return seq.replace("I", "J").replace("L", "J")


def clean_sequence(value) -> Optional[str]:
    if value is None:
        return None
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    seq = str(value).strip().upper()
    if not seq or seq in {"NAN", "NONE"}:
        return None

    # Handle common flanking-AA form K.PEPTIDE.R.
    m = re.fullmatch(r"[A-Z\-]\.([A-Z]+)\.[A-Z\-]", seq)
    if m:
        seq = m.group(1)

    # We intentionally do not "strip" arbitrary modification text, because that
    # can silently turn e.g. M(OXIDATION) into a false sequence. The manuscript
    # uses stripped/base sequence fields, so only plain AA strings are accepted.
    if not re.fullmatch(r"[A-Z]+", seq):
        return None
    return seq


def truthy(v) -> bool:
    if v is None:
        return False
    if isinstance(v, bytes):
        v = v.decode("utf-8", errors="replace")
    s = str(v).strip().lower()
    return s in TRUE_VALUES


def is_target_value(v) -> Optional[bool]:
    if v is None:
        return None
    if isinstance(v, bytes):
        v = v.decode("utf-8", errors="replace")
    s = str(v).strip().lower()
    if s in {"t", "target", "true", "1"}:
        return True
    if s in {"d", "decoy", "c", "contaminant", "false", "0"}:
        return False
    return None


def float_or_none(v) -> Optional[float]:
    try:
        if isinstance(v, bytes):
            v = v.decode("utf-8", errors="replace")
        return float(v)
    except Exception:
        return None


def detect_delimiter(path: Path) -> str:
    if path.suffix.lower() == ".csv":
        return ","
    return "\t"


def read_table_rows(path: Path) -> Iterator[dict]:
    delim = detect_delimiter(path)
    with path.open("rt", encoding="utf-8-sig", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        if not reader.fieldnames:
            raise ValueError(f"No header found in {path}")
        for row in reader:
            yield row


def find_col(fieldnames: Sequence[str], aliases: Sequence[str]) -> Optional[str]:
    lookup = {norm_col(x): x for x in fieldnames}
    for a in aliases:
        k = norm_col(a)
        if k in lookup:
            return lookup[k]
    return None


@dataclass
class Extraction:
    workflow: str
    source_files: List[str]
    sequences: Set[str]
    n_rows_seen: int
    n_rows_kept: int
    n_invalid_sequence: int
    notes: List[str]


# ---------------------------------------------------------------------
# FragPipe
# ---------------------------------------------------------------------

def parse_fragpipe(paths: Sequence[Path]) -> Extraction:
    seqs: Set[str] = set()
    seen = kept = invalid = 0
    notes: List[str] = []

    for path in paths:
        rows = read_table_rows(path)
        first = next(rows, None)
        if first is None:
            continue
        fields = list(first)
        seq_col = find_col(fields, ["Peptide", "Peptide Sequence", "Sequence"])
        if not seq_col:
            raise ValueError(f"FP: peptide sequence column not found in {path}; columns={fields}")

        decoy_col = find_col(fields, ["Is Decoy", "Decoy"])
        contam_col = find_col(fields, ["Is Contaminant", "Contaminant"])

        def handle(row):
            nonlocal seen, kept, invalid
            seen += 1
            if decoy_col and truthy(row.get(decoy_col)):
                return
            if contam_col and truthy(row.get(contam_col)):
                return
            seq = clean_sequence(row.get(seq_col))
            if not seq:
                invalid += 1
                return
            seqs.add(seq)
            kept += 1

        handle(first)
        for row in rows:
            handle(row)

    notes.append("psm.tsv treated as the already workflow-FDR-filtered FragPipe/Philosopher PSM report.")
    return Extraction("FP", [str(p) for p in paths], seqs, seen, kept, invalid, notes)


# ---------------------------------------------------------------------
# MetaMorpheus
# ---------------------------------------------------------------------

def parse_metamorpheus(paths: Sequence[Path]) -> Extraction:
    seqs: Set[str] = set()
    seen = kept = invalid = 0
    notes: List[str] = []

    for path in paths:
        rows = read_table_rows(path)
        first = next(rows, None)
        if first is None:
            continue
        fields = list(first)
        seq_col = find_col(fields, ["Base Sequence", "BaseSequence", "Peptide"])
        q_col = find_col(fields, ["QValue", "Q Value", "q-value", "qvalue"])
        tdc_col = find_col(fields, ["Decoy/Contaminant/Target", "Decoy Contaminant Target", "DCT"])

        if not seq_col:
            raise ValueError(f"MM: Base Sequence column not found in {path}; columns={fields}")
        if not q_col:
            raise ValueError(f"MM: QValue column not found in {path}; columns={fields}")
        if not tdc_col:
            raise ValueError(
                f"MM: Decoy/Contaminant/Target column not found in {path}; columns={fields}"
            )

        def handle(row):
            nonlocal seen, kept, invalid
            seen += 1
            q = float_or_none(row.get(q_col))
            if q is None or q > 0.01:
                return
            target = is_target_value(row.get(tdc_col))
            if target is not True:
                return
            raw = str(row.get(seq_col, "") or "")
            if "|" in raw:
                # Matches manuscript: pipe-delimited ambiguous Base Sequence excluded.
                return
            seq = clean_sequence(raw)
            if not seq:
                invalid += 1
                return
            seqs.add(seq)
            kept += 1

        handle(first)
        for row in rows:
            handle(row)

    notes.append("Applied QValue <= 0.01 and exact target classification; pipe-delimited Base Sequence excluded.")
    return Extraction("MM", [str(p) for p in paths], seqs, seen, kept, invalid, notes)


# ---------------------------------------------------------------------
# MaxQuant
# ---------------------------------------------------------------------

def parse_maxquant(paths: Sequence[Path]) -> Extraction:
    seqs: Set[str] = set()
    seen = kept = invalid = 0
    notes: List[str] = []

    for path in paths:
        rows = read_table_rows(path)
        first = next(rows, None)
        if first is None:
            continue
        fields = list(first)
        seq_col = find_col(fields, ["Sequence"])
        reverse_col = find_col(fields, ["Reverse"])
        contam_col = find_col(fields, ["Potential contaminant", "Potential Contaminant"])

        if not seq_col:
            raise ValueError(f"MQ: Sequence column not found in {path}; columns={fields}")

        def handle(row):
            nonlocal seen, kept, invalid
            seen += 1
            if reverse_col and truthy(row.get(reverse_col)):
                return
            if contam_col and truthy(row.get(contam_col)):
                return
            seq = clean_sequence(row.get(seq_col))
            if not seq:
                invalid += 1
                return
            seqs.add(seq)
            kept += 1

        handle(first)
        for row in rows:
            handle(row)

    notes.append("Removed rows marked Reverse or Potential contaminant from msms.txt.")
    return Extraction("MQ", [str(p) for p in paths], seqs, seen, kept, invalid, notes)


# ---------------------------------------------------------------------
# AlphaPept
# ---------------------------------------------------------------------

def _decode_array_value(v):
    if isinstance(v, bytes):
        return v.decode("utf-8", errors="replace")
    # numpy scalar
    try:
        return v.item()
    except Exception:
        return v


def _hdf_find_peptide_fdr_blocks(path: Path) -> List[Tuple[str, List[dict]]]:
    """
    Generic reader for AlphaPept HDF5 layouts.

    Supported:
    - compound/structured dataset called peptide_fdr
    - group called peptide_fdr containing same-length 1D datasets
    - pandas HDF table if pandas/PyTables are available

    Returns list of (HDF key/path, row dictionaries).
    """
    blocks: List[Tuple[str, List[dict]]] = []
    h5_errors = []

    try:
        import h5py  # type: ignore
        with h5py.File(path, "r") as h5:
            candidates = []

            def visitor(name, obj):
                if "peptide_fdr" in name.lower():
                    candidates.append((name, obj))
            h5.visititems(visitor)

            # Prefer exact/narrow peptide_fdr objects.
            candidates.sort(key=lambda x: (0 if x[0].lower().endswith("peptide_fdr") else 1, len(x[0])))

            consumed_groups = set()
            for name, obj in candidates:
                if isinstance(obj, h5py.Group):
                    if name in consumed_groups:
                        continue
                    cols = {}
                    lengths = []
                    for child_name, child in obj.items():
                        if isinstance(child, h5py.Dataset) and child.ndim == 1:
                            try:
                                arr = child[()]
                            except Exception:
                                continue
                            cols[child_name] = arr
                            lengths.append(len(arr))
                    if cols and lengths and len(set(lengths)) == 1:
                        n = lengths[0]
                        rows = []
                        keys = list(cols)
                        for i in range(n):
                            rows.append({k: _decode_array_value(cols[k][i]) for k in keys})
                        blocks.append((name, rows))
                        consumed_groups.add(name)

                elif isinstance(obj, h5py.Dataset):
                    dt_names = obj.dtype.names
                    if dt_names:
                        arr = obj[()]
                        rows = []
                        for rec in arr:
                            rows.append({k: _decode_array_value(rec[k]) for k in dt_names})
                        blocks.append((name, rows))
    except Exception as exc:
        h5_errors.append(f"h5py: {exc}")

    if blocks:
        return blocks

    # Fallback: pandas HDFStore / PyTables.
    try:
        import pandas as pd  # type: ignore
        with pd.HDFStore(str(path), mode="r") as store:
            keys = [k for k in store.keys() if "peptide_fdr" in k.lower()]
            for key in keys:
                df = store.get(key)
                blocks.append((key, df.to_dict(orient="records")))
    except Exception as exc:
        h5_errors.append(f"pandas/PyTables: {exc}")

    if not blocks:
        raise RuntimeError(
            f"Could not read a peptide_fdr block from {path}. "
            f"Reader diagnostics: {' | '.join(h5_errors)}"
        )
    return blocks


def _parse_ap_rows(rows: List[dict], source_label: str) -> Tuple[Set[str], int, int, int, List[str]]:
    if not rows:
        return set(), 0, 0, 0, [f"{source_label}: empty peptide_fdr block"]

    fields = list(rows[0])
    # AlphaPept provides both `sequence` (matched/modified representation)
    # and `sequence_naked` (unmodified amino-acid sequence). For sequence-space
    # remapping we MUST prefer sequence_naked.
    seq_col = find_col(fields, [
        "sequence_naked", "naked_sequence",
        "peptide_sequence", "peptide sequence", "peptide",
        "sequence"
    ])
    q_col = find_col(fields, [
        "q_value", "qvalue", "q value", "q", "fdr", "peptide_q_value"
    ])
    decoy_col = find_col(fields, [
        "decoy", "is_decoy", "is decoy", "decoy_flag"
    ])
    target_col = find_col(fields, [
        "target", "is_target", "is target", "target_decoy"
    ])

    if not seq_col:
        raise ValueError(
            f"AP: sequence column not detected in {source_label}. Columns={fields}"
        )
    if not q_col:
        raise ValueError(
            f"AP: q-value column not detected in {source_label}. Columns={fields}. "
            "The manuscript explicitly restricts AlphaPept peptide_fdr to q <= 0.01, "
            "so the script refuses to guess."
        )

    seqs = set()
    seen = kept = invalid = 0
    notes = [
        f"{source_label}: sequence={seq_col}; q={q_col}; "
        f"decoy={decoy_col or 'not detected'}; target={target_col or 'not detected'}"
    ]

    for row in rows:
        seen += 1
        q = float_or_none(row.get(q_col))
        if q is None or q > 0.01:
            continue

        if decoy_col and truthy(row.get(decoy_col)):
            continue

        if target_col:
            tv = row.get(target_col)
            parsed = is_target_value(tv)
            if parsed is False:
                continue
            # If value is an actual boolean/numeric target flag, truthy True is allowed.
            if parsed is None and not truthy(tv):
                continue

        seq = clean_sequence(row.get(seq_col))
        if not seq:
            invalid += 1
            continue
        seqs.add(seq)
        kept += 1

    return seqs, seen, kept, invalid, notes


def parse_alphapept(paths: Sequence[Path]) -> Extraction:
    seqs: Set[str] = set()
    seen = kept = invalid = 0
    notes: List[str] = []

    for path in paths:
        if path.suffix.lower() in {".tsv", ".txt", ".csv"}:
            rows_iter = read_table_rows(path)
            rows = list(rows_iter)
            s, a, b, c, n = _parse_ap_rows(rows, str(path))
            seqs |= s
            seen += a; kept += b; invalid += c; notes.extend(n)
        else:
            blocks = _hdf_find_peptide_fdr_blocks(path)
            # One AP ms_data.hdf should normally provide one relevant peptide_fdr block.
            # If several candidate blocks exist, unioning them is conservative at the
            # peptide-sequence level.
            for key, rows in blocks:
                s, a, b, c, n = _parse_ap_rows(rows, f"{path}::{key}")
                seqs |= s
                seen += a; kept += b; invalid += c; notes.extend(n)

    return Extraction("AP", [str(p) for p in paths], seqs, seen, kept, invalid, notes)


# ---------------------------------------------------------------------
# Auto-discovery
# ---------------------------------------------------------------------

def path_text(p: Path) -> str:
    return str(p).lower().replace("-", "_")


def score_candidate(path: Path, wf: str) -> int:
    s = path_text(path)
    score = 0
    if wf == "FP":
        if "fp_mbr_off" in s: score += 100
        if "fp_mbr_off_lfq" in s: score -= 40
        if "mbr_on" in s: score -= 100
        if path.name.lower() == "psm.tsv": score += 30
    elif wf == "MM":
        if "metamorpheus" in s or re.search(r"[\\/](mm)([\\/ _]|$)", s): score += 40
        if "mbr_off" in s: score += 20
        if "mbr_on" in s: score -= 30
        if "allpsms" in path.name.lower(): score += 80
    elif wf == "MQ":
        if "maxquant" in s or re.search(r"[\\/](mq)([\\/ _]|$)", s): score += 40
        if path.name.lower() == "msms.txt": score += 80
        if "mbr_on" in s: score -= 20
    return score


def discover_many(root: Path, wf: str) -> List[Path]:
    """
    Return all top-scoring files for a workflow.

    Proteomics programs commonly write one identification table per LC-MS/MS run.
    A tie among top-scoring files is therefore expected, not ambiguous. For this
    manuscript the primary FP identification set is explicitly the union of the
    nine MBR-OFF PSM reports, so all equally ranked FP_MBR_OFF psm.tsv files must
    be retained.
    """
    if wf == "FP":
        cands = list(root.rglob("psm.tsv"))
    elif wf == "MM":
        cands = [p for p in root.rglob("*") if p.is_file() and "allpsms" in p.name.lower()
                 and p.suffix.lower() in {".tsv", ".psmtsv", ".txt"}]
    elif wf == "MQ":
        cands = list(root.rglob("msms.txt"))
    else:
        raise ValueError(wf)

    if not cands:
        raise FileNotFoundError(f"No candidate file found for {wf} under {root}")

    ranked = sorted(((score_candidate(p, wf), p) for p in cands),
                    reverse=True, key=lambda x: x[0])
    best_score = ranked[0][0]
    best = sorted(p for sc, p in ranked if sc == best_score)

    eprint(f"[discover] {wf} candidates:")
    for sc, p in ranked[:20]:
        marker = "*" if sc == best_score else " "
        eprint(f"         {marker} score={sc:4d}  {p}")
    eprint(f"[discover] {wf}: selected {len(best)} top-scoring file(s); peptide sequences will be unioned.")

    return best


def discover_ap(root: Path) -> List[Path]:
    # Preferred: AlphaPept ms_data HDFs.
    hdfs = [p for p in root.rglob("*.hdf") if p.is_file() and p.name.lower().endswith("ms_data.hdf")]
    if hdfs:
        # Exclude directories that are clearly other workflow folders.
        filtered = []
        for p in hdfs:
            s = path_text(p)
            if any(tok in s for tok in ["fp_mbr", "fragpipe", "metamorpheus", "maxquant"]):
                continue
            filtered.append(p)
        if filtered:
            hdfs = filtered

        # Primary manuscript benchmark = AP_MBR_OFF.  Do not union MBR_ON.
        off = [p for p in hdfs if "ap_mbr_off" in path_text(p)]
        if off:
            hdfs = off

        eprint(f"[discover] AP: selected {len(hdfs)} primary *.ms_data.hdf files")
        for p in hdfs[:20]:
            eprint(f"           {p}")
        return sorted(hdfs)

    # Fallback: exported peptide_fdr tabular files.
    tabs = [p for p in root.rglob("*") if p.is_file()
            and "peptide_fdr" in p.name.lower()
            and p.suffix.lower() in {".tsv", ".txt", ".csv"}]
    if tabs:
        eprint(f"[discover] AP: found {len(tabs)} peptide_fdr table(s)")
        return sorted(tabs)

    raise FileNotFoundError(
        f"No AlphaPept *.ms_data.hdf or peptide_fdr table found under {root}. "
        "Supply --ap PATH (repeatable) or --ap-dir DIRECTORY."
    )


def resolve_ap_paths(args) -> List[Path]:
    if args.ap:
        return [Path(x) for x in args.ap]
    if args.ap_dir:
        d = Path(args.ap_dir)
        hdfs = sorted(d.rglob("*.ms_data.hdf"))
        if hdfs:
            return hdfs
        tabs = sorted(p for p in d.rglob("*") if p.is_file()
                      and "peptide_fdr" in p.name.lower()
                      and p.suffix.lower() in {".tsv", ".txt", ".csv"})
        if tabs:
            return tabs
        raise FileNotFoundError(f"No AP peptide_fdr source found under {d}")
    return discover_ap(Path(args.root))



# ---------------------------------------------------------------------
# Full common-reference FASTA remapping (manuscript reproduction layer)
# ---------------------------------------------------------------------

ISOFORM_SUFFIX_RE = re.compile(r"-(\d+)$")


@dataclass(frozen=True)
class RefProtein:
    accession: str
    base_accession: str
    gene: Optional[str]
    is_suffixed_isoform: bool


@dataclass
class RefClassification:
    category: str
    exact_entry: bool
    primary_isoform_discriminative: bool
    structural_discriminative: bool
    target_base_accession: Optional[str]
    target_gene: Optional[str]
    n_accessions: int
    n_base_accessions: int
    n_genes: int


def strip_isoform_suffix(accession: str) -> str:
    return ISOFORM_SUFFIX_RE.sub("", accession)


def parse_uniprot_header(header: str) -> Tuple[str, Optional[str]]:
    token = header.split()[0]
    if "|" in token:
        parts = token.split("|")
        accession = parts[1] if len(parts) >= 2 else token
    else:
        accession = token
    m = re.search(r"(?:^|\s)GN=([^\s]+)", header)
    gene = m.group(1) if m else None
    return accession, gene


def classify_reference_mapping(
    entry_indices: Set[int],
    proteins: Sequence[RefProtein],
    base_totals: Dict[str, int],
    gene_totals: Dict[str, int],
) -> RefClassification:
    if not entry_indices:
        return RefClassification(
            "unmapped", False, False, False, None, None, 0, 0, 0
        )

    entries = [proteins[i] for i in entry_indices]
    bases = {p.base_accession for p in entries}
    real_genes = {p.gene for p in entries if p.gene}
    has_missing_gene = any(p.gene is None for p in entries)
    n_genes_effective = len(real_genes) + (1 if has_missing_gene else 0)

    if len(entries) == 1:
        p = entries[0]
        if p.is_suffixed_isoform:
            return RefClassification(
                "single_isoform_unique", True, True, True,
                p.base_accession, p.gene, 1, 1, 1 if p.gene else 0
            )
        return RefClassification(
            "single_canonical_unique", True, False, True,
            p.base_accession, p.gene, 1, 1, 1 if p.gene else 0
        )

    if len(bases) == 1:
        base = next(iter(bases))
        gene = next(iter(real_genes)) if len(real_genes) == 1 and not has_missing_gene else None
        if len(entries) < base_totals[base]:
            return RefClassification(
                "within_family_subset_discriminative", False, True, True,
                base, gene, len(entries), 1, n_genes_effective
            )
        return RefClassification(
            "within_family_shared_all", False, False, False,
            base, gene, len(entries), 1, n_genes_effective
        )

    if len(real_genes) == 1 and not has_missing_gene:
        gene = next(iter(real_genes))
        if len(entries) < gene_totals[gene]:
            return RefClassification(
                "same_gene_subset_discriminative", False, False, True,
                None, gene, len(entries), len(bases), 1
            )
        return RefClassification(
            "same_gene_multi_entry_shared", False, False, False,
            None, gene, len(entries), len(bases), 1
        )

    if len(real_genes) >= 2:
        return RefClassification(
            "cross_gene_shared", False, False, False,
            None, None, len(entries), len(bases), n_genes_effective
        )

    return RefClassification(
        "multi_entry_gene_unresolved", False, False, False,
        None, None, len(entries), len(bases), n_genes_effective
    )


def remap_observed_to_full_fasta(
    fasta: Path,
    observed_keys: Set[str],
) -> Tuple[List[RefProtein], Dict[str, Set[int]], Dict[str, RefClassification]]:
    """
    Map every observed I/L-equivalent peptide key to every compatible protein
    entry in the full common UniProt FASTA by raw sequence containment.

    This is the manuscript-reproduction layer and is deliberately independent
    of enzyme/digestion assumptions.
    """
    try:
        import ahocorasick  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "Full common-reference remapping requires pyahocorasick.\n"
            "Install with: python -m pip install pyahocorasick"
        ) from exc

    A = ahocorasick.Automaton()
    for key in sorted(observed_keys):
        A.add_word(key, key)
    A.make_automaton()

    proteins: List[RefProtein] = []
    mappings: Dict[str, Set[int]] = defaultdict(set)
    base_totals: Counter = Counter()
    gene_totals: Counter = Counter()

    header = None
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

        accession, gene = parse_uniprot_header(header)
        base = strip_isoform_suffix(accession)
        idx = len(proteins)
        proteins.append(
            RefProtein(
                accession=accession,
                base_accession=base,
                gene=gene,
                is_suffixed_isoform=bool(ISOFORM_SUFFIX_RE.search(accession)),
            )
        )
        base_totals[base] += 1
        if gene:
            gene_totals[gene] += 1

        seq_key = normalize_il(seq)
        # X/U/non-standard residues stay in the protein sequence. Since peptide
        # patterns contain standard AA keys only, no match can span across a
        # non-standard character.
        for _, key in A.iter(seq_key):
            mappings[key].add(idx)

        header = None
        chunks = []

    with fasta.open("rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
            elif header is not None:
                chunks.append(line.strip())
        flush()

    classes: Dict[str, RefClassification] = {}
    for key in observed_keys:
        classes[key] = classify_reference_mapping(
            mappings.get(key, set()), proteins, base_totals, gene_totals
        )

    return proteins, mappings, classes


def summarize_reference_mapped_observed(
    wf_keys: Dict[str, Set[str]],
    ref_classes: Dict[str, RefClassification],
) -> Dict[str, dict]:
    out: Dict[str, dict] = {}
    for wf, keys in wf_keys.items():
        mapped = {k for k in keys if ref_classes[k].category != "unmapped"}
        cats = Counter(ref_classes[k].category for k in mapped)
        primary = {k for k in mapped if ref_classes[k].primary_isoform_discriminative}
        exact = {k for k in mapped if ref_classes[k].exact_entry}
        structural = {k for k in mapped if ref_classes[k].structural_discriminative}

        out[wf] = {
            "workflow": wf,
            "raw_extracted_IL_equivalent_keys": len(keys),
            "common_reference_mapped_keys": len(mapped),
            "unmapped_to_common_reference": len(keys - mapped),
            "single_isoform_unique": cats.get("single_isoform_unique", 0),
            "within_family_subset_discriminative": cats.get(
                "within_family_subset_discriminative", 0
            ),
            "primary_isoform_discriminative_total": len(primary),
            "exact_entry_total": len(exact),
            "structural_discriminative_total": len(structural),
            "single_canonical_unique": cats.get("single_canonical_unique", 0),
            "within_family_shared_all": cats.get("within_family_shared_all", 0),
            "same_gene_subset_discriminative": cats.get(
                "same_gene_subset_discriminative", 0
            ),
            "same_gene_multi_entry_shared": cats.get(
                "same_gene_multi_entry_shared", 0
            ),
            "cross_gene_shared": cats.get("cross_gene_shared", 0),
            "multi_entry_gene_unresolved": cats.get(
                "multi_entry_gene_unresolved", 0
            ),
        }
    return out


def write_reference_mapping_outputs(
    outdir: Path,
    wf_keys: Dict[str, Set[str]],
    reps: Dict[str, Dict[str, str]],
    ref_mappings: Dict[str, Set[int]],
    ref_classes: Dict[str, RefClassification],
    ref_summary: Dict[str, dict],
):
    rows = [ref_summary[wf] for wf in ["AP", "FP", "MM", "MQ"] if wf in ref_summary]
    write_tsv(outdir / "observed_common_reference_summary.tsv", rows)

    validation_rows = []
    for wf in ["AP", "FP", "MM", "MQ"]:
        s = ref_summary[wf]
        exp = EXPECTED[wf]
        checks = [
            ("common_reference", s["common_reference_mapped_keys"], exp["common_reference"]),
            ("single_isoform_unique", s["single_isoform_unique"], exp["single_isoform_unique"]),
            (
                "within_family_subset_discriminative",
                s["within_family_subset_discriminative"],
                exp["within_family_subset_discriminative"],
            ),
            (
                "primary_total",
                s["primary_isoform_discriminative_total"],
                exp["primary_total"],
            ),
        ]
        for metric, obs, expected in checks:
            validation_rows.append({
                "workflow": wf,
                "metric": metric,
                "observed_by_full_FASTA_remap": obs,
                "expected_from_manuscript": expected,
                "difference": obs - expected,
                "status": "MATCH" if obs == expected else "CHECK",
            })
    write_tsv(outdir / "validation_full_FASTA_remap_vs_manuscript.tsv", validation_rows)

    unmapped_rows = []
    for wf in ["AP", "FP", "MM", "MQ"]:
        for key in sorted(wf_keys[wf]):
            if ref_classes[key].category == "unmapped":
                unmapped_rows.append({
                    "workflow": wf,
                    "peptide_key_IL_equivalent": key,
                    "representative_sequence": reps[wf][key],
                })
    write_tsv(outdir / "observed_unmapped_to_common_FASTA.tsv", unmapped_rows)

    class_rows = []
    union_keys = sorted(set().union(*wf_keys.values()))
    for key in union_keys:
        c = ref_classes[key]
        class_rows.append({
            "peptide_key_IL_equivalent": key,
            "category": c.category,
            "exact_entry": int(c.exact_entry),
            "primary_isoform_discriminative": int(c.primary_isoform_discriminative),
            "structural_discriminative": int(c.structural_discriminative),
            "n_accessions": c.n_accessions,
            "n_base_accessions": c.n_base_accessions,
            "n_genes": c.n_genes,
            "target_base_accession": c.target_base_accession or "",
            "target_gene": c.target_gene or "",
            "AP": int(key in wf_keys.get("AP", set())),
            "FP": int(key in wf_keys.get("FP", set())),
            "MM": int(key in wf_keys.get("MM", set())),
            "MQ": int(key in wf_keys.get("MQ", set())),
        })
    write_tsv(outdir / "observed_full_FASTA_classification.tsv", class_rows)

# ---------------------------------------------------------------------
# Theoretical tryptic catalog integration
# ---------------------------------------------------------------------

def load_observed_keys(extractions: Dict[str, Extraction]) -> Tuple[Dict[str, Set[str]], Dict[str, Dict[str, str]]]:
    wf_keys: Dict[str, Set[str]] = {}
    representatives: Dict[str, Dict[str, str]] = {}
    for wf, ex in extractions.items():
        reps: Dict[str, str] = {}
        for seq in ex.sequences:
            key = normalize_il(seq)
            reps.setdefault(key, seq)
        wf_keys[wf] = set(reps)
        representatives[wf] = reps
    return wf_keys, representatives


def compare_to_catalog(
    catalog: Path,
    wf_keys: Dict[str, Set[str]],
) -> Tuple[Dict[str, dict], Dict[str, dict], Counter]:
    union_obs = set().union(*wf_keys.values()) if wf_keys else set()

    # Only hold classifications for observed keys; stream the ~3.16M-row catalog.
    observed_class: Dict[str, dict] = {}
    theoretical_counts = Counter()

    with gzip.open(catalog, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"peptide_key", "category", "exact_entry",
                    "primary_isoform_discriminative", "structural_discriminative"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(f"Theoretical catalog missing required columns: {catalog}")

        for row in reader:
            key = row["peptide_key"]
            cat = row["category"]
            theoretical_counts["all"] += 1
            theoretical_counts[f"class::{cat}"] += 1
            if row["exact_entry"] == "1":
                theoretical_counts["exact_entry"] += 1
            if row["primary_isoform_discriminative"] == "1":
                theoretical_counts["primary"] += 1
            if row["structural_discriminative"] == "1":
                theoretical_counts["structural"] += 1

            if key in union_obs:
                observed_class[key] = {
                    "category": cat,
                    "exact_entry": row["exact_entry"] == "1",
                    "primary": row["primary_isoform_discriminative"] == "1",
                    "structural": row["structural_discriminative"] == "1",
                    "target_base_accession": row.get("target_base_accession", ""),
                    "target_gene": row.get("target_gene", ""),
                }

    summaries: Dict[str, dict] = {}
    for wf in sorted(wf_keys):
        obs = wf_keys[wf]
        matched = obs & set(observed_class)
        cats = Counter(observed_class[k]["category"] for k in matched)
        n_exact = sum(1 for k in matched if observed_class[k]["exact_entry"])
        n_primary = sum(1 for k in matched if observed_class[k]["primary"])
        n_struct = sum(1 for k in matched if observed_class[k]["structural"])

        summaries[wf] = {
            "workflow": wf,
            "observed_unique_il_keys": len(obs),
            "observed_in_theoretical_tryptic_space": len(matched),
            "not_in_theoretical_tryptic_space": len(obs - matched),
            "strict_tryptic_catalog_overlap_pct_of_observed": 100 * len(matched) / len(obs)
                if obs else math.nan,
            "outside_strict_tryptic_catalog_pct_of_observed": 100 * len(obs - matched) / len(obs)
                if obs else math.nan,
            "observed_exact_entry": n_exact,
            "observed_primary_isoform_discriminative": n_primary,
            "observed_structural_discriminative": n_struct,
            "observed_single_isoform_unique": cats.get("single_isoform_unique", 0),
            "observed_within_family_subset_discriminative": cats.get(
                "within_family_subset_discriminative", 0
            ),
            "theoretical_tryptic_all": theoretical_counts["all"],
            "theoretical_tryptic_exact_entry": theoretical_counts["exact_entry"],
            "theoretical_tryptic_primary": theoretical_counts["primary"],
            "theoretical_tryptic_structural": theoretical_counts["structural"],
            "recovery_all_pct": 100 * len(matched) / theoretical_counts["all"]
                if theoretical_counts["all"] else math.nan,
            "recovery_exact_entry_pct": 100 * n_exact / theoretical_counts["exact_entry"]
                if theoretical_counts["exact_entry"] else math.nan,
            "recovery_primary_pct": 100 * n_primary / theoretical_counts["primary"]
                if theoretical_counts["primary"] else math.nan,
            "recovery_structural_pct": 100 * n_struct / theoretical_counts["structural"]
                if theoretical_counts["structural"] else math.nan,
        }

    return summaries, observed_class, theoretical_counts


# ---------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------

def write_tsv(path: Path, rows: List[dict], fieldnames: Optional[List[str]] = None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    if fieldnames is None:
        fieldnames = list(rows[0])
    with path.open("wt", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def write_outputs(
    outdir: Path,
    extractions: Dict[str, Extraction],
    wf_keys: Dict[str, Set[str]],
    reps: Dict[str, Dict[str, str]],
    summaries: Dict[str, dict],
    observed_class: Dict[str, dict],
    theoretical_counts: Counter,
):
    outdir.mkdir(parents=True, exist_ok=True)

    observed_rows = []
    for wf in ["AP", "FP", "MM", "MQ"]:
        if wf not in reps:
            continue
        for key in sorted(reps[wf]):
            observed_rows.append({
                "workflow": wf,
                "peptide_sequence": reps[wf][key],
                "peptide_key_IL_equivalent": key,
            })
    write_tsv(outdir / "observed_peptides.tsv", observed_rows)

    extraction_rows = []
    for wf in ["AP", "FP", "MM", "MQ"]:
        ex = extractions.get(wf)
        if not ex:
            continue
        extraction_rows.append({
            "workflow": wf,
            "n_source_files": len(ex.source_files),
            "n_rows_seen": ex.n_rows_seen,
            "n_rows_kept_before_sequence_dedup": ex.n_rows_kept,
            "n_invalid_or_nonplain_sequence_rows": ex.n_invalid_sequence,
            "n_unique_plain_sequences": len(ex.sequences),
            "n_unique_IL_equivalent_keys": len(wf_keys[wf]),
            "source_files": " | ".join(ex.source_files),
            "notes": " | ".join(ex.notes),
        })
    write_tsv(outdir / "observed_extraction_summary.tsv", extraction_rows)

    recovery_rows = [summaries[wf] for wf in ["AP", "FP", "MM", "MQ"] if wf in summaries]
    write_tsv(outdir / "observed_recovery_by_workflow.tsv", recovery_rows)

    # Manuscript extraction self-check.
    #
    # IMPORTANT: the manuscript's "common-reference peptide" count is an
    # observed peptide-set count after sequence normalization/remapping.  It is
    # NOT the count remaining after intersection with the strict in-silico
    # digestion catalog.  v1.0.0/v1.0.1 incorrectly conflated those metrics.
    val_rows = []
    for wf in ["AP", "FP", "MM", "MQ"]:
        if wf not in summaries:
            continue
        s = summaries[wf]
        exp = EXPECTED[wf]
        extracted_il = len(wf_keys[wf])
        val_rows.append({
            "workflow": wf,
            "metric": "common_reference_IL_equivalent_extracted",
            "observed_by_Step01B": extracted_il,
            "expected_from_manuscript": exp["common_reference"],
            "difference": extracted_il - exp["common_reference"],
            "status": "MATCH" if extracted_il == exp["common_reference"] else "CHECK",
        })

        # The following are intentionally labelled as theoretical-catalog
        # overlaps. They are NEW Step01B metrics, not manuscript reproduction
        # checks, because the strict digest model may exclude valid observed
        # sequences (e.g. search-engine terminal processing semantics).
        for metric, obs in [
            ("strict_tryptic_catalog_overlap", s["observed_in_theoretical_tryptic_space"]),
            ("strict_tryptic_single_isoform_unique", s["observed_single_isoform_unique"]),
            ("strict_tryptic_within_family_subset_discriminative",
             s["observed_within_family_subset_discriminative"]),
            ("strict_tryptic_primary_total", s["observed_primary_isoform_discriminative"]),
        ]:
            val_rows.append({
                "workflow": wf,
                "metric": metric,
                "observed_by_Step01B": obs,
                "expected_from_manuscript": "",
                "difference": "",
                "status": "NEW_METRIC",
            })
    write_tsv(outdir / "validation_vs_manuscript.tsv", val_rows)

    # Cross-workflow support among primary theoretical discriminative peptides.
    support = Counter()
    for wf, keys in wf_keys.items():
        for key in keys:
            c = observed_class.get(key)
            if c and c["primary"]:
                support[key] += 1

    support_rows = []
    for key, n in sorted(support.items(), key=lambda x: (-x[1], x[0])):
        c = observed_class[key]
        workflows = [wf for wf in ["AP", "FP", "MM", "MQ"] if key in wf_keys.get(wf, set())]
        support_rows.append({
            "peptide_key": key,
            "n_supporting_workflows": n,
            "supporting_workflows": ",".join(workflows),
            "category": c["category"],
            "target_base_accession": c["target_base_accession"],
            "target_gene": c["target_gene"],
        })
    write_tsv(outdir / "observed_primary_discriminative_support.tsv", support_rows)

    # Observed keys absent from theoretical tryptic space — useful QC.
    absent_rows = []
    for wf in ["AP", "FP", "MM", "MQ"]:
        for key in sorted(wf_keys.get(wf, set()) - set(observed_class)):
            absent_rows.append({
                "workflow": wf,
                "peptide_key_IL_equivalent": key,
                "representative_sequence": reps[wf][key],
            })
    write_tsv(outdir / "observed_not_in_theoretical_tryptic_space.tsv", absent_rows)

    # Manifest.
    manifest = {
        "script_version": SCRIPT_VERSION,
        "theoretical_counts": dict(theoretical_counts),
        "expected_manuscript_counts": EXPECTED,
        "interpretation_note": (
            "Observed recovery is calculated after I/L-equivalent normalization and "
            "intersection with the v1.0.3 theoretical tryptic catalog. A mismatch with "
            "the manuscript self-check should be resolved before biological interpretation."
        ),
    }
    (outdir / "Step01B_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    # Figures.
    try:
        import matplotlib.pyplot as plt  # type: ignore

        wfs = [wf for wf in ["AP", "FP", "MM", "MQ"] if wf in summaries]
        recovery = [summaries[wf]["recovery_primary_pct"] for wf in wfs]
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        ax.bar(wfs, recovery)
        ax.set_ylabel("Recovery of theoretical primary\nisoform-discriminative tryptic space (%)")
        ax.set_xlabel("Workflow")
        ax.set_title("Empirical recovery of the theoretical isoform-discriminative peptide space")
        ax.grid(axis="y", alpha=0.25)
        fig.tight_layout()
        fig.savefig(outdir / "Fig_Step01B_primary_recovery_pct.png", dpi=300)
        fig.savefig(outdir / "Fig_Step01B_primary_recovery_pct.pdf")
        plt.close(fig)

        counts = [theoretical_counts["primary"]] + [
            summaries[wf]["observed_primary_isoform_discriminative"] for wf in wfs
        ]
        labels = ["Theoretical"] + wfs
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        ax.bar(labels, counts)
        ax.set_yscale("log")
        ax.set_ylabel("Primary isoform-discriminative peptide keys (log scale)")
        ax.set_title("Theoretical versus empirically recovered isoform-discriminative space")
        ax.grid(axis="y", alpha=0.25)
        fig.tight_layout()
        fig.savefig(outdir / "Fig_Step01B_theoretical_vs_observed.png", dpi=300)
        fig.savefig(outdir / "Fig_Step01B_theoretical_vs_observed.pdf")
        plt.close(fig)
    except Exception as exc:
        eprint(f"[warning] Figure generation skipped: {exc}")


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description="Extract AP/FP/MM/MQ frozen peptide sets and compare them with Step01 theoretical tryptic space."
    )
    p.add_argument("--root", required=True, help="Benchmark_Program root directory.")
    p.add_argument("--theoretical-dir", required=True,
                   help="Step01_results_v103 directory containing peptide_catalog_trypsin.tsv.gz.")
    p.add_argument("--fasta", required=True,
                   help="Exact common UniProt FASTA used in the benchmark manuscript.")
    p.add_argument("--outdir", default="Step01B_observed_results")
    p.add_argument("--ap", action="append",
                   help="Explicit AlphaPept HDF/table path; repeat for multiple files.")
    p.add_argument("--ap-dir", help="Explicit AlphaPept results directory containing *.ms_data.hdf.")
    p.add_argument("--fp", action="append", help="Explicit FragPipe psm.tsv path; repeat for multiple run files.")
    p.add_argument("--mm", action="append", help="Explicit MetaMorpheus AllPSMs path; repeat for multiple run files.")
    p.add_argument("--mq", action="append", help="Explicit MaxQuant msms.txt path; repeat for multiple run files.")
    p.add_argument("--extract-only", action="store_true",
                   help="Only build observed_peptides.tsv; do not compare to theoretical catalog.")
    args = p.parse_args()

    root = Path(args.root)
    theoretical_dir = Path(args.theoretical_dir)
    fasta = Path(args.fasta)
    outdir = Path(args.outdir)

    if not root.exists():
        raise SystemExit(f"Benchmark root not found: {root}")
    if not fasta.exists():
        raise SystemExit(f"Common reference FASTA not found: {fasta}")

    eprint(f"[1/5] Step01B v{SCRIPT_VERSION}: discovering workflow outputs")
    ap_paths = resolve_ap_paths(args)
    fp_paths = [Path(x) for x in args.fp] if args.fp else discover_many(root, "FP")
    mm_paths = [Path(x) for x in args.mm] if args.mm else discover_many(root, "MM")
    mq_paths = [Path(x) for x in args.mq] if args.mq else discover_many(root, "MQ")

    for pth in [*ap_paths, *fp_paths, *mm_paths, *mq_paths]:
        if not pth.exists():
            raise SystemExit(f"Input file not found: {pth}")

    eprint("[2/5] Extracting frozen workflow peptide sets")
    extractions = {
        "AP": parse_alphapept(ap_paths),
        "FP": parse_fragpipe(fp_paths),
        "MM": parse_metamorpheus(mm_paths),
        "MQ": parse_maxquant(mq_paths),
    }

    for wf in ["AP", "FP", "MM", "MQ"]:
        ex = extractions[wf]
        eprint(
            f"      {wf}: {len(ex.sequences):,} unique plain sequences "
            f"from {len(ex.source_files)} source file(s)"
        )

    wf_keys, reps = load_observed_keys(extractions)

    eprint("[3/6] Full common-reference FASTA remapping (enzyme-independent)")
    union_observed = set().union(*wf_keys.values())
    ref_proteins, ref_mappings, ref_classes = remap_observed_to_full_fasta(
        fasta, union_observed
    )
    ref_summary = summarize_reference_mapped_observed(wf_keys, ref_classes)

    if args.extract_only:
        outdir.mkdir(parents=True, exist_ok=True)
        observed_rows = []
        for wf in ["AP", "FP", "MM", "MQ"]:
            for key in sorted(reps[wf]):
                observed_rows.append({
                    "workflow": wf,
                    "peptide_sequence": reps[wf][key],
                    "peptide_key_IL_equivalent": key,
                })
        write_tsv(outdir / "observed_peptides.tsv", observed_rows)
        write_reference_mapping_outputs(
            outdir, wf_keys, reps, ref_mappings, ref_classes, ref_summary
        )
        eprint(f"[done] Extract-only output: {outdir.resolve()}")
        return

    catalog = theoretical_dir / "peptide_catalog_trypsin.tsv.gz"
    if not catalog.exists():
        raise SystemExit(f"Theoretical trypsin catalog not found: {catalog}")

    eprint("[4/6] Streaming strict theoretical tryptic catalog and matching observed peptides")
    summaries, observed_class, theoretical_counts = compare_to_catalog(catalog, wf_keys)

    eprint("[5/6] Writing reference-remap validation, strict-catalog tables and figures")
    write_outputs(
        outdir, extractions, wf_keys, reps, summaries,
        observed_class, theoretical_counts
    )
    write_reference_mapping_outputs(
        outdir, wf_keys, reps, ref_mappings, ref_classes, ref_summary
    )

    eprint("[6/6] Full-FASTA manuscript reproduction self-check")
    all_ok = True
    for wf in ["AP", "FP", "MM", "MQ"]:
        s = ref_summary[wf]
        exp = EXPECTED[wf]
        vals = [
            ("common", s["common_reference_mapped_keys"], exp["common_reference"]),
            ("single", s["single_isoform_unique"], exp["single_isoform_unique"]),
            ("subset", s["within_family_subset_discriminative"],
             exp["within_family_subset_discriminative"]),
            ("primary", s["primary_isoform_discriminative_total"], exp["primary_total"]),
        ]
        ok = all(obs == expected for _, obs, expected in vals)
        all_ok &= ok
        eprint(
            f"      {wf}: {'MATCH' if ok else 'CHECK'}  "
            + ", ".join(f"{name}={obs}/{expected}" for name, obs, expected in vals)
            + f"; raw_extracted={s['raw_extracted_IL_equivalent_keys']}; "
              f"unmapped_to_FASTA={s['unmapped_to_common_reference']}"
        )

    if all_ok:
        eprint(
            "\n[success] Full FASTA sequence remapping reproduces the manuscript "
            "common-reference and primary isoform-discriminative counts for all workflows."
        )
    else:
        eprint(
            "\n[important] At least one workflow still differs after full FASTA remapping. "
            "Inspect validation_full_FASTA_remap_vs_manuscript.tsv and "
            "observed_unmapped_to_common_FASTA.tsv. Do not modify the theoretical "
            "denominator to force agreement."
        )

    eprint(
        "[note] The strict tryptic catalog remains a separate NEW analysis. "
        "Its overlap percentages should be interpreted as digestion-model/theoretical-space "
        "recovery, not as manuscript reproduction."
    )

    eprint(f"Done. Outputs: {outdir.resolve()}")


if __name__ == "__main__":
    main()
