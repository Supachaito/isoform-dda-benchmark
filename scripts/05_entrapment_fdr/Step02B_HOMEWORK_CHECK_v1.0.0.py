#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations
import argparse, csv, json, math, re, sqlite3
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


VERSION = "1.0.0"
WORKFLOWS = ["AP", "FP", "MM", "MQ"]

DEFAULT_ROOT = _public_project_root() / "ENTRAPMENT_FDR"

EXPECTED_STEP02A = {
    "target_proteins": 169637,
    "entrapment_proteins": 169637,
    "target_only_peptide_keys": 3158816,
    "entrapment_only_peptide_keys": 3169203,
    "shared_target_entrapment_peptide_keys": 5507,
}

def norm_col(x):
    return re.sub(r"[^a-z0-9]+", "", str(x or "").strip().lower())

def find_col(fields, aliases):
    by_norm = {norm_col(x): x for x in fields}
    for alias in aliases:
        n = norm_col(alias)
        if n in by_norm:
            return by_norm[n]
    for alias in aliases:
        n = norm_col(alias)
        for col in fields:
            if n and n in norm_col(col):
                return col
    return None

def normalize_il(seq):
    return seq.replace("I", "J").replace("L", "J")

def clean_sequence(value):
    if value is None:
        return None
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    s = str(value).strip().upper()
    m = re.fullmatch(r"[A-Z\-]\.([A-Z]+)\.[A-Z\-]", s)
    if m:
        s = m.group(1)
    return s if re.fullmatch(r"[A-Z]+", s or "") else None

def truthy(value):
    return str(value or "").strip().lower() in {
        "1", "true", "t", "yes", "y", "+",
        "decoy", "contaminant", "reverse"
    }

def float_or_none(value):
    try:
        if isinstance(value, bytes):
            value = value.decode("utf-8", errors="replace")
        return float(value)
    except Exception:
        return None

def target_value(value):
    s = str(value or "").strip().lower()
    if s in {"t", "target", "true", "1"}:
        return True
    if s in {"d", "decoy", "c", "contaminant", "false", "0"}:
        return False
    return None

def read_delimited(path):
    delim = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as fh:
        yield from csv.DictReader(fh, delimiter=delim)

def locate_step02a(root):
    direct = root / "Step02A_entrapment_db_v104"
    if direct.exists():
        return direct
    hits = [p for p in root.rglob("*") if p.is_dir() and "step02a_entrapment_db" in p.name.lower()]
    if not hits:
        raise FileNotFoundError(f"Could not find Step02A_entrapment_db_v104 under:\n{root}")
    return sorted(hits, key=lambda p: len(str(p)))[0]

def audit_step02a(step02a):
    required = [
        step02a / "Step02_target_plus_shuffled_entrapment_r1.fasta",
        step02a / "Step02_entrapment_peptide_space.sqlite",
        step02a / "Step02A_manifest.json",
        step02a / "Step02A_space_summary.tsv",
        step02a / "Step02A_reproduction_QC.tsv",
    ]
    missing = [p for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing Step02A files:\n" + "\n".join(map(str, missing)))

    manifest = json.loads((step02a / "Step02A_manifest.json").read_text(encoding="utf-8"))
    pc = manifest.get("protein_counts", {})
    pepc = manifest.get("peptide_space_counts", {})
    observed = {
        "target_proteins": int(pc.get("target", -1)),
        "entrapment_proteins": int(pc.get("entrapment", -1)),
        "target_only_peptide_keys": int(pepc.get("target_only", -1)),
        "entrapment_only_peptide_keys": int(pepc.get("entrapment_only", -1)),
        "shared_target_entrapment_peptide_keys": int(pepc.get("shared", -1)),
    }
    qc = []
    for metric, expected in EXPECTED_STEP02A.items():
        qc.append({
            "metric": metric,
            "expected_v103": expected,
            "observed": observed[metric],
            "pass": observed[metric] == expected,
        })
    r = float(manifest.get("effective_r_entrapment_over_target", float("nan")))
    return qc, r

def locate_workflow_dir(root, workflow):
    exact = [p for p in root.rglob("*") if p.is_dir() and p.name.lower() == f"{workflow}_entrap".lower()]
    if exact:
        return sorted(exact, key=lambda p: len(str(p)))[0]
    token = f"{workflow}_entrap".lower()
    hits = [p for p in root.rglob("*") if p.is_dir() and token in p.name.lower()]
    if not hits:
        raise FileNotFoundError(f"No {workflow}_ENTRAP directory found under:\n{root}")
    return sorted(hits, key=lambda p: len(str(p)))[0]

def _decode(v):
    if isinstance(v, bytes):
        return v.decode("utf-8", errors="replace")
    try:
        return v.item()
    except Exception:
        return v

def alpha_pept_blocks(path):
    import h5py
    blocks = []
    with h5py.File(path, "r") as h5:
        candidates = []
        def visitor(name, obj):
            if "peptide_fdr" in name.lower():
                candidates.append((name, obj))
        h5.visititems(visitor)
        for name, obj in candidates:
            if isinstance(obj, h5py.Group):
                cols, lengths = {}, []
                for col_name, child in obj.items():
                    if isinstance(child, h5py.Dataset) and child.ndim == 1:
                        arr = child[()]
                        cols[col_name] = arr
                        lengths.append(len(arr))
                if cols and len(set(lengths)) == 1:
                    n = lengths[0]
                    rows = [{k: _decode(v[i]) for k, v in cols.items()} for i in range(n)]
                    blocks.append((name, rows))
            elif isinstance(obj, h5py.Dataset) and obj.dtype.names:
                arr = obj[()]
                rows = [{k: _decode(rec[k]) for k in obj.dtype.names} for rec in arr]
                blocks.append((name, rows))
    return blocks

def parse_ap(workflow_dir, alpha=0.01):
    try:
        import h5py  # noqa
    except ImportError:
        raise RuntimeError("AP audit requires h5py. Install with: python -m pip install h5py")

    hdfs = sorted(workflow_dir.rglob("*.hdf"))
    if not hdfs:
        raise FileNotFoundError(f"[AP] No HDF files under {workflow_dir}")

    accepted, used = set(), []
    blocks_seen = 0
    for path in hdfs:
        try:
            blocks = alpha_pept_blocks(path)
        except Exception:
            continue
        if not blocks:
            continue
        used.append(path)
        for _, rows in blocks:
            blocks_seen += 1
            if not rows:
                continue
            fields = list(rows[0].keys())
            seq_col = find_col(fields, ["sequence_naked","naked_sequence","peptide_sequence","peptide","sequence"])
            q_col = find_col(fields, ["q_value","qvalue","q value","fdr","peptide_q_value"])
            decoy_col = find_col(fields, ["decoy","is_decoy","decoy_flag"])
            target_col = find_col(fields, ["target","is_target","target_decoy"])
            if seq_col is None:
                continue
            for row in rows:
                if q_col:
                    q = float_or_none(row.get(q_col))
                    if q is not None and q > alpha:
                        continue
                if decoy_col and truthy(row.get(decoy_col)):
                    continue
                if target_col and target_value(row.get(target_col)) is False:
                    continue
                seq = clean_sequence(row.get(seq_col))
                if seq:
                    accepted.add(normalize_il(seq))
    if not accepted:
        raise RuntimeError(f"[AP] No accepted peptides extracted; peptide_fdr blocks seen={blocks_seen}")
    return accepted, used, {"hdf_files_scanned": len(hdfs), "peptide_fdr_blocks_seen": blocks_seen}

def parse_fp(workflow_dir, alpha=0.01):
    paths = sorted(workflow_dir.rglob("psm.tsv"))
    if not paths:
        raise FileNotFoundError(f"[FP] No psm.tsv under {workflow_dir}")
    accepted, used = set(), []
    for path in paths:
        rows = read_delimited(path)
        first = next(rows, None)
        if first is None:
            continue
        fields = list(first.keys())
        seq_col = find_col(fields, ["Peptide","Peptide Sequence","Sequence"])
        decoy_col = find_col(fields, ["Is Decoy","Decoy"])
        contaminant_col = find_col(fields, ["Is Contaminant","Contaminant"])
        q_col = find_col(fields, ["Peptide QValue","QValue","q-value","q_value"])
        if seq_col is None:
            continue
        used.append(path)
        for row in [first, *rows]:
            if decoy_col and truthy(row.get(decoy_col)):
                continue
            if contaminant_col and truthy(row.get(contaminant_col)):
                continue
            if q_col:
                q = float_or_none(row.get(q_col))
                if q is not None and q > alpha:
                    continue
            seq = clean_sequence(row.get(seq_col))
            if seq:
                accepted.add(normalize_il(seq))
    if not accepted:
        raise RuntimeError("[FP] No accepted peptides extracted.")
    return accepted, used, {"psm_files": len(used)}

def parse_mm(workflow_dir, alpha=0.01):
    all_files = [p for p in workflow_dir.rglob("*") if p.is_file()]
    pep = sorted(p for p in all_files if "allpeptides" in p.name.lower() and p.suffix.lower() in {".tsv",".txt",".psmtsv"})
    psm = sorted(p for p in all_files if "allpsms" in p.name.lower() and p.suffix.lower() in {".tsv",".txt",".psmtsv"})
    paths = pep or psm
    if not paths:
        raise FileNotFoundError(f"[MM] No AllPeptides/AllPSMs under {workflow_dir}")
    accepted, used = set(), []
    for path in paths:
        rows = read_delimited(path)
        first = next(rows, None)
        if first is None:
            continue
        fields = list(first.keys())
        seq_col = find_col(fields, ["Base Sequence","BaseSequence","Peptide"])
        q_col = find_col(fields, ["QValue","Q Value","q-value","q_value","qvalue"])
        target_col = find_col(fields, ["Decoy/Contaminant/Target","Decoy Contaminant Target","DCT"])
        if seq_col is None:
            continue
        used.append(path)
        for row in [first, *rows]:
            if q_col:
                q = float_or_none(row.get(q_col))
                if q is not None and q > alpha:
                    continue
            if target_col and target_value(row.get(target_col)) is False:
                continue
            raw = str(row.get(seq_col, "") or "")
            if "|" in raw:
                continue
            seq = clean_sequence(raw)
            if seq:
                accepted.add(normalize_il(seq))
    if not accepted:
        raise RuntimeError("[MM] No accepted peptides extracted.")
    return accepted, used, {"source_type": "AllPeptides" if pep else "AllPSMs", "files": len(used)}

def parse_mq(workflow_dir, alpha=0.01):
    paths = sorted(workflow_dir.rglob("peptides.txt"))
    if not paths:
        raise FileNotFoundError(f"[MQ] No peptides.txt under {workflow_dir}")
    paths.sort(key=lambda p: (0 if p.parent.name.lower() == "txt" else 1, len(str(p))))
    accepted, used = set(), []
    for path in paths:
        rows = read_delimited(path)
        first = next(rows, None)
        if first is None:
            continue
        fields = list(first.keys())
        seq_col = find_col(fields, ["Sequence"])
        rev_col = find_col(fields, ["Reverse"])
        cont_col = find_col(fields, ["Potential contaminant","Potential Contaminant"])
        if seq_col is None:
            continue
        used.append(path)
        for row in [first, *rows]:
            if rev_col and truthy(row.get(rev_col)):
                continue
            if cont_col and truthy(row.get(cont_col)):
                continue
            seq = clean_sequence(row.get(seq_col))
            if seq:
                accepted.add(normalize_il(seq))
        if accepted:
            break
    if not accepted:
        raise RuntimeError("[MQ] No accepted peptides extracted.")
    return accepted, used, {"peptides_files_considered": len(paths), "peptides_files_used": len(used)}

def classify_keys(db_path, keys):
    con = sqlite3.connect(str(db_path))
    cur = con.cursor()
    found = {}
    kl = list(keys)
    for i in range(0, len(kl), 800):
        ch = kl[i:i+800]
        qs = ",".join("?" for _ in ch)
        cur.execute(
            f"SELECT peptide_key,is_target,is_entrapment FROM peptide_space WHERE peptide_key IN ({qs})",
            ch
        )
        for k, t, e in cur.fetchall():
            found[k] = (int(t), int(e))
    con.close()

    out = {"target_only": set(), "entrapment_only": set(), "shared": set(), "unmatched": set()}
    for k in keys:
        if k not in found:
            out["unmatched"].add(k)
        else:
            t, e = found[k]
            if t and e:
                out["shared"].add(k)
            elif t:
                out["target_only"].add(k)
            elif e:
                out["entrapment_only"].add(k)
            else:
                out["unmatched"].add(k)
    return out

def interpretation(lower, combined, nominal):
    if combined <= nominal:
        return "EVIDENCE_CONSISTENT_WITH_CONTROL"
    if lower > nominal:
        return "EVIDENCE_SUGGESTING_FAILURE"
    return "INCONCLUSIVE_BOUNDS_STRADDLE_NOMINAL"

def write_tsv(path, rows, fieldnames):
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(DEFAULT_ROOT))
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    root = Path(args.root)
    if not root.exists():
        raise SystemExit(f"Root does not exist:\n{root}")

    outdir = Path(args.outdir) if args.outdir else root / "Step02B_HOMEWORK_CHECK_v100"
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 100)
    print(f"STEP02B HOMEWORK CHECK v{VERSION}")
    print("=" * 100)
    print("Root:", root)
    print("Nominal peptide FDR:", f"{100*args.alpha:.2f}%")

    step02a = locate_step02a(root)
    qc_rows, effective_r = audit_step02a(step02a)
    db_path = step02a / "Step02_entrapment_peptide_space.sqlite"

    print("\n[1/4] Step02A frozen QC")
    print("Folder:", step02a)
    print("effective r =", f"{effective_r:.12f}")
    for row in qc_rows:
        print(("PASS" if row["pass"] else "FAIL"), row["metric"], row["observed"])
    if not all(r["pass"] for r in qc_rows):
        raise SystemExit("Step02A QC mismatch. Stop.")

    print("\n[2/4] Workflow folders")
    workflow_dirs = {}
    for wf in WORKFLOWS:
        workflow_dirs[wf] = locate_workflow_dir(root, wf)
        print(f"{wf}: {workflow_dirs[wf]}")

    parsers = {"AP": parse_ap, "FP": parse_fp, "MM": parse_mm, "MQ": parse_mq}

    summary_rows, discovery_rows, source_rows = [], [], []

    print("\n[3/4] Peptide-level fixed-cutoff evaluation")
    for wf in WORKFLOWS:
        print(f"\n[{wf}]")
        keys, used_files, meta = parsers[wf](workflow_dirs[wf], alpha=args.alpha)
        classes = classify_keys(db_path, keys)

        T = len(classes["target_only"])
        E = len(classes["entrapment_only"])
        S = len(classes["shared"])
        U = len(classes["unmatched"])
        denom = T + E

        lower = E / denom if denom else math.nan
        combined = E * (1.0 + 1.0 / effective_r) / denom if denom else math.nan
        status = interpretation(lower, combined, args.alpha)

        row = {
            "workflow": wf,
            "nominal_fdr": args.alpha,
            "effective_r": effective_r,
            "accepted_unique_IL_peptide_keys": len(keys),
            "target_only": T,
            "entrapment_only": E,
            "shared_excluded": S,
            "unmatched_qc": U,
            "evaluable_T_plus_E": denom,
            "lower_bound_fdp": lower,
            "combined_fdp": combined,
            "lower_bound_pct": 100 * lower,
            "combined_fdp_pct": 100 * combined,
            "interpretation": status,
        }
        summary_rows.append(row)

        print(f"accepted={len(keys):,} | T={T:,} | E={E:,} | shared={S:,} | unmatched={U:,}")
        print(f"lower FDP={100*lower:.4f}% | combined FDP={100*combined:.4f}% | {status}")

        for cls, seqs in classes.items():
            for seq in sorted(seqs):
                discovery_rows.append({
                    "workflow": wf,
                    "class": cls,
                    "peptide_key_IL_equivalent": seq
                })

        for p in used_files:
            source_rows.append({
                "workflow": wf,
                "source_file": str(p),
                "size_mb": round(p.stat().st_size / 1024 / 1024, 4)
            })
        source_rows.append({
            "workflow": wf,
            "source_file": "PARSER_METADATA=" + json.dumps(meta, sort_keys=True),
            "size_mb": ""
        })

    print("\n[4/4] Writing reports")
    summary_fields = [
        "workflow","nominal_fdr","effective_r","accepted_unique_IL_peptide_keys",
        "target_only","entrapment_only","shared_excluded","unmatched_qc",
        "evaluable_T_plus_E","lower_bound_fdp","combined_fdp",
        "lower_bound_pct","combined_fdp_pct","interpretation"
    ]
    write_tsv(outdir / "Step02B_entrapment_summary.tsv", summary_rows, summary_fields)
    write_tsv(outdir / "Step02B_discovery_classes.tsv", discovery_rows,
              ["workflow","class","peptide_key_IL_equivalent"])
    write_tsv(outdir / "Step02B_source_files.tsv", source_rows,
              ["workflow","source_file","size_mb"])
    write_tsv(outdir / "Step02A_reaudit.tsv", qc_rows,
              ["metric","expected_v103","observed","pass"])

    report = [
        f"Step02B HOMEWORK CHECK v{VERSION}",
        f"Root: {root}",
        f"Nominal peptide FDR: {100*args.alpha:.2f}%",
        f"Effective r: {effective_r:.12f}",
        "",
        "workflow\taccepted_unique_IL\tT\tE\tshared\tunmatched\tlower_FDP_pct\tcombined_FDP_pct\tinterpretation",
    ]
    for row in summary_rows:
        report.append(
            f"{row['workflow']}\t{row['accepted_unique_IL_peptide_keys']}\t"
            f"{row['target_only']}\t{row['entrapment_only']}\t"
            f"{row['shared_excluded']}\t{row['unmatched_qc']}\t"
            f"{row['lower_bound_pct']:.6f}\t{row['combined_fdp_pct']:.6f}\t"
            f"{row['interpretation']}"
        )

    report_path = outdir / "Step02B_HOMEWORK_REPORT.txt"
    report_path.write_text("\n".join(report), encoding="utf-8")

    print("\n" + "=" * 100)
    print("HOMEWORK CHECK COMPLETE")
    print("=" * 100)
    for row in summary_rows:
        print(
            f"{row['workflow']}: lower={row['lower_bound_pct']:.4f}% | "
            f"combined={row['combined_fdp_pct']:.4f}% | {row['interpretation']}"
        )
    print("\nPaste or upload:")
    print(report_path)
    print(outdir / "Step02B_entrapment_summary.tsv")
    print(outdir / "Step02B_source_files.tsv")
    print("=" * 100)

if __name__ == "__main__":
    main()
