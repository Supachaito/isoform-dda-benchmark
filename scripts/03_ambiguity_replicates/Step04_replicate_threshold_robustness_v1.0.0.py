#!/usr/bin/env python3
"""
Step04_replicate_threshold_robustness_v1.0.0.py

Downstream-only robustness analysis for the frozen MBR-OFF benchmark.

NO SEARCH ENGINE RERUN.
NO FASTA REMAPPING.
NO CHANGE TO FROZEN Step01B CLASSIFICATION.

Uses:
  AP_MBR_OFF/*.ms_data.hdf
  FP_MBR_OFF/<replicate>/psm.tsv
  MM_MBR_OFF/Task1-SearchTask/AllPSMs.psmtsv
  MQ_MBR_OFF/txt/msms.txt
  Step01_isoform_resolvability_package/Step01B_observed_results_v111/
      observed_full_FASTA_classification.tsv

Replicate-threshold definition
------------------------------
For each workflow and each biological cell-line group (C33A, HELA, SIHA),
a peptide is supported at:
  >=1/3 replicates
  >=2/3 replicates
  3/3 replicates

The workflow-level threshold set is the UNION of peptides satisfying that
threshold in at least one of the three cell lines.

This preserves condition-specific replicate support and avoids an arbitrary
pooled 9-run threshold.

Outputs:
  - validation against frozen Step01B workflow unions
  - threshold summary
  - structural ambiguity composition by threshold
  - cross-workflow support of primary isoform-discriminative peptides
  - leave-one-replicate-out robustness
  - figures

Stopping rule
-------------
If all validation anchors MATCH and outputs are generated, Step 4 is complete.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

VERSION = "1.0.0"
WORKFLOWS = ["AP", "FP", "MM", "MQ"]
CELL_LINES = ["C33A", "HELA", "SIHA"]
REPLICATES = [1, 2, 3]

EXPECTED_COMMON = {"AP":23767, "FP":22674, "MM":24543, "MQ":18962}
EXPECTED_PRIMARY = {"AP":180, "FP":101, "MM":149, "MQ":70}
EXPECTED_PRIMARY_UNION = 297
EXPECTED_PRIMARY_INCIDENCES = 500

VALID_CATEGORIES = {
    "single_canonical_unique",
    "single_isoform_unique",
    "within_family_subset_discriminative",
    "within_family_shared_all",
    "same_gene_subset_discriminative",
    "same_gene_multi_entry_shared",
    "multi_entry_gene_unresolved",
    "cross_gene_shared",
}

GROUP_ORDER = [
    "exact_entry_unique",
    "within_family_discriminative",
    "within_family_unresolved",
    "same_gene_discriminative",
    "same_gene_unresolved",
    "cross_gene_shared",
]

def iv(v):
    try:
        return int(float(str(v).strip()))
    except Exception:
        return 0

def fv(v):
    try:
        return float(str(v).strip())
    except Exception:
        return None

def normcol(s):
    return re.sub(r"[^a-z0-9]+","",str(s).strip().lower())

def findcol(fields, aliases):
    lookup={normcol(x):x for x in fields}
    for a in aliases:
        k=normcol(a)
        if k in lookup:
            return lookup[k]
    return None

def clean_seq(v):
    if v is None:
        return None
    if isinstance(v, bytes):
        v=v.decode("utf-8",errors="replace")
    s=str(v).strip().upper()
    if not s:
        return None
    # flanked peptide representation X.PEPTIDE.Y
    m=re.fullmatch(r"[A-Z\-]\.([A-Z]+)\.[A-Z\-]",s)
    if m:
        s=m.group(1)
    if not re.fullmatch(r"[A-Z]+",s):
        return None
    return s.replace("I","J").replace("L","J")

def truthy(v):
    s=str(v or "").strip().lower()
    return s in {"1","true","t","yes","y","+","decoy","contaminant"}

def target_value(v):
    s=str(v or "").strip().lower()
    if s in {"t","target","true","1"}:
        return True
    if s in {"d","decoy","c","contaminant","false","0"}:
        return False
    return None

def write_tsv(path, rows, fields=None):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    if not rows:
        path.write_text("",encoding="utf-8")
        return
    if fields is None:
        fields=list(rows[0].keys())
    with path.open("w",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

def read_table(path):
    path=Path(path)
    delim="," if path.suffix.lower()==".csv" else "\t"
    with path.open("r",encoding="utf-8-sig",errors="replace",newline="") as fh:
        yield from csv.DictReader(fh,delimiter=delim)

def rep_id_from_text(text):
    s=Path(str(text)).name
    s=re.sub(r"\.(raw|mzml|mzxml|mgf)$","",s,flags=re.I)
    s=s.upper().replace("-","_").replace(" ","_")
    # normalize HELA / C33A / SIHA
    for cell in CELL_LINES:
        m=re.search(rf"{cell}.*?([123])(?:$|_)",s)
        if m:
            return f"{cell}_{m.group(1)}"
    # fallback exact ending
    m=re.search(r"(C33A|HELA|SIHA)[^0-9]*([123])$",s)
    if m:
        return f"{m.group(1)}_{m.group(2)}"
    return None

def structural_group(cat):
    if cat in {"single_canonical_unique","single_isoform_unique"}:
        return "exact_entry_unique"
    if cat=="within_family_subset_discriminative":
        return "within_family_discriminative"
    if cat=="within_family_shared_all":
        return "within_family_unresolved"
    if cat=="same_gene_subset_discriminative":
        return "same_gene_discriminative"
    if cat in {"same_gene_multi_entry_shared","multi_entry_gene_unresolved"}:
        return "same_gene_unresolved"
    if cat=="cross_gene_shared":
        return "cross_gene_shared"
    return "other"

def load_frozen(step01b_dir):
    path=Path(step01b_dir)/"observed_full_FASTA_classification.tsv"
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r",encoding="utf-8-sig",errors="replace",newline="") as fh:
        r=csv.DictReader(fh,delimiter="\t")
        rows=[
            x for x in r
            if str(x.get("category","")).strip() in VALID_CATEGORIES
        ]
    by_key={}
    for x in rows:
        k=str(x.get("peptide_key_IL_equivalent","")).strip()
        if k:
            by_key[k]=x
    return path, rows, by_key

# --------------------------- AP ---------------------------

def _decode(v):
    if isinstance(v,bytes):
        return v.decode("utf-8",errors="replace")
    try:
        return v.item()
    except Exception:
        return v

def ap_peptide_fdr_rows(path):
    import h5py
    blocks=[]
    with h5py.File(path,"r") as h5:
        cands=[]
        def visitor(name,obj):
            if "peptide_fdr" in name.lower():
                cands.append((name,obj))
        h5.visititems(visitor)
        for name,obj in cands:
            if isinstance(obj,h5py.Group):
                cols={}
                lens=[]
                for cn,ch in obj.items():
                    if isinstance(ch,h5py.Dataset) and ch.ndim==1:
                        arr=ch[()]
                        cols[cn]=arr
                        lens.append(len(arr))
                if cols and len(set(lens))==1:
                    rows=[]
                    for i in range(lens[0]):
                        rows.append({k:_decode(v[i]) for k,v in cols.items()})
                    blocks.append(rows)
            elif isinstance(obj,h5py.Dataset) and obj.dtype.names:
                arr=obj[()]
                rows=[{k:_decode(rec[k]) for k in obj.dtype.names} for rec in arr]
                blocks.append(rows)
    if not blocks:
        raise RuntimeError(f"No peptide_fdr block found: {path}")
    return blocks

def extract_ap(root):
    d=Path(root)/"AP_MBR_OFF"
    out={}
    files=sorted(d.glob("*_*.ms_data.hdf"))
    if len(files)!=9:
        raise RuntimeError(f"Expected 9 AP MBR-OFF HDF files, found {len(files)}")
    for path in files:
        rid=rep_id_from_text(path.stem.replace(".ms_data",""))
        if not rid:
            raise RuntimeError(f"Cannot infer AP replicate from {path.name}")
        seqs=set()
        for rows in ap_peptide_fdr_rows(path):
            if not rows:
                continue
            fields=list(rows[0].keys())
            seq_col=findcol(fields,[
                "sequence_naked","naked_sequence","peptide_sequence","peptide","sequence"
            ])
            dec_col=findcol(fields,["decoy","is_decoy","decoy_flag"])
            tgt_col=findcol(fields,["target","is_target","target_decoy"])
            if not seq_col:
                continue
            for row in rows:
                if dec_col and truthy(row.get(dec_col)):
                    continue
                if tgt_col and target_value(row.get(tgt_col)) is False:
                    continue
                s=clean_seq(row.get(seq_col))
                if s:
                    seqs.add(s)
        out[rid]=seqs
    return out

# --------------------------- FP ---------------------------

def extract_fp(root):
    d=Path(root)/"FP_MBR_OFF"
    out={}
    for cell in CELL_LINES:
        for rep in REPLICATES:
            rid=f"{cell}_{rep}"
            path=d/rid/"psm.tsv"
            if not path.exists():
                raise FileNotFoundError(path)
            rows=read_table(path)
            first=next(rows,None)
            if first is None:
                out[rid]=set(); continue
            fields=list(first.keys())
            seq_col=findcol(fields,["Peptide","Peptide Sequence","Sequence"])
            dec_col=findcol(fields,["Is Decoy","Decoy"])
            con_col=findcol(fields,["Is Contaminant","Contaminant"])
            if not seq_col:
                raise RuntimeError(f"FP sequence column missing: {path}")
            seqs=set()
            for row in [first,*rows]:
                if dec_col and truthy(row.get(dec_col)): continue
                if con_col and truthy(row.get(con_col)): continue
                s=clean_seq(row.get(seq_col))
                if s: seqs.add(s)
            out[rid]=seqs
    return out

# --------------------------- MM ---------------------------

def extract_mm(root):
    path=Path(root)/"MM_MBR_OFF"/"Task1-SearchTask"/"AllPSMs.psmtsv"
    if not path.exists():
        raise FileNotFoundError(path)
    rows=read_table(path)
    first=next(rows,None)
    if first is None:
        raise RuntimeError("MM AllPSMs is empty")
    fields=list(first.keys())
    seq_col=findcol(fields,["Base Sequence","BaseSequence","Peptide"])
    q_col=findcol(fields,["QValue","Q Value","q-value","qvalue"])
    tdc_col=findcol(fields,["Decoy/Contaminant/Target","Decoy Contaminant Target","DCT"])
    file_col=findcol(fields,[
        "File Name","FileName","File","Filename","Spectra File","SpectraFile"
    ])
    if not all([seq_col,q_col,tdc_col,file_col]):
        raise RuntimeError(
            f"MM required columns missing. seq={seq_col}, q={q_col}, "
            f"tdc={tdc_col}, file={file_col}; fields={fields}"
        )
    out={f"{c}_{r}":set() for c in CELL_LINES for r in REPLICATES}
    for row in [first,*rows]:
        q=fv(row.get(q_col))
        if q is None or q>0.01: continue
        if target_value(row.get(tdc_col)) is not True: continue
        rid=rep_id_from_text(row.get(file_col))
        if rid not in out: continue
        raw=str(row.get(seq_col,"") or "")
        if "|" in raw: continue
        s=clean_seq(raw)
        if s: out[rid].add(s)
    return out

# --------------------------- MQ ---------------------------

def extract_mq(root):
    path=Path(root)/"MQ_MBR_OFF"/"txt"/"msms.txt"
    if not path.exists():
        raise FileNotFoundError(path)
    rows=read_table(path)
    first=next(rows,None)
    if first is None:
        raise RuntimeError("MQ msms.txt is empty")
    fields=list(first.keys())
    seq_col=findcol(fields,["Sequence"])
    raw_col=findcol(fields,["Raw file","Rawfile","Raw File"])
    rev_col=findcol(fields,["Reverse"])
    con_col=findcol(fields,["Potential contaminant","Potential Contaminant"])
    if not seq_col or not raw_col:
        raise RuntimeError(
            f"MQ required columns missing. seq={seq_col}, raw={raw_col}; fields={fields}"
        )
    out={f"{c}_{r}":set() for c in CELL_LINES for r in REPLICATES}
    for row in [first,*rows]:
        if rev_col and truthy(row.get(rev_col)): continue
        if con_col and truthy(row.get(con_col)): continue
        rid=rep_id_from_text(row.get(raw_col))
        if rid not in out: continue
        s=clean_seq(row.get(seq_col))
        if s: out[rid].add(s)
    return out

# --------------------------- analysis ---------------------------

def filter_to_frozen(rep_sets, frozen_keys):
    return {rid:(seqs & frozen_keys) for rid,seqs in rep_sets.items()}

def threshold_set(rep_sets, threshold):
    result=set()
    for cell in CELL_LINES:
        counts=Counter()
        for rep in REPLICATES:
            rid=f"{cell}_{rep}"
            counts.update(rep_sets[rid])
        result.update(k for k,n in counts.items() if n>=threshold)
    return result

def primary_set(keys, frozen):
    return {
        k for k in keys
        if iv(frozen[k].get("primary_isoform_discriminative"))==1
    }

def validation_rows(rep_by_wf, frozen):
    out=[]
    def add(metric,obs,exp):
        out.append({
            "metric":metric,
            "observed":obs,
            "expected":exp,
            "difference":obs-exp,
            "status":"MATCH" if obs==exp else "CHECK"
        })
    for wf in WORKFLOWS:
        add(f"replicate_count_{wf}",len(rep_by_wf[wf]),9)
        union=set().union(*rep_by_wf[wf].values())
        add(f"common_reference_union_{wf}",len(union),EXPECTED_COMMON[wf])
        add(f"primary_union_{wf}",len(primary_set(union,frozen)),EXPECTED_PRIMARY[wf])

    psets=[]
    incid=0
    for wf in WORKFLOWS:
        u=set().union(*rep_by_wf[wf].values())
        p=primary_set(u,frozen)
        psets.append(p); incid+=len(p)
    add("primary_crossworkflow_union",len(set().union(*psets)),EXPECTED_PRIMARY_UNION)
    add("primary_workflow_peptide_incidences",incid,EXPECTED_PRIMARY_INCIDENCES)
    return out

def replicate_summary(rep_by_wf, frozen):
    rows=[]
    for wf in WORKFLOWS:
        for rid in sorted(rep_by_wf[wf]):
            keys=rep_by_wf[wf][rid]
            p=primary_set(keys,frozen)
            cell,rep=rid.split("_")
            rows.append({
                "workflow":wf,
                "cell_line":cell,
                "replicate":int(rep),
                "n_common_reference_peptide_keys":len(keys),
                "n_primary_isoform_discriminative":len(p),
                "primary_pct":100*len(p)/len(keys) if keys else math.nan
            })
    return rows

def threshold_summary(rep_by_wf, frozen):
    rows=[]
    sets={}
    for wf in WORKFLOWS:
        for th in [1,2,3]:
            keys=threshold_set(rep_by_wf[wf],th)
            p=primary_set(keys,frozen)
            sets[(wf,th)]=keys
            groups=Counter(structural_group(frozen[k]["category"]) for k in keys)
            rows.append({
                "workflow":wf,
                "replicate_threshold":th,
                "threshold_label":f">={th}/3 in >=1 cell line" if th<3 else "3/3 in >=1 cell line",
                "n_common_reference_peptide_keys":len(keys),
                "pct_retained_vs_ge1":
                    math.nan, # filled below
                "n_primary_isoform_discriminative":len(p),
                "primary_pct_of_common":100*len(p)/len(keys) if keys else math.nan,
                "n_exact_entry_unique":groups["exact_entry_unique"],
                "n_within_family_discriminative":groups["within_family_discriminative"],
                "n_within_family_unresolved":groups["within_family_unresolved"],
                "n_same_gene_discriminative":groups["same_gene_discriminative"],
                "n_same_gene_unresolved":groups["same_gene_unresolved"],
                "n_cross_gene_shared":groups["cross_gene_shared"],
            })
    base={(r["workflow"]):r["n_common_reference_peptide_keys"]
          for r in rows if r["replicate_threshold"]==1}
    for r in rows:
        b=base[r["workflow"]]
        r["pct_retained_vs_ge1"]=100*r["n_common_reference_peptide_keys"]/b if b else math.nan
    return rows,sets

def group_threshold_rows(threshold_sets,frozen):
    rows=[]
    for (wf,th),keys in threshold_sets.items():
        denom=len(keys)
        cnt=Counter(structural_group(frozen[k]["category"]) for k in keys)
        for grp in GROUP_ORDER:
            n=cnt[grp]
            rows.append({
                "workflow":wf,
                "replicate_threshold":th,
                "ambiguity_group":grp,
                "n_peptide_keys":n,
                "total_threshold_peptides":denom,
                "pct_of_threshold_set":100*n/denom if denom else math.nan
            })
    return rows

def crossworkflow_primary_support(threshold_sets,frozen):
    rows=[]
    summary=[]
    for th in [1,2,3]:
        psets={wf:primary_set(threshold_sets[(wf,th)],frozen) for wf in WORKFLOWS}
        union=set().union(*psets.values())
        support=Counter(sum(k in psets[wf] for wf in WORKFLOWS) for k in union)
        for n in [1,2,3,4]:
            c=support[n]
            rows.append({
                "replicate_threshold":th,
                "n_supporting_workflows":n,
                "n_unique_primary_peptide_keys":c,
                "threshold_primary_union":len(union),
                "pct_of_threshold_primary_union":100*c/len(union) if union else math.nan
            })
        summary.append({
            "replicate_threshold":th,
            "primary_union_unique_keys":len(union),
            "workflow_specific_n":support[1],
            "workflow_specific_pct":100*support[1]/len(union) if union else math.nan,
            "supported_ge2_n":support[2]+support[3]+support[4],
            "supported_ge2_pct":
                100*(support[2]+support[3]+support[4])/len(union) if union else math.nan,
            "supported_all4_n":support[4],
            "supported_all4_pct":100*support[4]/len(union) if union else math.nan,
            "AP_n":len(psets["AP"]),
            "FP_n":len(psets["FP"]),
            "MM_n":len(psets["MM"]),
            "MQ_n":len(psets["MQ"]),
        })
    return rows,summary

def leave_one_out(rep_by_wf,frozen):
    rows=[]
    for wf in WORKFLOWS:
        all_union=set().union(*rep_by_wf[wf].values())
        all_primary=primary_set(all_union,frozen)
        for rid in sorted(rep_by_wf[wf]):
            kept=[s for rr,s in rep_by_wf[wf].items() if rr!=rid]
            u=set().union(*kept)
            p=primary_set(u,frozen)
            rows.append({
                "workflow":wf,
                "removed_replicate":rid,
                "baseline_common_reference":len(all_union),
                "leave_one_out_common_reference":len(u),
                "common_retained_pct":
                    100*len(u)/len(all_union) if all_union else math.nan,
                "baseline_primary":len(all_primary),
                "leave_one_out_primary":len(p),
                "primary_retained_pct":
                    100*len(p)/len(all_primary) if all_primary else math.nan,
                "primary_lost_n":len(all_primary-p),
            })
    return rows

def pairwise_primary_jaccard(threshold_sets,frozen):
    rows=[]
    for th in [1,2,3]:
        p={wf:primary_set(threshold_sets[(wf,th)],frozen) for wf in WORKFLOWS}
        for a,b in combinations(WORKFLOWS,2):
            inter=len(p[a]&p[b]); union=len(p[a]|p[b])
            rows.append({
                "replicate_threshold":th,
                "workflow_A":a,
                "workflow_B":b,
                "n_A":len(p[a]),
                "n_B":len(p[b]),
                "intersection":inter,
                "union":union,
                "jaccard":inter/union if union else math.nan,
                "jaccard_pct":100*inter/union if union else math.nan,
            })
    return rows

def make_figures(outdir,threshold_rows,cross_summary,loro_rows):
    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[warning] matplotlib unavailable: {e}")
        return

    # Fig 1 primary counts vs replicate threshold
    fig,ax=plt.subplots(figsize=(7.2,4.8))
    for wf in WORKFLOWS:
        rr=sorted([r for r in threshold_rows if r["workflow"]==wf],
                  key=lambda x:x["replicate_threshold"])
        ax.plot(
            [r["replicate_threshold"] for r in rr],
            [r["n_primary_isoform_discriminative"] for r in rr],
            marker="o",label=wf
        )
    ax.set_xticks([1,2,3])
    ax.set_xticklabels([">=1/3",">=2/3","3/3"])
    ax.set_xlabel("Replicate support within at least one cell line")
    ax.set_ylabel("Primary isoform-discriminative peptide keys")
    ax.set_title("Replicate-threshold robustness of isoform-discriminative evidence")
    ax.legend(frameon=False); ax.grid(axis="y",alpha=.25)
    fig.tight_layout()
    fig.savefig(outdir/"Fig_Step04_primary_by_replicate_threshold.png",dpi=300)
    fig.savefig(outdir/"Fig_Step04_primary_by_replicate_threshold.pdf")
    plt.close(fig)

    # Fig 2 workflow-specific fraction across thresholds
    fig,ax=plt.subplots(figsize=(6.5,4.6))
    xs=[r["replicate_threshold"] for r in cross_summary]
    ys=[r["workflow_specific_pct"] for r in cross_summary]
    ax.plot(xs,ys,marker="o")
    ax.set_xticks([1,2,3]); ax.set_xticklabels([">=1/3",">=2/3","3/3"])
    ax.set_xlabel("Replicate support within at least one cell line")
    ax.set_ylabel("Workflow-specific primary peptides (%)")
    ax.set_title("Workflow specificity under stricter replicate support")
    ax.grid(axis="y",alpha=.25)
    fig.tight_layout()
    fig.savefig(outdir/"Fig_Step04_workflow_specificity_by_threshold.png",dpi=300)
    fig.savefig(outdir/"Fig_Step04_workflow_specificity_by_threshold.pdf")
    plt.close(fig)

    # Fig 3 leave-one-out primary retention
    fig,ax=plt.subplots(figsize=(8.0,4.8))
    data=[]
    labels=[]
    for wf in WORKFLOWS:
        vals=[r["primary_retained_pct"] for r in loro_rows if r["workflow"]==wf]
        data.append(vals); labels.append(wf)
    ax.boxplot(data,labels=labels,showfliers=True)
    ax.set_ylabel("Primary evidence retained after removing one run (%)")
    ax.set_title("Leave-one-replicate-out robustness")
    ax.grid(axis="y",alpha=.25)
    fig.tight_layout()
    fig.savefig(outdir/"Fig_Step04_leave_one_out_primary_retention.png",dpi=300)
    fig.savefig(outdir/"Fig_Step04_leave_one_out_primary_retention.pdf")
    plt.close(fig)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",required=True)
    ap.add_argument("--step01b-dir",default=None)
    ap.add_argument("--outdir",default=None)
    args=ap.parse_args()

    root=Path(args.root).resolve()
    step01b=(
        Path(args.step01b_dir).resolve()
        if args.step01b_dir
        else root/"Step01_isoform_resolvability_package"/"Step01B_observed_results_v111"
    )
    outdir=(
        Path(args.outdir).resolve()
        if args.outdir
        else root/"Step04_replicate_threshold_robustness_v100"
    )
    outdir.mkdir(parents=True,exist_ok=True)

    print(f"[1/7] Step04 v{VERSION}: loading frozen structural classification")
    frozen_path,frozen_rows,frozen=load_frozen(step01b)
    frozen_keys=set(frozen)
    print(f"      frozen common-reference keys: {len(frozen_keys):,}")

    print("[2/7] Extracting replicate-level MBR-OFF peptide evidence (no rerun)")
    raw={}
    raw["AP"]=extract_ap(root)
    raw["FP"]=extract_fp(root)
    raw["MM"]=extract_mm(root)
    raw["MQ"]=extract_mq(root)

    print("[3/7] Intersecting replicate evidence with frozen full-FASTA common reference")
    rep_by_wf={wf:filter_to_frozen(raw[wf],frozen_keys) for wf in WORKFLOWS}

    print("[4/7] Validating workflow unions against frozen Step01B anchors")
    val=validation_rows(rep_by_wf,frozen)
    write_tsv(outdir/"Step04_validation.tsv",val)
    bad=[r for r in val if r["status"]!="MATCH"]
    if bad:
        for r in bad: print("CHECK",r)
        raise RuntimeError(
            "Step04 validation failed. Replicate extraction does not reproduce "
            "the frozen Step01B workflow unions; do not interpret outputs."
        )
    print("      All anchors MATCH.")

    print("[5/7] Calculating replicate-threshold and leave-one-out robustness")
    rep_sum=replicate_summary(rep_by_wf,frozen)
    th_sum,th_sets=threshold_summary(rep_by_wf,frozen)
    grp=group_threshold_rows(th_sets,frozen)
    support,support_sum=crossworkflow_primary_support(th_sets,frozen)
    loro=leave_one_out(rep_by_wf,frozen)
    jac=pairwise_primary_jaccard(th_sets,frozen)

    print("[6/7] Writing outputs")
    write_tsv(outdir/"Step04_replicate_summary.tsv",rep_sum)
    write_tsv(outdir/"Step04_threshold_summary.tsv",th_sum)
    write_tsv(outdir/"Step04_structural_groups_by_threshold.tsv",grp)
    write_tsv(outdir/"Step04_crossworkflow_primary_support.tsv",support)
    write_tsv(outdir/"Step04_crossworkflow_primary_summary.tsv",support_sum)
    write_tsv(outdir/"Step04_leave_one_replicate_out.tsv",loro)
    write_tsv(outdir/"Step04_primary_pairwise_jaccard_by_threshold.tsv",jac)

    manifest={
        "script_version":VERSION,
        "root":str(root),
        "frozen_classification":str(frozen_path),
        "no_search_rerun":True,
        "no_fasta_remapping":True,
        "replicate_threshold_definition":
            "Within each cell line separately, peptide support is >=1/3, >=2/3, or 3/3; workflow-level set is union across C33A/HELA/SIHA.",
        "leave_one_out_definition":
            "Remove one of the nine MBR-OFF runs from a workflow and recompute the workflow union.",
        "stopping_rule":
            "If all validation anchors MATCH and outputs are generated, Step 4 is complete."
    }
    (outdir/"Step04_manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")

    print("[7/7] Generating figures")
    make_figures(outdir,th_sum,support_sum,loro)

    print("\nTHRESHOLD SUMMARY")
    for r in th_sum:
        print(
            f"  {r['workflow']} >= {r['replicate_threshold']}/3: "
            f"common={r['n_common_reference_peptide_keys']}, "
            f"primary={r['n_primary_isoform_discriminative']} "
            f"({r['primary_pct_of_common']:.3f}%)"
        )

    print("\nCROSS-WORKFLOW PRIMARY SUPPORT")
    for r in support_sum:
        print(
            f"  >= {r['replicate_threshold']}/3: union={r['primary_union_unique_keys']}, "
            f"workflow-specific={r['workflow_specific_pct']:.1f}%, "
            f">=2 workflows={r['supported_ge2_pct']:.1f}%, "
            f"all4={r['supported_all4_pct']:.1f}%"
        )

    print(f"\nSTEP 4 COMPLETE: {outdir}")

if __name__=="__main__":
    main()
