#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path
from collections import Counter, defaultdict
from itertools import combinations

VERSION='1.0.1'
WF=['AP','FP','MM','MQ']
EXPECTED_COMMON={'AP':23767,'FP':22674,'MM':24543,'MQ':18962}
EXPECTED_PRIMARY={'AP':180,'FP':101,'MM':149,'MQ':70}
EXPECTED_UNION=297
EXPECTED_INCID=500

VALID_CATEGORIES={
 'single_canonical_unique',
 'single_isoform_unique',
 'within_family_subset_discriminative',
 'within_family_shared_all',
 'same_gene_subset_discriminative',
 'same_gene_multi_entry_shared',
 'multi_entry_gene_unresolved',
 'cross_gene_shared'
}

GROUPS=[
 'exact_entry_unique',
 'within_family_discriminative',
 'within_family_unresolved',
 'same_gene_discriminative',
 'same_gene_unresolved',
 'cross_gene_shared',
 'other'
]

def iv(x):
    try:return int(float(str(x).strip()))
    except:return 0

def grp(cat):
    if cat in {'single_canonical_unique','single_isoform_unique'}: return 'exact_entry_unique'
    if cat=='within_family_subset_discriminative': return 'within_family_discriminative'
    if cat=='within_family_shared_all': return 'within_family_unresolved'
    if cat=='same_gene_subset_discriminative': return 'same_gene_discriminative'
    if cat in {'same_gene_multi_entry_shared','multi_entry_gene_unresolved'}: return 'same_gene_unresolved'
    if cat=='cross_gene_shared': return 'cross_gene_shared'
    return 'other'

def write_tsv(p,rows,fields=None):
    p.parent.mkdir(parents=True,exist_ok=True)
    if not rows:
        p.write_text('',encoding='utf-8'); return
    if fields is None: fields=list(rows[0])
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(rows)

def load(path):
    with path.open('r',encoding='utf-8-sig',errors='replace',newline='') as f:
        r=csv.DictReader(f,delimiter='\t')
        need={'peptide_key_IL_equivalent','category','primary_isoform_discriminative',*WF}
        miss=need-set(r.fieldnames or [])
        if miss: raise RuntimeError(f'Missing required columns: {sorted(miss)}')
        return list(r)

def validate(rows):
    out=[]
    def add(m,o,e): out.append({'metric':m,'observed':o,'expected':e,'difference':o-e,'status':'MATCH' if o==e else 'CHECK'})
    for w in WF:
        add(f'common_reference_{w}',sum(iv(r[w])==1 for r in rows),EXPECTED_COMMON[w])
        add(f'primary_{w}',sum(iv(r[w])==1 and iv(r['primary_isoform_discriminative'])==1 for r in rows),EXPECTED_PRIMARY[w])
    prim=[r for r in rows if iv(r['primary_isoform_discriminative'])==1 and any(iv(r[w])==1 for w in WF)]
    add('primary_union',len(prim),EXPECTED_UNION)
    add('primary_workflow_peptide_incidences',sum(sum(iv(r[w]) for w in WF) for r in prim),EXPECTED_INCID)
    return out

def workflow_summary(rows):
    out=[]
    for w in WF:
        wr=[r for r in rows if iv(r[w])==1]; n=len(wr)
        c=Counter(grp(str(r['category']).strip()) for r in wr)
        primary=sum(iv(r['primary_isoform_discriminative'])==1 for r in wr)
        exact=c['exact_entry_unique']
        out.append({
            'workflow':w,'total_common_reference_peptides':n,
            'exact_entry_unique_n':exact,'exact_entry_unique_pct':100*exact/n,
            'primary_isoform_discriminative_n':primary,'primary_isoform_discriminative_pct':100*primary/n,
            'within_family_unresolved_n':c['within_family_unresolved'],'within_family_unresolved_pct':100*c['within_family_unresolved']/n,
            'same_gene_unresolved_n':c['same_gene_unresolved'],'same_gene_unresolved_pct':100*c['same_gene_unresolved']/n,
            'cross_gene_shared_n':c['cross_gene_shared'],'cross_gene_shared_pct':100*c['cross_gene_shared']/n,
            'all_non_exact_entry_n':n-exact,'all_non_exact_entry_pct':100*(n-exact)/n
        })
    return out

def grouped_by_workflow(rows):
    out=[]
    for w in WF:
        wr=[r for r in rows if iv(r[w])==1]; n=len(wr); c=Counter(grp(str(r['category']).strip()) for r in wr)
        for g in GROUPS:
            out.append({'workflow':w,'ambiguity_group':g,'n_peptide_keys':c[g],'workflow_total_peptide_keys':n,'pct_of_workflow_peptides':100*c[g]/n if n else math.nan})
    return out

def category_by_workflow(rows):
    cats=sorted(set(str(r['category']).strip() or 'other' for r in rows))
    out=[]
    for w in WF:
        wr=[r for r in rows if iv(r[w])==1]; n=len(wr); c=Counter(str(r['category']).strip() or 'other' for r in wr)
        for cat in cats:
            out.append({'workflow':w,'category':cat,'n_peptide_keys':c[cat],'workflow_total_peptide_keys':n,'pct_of_workflow_peptides':100*c[cat]/n if n else math.nan})
    return out

def support_by_group(rows):
    d=defaultdict(Counter); totals=Counter()
    for r in rows:
        s=sum(iv(r[w]) for w in WF)
        if s<1: continue
        g=grp(str(r['category']).strip()); d[g][s]+=1; totals[g]+=1
    out=[]
    for g in GROUPS:
        if totals[g]==0: continue
        for s in (1,2,3,4):
            n=d[g][s]
            out.append({'ambiguity_group':g,'n_supporting_workflows':s,'n_unique_peptide_keys':n,'group_total_unique_peptide_keys':totals[g],'pct_within_group':100*n/totals[g]})
    return out

def pairwise(rows):
    sets={}
    for w in WF:
        sets[(w,'ALL')]={r['peptide_key_IL_equivalent'] for r in rows if iv(r[w])==1}
        sets[(w,'PRIMARY')]={r['peptide_key_IL_equivalent'] for r in rows if iv(r[w])==1 and iv(r['primary_isoform_discriminative'])==1}
        for g in GROUPS:
            sets[(w,g)]={r['peptide_key_IL_equivalent'] for r in rows if iv(r[w])==1 and grp(str(r['category']).strip())==g}
    out=[]
    for scope in ['ALL','PRIMARY',*GROUPS]:
        for a,b in combinations(WF,2):
            A,B=sets[(a,scope)],sets[(b,scope)]; inter=len(A&B); union=len(A|B)
            out.append({'scope':scope,'workflow_A':a,'workflow_B':b,'n_A':len(A),'n_B':len(B),'intersection':inter,'union':union,'jaccard':inter/union if union else math.nan,'jaccard_pct':100*inter/union if union else math.nan})
    return out

def figures(outdir,grouped,support,pairs):
    try: import matplotlib.pyplot as plt
    except Exception as e:
        print('[warning] figures skipped:',e); return
    data={w:{g:0 for g in GROUPS} for w in WF}
    for r in grouped:data[r['workflow']][r['ambiguity_group']]=float(r['pct_of_workflow_peptides'])
    fig,ax=plt.subplots(figsize=(8.8,5.2)); bottom=[0]*4
    for g in GROUPS:
        vals=[data[w][g] for w in WF]; ax.bar(WF,vals,bottom=bottom,label=g.replace('_',' ')); bottom=[a+b for a,b in zip(bottom,vals)]
    ax.set_ylabel('Observed common-reference peptides (%)'); ax.set_title('Structural ambiguity composition across workflows'); ax.legend(frameon=False,bbox_to_anchor=(1.02,1),loc='upper left'); fig.tight_layout(); fig.savefig(outdir/'Fig_Step03_workflow_ambiguity_composition.png',dpi=300); fig.savefig(outdir/'Fig_Step03_workflow_ambiguity_composition.pdf'); plt.close(fig)

    groups=[g for g in GROUPS if any(r['ambiguity_group']==g for r in support)]; x=list(range(len(groups))); width=.18
    fig,ax=plt.subplots(figsize=(10,5.3))
    for idx,s in enumerate((1,2,3,4)):
        vals=[]
        for g in groups:
            m=[r for r in support if r['ambiguity_group']==g and iv(r['n_supporting_workflows'])==s]
            vals.append(float(m[0]['pct_within_group']) if m else 0)
        ax.bar([i+(idx-1.5)*width for i in x],vals,width=width,label=f'{s} workflow(s)')
    ax.set_xticks(x); ax.set_xticklabels([g.replace('_',' ') for g in groups],rotation=25,ha='right'); ax.set_ylabel('Unique peptide keys within class (%)'); ax.set_title('Cross-workflow reproducibility by structural class'); ax.legend(frameon=False); fig.tight_layout(); fig.savefig(outdir/'Fig_Step03_support_by_ambiguity_class.png',dpi=300); fig.savefig(outdir/'Fig_Step03_support_by_ambiguity_class.pdf'); plt.close(fig)

    mat=[[math.nan]*4 for _ in range(4)]
    for i in range(4):mat[i][i]=1
    for r in pairs:
        if r['scope']!='PRIMARY':continue
        i,j=WF.index(r['workflow_A']),WF.index(r['workflow_B']); v=float(r['jaccard']); mat[i][j]=mat[j][i]=v
    fig,ax=plt.subplots(figsize=(5.4,4.8)); im=ax.imshow(mat,vmin=0,vmax=1); ax.set_xticks(range(4)); ax.set_xticklabels(WF); ax.set_yticks(range(4)); ax.set_yticklabels(WF); ax.set_title('Primary isoform-discriminative peptide Jaccard')
    for i in range(4):
        for j in range(4):ax.text(j,i,f'{mat[i][j]:.2f}',ha='center',va='center')
    fig.colorbar(im,ax=ax,label='Jaccard'); fig.tight_layout(); fig.savefig(outdir/'Fig_Step03_primary_pairwise_jaccard.png',dpi=300); fig.savefig(outdir/'Fig_Step03_primary_pairwise_jaccard.pdf'); plt.close(fig)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',required=True)
    ap.add_argument('--step01b-dir',default=None)
    ap.add_argument('--outdir',default=None)
    a=ap.parse_args()
    root=Path(a.root).resolve()
    step01b=Path(a.step01b_dir).resolve() if a.step01b_dir else root/'Step01_isoform_resolvability_package'/'Step01B_observed_results_v111'
    inp=step01b/'observed_full_FASTA_classification.tsv'
    if not inp.exists(): raise SystemExit(f'Frozen Step01B file not found: {inp}')
    outdir=Path(a.outdir).resolve() if a.outdir else root/'Step03_ambiguity_decomposition_v101'; outdir.mkdir(parents=True,exist_ok=True)

    print('[1/6] Loading frozen Step01B')
    raw_rows=load(inp)
    print('      raw rows:',f'{len(raw_rows):,}')

    # Structural decomposition is defined only for the eight frozen full-FASTA
    # mapping categories. This intentionally removes raw observed keys that did
    # not map back to the common full-FASTA reference space (known MQ: 48).
    rows=[r for r in raw_rows if str(r.get('category','')).strip() in VALID_CATEGORIES]
    excluded=[r for r in raw_rows if str(r.get('category','')).strip() not in VALID_CATEGORIES]
    print('      structurally classified rows:',f'{len(rows):,}')
    print('      excluded unclassified/unmapped rows:',f'{len(excluded):,}')

    diag=[]
    for w in WF:
        n=sum(iv(r.get(w))==1 for r in excluded)
        diag.append({'workflow':w,'excluded_unclassified_or_unmapped_peptide_keys':n})
    write_tsv(outdir/'Step03_excluded_unclassified_or_unmapped.tsv',diag)

    print('[2/6] Validating frozen common-reference anchors')
    val=validate(rows); write_tsv(outdir/'Step03_validation.tsv',val)
    bad=[r for r in val if r['status']!='MATCH']
    if bad:
        for r in bad: print('CHECK',r)
        raise RuntimeError('Validation failed. Do not interpret Step03.')
    print('      All anchors MATCH.')
    print('[3/6] Decomposing ambiguity'); ws=workflow_summary(rows); gbw=grouped_by_workflow(rows); cbw=category_by_workflow(rows)
    print('[4/6] Cross-workflow support'); sbg=support_by_group(rows); pw=pairwise(rows)
    print('[5/6] Writing outputs'); write_tsv(outdir/'Step03_workflow_summary.tsv',ws); write_tsv(outdir/'Step03_grouped_ambiguity_by_workflow.tsv',gbw); write_tsv(outdir/'Step03_category_by_workflow.tsv',cbw); write_tsv(outdir/'Step03_crossworkflow_support_by_group.tsv',sbg); write_tsv(outdir/'Step03_pairwise_jaccard.tsv',pw)
    (outdir/'Step03_manifest.json').write_text(json.dumps({'version':VERSION,'input':str(inp),'no_rerun':True,'no_remapping':True,'input_filter':'Only the eight frozen full-FASTA structural categories; unclassified/unmapped raw keys excluded','groups':GROUPS,'stopping_rule':'If validation MATCH, Step 3 is complete.'},indent=2),encoding='utf-8')
    print('[6/6] Figures'); figures(outdir,gbw,sbg,pw)
    print('\nWORKFLOW SUMMARY')
    for r in ws:
        print(f"  {r['workflow']}: total={r['total_common_reference_peptides']}, exact={r['exact_entry_unique_pct']:.2f}%, primary={r['primary_isoform_discriminative_pct']:.3f}%, within-family-unresolved={r['within_family_unresolved_pct']:.2f}%, same-gene-unresolved={r['same_gene_unresolved_pct']:.2f}%, cross-gene={r['cross_gene_shared_pct']:.2f}%")
    print('\nSTEP 3 COMPLETE:',outdir)

if __name__=='__main__': main()
