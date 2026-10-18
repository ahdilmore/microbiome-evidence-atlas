"""Read-only compatibility audit; no cross-cohort disease or distance estimates.
Run with Python plus numpy/pandas. Inputs and checksums are in sources.json.
"""
from pathlib import Path
import bz2, csv, hashlib, json, re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
SID = ROOT.parent / 'site/dist/sid-data.json'
TREE = ROOT.parent / 'site/source-data/sid/mpa_vJun23_CHOCOPhlAnSGB_202307.nwk'
def stats(x):
    x=np.asarray(x,dtype=float)
    return dict(zip(['min','q25','median','q75','max'],map(float,np.quantile(x,[0,.25,.5,.75,1]))))

def main():
    raw=(ROOT/'merged_table.tsv').read_text()
    header=next(csv.reader(raw.splitlines(),delimiter='\t'))
    assert len(header)==len(set(header)), 'Duplicate headers'
    frame=pd.read_csv(ROOT/'merged_table.tsv',sep='\t')
    columns=[c for c in frame if c.startswith('t__')]
    x=frame[columns].to_numpy(float)
    assert np.isfinite(x).all() and (x>=0).all()
    sums=x.sum(axis=1)
    assert (sums>0).all()
    sid=json.loads(SID.read_text())
    sidids=[str(f['metaphlan_id']) for f in sid['features']]
    sv=np.array([s['shares'] for s in sid['samples']])
    assert np.allclose(sv.sum(axis=1),1)
    maps={}
    for f in sorted(ROOT.glob('*_species.txt.bz2')):
        maps[f.name.split('_')[1]]={line.split('\t',1)[0]:line.split('\t',1)[1] for line in bz2.decompress(f.read_bytes()).decode().splitlines() if '\t' in line}
    features=[c.removeprefix('t__') for c in columns]
    numeric=[re.fullmatch(r'SGB(\d+)',f).group(1) if re.fullmatch(r'SGB(\d+)',f) else None for f in features]
    grouped={re.fullmatch(r'SGB(\d+)_group',f).group(1) for f in features if re.fullmatch(r'SGB(\d+)_group',f)}
    tips=set(re.findall(r'(?:\(|,)(\d+):',TREE.read_text()))
    assert len(tips)>30000
    # Mapping is purely structural; it is NOT evidence of the generating DB or identical marker sets.
    keep=np.array([f is not None and f in tips and f not in grouped for f in numeric])
    sidkeep=np.array([f in tips and f not in grouped for f in sidids])
    retained=x[:,keep].sum(axis=1)/sums
    cohort={}
    for group,g in frame.groupby('cohortID'):
        i=g.index.to_numpy()
        cohort[group]={'profiles':len(i),'sum_percent':stats(sums[i]),'nonzero_features':stats((x[i]>0).sum(axis=1)),
                       'tree_mapping_retained_percent':stats(retained[i]*100),'profiles_below_95_percent_retained':int((retained[i]<.95).sum())}
    coverage=[]
    for i,s in enumerate(sid['samples']):
        coverage.append({'sample':s['id'],'date':s['date'],'tree_mapping_retained_percent':float(sv[i,sidkeep].sum()*100),
                         'share_in_exact_numeric_columns_observed_in_reference_percent':float(sv[i,[f in set(numeric) for f in sidids]].sum()*100)})
    memberships={}
    for version,taxonomy in maps.items():
        missing=[f for f in features if f not in taxonomy]
        memberships[version]={'exact_feature_matches':len(features)-len(missing),'unmatched_count':len(missing),'unmatched':missing,
                              'sid_numeric_matches':sum('SGB'+f in taxonomy for f in sidids)}
    # Hypothetical name repair is a diagnostic only; source columns are never modified.
    bases=[re.sub(r'\.\.\.\d+$','',f) for f in features]
    collisions={b:[j for j,v in enumerate(bases) if v==b] for b in set(bases) if bases.count(b)>1}
    cooccurring={b:int(((x[:,j]>0).sum(axis=1)>1).sum()) for b,j in collisions.items()}
    fingerprint={}
    ambiguous_union=set()
    for version,taxonomy in maps.items():
        normalized={k.removesuffix('_group') for k in taxonomy}
        group_alias={b for b in set(bases) if b not in taxonomy and b+'_group' in taxonomy}
        ambiguous_union.update(group_alias)
        unsupported=x[:,[b not in normalized for b in bases]].sum(axis=1)
        fingerprint[version]={'absent_unique_ids_after_suffix_and_group_normalization':sorted(set(bases)-normalized),
          'bare_labels_matching_only_group_entries':sorted(group_alias),
          'profiles_whose_nonzero_ids_all_exist_in_catalog':int((unsupported==0).sum()),
          'unsupported_percent_by_cohort':{g:stats(unsupported[i]) for g,i in frame.groupby('cohortID').groups.items()}}
        sample_col='unsupported_'+version+'_percent'
        frame[sample_col]=unsupported
    alias_keep=np.array([b.startswith('SGB') and b[3:] in tips for b in bases])
    conservative_keep=np.array([k and b not in ambiguous_union for k,b in zip(alias_keep,bases)])
    conservative_sid=np.array(['SGB'+f not in ambiguous_union and f in tips for f in sidids])
    scenarios={}
    for name,mask in [('hypothetical_suffix_repair_tree_only',alias_keep),('hypothetical_suffix_repair_excluding_group_ambiguity_across_catalogs',conservative_keep)]:
        c=x[:,mask].sum(axis=1)/sums*100
        scenarios[name]={'retained_percent':stats(c),'profiles_below_95_percent':int((c<95).sum()),
                         'by_cohort':{g:stats(c[i]) for g,i in frame.groupby('cohortID').groups.items()}}
    with (ROOT/'duplicate-base-labels.tsv').open('w') as out:
        w=csv.writer(out,delimiter='\t');w.writerow(['base_label','source_columns','profiles_with_multiple_nonzero_columns'])
        for b,j in sorted(collisions.items()):w.writerow([b,';'.join(columns[i] for i in j),cooccurring[b]])
    # Per-feature diagnostics allow future mapping work without silently rewriting IDs.
    with (ROOT/'feature-audit.tsv').open('w') as out:
        w=csv.writer(out,delimiter='\t');w.writerow(['feature','nonzero_profiles','mean_percent','max_percent','numeric_tree_match',*maps])
        for j,f in enumerate(features):w.writerow([f,int((x[:,j]>0).sum()),float(x[:,j].mean()),float(x[:,j].max()),bool(keep[j]),*[f in m for m in maps.values()]])
    sample=frame[['sampleID','cohortID']].copy()
    sample['sum_percent']=sums;sample['retained_percent']=retained*100
    sample['nonzero_features']=(x>0).sum(axis=1)
    sample.to_csv(ROOT/'sample-audit.tsv',sep='\t',index=False)
    result={'scope':'Structural audit only. No disease classification, health score or clinical similarity estimated.',
      'source_sha256':hashlib.sha256((ROOT/'merged_table.tsv').read_bytes()).hexdigest(),
      'sid_sha256':hashlib.sha256(SID.read_bytes()).hexdigest(),'tree_sha256':hashlib.sha256(TREE.read_bytes()).hexdigest(),
      'profiles':len(frame),'unique_sample_ids':int(frame.sampleID.nunique()),'microbial_features':len(columns),
      'metadata_columns':[c for c in header if c not in columns],'metadata_missing':frame[[c for c in header if c not in columns]].isna().sum().to_dict(),
      'nonfinite_or_negative_values':False,'sum_percent':stats(sums),'rows_with_sum_outside_99_9_to_100_1':int(((sums<99.9)|(sums>100.1)).sum()),
      'all_zero_features':int((x.sum(axis=0)==0).sum()),'duplicate_abundance_profiles':int(frame[columns].duplicated().sum()),
      'cohorts':cohort,'database_membership':memberships,
      'suffix_diagnostics':{'suffixed_columns':sum(b!=f for b,f in zip(bases,features)),'unique_base_labels':len(set(bases)),
        'duplicated_base_labels':len(collisions),'base_labels_with_cooccurring_nonzeros':sum(v>0 for v in cooccurring.values()),
        'interpretation':'Looks consistent with name-repair/merge artifacts, but origin and equivalence are not proven. Do not strip and pool silently.'},
      'catalog_fingerprint_diagnostic':fingerprint,
      'hypothetical_mapping_scenarios_not_approved_for_analysis':scenarios,
      'sid_conservative_mapping_retained_percent':(sv[:,conservative_sid].sum(axis=1)*100).tolist(),
      'mapping_policy':'Exact numeric SGB IDs in the existing Jun23 tree, excluding grouped labels and their numeric base IDs symmetrically. Coverage is relative to each supplied profile total, not total reads. This does NOT validate cross-version marker equivalence.',
      'grouped_labels':[f for f in features if f.endswith('_group')],
      'other_non_numeric':[f for f,n in zip(features,numeric) if n is None and not f.endswith('_group')],
      'numeric_ids_not_in_sid_tree':[f for f in numeric if f is not None and f not in tips],
      'reference_features_retained':int(keep.sum()),'sid_features_retained':int(sidkeep.sum()),
      'sid_coverage':coverage,'exact_numeric_observed_overlap':len(set(numeric)&set(sidids)),
      'figure5_csv_identical_to_deposit':(ROOT/'merged_table.tsv').read_bytes()==(ROOT/'figure5.csv').read_bytes(),
      'database_verdict':'Unconfirmed: deposited merged table lacks native MetaPhlAn version/database/command headers. Feature membership cannot identify the generating database.',
      'analysis_verdict':'Hold joint PCoA/weighted UniFrac and cancer-group ranking. Resolve duplicate-like suffixes, stripped group semantics, per-sample database provenance and metadata first. Catalog overlap is not marker equivalence.'}
    (ROOT/'audit-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('profiles','microbial_features','sum_percent','metadata_columns','suffix_diagnostics','hypothetical_mapping_scenarios_not_approved_for_analysis','sid_conservative_mapping_retained_percent','analysis_verdict')},indent=2))
    print('DB membership:',{k:{a:b for a,b in v.items() if a!='unmatched'} for k,v in memberships.items()})
    print('Grouped labels:',len(grouped))

if __name__=='__main__':main()
