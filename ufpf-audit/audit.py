"""Reproduce the UFPF/Sid compatibility audit without modifying source data."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
m = pd.read_csv(ROOT/'Metadata All Samples.csv')
m.columns = m.columns.str.strip()
a = pd.read_csv(ROOT/'metaphlan_all_combined.tsv', sep='\t', comment='#', index_col=0)
a.columns = a.columns.str.replace('_metaphlan_rel_ab','',regex=False)
b = pd.read_csv(ROOT/'metaphlan_unknown_all_combined.tsv', sep='\t', comment='#', index_col=0)
b.columns = b.columns.str.replace('_metaphlan_rel_ab_unknown','',regex=False)
sid = json.loads((ROOT.parent/'site/dist/sid-data.json').read_text())
assert not a.index.duplicated().any() and not a.columns.duplicated().any()
assert not m.Sample.duplicated().any()
assert set(a.columns) == set(b.columns)
assert np.isfinite(a.to_numpy()).all() and (a.to_numpy() >= 0).all()
assert np.isfinite(b.to_numpy()).all() and (b.to_numpy() >= 0).all()
assert a.index.str.contains(r'\|t__').any()
t = a.loc[a.index.str.contains(r'\|t__SGB')].copy()
t.index = t.index.str.extract(r't__SGB(\d+(?:_group)?)$',expand=False)
assert not t.index.isna().any() and not t.index.duplicated().any()
leaves = a.loc[a.index.str.contains(r'\|t__')]
assert np.allclose(leaves.sum(),100,atol=.001)
overlap = set(t.index) & {x['metaphlan_id'] for x in sid['features']}
joined = m[m.Sample.isin(a.columns)]
normalized_diagnosis = m.Diagnosis.replace({'CD':'IBD','UC':'IBD','CD/UC':'IBD'})
conflict = m.loc[normalized_diagnosis != m.Diagnosis2, ['Sample','Diagnosis','Diagnosis2']]
groups={}
for label,g in joined.groupby('Diagnosis2'):
    groups[label] = dict(n=len(g), age_min=float(g.Age.min()),age_median=float(g.Age.median()),
      age_max=float(g.Age.max()),sex=g.Sex.value_counts().to_dict(),
      reads_min=int(g.Reads.min()),reads_median=float(g.Reads.median()),reads_max=int(g.Reads.max()))
common = a.index.intersection(b.index)
error = (a.loc[common]-b.loc[common,a.columns].div(1-b.loc['UNCLASSIFIED',a.columns]/100,axis=1)).abs().max().max()
assert error < .0001
issues=[]
for sample in sorted(set(a.columns)-set(m.Sample)):
    issues.append(dict(sample=sample,issue='profile_without_exact_metadata',action='Exclude pending verified ID resolution'))
for sample in sorted(set(m.Sample)-set(a.columns)):
    issues.append(dict(sample=sample,issue='metadata_without_exact_profile',action='Do not infer or silently repair ID'))
for row in conflict.to_dict('records'):
    issues.append(dict(sample=row['Sample'],issue=f"diagnosis_conflict:{row['Diagnosis']} vs {row['Diagnosis2']}",action='Exclude from diagnosis-stratified primary comparisons; include only as ambiguous in all-profile exploration'))
with (ROOT/'sample-issues.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=['sample','issue','action']);w.writeheader();w.writerows(issues)
summary=dict(repository_commit=json.loads((ROOT/'repository-tree.json').read_text())['sha'],
    profiles=a.shape[1],metadata_rows=len(m),exact_matches=len(joined),
    groups_as_Diagnosis2=groups,diagnosis_conflicts=conflict.to_dict('records'),
    conservative_labeled_profiles=len(joined)-len(set(conflict.Sample)&set(joined.Sample)),
    profiles_without_metadata=sorted(set(a.columns)-set(m.Sample)),
    metadata_without_profiles=sorted(set(m.Sample)-set(a.columns)),
    database_header=(ROOT/'metaphlan_all_combined.tsv').read_text().splitlines()[0],
    taxonomic_rows=len(a),sgb_rows=len(t),non_sgb_leaf_rows=len(leaves)-len(t),
    grouped_sgb_rows=int(t.index.str.endswith('_group').sum()),
    grouped_ids_also_seen_as_numeric_sid_ids=sorted(set(x.removesuffix('_group') for x in t.index if x.endswith('_group')) & {x['metaphlan_id'] for x in sid['features']}),
    leaf_sums_percent=[float(leaves.sum().min()),float(leaves.sum().max())],
    sgb_sums_percent=[float(t.sum().min()),float(t.sum().max())],
    unknown_percent_range=[float(b.loc['UNCLASSIFIED'].min()),float(b.loc['UNCLASSIFIED'].max())],
    unknown_normalization_max_error_percentage_points=float(error),
    sid_sgbs=len(sid['features']),shared_sgbs=len(overlap),
    sid_abundance_fraction_in_shared_sgbs=[dict(sample=s['id'],date=s['date'],fraction=sum(v for f,v in zip(sid['features'],s['shares']) if f['metaphlan_id'] in overlap)) for s in sid['samples']],
    missing_sid_sgbs=[f for f in sid['features'] if f['metaphlan_id'] not in overlap],
    missing_metadata_counts=m.isna().sum().to_dict(),
    interpretation='Same database confirmed in actual headers and commands. Profiling minor version and Sid preprocessing not fully verified. Overlap measures detection across supplied tables, not sensitivity or proof of biological absence. Cross-sectional exploratory reference only.')
z=pd.read_excel(ROOT/'Source Data.xlsx',sheet_name='UFPF_metadata')
za=pd.read_excel(ROOT/'Source Data.xlsx',sheet_name='UFPF_metaphlan_relab',index_col=0)
z.Sample_ID=z.Sample_ID.str.replace('.','-',regex=False)
za.columns=za.columns.str.replace('.','-',regex=False)
assert not z.Sample_ID.duplicated().any() and set(za.columns)==set(z.Sample_ID)
common=a.columns.intersection(za.columns)
archive_error=float(np.max(np.abs(a.loc[za.index,common].values-za.loc[:,common].values)))
assert archive_error==0
alias_error=float(np.max(np.abs(a.loc[za.index,'UF-PF-2021-023-1N']-za['UF-PF-2020-023-1N'])))
assert alias_error==0
zt=za.loc[za.index.str.contains(r'\|t__SGB')].copy()
zt.index=zt.index.str.extract(r't__SGB(\d+(?:_group)?)$',expand=False)
archive_observed=set(zt.index[zt.sum(axis=1)>0])
archive_overlap=archive_observed & {x['metaphlan_id'] for x in sid['features']}
summary['archive'] = dict(metadata_rows=len(z),profiles=za.shape[1],exact_matches=len(z),
    groups=z.Diagnosis2.value_counts().to_dict(),
    conservative_diagnosis_labeled_profiles=len(z)-len(conflict),
    source='https://zenodo.org/records/10912505',
    common_profile_max_error=archive_error,year_alias_profile_max_error=alias_error,
    verified_profile_alias={'github':'UF-PF-2021-023-1N','archive':'UF-PF-2020-023-1N'},
    github_extra_profiles=sorted(set(a.columns)-set(za.columns)-{'UF-PF-2021-023-1N'}),
    github_unknown_probiotics=int(m['Taking Probiotics'].isna().sum()),
    archive_unknown_probiotics=int(z['Taking.Probiotics'].isna().sum()),
    primary_metadata_policy='Use archive IDs and profiles; retain disputed diagnosis as ambiguous; preserve GitHub missing probiotics as unknown, not No.',
    observed_sgb_rows=len(archive_observed),shared_observed_sgbs=len(archive_overlap),
    sid_abundance_fraction_in_shared_observed_sgbs=[sum(v for f,v in zip(sid['features'],s['shares']) if f['metaphlan_id'] in archive_overlap) for s in sid['samples']])
summary['file_sha256']={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in ROOT.iterdir() if f.name in ['Source Data.xlsx','Metadata All Samples.csv','metaphlan_all_combined.tsv','metaphlan_unknown_all_combined.tsv','paper.xml','Key Resource Table.xlsx','Bioinformatic Processing of Sequences']}
(ROOT/'audit-results.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:summary[k] for k in ['profiles','metadata_rows','exact_matches','conservative_labeled_profiles','shared_sgbs','leaf_sums_percent','unknown_normalization_max_error_percentage_points']},indent=2))
