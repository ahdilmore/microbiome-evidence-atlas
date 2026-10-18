"""Reconcile public HMP2Data counts with the Stanford abundance table.
Run with Python, rdata, numpy and pandas. Does not modify the app.
"""
from pathlib import Path
import json, re, hashlib, html
import numpy as np
import pandas as pd
import rdata

ROOT = Path(__file__).resolve().parent
SRC = ROOT / 'source'
out = ROOT / 'output'
out.mkdir(exist_ok=True)
mtx = rdata.read_rda(SRC / 'T2D16S_mtx.rda')['T2D16S_mtx']
meta = rdata.read_rda(SRC / 'T2D16S_samp.rda')['T2D16S_samp']
parsed = rdata.parser.parse_file(SRC / 'T2D16S_tax.rda')
attrs = rdata.conversion.SimpleConverter().convert(parsed.object.value[0].attributes)
tax_values = rdata.read_rda(SRC / 'T2D16S_tax.rda')['T2D16S_tax']
# The reader returns an R character matrix as a flat vector; restore its declared dimensions.
tax = pd.DataFrame(tax_values.reshape(tuple(attrs['dim']), order='F'),
                   index=attrs['dimnames'][0], columns=attrs['dimnames'][1])
counts = pd.DataFrame(mtx.values, index=mtx.coords['dim_0'].values,
                      columns=mtx.coords['dim_1'].values)
assert counts.index.equals(tax.index)
assert counts.columns.equals(meta.index)
assert np.isfinite(counts.values).all()
assert (counts.values >= 0).all()
assert (counts.values == np.floor(counts.values)).all()
counts = counts.astype('int64')
meta = meta.copy()
meta['visit_token'] = meta['file_name'].str.extract(r'_(Z[A-Z0-9]+-[^_]+)_[^_]+$', expand=False)
assert meta['visit_token'].notna().all()
fecal = meta.loc[meta['sample_body_site'] == 'feces'].copy()
fecal['library_size'] = counts[fecal.index].sum(axis=0)
assert (fecal.library_size > 0).all()
stanford = pd.read_csv(ROOT.parent / 'data/gut_16s_abundance.tsv', sep='\t', index_col=0)
assert stanford.index.is_unique
groups = fecal.groupby('visit_token').size()
matches = pd.DataFrame(index=stanford.index)
matches['count_records'] = groups.reindex(stanford.index).fillna(0).astype(int)
matches['match_status'] = np.select([matches.count_records == 0, matches.count_records == 1],
                                   ['unmatched', 'unique_exact_visit_token'], default='multiple_count_records')
matches['count_file_names'] = [';'.join(fecal.index[fecal.visit_token == s]) for s in matches.index]
matches.to_csv(out/'sample_match_audit.csv', index_label='StanfordSampleID')
fecal.to_csv(out/'fecal_sample_metadata.csv', index_label='CountSampleID')
tax.to_csv(out/'otu_taxonomy.csv', index_label='OTU_ID')
counts[fecal.index].to_csv(out/'fecal_otu_counts.csv.gz', index_label='OTU_ID', compression='gzip')

unique_ids = matches.index[matches.count_records == 1]
sample_lookup = fecal.reset_index(names='count_name').set_index('visit_token')
# Select only one-to-one exact tokens: no suffix edits, pooling or replicate assumptions.
count_names = [fecal.index[fecal.visit_token == s][0] for s in unique_ids]
c = counts[count_names].copy()
c.columns = unique_ids
denom = c.sum(axis=0)
known = tax.Genus.notna() & (tax.Genus != '')
genus_counts = c.loc[known].groupby(tax.loc[known, 'Genus']).sum()
relative = genus_counts.divide(denom, axis=1)
original_names = [x[6:] for x in stanford.columns if x.startswith('genus_')]
overlap = sorted(set(original_names) & set(genus_counts.index))
tax_audit = pd.DataFrame({'StanfordGenus': original_names})
tax_audit['mapping'] = ['exact_name' if x in overlap else 'no_exact_genus_name' for x in original_names]
tax_audit.to_csv(out/'genus_mapping_audit.csv', index=False)
records=[]
for genus in overlap:
 a=stanford.loc[unique_ids,'genus_'+genus].to_numpy(float)
 b=relative.loc[genus,unique_ids].to_numpy(float)
 d=np.abs(a-b)
 records.append({'genus':genus,'samples':len(a),'mean_absolute_difference_pp':float(d.mean()*100),
  'median_absolute_difference_pp':float(np.median(d)*100),'max_absolute_difference_pp':float(d.max()*100),
  'fraction_within_1e-6':float(np.mean(d<=1e-6)),
  'pearson_correlation':float(np.corrcoef(a,b)[0,1]) if np.std(a)>0 and np.std(b)>0 else None})
comparison=pd.DataFrame(records)
comparison.to_csv(out/'genus_comparison.csv',index=False)
relative.loc[overlap].T.to_csv(out/'matched_count_derived_relative_abundance.csv',index_label='SampleID')
genus_counts.T.to_csv(out/'matched_classified_genus_counts.csv.gz', index_label='SampleID',compression='gzip')
pd.DataFrame({'SampleID':unique_ids,'CountSampleID':count_names,'total_retained_OTU_reads':denom.values,
              'genus_unassigned_reads':c.loc[~known].sum(axis=0).values}).to_csv(out/'matched_library_sizes.csv',index=False)
all_a=stanford.loc[unique_ids,['genus_'+x for x in overlap]].to_numpy(float)
all_b=relative.loc[overlap,unique_ids].T.to_numpy(float)
diff=np.abs(all_a-all_b)
sample_metadata=pd.read_csv(ROOT.parent/'data/sample_metadata.csv',keep_default_na=False)
eligible=sample_metadata[(sample_metadata.visit_match=='exact') & (sample_metadata.visit_Event=='Healthy') & (sample_metadata.subject_IR_IS_classification.isin(['IR','IS']))]
eligible_matched=eligible[eligible.SampleID.isin(unique_ids)]
fully_reproduced=bool(np.all(diff<=1e-6))
report={
 'audit_date':'2026-10-17',
 'sources':{f.name:{'url':'https://raw.githubusercontent.com/dozmorovlab/HMP2Data/master/data/'+f.name,
                    'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in SRC.glob('*.rda')},
 'count_shape':list(counts.shape),'nonnegative_integer_counts':True,
 'body_sites':{str(k):int(v) for k,v in meta.sample_body_site.value_counts().items()},
 'fecal_unique_visit_tokens':int(fecal.visit_token.nunique()),
 'fecal_duplicated_visit_tokens':int((groups>1).sum()),
 'fecal_library_size_quantiles':{str(q):float(fecal.library_size.quantile(q)) for q in [0,.25,.5,.75,1]},
 'stanford_samples':len(stanford),'match_status_counts':{str(k):int(v) for k,v in matches.match_status.value_counts().items()},
 'stanford_named_and_unclassified_genus_features':len(original_names),'exact_shared_genus_names':len(overlap),
 'unmapped_stanford_genus_names':sorted(set(original_names)-set(overlap)),
 'comparison':{'unique_matched_samples':len(unique_ids),'sample_genus_pairs':int(diff.size),
  'mean_absolute_difference_percentage_points':float(diff.mean()*100),
  'median_absolute_difference_percentage_points':float(np.median(diff)*100),
  'max_absolute_difference_percentage_points':float(diff.max()*100),
  'fraction_pairs_within_1e-6':float(np.mean(diff<=1e-6)),
  'samples_all_shared_genera_within_1e-6':int(np.all(diff<=1e-6,axis=1).sum())},
 'default_IR_IS_cohort':{'original_healthy_samples':len(eligible),'original_participants':int(eligible.SubjectID.nunique()),'uniquely_matched_healthy_samples':len(eligible_matched),'uniquely_matched_participants':int(eligible_matched.SubjectID.nunique())},
 'matched_libraries_below_1000_retained_reads':int((denom<1000).sum()),
 'decision':('Shared genera reproduce within tolerance; still verify provenance and replicate identity before replacement.' if fully_reproduced else 'Do not replace the app input: these counts do not reproduce the Stanford abundance table on shared exact-name genera.'),
 'limitations':['Exact visit-token agreement supports matching but does not prove aliquot or sequencing-run identity.',
 'Multiple records per visit are excluded from quantitative comparison, not pooled.',
 'Genus mapping uses exact strings only; Greengenes taxonomy can split or merge taxa differently from Stanford labels.',
 'Relative abundance denominator is all retained OTU counts in each package sample, including OTUs without genus assignments; it is not verified to equal raw sequencing depth.',
 'Unassigned-genus reads are preserved separately and are not distributed among named genera.',
 'The package provides processed OTU read counts, not bacterial cell counts or unprocessed FASTQ reads.',
 'Differences may reflect processing, taxonomy, filtering or sample handling; this audit cannot uniquely attribute the cause.']}
(out/'audit_summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
