"""Build the browser dataset from the separately audited HMP2Data count files."""
from pathlib import Path
import json, hashlib, shutil
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent
OUT=ROOT/'dist'
SRC=BASE/'count-audit/output'
counts=pd.read_csv(SRC/'fecal_otu_counts.csv.gz',index_col=0)
tax=pd.read_csv(SRC/'otu_taxonomy.csv',index_col=0,keep_default_na=False)
meta=pd.read_csv(SRC/'fecal_sample_metadata.csv',index_col=0,keep_default_na=False)
visits=pd.read_csv(BASE/'data/visit_metadata.csv',index_col=0,keep_default_na=False)
subjects=pd.read_csv(BASE/'data/subject_metadata.csv',index_col=0,keep_default_na=False)
assert counts.index.equals(tax.index) and counts.columns.equals(meta.index)
assert visits.index.is_unique and subjects.index.is_unique
assert np.isfinite(counts.values).all() and (counts.values>=0).all()
assert (counts.values==np.floor(counts.values)).all()
totals=counts.sum(axis=0)
assert np.array_equal(totals.to_numpy(),meta.library_size.to_numpy())
named=tax.Genus!=''
genus=counts.loc[named].groupby(tax.loc[named,'Genus']).sum()
genus=genus.loc[genus.sum(axis=1)>0].sort_index()
relative=genus.divide(totals,axis=1)
unknown=counts.loc[~named].sum(axis=0)
assert np.array_equal(genus.sum(axis=0).to_numpy()+unknown.to_numpy(),totals.to_numpy())
duplicates=set(meta.visit_token[meta.visit_token.duplicated(keep=False)])
records=[];audit=[]
for name,m in meta.iterrows():
 token=m.visit_token
 reason='eligible'
 if '_ST_' not in name:reason='body_site_filename_conflict'
 elif token in duplicates:reason='ambiguous_duplicate_visit'
 elif totals[name]<1000:reason='below_1000_retained_reads'
 elif token not in visits.index:reason='no_exact_visit_metadata'
 elif visits.loc[token,'SubjectID'] not in subjects.index:reason='no_subject_metadata'
 audit.append({'CountSampleID':name,'VisitID':token,'retained_OTU_reads':int(totals[name]),'status':reason})
 if reason!='eligible':continue
 v=visits.loc[token];sid=v.SubjectID
 assert sid==token.split('-')[0]
 s=subjects.loc[sid]
 records.append({'id':token,'countSampleId':name,'subject':sid,'match':'exact','event':v.Event,
  'race':s.Race,'study':s.Study,'visitSubstudy':v.SubStudy,'collectionDate':v.CollectionDate,
  'ir':s.IR_IS_classification,'sex':s.Sex,'age':s.Age,'bmi':s.BMI,'sspg':s.SSPG,
  'librarySize':int(totals[name]),'unassignedGenusReads':int(unknown[name]),
  'counts':[int(x) for x in genus[name]],'values':[float(x) for x in relative[name]]})
quality=pd.DataFrame(audit)
healthy=[r for r in records if r['event']=='Healthy']
default=[r for r in healthy if r['ir'] in ['IR','IS']]
groups={r['subject']:r['ir'] for r in default}
provenance={'dataset':'HMP2Data T2D16S','taxonomy':'Greengenes 13_5','version':'greengenes-v1',
 'source':'https://github.com/dozmorovlab/HMP2Data','retrieved':'2026-10-17',
 'rawFecalSamples':len(meta),'otuFeatures':len(counts),'namedGenusFeatures':len(genus),
 'minimumRetainedReads':1000,'exclusions':{str(k):int(v) for k,v in quality.status.value_counts().items()},
 'eligibleSamples':len(records),'eligibleParticipants':len({r['subject'] for r in records}),
 'healthySamples':len(healthy),'healthyParticipants':len({r['subject'] for r in healthy}),
 'defaultSamples':len(default),'defaultParticipants':len(groups),'defaultIR':sum(v=='IR' for v in groups.values()),'defaultIS':sum(v=='IS' for v in groups.values()),
 'denominator':'All retained OTU counts per fecal sample, including OTUs without a genus assignment.',
 'duplicatePolicy':'Exclude every count record with a duplicated fecal visit token; do not pool or choose replicates.',
 'analysis':'Median per participant across Healthy visits; unadjusted rank tests on count-derived relative abundances; not CLR or DREAM.',
 'thresholdNote':'1000 retained reads is a transparent exploratory QC choice, not a threshold reproduced from the source paper.',
 'sourceSHA256':{f:hashlib.sha256((SRC/f).read_bytes()).hexdigest() for f in ['fecal_otu_counts.csv.gz','otu_taxonomy.csv','fecal_sample_metadata.csv']}}
(OUT/'data.json').write_text(json.dumps({'taxa':genus.index.tolist(),'samples':records,'provenance':provenance},separators=(',',':')))
(OUT/'data-provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
quality.to_csv(OUT/'sample-qc.csv',index=False)
shutil.copyfile(SRC/'fecal_otu_counts.csv.gz',OUT/'fecal-otu-counts.csv.gz')
shutil.copyfile(SRC/'otu_taxonomy.csv',OUT/'otu-taxonomy.csv')
print(json.dumps(provenance,indent=2))
