from pathlib import Path
import csv,json,hashlib,math,shutil
p=Path(__file__).resolve().parent;src=p/'source-data/sid';out=p/'dist'
rows=list(csv.DictReader(open(src/'microbiome_long.tsv'),delimiter='\t'))
rows=[r for r in rows if '-stool-' in r['sample_id']]
ids=sorted({r['sample_id'] for r in rows});features={}
for r in rows:
 sid=r['metaphlan_id'];f={k:r[k] for k in ['metaphlan_id','genus','species','family','phylum']}
 if sid in features:assert features[sid]==f
 features[sid]=f
features=sorted(features.values(),key=lambda r:int(r['metaphlan_id']))
dates=['2025-10-22','2025-12-02','2026-02-13','2026-03-10','2026-03-27']
notes=['Date agrees between the public summary and lab report.','Date agrees between the public summary and lab report.','Provisional: relative-day metadata and folder give Feb 13; the lab report and abundance filename give Feb 19.','Receipt date in the lab report; actual collection date unavailable.','Receipt date in the lab report; actual collection date unavailable.']
samples=[]
for i,s in enumerate(ids):
 rr={r['metaphlan_id']:float(r['sgb_abundance']) for r in rows if r['sample_id']==s};assert len(rr)==sum(r['sample_id']==s for r in rows)
 weights=[rr.get(f['metaphlan_id'],0) for f in features];total=sum(weights);assert total>0 and all(v>=0 and math.isfinite(v) for v in weights)
 shares=[v/total for v in weights];shannon=-sum(v*math.log(v) for v in shares if v);simpson=1/sum(v*v for v in shares)
 assert round(shannon,2)==[3.02,3.75,3.65,3.84,3.58][i]
 assert round(simpson,2)==[11.14,24.78,13.92,20.81,21.65][i]
 samples.append({'id':s,'date':dates[i],'dateNote':notes[i],'dateStatus':['reported','reported','conflicting','received','received'][i],'alternateDate':'2026-02-19' if i==2 else None,'weights':weights,'total':total,'shares':shares,'observed':len(rr),'shannon':shannon})
# Validate downloaded summary against original lab CSV for first visit.
orig=list(csv.DictReader(open(src/'original-stool-001.csv')));original={r['SGB_ID']:float(r['RCRF-H37004-Collection-1']) for r in orig}
assert all(original.get(f['metaphlan_id'],0)==samples[0]['weights'][j] for j,f in enumerate(features))
prov={'source':'https://osteosarc.com/microbiome/','sourceTable':'https://osteosarc.com/data/microbiome/microbiome_long.tsv','report':'https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/microbiome/NGS/RCRF-H37004-20260422-Stool/processed_data/Stool_report_4_22_26.pdf','retrieved':'2026-10-17','software':'MetaPhlAn 4.1.1','database':'mpa_vJun23_CHOCOPhlAnSGB_202307','pipelineSource':'June 23, 2026 lab report covering the first six reported visits. Full preprocessing parameters are not documented here.','measure':'Share of reported SGB abundance weights; each supplied weight divided by the sum across all reported stool SGBs for that visit. Underlying integer units are unverified; these are not asserted to be raw reads, cells, or absolute microbial load.','validation':'First visit matches original lab CSV. Normalized weights reproduce all five lab-report Shannon and inverse-Simpson values at two decimals.','missing':'An absent SGB row is treated as not detected (zero in the supplied profile), not proven biological absence.','scope':'Five stool profiles in the downloadable public summary. Later archive collections have not been imported. No tumor matching, external reference, p-values or clinical scoring.','sourceSHA256':{name:hashlib.sha256((src/name).read_bytes()).hexdigest() for name in ['microbiome_long.tsv','samples.tsv','original-stool-001.csv','april-report.pdf']}}
(out/'sid-data.json').write_text(json.dumps({'features':features,'samples':samples,'provenance':prov},separators=(',',':')))
(out/'sid-provenance.json').write_text(json.dumps(prov,indent=2))
with open(out/'sid-stool-profiles.tsv','w') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
print('Prepared',len(features),'SGBs across',len(samples),'visits; lab summary checks passed.')
