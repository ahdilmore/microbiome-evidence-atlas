"""Summarize completed HUMAnN runs; never substitute taxonomic predictions."""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parent
SPECIAL = {'UNMAPPED', 'UNINTEGRATED', 'UNGROUPED'}

def read_table(path):
    values = {}
    with path.open() as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)
        for row in reader:
            if not row: continue
            if len(row) != 2 or row[0] in values:
                raise ValueError('Expected unique single-sample rows in ' + str(path))
            val = float(row[1])
            if not math.isfinite(val) or val < 0:
                raise ValueError('Invalid abundance')
            values[row[0]] = val
    return values

def summarize(samples, tables, minimum):
    features = sorted(set().union(*(set(t) for t in tables)))
    rows = []
    for feature in features:
        if '|' in feature or feature in SPECIAL: continue
        vals = [t.get(feature, 0.) for t in tables]
        med = statistics.median(vals)
        deltas = [b-a for a,b in zip(vals, vals[1:])]
        rows.append(dict(pathway=feature, abundance_cpm=vals,
            visits_above_threshold=sum(v > 0 and v >= minimum for v in vals),
            median_cpm=med, range_cpm=max(vals)-min(vals),
            range_over_median=(max(vals)-min(vals))/med if med else None,
            first_to_last_delta_cpm=vals[-1]-vals[0], adjacent_deltas_cpm=deltas,
            largest_adjacent_absolute_delta_cpm=max(map(abs, deltas))))
    persistent = [r for r in rows if r['visits_above_threshold']==len(samples)]
    persistent.sort(key=lambda r: (r['range_over_median'], -r['median_cpm']))
    return dict(samples=samples, units='HUMAnN copies per million, special features included in denominator',
        minimum_cpm=minimum, pathways=rows,
        most_changed=sorted(rows,key=lambda r:r['largest_adjacent_absolute_delta_cpm'],reverse=True)[:30],
        persistent_least_variable=persistent[:30],
        special_features=[{k:v for k,v in t.items() if k in SPECIAL} for t in tables],
        interpretation='Descriptive DNA-based functional potential. Persistence is not proof of biological stability. Missing pathway rows represent non-detection, not proven absence. No causal inference.')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--minimum-cpm',type=float,default=10)
    p.add_argument('--output',type=Path,default=ROOT/'pathway-summary.json')
    a=p.parse_args()
    if not math.isfinite(a.minimum_cpm) or a.minimum_cpm<0: p.error('Invalid threshold')
    samples=json.loads((ROOT/'samples.json').read_text()); tables=[]; settings=[]
    for s in samples:
        folder=ROOT/'work'/s['id']
        provenance=json.loads((folder/'run.json').read_text())
        if provenance['status']!='complete': raise ValueError('Incomplete sample '+s['id'])
        settings.append({k:v for k,v in provenance['settings'].items() if k not in ['sample','threads']})
        tables.append(read_table(folder/'humann'/f"{s['id']}_pathabundance_cpm.tsv"))
    if any(x!=settings[0] for x in settings): raise ValueError('Inconsistent profiling settings')
    a.output.write_text(json.dumps(summarize(samples,tables,a.minimum_cpm),indent=2,allow_nan=False)+'\n')
    print(a.output)

if __name__=='__main__': main()
