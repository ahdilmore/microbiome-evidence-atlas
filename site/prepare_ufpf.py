"""Exploratory cross-study ordination; no hypothesis tests or health classification.

Inputs: archived UFPF Source Data.xlsx (Zenodo 10912505), existing Sid shares,
and the matching MetaPhlAn Newick. Requires numpy, pandas, openpyxl, Biopython.
Standalone use: python prepare_ufpf.py --archive 'Source Data.xlsx'
  --sid sid-data.json --tree mpa_vJun23_CHOCOPhlAnSGB_202307.nwk --output-dir results
Download source files from the URLs and compare SHA256 hashes in sid-ufpf.json.
"""
from pathlib import Path
import argparse, csv, hashlib, json, re
import numpy as np
import pandas as pd
from Bio import Phylo

ROOT = Path(__file__).resolve().parent

def weighted_distance(tree, features, values):
    lookup = {f:i for i,f in enumerate(features)}
    assert set(features) <= {t.name for t in tree.get_terminals()}
    assert np.isfinite(values).all() and (values >= 0).all()
    assert np.allclose(values.sum(axis=1), 1)
    masses = {}; numerator = np.zeros((len(values),len(values))); depth = np.zeros(len(values))
    for node in tree.find_clades(order='postorder'):
        mass = (values[:,lookup[node.name]] if node.name in lookup else np.zeros(len(values))) if node.is_terminal() else sum((masses.pop(c) for c in node.clades),np.zeros(len(values)))
        masses[node] = mass
        if node is tree.root: continue
        length = node.branch_length
        assert length is not None and np.isfinite(length) and length >= 0
        if np.any(mass):
            numerator += length * np.abs(mass[:,None]-mass[None,:])
            depth += length * mass
    den = depth[:,None]+depth[None,:]
    return np.divide(numerator, den, out=np.zeros_like(numerator), where=den>0)

def pcoa(distance):
    n = len(distance); j = np.eye(n)-np.ones((n,n))/n
    eig = np.linalg.eigvalsh(-.5*j@(distance**2)@j)
    correction = max(0.,-float(eig.min())) if eig.min() < -1e-12 else 0.
    squared = distance**2 + 2*correction
    np.fill_diagonal(squared,0)
    vals, vecs = np.linalg.eigh(-.5*j@squared@j)
    order = np.argsort(vals)[::-1]; vals = vals[order]; vecs = vecs[:,order]
    assert vals.min() > -1e-9
    positive = vals > 1e-12
    coordinates = vecs[:,positive] * np.sqrt(vals[positive])
    # Stable signs for reproducible plots; signs have no scientific meaning.
    for k in range(coordinates.shape[1]):
        if coordinates[np.argmax(np.abs(coordinates[:,k])),k] < 0: coordinates[:,k] *= -1
    reconstructed = np.sum((coordinates[:,None,:]-coordinates[None,:,:])**2,axis=2)
    assert np.allclose(reconstructed,squared,atol=1e-9)
    return dict(coordinates=coordinates[:,:3].tolist(),eigenvalues=vals.tolist(),
        axisFractions=(vals[:3]/vals[positive].sum()).tolist(),lingoesConstant=correction,
        rawEigenvalues=eig[::-1].tolist(),negativeInertiaFraction=float(-eig[eig<0].sum()/np.abs(eig).sum()))

def hull(points):
    pts=sorted(set(map(tuple,points)))
    def cross(o,a,b): return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    lower=[]; upper=[]
    for p in pts:
        while len(lower)>=2 and cross(lower[-2],lower[-1],p)<=0: lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper)>=2 and cross(upper[-2],upper[-1],p)<=0: upper.pop()
        upper.append(p)
    return np.array(lower[:-1]+upper[:-1])

def hulls_disjoint(a,b):
    # Separating-axis theorem for convex polygons, touching counts as overlap.
    for polygon in (a,b):
        for p,q in zip(polygon,np.roll(polygon,-1,axis=0)):
            axis=np.array([q[1]-p[1],p[0]-q[0]])
            x=a@axis; y=b@axis
            if x.max()<y.min()-1e-12 or y.max()<x.min()-1e-12: return True
    return False

def summary(distance):
    cross = distance[:5,5:]; within_sid = distance[:5,:5][np.triu_indices(5,1)]
    within_uf = distance[5:,5:][np.triu_indices(len(distance)-5,1)]
    stats=lambda v:dict(min=float(v.min()),median=float(np.median(v)),max=float(v.max()))
    return dict(between=stats(cross),withinSid=stats(within_sid),withinUF=stats(within_uf),
        strictDistanceGap=bool(cross.min()>max(within_sid.max(),within_uf.max())),
        nearestUFIndices=(cross.argmin(axis=1)+5).tolist(),nearestUFDistances=cross.min(axis=1).tolist())

def build(archive, sidpath=None, treepath=None, outputdir=None):
    sidpath=sidpath or ROOT/'dist/sid-data.json'
    treepath=treepath or ROOT/'source-data/sid/mpa_vJun23_CHOCOPhlAnSGB_202307.nwk'
    outputdir=outputdir or ROOT/'dist'; outputdir.mkdir(parents=True,exist_ok=True)
    sid=json.loads(sidpath.read_text()); tree=Phylo.read(treepath,'newick')
    tips=[t.name for t in tree.get_terminals()]; assert len(tips)==len(set(tips))
    metadata=pd.read_excel(archive,sheet_name='UFPF_metadata')
    metadata.Sample_ID=metadata.Sample_ID.str.replace('.','-',regex=False)
    abundance=pd.read_excel(archive,sheet_name='UFPF_metaphlan_relab',index_col=0)
    abundance.columns=abundance.columns.str.replace('.','-',regex=False)
    assert not metadata.Sample_ID.duplicated().any() and not abundance.index.duplicated().any()
    assert set(metadata.Sample_ID)==set(abundance.columns) and len(metadata)==96
    leaf=abundance.loc[abundance.index.str.contains(r'\|t__')]
    assert np.allclose(leaf.sum(),100,atol=.001)
    grouped=leaf.loc[leaf.index.str.endswith('_group')]
    ambiguous={re.search(r't__SGB(\d+)_group$',s).group(1) for s in grouped.index}
    uf=leaf.loc[leaf.index.str.contains(r'\|t__SGB\d+$')].copy()
    uf.index=uf.index.str.extract(r't__SGB(\d+)$',expand=False)
    sidids=[f['metaphlan_id'] for f in sid['features']]
    # Keep union, NOT just observed overlap. Symmetrically exclude ambiguous IDs.
    eligible=(set(uf.index)|set(sidids)) & set(tips) - ambiguous
    features=sorted(eligible,key=int); lookup={f:i for i,f in enumerate(features)}
    rows=[]; samples=[]
    for s in sid['samples']:
        row=np.zeros(len(features))
        for f,v in zip(sidids,s['shares']):
            if f in lookup: row[lookup[f]]=v
        rows.append(row)
        samples.append(dict(id=s['id'],cohort='Sid',group='Sid',date=s['date'],dateStatus=s['dateStatus'],dateNote=s['dateNote'],coverage=float(row.sum()),coverageDenominator='Supplied Sid SGB weights'))
    for m in metadata.to_dict('records'):
        row=np.zeros(len(features))
        for f,v in uf[m['Sample_ID']].items():
            if f in lookup: row[lookup[f]]=v/100
        rows.append(row)
        group='Uncertain diagnosis' if m['Sample_ID']=='UF-PF-2022-074-1N' else m['Diagnosis2']
        samples.append(dict(id=m['Sample_ID'],cohort='UFPF',group=group,age=int(m['Age']),sex=str(m['Sex']),reads=int(m['Reads']),coverage=float(row.sum()),coverageDenominator='Classified MetaPhlAn leaf abundance',bowelPrep=group=='IBD'))
    values=np.array(rows); coverage=values.sum(axis=1)
    assert (coverage>0).all() and (coverage<=1.00001).all()
    values=values/coverage[:,None]
    nonzero=values.sum(axis=0)>0; values=values[:,nonzero]; features=[f for f,keep in zip(features,nonzero) if keep]
    matrices={'unifrac':weighted_distance(tree,features,values)}
    matrices['braycurtis']=np.array([[np.abs(a-b).sum()/2 for b in values] for a in values])
    scenarios={}
    for group in ('all','nonIBD','Control','PD','IBD'):
        for floor in (0,95):
            indices=list(range(5))+[i for i,s in enumerate(samples) if i>=5 and (group=='all' or group=='nonIBD' and s['group']!='IBD' or s['group']==group) and s['coverage']>=floor/100]
            assert len(indices)>7
            result=dict(indices=indices,excludedReferenceCount=96-(len(indices)-5),metrics={})
            for key,d in matrices.items():
                matrix=d[np.ix_(indices,indices)]; ordination=pcoa(matrix); coords=np.array(ordination['coordinates'])
                ordination['hullsDisjoint']={f'0,{axis}':hulls_disjoint(hull(coords[:5,[0,axis]]),hull(coords[5:,[0,axis]])) for axis in (1,2)}
                result['metrics'][key]=dict(**ordination,summary=summary(matrix))
            scenarios[f'{group}-{floor}']=result
    hashes={str(p.name):hashlib.sha256(p.read_bytes()).hexdigest() for p in (archive,sidpath,treepath)}
    output=dict(generated='2026-10-17',database='mpa_vJun23_CHOCOPhlAnSGB_202307',samples=samples,features=features,matrices={k:v.tolist() for k,v in matrices.items()},scenarios=scenarios,
        provenance=dict(sha256=hashes,referenceSource='https://zenodo.org/records/10912505',paper='https://pmc.ncbi.nlm.nih.gov/articles/PMC11923181/',treeSource='https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202307.nwk',
        groupedRowsExcluded=grouped.index.tolist(),numericIdsExcludedSymmetrically=sorted(ambiguous,key=int),unmappedNumericIds=sorted((set(uf.index)|set(sidids))-set(tips),key=int),
        nonSGBrankRowsExcluded=leaf.loc[~leaf.index.str.contains(r'\|t__SGB')].index.tolist(),
        featurePolicy='Union of observed exact numeric SGB tip IDs. Exclude grouped SGB rows and corresponding numeric IDs in both cohorts; exclude non-SGB leaves and unmapped tips. Renormalize each retained profile to one. No pseudocounts, rarefaction, or shared-observation intersection.',
        metric='Normalized weighted UniFrac: sum of branch length times absolute difference of descendant shares / sum of branch length times sum of descendant shares. Reference root retained, root stem excluded.',
        ordination='Classical PCoA; Lingoes correction sqrt(d^2+2c) off diagonal, c=-minimum eigenvalue when negative. Percentages use corrected positive inertia. Raw distances remain available and used in summaries.',
        limitations='Same reference database, not verified identical processing. Sid weight units and preprocessing remain unverified. UFPF cross-sectional, five Sid profiles from one treated person. IBD bowel preparation confound. No inference, health score, or disease classification. Reference coverage denominators differ by source.',
        diagnosisConflict='UF-PF-2022-074-1N: Diagnosis Control vs Diagnosis2 PD; labeled Uncertain diagnosis and excluded from diagnosis-specific views.'))
    (outputdir/'sid-ufpf.json').write_text(json.dumps(output,separators=(',',':'),allow_nan=False)+'\n')
    with (outputdir/'sid-ufpf-profiles.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['sample','cohort','group','retained_fraction',*['SGB'+v for v in features]])
        for s,row in zip(samples,values): w.writerow([s['id'],s['cohort'],s['group'],s['coverage'],*row])
    (outputdir/'prepare_ufpf.py').write_text(Path(__file__).read_text())
    print(json.dumps(dict(features=len(features),coverageSid=coverage[:5].tolist(),coverageUF=[float(coverage[5:].min()),float(np.median(coverage[5:])),float(coverage[5:].max())],
        results={k:dict(n=len(s['indices']),axes=s['metrics']['unifrac']['axisFractions'],hulls=s['metrics']['unifrac']['hullsDisjoint'],summary=s['metrics']['unifrac']['summary']) for k,s in scenarios.items()}),indent=2))
    return output

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--archive',type=Path,default=ROOT.parent/'ufpf-audit/Source Data.xlsx')
    parser.add_argument('--sid',type=Path); parser.add_argument('--tree',type=Path); parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args(); build(args.archive,args.sid,args.tree,args.output_dir)
