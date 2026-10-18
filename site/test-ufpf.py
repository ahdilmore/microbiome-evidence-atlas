"""Independent mathematical and data-integrity checks for the comparison."""
from io import StringIO
from pathlib import Path
import json
import numpy as np
import pandas as pd
from Bio import Phylo
from prepare_ufpf import weighted_distance, pcoa, hull, hulls_disjoint

ROOT=Path(__file__).resolve().parent
tree=Phylo.read(StringIO('(a:1,(b:1,c:1):1);'),'newick')
v=np.array([[1.,0,0],[0,1.,0],[0,0,1.],[0,.5,.5]])
d=weighted_distance(tree,['a','b','c'],v)
assert np.allclose(d,[[0,1,1,1],[1,0,.5,.25],[1,.5,0,.25],[1,.25,.25,0]])

# Example from scikit-bio's official weighted_unifrac documentation.
t=Phylo.read(StringIO('(((((U1:0.5,U2:0.5):0.5,U3:1.0):1.0):0.0,(U4:0.75,(U5:0.5,((U6:0.33,U7:0.62):0.5,U8:0.5):0.5):0.5):1.25):0.0)root;'),'newick')
v=np.array([[1,0,0,4,1,2,3,0],[0,1,1,6,0,1,0,0]],float);v/=v.sum(axis=1)[:,None]
assert round(weighted_distance(t,[f'U{i}' for i in range(1,9)],v)[0,1],2)==.33

data=json.loads((ROOT/'dist/sid-ufpf.json').read_text())
profiles=pd.read_csv(ROOT/'dist/sid-ufpf-profiles.tsv',sep='\t')
v=profiles.iloc[:,4:].values
assert np.isfinite(v).all() and (v>=0).all() and np.allclose(v.sum(axis=1),1)
assert len(data['samples'])==101 and len(data['features'])==1276
assert profiles['sample'].tolist()==[s['id'] for s in data['samples']]
assert profiles['group'].value_counts().to_dict()=={'PD':53,'IBD':26,'Control':16,'Sid':5,'Uncertain diagnosis':1}
assert not set(data['features']) & set(data['provenance']['numericIdsExcludedSymmetrically'])
assert any((v[:5].sum(axis=0)>0)&(v[5:].sum(axis=0)==0)), 'Sid-only taxa must remain'
assert any((v[:5].sum(axis=0)==0)&(v[5:].sum(axis=0)>0)), 'UF-only taxa must remain'
assert len(data['scenarios']['all-95']['indices'])==80
assert min(s['coverage'] for s in data['samples'][:5])>.997

# Independent branch reconstruction: walk each tip toward the supplied root.
t=Phylo.read(ROOT/'source-data/sid/mpa_vJun23_CHOCOPhlAnSGB_202307.nwk','newick')
parents={c:p for p in t.find_clades() for c in p.clades};tips={tip.name:tip for tip in t.get_terminals()}
selected=[0,1,2,3,4,5,30,60,90,100]; masses={}
for k,feature in enumerate(data['features']):
    node=tips[feature]
    while node is not t.root:
        masses.setdefault(node,np.zeros(len(selected)))
        masses[node]+=v[selected,k]
        node=parents[node]
actual=np.array(data['matrices']['unifrac'])
for i,ii in enumerate(selected):
    for j,jj in enumerate(selected):
        numerator=sum(node.branch_length*abs(m[i]-m[j]) for node,m in masses.items())
        denominator=sum(node.branch_length*(m[i]+m[j]) for node,m in masses.items())
        assert abs(actual[ii,jj]-numerator/denominator)<1e-12
for metric,matrix in data['matrices'].items():
    matrix=np.array(matrix)
    assert np.allclose(matrix,matrix.T) and np.allclose(np.diag(matrix),0)
    assert (matrix>=0).all() and (matrix<=1+1e-12).all()
    for scenario in data['scenarios'].values():
        idx=scenario['indices']; result=scenario['metrics'][metric]
        repeat=pcoa(matrix[np.ix_(idx,idx)]) # also asserts reconstruction in all dimensions
        assert np.allclose(repeat['coordinates'],result['coordinates'])
        for i,nearest in enumerate(result['summary']['nearestUFIndices']):
            distance=matrix[idx[i],idx[nearest]]
            assert abs(distance-matrix[idx[i],idx[5:]].min())<1e-12
        assert not result['summary']['strictDistanceGap']

square=np.array([[0,0],[1,0],[1,1],[0,1]])
assert not hulls_disjoint(hull(square),hull(square+.5))
assert hulls_disjoint(hull(square),hull(square+2))
assert not hulls_disjoint(hull(square),hull(square+[1,0]))
assert not data['scenarios']['all-0']['metrics']['unifrac']['hullsDisjoint']['0,1']
assert data['scenarios']['Control-0']['metrics']['unifrac']['hullsDisjoint']['0,1']
assert not data['scenarios']['Control-95']['metrics']['unifrac']['hullsDisjoint']['0,1']
print('PASS: published example, independent leaf-path distances, data integrity, 20 PCoA reconstructions, nearest profiles and hull geometry.')
