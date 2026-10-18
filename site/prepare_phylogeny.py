"""Generate normalized weighted UniFrac from the matching reference tree.
Requires Biopython 1.85 and numpy. No inferred tips, re-rooting, pruning or pseudocounts.
"""
from pathlib import Path
import json, hashlib, math
from Bio import Phylo
import numpy as np
ROOT=Path(__file__).resolve().parent

def weighted_matrix(tree, features, profiles):
    tips=tree.get_terminals(); names=[t.name for t in tips]
    assert len(names)==len(set(names)), 'Duplicate tips'
    assert set(features)<=set(names), 'Unmapped SGBs'
    values=np.asarray(profiles,dtype=float)
    assert np.isfinite(values).all() and (values>=0).all()
    assert np.allclose(values.sum(axis=1),1)
    mapping={name:i for i,name in enumerate(features)}
    masses={}; edges=[]
    for node in tree.find_clades(order='postorder'):
        if node.is_terminal(): mass=values[:,mapping[node.name]] if node.name in mapping else np.zeros(len(values))
        else: mass=sum((masses[c] for c in node.clades),np.zeros(len(values)))
        masses[node]=mass
        if node is tree.root: continue
        length=node.branch_length
        assert length is not None and math.isfinite(length) and length>=0, 'Invalid branch length'
        if np.any(mass): edges.append((length,mass))
    out=np.zeros((len(values),len(values)))
    for i in range(len(values)):
        for j in range(i+1,len(values)):
            numerator=sum(length*abs(m[i]-m[j]) for length,m in edges)
            denominator=sum(length*(m[i]+m[j]) for length,m in edges)
            out[i,j]=out[j,i]=numerator/denominator if denominator else 0
    return out,edges

if __name__=='__main__':
    path=ROOT/'source-data/sid/mpa_vJun23_CHOCOPhlAnSGB_202307.nwk'
    tree=Phylo.read(path,'newick');data=json.loads((ROOT/'dist/sid-data.json').read_text())
    features=[f['metaphlan_id'] for f in data['features']]
    matrix,edges=weighted_matrix(tree,features,[s['shares'] for s in data['samples']])
    # Independently reconstruct each sample's mass on branches by walking leaf-to-root paths.
    parents={c:p for p in tree.find_clades() for c in p.clades}; tips={t.name:t for t in tree.get_terminals()}
    independent={}
    for j,name in enumerate(features):
        node=tips[name]
        while node is not tree.root:
            independent.setdefault(node,np.zeros(len(data['samples'])))
            independent[node]+=np.array([s['shares'][j] for s in data['samples']]);node=parents[node]
    for i in range(5):
        for j in range(5):
            num=sum(n.branch_length*abs(m[i]-m[j]) for n,m in independent.items())
            den=sum(n.branch_length*(m[i]+m[j]) for n,m in independent.items())
            assert abs(matrix[i,j]-(num/den if den else 0))<1e-12
    result={'metric':'normalized weighted UniFrac','database':'mpa_vJun23_CHOCOPhlAnSGB_202307','treeSource':'https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202307.nwk','retrieved':'2026-10-17','treeSha256':hashlib.sha256(path.read_bytes()).hexdigest(),'dataSha256':hashlib.sha256((ROOT/'dist/sid-data.json').read_bytes()).hexdigest(),'referenceTips':len(tips),'mappedFeatures':len(features),'totalFeatures':len(features),'mappedAbundanceFraction':[1]*5,'sampleIds':[s['id'] for s in data['samples']],'matrix':matrix.tolist(),'method':'sum(branch length × absolute difference in descendant shares) / sum(branch length × sum of descendant shares). Reference root retained; root stem excluded. All supplied SGB shares included.','validation':'Biopython 1.85 parser; branch-sum and independent leaf-path calculations agree to 1e-12.'}
    (ROOT/'dist/sid-phylogeny.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({'mapped':len(features),'baseline':matrix[0].tolist()}))
