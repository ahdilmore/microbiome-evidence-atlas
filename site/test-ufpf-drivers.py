from pathlib import Path
import json
import numpy as np
from prepare_ufpf_drivers import feature_stats

ROOT=Path(__file__).resolve().parent
# Distinguish mean pairwise absolute difference from distance to a median.
ref=np.array([[0.,1.],[0.,1.],[.9,.1]])
sid=np.array([[0.,1.],[.6,.4]])
s=feature_stats(sid,ref)
assert np.allclose(s['median'],[0,1])
assert np.allclose(s['contribution'],[[.15,.15],[.25,.25]])
assert np.allclose(s['delta'],[[0,0],[.6,-.6]])
assert s['above'].tolist()==[1,0] and s['below'].tolist()==[0,1]

d=json.loads((ROOT/'dist/sid-ufpf-drivers.json').read_text())
original=json.loads((ROOT/'dist/sid-ufpf.json').read_text())
for floor in ['0','95']:
    scenario=d['scenarios'][floor]
    for metric,groupstats in scenario['similarity'].items():
        matrix=np.array(original['matrices'][metric])
        for group,stats in groupstats.items():
            idx=[i for i,s in enumerate(original['samples']) if s['group']==group and s['coverage']>=int(floor)/100]
            assert len(idx)==stats['n']
            assert all(original['samples'][i]['group']!='Uncertain diagnosis' for i in idx)
            assert np.allclose(stats['median'],np.median(matrix[:5,idx],axis=1))
            # Independently reconstruct leave-one-out medians, without NaN masking.
            baseline=[np.median([matrix[i,j] for j in idx if i!=j]) for i in idx]
            assert abs(stats['typicalMedian']-np.median(baseline))<1e-12
    for group,result in scenario['drivers'].items():
        for level,x in result['levels'].items():
            assert len(x['rows'])==len(d['levels'][level]['features'])
            contributions=np.array([r['contribution'] for r in x['rows']])
            assert np.allclose(contributions.sum(axis=0),x['totalMeanDistance'],atol=1e-12)
            assert (contributions>=0).all()
            for i,r in enumerate(x['rows']):
                v=np.array(d['levels'][level]['sidValues'][i])
                assert np.allclose(r['delta'],v-r['median'])
                assert r['above']==np.sum(v-r['median']>1e-12)
                assert r['below']==np.sum(v-r['median'] < -1e-12)
                assert r['high']==np.sum(v-r['q90']>1e-12)
                assert r['low']==np.sum(v-r['q10'] < -1e-12)
                assert r['q10']<=r['q25']<=r['median']<=r['q75']<=r['q90']
                assert r['above']+r['below']<=5
    stats=scenario['similarity']['unifrac']
    closest=[min(d['groups'],key=lambda g:stats[g]['median'][i]) for i in range(5)]
    assert closest==['PD','Control','Control','Control','Control']
    taxa=d['levels']['genus']['features'];lookup={t['name']:i for i,t in enumerate(taxa)}
    controls=scenario['drivers']['Control']['levels']['genus']['rows']
    for name in ['Blautia','Faecalibacterium']:assert controls[lookup[name]]['below']==5
    assert controls[lookup['Alistipes']]['above']==5
print('PASS: mean-pairwise contribution fixture; 16 decompositions; visit-direction and reference quantile checks; per-group medians and leave-one-out baselines; sensitivity findings.')
