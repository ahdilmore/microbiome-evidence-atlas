"""Descriptive per-group similarity and exact Bray-Curtis taxon contributions.

Uses the already audited, retained and normalized Sid/UFPF profiles. No additional
normalization, disease classifier, differential-abundance testing, or causal model.
Requirements: numpy, pandas, openpyxl. Run from the Site checkout after prepare_ufpf.py.
"""
from pathlib import Path
import hashlib, json, re
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent

def feature_stats(sid,ref):
    median=np.median(ref,axis=0)
    delta=sid-median
    # Each feature's exact contribution to the MEAN pairwise Bray-Curtis distance.
    contribution=np.abs(sid[:,None,:]-ref[None,:,:]).mean(axis=1)/2
    q10,q25,q75,q90=np.quantile(ref,[.1,.25,.75,.9],axis=0)
    return dict(median=median,mean=ref.mean(axis=0),q10=q10,q25=q25,q75=q75,q90=q90,
        delta=delta,contribution=contribution,above=(delta>1e-12).sum(axis=0),below=(delta < -1e-12).sum(axis=0),
        high=(sid>q90+1e-12).sum(axis=0),low=(sid<q10-1e-12).sum(axis=0))

def build():
    inputs=[ROOT/'dist/sid-ufpf.json',ROOT/'dist/sid-ufpf-profiles.tsv',ROOT/'dist/sid-data.json',ROOT.parent/'ufpf-audit/Source Data.xlsx']
    reference=json.loads(inputs[0].read_text()); profiles=pd.read_csv(inputs[1],sep='\t');sid=json.loads(inputs[2].read_text())
    workbook=pd.read_excel(inputs[3],sheet_name='UFPF_metaphlan_relab',index_col=0)
    values=profiles.iloc[:,4:].values;features=reference['features'];samples=reference['samples']
    assert profiles['sample'].tolist()==[s['id'] for s in samples]
    assert profiles.columns[4:].tolist()==['SGB'+f for f in features]
    assert np.allclose(values.sum(axis=1),1) and np.isfinite(values).all() and (values>=0).all()
    sidmap={f['metaphlan_id']:f for f in sid['features']}
    ufmap={re.search(r't__SGB(\d+)$',t).group(1):t.split('|') for t in workbook.index if re.search(r't__SGB(\d+)$',t)}
    sgb=[];groups={};conflicts=[]
    for j,f in enumerate(features):
        original=sidmap.get(f,{})
        if f in ufmap:
            genus=next((v for v in ufmap[f] if v.startswith('g__')),'')
            species=next((v for v in ufmap[f] if v.startswith('s__')),'')
            if original.get('genus') and original['genus']!=genus:
                conflicts.append(dict(sgb=f,sid=original['genus'],ufpf=genus))
                assert not original['genus'].startswith('g__'), 'Conflicting named genera require review'
        else: genus=original.get('genus','');species=original.get('species','')
        named=bool(genus.startswith('g__') and not re.search(r'(GGB\d|unclassified)',genus))
        groupkey=genus if genus.startswith('g__') and 'unclassified' not in genus else 'unresolved-SGB'+f
        groupname=genus[3:].replace('_',' ') if named else 'Unresolved '+(genus[3:] if groupkey==genus else 'SGB '+f)
        speciesnamed=species.startswith('s__') and not re.search(r'(GGB\d|unclassified)',species)
        speciesname=species[3:].replace('_',' ') if speciesnamed else 'Unresolved SGB '+f
        sgb.append(dict(id=f,name=speciesname+' · SGB '+f if speciesnamed else speciesname,resolved=bool(speciesnamed),members=[f]))
        groups.setdefault(groupkey,dict(id=groupkey,name=groupname,resolved=named,members=[],columns=[]))
        groups[groupkey]['members'].append(f);groups[groupkey]['columns'].append(j)
    genusvalues=np.stack([values[:,g['columns']].sum(axis=1) for g in groups.values()],axis=1)
    genera=[{k:v for k,v in g.items() if k!='columns'} for g in groups.values()]
    assert np.allclose(genusvalues.sum(axis=1),1)
    levels={'genus':(genera,genusvalues),'sgb':(sgb,values)}
    output=dict(generated='2026-10-17',samples=samples[:5],groups=['Control','PD','IBD'],levels={},scenarios={},
        provenance=dict(inputSHA256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        taxonomy='Exact SGB IDs use UFPF taxonomy when available, otherwise Sid taxonomy. Named shared genera agree. UF GGB placeholders remain unresolved; unresolved Sid-only SGBs remain separate.',
        labelDifferences=conflicts,contribution='For each Sid visit and taxon: mean over reference profiles of abs(Sid share - reference share) / 2. Sum over all taxa equals mean pairwise Bray-Curtis at that rank. Across-visit ranking averages the five visits equally. These are descriptive contributions to abundance distance, not contributions to weighted UniFrac or causes of biological change.',
        direction='Sid share minus reference-group median share at each visit. Direction uses a 1e-12 tolerance. Reference 10th–90th percentiles describe between-person spread, not a confidence interval. No pseudocounts or p-values.',
        similarity='Median pairwise distance from each Sid profile to every eligible reference profile in a group. Lower means more similar for this metric. Reference baseline: median distance from each reference profile to all other profiles in its own group, excluding self; show median and IQR of these per-profile medians.',
        scope='Same retained SGB profiles and coverage screen as the existing comparison. One diagnosis-conflicted profile remains only in the all-UFPF driver view, never in named diagnosis groups. IBD samples followed initiation of bowel preparation. Repeated Sid visits are not independent people. Study effects, Sid input-unit semantics and full processing equivalence remain unresolved.'))
    for level,(taxa,v) in levels.items():
        output['levels'][level]=dict(features=taxa,sidValues=v[:5].T.tolist())
    for floor in (0,95):
        scenario=dict(similarity={},drivers={})
        for metric,matrix in reference['matrices'].items():
            d=np.asarray(matrix);scenario['similarity'][metric]={}
            for group in output['groups']:
                idx=[i for i,s in enumerate(samples) if s['group']==group and s['coverage']>=floor/100]
                cross=d[:5,idx];within=d[np.ix_(idx,idx)].copy();np.fill_diagonal(within,np.nan)
                baseline=np.nanmedian(within,axis=1)
                scenario['similarity'][metric][group]=dict(n=len(idx),median=np.median(cross,axis=1).tolist(),
                    q25=np.quantile(cross,.25,axis=1).tolist(),q75=np.quantile(cross,.75,axis=1).tolist(),
                    typicalMedian=float(np.median(baseline)),typicalQ25=float(np.quantile(baseline,.25)),typicalQ75=float(np.quantile(baseline,.75)))
        for group in ('all','Control','PD','IBD'):
            idx=[i for i,s in enumerate(samples) if i>=5 and (group=='all' or s['group']==group) and s['coverage']>=floor/100]
            result=dict(n=len(idx),sampleIds=[samples[i]['id'] for i in idx],levels={})
            for level,(taxa,v) in levels.items():
                stats=feature_stats(v[:5],v[idx]);contribution=stats['contribution']
                direct=np.array([np.abs(s-v[idx]).sum(axis=1).mean()/2 for s in v[:5]])
                assert np.allclose(contribution.sum(axis=1),direct,atol=1e-12)
                if level=='sgb':assert np.allclose(direct,np.array(reference['matrices']['braycurtis'])[:5,idx].mean(axis=1),atol=1e-12)
                result['levels'][level]=dict(totalMeanDistance=direct.tolist(),rows=[{k:(a[:,j].tolist() if a.ndim==2 else float(a[j])) for k,a in stats.items()} for j in range(len(taxa))])
            scenario['drivers'][group]=result
        output['scenarios'][str(floor)]=scenario
    (ROOT/'dist/sid-ufpf-drivers.json').write_text(json.dumps(output,separators=(',',':'),allow_nan=False)+'\n')
    (ROOT/'dist/prepare_ufpf_drivers.py').write_text(Path(__file__).read_text())
    for group in ('all','Control','PD','IBD'):
        rows=output['scenarios']['0']['drivers'][group]['levels']['genus']['rows'];taxa=output['levels']['genus']['features']
        rank=sorted(range(len(rows)),key=lambda i:np.mean(rows[i]['contribution']),reverse=True)
        print(group,[(taxa[i]['name'],round(np.mean(rows[i]['contribution'])*100,2),int(rows[i]['above']),int(rows[i]['below'])) for i in rank if taxa[i]['resolved']][:10])
    print('Generated',len(genera),'genus bins;',len(sgb),'SGBs;',len(conflicts),'unresolved label differences; named shared genera agree.')
    return output

if __name__=='__main__':build()
