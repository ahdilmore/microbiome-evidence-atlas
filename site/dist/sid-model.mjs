export function featuresAtLevel(data,level='genus'){
 if(level==='sgb')return data.features.map((f,j)=>({id:f.metaphlan_id,name:clean(f.species)||`SGB ${f.metaphlan_id}`,genus:clean(f.genus),species:clean(f.species),resolved:f.species.startsWith('s__')&&!f.species.includes('unclassified'),sgbs:[f.metaphlan_id],values:data.samples.map(s=>s.shares[j])}));
 const rank=level==='family'?'family':'genus',prefix=rank==='family'?'f__':'g__';
 const groups=new Map();data.features.forEach((f,j)=>{const label=f[rank]||'unassigned',resolved=label.startsWith(prefix)&&!label.includes('unclassified');const key=label;if(!groups.has(key))groups.set(key,{id:key,name:resolved?clean(key):`Unresolved · ${clean(key)}`,genus:rank==='genus'&&resolved?clean(key):'',family:rank==='family'&&resolved?clean(key):'',species:'',resolved,sgbs:[],values:data.samples.map(()=>0)});const g=groups.get(key);g.sgbs.push(f.metaphlan_id);g.values.forEach((_,i)=>g.values[i]+=data.samples[i].shares[j]);});return [...groups.values()];
}
export function clean(s){return (s||'').replace(/^LKT_/,'').replace(/^[a-z]__/,'').replaceAll('_',' ')}
export function changes(features,start,end){return features.map(f=>{const a=f.values[start],b=f.values[end];return{...f,start:a,end:b,delta:(b-a)*100,fold:a>0?b/a:null,status:a===0&&b>0?'Newly detected':a>0&&b===0?'Not detected at end':a===0&&b===0?'Not detected at either visit':a===b?'Unchanged':b>a?'Increased':'Decreased'}}).sort((a,b)=>Math.abs(b.delta)-Math.abs(a.delta)||a.name.localeCompare(b.name))}
export function evidenceMatch(card,feature,level){if(!feature.resolved||level==='family')return null;if(level==='sgb'&&card.species?.includes(feature.species))return 'Species name matches; strain/function unverified';if(card.genera.includes(feature.genus))return level==='genus'?'Genus context; species/strain may differ':'Genus context only; species match not established';return null}

// Thresholds are display filters on the original, unrenormalized shares.
export function filterChanges(rows,{abundance=0.1,minSamples=2,minDelta=0,query='',direction='all'}={}){
 return rows.map(r=>({...r,present:r.values.filter(v=>v>0&&v*100>=abundance).length})).filter(r=>r.present>=minSamples&&Math.abs(r.delta)>=minDelta&&(r.name.toLowerCase().includes(query.toLowerCase())||r.sgbs.some(id=>id.includes(query)))&&(direction==='all'||(direction==='up'?r.delta>0:r.delta<0)));
}
export function contextEvents(events,start,end,padding=0,category='Treatments'){
 const day=86400000,lo=Date.parse(start),hi=Date.parse(end);
 return events.filter(e=>category==='all'||e.category===category).filter(e=>Date.parse(e.date)<=hi+padding*day&&Date.parse(e.end_date||e.date)>=lo-padding*day).map(e=>({...e,relation:Date.parse(e.end_date||e.date)<lo?'Before window':Date.parse(e.date)>hi?'After window':e.end_date?'Overlaps window':'Within window'})).sort((a,b)=>a.date.localeCompare(b.date)||a.title.localeCompare(b.title));
}
