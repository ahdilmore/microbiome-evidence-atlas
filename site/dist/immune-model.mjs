export function quantile(values,p){
 const a=values.filter(Number.isFinite).sort((x,y)=>x-y);
 if(!a.length)return null;
 const n=(a.length-1)*p,i=Math.floor(n),f=n-i;
 return a[i]+(a[i+1]===undefined?0:(a[i+1]-a[i])*f);
}
export function parseProfiles(text){
 const [head,...lines]=text.trim().split(/\r?\n/),keys=head.split('\t');
 return lines.filter(Boolean).map(line=>{const cells=line.split('\t');return Object.fromEntries(keys.map((k,i)=>[k,cells[i]]));});
}
export const featureMatches=(candidate,features)=>candidate.rank==='genus'?features.filter(f=>f.id===candidate.taxon):features.filter(f=>String(f.name).replace(/ · SGB \d+$/,'')===candidate.taxon);
export function aggregate(row,matched){
 if(!row||!matched.length)return null;
 const members=[...new Set(matched.flatMap(f=>f.members||[f.id]))];
 let total=0;
 for(const id of members){const v=row[`SGB${id}`];if(v===undefined||v===null||v===''||!Number.isFinite(+v)||+v<0)return null;total+=+v;}
 return total;
}
export function immuneStories(candidates,profiles,features,sampleIds){
 const bySample=new Map(profiles.map(r=>[r.sample,r]));
 return candidates.map(c=>{
  const matched=featureMatches(c,features),sid=sampleIds.map(id=>aggregate(bySample.get(id),matched));
  const controls=profiles.filter(r=>r.group==='Control').map(r=>aggregate(r,matched)).filter(Number.isFinite);
  const available=sid.filter(Number.isFinite),positive=available.filter(v=>v>0).length;
  const availability=!available.length?'unavailable':positive?`Reported at ${positive}/${sampleIds.length} visits`:'Not detected in supplied profiles';
  return {...c,matched,sid,controls,availability,availableVisits:available.length,reference:{q25:quantile(controls,.25),median:quantile(controls,.5),q75:quantile(controls,.75),n:controls.length}};
 });
}
