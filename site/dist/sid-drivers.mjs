const $=id=>document.getElementById(id);
const groups=['Control','PD','IBD'];
const labels={all:'All UFPF profiles',Control:'Controls',PD:'Parkinson’s disease',IBD:'IBD · bowel prep'};
const colors={Control:'#416eb0',PD:'#25786c',IBD:'#bb7b26'};
const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const mean=v=>v.reduce((a,b)=>a+b,0)/v.length;
const fixed=v=>v.toFixed(3),pct=v=>(100*v).toFixed(2)+'%';
const pp=v=>(v>1e-12?'+':'')+(100*v).toFixed(2);
let data,selected=null,currentRows=[];
const shortDate=(s,i)=>new Intl.DateTimeFormat('en-US',{month:'short',day:'numeric',timeZone:'UTC'}).format(new Date(s.date))+(i===2?'*':i>2?'†':'');
function download(name,rows){const csv=rows.map(r=>r.map(v=>'"'+String(v??'').replaceAll('"','""')+'"').join(',')).join('\r\n');const url=URL.createObjectURL(new Blob([csv],{type:'text/csv'}));const a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
function state(){const coverage=$('groupCoverage').value,group=$('driverGroup').value,level=$('driverLevel').value;return{coverage,group,level,scenario:data.scenarios[coverage],taxa:data.levels[level].features,values:data.levels[level].sidValues,result:data.scenarios[coverage].drivers[group].levels[level]};}

function chartBase(title,maxY,percent=false){
 const W=900,H=330,L=72,R=846,T=45,B=257;
 const times=data.samples.map(s=>Date.parse(s.date)),lo=times[0],span=times.at(-1)-lo;
 const x=i=>L+(R-L)*(times[i]-lo)/span,y=v=>B-(B-T)*v/maxY;
 let svg=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(title)}" style="font-family:Arial,sans-serif"><title>${esc(title)}</title>`;
 for(let tick=0;tick<=4;tick++){const v=maxY*tick/4;svg+=`<path d="M${L},${y(v)}H${R}" stroke="#dce3db"/><text x="${L-12}" y="${y(v)+4}" text-anchor="end" font-size="12" fill="#677b77">${percent?(v*100).toFixed(1)+'%':v.toFixed(2)}</text>`;}
 for(let i=0;i<5;i++)svg+=`<text transform="translate(${x(i)},${B+25}) rotate(-20)" text-anchor="end" font-size="12" fill="#677b77">${shortDate(data.samples[i],i)}</text>`;
 svg+=`<text x="${L}" y="21" font-size="13" fill="#193a38">${esc(title)}</text>`;
 return {svg,x,y,L,R};
}

function renderGroups(){
 const metric=$('groupMetric').value,coverage=$('groupCoverage').value;
 const stats=data.scenarios[coverage].similarity[metric];
 const closest=data.samples.map((_,i)=>[...groups].sort((a,b)=>stats[a].median[i]-stats[b].median[i])[0]);
 const counts=groups.map(g=>[g,closest.filter(x=>x===g).length]).filter(([,n])=>n);
 const metricName=metric==='unifrac'?'weighted UniFrac':'SGB Bray–Curtis';
 $('groupFinding').textContent=`By median ${metricName} distance, the closest group is ${counts.map(([g,n])=>`${labels[g]} at ${n} of 5 visits`).join('; ')}. This ranks these reference groups only; it does not assign Sid to a disease category. Small differences and overlapping distance ranges should not be treated as distinct biological states.`;
 const chart=chartBase(`Median ${metricName} distance to each group · lower = more similar`,1);
 for(const g of groups){const s=stats[g];chart.svg+=`<polyline points="${s.median.map((v,i)=>`${chart.x(i)},${chart.y(v)}`).join(' ')}" fill="none" stroke="${colors[g]}" stroke-width="2.5"/>`;for(let i=0;i<5;i++)chart.svg+=`<circle cx="${chart.x(i)}" cy="${chart.y(s.median[i])}" r="5" fill="${colors[g]}" stroke="white"><title>${esc(labels[g])} · ${data.samples[i].date}: ${fixed(s.median[i])}; reference-distance IQR ${fixed(s.q25[i])}–${fixed(s.q75[i])}</title></circle>`;}
 $('groupPlot').innerHTML=chart.svg+'</svg>';
 $('groupTable').innerHTML=`<table><thead><tr><th>Sid visit</th>${groups.map(g=>`<th>${labels[g]}<small>n = ${stats[g].n}</small></th>`).join('')}<th>Smallest median</th></tr></thead><tbody>${data.samples.map((s,i)=>`<tr><td>${s.date}${i===2?'*':i>2?'†':''}</td>${groups.map(g=>`<td class="${closest[i]===g?'closest':''}">${fixed(stats[g].median[i])}<small>IQR ${fixed(stats[g].q25[i])}–${fixed(stats[g].q75[i])}</small></td>`).join('')}<td>${labels[closest[i]]}</td></tr>`).join('')}</tbody></table>`;
 $('groupBaseline').innerHTML=`<table><thead><tr><th>Reference group</th><th>Profiles</th><th>Typical within-group distance</th><th>Middle 50% of per-profile medians</th></tr></thead><tbody>${groups.map(g=>`<tr><td>${labels[g]}</td><td>${stats[g].n}</td><td>${fixed(stats[g].typicalMedian)}</td><td>${fixed(stats[g].typicalQ25)}–${fixed(stats[g].typicalQ75)}</td></tr>`).join('')}</tbody></table>`;
 $('groupExport').disabled=false;
}

function direction(r){return r.above===5?'Higher 5/5':r.below===5?'Lower 5/5':r.above===0&&r.below===0?'At median 5/5':`${r.above} higher · ${r.below} lower`;}
function eligibleRows(){
 const {taxa,result}=state(),query=$('driverSearch').value.toLowerCase().trim(),filter=$('driverDirection').value;
 return taxa.map((f,i)=>({...f,index:i,...result.rows[i],score:mean(result.rows[i].contribution)})).filter(r=>
  ($('driverUnresolved').checked||r.resolved)&&(!query||r.name.toLowerCase().includes(query)||r.members.some(id=>id.includes(query)))&&
  (filter==='all'||filter==='above'&&r.above===5||filter==='below'&&r.below===5||filter==='outside'&&(r.high===5||r.low===5)))
  .sort((a,b)=>b.score-a.score||a.name.localeCompare(b.name));
}
function renderDrivers(){
 const {coverage,group,level,result,taxa}=state(),total=mean(result.totalMeanDistance),n=data.scenarios[coverage].drivers[group].n;
 currentRows=eligibleRows();const shown=currentRows.slice(0,20);
 if(!shown.some(r=>r.id===selected))selected=shown[0]?.id??null;
 $('driverStatus').textContent=`${labels[group]} · ${n} reference profiles · showing ${shown.length} of ${currentRows.length} matching ${level==='genus'?'genus groups':'SGBs'}, ranked across all five visits. Mean ${level}-level Bray–Curtis distance: ${fixed(total)}. Hidden features remain in the distance denominator (${taxa.length} total features).`;
 $('driverTable').innerHTML=`<table><thead><tr><th>Microbe</th><th>Share of mean distance</th>${data.samples.map((s,i)=>`<th>${shortDate(s,i)} Δ pp</th>`).join('')}<th>Vs reference median</th></tr></thead><tbody>${shown.map(r=>`<tr class="${r.id===selected?'selected':''}"><td><button data-driver-id="${esc(r.id)}" aria-pressed="${r.id===selected}">${esc(r.name)}</button>${!r.resolved?'<small> · unresolved</small>':''}</td><td>${(100*r.score/total).toFixed(1)}%</td>${r.delta.map(v=>`<td style="color:${v>1e-12?'#a84349':v < -1e-12?'#25786c':'#677b77'};background:${v>1e-12?'#fbf0ed':v < -1e-12?'#edf5ef':'transparent'}">${Math.abs(v)<.00005&&Math.abs(v)>1e-12?(v>0?'+<0.01':'−<0.01'):pp(v)}</td>`).join('')}<td>${direction(r)}</td></tr>`).join('')}</tbody></table>`;
 if(!shown.length)$('driverTable').innerHTML='<p class="empty">No features match these filters. Try another pattern or include unresolved bins.</p>';
 $('driverExport').disabled=!currentRows.length;
 renderDetail();
}
function renderDetail(){
 const {values,group}=state(),r=currentRows.find(r=>r.id===selected);
 $('driverDetail').hidden=!r;if(!r)return;
 const v=values[r.index];$('driverTitle').textContent=r.name;
 $('driverPattern').textContent=`${direction(r)} relative to the ${labels[group]} median. Above the reference 90th percentile at ${r.high}/5 visits; below its 10th percentile at ${r.low}/5. Consistency describes the direction of a relative-abundance difference, not a stable absolute amount or treatment effect.`;
 const chart=chartBase(`${r.name} · Sid share vs ${labels[group]}`,Math.max(...v,r.q90,.001)*1.15,true);
 chart.svg+=`<rect x="${chart.L}" y="${chart.y(r.q90)}" width="${chart.R-chart.L}" height="${chart.y(r.q10)-chart.y(r.q90)}" fill="#416eb0" fill-opacity=".12"/><path d="M${chart.L},${chart.y(r.median)}H${chart.R}" stroke="#416eb0" stroke-width="2" stroke-dasharray="6 4"/><polyline points="${v.map((s,i)=>`${chart.x(i)},${chart.y(s)}`).join(' ')}" fill="none" stroke="#a84349" stroke-width="2.5"/>`;
 for(let i=0;i<5;i++)chart.svg+=`<circle cx="${chart.x(i)}" cy="${chart.y(v[i])}" r="5" fill="#a84349" stroke="white"><title>${data.samples[i].date}: ${pct(v[i])}; reference median ${pct(r.median)}</title></circle>`;
 $('driverPlot').innerHTML=chart.svg+'</svg>';
 $('driverValues').innerHTML=`<table><thead><tr><th>Visit</th><th>Sid share</th><th>Reference median</th><th>Reference 10th–90th</th><th>Δ pp</th></tr></thead><tbody>${data.samples.map((s,i)=>`<tr><td>${s.date}</td><td>${pct(v[i])}</td><td>${pct(r.median)}</td><td>${pct(r.q10)}–${pct(r.q90)}</td><td>${pp(r.delta[i])}</td></tr>`).join('')}</tbody></table>`;
}

for(const id of ['groupMetric','groupCoverage'])$(id).addEventListener('change',()=>{if(data){renderGroups();renderDrivers();}});
for(const id of ['driverGroup','driverLevel','driverDirection','driverUnresolved'])$(id).addEventListener('change',()=>{if(data)renderDrivers();});
$('driverSearch').addEventListener('input',()=>{if(data)renderDrivers();});
$('driverTable').addEventListener('click',e=>{const button=e.target.closest('[data-driver-id]');if(!button)return;selected=button.dataset.driverId;renderDrivers();});
$('groupExport').addEventListener('click',()=>{const metric=$('groupMetric').value,coverage=$('groupCoverage').value,s=data.scenarios[coverage].similarity[metric];download(`sid-uf-groups-${metric}-coverage-${coverage}.csv`,[['sid_sample','reported_date','date_note','metric','minimum_reference_coverage_percent','group','n','median_distance','q25_distance','q75_distance','typical_within_group_median','within_group_q25','within_group_q75'],...data.samples.flatMap((sample,i)=>groups.map(g=>[sample.id,sample.date,sample.dateNote,metric,coverage,g,s[g].n,s[g].median[i],s[g].q25[i],s[g].q75[i],s[g].typicalMedian,s[g].typicalQ25,s[g].typicalQ75]))]);});
$('driverExport').addEventListener('click',()=>{const {coverage,group,level,result,values}=state(),total=mean(result.totalMeanDistance);download(`sid-uf-microbial-differences-${group}-${level}-coverage-${coverage}.csv`,[['taxon_id','taxon','resolved','reference_group','rank','minimum_reference_coverage_percent','pattern_filter','search_filter','include_unresolved','mean_Bray_Curtis_contribution','fraction_of_total_mean_distance','reference_median_share','reference_q10_share','reference_q90_share','visits_above_median','visits_below_median','visits_above_q90','visits_below_q10',...data.samples.flatMap(s=>[s.id+'_reported_date',s.id+'_date_note',s.id+'_share',s.id+'_delta_pp',s.id+'_BC_contribution'])],...currentRows.map(r=>[r.id,r.name,r.resolved,group,level,coverage,$('driverDirection').value,$('driverSearch').value,$('driverUnresolved').checked,r.score,r.score/total,r.median,r.q10,r.q90,r.above,r.below,r.high,r.low,...data.samples.flatMap((s,i)=>[s.date,s.dateNote,values[r.index][i],r.delta[i]*100,r.contribution[i]])])]);});

try{const response=await fetch('sid-ufpf-drivers.json');if(!response.ok)throw Error(`HTTP ${response.status}`);data=await response.json();renderGroups();renderDrivers();}catch(error){$('groupFinding').textContent='Group and microbe comparisons could not load. Please reload; the other timeline analyses remain available.';console.error(error);}
