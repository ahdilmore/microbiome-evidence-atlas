const $=id=>document.getElementById(id);
const colors={Sid:'#a84349',Control:'#416eb0',PD:'#25786c',IBD:'#bb7b26','Uncertain diagnosis':'#737781'};
const labels={Sid:'Sid · repeated visits',Control:'UFPF · control',PD:'UFPF · Parkinson’s disease',IBD:'UFPF · IBD / bowel prep','Uncertain diagnosis':'UFPF · uncertain diagnosis'};
const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const percent=v=>(100*v).toFixed(1)+'%';
const decimal=v=>v.toFixed(3);
let data,current;

function download(name,body,type){const url=URL.createObjectURL(new Blob([body],{type}));const a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
function csv(rows){return rows.map(row=>row.map(v=>'"'+String(v??'').replaceAll('"','""')+'"').join(',')).join('\r\n');}
function hull(points){const p=[...points].sort((a,b)=>a[0]-b[0]||a[1]-b[1]);const cross=(o,a,b)=>(a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0]);const half=list=>{const h=[];for(const v of list){while(h.length>=2&&cross(h.at(-2),h.at(-1),v)<=0)h.pop();h.push(v);}return h.slice(0,-1);};return [...half(p),...half([...p].reverse())];}
function description(s,i){return s.cohort==='Sid'?`Sid visit ${i+1} · ${s.date} · ${s.dateNote} · ${percent(s.coverage)} of supplied SGB abundance retained`:`${s.id} · ${labels[s.group]} · age ${s.age} · ${s.reads.toLocaleString()} reported reads · ${percent(s.coverage)} of classified abundance retained`;}

function render(){
  const group=$('ufGroup').value,coverage=$('ufCoverage').value,metric=$('ufMetric').value,axis=Number($('ufAxes').value);
  const scenario=data.scenarios[`${group}-${coverage}`],result=scenario.metrics[metric];
  const samples=scenario.indices.map(i=>data.samples[i]);
  const metricName=metric==='unifrac'?'Normalized weighted UniFrac':'Bray–Curtis';
  current={group,coverage,metric,axis,scenario,result,samples,metricName};
  const refs=samples.slice(5),counts={};for(const s of samples)counts[s.group]=(counts[s.group]||0)+1;
  $('ufStatus').textContent=`5 Sid visits + ${refs.length} UFPF profiles · ${data.features.length.toLocaleString()} SGB tips · ordination recalculated for this reference subset. ${coverage==='95'?'All five Sid profiles also retain ≥ 95%.':''}`;
  const disjoint=result.hullsDisjoint[`0,${axis}`];
  $('ufCards').innerHTML=`<div><span>Projected cohort envelopes</span><b>${disjoint?'Separate in this view':'Overlap in this view'}</b><small>Convex hulls on axes 1 & ${axis+1}; not a significance test.</small></div><div><span>Displayed axes explain</span><b>${percent(result.axisFractions[0]+result.axisFractions[axis])}</b><small>Of corrected positive inertia. Other axes contain the rest.</small></div><div><span>Reference abundance retained</span><b>${percent(Math.min(...refs.map(s=>s.coverage)))}–${percent(Math.max(...refs.map(s=>s.coverage)))}</b><small>Before renormalization; source-specific denominators.</small></div>`;
  $('ufLegend').innerHTML=Object.entries(counts).map(([g,n])=>`<span><i style="background:${colors[g]}" class="uf-swatch ${g==='Sid'?'uf-diamond':''}"></i>${esc(labels[g])} (${n})</span>`).join('');
  plot(samples,result,axis,metricName);
  $('ufPoint').textContent='Hover, focus, or tap a point for sample details. Sid’s numbered visits are connected in chronological order.';
  $('ufProjectionNote').textContent=`${disjoint?'The two cohort envelopes do not overlap in this projection. This does not establish complete separation in the full distance matrix.':'Sid’s and UFPF’s convex hulls overlap in this projection; this view does not show complete separation.'} ${result.lingoesConstant>0?`Lingoes correction applied (c = ${result.lingoesConstant.toFixed(5)}); negative inertia before correction: ${percent(result.negativeInertiaFraction)}.`:'No negative-eigenvalue correction needed.'} Changing the reference subset rotates/recomputes the axes; plot positions cannot be compared directly across subsets.`;
  const summary=result.summary;
  $('ufDistanceSummary').innerHTML=`<table><thead><tr><th>Original ${metricName} distances</th><th>Minimum</th><th>Median</th><th>Maximum</th></tr></thead><tbody>${[['Sid ↔ UFPF',summary.between],['Within Sid',summary.withinSid],['Within UFPF',summary.withinUF]].map(([label,s])=>`<tr><td>${label}</td><td>${decimal(s.min)}</td><td>${decimal(s.median)}</td><td>${decimal(s.max)}</td></tr>`).join('')}</tbody></table><p class="small">Strict distance-gap check: ${summary.strictDistanceGap?'every between-cohort distance exceeds every within-cohort distance.':'not met; between-cohort and within-cohort distance ranges overlap.'} This is a descriptive check, not a classifier or formal hypothesis test.</p>`;
  $('ufNearest').innerHTML=`<table><thead><tr><th>Sid visit</th><th>Nearest UFPF profile</th><th>Source group</th><th>Distance</th><th>Reference retained</th></tr></thead><tbody>${samples.slice(0,5).map((s,i)=>{const r=samples[summary.nearestUFIndices[i]];return `<tr><td>${i+1} · ${s.date}</td><td>${r.id}</td><td>${esc(labels[r.group])}</td><td>${decimal(summary.nearestUFDistances[i])}</td><td>${percent(r.coverage)}</td></tr>`;}).join('')}</tbody></table>`;
  $('ufCoverageNote').textContent=`Across all input profiles, Sid retains ${percent(Math.min(...data.samples.slice(0,5).map(s=>s.coverage)))}–100% of supplied SGB weight. UFPF retains ${percent(Math.min(...data.samples.slice(5).map(s=>s.coverage)))}–100% of classified leaf abundance. ${data.provenance.groupedRowsExcluded.length} grouped SGB rows and ${data.provenance.nonSGBrankRowsExcluded.length} non-SGB leaf rows are excluded. ${data.provenance.unmappedNumericIds.length} numeric tip IDs are unmapped. Coverage is not the proportion of sequencing reads classified; the sources have different original denominators. The ≥ 95% screen excludes 21 of 96 UFPF profiles before diagnosis filtering. Excluded abundance can materially change composition after renormalization.`;
  for(const id of ['ufSvg','ufCsv','ufDistances'])$(id).disabled=false;
}

function plot(samples,result,axis,metricName){
  const coords=result.coordinates.map(c=>[c[0],c[axis]]);const W=920,H=560,left=92,top=75,pw=766,ph=390;
  const xs=coords.map(c=>c[0]),ys=coords.map(c=>c[1]);
  const cx=(Math.max(...xs)+Math.min(...xs))/2,cy=(Math.max(...ys)+Math.min(...ys))/2;
  const unit=Math.max((Math.max(...xs)-Math.min(...xs))/pw,(Math.max(...ys)-Math.min(...ys))/ph)*1.22||1;
  const x=v=>left+pw/2+(v-cx)/unit,y=v=>top+ph/2-(v-cy)/unit;
  const pixels=coords.map(([a,b])=>[x(a),y(b)]);
  let svg=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" role="img" aria-labelledby="ufSvgTitle ufSvgDesc" style="font-family:Arial,sans-serif;background:#fffefa"><title id="ufSvgTitle">${metricName} PCoA, Sid and UFPF</title><desc id="ufSvgDesc">${samples.length-5} UFPF profiles; Sid’s five repeated profiles are numbered diamonds joined chronologically. The two axes use equal distance scales. ${result.hullsDisjoint[`0,${axis}`]?'Cohort envelopes separate':'Cohort envelopes overlap'} in this projection.</desc><rect width="920" height="560" fill="#fffefa"/><text x="${left}" y="25" fill="#29392f" font-size="16" font-weight="bold">${metricName} · PCoA</text><text x="${left}" y="46" fill="#657267" font-size="12">${samples.length-5} UFPF profiles · 5 Sid visits · ${current.coverage==='95'?'≥95% reference abundance retained':'all reference coverage levels'} · Lingoes-corrected</text>`;
  for(let tick=0;tick<=4;tick++){
    const a=left+pw*tick/4,b=top+ph*tick/4;
    svg+=`<path d="M${a},${top}V${top+ph} M${left},${b}H${left+pw}" stroke="#e6e8df" stroke-width="1"/><text x="${a}" y="${top+ph+24}" text-anchor="middle" font-size="12" fill="#657267">${((a-left-pw/2)*unit+cx).toFixed(2)}</text><text x="${left-12}" y="${b+4}" text-anchor="end" font-size="12" fill="#657267">${(cy-(b-top-ph/2)*unit).toFixed(2)}</text>`;
  }
  if($('ufHulls').checked)for(const [p,color] of [[pixels.slice(5),'#3e766c'],[pixels.slice(0,5),colors.Sid]])svg+=`<polygon points="${hull(p).map(c=>c.join(',')).join(' ')}" fill="${color}" fill-opacity=".06" stroke="${color}" stroke-opacity=".6" stroke-dasharray="5 4"/>`;
  svg+=`<polyline points="${pixels.slice(0,5).map(c=>c.join(',')).join(' ')}" fill="none" stroke="${colors.Sid}" stroke-width="2" stroke-opacity=".7"/>`;
  for(const i of [...Array(samples.length-5)].map((_,i)=>i+5).concat([0,1,2,3,4])){
    const s=samples[i],[px,py]=pixels[i],label=esc(description(s,i));
    svg+=`<g data-uf-index="${i}" tabindex="0" role="button" aria-label="${label}" style="cursor:pointer"><title>${label}</title><circle cx="${px}" cy="${py}" r="14" fill="transparent" pointer-events="all"/>`;
    svg+=i<5?`<path d="M${px},${py-10}l10,10 -10,10 -10,-10Z" fill="${colors.Sid}" stroke="white" stroke-width="1.5"/><text x="${px+13}" y="${py-10}" font-size="13" font-weight="bold" fill="${colors.Sid}" stroke="#fffefa" stroke-width="3" paint-order="stroke">${i+1}</text>`:`<circle cx="${px}" cy="${py}" r="5" fill="${colors[s.group]}" fill-opacity=".75" stroke="white" stroke-width=".7"/>`;
    svg+='</g>';
  }
  svg+=`<text x="${left+pw/2}" y="${H-28}" text-anchor="middle" font-size="14" fill="#29392f">PCoA 1 · ${percent(result.axisFractions[0])}</text><text transform="translate(25,${top+ph/2}) rotate(-90)" text-anchor="middle" font-size="14" fill="#29392f">PCoA ${axis+1} · ${percent(result.axisFractions[axis])}</text></svg>`;
  $('ufPlot').innerHTML=svg;
  for(const point of $('ufPlot').querySelectorAll('[data-uf-index]')){
    const show=()=>{const i=Number(point.dataset.ufIndex);$('ufPoint').textContent=description(samples[i],i);};
    for(const event of ['pointerenter','focus','click'])point.addEventListener(event,show);
    point.addEventListener('keydown',event=>{if(event.key==='Enter'||event.key===' '){event.preventDefault();show();}});
  }
}

for(const id of ['ufGroup','ufMetric','ufCoverage','ufAxes','ufHulls'])$(id).addEventListener('change',()=>{if(data)render();});
$('ufCsv').addEventListener('click',()=>{const {samples,result,metric,group,coverage}=current;download(`sid-ufpf-${metric}-${group}-${coverage}-coordinates.csv`,csv([['sample','cohort','group','reported_date','date_note','retained_fraction','PCoA1','PCoA2','PCoA3','axis1_fraction','axis2_fraction','axis3_fraction','lingoes_constant','metric','reference_group','minimum_reference_coverage_percent'],...samples.map((s,i)=>[s.id,s.cohort,s.group,s.date,s.dateNote,s.coverage,...result.coordinates[i],...result.axisFractions,result.lingoesConstant,metric,group,coverage])]),'text/csv');});
$('ufDistances').addEventListener('click',()=>{const {samples,scenario,metric,group,coverage}=current;download(`sid-ufpf-${metric}-${group}-${coverage}-raw-distances.csv`,csv([['sample',...samples.map(s=>s.id)],...scenario.indices.map((i,j)=>[samples[j].id,...scenario.indices.map(k=>data.matrices[metric][i][k])])]),'text/csv');});
$('ufSvg').addEventListener('click',()=>{const svg=$('ufPlot').querySelector('svg').cloneNode(true);svg.setAttribute('viewBox','0 0 920 680');svg.setAttribute('height','680');svg.setAttribute('width','920');const footer=document.createElementNS('http://www.w3.org/2000/svg','g');footer.innerHTML=`<rect x="0" y="560" width="920" height="120" fill="#fffefa"/>`+Object.entries(colors).map(([g,color],i)=>`<text x="92" y="${578+i*17}" font-size="12" fill="${color}">${esc(labels[g])}</text>`).join('')+'<text x="380" y="595" font-size="12" fill="#657267">Exploratory comparison; same database, processing equivalence unverified.</text><text x="380" y="615" font-size="12" fill="#657267">Sid: 1 Oct 22; 2 Dec 2; 3 Feb 13*; 4 Mar 10†; 5 Mar 27†.</text><text x="380" y="635" font-size="12" fill="#657267">* Conflicting February date. † Receipt dates. 2025–2026.</text>';svg.append(footer);download(`sid-ufpf-${current.metric}-${current.group}-${current.coverage}.svg`,new XMLSerializer().serializeToString(svg),'image/svg+xml');});

try{const response=await fetch('sid-ufpf.json');if(!response.ok)throw Error(`HTTP ${response.status}`);data=await response.json();render();}catch(error){$('ufStatus').textContent='The UFPF comparison could not load. Please reload the page; the other timeline tools remain available.';console.error(error);}
