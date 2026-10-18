import {featuresAtLevel} from './sid-model.mjs';
import {bacteria} from './research-story-bacteria.mjs';
import {immuneStories,parseProfiles} from './immune-model.mjs';
import {eventsForInterval,ongoingAt,largestObservedChange,relevantTreatments,treatmentUncertainty} from './treatment-context.mjs';

const $=s=>document.querySelector(s), DAY=86400000;
export const uniqueEvents=events=>[...new Map(events.map(e=>[`${e.date}|${e.end_date||''}|${e.category}|${e.subcategory}|${e.title}|${e.dose||''}`,e])).values()];
export const overlaps=(e,start,end)=>Date.parse(e.date)<=Date.parse(end)&&Date.parse(e.end_date||e.date)>=Date.parse(start);
export const clipInterval=(e,start,end)=>({start:e.date<start?start:e.date,end:(e.end_date||e.date)>end?end:(e.end_date||e.date),point:!e.end_date});
export const eventsSincePrevious=(events,samples,index)=>index===0?ongoingAt(events,samples[0].date):eventsForInterval(events,samples[index-1].date,samples[index].date);
export function bacteroidesSummary(species,genusValues,from=3,to=4){
 const uniformis=species[0]?.sid?.map(v=>v*100)||[], theta=species[1]?.sid?.map(v=>v*100)||[];
 const genus=genusValues||[], net=genus[from]-genus[to], uniformisDrop=uniformis[from]-uniformis[to];
 return {uniformis,theta,genus,net,uniformisDrop,contribution:net>0?uniformisDrop/net*100:null};
}
const fmtDate=iso=>new Intl.DateTimeFormat('en-US',{month:'short',day:'numeric',year:'numeric',timeZone:'UTC'}).format(new Date(`${iso}T00:00:00Z`));
const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const timestamp=s=>Date.parse(`${s}T00:00:00Z`), label=s=>s.dateStatus==='received'?`${fmtDate(s.date)}†`:s.alternateDate?`${fmtDate(s.date)}*`:fmtDate(s.date);
const dateRange=s=>s.alternateDate?`${label(s)} (lab report: ${fmtDate(s.alternateDate)})`:label(s);
const shortLabel=s=>new Intl.DateTimeFormat('en-US',{month:'short',day:'numeric',timeZone:'UTC'}).format(new Date(`${s.date}T00:00:00Z`))+(s.dateStatus==='received'?'†':s.alternateDate?'*':'');

export function sourceStories(data,drivers){
 const genera=featuresAtLevel(data), byName=name=>genera.find(f=>f.name===name);
 const shared=drivers.levels.genus, control=drivers.scenarios['0'].drivers.Control.levels.genus;
 const ref=name=>{const i=shared.features.findIndex(f=>f.name===name);return {values:shared.sidValues[i],delta:control.rows[i].delta};};
 const akk=byName('Akkermansia'), bact=ref('Bacteroides'), bla=ref('Blautia'), fae=ref('Faecalibacterium');
 return {akk:{title:'Akkermansia',sub:'A late rise that raises a question, not an answer.',name:'Akkermansia (genus)',note:'Original Sid genus shares · the group includes the displayed Akkermansia sp KLE1798 annotation',guide:'<b>Akkermansia</b> is a group of gut bacteria that lives in the intestine’s mucus layer. The best-studied species, <em>A. muciniphila</em>, uses mucus as food. Cancer studies in lung and kidney cancers associated it with PD-1 immunotherapy response, but this Atlas label is KLE1798 and species equivalence is not established. <a target="_blank" rel="noopener" href="https://pubmed.ncbi.nlm.nih.gov/29097494/">Cancer study ↗</a>',values:akk.values.map(v=>v*100),color:'#b58b48',observed:`The final two reported dates show ${(akk.values[3]*100).toFixed(3)}% → ${(akk.values[4]*100).toFixed(3)}%.`,question:'Does the final Akkermansia increase coincide with a change in treatment experience?',key:'The line only joins five date-positioned observations; it does not show continuous sampling. * disputed February date; † lab receipt date.',conclusion:'A late sharp rise is observed. Symptom and outcome data are not linked here, and March collection dates are uncertain; this cannot establish cause or benefit.',next:'Confirm collection dates and processing quality, link symptom history to the intervals, and see whether the increase persists in the next sample.'},bact:{title:'Bacteroides',sub:'A first-to-last comparison misses the variation.',name:'Bacteroides',note:'Relative share retained from the separate REFERENCE-COMPARISON scale',guide:'<b>Bacteroides</b> is a large group of common gut bacteria; some members break down complex carbohydrates in food. Different species have been linked to different immunotherapy outcomes, and associations varied across melanoma cohorts. A whole-group rise or fall cannot be labeled good or bad. <a target="_blank" rel="noopener" href="https://pmc.ncbi.nlm.nih.gov/articles/PMC8938272/">Melanoma study ↗</a>',values:bact.values.map(v=>v*100),color:'#718166',observed:`Repeated swings are visible, ending at approximately ${(bact.values[4]*100).toFixed(2)}%.`,question:'Do the swings line up with documented exposures or symptoms?',key:'This series retains the reference-comparison scale and is not combined with the original within-Sid genus summary. * disputed February date; † lab receipt date.',conclusion:'The displayed shares rise and fall across five observations. There is no linked symptom record here, and no recovery or benefit judgment follows from this pattern.',next:'Validate comparable measurement, dates, and QC, then link interval-specific clinical observations before interpreting the swings.'},ref:{title:'Blautia + Faecalibacterium',sub:'Repeated differences need comparable evidence.',name:'Deviation from Control median',note:'Percentage points · Control reference n=16',guide:'<b>Blautia</b> includes bacteria that ferment food; some make acetate. <b>Faecalibacterium</b> includes bacteria that make butyrate. Studies in melanoma connect both groups with immune-treatment outcomes, but the findings are other-cancer research, not recurrence measurements or treatment targets for Sid. <a target="_blank" rel="noopener" href="https://pmc.ncbi.nlm.nih.gov/articles/PMC13384208/">Blautia study ↗</a> · <a target="_blank" rel="noopener" href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5827966/">Faecalibacterium study ↗</a>',values:[bla.delta.map(v=>v*100),fae.delta.map(v=>v*100)],labels:['Blautia','Faecalibacterium'],colors:['#57766b','#9c7460'],observed:'Blautia and Faecalibacterium are both below the Control reference median at 5 of 5 visits.',question:'Are recurring differences meaningful after accounting for study and processing differences?',key:'Both lines are Sid share minus the Control median. Processing equivalence is unverified; no health or disease classification is implied.',conclusion:'Both displayed taxa are below the reference median at 5 of 5 observations. That recurring difference needs compatible study and processing evidence before it can be interpreted clinically.',next:'Check comparison compatibility and processing details, then seek replication rather than treating a reference difference as a classification.'}};
}

function chart(story,samples,selected,onPick){
 const multi=Array.isArray(story.values[0]),values=multi?story.values.flat():story.values,known=values.filter(Number.isFinite);
 const percentSeries=multi&&story.unit==='%',min=percentSeries?0:multi?Math.min(-10,...known):0,ref=story.reference;
 const max=multi?Math.max(1,...known)*1.14:Math.max(.001,...known,ref?.q75||0,ref?.median||0)*1.14;
 const W=720,H=245,L=48,R=16,T=22,B=37,t0=timestamp(samples[0].date),t1=timestamp(samples.at(-1).date);
 const x=i=>L+(timestamp(samples[i].date)-t0)*(W-L-R)/(t1-t0),y=v=>T+(max-v)*(H-T-B)/(max-min);
 const decimals=multi?0:max<.01?4:max<1?3:1;
 let out=`<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(story.name)} as percent of the retained microbial profile across five reported dates; highlighted selected interval${selected?` from ${esc(shortLabel(samples[selected-1]))} to ${esc(shortLabel(samples[selected]))}`:''}">`;
 if(selected){const focusStart=x(selected-1),focusWidth=x(selected)-x(selected-1);out+=`<rect class="focus-window" x="${focusStart}" y="${T}" width="${focusWidth}" height="${H-T-B}"/>`;}
 if(known.length){
  if(ref&&Number.isFinite(ref.q25)&&Number.isFinite(ref.q75))out+=`<rect class="control-band" x="${L}" y="${y(ref.q75)}" width="${W-L-R}" height="${y(ref.q25)-y(ref.q75)}" fill="#aeb6b1" opacity=".28"/>`;
  [min,(min+max)/2,max].forEach(v=>out+=`<line x1="${L}" x2="${W-R}" y1="${y(v)}" y2="${y(v)}" stroke="#d9d4c6"/><text x="${L-7}" y="${y(v)+4}" text-anchor="end">${v.toFixed(decimals)}${story.unit||(multi?' pp':'%')}</text>`);
  if(ref&&Number.isFinite(ref.median))out+=`<line class="control-median" x1="${L}" x2="${W-R}" y1="${y(ref.median)}" y2="${y(ref.median)}" stroke="#66736e" stroke-dasharray="5 4"/>`;
  if(multi&&!percentSeries)out+=`<line x1="${L}" x2="${W-R}" y1="${y(0)}" y2="${y(0)}" stroke="#718166"/>`;
  (multi?story.values:[story.values]).forEach((series,j)=>{
   const c=multi?story.colors[j]:story.color;
   const path=series.map((v,i)=>Number.isFinite(v)?`${i&&Number.isFinite(series[i-1])?'L':'M'}${x(i)},${y(v)}`:'').join(' ');
   out+=`<path d="${path}" fill="none" stroke="${c}" stroke-width="${multi?2.4:3}"/>`;
   series.forEach((v,i)=>{if(Number.isFinite(v))out+=`<g class="point" data-visit="${i}" data-series="${j}" tabindex="0" role="button" aria-label="${esc(label(samples[i]))}: ${multi?story.labels[j]+' ':''}${v.toFixed(multi?2:3)}${percentSeries||!multi?' percent of the retained microbial profile':' percentage points'}"><circle cx="${x(i)}" cy="${y(v)}" r="${i===selected?6:4.5}" fill="${c}" stroke="${i===selected?'#163b36':'#fcfbf6'}" stroke-width="${i===selected?2:1}"/></g>`;});
  });
 }else out+=`<text x="${L}" y="100">No matched measurements are available for this candidate.</text>`;
 samples.forEach((s,i)=>out+=`<text x="${x(i)}" y="${H-12}" text-anchor="${i===samples.length-1?'end':'middle'}">${esc(shortLabel(s))}</text>`);
 if(ref&&ref.n)out+=`<text x="${L}" y="13">Sid: gold line · grey: Control middle 50% · dashed: median</text>`;
 out+='</svg>';$('#chart').innerHTML=out;
 $('#chart').querySelectorAll('.point').forEach(el=>{const choose=()=>onPick(+el.dataset.visit);el.onclick=choose;el.onkeydown=e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();choose();}}});
}

function timeline(events,samples,enabled,selected){
 const start=samples[0].date,end=samples.at(-1).date,W=720,L=48,R=16;
 const t0=timestamp(start),t1=timestamp(end),x=d=>L+(timestamp(d)-t0)*(W-L-R)/(t1-t0);
 const rows=[];
 for(const type of enabled){
  rows.push({type,heading:true});
  for(const title of [...new Set(events[type].map(e=>e.title))]) rows.push({type,title,events:events[type].filter(e=>e.title===title)});
 }
 const H=34+rows.reduce((n,r)=>n+(r.heading?30:42),0);
 const interval=selected?eventsForInterval(Object.values(events).flat(),samples[selected-1].date,samples[selected].date):ongoingAt(Object.values(events).flat(),samples[0].date);
 const focusKeys=new Set(interval.map(e=>`${e.date}|${e.end_date||''}|${e.title}`));
 let svg=`<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Documented treatment periods and individual doses, aligned to reported stool dates${selected?`; highlighted selected interval from ${esc(shortLabel(samples[selected-1]))} to ${esc(shortLabel(samples[selected]))}`:''}">${selected?`<rect class="focus-window" x="${x(samples[selected-1].date)}" y="19" width="${x(samples[selected].date)-x(samples[selected-1].date)}" height="${H-27}"/>`:''}`;
 samples.forEach((s,i)=>svg+=`<line data-sample="${i}" x1="${x(s.date)}" x2="${x(s.date)}" y1="19" y2="${H-8}" stroke="#d9d4c6"/><text x="${x(s.date)}" y="12" text-anchor="${i===samples.length-1?'end':'middle'}">${esc(shortLabel(s))}</text>`);
 let y=25;
 for(const row of rows){
  if(row.heading){svg+=`<text class="lane-category" x="${L}" y="${y+17}">${esc(row.type)}</text>`;y+=30;continue;}
  const title=row.title==='mRNA neoantigen vaccine'?'mRNA neoantigen vaccine':row.title;
  svg+=`<text class="lane-drug" x="${L}" y="${y+10}">${esc(title)}</text><line x1="${L}" x2="${W-R}" y1="${y+26}" y2="${y+26}" stroke="#d9d4c6"/>`;
  for(const e of row.events){
   const c=clipInterval(e,start,end),tip=esc(e.title+' · '+fmtDate(e.date)+(e.end_date?'–'+fmtDate(e.end_date):'')+(e.dose?' · '+e.dose:''));
   const dim=focusKeys.has(`${e.date}|${e.end_date||''}|${e.title}`)?'':' class="event-dim"';
   svg+=c.point?`<g${dim}><title>${tip}</title><line x1="${x(c.start)}" x2="${x(c.start)}" y1="${y+19}" y2="${y+33}" stroke="#b58b48"/><circle cx="${x(c.start)}" cy="${y+26}" r="4" fill="#b58b48"/></g>`:`<g${dim}><title>${tip}</title><line x1="${x(c.start)}" x2="${x(c.end)}" y1="${y+26}" y2="${y+26}" stroke="#b58b48" stroke-width="6" stroke-linecap="round"/></g>`;
  }
  y+=42;
 }
 return `<div class="timeline-scroll" tabindex="0" role="region" aria-label="Treatment chart; scroll horizontally on a narrow screen">${svg}</svg></div>`;
}

async function init(){const [data,drivers,eventData,registry,treatmentEvidence,bactEvidence,antibioticContext,tsv]=await Promise.all(['sid-data.json','sid-ufpf-drivers.json','sid-events.json','immune-evidence.json','treatment-evidence.json','bacteroides-species.json','antibiotic-context.json'].map(u=>fetch(u).then(r=>r.json())).concat(fetch('sid-ufpf-profiles.tsv').then(r=>r.text()))),samples=data.samples,all=uniqueEvents(eventData.events),start=samples[0].date,end=samples.at(-1).date;
 const profiles=parseProfiles(tsv),allFeatures=[...drivers.levels.genus.features,...drivers.levels.sgb.features];
 const immune=immuneStories(registry.candidates,profiles,allFeatures,samples.map(s=>s.id));
 const immuneMap=Object.fromEntries(immune.map(r=>[r.id,{
  ...r,title:r.name,sub:`${r.availability} · ${r.direction}`,name:r.name,
  note:'Reference-comparison profiles · share (%)',values:r.sid.map(v=>v===null?null:v*100),color:'#b58b48',
  observed:r.availability,question:`What does ${r.name} research let us ask about Sid’s profile?`,
  key:'',
  conclusion:'These are research leads. Neither a reference difference nor a non-detection establishes Sid’s immune strength, treatment response, or recurrence risk.',
  next:'Confirm the species and processing match, then link the measured intervals to clinical outcomes or direct immune measurements.',
  reference:{n:r.reference.n,...Object.fromEntries(['q25','median','q75'].map(k=>[k,r.reference[k]===null?null:r.reference[k]*100]))}
 }]));
 const speciesNames=['Bifidobacterium bifidum','Bifidobacterium adolescentis','Bifidobacterium animalis','Bifidobacterium longum','Bifidobacterium breve'];
 const bifidSpecies=immuneStories(speciesNames.map((name,i)=>({id:`bifid-${i}`,name,rank:'species',taxon:name})),profiles,allFeatures,samples.map(s=>s.id));
 const bactSpecies=immuneStories(bactEvidence.entries.map((entry,index)=>({id:`bact-species-${index}`,name:entry.name,rank:'species',taxon:entry.name})),profiles,allFeatures,samples.map(s=>s.id));
 immuneMap.bact.genusValues=immuneMap.bact.values;immuneMap.bact.values=[bactSpecies[0].sid.map(v=>v*100),bactSpecies[1].sid.map(v=>v*100)];immuneMap.bact.labels=['B. uniformis','B. thetaiotaomicron'];immuneMap.bact.colors=['#718166','#a36b52'];immuneMap.bact.name='Bacteroides species shares';immuneMap.bact.summary='Most of the final Bacteroides decline comes from B. uniformis, while B. thetaiotaomicron increases. Cefazolin and a hip procedure are recorded within this interval, making a treatment-related shift plausible; uncertain collection dates prevent confirming the cause. This lower relative share alone does not show worsening health or a fall in total bacterial count.';immuneMap.bact.reference=null;immuneMap.bact.unit='%';
 const bactSummary=bacteroidesSummary(bactSpecies,immuneMap.bact.genusValues);
 const mode='immune';
 let stories=immuneMap;

 const lanes={"Weight / metabolic medicines":uniqueEvents(all.filter(e=>['Zepbound','Wegovy','Metformin'].includes(e.title)&&overlaps(e,start,end))),Antibiotics:uniqueEvents(all.filter(e=>e.subcategory==='Antibiotics'&&overlaps(e,start,end))),Procedures:uniqueEvents(all.filter(e=>e.subcategory==='Surgery'&&overlaps(e,start,end))),"Cancer immunotherapy":uniqueEvents(all.filter(e=>(e.title==='Dostarlimab'||/neoantigen vaccine/i.test(e.title))&&overlaps(e,start,end)))};
 const treatmentByTitle=Object.fromEntries(treatmentEvidence.treatments.map(item=>[item.title,item]));
 const treatmentSection=$('.lanes');treatmentSection.id='treatments';$('#story').after(treatmentSection);
 let active='bact',selected=4,treatmentSelected=4;let tabOrder=['bact','akk','bifido','blautia','faecali'];
 const renderLanes=()=>{
  const focusedTreatment=document.activeElement?.dataset?.treatmentVisit,timelineLeft=$('.timeline-scroll')?.scrollLeft||0;
  const prior=samples[treatmentSelected-1],current=samples[treatmentSelected],events=eventsForInterval(all,prior.date,current.date),choices=relevantTreatments(events,treatmentByTitle,'shared');
  if(!$('#treatmentIntervals')){$('#intervalSummary').before(Object.assign(document.createElement('div'),{id:'treatmentIntervals',className:'visit-buttons'}));}
  $('#intervalSummary').textContent=`${dateRange(prior)} → ${dateRange(current)} · shared treatment interval`;
  $('#treatmentIntervals').innerHTML=samples.slice(1).map((sample,index)=>`<button data-treatment-visit="${index+1}" aria-pressed="${index+1===treatmentSelected}">${esc(shortLabel(samples[index]))} → ${esc(shortLabel(sample))}</button>`).join('');
  $('#treatmentIntervals').querySelectorAll('button').forEach(button=>button.onclick=()=>{treatmentSelected=+button.dataset.treatmentVisit;renderLanes();});
  const titles=new Set(choices.map(item=>item.title)),sourceDates=events.filter(event=>titles.has(event.title)),sourceDate=event=>event.end_date?`${esc(event.title)}: ${event.date<=prior.date?'already recorded before this interval':'recorded dates'} (${esc(fmtDate(event.date))}–${esc(fmtDate(event.end_date))})`:`${esc(event.title)} recorded ${esc(fmtDate(event.date))}`;
  const uncertainty=treatmentUncertainty(prior,current);
  const antibioticNote=antibioticContext.intervals[String(treatmentSelected)];
  $('#treatmentCallout').innerHTML=`${choices.map(item=>`<section class="treatment-explanation"><strong>${esc(item.title)}</strong><p>${esc(item.explanation)} <a href="${esc(item.source.url)}" target="_blank" rel="noopener">${esc(item.source.label)} ↗</a></p></section>`).join('')||'<p>No selected treatment explanation is available for this interval.</p>'}${antibioticNote?`<section class="treatment-explanation antibiotic-context"><strong>${esc(antibioticNote.title)}</strong><p>${esc(antibioticNote.text)} <a href="${esc(antibioticContext.source.url)}" target="_blank" rel="noopener">${esc(antibioticContext.source.label)} ↗</a></p></section>`:''}${sourceDates.length?`<p class="small">Source-record dates: ${sourceDates.map(sourceDate).join('; ')}. Periods are reported periods, not daily dosing.</p>`:''}${uncertainty?`<p class="small">${esc(uncertainty)}</p>`:''}`;
  $('#timelineLanes').innerHTML=timeline(lanes,samples,Object.keys(lanes),treatmentSelected)+`<details class="event-details"><summary>Read all source events by treatment group</summary><p class="small">Source snapshot retrieved Oct 17, 2026 from <a href="https://osteosarc.com/timeline/" target="_blank" rel="noopener">the public timeline ↗</a>.</p>${Object.entries(lanes).map(([name,items])=>`<h4>${esc(name)}</h4><ul class="small">${items.map(event=>`<li>${esc(fmtDate(event.date))}${event.end_date?`–${esc(fmtDate(event.end_date))}`:''} · ${esc(event.title)}</li>`).join('')}</ul>`).join('')}</details>`;
  const laneScroll=$('.timeline-scroll');
  const keepLabelsVisible=()=>{const width=laneScroll.querySelector('svg').getBoundingClientRect().width;laneScroll.querySelectorAll('.lane-category,.lane-drug').forEach(text=>text.setAttribute('x',48+laneScroll.scrollLeft*720/width));};
  laneScroll.scrollLeft=timelineLeft;laneScroll.onscroll=keepLabelsVisible;keepLabelsVisible();
  if(focusedTreatment)$('#treatmentIntervals').querySelector(`[data-treatment-visit="${focusedTreatment}"]`)?.focus({preventScroll:true});
 };
 const renderVisit=()=>{
  const focusedVisit=document.activeElement?.dataset?.visit,s=samples[selected],previous=samples[selected-1],story=stories[active],multi=Array.isArray(story.values[0]);
  const values=multi?story.values.map((row,index)=>({label:story.labels[index],before:row[selected-1],after:row[selected]})):[{label:story.name,before:story.values[selected-1],after:story.values[selected]}];
  let bacteriaSummary=$('#bacteriaIntervalSummary');if(!bacteriaSummary){bacteriaSummary=document.createElement('p');bacteriaSummary.id='bacteriaIntervalSummary';bacteriaSummary.className='small';$('#visits').before(bacteriaSummary);}
  bacteriaSummary.textContent=selected?`${dateRange(previous)} → ${dateRange(s)} · ${values.map(v=>`${v.label}: ${Number.isFinite(v.before)?v.before.toFixed(3)+'%':'unavailable'} → ${Number.isFinite(v.after)?v.after.toFixed(3)+'%':'unavailable'}`).join(' · ')}`:`No preceding reported sample before ${dateRange(s)}.`;
  const candidateIsFlat=active==='a-muciniphila'&&values.every(value=>value.before===0&&value.after===0);
  if(selected===0){
   const ongoing=ongoingAt(all,s.date);
   bacteriaSummary.textContent=`No preceding reported sample before ${dateRange(s)}; no measured change is available.`;
  }else if(candidateIsFlat){
   bacteriaSummary.textContent=`${dateRange(previous)} → ${dateRange(s)}. No measured change: A. muciniphila is 0.000% at both reported dates.`;
  }else{
   const unit=story.unit||(multi?' pp':'%'),digits=multi?2:3;
   bacteriaSummary.textContent=`${dateRange(previous)} → ${dateRange(s)} · `+values.map(value=>`${value.label}: ${Number.isFinite(value.before)?value.before.toFixed(digits):'not available'}${unit} → ${Number.isFinite(value.after)?value.after.toFixed(digits):'not available'}${unit}${Number.isFinite(value.before)&&Number.isFinite(value.after)?` (${(value.after-value.before)>=0?'+':''}${(value.after-value.before).toFixed(digits)} pp)`:''}`).join(' · ');
   const inWindow=eventsForInterval(all,previous.date,s.date),choices=relevantTreatments(inWindow,treatmentByTitle,active);
   const uncertainty=treatmentUncertainty(previous,s);
   const choiceTitles=new Set(choices.map(item=>item.title));
   const sourceDates=inWindow.filter(event=>choiceTitles.has(event.title));
   const sourceDateLabel=event=>event.end_date?`${esc(event.title)}: ${event.date<=previous.date?'already recorded before this interval':'recorded dates'} (${esc(fmtDate(event.date))}–${esc(fmtDate(event.end_date))})`:`${esc(event.title)} recorded ${esc(fmtDate(event.date))}`;
  }
  $('#visits').innerHTML=samples.map((sample,index)=>`<button data-visit="${index}" aria-pressed="${index===selected}">${esc(shortLabel(sample))}</button>`).join('');
  $('#visits').querySelectorAll('button').forEach(button=>button.onclick=()=>{selected=+button.dataset.visit;render();});
  if(focusedVisit)$('#visits').querySelector(`[data-visit="${focusedVisit}"]`)?.focus({preventScroll:true});
  const largest=largestObservedChange(story.values),jump=$('#largestChange');
  jump.hidden=!largest;jump.disabled=!largest;
  if(largest)jump.onclick=()=>{selected=largest.index;render();};
 };
 function render(){const focus=document.activeElement?.dataset?.story,s=stories[active];$('#question').textContent=s.question;let guide=$('#bacteriaGuide');if(!guide){guide=document.createElement('aside');guide.id='bacteriaGuide';guide.className='bacteria-guide';$('#question').after(guide)}if(mode==='immune'){
 guide.innerHTML=`<section class="research-card"><p class="eyebrow">${esc(s.direction)}</p><h3>${esc(s.type)}</h3><p class="research-context">${esc(s.population)} · ${esc(s.treatment)}</p><p>${esc(s.claim)}</p><dl><dt>What was studied</dt><dd>${esc(s.studied)}</dd><dt>How closely does Sid’s measurement match?</dt><dd>${esc(s.matchNote)}</dd></dl><p class="research-limit">${esc(s.limits)}</p><p class="research-sources">${s.sources.map(source=>`<a href="${esc(source.url)}" target="_blank" rel="noopener">${esc(source.label)} ↗</a>`).join(' · ')}</p></section>`;
 }else guide.innerHTML=bacteria[active].map(b=>`<section class="bacterium"><h3>${b.name}, in plain language</h3><div class="bio-columns"><div><span class="eyebrow">What it is</span><p>${b.what} <a href="${b.definitionUrl}" target="_blank" rel="noopener">Biology source ↗</a></p></div><div><span class="eyebrow">Why cancer researchers care</span><p>${b.why}</p><span class="bio-evidence">${b.evidence}</span><a href="${b.researchUrl}" target="_blank" rel="noopener">${b.researchLabel} ↗</a></div></div><p class="sid-meaning"><strong>What this means for Sid:</strong> ${b.sid}</p></section>`).join('');
 $('#observed').textContent=s.observed;$('#chartName').textContent=s.name;$('#chartNote').textContent=s.note;$('#chartKey').textContent=s.key;$('#storyPanel').setAttribute('aria-labelledby','tab-'+active);chart(s,samples,selected,i=>{const focus=document.activeElement?.closest('.point'),series=focus?.dataset.series;selected=i;render();if(focus)$('#chart').querySelector(`[data-visit="${i}"][data-series="${series}"]`)?.focus({preventScroll:true});});const legend=$('#chartLegend'),multi=Array.isArray(s.values[0]);legend.hidden=!multi;legend.innerHTML=multi?s.labels.map((name,index)=>`<span><i style="--series-color:${esc(s.colors[index])}" aria-hidden="true"></i>${esc(name)}</span>`).join(''):'';renderVisit();const chartScroll=$('#chart'),plot=$('#chart svg'),fraction=(timestamp(samples[selected].date)-timestamp(samples[0].date))/(timestamp(samples.at(-1).date)-timestamp(samples[0].date));chartScroll.scrollLeft=fraction*plot.getBoundingClientRect().width-chartScroll.clientWidth/2;$('#storyTabs').innerHTML=tabOrder.map((id,i)=>`<button class="story-tab" id="tab-${id}" role="tab" aria-controls="storyPanel" tabindex="${id===active?0:-1}" aria-selected="${id===active}" data-story="${id}"><span>0${i+1} / ${mode==='immune'?'research candidate':stories[id].title==='Akkermansia'?'change at the end':stories[id].title==='Bacteroides'?'swings over time':'recurring reference gap'}</span><strong>${esc(stories[id].title)}</strong><small>${esc(stories[id].sub)}</small></button>`).join('');[...$('#storyTabs').children].forEach((b,i)=>{const choose=()=>{active=b.dataset.story;selected=4;$('#storyPanel').scrollTop=0;render()};b.onclick=choose;b.onkeydown=e=>{let n;if(['ArrowRight','ArrowDown'].includes(e.key))n=(i+1)%tabOrder.length;if(['ArrowLeft','ArrowUp'].includes(e.key))n=(i+tabOrder.length-1)%tabOrder.length;if(e.key==='Home')n=0;if(e.key==='End')n=tabOrder.length-1;if(n!==undefined){e.preventDefault();active=tabOrder[n];selected=4;$('#storyPanel').scrollTop=0;render();document.querySelectorAll('.story-tab')[n].focus()}}});if(focus)$('#storyTabs').querySelector(`[data-story="${focus}"]`)?.focus();}
 const baseRender=render;
 render=function(){
  baseRender();const story=stories[active],panel=$('#storyPanel'),chartCard=$('.chart-card'),guide=$('#bacteriaGuide');
  const immuneMode=mode==='immune';
  $('.chart-top').hidden=true;$('#chartKey').hidden=immuneMode;
  $('#question').hidden=immuneMode;$('#question').previousElementSibling.hidden=immuneMode;
  $('#observed').parentElement.hidden=immuneMode;
  let heading=$('#candidateTitle'),lead=$('#researchSummary');
  if(!heading){heading=document.createElement('h2');heading.id='candidateTitle';lead=document.createElement('p');lead.id='researchSummary';lead.className='research-summary';chartCard.before(heading,lead);}
  heading.hidden=lead.hidden=!immuneMode;
  let detail=$('#researchDetails');
  if(!detail){detail=document.createElement('details');detail.id='researchDetails';detail.innerHTML='<summary>Study details and sources</summary>';chartCard.after(detail);}
  detail.hidden=!immuneMode;
  if(immuneMode){heading.textContent=story.title;lead.textContent=story.summary;detail.append(guide);}
  else $('#question').after(guide);

  let summary=$('#controlSummary');if(!summary){summary=document.createElement('p');summary.id='controlSummary';summary.className='small';$('#chartKey').before(summary);}
  summary.hidden=mode!=='immune'||!story.reference;
  if(mode==='immune'&&story.reference){
   const f=v=>Number.isFinite(v)?v.toFixed(3)+'%':'unavailable';
   summary.textContent=`Control median ${f(story.reference.median)} · middle 50% ${f(story.reference.q25)}–${f(story.reference.q75)}`;
  }
  let table=$('#speciesBreakdown');
  if(mode==='immune'&&active==='bact'){
   if(!table){table=document.createElement('details');table.id='speciesBreakdown';table.className='species-breakdown';$('.chart-card').after(table);}
   const wholeFrom=Number.isFinite(bactSummary.genus[selected])?bactSummary.genus[selected].toFixed(3)+'%':'unavailable';
   const allFive=bactSpecies.reduce((total,item)=>total+(Number.isFinite(item.sid[selected])?item.sid[selected]*100:0),0);
   table.innerHTML=`<summary>Species roles, cautions, and selected-date breakdown</summary><div class="species-content"><p class="eyebrow">Bacteroides species · ${esc(shortLabel(samples[selected]))}</p><p class="small">Whole Bacteroides group: ${wholeFrom} of the retained profile. The five displayed species total ${allFive.toFixed(3)}%; they are included within the whole group, not additional compositional groups. The final whole-group change is ${bactSummary.genus[3].toFixed(3)}% → ${bactSummary.genus[4].toFixed(3)}%; B. uniformis accounts for about ${bactSummary.contribution.toFixed(0)}% of that net decline.</p><p class="small">${esc(bactEvidence.scope)}</p><table class="bact-species-table"><thead><tr><th>Species and context</th><th>Profile share</th><th>Control median</th></tr></thead><tbody>${bactEvidence.entries.map((entry,index)=>{const record=bactSpecies[index],share=record.sid[selected],median=record.reference.median;return `<tr><td class="species-context"><strong>${esc(entry.name)}</strong><small>${esc(entry.usualRole)}</small><small>${esc(entry.caution)}</small><small>${entry.sources.map(source=>`<a href="${esc(source.url)}" target="_blank" rel="noopener">${esc(source.label)} ↗</a>`).join(' · ')}</small></td><td>${Number.isFinite(share)?(share*100).toFixed(3)+'%':'Not available'}</td><td>${Number.isFinite(median)?(median*100).toFixed(3)+'%':'Not available'}</td></tr>`;}).join('')}</tbody></table></div>`;
  }else if(mode==='immune'&&active==='bifido'){
   if(!table){table=document.createElement('details');table.id='speciesBreakdown';table.className='species-breakdown';$('.chart-card').after(table);}
   table.innerHTML=`<summary>Detailed species distribution</summary><div class="species-content"><p class="eyebrow">Bifidobacterium species · ${esc(shortLabel(samples[selected]))}</p><p class="small">The group total can be dominated by other species. B. longum and B. breve are among those studied in the linked cancer research. Non-detection does not prove biological absence.</p><table><thead><tr><th>Reported species</th><th>Profile share</th></tr></thead><tbody>${bifidSpecies.map(r=>`<tr><td>${esc(r.name)}${/longum|breve/.test(r.name)?'<small>Among the studied species</small>':''}</td><td>${Number.isFinite(r.sid[selected])?(r.sid[selected]*100).toFixed(3)+'%':'Not available'}</td></tr>`).join('')}</tbody></table></div>`;
  }else table?.remove();
  if(table&&table.isConnected)detail.after(table);
 };
 document.querySelector('#methods .details-grid').insertAdjacentHTML('beforeend',`<details><summary>How immune-response research is selected</summary><p class="small">${esc(registry.scope)}</p><p class="small">Reviewed ${esc(registry.reviewed)}. New studies can be added to the reviewed list and matched to the measured taxon; no automatic response prediction is run.</p><p class="small"><a href="immune-evidence.json" download>Reviewed evidence list</a> · <a href="sid-ufpf-profiles.tsv" download>Comparison profiles</a></p></details>`);
 renderLanes();render();}
if(typeof document!=='undefined')init().catch(e=>{$('#storyPanel').innerHTML=`<p>Unable to load the local research-story data: ${esc(e.message)}</p>`});
