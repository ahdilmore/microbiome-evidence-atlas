import {featuresAtLevel,changes} from './sid-model.mjs';
import {stability,adjacentChanges} from './sid-longitudinal-model.mjs';

const $=id=>document.getElementById(id);
const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const pct=n=>(100*n).toFixed(1)+'%';
const date=s=>new Intl.DateTimeFormat('en-US',{month:'short',day:'numeric',year:'numeric',timeZone:'UTC'}).format(new Date(s));
async function read(path){const r=await fetch(path);if(!r.ok)throw Error(`Unable to load ${path}`);return r.json();}

async function sidSummary(){
 try{
  const data=await read('sid-data.json'),features=featuresAtLevel(data,'genus');
  const biggest=adjacentChanges(features,data.samples).sort((a,b)=>b.brayCurtis-a.brayCurtis)[0];
  const named=changes(features,0,data.samples.length-1).filter(f=>f.resolved);
  const fall=named.find(f=>f.delta<0),rise=named.find(f=>f.delta>0);
  $('overviewChange').innerHTML=`<strong class="overview-value">${esc(date(data.samples[biggest.start].date))} → ${esc(date(data.samples[biggest.end].date))}</strong><p>The largest observed between-visit shift: <b>${biggest.brayCurtis.toFixed(3)} Bray–Curtis dissimilarity</b> across genus groups (0 = identical; 1 = disjoint).</p><p>Across the full ${esc(date(data.samples[0].date))} → ${esc(date(data.samples.at(-1).date))} window, the largest named-genus decrease was <i>${esc(fall.name)}</i> (${pct(fall.start)} → ${pct(fall.end)}); the largest increase was <i>${esc(rise.name)}</i> (${pct(rise.start)} → ${pct(rise.end)}).</p>`;
  const persistent=stability(features,.1),steady=persistent[0];
  $('overviewPersistence').innerHTML=`<strong class="overview-value">${persistent.length} named genera at every visit</strong><p>Each stayed at or above <b>0.1% share in all ${data.samples.length} profiles</b>, including ${persistent.slice(0,2).map(f=>`<i>${esc(f.name)}</i>`).join(' and ')}.</p><p>Presence persisted, but abundance varied. <i>${esc(steady.name)}</i>, the least variable by range relative to median among these genera, still ranged from <b>${pct(steady.min)} to ${pct(steady.max)}</b>.</p>`;
 }catch(error){for(const id of ['overviewChange','overviewPersistence'])$(id).textContent='Summary unavailable. The detailed sections below provide the source data and methods.';console.error(error);}
}

async function referenceSummary(){
 try{
  const data=await read('sid-ufpf.json'),full=data.scenarios['all-0'],metric=full.metrics.unifrac;
  const disjoint=metric.hullsDisjoint['0,1'];
  const controls=data.scenarios['Control-0'].metrics.unifrac.hullsDisjoint['0,1'];
  const screened=data.scenarios['Control-95'].metrics.unifrac.hullsDisjoint['0,1'];
  const sensitive=controls&&!screened?'Controls alone appear separate on these axes, but that separation disappears with the ≥95% retained-abundance screen.':'Reference selection and retained-abundance screening can change the projection; inspect the sensitivity views below.';
  $('overviewReference').innerHTML=`<strong class="overview-value">${disjoint?'Separate in this projection':'Overlap with the UF cohort'}</strong><p>Sid’s five visits ${disjoint?'have a separate envelope from':'overlap the envelope of'} <b>${full.indices.length-5} UFPF profiles</b> on weighted UniFrac PCoA axes 1–2. Those axes show <b>${pct(metric.axisFractions[0]+metric.axisFractions[1])}</b> of corrected inertia.</p><p>${sensitive}</p><p class="overview-caveat">Same reference database; full processing equivalence is unverified. This comparison does not establish health or disease status.</p>`;
 }catch(error){$('overviewReference').textContent='Reference summary unavailable. Open the UFPF comparison below for its data and methods.';console.error(error);}
}

await Promise.all([sidSummary(),referenceSummary()]);
$('sid-overview').setAttribute('aria-busy','false');
