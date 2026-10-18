import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {featuresAtLevel,changes,evidenceMatch} from './dist/sid-model.mjs';
const data=JSON.parse(readFileSync(new URL('./dist/sid-data.json',import.meta.url))),cards=JSON.parse(readFileSync(new URL('./dist/sid-evidence.json',import.meta.url)));
assert.equal(data.samples.length,5);assert.equal(data.features.length,236);
for(const s of data.samples){assert.ok(Math.abs(s.shares.reduce((a,b)=>a+b,0)-1)<1e-12);assert.equal(s.weights.reduce((a,b)=>a+b,0),s.total);s.shares.forEach((v,i)=>assert.equal(v,s.weights[i]/s.total));}
const genera=featuresAtLevel(data,'genus');for(let i=0;i<5;i++)assert.ok(Math.abs(genera.reduce((a,g)=>a+g.values[i],0)-1)<1e-12);
const deltas=changes(genera,0,4);assert.ok(Math.abs(deltas.reduce((a,g)=>a+g.delta,0))<1e-10);
const akk=deltas.find(r=>r.name==='Akkermansia');assert.ok(akk.end>.13&&akk.end<.14);assert.ok(Math.abs(akk.delta-(10940293/82887558-232927/12171318)*100)<1e-10);
const synthetic=changes([{id:'a',name:'a',values:[0,.1]},{id:'b',name:'b',values:[.2,0]},{id:'c',name:'c',values:[.01,.04]}],0,1);assert.equal(synthetic.find(r=>r.id==='a').fold,null);assert.equal(synthetic.find(r=>r.id==='a').status,'Newly detected');assert.equal(synthetic.find(r=>r.id==='b').status,'Not detected at end');assert.equal(synthetic.find(r=>r.id==='c').fold,4);
const sgb=featuresAtLevel(data,'sgb'),kle=sgb.find(r=>r.name==='Akkermansia sp KLE1798'),card=cards.find(c=>c.id==='depommier2019');assert.ok(kle);assert.match(evidenceMatch(card,kle,'sgb'),/Genus context only/);for(const f of sgb.filter(f=>!f.resolved))for(const c of cards)assert.equal(evidenceMatch(c,f,'sgb'),null);
assert.equal(data.samples[2].alternateDate,'2026-02-19');assert.equal(data.samples[3].dateStatus,'received');assert.equal(data.samples[4].dateStatus,'received');
for(const c of cards){assert.ok(c.claim&&c.limits&&c.deltaMeaning);assert.equal(new URL(c.url).protocol,'https:')}
console.log('Passed: source normalization, total-preserving aggregation, delta arithmetic, zero handling, Akkermansia profile, taxonomic evidence mismatch, unresolved exclusions, and uncertain-date provenance.');
