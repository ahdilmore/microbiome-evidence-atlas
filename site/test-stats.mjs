import fs from 'node:fs';
import assert from 'node:assert/strict';
import {ranks,bh,analyze,bootstrapDelta} from './dist/stats.mjs';
assert.deepEqual(ranks([4,1,1,3]),[4,1.5,1.5,3]);
assert.deepEqual(bh([.01,.04,.03]),[.03,.04,.04]);
assert.deepEqual(bootstrapDelta([2,2],[1,1]),[100,100]);
const data=JSON.parse(fs.readFileSync(new URL('./dist/data.json',import.meta.url)));
assert.equal(data.provenance.taxonomy,'Greengenes 13_5');
for(const sample of data.samples){assert(sample.librarySize>=1000);assert(sample.counts.every(Number.isInteger));assert.equal(sample.counts.reduce((a,b)=>a+b,0)+sample.unassignedGenusReads,sample.librarySize);sample.values.forEach((v,i)=>assert(Math.abs(v-sample.counts[i]/sample.librarySize)<1e-12));}
for(const field of ['ir','sex','age','bmi','sspg']){
 const a=analyze(data,field);
 assert(a.results.every(r=>r.p>0&&r.p<=r.q+1e-12&&r.q<=1&&Math.abs(r.effect)<=1));
 assert.equal(new Set(a.people.map(p=>p.id)).size,a.n);
 if(field==='ir'){
  assert.equal(a.n,50);assert.equal(a.n1,27);assert.equal(a.n0,23);
  assert.deepEqual(analyze(data,field),a);
  for(const result of a.results){
   const x=a.people.filter(p=>p.value==='IR').map(p=>p.values[result.index]);
   const y=a.people.filter(p=>p.value==='IS').map(p=>p.values[result.index]);
   let pairwise=0;for(const u of x)for(const v of y)pairwise+=(u>v?1:u<v?-1:0);
   assert(Math.abs(pairwise/(x.length*y.length)-result.effect)<1e-12);
  }
 }
}
const fixture={taxa:['separated','constant'],samples:Array.from({length:20},(_,i)=>({id:'sample'+i,subject:'p'+i,match:'exact',event:'Healthy',ir:i<10?'IS':'IR',values:[i<10?0:1,.1]}))};
const result=analyze(fixture,'ir');assert.equal(result.results.length,1);assert.equal(result.results[0].effect,1);assert(result.results[0].p<.01);
const evidence=JSON.parse(fs.readFileSync(new URL('./dist/literature.json',import.meta.url)));
for(const e of evidence){for(const key of ['genus','title','claim','scope','studyType','population','limitation','url','doi','year'])assert(e[key]);assert(e.url.startsWith('https://'));}
console.log('Passed: cohort inclusion, average ranks, BH correction, pairwise effect verification, reproducibility, complete separation, constant exclusion, bootstrap and evidence schema.');
