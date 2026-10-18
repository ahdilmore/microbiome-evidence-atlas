import assert from 'node:assert/strict';import fs from 'node:fs';
import {featuresAtLevel,evidenceMatch,changes} from './dist/sid-model.mjs';
import {distanceFromFirst} from './dist/sid-longitudinal-model.mjs';
import {childFeatures} from './dist/sid-scale.mjs';
const d=JSON.parse(fs.readFileSync('dist/sid-data.json')),p=JSON.parse(fs.readFileSync('dist/sid-phylogeny.json'));
for(const rank of ['sgb','genus','family']){const f=featuresAtLevel(d,rank);for(let i=0;i<5;i++)assert.ok(Math.abs(f.reduce((s,r)=>s+r.values[i],0)-1)<1e-12);assert.equal(distanceFromFirst(f,d.samples)[0].brayCurtis,0);if(rank!=='sgb')for(const parent of f){const children=childFeatures(d,parent,rank);for(let i=0;i<5;i++)assert.ok(Math.abs(children.reduce((s,c)=>s+c.values[i],0)-parent.values[i])<1e-12)}}
const species=distanceFromFirst(featuresAtLevel(d,'sgb'),d.samples),genus=distanceFromFirst(featuresAtLevel(d,'genus'),d.samples),family=distanceFromFirst(featuresAtLevel(d,'family'),d.samples);
for(let i=0;i<5;i++){assert.ok(family[i].brayCurtis<=genus[i].brayCurtis+1e-12);assert.ok(genus[i].brayCurtis<=species[i].brayCurtis+1e-12)}
assert.ok(Math.abs(genus[1].brayCurtis-.599492883089399)<1e-12);
assert.equal(evidenceMatch({genera:['Bacteroides']},{resolved:true,genus:'Bacteroides'},'family'),null);
assert.equal(p.mappedFeatures,236);assert.equal(p.totalFeatures,236);assert.deepEqual(p.sampleIds,d.samples.map(s=>s.id));
for(let i=0;i<5;i++)for(let j=0;j<5;j++){assert.equal(p.matrix[i][j],p.matrix[j][i]);assert.ok(p.matrix[i][j]>=0&&p.matrix[i][j]<=1);if(i===j)assert.equal(p.matrix[i][j],0)}
console.log('Passed: rank totals, child aggregation, baseline reference, coarsening distances, family evidence scope, and phylogenetic matrix integrity.');
