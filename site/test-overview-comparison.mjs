import assert from 'node:assert/strict';
import {compositionShares,sharedPercentDomain} from './dist/overview-comparison.mjs';
const composition=compositionShares([21.646,.1],[.014,13.202]);
assert.deepEqual(composition[0],{bact:21.646,akk:.014,other:78.34},'other reconstructs the retained whole profile');
assert(composition.every(row=>row&&row.bact+row.akk+row.other===100),'each composition totals 100%');
assert.equal(sharedPercentDomain([[21.646,.1],[.014,13.202]]),30,'both panels share a rounded percentage domain');
assert.equal(sharedPercentDomain([[1,null],[2,3]]),null,'missing measurements are not replaced with zero');
assert.deepEqual(compositionShares([80],[21]),[null],'invalid pair values do not invent a negative remainder');
console.log('overview comparison checks passed');
