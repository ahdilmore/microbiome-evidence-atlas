import assert from 'node:assert/strict';
import {eventsForInterval,largestObservedChange,relevantTreatments,treatmentUncertainty} from './dist/treatment-context.mjs';

const events=[
 {title:'past',date:'2026-01-01',end_date:'2026-01-04'},
 {title:'period',date:'2026-01-05',end_date:'2026-01-15'},
 {title:'boundary point',date:'2026-01-10',end_date:''},
 {title:'previous point',date:'2026-01-05',end_date:''},
 {title:'future',date:'2026-01-11',end_date:''}
];
assert.deepEqual(eventsForInterval(events,'2026-01-05','2026-01-10').map(e=>e.title),['period','boundary point'],'periods overlap the interval; points are after the previous date and include the endpoint');
assert.equal(largestObservedChange([1,1,1]),null,'flat series has no jump target');
assert.deepEqual(largestObservedChange([[1,3,2],[10,4,5]]),{index:1,series:1,delta:-6},'multi-series comparison finds the largest absolute consecutive change');
assert.match(treatmentUncertainty({alternateDate:'2026-02-19'},{dateStatus:'received'}),/February.*receipt timing/,'both documented date ambiguities are retained');
const evidence={Metformin:{explanation:'m',source:{}},Cefazolin:{explanation:'c',source:{}},Dalbavancin:{explanation:'d',source:{}}};
assert.deepEqual(relevantTreatments([{title:'Cefazolin'},{title:'Metformin'},{title:'Dalbavancin'}],evidence,'akk').map(x=>x.title),['Metformin','Cefazolin'],'Akkermansia ordering prefers Metformin');
console.log('treatment context checks passed');
