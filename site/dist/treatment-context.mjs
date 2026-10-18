const endOf=e=>e.end_date||e.date;

// A point is only in the interval after the previous reported date. A reported
// period can overlap either boundary, because its daily dosing is not known.
export const eventsForInterval=(events,previousDate,selectedDate)=>events.filter(event=>
 event.end_date
  ? event.date<=selectedDate&&endOf(event)>=previousDate
  : event.date>previousDate&&event.date<=selectedDate
);

export const ongoingAt=(events,date)=>events.filter(event=>event.end_date&&event.date<=date&&endOf(event)>=date);

export function largestObservedChange(values){
 const series=Array.isArray(values[0])?values:[values];
 let largest=null;
 series.forEach((row,seriesIndex)=>row.forEach((value,index)=>{
  if(index===0||!Number.isFinite(value)||!Number.isFinite(row[index-1]))return;
  const delta=value-row[index-1];
  if(delta!==0&&(!largest||Math.abs(delta)>Math.abs(largest.delta)))largest={index,series:seriesIndex,delta};
 }));
 return largest;
}

export const treatmentUncertainty=(previous,selected)=>{
 const notes=[];
 if(previous?.alternateDate||selected?.alternateDate)notes.push('February has conflicting Feb 13 and Feb 19 boundaries.');
 if(previous?.dateStatus==='received'||selected?.dateStatus==='received')notes.push('March dates are lab receipt dates, so receipt timing cannot confirm collection ordering.');
 return notes.join(' ');
};

export function relevantTreatments(events,evidence,active){
 const preferred=(active==='bifido'||active==='akk')
  ? ['Metformin','Cefazolin','Dalbavancin','Wegovy','Zepbound']
  : ['Cefazolin','Dalbavancin','Metformin','Wegovy','Zepbound'];
 const titles=new Set(events.map(event=>event.title));
 return preferred.filter(title=>titles.has(title)&&evidence[title]).slice(0,2).map(title=>({title,...evidence[title]}));
}
