export const sharedPercentDomain=series=>{
 const values=series.flat();
 return values.every(Number.isFinite)?Math.max(10,Math.ceil(Math.max(...values)/10)*10):null;
};
export const compositionShares=(bact,akk)=>bact.map((value,index)=>{
 const other=100-value-akk[index];
 return Number.isFinite(value)&&Number.isFinite(akk[index])&&other>=0?{bact:value,akk:akk[index],other}:null;
});
