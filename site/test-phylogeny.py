from io import StringIO
from Bio import Phylo
from prepare_phylogeny import weighted_matrix
# Equal-height tree: sister replacement is smaller than replacement across the root.
t=Phylo.read(StringIO('((a:1,b:1):2,c:3);'),'newick')
m,_=weighted_matrix(t,['a','b','c'],[[1,0,0],[0,1,0],[0,0,1]])
assert abs(m[0,1]-1/3)<1e-12
assert abs(m[0,2]-1)<1e-12
assert all(m[i,i]==0 for i in range(3))
try: weighted_matrix(t,['missing'],[[1]])
except AssertionError: pass
else: raise AssertionError('Missing tips must not be dropped silently')
print('Passed: hand-calculated weighted UniFrac, sister-vs-distant replacement, identity and missing-tip rejection.')
