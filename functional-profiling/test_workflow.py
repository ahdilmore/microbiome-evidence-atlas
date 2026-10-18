import gzip
from pathlib import Path
import tempfile
import unittest
from run import concatenate_fastq
from summarize import summarize, read_table

class Checks(unittest.TestCase):
    def test_mates_get_unique_ids_and_bad_fastq_fails(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d); files=[d/'r1.gz',d/'r2.gz']
            for f in files:
                with gzip.open(f,'wt') as g: g.write('@same\nACGT\n+\nIIII\n')
            self.assertEqual(concatenate_fastq(files,d/'combined.gz'),2)
            with gzip.open(d/'combined.gz','rt') as f:
                content=f.read(); self.assertIn('@read_1',content); self.assertIn('@read_2',content)
            with gzip.open(files[0],'wt') as g: g.write('@bad\nACGT\n+\nI\n')
            with self.assertRaises(ValueError): concatenate_fastq(files,d/'bad.gz')

    def test_special_and_stratified_rows_not_double_counted(self):
        tables=[{'PWY-A':20,'PWY-A|taxon':20,'UNMAPPED':900000,'PWY-B':30},
                {'PWY-A':22,'PWY-A|taxon':22,'UNMAPPED':800000}]
        out=summarize([{'id':'a'},{'id':'b'}],tables,10)
        self.assertEqual(len(out['pathways']),2)
        self.assertEqual(out['persistent_least_variable'][0]['pathway'],'PWY-A')
        self.assertEqual(out['most_changed'][0]['pathway'],'PWY-B')
        self.assertEqual(out['most_changed'][0]['adjacent_deltas_cpm'],[-30])
        self.assertEqual(len(summarize([{},{}],[{'A':0},{'A':0}],0)['persistent_least_variable']),0)

    def test_nonfinite_values_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'table.tsv';p.write_text('# pathway\tsample\nA\tnan\n')
            with self.assertRaises(ValueError): read_table(p)

if __name__=='__main__': unittest.main()
