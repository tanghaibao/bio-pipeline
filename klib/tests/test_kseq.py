import os.path as op

import unittest
import klib
import klib.kseq


DATADIR = "data"


class TestKseqMethods(unittest.TestCase):

    def setUp(self):
        self.fastafile = op.join(DATADIR, "test.fasta")
        self.reader = klib.kseq.Kseq(self.fastafile)

    def test_get_name_length_tuples(self):
        results = self.reader.get_name_length_tuples()
        self.assertEqual(results, [('HSBGPG', 1231), ('HSGLTH1', 1020)])


if __name__ == '__main__':
    unittest.main()
