#! /usr/bin/env python

import sys
from centr_utils.utils import readfq


input_file = sys.argv[1]

seq_db = {}
with open(input_file, 'r') as hin:
    for name, seq, qual in readfq(hin):
        seq_db[name] = seq

for name in sorted(seq_db):
    seq_ext = seq_db[name] + seq_db[name][:170]
    print(">%s\n%s" % (name, seq_ext))

