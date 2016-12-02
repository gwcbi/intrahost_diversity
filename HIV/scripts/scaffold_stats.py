#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"

for s in SeqIO.parse(sys.stdin, 'fasta'):
    subtype = s.id.split('.')[0]
    total_length = len(s)
    num_call = total_length - sum(_ not in 'ACGT' for _ in str(s.seq).upper())
    pct_call = float(num_call) / total_length
    print >>sys.stdout, '%s\t%d\t%d\t%.1f' % (subtype, total_length, num_call, pct_call*100)
