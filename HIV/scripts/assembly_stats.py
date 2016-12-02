#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"

def N50(l):
    slen = sorted(l)
    cumsum = [sum(slen[:i+1]) for i in range(len(slen))]
    i = len([c for c in cumsum if c < sum(slen)*0.5])
    return len(slen)-i, slen[i]

contig_lengths = sorted(len(s) for s in SeqIO.parse(sys.stdin, 'fasta'))
print >>sys.stdout, 'num_contigs\t%d' % len(contig_lengths)
print >>sys.stdout, 'max_contig\t%d' % contig_lengths[-1]
print >>sys.stdout, 'contig>1kb\t%d' % sum(l>=1000 for l in contig_lengths)
l50,n50 = N50(contig_lengths)
print >>sys.stdout, 'contig_N50\t%d' % n50
print >>sys.stdout, 'contig_L50\t%d' % l50
contig_lengths.sort(reverse=True)
for i,l in enumerate(contig_lengths[1:5]):
    print >>sys.stdout, 'rank%d\t%d' % (i+2, l)
