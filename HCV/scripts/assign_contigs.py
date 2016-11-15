#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
from collections import defaultdict
from subprocess import Popen, PIPE
from Bio import SeqIO

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"

fmt6_cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
             "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

def main(args):
    p = Popen(['blastn',
               '-db', args.ref_prefix,
               '-outfmt', '6',
               '-num_alignments', '1'
               ], 
               stdin=open(args.contigs, 'r'), stdout=PIPE, stderr=PIPE)
    o, e = p.communicate()
    with open(os.path.join(args.outdir, 'blast.hits'), 'w') as outh:
        print >>outh, o.strip('\n')
    
    # Create dictionary with hits
    hits = [l.split('\t') for l in o.strip('\n').split('\n')]
    hits = {h[0]:h for h in hits}

    bysubtype = defaultdict(list)
    for s in SeqIO.parse(args.contigs, 'fasta'):
        if s.id in hits:
            hit = dict(zip(fmt6_cols, hits[s.id]))
            bysubtype[hit["sseqid"].split('.')[0]].append(s)
        else:
            bysubtype['unassigned'].append(s)

    for subtype,seqs in bysubtype.iteritems():
        with open(os.path.join(args.outdir, '%s.contigs.fasta' % subtype), 'w') as outh:
            for s in seqs:
                print >>outh, '>%s\n%s' % (s.id, str(s.seq))

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Assign contigs to subtypes')
    parser.add_argument('contigs',
                        help='''File containing contigs, fasta format''')
    parser.add_argument('ref_prefix',
                        help='''Prefix for reference blast database''')
    parser.add_argument('outdir',
                        help='''Output directory''')
    main(parser.parse_args())
