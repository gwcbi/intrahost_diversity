#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
from collections import defaultdict, Counter
from subprocess import Popen, PIPE
from Bio import SeqIO

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"

fmt6_cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
             "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

def merge_interval_list(ivs, dist=1):
    """ Merge intervals

    Args:
        ivs (list): List of intervals. Each interval is represented by a tuple of
            integers (start, end) where end > start.
        dist (int): Distance between intervals to be merged. Setting dist=1 (default) will
            merge adjacent ("bookended") intervals

    Returns:
        list: Merged list of intervals

    Examples:
        >>> merge_interval_list([])
        []
        >>> merge_interval_list([(1,10)])
        [(1, 10)]
        >>> merge_interval_list([(4, 9), (10, 14), (1, 3)])
        [(1, 3), (4, 9), (10, 14)]
        >>> merge_interval_list([(4, 9), (10, 14), (1, 3)], dist=1)
        [(1, 14)]
    """
    if len(ivs)<= 1: return ivs
    ivs.sort(key=lambda x:x[0])
    ret = [ivs[0]]
    for iv in ivs[1:]:
        if iv[0] - ret[-1][1] > dist:
            ret.append(iv)
        else:
           ret[-1] = (ret[-1][0], max(iv[1],ret[-1][1]))
    return ret

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
    
    # Parse hits
    hits = [l.split('\t') for l in o.strip('\n').split('\n')]
    
    # Find intervals
    intervals = defaultdict(list)
    for h in hits:
        sstart, send = (int(h[8]), int(h[9])) if int(h[8]) < int(h[9]) else (int(h[9]), int(h[8]))
        intervals[h[1]] = merge_interval_list(intervals[h[1]] + [(sstart, send)])

    # Determine isolate with best mapping
    subisolates = defaultdict(Counter)
    for ref, ivs in intervals.iteritems():
        subisolates[ref.split('.')[0]][ref] = sum([t[1]-t[0] for t in ivs])
        
    with open(os.path.join(args.outdir, 'subtypes.config'), 'w') as outc:
        for subtype in sorted(subisolates.keys()):
            mc = subisolates[subtype].most_common()
            if mc[0][1] >= args.cutoff:
                print >>sys.stdout, '--> Subtype: %s' % subtype
                print >>sys.stdout, '--> \t%s\t%s' % mc[0]
                print >>sys.stdout, '\n'.join('[X]\t%s\t%s' % _ for _ in mc[1:])
                print >>outc, "%s\t%s" % (subtype, mc[0][0])
            else:
                print >>sys.stdout, '[X] Subtype: %s' % subtype
                print >>sys.stdout, '\n'.join('[X]\t%s\t%s' % _ for _ in mc)
    
    # Create dictionary with hits
    hits = {h[0]:h for h in hits}

    bysubtype = defaultdict(list)
    for s in SeqIO.parse(args.contigs, 'fasta'):
        if s.id in hits:
            hit = dict(zip(fmt6_cols, hits[s.id]))
            bysubtype[hit["sseqid"].split('.')[0]].append(s)
        else:
            bysubtype['unassigned'].append(s)

    for subtype,seqs in bysubtype.iteritems():
        with open(os.path.join(args.outdir, '%s.fasta' % subtype), 'w') as outh:
            for s in seqs:
                print >>outh, '>%s\n%s' % (s.id, str(s.seq))

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Assign contigs to subtypes')
    parser.add_argument('--cutoff', default=0, type=int,
                        help='''Minimum coverage to output sequences''')
    parser.add_argument('contigs',
                        help='''File containing contigs, fasta format''')
    parser.add_argument('ref_prefix',
                        help='''Prefix for reference blast database''')
    parser.add_argument('outdir',
                        help='''Output directory''')
    main(parser.parse_args())
