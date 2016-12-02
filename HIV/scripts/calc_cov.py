#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import re
import math
import gzip

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"

vcf = sys.argv[1]
pct = float(sys.argv[2])

if pct > 1:
    print >>sys.stdout, '%d' % int(pct)
    sys.exit()

lines = (l.strip('\n').split('\t') for l in gzip.open(vcf, 'rb') if not l.startswith('#'))
covs = []
for l in lines:
    m = re.search('DP=(\d+);', l[7])
    if m: covs.append(int(m.group(1)))

covs.sort()
med = covs[len(covs)/2]
print >>sys.stdout, '%d' % (med * pct)
# print >>sys.stdout, '%d' % covs[int(math.ceil(len(covs) * pct))]
