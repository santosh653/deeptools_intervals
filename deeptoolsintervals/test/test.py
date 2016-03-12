#!/usr/bin/env python
from deeptoolsintervals import parse
from deeptoolsintervals.tree import printGTFtree
import os.path as op

gtf = parse.GTF()
gtf.parse(["{}/EColi.gtf.gz".format(op.dirname(op.realpath(__file__)))])
gtf.tree.printGTFtree()
