#!/usr/bin/env python2.7
import numpy as np
import pandas as pd
import argparse, sys
import datetime as dt

parser = argparse.ArgumentParser(description="""

    Given number of reads and genome size,
    Returns bin sizes for a range of expected numbers of reads per bin.

    
    """, formatter_class= argparse.RawTextHelpFormatter)



parser.add_argument('-g','--genomesize',
                   type=str, required=True, 
                   help='''Genome size integer or K, M, G format. ''')

parser.add_argument('-r', '--numreads', 
                   type=str, required=True, 
                   help='''number of reads -- integer or K, M, G format. ''')

parser.add_argument('-s','--stranded', 
                   action='store_true', default=False, 
                   help='''Whether or not reads will be counted per strand. ''')

parser.add_argument('-m','--minreads', 
                   type=int, default=10, 
                   help='''Min number of reads per bin. Always reports bin size for 1. So min other than 1. Default = 10. ''')

parser.add_argument('-M','--maxreads', 
                   type=int, default=100, 
                   help='''Max number reads per bin. Default 100. ''')

parser.add_argument('-S','--step', 
                   type=int, default=5, 
                   help='''Step size. ''')


args = parser.parse_args()



def convert(x):
    d = {'K':1e3, 'M':1e6, 'G':1e9}
    if x[-1].upper() in 'KMG':
        x = d[x[-1].upper()] * float(x[:-1])
    else:
        x = float(x)
    return x

def run(args):
    args.genomesize = convert(args.genomesize)
    args.numreads = convert(args.numreads)

    if args.stranded:
        args.genomesize = args.genomesize*2
    
    ## 1 read
    print(args.genomesize, args.numreads)
    B = args.genomesize/args.numreads
    print ('\t'.join(['expNumReads','binsize']))
    print ('\t'.join(str(e) for e in [1,round(B)]))

    ## Iter
    for i in range(args.minreads, args.maxreads+args.step, args.step):
        print( '\t'.join(
            [str(e) for e in [i, round(i*B)]]
            ))
    



#######################################################################################################
run(args)
    

