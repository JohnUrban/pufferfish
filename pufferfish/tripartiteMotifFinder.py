#!/usr/bin/env python2.7
import argparse
import itertools
import numpy as np
import pandas as pd
import time
from collections import defaultdict
from scipy.stats import norm, expon, poisson, geom
from pfmTools import *
from Bio import SeqIO

desc='''
        Designed for EcRE project 2015-2020.
        Can be used for palindromes generally.
        Can be used for any tripartite motif that has a LHS, RHS, and middle...
        ...so long as the way the middle is modeled here works.

        Requires a tab-separated file path for examples of the palindrome with the following columns:
        1 = Group / species = a string name assigning each row to a group. This is not used at the moment.
        2 = seq name = a string identifier for the sequence in this row. These need to be unique since sequences may not be.
        3 = left hand side of palindrome. All sequences need to be the same strand and the same length. Add Ns to some if need be.
        4 = middle of palindrome - sometimes there is a central bp or subseq. If there is none use "-".
        5 = right hand side of palindrome. Same rules as LHS.

        Assumes header is present: group   name    left    middle  right
        
        M = number of training examples.
        
        For the LHS (same rules for RHS):
            Seq len = L (should be >= 1)

            There are 4 states per position : A C G T
            Thus there are 4*L states overall: A1 C1 G1 T1 ... AL CL GL TL
        
            Init/Emit probs
                - prob(b_i) = freq_b_i/sum(freq_b_i for b in ACGT) = freq_b_i/M of base b at position i. (when N is encounter 0.25 added for all)
                - L x 4 pandas matrix index is range(0,L), columns are bases
            
            Trans probs:
                A dictionary of L (4 x 4) pandas matrices ACGT x ACGT
                    - dictionary is index:matrix pairs (indexes range(0,L))
                    - T1 matrix is all ones for LHS
                    - For RHS, T1 is prob of transitiong from the last base in Middle to the first base in RHS
                    - T2 - TL matrices are defined by the freq of seeing base bi after base bi-1


            p(LHS) = Product of Ei * T_i_j for i,j in range(0,L)
            
        For the middle sequence:
            Seq len = Lm
                - seqlens >= 0
                - string should be >= 1
                - "-" character represents seqlen of 0

            LHS to Middle trans probs
                = the prob of transitiong from the last base in LHS to the first base in middle
                = since the middle can be empty, the prob of transitioning to '-' is included.
            Prob of sequence length:
                = geometric prob of middle seq len given the average.
            Prob of sequence composition:
                - seq composition probs is learned by the frequencies of seeing base b in the middle sequence
                - the prob of a sequence with a given composition is the product of base probabilities

            p(Middle) = p(LHS to Mid) * p(seqlen) * p(seqcomp)

            An alternative would just be treating all seqs as 1 possible object to encounter.




    INPUTS:
    1. Motif table as described above.
    2. Scale factor for likelihood cutoff to call bipartite motif (just LHS and RHS). Generally larger than or equal to tripartite scale factor.
    3. Scale factor for likelihood cutoff to call tripartite motif.
        - 0 to Inf
        - This scales the minimum likelihood of the training data to set as the cutoff for the test data.
        
        
    
    OUTPUTS:
    Wig files:
    fixedStep  chrom=chrN
    start=1  step=1

    BED files:
    chr start end name score strand
    

    prefix.positive_strand.lhs_likelihood.wig = lhs_emit + lhs_trans

    prefix.positive_strand.rhs_likelihood.wig = rhs_emit + rhs_trans

    prefix.tripartite-motif-locations.bed (score = total likelihood)


    Algorithm:
    For each sequence given:
        1. Calculate LHS likelihood for each position.
        2. Calculate RHS likelihood for each position.
        3. For each position, check if sum of LHS + candidate RHS's based on possible mid seq lengths
            - is higher than some bipartite cutoff
            - if so check of its higher than some tripartite cutoff
            -if so, keep it.


    TODO:
    Process the complement strand.


    Make way faster.
        - store and use results for previously seen kmers instead of repeagt computations

    Iterate:
        - Take outputs, retrain, recall, etc

    Weight:
        - For training sequences, add weight that each seq should have
        - When calling spots, treat every spot as a training seq weighted by its likelihood
        - And/or when printing BED of coordinates, also include LHS/Mid/RHS seqs
            - Then BED can be converted into a training table for the next round.
    

'''

parser = argparse.ArgumentParser(description=desc,
    formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('-i', '--trainingTable',
                    type=str,
                    help='''Path to the training table as described above..
                    ''',
                    required=True)

parser.add_argument('-f', '--fasta',
                    type=str,
                    help='''Path to fasta to find motifs..
                    ''',
                    required=True)

parser.add_argument('-rc', '--revcomp',
                    action='store_true',
                    default=False,
                    help='''Compute on reverse complement of given sequence.
                    Note that this still returns Wigs and BEDs w/r/t to the original strand.
                    Thus you can compute likelihoods on the fwd strand, use this to compute
                    on the rev strand, and visualize both in separate tracks w/r.t to the fwd strand.
                    ''')

parser.add_argument('-adj', '--adjust',
                    action='store_true',
                    default=False,
                    help='''Adjusts the bipartite and tripartite wig values to be original likelihood + cutoff.
                    This results in all passing values to be >=0 and failing values to remain neagtive.
                    It is useful for viz in IGV -- easy to understand >0 and <0.
                    ''')

parser.add_argument('-c', '--bipcutoff',
                    type=str,
                    help='''cutoffs can be any number or the following strings: auto, min, max, mean, median.
                    The statistics are computed on the likelihoods of the input training set.
                    Default: auto/min''',
                    default="auto")

parser.add_argument('-C', '--tripcutoff',
                    type=str,
                    help='''cutoffs can be any number or the following strings: auto, min, max, mean, median.
                    The statistics are computed on the likelihoods of the input training set.
                    Default: auto/min''',
                    default="auto")


parser.add_argument('-p', '--prefix',
                    type=str,
                    help='''prefix for files. E.g. ecre_genome or output/ecre_genome''',
                    default="auto")

parser.add_argument('-ps', '--pseudo',
                    type=float,
                    help='''Pseudocount to avoid multiplying by 0 (or log of 0). Default 1e-2. Recommend much smaller if you want to be strict.''',
                    default=1e-2)

parser.add_argument('-sh', '--startshift',
                    type=int, default=0,
                    help='''If analyzing a sub-sequence extracted from a larger sequence,
                    give the 0-based start position within the longer sequence to report Wigs and BEDs accordingly.''',)

parser.add_argument('-nsh', '--nameshift',
                    action='store_true', default=False,
                    help='''If analyzing a sub-sequence extracted from a larger sequence,
                    this assumes the sequence has a name structure chr:start=end, and that you want to only report the "chr" component in output Wigs and BEDs.''',)

parser.add_argument('-msml', '--midSeqMinLen',
                    type=str, default='0',
                    help='''By default, this script checks all possible middle sequences from 0 bp to the max middle sequence length in the training table.
                            This tells it to only search starting from given minimum middle seq len.
                            For example, if there is only 1 middle sequence length and you dont want to look for other lengths, just provide that length.
                            Use "auto" to automatically use the shortest middle seq len in the training data.''',)


args = parser.parse_args()
                                                               



def run(args):
    # Adjust args as needed
    args.bipcutoff = float(args.bipcutoff) if args.bipcutoff not in ('auto', 'min', 'max', 'mean', 'median') else args.bipcutoff
    args.tripcutoff = float(args.tripcutoff) if args.tripcutoff not in ('auto', 'min', 'max', 'mean', 'median') else args.tripcutoff
    
    # Train
    #model = tripartiteProfileModel('/Users/johnurban/Documents/data/sciara/EcRE/Table1-ver5-April16-MikeF-or-Yutaka-EcRE-analysis.transcribed.txt', pseudo=args.pseudo)
    model = tripartiteProfileModel(args.trainingTable, pseudo=args.pseudo)
    model.train.to_csv(args.prefix + '.training-results.txt', sep="\t", index=False)
    
    # Open output files
    lhswig = open(args.prefix + '.lhs_likelihood.wig', 'w')
    rhswig = open(args.prefix + '.rhs_likelihood.wig', 'w')
    bipwig = open(args.prefix + '.bipartite_likelihood.wig', 'w')
    tripwig = open(args.prefix + '.tripartite_likelihood.wig', 'w')
    bipbed = open(args.prefix + '.bipartite.bedGraph', 'w')
    tripbed = open(args.prefix + '.tripartite.bedGraph', 'w')
    
    # Iter through seqs
    for fa in SeqIO.parse(args.fasta, 'fasta'):
        if args.nameshift:
            fa.id = str(fa.id).split(":")[0]
            
        if args.revcomp:
            seq = triPartiteSequenceSearch(sequence = str(fa.seq.reverse_complement()).upper(),
                                           name = str(fa.id),
                                           model = model,
                                           bipcutoff = args.bipcutoff,
                                           tripcutoff = args.tripcutoff)
        else:
            seq = triPartiteSequenceSearch(sequence = str(fa.seq).upper(),
                                           name = str(fa.id),
                                           model = model,
                                           bipcutoff = args.bipcutoff,
                                           tripcutoff = args.tripcutoff,
                                           minMidSeqLen = args.midSeqMinLen)

        ##LHS WIG
        seq.to_stderr("Writing LHS likelihood wig.")
        wig = seq.lhs_wig(reverse=args.revcomp,
                                  start=args.startshift+1)
        if wig:
            lhswig.write( wig + '\n' )

        ## RHS WIG
        seq.to_stderr("Writing RHS likelihood wig.")
        wig = seq.rhs_wig(reverse=args.revcomp,
                                  start=args.startshift+1)
        if wig:
            rhswig.write( wig + '\n' )

        ## BIP WIG
        seq.to_stderr("Writing bipartite likelihood wig.")
        wig = seq.bipartite_wig(adjust=args.adjust,
                                        reverse=args.revcomp,
                                        start=args.startshift+1)
        if wig:
            bipwig.write( wig + '\n' )

        ## TRIP WIG
        seq.to_stderr("Writing tripartite likelihood wig.")
        wig = seq.tripartite_wig(adjust=args.adjust,
                                          reverse=args.revcomp,
                                          start=args.startshift+1)
        if wig:
            tripwig.write( wig + '\n' )

        ## BIP BDG
        seq.to_stderr("Writing bipartite BEDGRAPH.")
        bdg = seq.bipartite_bedgraph(adjust=args.adjust,
                                        reverse=args.revcomp,
                                        start=args.startshift)
        if bdg:
            bipbed.write( bdg + '\n' )

        ## TRIP BDG
        seq.to_stderr("Writing tripartite BEDGRAPH.")
        bdg = seq.tripartite_bedgraph(adjust=args.adjust,
                                        reverse=args.revcomp,
                                        start=args.startshift)
        if bdg:
            tripbed.write( bdg + '\n' )

        seq.to_stderr("Done!")
        
    # Close files
    lhswig.close()
    rhswig.close()
    bipwig.close()
    tripwig.close()
    bipbed.close()
    tripbed.close()
        



#### EXECUTE
if __name__ == "__main__":
    run(args)

