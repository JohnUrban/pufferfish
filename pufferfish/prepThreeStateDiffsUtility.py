#!/usr/bin/env python2.7
import numpy as np
import pandas as pd
import argparse, sys
import datetime as dt

parser = argparse.ArgumentParser(description="""

    Usage case - chromatin states on two different cell types.

    E.g. Nurse cells (NCs) vs Germline Stem Cells (GSCs)

    For each, prepare/process the data identically.
    - same mapping program and parameters
    - same filtering
    - same windows/bins for counting
    - same normalization procedure (e.g. eliminate 0 bins, median ratio norm --protocol12)
    - Use puffcn HMM to segment into 3-states (1,2,3) representing depleted, no difference, and enriched.
    - Use BEDtools as follows to pair identical bins in both samples with states from both samples
        :: bedtools intersect -wo -a GSC_Q10_iter22.bedGraph -b St5NC_Q10_iter43.bedGraph | awk '$1==$5 && $2==$6 && $3==$7 && $4>0 && $8>0 {OFS="\t"; print $1,$2,$3,$4,$8}' > GSC-NC-states.tab
    - That returns a 5-column table -- like a bedGraph with two score columns.
    - Give that as input to this utility script.
    - This script:
        - Converts states 1,2,3 to -1,0,1.
        - Subtracts column4 from column5: difference = C5-C4
        - This converts the state pairs into 3 states: -2,-1,0,1,2
        - It then converts -2,-1,0,1,2 to -1,0,1 to to 1,2,3
        - It returns bins as bedGraph with the new 1,2,3 states
        - This in turn can be input to puffcn and segmented using --emodel discrete
            - and setting emission priors as something like :
                :: --mu "0.98,0.01,0.01;0.01,0.98,0.01;0.01,0.01,0.98"

    - 
    # strategy
    HMM models 3 3-sided dice (or 5 5-sided) with different biases to model dependencies and call larger domains of "no change", "change up", "change down". For the latter, the states from the former (1,2,3) would be converted to (-1,0,1). Then we could get the state differences as 3 states (-2, 0, 2) as shown below:

    # state_gsc     state_nc      difference(nc-gsc)             interpretation         collapsed_interp
    # -1             -1             0                              no change             no change
    # -1              1             2                              2 change up           up
    # -1              0             1                              1 change up           up
    #  1             -1            -2                              2 change down         down
    #  1              1             0                              no change             no change
    #  1              0            -1                              1 change down         down
    #  0              0             0                              no change             no change
    #  0              1             1                              1 change up           up
    #  0             -1            -1                              1 change down         down

    ## Collapse -1 and -2 to -1; 1 and 2 to 1


    # After getting the difference (-1,0,1), I'd use an HMM and automated learning to get the larger domains. 
    # Each domain type would be modeled as a biased 3-sided dice that can turn up -2, 0, 2 ... or 1,2,3 (doesn't matter), 
    #	preferentially turning up only one of those. It is a very simple approach and is similar to how CpG island domains are called. 
    # In that example, there are two dice that turn up A,C,G,T or (dimers) with different biases.
        
    """, formatter_class= argparse.RawTextHelpFormatter)



parser.add_argument('-t','--table',
                   type=str, required=True, 
                   help='''Path to input table.''')


args = parser.parse_args()

def run(args):
    ## VARIABLES/CONSTANTS
    expected_states = (1,2,3)
    in_convert = {1:-1, 2:0, 3:1}
    mid_convert = {-2:-1, -1:-1, 0:0, 1:1, 2:1}
    out_convert = {-1:1, 0:2, 1:3}

    ## EXECUTE
    with open(args.table) as tab:
        for line in tab:
            line = line.strip().split()
            A = int(float(line[3])) 
            B = int(float(line[4]))
            try:
                assert A in expected_states and B in expected_states
            except AssertionError:
                sys.stderr.write("Encountered unexpected value in table: A:" + str(A) + ", B:" + str(B) + "\n")
                quit()
            Ai = in_convert[A]
            Bi = in_convert[B]
            DIFF = Bi-Ai
            DIFFm = mid_convert[DIFF]
            DIFFo = out_convert[DIFFm] ## could of course do mid and out in 1 step
            ## out = line[:3] + [A, B, Ai, Bi, DIFF, DIFFm, DIFFo]
            out = line[:3] + [DIFFo]
            out = '\t'.join(str(e) for e in out)
            print( out )



#########################################################################


run(args)

        
        
