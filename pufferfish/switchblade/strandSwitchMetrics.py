######## DEPRECATED
######## LIVE VERSION @ https://github.com/JohnUrban/switchblade

#!/usr/bin/env python2.7

import sys
from collections import defaultdict, deque
from scipy.stats.stats import pearsonr
from math import log
import argparse


## THIS PARTICULAR SCRIPT IS DEPRECATED
## It is now under the name strandSwitchMetrics.py in PufferFish suite.

## Created circa June 11, 2014.
##      Former scriptname = rawPileUpToOEMinput.py
##      - original use: convert bedtools stranded depth files into input format for OEM.py.
##      - evolved use : instead of needing 2 scripts and 2 steps, just perform OEM steps here as well.
##      old examples
##      ./rawPileUpToOEMinput.py loessAsOEMinput/dbf4loessSpan0.3.imputedMissingData.txt 300 0.000000001 > callPeaksFromOEMCorrectedHeightSignal/ddbf4loessSpan0.3.imputedMissingData.300bpWin.1e-9pseudo.oemFile
## September 2021
##      - added argument parser, and re-structured accordingly.
##      - now can use scriptname.py -h for more usage info.
##      - legacy usage still possible with input format 2 and output format 1 (see help options for more info on formats).



##############################################################################
''' ARGUMENT PARSER '''
##############################################################################

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    

    Usage:


          

              """, formatter_class= argparse.RawTextHelpFormatter)

parser_input = parser.add_mutually_exclusive_group()
parser_input.add_argument('-i', "--input",
                   type= str, default="-",
                   help='''Text file (or stdin) of variety of formats to be specified with -f. Use "-" or "stdin" or leave blank for piping into stdin.''')

parser.add_argument('-f', "--inputformat",
                   type=int, required=True,
                   help='''What input format to expect?
1   =   CSV: Complete position records: Chr, pos, fwd strand coverage, reverse strand coverage. (pos is 1 -based)
2   =   TSV: Stranded position records format1: Chr, pos, strand, coverage.
            Forward strand can be specified as 0, +, plus, f, fwd, forward, pos, positive.
            Reverse strand can be specified as 16, -, minus, r, rev, reverse, neg, negative.
3   =   CSV: Loess smoothed position records: chr, pos, smoothed value*.
        The smoothed value is obtained first by getting counts from fwd and reverse strands independently.
        Then counts on negative strand are moved as needed (e.g. to represent 5' end) and made negative (e.g. 10 becomes -10).
        Then positive and negative counts are combined. Then smoothing is performed to give smoothed averages over each position weighted toward the strand with most coverage.
        Thus, smoothed values range from negative to positive real values.''')

parser.add_argument('-F', "--outputformat",
                   type=int, default=2,
                   help='''What input format to expect? Default is format 2.
1   =   CSV: OEM input as output (same as input format 1): Chr, pos, fwd strand coverage, reverse strand coverage. (pos is 1 -based)
        Only possible with input format 2.
2   =   TSV: OEM expanded output: Skip the 2 step process that pipes into OEM.py and get OEM output here.''')


parser.add_argument('-w', '--windowsize', type=int, default=10000,
                    help='''Window size to use. Default = 10000 (10 kb). windowSize is how big the window is to left of a position (for WL and CL calculations) and to the right (for WR and CR). This should be set to '10000' if unsure.''')

parser.add_argument('-p', '--pseudocount', type=float, default=0,
                    help='''Pseudocount to use. Default = 0. The pseudocount is used to prevent division by 0.
              In 10kb bins on yeast genome it is unlikely to get a division-by-zero error.
              This is not necessarily true for the human genome or smaller window sizes, thus a pseudocount can help.
              Suggested pseudocount if needed is 1, but using sequencing depth and other knowledge to decide is better.''')

parser.add_argument('-a', '--negStrandPositionAdjustment', type=int, default=0,
                    help='''In some applications, one needs to adjust the reported positions on the negative strands.
For example, if you create an input file (format 1) from SAM records (Chr, pos, strand, coverage), and you want to assign the counts to 5' ends, you might adjust negative strand position by adding the readlength -1.
This is not perfect. It is advised to do all adjustments, including this type, BEFORE using this script with utilities such as AWK.
For the case given, BEDtools might be all you need to get actual coverage values of 5 prime ends on each strand.''')
 
#parser.add_argument('-H', '--hasheader', action='store_true', default=False,
#                    help='''Use this flag if the inputCSV has the header line 'chr,pos,fwd_str,rev_str'.''')

#parser.add_argument('-c', '--csvformat', action='store_true', default=False,
#                    help='''Use this flag for CSV output formatting instead of bedGraph formatting.''')



args = parser.parse_args()

##############################################################################
''' ASSERTIONS '''
##############################################################################
if args.outputformat == 1:
    assert args.inputformat == 2


##############################################################################
''' FUNCTIONS '''
##############################################################################



def log_b(b):
    ''' returns function that does log_b(x)'''
    logOfB = log(b)
    def logb(x):
        return log(x)/logOfB
    return logb
    
log2 = log_b(b=2)


class PositionCounts(object):
    def __init__(self, negStrandPositionAdjustment=0, windowsize=500, pseudocount=0.1):
        self.posdict = dict()
        self.negStrandPositionAdjustment = negStrandPositionAdjustment
        self.chrom = None
        self.windowsize = windowsize
        self.pseudocount = pseudocount
    def sameChrom(self, chrom):
        if self.chrom == None:
            self.chrom = chrom
            return True
        elif self.chrom != chrom:
            return False
        else: #self.chrom == chrom
            return True
    def changeChrom(self, chrom):
        self.chrom = chrom
            
##    def addCount(self, pos, strand, count):
##        '''0 = fwd strand; 16 = rev strand. Integers match SAM flags - otherwise arbitrary.'''
##        if strand == 16:
##            #pos += self.readlen - 1
##            pos += self.negStrandPositionAdjustment
##        try:
##            self.posdict[pos][strand] += float(count)
##        except KeyError:
##            try:
##                self.posdict[pos][strand] = float(count)
##            except:
##                self.posdict[pos] = {0:0,16:0}
##                self.posdict[pos][strand] = float(count)
    def addCount(self, chrom, pos, strand, count):
        '''0 = fwd strand; 16 = rev strand. Integers match SAM flags - otherwise arbitrary.'''
            
        if strand == 16:
            #pos += self.readlen - 1
            pos += self.negStrandPositionAdjustment

        # Ensure/initialize chrom key
        try:
            self.posdict[chrom]
        except:
            self.posdict[chrom] = {}

        # Ensure/initialize chrom pos key
        try:
            self.posdict[chrom][pos]
        except:
            self.posdict[chrom][pos] = {0:0,16:0}

        # Add count to strand
        self.posdict[chrom][pos][strand] += float(count)

            
##    def outputChrom(self, newchrom):
##        for pos in self.posdict.keys():
##            fwd = self.posdict[pos][0]
##            rev = self.posdict[pos][16]
##            print ",".join(str(e) for e in [self.chrom,pos,fwd,rev])

    def outputChrom(self, chrom=None):
        if chrom is None:
            chrom=self.chrom
        for pos in sorted(self.posdict[chrom].keys()):
            fwd = self.posdict[chrom][pos][0]
            rev = self.posdict[chrom][pos][16]
            print ",".join(str(e) for e in [chrom,pos,fwd,rev])

    def outputAllChrom(self):
        for chrom in sorted(self.posdict.keys()):
            self.outputChrom(chrom)
            
    def outputChromOEM(self, chrom=None):
        if chrom is None:
            chrom=self.chrom
        OEM(chrom, self.posdict[chrom], self.windowsize, self.pseudocount)
        #self.chrom = newchrom
        #self.posdict = dict()

    def outputAllChromOEM(self):
        for chrom in sorted(self.posdict.keys()):
            self.outputChromOEM(chrom)
        

    def outputChromLocalCorrelation(self,newchrom):
        localCor(self.chrom)

    def reset(self, newchrom):
        self.chrom = newchrom
        self.posdict = dict()        


def updateWC(Pos, posdict, WL, CL, WR, CR, windowsize, pseudocount,fwd=0,rev=16):
    ''' part of OEM below'''
    try:
        WL += posdict[Pos][fwd] + pseudocount
        CL += posdict[Pos][rev] + pseudocount
        WR -= posdict[Pos][fwd]
        CR -= posdict[Pos][rev]
    except KeyError:
        WL += pseudocount
        CL += pseudocount
    try:
        WL -= posdict[Pos-windowsize][fwd]
        CL -= posdict[Pos-windowsize][rev]
    except KeyError:
        pass
    try:
        WR += posdict[Pos+windowsize][fwd] + pseudocount
        CR += posdict[Pos+windowsize][rev] + pseudocount
    except KeyError:
        WR += pseudocount
        CR += pseudocount
    return WL, CL, WR, CR


def OEM(chrom, posdict, windowsize=500, pseudocount=0.1):
    """Takes in posDict (of single chr), returns oemDict"""
    def findbalance(WL, CR, log2Dict):
        try: log2Dict[WL] ## if it sets of KeyError put it in
        except KeyError: log2Dict[WL] = log2(WL)
        try: log2Dict[CR]
        except KeyError: log2Dict[CR] = log2(CR)
        return log2Dict[WL] - log2Dict[CR], log2Dict
    fwd=0
    rev=16
    # initialize oemDict
    oemDict = {}
    positions = deque(sorted(posdict.keys()))
    chrStart = positions[0]
    chrLen = positions[-1]
    log2Dict = {2:log2(2)}
    # Set the start and end position (end position is 10,000 bp less than last position) (JB)
    start = chrStart+windowsize-1
    ##end = chrLen - windowsize ## make sure this works
    end = chrLen - windowsize + 1 ## make sure this works

    ## GIVE DUMMY LINE AND STOP IF CHRLEN < MinLen == (w+1)*2
    if chrLen < (windowsize+1)*2:
        Pos = 1
        oem, ratio, height, balance, bias, skew, summit = 0, 0, 0, 0, 0, 0, 0
        sys.stdout.write("\t".join([str(e) for e in [chrom, Pos, oem, ratio, height, balance, bias, skew, summit]])+"\n")
        return

    ### CONTINUE IF CHRLEN >= MinLen == (w+1)*2
    # Calculate the first window (JB) -- float() used to avoid default 'int division'
    WL, CL, WR, CR = pseudocount, pseudocount, pseudocount, pseudocount
    L = [chrStart,start+1]
    R = [start+1,start+1+windowsize]
    while positions[0] >= L[0] and positions[0] < L[1]:
        WL += posdict[positions[0]][fwd]
        CL += posdict[positions[0]][rev]
        positions.popleft()
    while positions[0] >= R[0] and positions[0] < R[1]:
        WR += posdict[positions[0]][fwd]
        CR += posdict[positions[0]][rev]
        positions.popleft()

    # Set up the storage oem (JB)
    L = (WL / (WL + CL))
    R = (WR / (WR + CR))
##    print L, R
    oem = L - R
    ratio = L/R
    height = sum([WL,CL,WR,CR])
##    oemCorrectedHeight = oem*height ## really oem can be used to correct FE signal not just treatment read signal
    balance, log2Dict = findbalance(WL, CR, log2Dict)  ## divind by max(balance) in R puts everything between 1 to -1
    bias = (WL-CR)/(WL+CR)
    ## bias is related to balance -- +1 means the balance is skewed toward left peak, -1 means skewed towad right peak, 0 is noskew (balance) what you want
    ## you can have no-skew and balance in non-double peak situations though
##    biasInProportion = bias*(WL+CR) ## this gives a bias proportional to height of signal looked at
    ##  -- with shorter double peaks, smaller fluxuations look more biased so multiplying by the smaller number sends it less far away from 0 than if they were taller double peaks
    skew = (WL-CL)/(WL+CL) - (WR-CR)/(WR+CR)
    ## skew is related to oem. it can go from +2 to -2.
    summit = 2 * (WL * CR)**0.5 - WR - CL
    #oem = (WL / (WL + CL)) - (WR / (WR + CR))
    sys.stdout.write("\t".join([str(e) for e in [chrom, start, oem, ratio, height, balance, bias, skew, summit]])+"\n")

    ## Iterative steps
    numPositions = len(range(start,end))
    Pos = start
    for i in range(1, numPositions):
        # Update the pos
        Pos += 1
        WL, CL, WR, CR = updateWC(Pos, posdict, WL, CL, WR, CR, windowsize, pseudocount)
        L = (WL / (WL + CL))
        R = (WR / (WR + CR))
        oem = L - R
##        heightWC = WL - CR
        ratio = L/R
        height = sum([WL,CL,WR,CR])
##        oemCorrectedHeight = oem*height ## can correct height in R by oem*height or oem*skew
        balance, log2Dict = findbalance(WL, CR, log2Dict)
        bias = (WL-CR)/(WL+CR)
##        biasInProportion = bias*(WL+CR)
        skew = (WL-CL)/(WL+CL) - (WR-CR)/(WR+CR)
        summit = 2 * (WL * CR)**0.5 - WR - CL
            ## in practice, for positive values macs2RefinePeak(summit) is proportional to oemCorrectedHeight
            ## and for negative values just follows the of the loess smoothed curve whether the values were negative or, if positive, reflected over 0
            ## in some ways its is just reflecting height over 0 so all height is negative and"dips/valleys" become peaks (due to the reflection), some of which summit over 0
            ## oemCorrectedHeight is a little better for the following reasons:
            ##  1. it sees both WL-->CR switches (positive peaks, oris) and CL --> WR (neg peaks, fork merge spots) as oem is able to do
            ##  2. everything else is essentially 0 so it just looks like +/- spikes occurring along a line
        
        sys.stdout.write("\t".join([str(e) for e in [chrom, Pos, oem, ratio, height, balance, bias, skew, summit]])+"\n")

### What needs to be done: all that should be calculated here is oem
###     Then there needs to be a function that goes through oem signal and marks peak boundaries give some cutoff.
###         Could have separate ones that call pos OEMs and neg OEMs
###     Then something would take peak boundaries in post-hoc and calculate these metrics (height, ormCorrHeight, balance, bias, biasInPro, skew)
###     One could do it across the peak in sliding windows, or just from the peak OEM summit to left and right some fixed distance (so result is not influenced by peak width)
###     Below I pasted macs2 refine peak:
##        It starts at peak start and takes WL,WR,CL,CR in 100 bp windows to each side
##            It applies/records this calc: 2 * (watson_left * crick_right)**0.5 - watson_right - crick_left
##        It then iterates moving up 1 bp and does it all again up intil peak end
##        It then takes the index of the max score (from calculation) as the summit
##        Thus, this calculation behaves in some proportional way to the OEM calc

##def macs2_find_summit(chrom, plus, minus, peak_start, peak_end, name = "peak", window_size=100, cutoff = 5):
##    
##    left_sum = lambda strand, pos, width = window_size: sum([strand[x] for x in strand if x <= pos and x >= pos - width])
##    right_sum = lambda strand, pos, width = window_size: sum([strand[x] for x in strand if x >= pos and x <= pos + width])
##    left_forward = lambda strand, pos: strand.get(pos,0) - strand.get(pos-window_size, 0)
##    right_forward = lambda strand, pos: strand.get(pos + window_size, 0) - strand.get(pos, 0)
##
##    watson, crick = (Counter(plus), Counter(minus))
##    watson_left = left_sum(watson, peak_start)
##    crick_left = left_sum(crick, peak_start)
##    watson_right = right_sum(watson, peak_start)
##    crick_right = right_sum(crick, peak_start)
##
##    wtd_list = []
##    for j in range(peak_start, peak_end+1):
##        wtd_list.append(2 * (watson_left * crick_right)**0.5 - watson_right - crick_left)
##        watson_left += left_forward(watson, j)
##        watson_right += right_forward(watson, j)
##        crick_left += left_forward(crick, j)
##        crick_right += right_forward(crick,j)
##
##    wtd_max_val = max(wtd_list)
##    wtd_max_pos = wtd_list.index(wtd_max_val) + peak_start
##
##    #return (chrom, wtd_max_pos, wtd_max_pos+1, wtd_max_val)
##
##    if wtd_max_val > cutoff:
##        return (chrom, wtd_max_pos, wtd_max_pos+1, name+"_R" , wtd_max_val) # 'R'efined
##    else:
##        return (chrom, wtd_max_pos, wtd_max_pos+1, name+"_F" , wtd_max_val) # 'F'ailed

def localCorr():
    
    pass

def getstrand(strand):
    if strand in ('0', '+', 'plus', 'fwd', 'f', 'forward', 'pos', 'positive'): 
        return 0
    elif strand in ('16', '-', 'minus', 'rev', 'r', 'reverse', 'neg', 'negative'):
        return 16
    else:
        sys.stderr("Unexpected strand format encountered.... Exiting....")
        quit()
        
##def processCounts(fileconnection,
##                  negStrandPositionAdjustment,
##                  windowsize=500,
##                  pseudocount=0.1,
##                  delim="\t"):
##    ## initial
##    countsOf5pEnds = PositionCounts(negStrandPositionAdjustment=negStrandPositionAdjustment,
##                                    windowsize=windowsize,
##                                    pseudocount=pseudocount)
##    
##    for line in fileconnection:
##        line = line.strip().split( delim )
##        chrom = line[0]
##        pos = int(line[1])
##        strand = getstrand(line[2])
##        count = int(line[3])
##        countsOf5pEnds.addCount(chrom=chrom,
##                                pos=pos,
##                                strand=strand,
##                                count=count)
####        if countsOf5pEnds.sameChrom(chrom):
####            countsOf5pEnds.addCount(chrom=chrom,
####                                    pos=pos,
####                                    strand=strand,
####                                    count=count)
####        else:
####            countsOf5pEnds.outputChrom()
####            countsOf5pEnds.changeChrom(chrom)
####            countsOf5pEnds.addCount(chrom=chrom,
####                                    pos=pos,
####                                    strand=strand,
####                                    count=count)
####    ## print last chr
####    countsOf5pEnds.outputChrom()
##    countsOf5pEnds.outputAllChrom()


##def processOEMs(fileconnection,
##                negStrandPositionAdjustment,
##                windowsize=500,
##                pseudocount=0.1):
##    ## initial
##    countsOf5pEnds = PositionCounts(negStrandPositionAdjustment=negStrandPositionAdjustment,
##                                    windowsize=windowsize,
##                                    pseudocount=pseudocount)
##    for line in fileconnection:
##        line = line[:-1].split()
##        chrom = line[0]
##        pos = int(line[1])
##        strand = getstrand(line[2])
##        count = int(line[3])
##        countsOf5pEnds.addCount(chrom=chrom,
##                                pos=pos,
##                                strand=strand,
##                                count=count)
####        if countsOf5pEnds.sameChrom(chrom):
####            countsOf5pEnds.addCount(chrom=chrom,
####                                    pos=pos,
####                                    strand=strand,
####                                    count=count)
####        else:
####            countsOf5pEnds.outputChromOEM()
####            countsOf5pEnds.changeChrom(chrom)
####            countsOf5pEnds.addCount(chrom=chrom,
####                                    pos=pos,
####                                    strand=strand,
####                                    count=count)
####    ## print last chr
####    countsOf5pEnds.outputChromOEM()
##    countsOf5pEnds.outputAllChromOEM()

    

def processCounts(fileconnection,
                  negStrandPositionAdjustment,
                  windowsize=500,
                  pseudocount=0.1,
                  delim="\t"):
    ## initial
    countsOf5pEnds = PositionCounts(negStrandPositionAdjustment=negStrandPositionAdjustment,
                                    windowsize=windowsize,
                                    pseudocount=pseudocount)
    
    for line in fileconnection:
        line = line.strip().split( delim )
        chrom = line[0]
        pos = int(line[1])
        strand = getstrand(line[2])
        count = int(line[3])
        countsOf5pEnds.addCount(chrom=chrom,
                                pos=pos,
                                strand=strand,
                                count=count)
    return countsOf5pEnds

def printOEMstyleInput(fileconnection,
                  negStrandPositionAdjustment,
                  windowsize=500,
                  pseudocount=0.1,
                  delim="\t"):
    counts = processCounts(fileconnection,
                  negStrandPositionAdjustment,
                  windowsize,
                  pseudocount,
                  delim)
    counts.outputAllChrom()


         
##def processOEMs(fileconnection,
##                negStrandPositionAdjustment,
##                windowsize=500,
##                pseudocount=0.1,
##                  delim="\t"):
##    ## initial
##    countsOf5pEnds = PositionCounts(negStrandPositionAdjustment=negStrandPositionAdjustment,
##                                    windowsize=windowsize,
##                                    pseudocount=pseudocount)
##    for line in fileconnection:
##        line = line.strip().split( delim )
##        chrom = line[0]
##        pos = int(line[1])
##        strand = getstrand(line[2])
##        count = int(line[3])
##        countsOf5pEnds.addCount(chrom=chrom,
##                                pos=pos,
##                                strand=strand,
##                                count=count)
##
##    countsOf5pEnds.outputAllChromOEM()

def printStrandSwitchMetrics(fileconnection,
                negStrandPositionAdjustment,
                windowsize=500,
                pseudocount=0.1,
                  delim="\t"):
    counts = processCounts(fileconnection,
                  negStrandPositionAdjustment,
                  windowsize,
                  pseudocount,
                  delim)
    counts.outputAllChromOEM()

# def processOEMsFromLoessCSV
def printStrandSwitchMetricsFromLoessCSV(fileconnection,
                            negStrandPositionAdjustment,
                            windowsize=300,
                            pseudocount=1e-9):
    ## initial
    countsOf5pEnds = PositionCounts(negStrandPositionAdjustment=negStrandPositionAdjustment,
                                    windowsize=windowsize,
                                    pseudocount=pseudocount)
    for line in fileconnection:
        line = line[:-1].split(",")
        chrom = line[0]
        pos = int(line[1])
        count = float(line[2])
        if count < 0:
            strand = 16
            count = count*-1
        else:
            strand = 0
        countsOf5pEnds.addCount(chrom,pos,strand,count)
##        if countsOf5pEnds.sameChrom(chrom):
##            countsOf5pEnds.addCount(pos,strand,count)
##        else:
##            countsOf5pEnds.outputChromOEM(chrom)
##            countsOf5pEnds.addCount(pos,strand,count)
##    ## print last chr
##    countsOf5pEnds.outputChromOEM(chrom)
    countsOf5pEnds.outputAllChromOEM()

#def processOEMsFromOEMinput
def printStrandSwitchMetricsFromOEMstyleInput(fileconnection,
                            negStrandPositionAdjustment,
                            windowsize=300,
                            pseudocount=1e-9):
    ## initial
    countsOf5pEnds = PositionCounts(negStrandPositionAdjustment=negStrandPositionAdjustment,
                                    windowsize=windowsize,
                                    pseudocount=pseudocount)
    for line in fileconnection:
        line = line[:-1].split(",")
        chrom = line[0]
        pos = int(line[1])
        fwd = float(line[2])
        rev = float(line[3])
        countsOf5pEnds.addCount(chrom,pos,0,fwd)
        countsOf5pEnds.addCount(chrom,pos,16,rev)
    ## print 
    countsOf5pEnds.outputAllChromOEM()


##############################################################################
''' EXECUTION '''
##############################################################################

if __name__ == '__main__':
    if args.input in ("-", "stdin"):
        fileconnection = sys.stdin
    else:
        fileconnection = open(args.input, 'r')


    if args.inputformat == 1:
        printStrandSwitchMetricsFromOEMstyleInput(fileconnection,
                                negStrandPositionAdjustment=args.negStrandPositionAdjustment,
                                windowsize=args.windowsize,
                                pseudocount = args.pseudocount)
    elif args.inputformat == 2:
        if args.outputformat == 1:
            printOEMstyleInput(fileconnection,
                          negStrandPositionAdjustment = args.negStrandPositionAdjustment,
                          windowsize = args.windowsize,
                          pseudocount = args.pseudocount)
        else: #2
            printStrandSwitchMetrics(fileconnection,
                        negStrandPositionAdjustment = args.negStrandPositionAdjustment,
                        windowsize = args.windowsize,
                        pseudocount = args.pseudocount)
    elif args.inputformat == 3:
        printStrandSwitchMetricsFromLoessCSV(fileconnection,
                                negStrandPositionAdjustment = args.negStrandPositionAdjustment,
                                windowsize = args.windowsize,
                                pseudocount = args.pseudocount)


