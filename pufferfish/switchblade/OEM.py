######## DEPRECATED
######## LIVE VERSION @ https://github.com/JohnUrban/switchblade

#!/usr/bin/env python2.7
import sys
import numpy as np
import argparse
## 

## Date Created: October 02, 2013
## Author: John Urban
##
## June 11,2014::  modified to take variable step csv -- 
##
##	Future improvements to this script might include
##              - writing the output file directly instead of building up and storing the oemDict
##		    I imagine this will help on the human genome
##              - making it write out in correct order
##                    as of now, the order is scrambled b/c it loops over keys (chromosome names) of a dictionary
##
## September 2021::
##      - modified to use argparse.




##############################################################################
''' ARGUMENT PARSER '''
##############################################################################

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Origin Efficiency Metric.

    Usage:

        OEM.py -i inputCSV [-w W -p pseudocount -o outname.ext -H --csvformat]

    INPUT:
    CSV containing chromosome position and coverage for both strands.

    OUTPUT:
    bedGraph or CSV of OEM for each position.


        
          input CSV: chr,pos,fwd,rev

          outputFormat = either 'bdg' or 'csv'
              bdg will return a bedGraph, a tab-delimited BED3 + score format: chr start end oem
                  Note that bedGraphs (and BED files) are 0-based, half open intervals
                  Thus, the first base of a sequence would be identified as:  chr  0  1  score
                   i.e. start = 0, end = 1 -- but the end is not included, thus only base 0 is represented here
              csv will return a CSV file similar to the input of format: chr,pos,oem
                  The CSV output (and input) is 1-based. Thus, base 1 of a sequence would be represented as: chr,1,oem

          

              """, formatter_class= argparse.RawTextHelpFormatter)

parser_input = parser.add_mutually_exclusive_group()
parser_input.add_argument('-i', "--input",
                   type= str, default="-",
                   help='''input CSV: chr,pos,fwdcov,revcov. Use "-" or "stdin" or leave blank during piping.''')

parser.add_argument('-o', "--outname",
                   type= str, default="-",
                   help='''Defaults to stdout. Otherwise, specify a filename for output: e.g. output.bedGraph or output.csv. Use "-" or "stdout" or leave blank during piping. ''')

parser.add_argument('-H', '--hasheader', action='store_true', default=False,
                    help='''Use this flag if the inputCSV has the header line 'chr,pos,fwd_str,rev_str'.''')

parser.add_argument('-c', '--csvformat', action='store_true', default=False,
                    help='''Use this flag for CSV output formatting instead of bedGraph formatting.''')

parser.add_argument('-w', '--windowsize', type=int, default=10000,
                    help='''Window size to use. Default = 10000 (10 kb). windowSize is how big the window is to left of a position (for WL and CL calculations) and to the right (for WR and CR). This should be set to '10000' if unsure.''')
parser.add_argument('-p', '--pseudocount', type=int, default=0,
                    help='''Pseudocount to use. Default = 0. The pseudocount is used to prevent division by 0.
              In 10kb bins on yeast genome it is unlikely to get a division-by-zero error.
              This is not necessarily true for the human genome or smaller window sizes, thus a pseudocount can help.
              Suggested pseudocount if needed is 1, but using sequencing depth and other knowledge to decide is better.''')


args = parser.parse_args()




##############################################################################
''' FUNCTIONS '''
##############################################################################

def okazakiFileReader(inputfile, header=True):
    """Takes in .csv a la Jason Belsky, returns Dict with that info """
    # Open connection to input.csv
    if inputfile in ('-','stdin'):
        okazaki = sys.stdin
    else:
        okazaki = open(inputfile, "r")

    # If header line present, store it
    if header:
        headerLine = okazaki.readline()[:-1].split(",")

    # initialize dictionary for storing information from file
    okazakiDict = {}

    # Parse file line by line
    for line in okazaki:
        line = line[:-1].split(",")
        chromosome, pos, fwd, rev = line[0], int(line[1]), float(line[2]), float(line[3])
        try:
            #okazakiDict[chromosome]
            okazakiDict[chromosome]["pos"] += [pos]
            okazakiDict[chromosome]["fwd"] += [fwd]
            okazakiDict[chromosome]["rev"] += [rev]
        except KeyError:
            okazakiDict[chromosome] = {}
            okazakiDict[chromosome]["pos"] = [pos]
            okazakiDict[chromosome]["fwd"] = [fwd]
            okazakiDict[chromosome]["rev"] = [rev]

    if inputfile in ('-','stdin'):
        okazaki.close()
    return okazakiDict


def OEMDict(okazakiDict, windowSize=500, pseudocount=0):
    """Takes in okazakiDict, returns oemDict"""
    # chr names
    chromosomes = okazakiDict.keys()

    # initialize oemDict
    oemDict = {}

    # start loop over all chromosomes
    for chromosome in chromosomes:
        sys.stderr.write( "Processing chromosome: " + str(chromosome) + "\n")
        chrLen = len(okazakiDict[chromosome]["pos"])
        sys.stderr.write( "Chromosome Length = " + str(chrLen) + "\n" )

        ## GIVE DUMMY LINE AND STOP IF CHRLEN < MinLen == (w+1)*2
        if chrLen < (windowSize+1)*2:
            oemDict[chromosome] = {'pos':[1], 'oem':[0]}
            continue


        ### CONTINUE IF CHRLEN >= MinLen == (w+1)*2
        # Set the start and end position (end position is 10,000 bp less than last position) (JB)
        start = windowSize-1
        end = chrLen - windowSize ## make sure this works
        
        # Calculate the first window (JB) -- float() used to avoid default 'int division'
        WL = float(sum(okazakiDict[chromosome]['fwd'][0:start+1])) + pseudocount ## goes up to and includes start as in R
        CL = float(sum(okazakiDict[chromosome]['rev'][0:start+1])) + pseudocount
        WR = float(sum(okazakiDict[chromosome]['fwd'][start+1:start+1+windowSize])) + pseudocount ## starts at 1 above start as in R
        CR = float(sum(okazakiDict[chromosome]['rev'][start+1:start+1+windowSize])) + pseudocount

        # Set up the storage oem (JB)
        pos = okazakiDict[chromosome]['pos'][start:end] ##+!
        oemDict[chromosome] = {'pos':pos, 'oem':[0]*len(pos)}

        oem = (WL / (WL + CL)) - (WR / (WR + CR))

        oemDict[chromosome]['oem'][0] = oem


        ## Iterative steps
        numPositions = len(range(start,end))
        percentsToReport = range(0,101,10)

        Pos = start
        for i in range(1, numPositions):
            # report progress
            ## -- since this script goes so fast, progress other than which chr is currently being worked on was not necessary to report
            ##percentComplete = int(100.0*(i+1)/numPositions) ## i+1 b/c initialization added in

            ##if percentComplete in percentsToReport:
                ##print str(percentComplete) + "% complete for chromosome: " + str(chromosome)
                ##percentsToReport[int(percentComplete/10.0)] = -1

            # Update the pos
            Pos += 1

            # Update WL, CL, WR, CR
            WL = WL + okazakiDict[chromosome]['fwd'][Pos] - okazakiDict[chromosome]['fwd'][Pos - windowSize]
            CL = CL + okazakiDict[chromosome]['rev'][Pos] - okazakiDict[chromosome]['rev'][Pos - windowSize]
            #WR = WR + okazakiDict[chromosome]['fwd'][Pos+windowSize-1] - okazakiDict[chromosome]['fwd'][Pos]
            #CR = CR + okazakiDict[chromosome]['rev'][Pos+windowSize-1] - okazakiDict[chromosome]['rev'][Pos]
            WR = WR + okazakiDict[chromosome]['fwd'][Pos+windowSize] - okazakiDict[chromosome]['fwd'][Pos]
            CR = CR + okazakiDict[chromosome]['rev'][Pos+windowSize] - okazakiDict[chromosome]['rev'][Pos]
            ## I am considering making it directly centered on pos by having the first window end on pos and 2nd window start on pos
            ## No reason that shouldn't/cannot be done and it is more reflective of the pos
            ## This will require some changes in the other parts of script as well such as end position

            # Store the oem
            oem = (WL / (WL + CL)) - (WR / (WR + CR))

            oemDict[chromosome]['oem'][i] = oem

    return oemDict


def writeOEM2Bdg(oemDict, out):
    """Takes in oemDict and writes out a bedgraph
    This will equate to chromosome, oemDict[chromosome]['pos'][i]-1, oemDict[chromosome]['pos'][i], oemDict[chromosome]['oem'][i]
    Note that bedGraphs are 0-based, half-open intervals"""
    chromosomes = oemDict.keys()
    for chromosome in chromosomes:
        numPositions = len(oemDict[chromosome]['pos'])
        for i in range(numPositions):
            chrom = str(chromosome)
            start = str(oemDict[chromosome]['pos'][i]-1)
            end = str(oemDict[chromosome]['pos'][i])
            oem = str(oemDict[chromosome]['oem'][i])
            out.writelines('\t'.join([str(e) for e in [chrom, start, end, oem]]) + "\n")

def writeOEM2CSV(oemDict, out):
    """Takes in oemDict and writes out a CSV similar to input CSV: chr,pos,oem
    This will equate to chromosome, oemDict[chromosome]['pos'][i], oemDict[chromosome]['oem'][i]
    Note that this outputCSV is a 1-based position as was the inputCSV
    If you are to make a bedGraph from the commandline, it will require something like:
        awk '{gsub(/,/,"\t"); print}' output.csv | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $3}' > output.bedGraph"""
    chromosomes = oemDict.keys()
    for chromosome in chromosomes:
        numPositions = len(oemDict[chromosome]['pos'])
        for i in range(numPositions):
            chrom = str(chromosome)
            pos = str(oemDict[chromosome]['pos'][i])
            oem = str(oemDict[chromosome]['oem'][i])
            out.writelines(','.join([str(e) for e in [chrom, pos, oem]]) + "\n")


    
##############################################################################
''' EXECUTION '''
##############################################################################

##If run from the linux commandline then execute the following
#### See 6.1.1 here for more info on how this works: http://docs.python.org/2/tutorial/modules.html
if __name__ == "__main__":

    #EXECUTE
    okazakiDict = okazakiFileReader(args.input, args.hasheader)
    oemDict = OEMDict(okazakiDict, args.windowsize, args.pseudocount)


    ## OUT
    if args.outname in ("-","stdout"):
        out = sys.stdout
    else:
        out = open(args.outname, 'w')
        
    ## WRITING OUT
    if args.csvformat:
        writeOEM2CSV(oemDict, out)
    else:
        writeOEM2Bdg(oemDict, out)

    ## QUIT
    if args.outname in ("-","stdout"):
        out.close()
    quit()


