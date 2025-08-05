#!/usr/bin/env python2.7
import os.path
import sys
import argparse

#logger
import logging
logger = logging.getLogger('pufferfish')

# pufferfish imports
#import pufferfish.version
pfv = '0.0.0'


## FUNCTIONS
def run_subtool(parser, args):
    ## the function to be used in subparser.set_defaults(func=func)
    ## it selects the module that the sub-parser is made for
    ## then uses the "run" function in that module
    if args.command == 'mapreads':
        import mapreads as submodule
    elif args.command == 'getcov':
        import getcov as submodule
    elif args.command == 'findpuffs':
        import findpuffs as submodule
    elif args.command == 'dump':
        import pk2txt as submodule
    elif args.command == 'puffcn':
        import cn as submodule
    elif args.command == 'puffcnpy':
        import cnpy as submodule
    elif args.command == 'summits':
        import findsummits as submodule
    elif args.command == 'normalize':
        import normalize as submodule
    elif args.command == 'hmm':
        import generalhmm as submodule
    elif args.command == 'generate':
        import generate as submodule
    elif args.command == 'filter':
        import filterfish as submodule
    elif args.command == 'help':
        import helper as submodule
    # run the chosen submodule
    submodule.run(parser, args)


## This subclass is used as parser_class for submodule sub-parsers
class ArgumentParserWithDefaults(argparse.ArgumentParser):
    ## Child/sub of argparse.ArgumentParser class (the super)
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        ## Add arguments that can be used by all sub-commands
        self.add_argument("-q", "--quiet", help='''Do not output warnings to stderr.''',
                          action="store_true",
                          dest="quiet")



    
    
def main():
    logging.basicConfig()

    ## Create top-level parser
    parser = argparse.ArgumentParser(prog="pufferfish", description=''' PufferFish - HMM-based approach(es) to finding and analyzing developmentally regulated amplicons (genomic sites that are programmed to increase in copy number over time).''',
                                     formatter_class=argparse.RawTextHelpFormatter)#ArgumentDefaultsHelpFormatter
    parser.add_argument('-v', '--version', help='''Installed pufferfish version.''',
                        action='version',
                        version='%(prog)s ' + pfv)#str(pufferfish.version.__version__))
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command',
                                       parser_class=ArgumentParserWithDefaults)


    ## Create a sub-command parser for mapreads
    parser_mapreads = subparsers.add_parser('mapreads',
                                            help='''Depends on Bowtie2 and SAMtools.
Maps fasta/fastq files to genome (can provide bt2 index or fasta reference (which will first be converted to bt2 index).
Maps reads to genome and filters out unmapped reads before sorting and indexing.''')
    parser_mapreads.add_argument('fastxfiles', metavar='fastxfiles', nargs='+', type=str,
                                 help='''Paths to as many fasta/fastq files as youd like to map. Can be gz or bz2.''')
    parser_mapreads_reftype = parser_mapreads.add_mutually_exclusive_group(required=True)
    parser_mapreads_reftype.add_argument('-b', '--bt2', type=str,
                               help='''Path to bt2 index prefix.''')
    parser_mapreads_reftype.add_argument('-r', '--ref_fasta', type=str,
                               help='''Path to reference fasta. A bt2 index will be generated in same dir.''')
    parser_mapreads.add_argument('--threads', '-p', type=int, default=1, help='''Number of parallel threads for bowtie2 to use. Default: 1.''')

    parser_mapreads.add_argument('--dry', action='store_true', default=False, help='''Only writes out the commands that will be used if set.''')

    parser_mapreads.set_defaults(func=run_subtool)


##    ## Create a sub-command parser for filterdup
##    parser_mapreads = subparsers.add_parser('filterdup',
##                                            help='''Depends on Picard Tools 2.1.1, BEDtools, pybedtools, pysam.
##Remove optical duplicates and marks PCR duplicates.
##All PCR duplicates except K at a given site are removed.
##K is determined by a binomial calculuation using a bin size and number of reads in a given bin.
##Then any duplicates in that bin are subject to filtering down to K.
##1. Remove optical duplicates and mark PCR duplicates.
##2. Make bins
##3. Get read count in those bins
##4. For each bin, check if there are marked reads. If so, calculate K and filter.
##5. write remaining reads as you go...
##''')
##    parser_filterdup.add_argument('bams', metavar='bams', nargs='+',
##                               type=str,
##                               help=''' Paths to BAM files that need duplicate filtering.''')
##    parser_filterdup.add_argument('-g', '--genome', type=str,
##                               help='''Path to BEDtools genome file describing reference reads were mapped to.''')
##    parser_filterdup.add_argument()
##    parser_filterdup.add_argument('--dry', action='store_true', default=False, help='''Only writes out the commands that will be used if set.''')
##
##    parser_filterdup.set_defaults(func=run_subtool)
##    

    ## Create sub-command parser for getcov
    ## TODO add filterdup possibility from macs2... rm pcr dups
    parser_getcov = subparsers.add_parser('getcov',
                                          help=''' Depends on SAMtools and BEDtools.''')
    parser_getcov.add_argument('bams', metavar='bams', nargs='+',
                               type=str,
                               help=''' Paths to as many bam files as you need to get coverage for.
Can include a file-of-filenames (FOFN) and tarballs as well.''')
    
    parser_getcov.add_argument('-f','--filterdup', type=str,
                               help='''Provide /path/to/picard.jar  (picard v2.1.1 or higher)''')
    parser_getcov.add_argument('-g', '--genome', type=str, required=True,
                               help='''Path to file.genome as needed and defined by BEDTools. See "bedtools makewindows" or "bedtools coverage"''')
    parser_getcov.add_argument('-w', '--window', type=str, default='500',
                               help='''Integer window size - will be counting in windows of this size. Default: 500.''')
    parser_getcov.add_argument('-s', '--step', type=str, default='500',
                               help='''Integer step size - will slide window over this much. Default: 500.''')
    parser_getcov.add_argument('-Q', '--mapq', type=str, default='0',
                               help='''Integer mapq cut-off - only include reads with mapping quality >= INT. Default: 0.''')
    parser_getcov.add_argument('-m', '--picardmem', type=str, default='4g',
                               help='''Provide memory needed/available to Picard MarkDuplicates as integer_letter string, such as 500m, 1g, 2g, 64g, etc. Default: 4g.''')
    parser_getcov.add_argument('--keepopt',action='store_true', default=False, help='''Optical duplicates are removed by default. This flag says to mark them instead.''')
    parser_getcov.add_argument('--rmdup',action='store_true', default=False, help='''PCR duplicates are marked by default. This flag will result in removing them (as well as optical duplicates).''')
    parser_getcov.add_argument('--dry',action='store_true', default=False, help='''Only writes out the commands that will be used if set.''')
    parser_getcov.add_argument('--clean',action='store_true',default=False,help='''Remove intermediate files... Default: False.''')
    parser_getcov.add_argument('--force',action='store_true',default=False,help='''Ignore assertions. Good for made-up filenames when debugging in dry-runs. Do not use this for real run. Default: False.''')
    parser_getcov.set_defaults(func=run_subtool)











    ## Create sub-command parser for findpuffs
    parser_findpuffs = subparsers.add_parser('findpuffs',
                                             help='''Take in getcov bedGraphs, do stuff.''')
    parser_findpuffs_input = parser_findpuffs.add_mutually_exclusive_group(required=True)
    parser_findpuffs_input.add_argument('-i','--input', type=str,
                                  help='''Input file -- a tab-sep file with 2 columns (stage number, filename) with a single line for all getcov bedGraph files you wish to include.
Example:
1\tstage1.bedGraph''')
    parser_findpuffs_input.add_argument('-ip','--inpickle', type=str,
                                        help='''Pickle file (e.g. data.pk) containing already processed getcov bedGraphs as MultiCovBed object.''')

    parser_findpuffs.add_argument('-op','--outpickle', type=str, default='data.fp.pk',
                                        help='''Name for output pickle file (e.g. data.fp.pk.gz) that will contain the MultiCovBed object made when this is run.
Pickled data is automatically gzipped. If .gz not at end of given filename, it will be added.
If the filename exists it will be erased and written over.
Default: data.fp.pk.gz''')
    
    parser_findpuffs.add_argument('-s1', '--smoothbeforemedian', default=False, type=int,
                                  help='''Smooth counts in bins before finding the median bin counts (and before median normalization).
Must provide integer window size to smooth in (should be longer than bin size)
-- e.g. if bin size = 500, smoothing bandwidth could be 10000.
One probably should not do both --smoothbeforemedian and --smoothbeforecor. Pick one or none.''')
    parser_findpuffs.add_argument('-s2', '--smoothbeforecor', default=False, type=int,
                                  help='''Smooth counts in bins after median normalization, but before finding correlations in each bin.
Must provide integer window size to smooth in (should be longer than bin size)
-- e.g. if bin size = 500, smoothing bandwidth could be 10000.
One probably should not do both --smoothbeforemedian and --smoothbeforecor. Pick one or none.''')

    parser_findpuffs.add_argument('-bw', '--corsmoothbandwidth', type=int, default=15000,
                                  help='''For smoothing correlation scores. Provide integer window size to smooth in (should be longer than bin size)
-- e.g. if bin size = 500, smoothing bandwidth could be 10000. Default: 15000.''')

    parser_findpuffs.add_argument('-mep', '--max_empty_bin_pct', type=float, default=0.4,
                                  help='''For filtering contigs out, contig is allowed to have up to X (proportion between 0 and 1) bins with 0 coverage. Default: 0.4.''')

    parser_findpuffs.add_argument('-mop', '--max_offending_samples_pct', type=float, default=0.4,
                                  help='''For filtering contigs out, contig is allowed to exceed max_empty_bin_pct in Y (proportion between 0 and 1) of the samples. Default: 0.4.''')


    parser_findpuffs.set_defaults(func=run_subtool)
    














    ## Create sub-command for dump
    parser_dump = subparsers.add_parser('dump',
                                                help=''' Take in pickled object containing MultiCovBed object where statepath has been found.
Output bedGraph of statepath and/or BED file containing coordinates of states.''')
    parser_dump.add_argument('-ip', '--inpickle', type=str, required=True,
                                    help='''Path to input pickle.''')
    parser_dump.add_argument('-p','--prefix', type=str, required=True,
                                     help='''Output files with provided --prefix''')
    parser_dump.add_argument('-c', '--counts', default=False, action='store_true',
                             help='''Output counts as expanded single-step, bedGraph-like file. Counts will be normalized.
Columns 4+ will have counts from files in order files were given.
Can use this with awk to create bedGraphs for each -- e.g. for i in {4..8}; do awk -v "i=$i" 'OFS="\t" {print $1,$2,$3,$i}' q30.w500.s500.default.counts.bedGraph > q30.w500.s500.default.counts.$i.bedGraph; done ''')
    parser_dump.add_argument('-fc', '--filtered_counts', default=False, action='store_true',
                             help='''Output counts from filtered contigs as expanded single-step, bedGraph-like file. Counts will likely NOT be normalized (as filtereing is done prior to median normalization).
Columns 4+ will have counts from files in order files were given.
Can use this with awk to create bedGraphs for each -- e.g. awk 'OFS="\t" {print $1,$2,$3,$4}' filtered_counts.bedGraph > file.filtered_counts.bedGraph ''')
    parser_dump.add_argument('-vc','--viterbi_collapsed', default=False, action='store_true',
                                     help='''Output viterbi statepath expanded single-step bedGraph with provided --prefix''')
    parser_dump.add_argument('-ve','--viterbi_expanded', default=False, action='store_true',
                                     help='''Output viterbi statepath collapsed varstep bedGraph with provided --prefix''')
    parser_dump.add_argument('-cc','--correlations_collapsed', default=False, action='store_true',
                                     help='''Output bin correlation scores as collapsed varstep bedGraph with provided --prefix''')
    parser_dump.add_argument('-ce','--correlations_expanded', default=False, action='store_true',
                                     help='''Output bin correlation scores as expanded single-step bedGraph with provided --prefix''')
    parser_dump.add_argument('-scc','--smoothed_correlations_collapsed', default=False, action='store_true',
                                     help='''Output smoothed bin correlation scores as collapsed varstep bedGraph with provided --prefix''')
    parser_dump.add_argument('-sce','--smoothed_correlations_expanded', default=False, action='store_true',
                                     help='''Output smoothed bin correlation scores as expanded single-step bedGraph with provided --prefix''')
    parser_dump.add_argument('-sc','--slopes_collapsed', default=False, action='store_true',
                                     help='''Output bin slopes as collapsed varstep bedGraph with provided --prefix''')
    parser_dump.add_argument('-se','--slopes_expanded', default=False, action='store_true',
                                     help='''Output bin slopes as expanded single-step bedGraph with provided --prefix''')

    parser_dump.set_defaults(func=run_subtool)














    ## create sub-sommand for puffcn (puff copy number)
    parser_puffcn = subparsers.add_parser('puffcn',
                                          help = '''Given a latest-stage sample (where all or most puffs have grown) and an optional earliest stage sample
(for additional Fold-enrichment normalization), define whether a region is best explained by cn=1,2,4,8,16,32,64.
Ideally, this can give an idea of where replication forks approximately reach from each firing.''')
    parser_puffcn.add_argument('-l','--latestage', type=str, required=True,
                                help='''Provide path to bedGraph (e.g. made from getcov) for a late stage sample.''')
    parser_puffcn.add_argument('-e','--earlystage', type=str, required=False, default=False,
                               help=''' Optional: Provide path to bedGraph (e.g. made from getcov) for an early stage sample. This is used after smoothing and median normalization to further normalize the late-stage sample (e.g. can correct for sequencing biases)''')
    parser_puffcn.add_argument('--replace', action='store_true', default=False,
                                        help='''Turn on "replace" functionality. By default this will replace '.' in the count column of bedGraphs with '0'.
Use --replace_with and --replace_this to change.''')
    parser_puffcn.add_argument('--replace_this', type=str, default='.',
                                        help='''Used with --replace. Specify the character in count column to replace. Default = '.' ''')
    parser_puffcn.add_argument('--replace_with', type=str, default='0',
                                        help='''Used with --replace. Specify the character to replace the --replace_this character with.
Must be a string that can be converted to a float. Default = '0' ''')


    ## PROTOCOLS
    parser_puffcn_protocol = parser_puffcn.add_mutually_exclusive_group(required=True)
    parser_puffcn_protocol.add_argument('-1', '--protocol1', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are median normalized.
Then late stage is normalized to early stage if available.
Then the HMM is run.''')
    parser_puffcn_protocol.add_argument('-2', '--protocol2', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first smoothed with bandwidth given by --bandwidth.
Then they are median normalized.
Then late stage is normalized to early stage if available.
Then the HMM is run.''')
    parser_puffcn_protocol.add_argument('-3', '--protocol3', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first median normalized.
Then they are smoothed with bandwidth given by --bandwidth.
Then late stage is normalized to early stage if available.
Then the HMM is run.
Note: if early is not present, this is same as protocol 4.''')
    parser_puffcn_protocol.add_argument('-4', '--protocol4', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first median normalized.
Then late stage is normalized to early stage if available.
Then late/early is smoothed with bandwidth given by --bandwidth.
Then the HMM is run.
Note: if early is not present, this is same as protocol 3.''')

    parser_puffcn_protocol.add_argument('-5', '--protocol5', action='store_true', default=False,
                                        help='''Late stage is normalized to early stage if available.
Then late/early is smoothed with bandwidth given by --bandwidth. (i.e. L/E -> smooth)
Then the HMM is run.
Note: if early is not present, this is same as protocol 6.''')

    parser_puffcn_protocol.add_argument('-6', '--protocol6', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first smoothed with bandwidth given by --bandwidth.
Then late stage is normalized to early stage if available. (i.e. smooth -> L/E)
Then the HMM is run.
Note: if early is not present, this is same as protocol 5.''')
    parser_puffcn_protocol.add_argument('-7', '--protocol7', action='store_true', default=False,
                                        help='''Late stage is normalized to early stage if available.
Then the HMM is run.
Note: No median normalization or smoothing is performed. If only late is given, then this is just an identity/pass-through function.''')

    parser_puffcn_protocol.add_argument('-8', '--protocol8', action='store_true', default=False,
                                        help='''SKEW. Only accepts one file (e.g. latestage).
                                                For file with N rows, returns N-1 rows.
                                                For each value, V[i], in 4th column, for i in 2:N, it computes Skew = (V[i]-V[i-1]) / (V[i]+V[i-1])''')

    parser_puffcn_protocol.add_argument('-9', '--protocol9', action='store_true', default=False,
                                        help='''PERCENT CHANGE. Only accepts one file (e.g. latestage).
                                                For file with N rows, returns N-1 rows.
                                                For each value, V[i], in 4th column, for i in 2:N, it computes Skew = 100*(V[i]-V[i-1]) / V[i-1]''')

    parser_puffcn_protocol.add_argument('-10', '--protocol10', action='store_true', default=False,
                                        help='''SKEW CHANGE or SKEW DERIVATIVE. Only accepts one file (e.g. latestage).
                                                For file with N rows, returns N-1 rows.
                                                For each value, V[i], in 4th column, for i in 2:N, it computes Skew = 100*(V[i]-V[i-1]) / V[i-1]''')

    parser_puffcn_protocol.add_argument('-11', '--protocol11', action='store_true', default=False,
                                        help='''PERCENT CHANGE DERIVATIVE. Only accepts one file (e.g. latestage).
                                                For file with N rows, returns N-1 rows.
                                                For each value, V[i], in 4th column, for i in 2:N, it computes Skew = 100*(V[i]-V[i-1]) / V[i-1]''')

    parser_puffcn_protocol.add_argument('-12', '--protocol12', action='store_true', default=False,
                                        help=''' Median ratio normalization. Late is normalized to early. Then those ratios are median normalized.
                                                This is similar to a global version of what's used in DEseq2 (or TMM for EdgeR).''')

    parser_puffcn_protocol.add_argument('-13', '--protocol13', action='store_true', default=False,
                                        help=''' Median ratio normalization with pre-smoothing. Late (and early if present) is smoothed. Late is normalized to early. Then those ratios are median normalized.
                                                This is similar to a global version of what's used in DEseq2 (or TMM for EdgeR).''')

    parser_puffcn_protocol.add_argument('-14', '--protocol14', action='store_true', default=False,
                                        help=''' Median ratio normalization with post-smoothing. Late is normalized to early. Then those ratios are smoothed. Then median normalized.
                                                This is similar to a global version of what's used in DEseq2 (or TMM for EdgeR).''')

    parser_puffcn_protocol.add_argument('-15', '--protocol15', action='store_true', default=False,
                                        help=''' Median ratio normalization with end-smoothing. Late is normalized to early. Then those ratios are median normalized. Then smoothed.
                                                This is similar to a global version of what's used in DEseq2 (or TMM for EdgeR).''')
    parser_puffcn_protocol.add_argument('-16', '--protocol16', action='store_true', default=False,
                                        help=''' Glocal Median ratio normalization. Late is normalized to early. Then those ratios are median normalized based on coverage value in late.
                                                This is similar to a coverage-based local version of what's used in DEseq2 (or TMM for EdgeR).''')
    
    parser_puffcn_protocol.add_argument('-17', '--protocol17', action='store_true', default=False,
                                        help=''' First normalizes late stage to median (X/Med), and scales it with --scalecov to make the median equal scalecov.
                                                It does the same to early if present.
                                                Then Late is normalized to Early (FE) using given pseudocount (not pseudocount given has same ratio to median in both samples).
                                                It then does median ratio normalization.
                                                It then returns it as is or logged if --log10 or --log2 specified''')

    parser_puffcn_protocol.add_argument('-18', '--protocol18', action='store_true', default=False,
                                        help='''Robust Z scores. If early (control) sample given, first is control norm (late/early, test/control) followed by robust Z. The fold-change calculation is not treated specially as it might be with other options. However, the robust Z should be the same or similar even if the t/c was not then median ratio normalized, for example, since that is just a scaling factor.''')

    parser_puffcn_protocol.add_argument('-19', '--protocol19', action='store_true', default=False,
                                        help='''Rank scores. If early (control) sample given, first is control norm (late/early, test/control) followed by ranking. The fold-change calculation is not treated specially as it might be with other options. However, the rank should be the same or similar even if the t/c was not then median ratio normalized, for example, since that is just a scaling factor.''')

    parser_puffcn_protocol.add_argument('-20', '--protocol20', action='store_true', default=False,
                                        help='''Robust Z score differences between two samples. Requires both late (test) and early (control) samples. Robust Z scores are calculated independently for each sample. Then R_control is subtracted from R_test = R_t - R_c.''')

    parser_puffcn_protocol.add_argument('-21', '--protocol21', action='store_true', default=False,
                                        help='''Rank score differences between two samples. Requires both late (test) and early (control) samples. Rank scores are calculated independently for each sample. Then R_control is subtracted from R_test = R_t - R_c.''')

    parser_puffcn_protocol.add_argument('-22', '--protocol22', action='store_true', default=False,
                                        help='''Signal Per Million Reads (Or Counts or Per Million whatever can be summed up in the 4th column). Bin_spmr = 1e6*Bin/Sum(Bins).
Use --SPXR to change scaling factor from 1e6 to X. If early (control) given, both are SPMR'd independently, then late/early (test/control).
When an early (control) sample is provided, you may also want to check the default pseudocount applied.''')

    parser_puffcn_protocol.add_argument('-23', '--protocol23', action='store_true', default=False,
                                        help='''Rank standardize scores. First rank, then subtract and divide by middle: (r-M)/M, where r is a bins rank, and M is the theoretical middle rank: M=(min+max)/2. If early (control) sample given, first is control norm (late/early, test/control) followed by ranking. The fold-change calculation is not treated specially as it might be with other options. However, the rank should be the same or similar even if the t/c was not then median ratio normalized, for example, since that is just a scaling factor.''')


    parser_puffcn.add_argument('--stringcols', action='store_true', default=False,
                               help='''Just treat columns other than 4 as strings...''')
    parser_puffcn.add_argument('--log2', action='store_true', default=False,
                               help='''Return log2 values. Default = False.''')
    parser_puffcn.add_argument('--log10', action='store_true', default=False,
                               help='''Return log10 values. Default = False.''')
    parser_puffcn.add_argument('--scalecov', type=float, default=1,
                               help='''Multiply coverage by this as part of protocol 17.''')
    parser_puffcn.add_argument('--SPXR', type=float, default=1e6,
                               help='''In essence, this is like --scalecov with a different default: 1e6.''')    
    parser_puffcn.add_argument('--pseudoZeroBins', action='store_true', default=False,
                               help='''Not to be confused with --pseudo. This option applies only to protocols 24-27 right now. It only need be used when there are zeros in the control (early) sample. In protocols 26 and 27, this is likely to happen from the robust z-score pre-processing. If an error is thrown, try --pseudoZeroBins or --addMinOtherPlusOneToBoth. --pseudoZeroBins adds min(abs(values))/10 to bins in both samples that were 0 in contorl (early). --addMinOtherPlusOneToBoth shifts both distributions up by min(control values)+1, setting the min control value to 1.  These are not meant to be used together, but won't throw an error if they are.''')    
    parser_puffcn.add_argument('--addMinOtherPlusOneToBoth', action='store_true', default=False,
                               help='''This option applies only to protocols 24-27 right now. It only need be used when there are zeros in the control (early) sample. In protocols 26 and 27, this is likely to happen from the robust z-score pre-processing. If an error is thrown, try --pseudoZeroBins or --addMinOtherPlusOneToBoth. --pseudoZeroBins adds min(abs(values))/10 to bins in both samples that were 0 in contorl (early). --addMinOtherPlusOneToBoth shifts both distributions up by min(control values)+1, setting the min control value to 1. These are not meant to be used together, but won't throw an error if they are.''')    


    
    parser_puffcn.add_argument('-m', '--emodel', type=str, default='normal',
                               help='''Specify emissions model to assume for HMM. Options: normal, exponential, poisson, geometric, gamma, and discrete. Default: normal.
Note that all you ever need to do is given the expected means and standard deviations for each state from sampled data.
The normal model will use those directly. Poisson will use the means as lambda. Exponential will use the means as B w/ rate 1/B. Geometric will also use the means as 1/mu.
Gamma will estimate alpha and beta (shape and scale) parameters from the means and standard deviations - if you have A/B in mind, convert to mu and sigma by A*B and (A*B^2)^0.5.
Note that the exponential is same as the gamma when shape is set to 1, making mu and sigma equal B. Thus, if you do gamma w/ muscale=1, then you should get same as exponential.
"Discrete" is when you have a finite number of categories that occur at various frequencies like the sides of a coin or dice. This expects a nState X nSymbol matrix in --mu.
"Discrete" assumes the symbols emitted are sequential integers from 1:N where N = number states.''')
    parser_puffcn.add_argument('-p', '--path', type=str, default='viterbi',
                               help='''Specify whether to take state path defined by viterbi or posterior decoding. Options: viterbi, posterior. Default: viterbi.''')
    parser_puffcn.add_argument('-s', '--scale', action='store_true', default=False,
                               help='''Before going back into the HMM, re-scale counts back towards original magnitude by multiplying by median.
This will also scale the HMM state means by the median to cooperate.''')
    parser_puffcn.add_argument('-c', '--collapsed', action='store_true', default=False,
                               help='''Return collapsed variable-step bedGraph instead of expanded single-step bedGraph.
This is often a much smaller file.''')
    parser_puffcn.add_argument('-ps', '--pseudo', type=float, default=0.1,
                               help=''' Before normalizing late to early, add this pseudocount to all counts in order to avoid division by zero.
Should be between 0 and 1.
Should be small enough to not change other values much,
but big enough such that numbers divided by 0+pseudo do not become massive.
Default: 0.1.''')
    parser_puffcn.add_argument('-bw', '--bandwidth', type=int, default=2500,
                               help=''' If kernel smoothing, specify bandwidth (int).
Bandwidth should be bigger when no early stage normalization to try to smooth out sequencing biases, mappability biases, etc.
Default: 2500.''')
    parser_puffcn.add_argument('--endsmoothing', action='store_true', default=False,
                               help=''' Add smoothing to the absolute end of any of the protocols for more flexibility here. This comes after log-transformation steps, for example, which optionally comes at the end of any protocol.''')
    parser_puffcn.add_argument('--replaceNaN', action='store_true', default=False,
                               help=''' If kernel smoothing, NaNs can be generated. This option replaces those with local averages (see --localwinsize, default=5 bins). In cases where localaverages return NaN (very rare), it fills NaN with the global average for the given chrom/sequence (not whole geneome, so still local-ish).''')
    parser_puffcn.add_argument('--localwinsize', type=int, default=5,
                               help=''' If kernel smoothing and/or using --replaceNan, this specifies the number of bins to use (centered on this so odd numbers preferable). Default = 5.''')
    parser_puffcn.add_argument('--impute', type=int, default=False,
                               help=''' If imputing, specify bandwidth (int) for  kernel smoothing.
This bandwidth is generally longer than the one you would provide for regular smoothing.
Only bins with a count of 0 will take on smoothed (imputed) values.
Try: 10000.
NOTE: In practice, this lead to its own set of problems and I do not recommend using it in its current form.''')
    parser_puffcn.add_argument('--counts', type=str, default=False,
                               help=''' Use this flag and specify an output prefix for the final normalized late stage bin counts bedGraph.''')
    parser_puffcn.set_defaults(func=run_subtool)
    parser_puffcn.add_argument('--levels', action='store_true', default=False,
                               help=''' Use this flag to output levels bdg instead of state bdg.
Levels are the means 1,2,4,8,16,32,64.
These are obtained by 2**(state-1).''')


    parser_puffcn.add_argument('--mu', '--discreteEmat', type=str, default='1,2,4,8,16,32,64',
                               help=''' PuffCN has been optimized for mapping DNA puffs in the fungus fly.
The default state means were previously hard-coded.
This option allows some flexibility from the command-line to change the state means.
Default: 1,2,4,8,16,32,64
To change: Provide comma-seprated list of state means.
The number of states will be calculated from this list.
If changing state sigmas (used in normal model), it must have same number of states represented.

NOTE: If using exponential or geometric distribution, provide the expected mean RCN values of the states
    as you would for normal or poisson models. This script will automatically take their inverses to work
    in the exponential and geometric models.

NOTE2: For "--emodel discrete" can call --mu as --discreteEmat for better-readability at commandline.
    Instead of a comma-sep list of means, provide comma-/semicolon-separated values to make up a nState X nSymbol matrix.
    Example of a 3-state x 4 symbol matrix: "0.97,0.01,0.01,0.01;0.01,0.97,0.01,0.01;0.01,0.01,0.01,0.97".
    That can be thought of as 3 4-sided dice.''')

    parser_puffcn.add_argument('--sigma', type=str, default=None,
                               help=''' PuffCN has been optimized for mapping DNA puffs in the fungus fly.
The default state sigmas (stdevs) were previously hard-coded.
This option allows some flexibility from the command-line to change the state sigmas.
Default: if not changed, defaults to square root of state means (Poisson-like).
To change: Provide comma-seprated list of state sigmas.
Alternatively: Use --mu_scale (default False) with a scaling factor multiplied against the MUs.
The number of states is calculated from this state mean list, which defaults to 7.
If changing state sigmas (used in normal model), it must have same number of states represented as state means.''')

    parser_puffcn.add_argument('--mu_scale', type=float, default=None,
                               help=''' See --sigma for more details on sigmas.
Use this to scale means (--mu) to use as stdevs (sigma) instead of taking square roots of means.
For example, --mu_scale 0.5 will use mu*0.5 as the stdev.''')

##    parser_puffcn.add_argument('--changestate', type=float, default=0.001,
##                               help=''' PuffCN has been optimized for mapping DNA puffs in the fungus fly.
##The default transition probabilities were previously hard-coded.
##For now, there are two parameters for transition probabilities: changing states or staying in a state.
##They are shard by all states.
##This option allows some flexibility from the command-line to change the probability of changing states.
##Default: 0.001.
##NOTE: To ensure transition probabilities from state i to j sum to 1,
##    the final transition probabilities will actually be X/sum(all trans i to j),
##    where X is 0.001 by default or what user has given.
##''')
##
##    parser_puffcn.add_argument('--samestate', type=float, default=0.999,
##                               help=''' PuffCN has been optimized for mapping DNA puffs in the fungus fly.
##The default transition probabilities were previously hard-coded.
##For now, there are two parameters for transition probabilities: changing states or staying in a state.
##They are shard by all states.
##This option allows some flexibility from the command-line to change the probability of staying in the current state (not changing).
##Default: 0.999.
##NOTE: To ensure transition probabilities from state i to j sum to 1,
##    the final transition probabilities will actually be X/sum(all trans i to j),
##    where X is 0.999 by default or what user has given.
##''')


    parser_puffcn.add_argument('--special_idx', type=int, default=0,
                               help='''Only for use if you're very familiar with the program (and change defaults).
The default state means is 1,2,4,8,16,32,64.
The default index for the mean that represents copy number 1 is 0.
In this lingo - CN=1 is the special state, and the 0-based index of the special state in that list is 0.
If you were to change parameters that affect where the special state is in the list, make sure to change this index.
This index is only used to help construct initial probabilities and transition probabilies.
If understood, it can be used to designate any single special state (not necessarily the one that corresponds to CN=1).
The other parameters to use with this are:
--init_special (probability of starting in the special state (usually CN=1).
    The probabity of starting in a another state (usually copy number variant states) defaults to (1-init_special)/(nstates-1).
--prob_leave_special
--prob_stay_special
--prob_other_to_special
--prob_other_to_other
--prob_other_to_self

Alternative to, an initial probability vector can be given with --initialprobs 
''')

    parser_puffcn.add_argument('--init_special', type=float, default=0.997,
                               help='''Probability of starting in the 'special state' (usually copy number = 1). Default: 0.997.
The probabity of starting in a another state (usually copy number variant states) defaults to (1-init_special)/(nstates-1).
''')



    parser_puffcn.add_argument('--leave_special_state', type=float, default=0.001,
                               help='''Probability of leaving the 'special state' (usually copy number = 1).
Default: 0.001.
If number is betwen 0 and 1, it will be assumed a probability.
If number given is > 1, then it will be treated as the average length (number of bins) of the special state.
For example, if 1000 is given, it will be 1/1000 = 0.001.
In terms of bp lengths, one would need to multiply n_bins * bin_length OR divide bp length by bin_length
Thus, if you want to see a change every 500 kb w/ 500 bp bins, then 500kb/500 = 1 kb = 1000 -- which will be interpreted as 0.001.
Or as another example, if you expect to see a change every 2 Mb with 100 bp bins, then 2e6/1e2 = 1e4 = 10 kb = 10000, interpreted as 0.0001.

The probability of staying in this state is the complement: 1-p
''')

    parser_puffcn.add_argument('--leave_other', type=str, default=None,
                               help='''Probability of leaving one of the other states.
This defaults to --leave_special_state making all transition probabilities out of states the same (0.001 by default).

To change, provide a probability of leaving (p).

If the  first number is betwen 0 and 1, it will be assumed a probability.
If the first number given is > 1, then it will be treated as the average length (number of bins).
For example, if 1000 is given, it will be 1/1000 = 0.001.

If only 1 number is given, then that is assumed to be the probability of transitioning to all the other states.

You can also give a comma-separated pair of 2 probabilities:
    prob of leaving to special state
    prob of leaving to another 'non-special' state.
Make sure the probabilities sum to what you expect the overall probability of leaving the state is...
    which should be p_to_special + p_to_nonspecial * (num_non_special-1) = p_to_special + p_to_nonspecial (nstates-2)

For example, in a 7-state model:
    0.001,0.0002 --> 0.001 + 0.0002 * 5 = 0.001 + 0.001 = 0.002
    OR
    0.001,0.001 --> 0.001 + 0.001 * 5 = 0.006

If the second number is > 1, the same rules apply as to the first number.

For other analyses, I've used:
0.00001,0.000000000001
OR
0.001,0.0000000001

The probability of staying in these states is the complement: 1-p1-p2

NOTE: the program forces the transition probabilities of a given state to sum to 1.
''')

    parser_puffcn.add_argument('--initialprobs', type=str, default=None,
                               help='''PuffCN has been optimized for mapping DNA puffs in the fungus fly.
The default state means were previously hard-coded.
This option allows some flexibility from the command-line to change the state means.
Default: [0.997, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005]
The default will change with more or less states described w/ --mu and --sigma.
By default, the first state will start out as 0.997 as above, all other states will be (1-0.997)/n_other_states.
That behavior also changes with following parameters:
--special_idx -- determines which state (not necessarily first) will be given default 0.997 (OR other with --initcn1)
--init_special (probability of starting in the special state (usually CN=1).
    The probabity of starting in a another state (usually copy number variant states) defaults to (1-init_special)/(nstates-1).
--leave_special_state

--prob_other_to_special
--prob_other_to_other
--prob_other_to_self
To change the initial probs manually: Provide comma-separated list of initial probs -e.g.: '0.997,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005'
This must have same number of states represented as state means (--mu; default 7).

''')


    parser_puffcn.add_argument('--transprobs', type=str, default=None,
                               help='''Provide trans probs with semi-colon separated rows that have comma-separated values.
E.g. a 2-state t matrix where self-self is 0.99 and self-other is 0.01 looks like: 0.99,0.01;0.01,0.99

''')
##'0.997,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005'  [0.997, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005]

    parser_puffcn.add_argument('--kmeans', type=int, default=None,
                               help='''PuffCN has been optimized for mapping DNA puffs in the fungus fly.
The default HMM parameters were previously hard-coded.
There are now other options available to tune the paramters.
This option is the first way to learn parameters.
It uses kmeans clustering of the data to estimate initial, transition, and emission probs.
For now, you need to make an assumption about k by providing an integer.
This model probably works best when you expect 2-3 states...
This option over-rides all other parameter options (which will be ignored).
''')

    parser_puffcn.add_argument('--iters', type=int, default=1,
                               help='''Number of iterations to run for updating parameters. Default: 1 (no updates).
''')
    parser_puffcn.add_argument('--converge', type=float, default=1e-9,
                               help='''When multiple iterations, stop iterating if difference between the log likelihood of the current state path is less than this much different than previous. Default = 1e-9.
''')
    parser_puffcn.add_argument('--learnpseudo', type=float, default=1e-323,
                               help='''When learning state transitions from previous state path, zero counts can throw a wrench into the spokes.
This prevents that. Can be 1e-323 to Inf, but recommend <= 1.
Default pseudocount is basically the lowest non-zero number possible: 1e-323..
''')
    parser_puffcn.add_argument('--emitpseudo', type=float, default=1e-7,
                               help='''Specifically for --model discrete.
When learning emission fequencies from previous state path, zero counts can throw a wrench into the spokes.
This prevents that. Can be 1e-323 to Inf, but recommend <= 1.
Default pseudocount: 1e-7 (1 in 10 million).
''')
    parser_puffcn.add_argument('--constrainEmit', action='store_true', default=False,
                               help='''When iterating, do not update the emission probabilities: constrain them to what was given to initialize.
''')
    parser_puffcn.add_argument('--outpfx', type=str, default=None,
                               help='''Prefix for output. If not used, all output goes to stdout.
This is particularly useful when iters > 1, which outputs the statepath from each round.
If an out prefix is not specified and iters > 1, you will be warned/reminded in a stderr message that can be ignored if purposeful.
''')

    parser_puffcn.set_defaults(func=run_subtool)






















    ## create sub-command for summits
    parser_summits = subparsers.add_parser('summits',
                                           help=''' Find summits...''')
    parser_summits.add_argument('-l','--latestage', type=str, required=True,
                                help='''Provide path to bedGraph (e.g. made from getcov) for a late stage sample.''')
    parser_summits.add_argument('-e','--earlystage', type=str, required=False, default=False,
                               help=''' Optional: Provide path to bedGraph (e.g. made from getcov) for an early stage sample. This is used after smoothing and median normalization to further normalize the late-stage sample (e.g. can correct for sequencing biases)''')

    parser_summits.add_argument('--replace', action='store_true', default=False,
                                        help='''Turn on "replace" functionality. By default this will replace '.' in the count column of bedGraphs with '0'.
Use --replace_with and --replace_this to change.''')
    parser_summits.add_argument('--replace_this', type=str, default='.',
                                        help='''Used with --replace. Specify the character in count column to replace. Default = '.' ''')
    parser_summits.add_argument('--replace_with', type=str, default='0',
                                        help='''Used with --replace. Specify the character to replace the --replace_this character with.
Must be a string that can be converted to a float. Default = '0' ''')



    parser_summits_protocol = parser_summits.add_mutually_exclusive_group(required=True)
    parser_summits_protocol.add_argument('-1', '--protocol1', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are median normalized.
Then late stage is normalized to early stage if available.
Then the HMM is run.''')
    parser_summits_protocol.add_argument('-2', '--protocol2', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first smoothed with bandwidth given by --bandwidth.
Then they are median normalized.
Then late stage is normalized to early stage if available.
Then the HMM is run.''')
    parser_summits_protocol.add_argument('-3', '--protocol3', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first median normalized.
Then they are smoothed with bandwidth given by --bandwidth.
Then late stage is normalized to early stage if available.
Then the HMM is run.
Note: if early is not present, this is same as protocol 4.''')
    parser_summits_protocol.add_argument('-4', '--protocol4', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first median normalized.
Then late stage is normalized to early stage if available.
Then late/early is smoothed with bandwidth given by --bandwidth.
Then the HMM is run.
Note: if early is not present, this is same as protocol 3.''')

    parser_summits_protocol.add_argument('-5', '--protocol5', action='store_true', default=False,
                                        help='''Late stage is normalized to early stage if available.
Then late/early is smoothed with bandwidth given by --bandwidth. (i.e. L/E -> smooth)
Then the HMM is run.
Note: if early is not present, this is same as protocol 6.''')

    parser_summits_protocol.add_argument('-6', '--protocol6', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first smoothed with bandwidth given by --bandwidth.
Then late stage is normalized to early stage if available. (i.e. smooth -> L/E)
Then the HMM is run.
Note: if early is not present, this is same as protocol 5.''')

    parser_summits_protocol.add_argument('-7', '--protocol7', action='store_true', default=False,
                                        help='''Late stage is normalized to early stage if available.
Note: No median normalization or smoothing is performed. If only late is given, then this is just an identity/pass-through function.''')

    parser_summits_protocol.add_argument('-8', '--protocol8', action='store_true', default=False,
                                        help='''SKEW. Only accepts one file (e.g. latestage).
                                                For file with N rows, returns N-1 rows.
                                                For each value, V[i], in 4th column, for i in 2:N, it computes Skew = (V[i]-V[i-1]) / (V[i]+V[i-1])''')

    parser_summits_protocol.add_argument('-9', '--protocol9', action='store_true', default=False,
                                        help='''PERCENT CHANGE. Only accepts one file (e.g. latestage).
                                                For file with N rows, returns N-1 rows.
                                                For each value, V[i], in 4th column, for i in 2:N, it computes Skew = 100*(V[i]-V[i-1]) / V[i-1]''')

    parser_summits_protocol.add_argument('-10', '--protocol10', action='store_true', default=False,
                                        help='''SKEW CHANGE or SKEW DERIVATIVE. Only accepts one file (e.g. latestage).
                                                For file with N rows, returns N-1 rows.
                                                For each value, V[i], in 4th column, for i in 2:N, it computes Skew = 100*(V[i]-V[i-1]) / V[i-1]''')

    parser_summits_protocol.add_argument('-11', '--protocol11', action='store_true', default=False,
                                        help='''PERCENT CHANGE DERIVATIVE. Only accepts one file (e.g. latestage).
                                                For file with N rows, returns N-1 rows.
                                                For each value, V[i], in 4th column, for i in 2:N, it computes Skew = 100*(V[i]-V[i-1]) / V[i-1]''')


    parser_summits.add_argument('--stringcols', action='store_true', default=False,
                               help='''Just treat columns other than 4 as strings...''')
    
    parser_summits.add_argument('--log2', action='store_true', default=False,
                               help='''Return log2 values. Default = False.''')
    parser_summits.add_argument('--log10', action='store_true', default=False,
                               help='''Return log10 values. Default = False.''')
                              
    parser_summits.add_argument('--endsmoothing', action='store_true', default=False,
                              help=''' Add smoothing to the absolute end of any of the protocols for more flexibility here. This comes after log-transformation steps, for example, which optionally comes at the end of any protocol.''')
   
   
    parser_summits.add_argument('-ps', '--pseudo', type=float, default=0.1,
                               help=''' Before normalizing late to early, add this pseudocount to all counts in order to avoid division by zero.
Should be between 0 and 1.
Should be small enough to not change other values much,
but big enough such that numbers divided by 0+pseudo do not become massive.
Default: 0.1.''')
    parser_summits.add_argument('-bw', '--bandwidth', type=int, default=2500,
                               help=''' If kernel smoothing, specify bandwidth (int).
Bandwidth should be bigger when no early stage normalization to try to smooth out sequencing biases, mappability biases, etc.
Default: 2500.''')
    parser_summits.add_argument('--impute', type=int, default=False,
                               help=''' If imputing, specify bandwidth (int) for  kernel smoothing.
This bandwidth is generally longer than the one you would provide for regular smoothing.
Only bins with a count of 0 will take on smoothed (imputed) values.
Try: 10000.''')
    
    parser_summits_regions = parser_summits.add_mutually_exclusive_group(required=True)

    parser_summits_regions.add_argument('--regions', type=str, default=False,
                                help = ''' Find summits in these regions - provide BED file.''')
    parser_summits_regions.add_argument('--states', type=str, default=False,
                                        help=''' Provide statepath bedGraph output by "cn" sub-command. Peak regions will be found automatically.''')
    parser_summits.add_argument('--thresh_state', type=int, default=1,
                                help=''' Used with --states. Only consider regions with states higher than state given. Default: 1.''')
    parser_summits.add_argument('--merge1', type=float, default=10e3,
                                help = '''Used with --states. After extracting only bins with higher state value than --thresh_state, merge bins if they are with in --merge1 bp from each other. Default: 10e3.''')
    parser_summits.add_argument('--minwidth', type=float, default=50e3,
                                help = '''After extracting bins with states > --thresh_state and merging remaining bins that are within --merge1 bp of each other,
only keep merged regions > --minwidth.''')
    parser_summits.add_argument('--merge2', type=float, default=40e3,
                                help = '''After (i) extracting bins with states > --thresh_state, (ii) merging remaining bins that are within --merge1 bp of each other,
(iii) retaining only merged regions > --minwidth, merge regions that are within --merge2 bp of each other.''')
    parser_summits.add_argument('--max_state_thresh', type=int, default=2,
                                help = ''' After (i) extracting bins with states > --thresh_state, (ii) merging remaining bins that are within --merge1 bp of each other,
(iii) retaining only merged regions > --minwidth, (iv) merging filtered regions that are within --merge2 bp of each other,
only retain the remaining merged regions if their maximum state is > --max_State_thresh''')

    
    parser_summits.set_defaults(func=run_subtool)
















    ## create sub-command for normalize
    parser_normalize = subparsers.add_parser('normalize',
                                          help = '''Given a latest-stage sample (where all or most puffs have grown) and an optional earliest stage sample
(for additional Fold-enrichment normalization), just return the late-stage sample with normalized values as specified by protocol options below.''')

    parser_normalize.add_argument('-l','--latestage', type=str, required=True,
                                help='''Provide path to bedGraph (e.g. made from getcov) for a late stage sample.''')
    parser_normalize.add_argument('-e','--earlystage', type=str, required=False, default=False,
                               help=''' Optional: Provide path to bedGraph (e.g. made from getcov) for an early stage sample. This is used after smoothing and median normalization to further normalize the late-stage sample (e.g. can correct for sequencing biases)''')
    parser_normalize.add_argument('--replace', action='store_true', default=False,
                                        help='''Turn on "replace" functionality. By default this will replace '.' in the count column of bedGraphs with '0'.
Use --replace_with and --replace_this to change.''')
    parser_normalize.add_argument('--replace_this', type=str, default='.',
                                        help='''Used with --replace. Specify the character in count column to replace. Default = '.' ''')
    parser_normalize.add_argument('--replace_with', type=str, default='0',
                                        help='''Used with --replace. Specify the character to replace the --replace_this character with.
Must be a string that can be converted to a float. Default = '0' ''')

    parser_normalize_protocol = parser_normalize.add_mutually_exclusive_group(required=True)
    parser_normalize_protocol.add_argument('-1', '--protocol1', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are median normalized.
Then late stage is normalized to early stage if available.
Then the HMM is run.''')
    parser_normalize_protocol.add_argument('-2', '--protocol2', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first smoothed with bandwidth given by --bandwidth.
Then they are median normalized.
Then late stage is normalized to early stage if available.
Then the HMM is run.''')
    parser_normalize_protocol.add_argument('-3', '--protocol3', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first median normalized.
Then they are smoothed with bandwidth given by --bandwidth.
Then late stage is normalized to early stage if available.
Then the HMM is run.
Note: if early is not present, this is same as protocol 4.''')
    parser_normalize_protocol.add_argument('-4', '--protocol4', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first median normalized.
Then late stage is normalized to early stage if available.
Then late/early is smoothed with bandwidth given by --bandwidth.
Then the HMM is run.
Note: if early is not present, this is same as protocol 3.''')

    parser_normalize_protocol.add_argument('-5', '--protocol5', action='store_true', default=False,
                                        help='''Late stage is normalized to early stage if available.
Then late/early is smoothed with bandwidth given by --bandwidth. (i.e. L/E -> smooth)
Then the HMM is run.
Note: if early is not present, this is same as protocol 6.''')

    parser_normalize_protocol.add_argument('-6', '--protocol6', action='store_true', default=False,
                                        help='''Late stage (and early stage if present) bin counts are first smoothed with bandwidth given by --bandwidth.
Then late stage is normalized to early stage if available. (i.e. smooth -> L/E)
Then the HMM is run.
Note: if early is not present, this is same as protocol 5.''')

    parser_normalize_protocol.add_argument('-7', '--protocol7', action='store_true', default=False,
                                        help='''Late stage is normalized to early stage if available.
Note: No median normalization or smoothing is performed. If only late is given, then this is just an identity/pass-through function.''')

    parser_normalize_protocol.add_argument('-8', '--protocol8', action='store_true', default=False,
                                        help='''SKEW. Only accepts one file (e.g. latestage).
                                                For file with N rows, returns N-1 rows.
                                                For each value, V[i], in 4th column, for i in 2:N, it computes Skew = (V[i]-V[i-1]) / (V[i]+V[i-1])''')

    parser_normalize_protocol.add_argument('-9', '--protocol9', action='store_true', default=False,
                                        help='''PERCENT CHANGE. Only accepts one file (e.g. latestage).
                                                For file with N rows, returns N-1 rows.
                                                For each value, V[i], in 4th column, for i in 2:N, it computes Skew = 100*(V[i]-V[i-1]) / V[i-1]''')

    parser_normalize_protocol.add_argument('-10', '--protocol10', action='store_true', default=False,
                                        help='''SKEW CHANGE or SKEW DERIVATIVE. Only accepts one file (e.g. latestage).
                                                For file with N rows, returns N-1 rows.
                                                For each value, V[i], in 4th column, for i in 2:N, it computes Skew = 100*(V[i]-V[i-1]) / V[i-1]''')

    parser_normalize_protocol.add_argument('-11', '--protocol11', action='store_true', default=False,
                                        help='''PERCENT CHANGE DERIVATIVE. Only accepts one file (e.g. latestage).
                                                For file with N rows, returns N-1 rows.
                                                For each value, V[i], in 4th column, for i in 2:N, it computes Skew = 100*(V[i]-V[i-1]) / V[i-1]''')

    parser_normalize_protocol.add_argument('-12', '--protocol12', action='store_true', default=False,
                                        help=''' Median ratio normalization. Late is normalized to early. Then those ratios are median normalized.
                                                This is similar to a global version of what's used in DEseq2 (or TMM for EdgeR).''')

    parser_normalize_protocol.add_argument('-13', '--protocol13', action='store_true', default=False,
                                        help=''' Median ratio normalization with pre-smoothing. Late (and early if present) is smoothed. Late is normalized to early. Then those ratios are median normalized.
                                                This is similar to a global version of what's used in DEseq2 (or TMM for EdgeR).''')

    parser_normalize_protocol.add_argument('-14', '--protocol14', action='store_true', default=False,
                                        help=''' Median ratio normalization with post-smoothing. Late is normalized to early. Then those ratios are smoothed. Then median normalized.
                                                This is similar to a global version of what's used in DEseq2 (or TMM for EdgeR).''')

    parser_normalize_protocol.add_argument('-15', '--protocol15', action='store_true', default=False,
                                        help=''' Median ratio normalization with end-smoothing. Late is normalized to early. Then those ratios are median normalized. Then smoothed.
                                                This is similar to a global version of what's used in DEseq2 (or TMM for EdgeR).''')

    parser_normalize_protocol.add_argument('-16', '--protocol16', action='store_true', default=False,
                                        help=''' Glocal Median ratio normalization. Late is normalized to early. Then those ratios are median normalized based on coverage value in late.
                                                This is similar to a coverage-based local version of what's used in DEseq2 (or TMM for EdgeR).''')
    

    parser_normalize_protocol.add_argument('-17', '--protocol17', action='store_true', default=False,
                                        help=''' First normalizes late stage to median (X/Med), and scales it with --scalecov to make the median equal scalecov.
                                                It does the same to early if present.
                                                Then Late is normalized to Early (FE) using given pseudocount (not pseudocount given has same ratio to median in both samples).
                                                It then does median ratio normalization.
                                                It then returns it as is or logged if --log10 or --log2 specified''')

    parser_normalize_protocol.add_argument('-18', '--protocol18', action='store_true', default=False,
                                        help='''Robust Z scores. If early (control) sample given, first is control norm (late/early, test/control) followed by robust Z. The fold-change calculation is not treated specially as it might be with other options. However, the robust Z should be the same or similar even if the t/c was not then median ratio normalized, for example, since that is just a scaling factor.''')

    parser_normalize_protocol.add_argument('-19', '--protocol19', action='store_true', default=False,
                                        help='''Rank scores. If early (control) sample given, first is control norm (late/early, test/control) followed by ranking. The fold-change calculation is not treated specially as it might be with other options. However, the rank should be the same or similar even if the t/c was not then median ratio normalized, for example, since that is just a scaling factor.''')


    parser_normalize_protocol.add_argument('-20', '--protocol20', action='store_true', default=False,
                                        help='''Robust Z score differences between two samples. Requires both late (test) and early (control) samples. Robust Z scores are calculated independently for each sample. Then R_control is subtracted from R_test = R_t - R_c.''')

    parser_normalize_protocol.add_argument('-21', '--protocol21', action='store_true', default=False,
                                        help='''Rank score differences between two samples. Requires both late (test) and early (control) samples. Rank scores are calculated independently for each sample. Then R_control is subtracted from R_test = R_t - R_c.''')

    parser_normalize_protocol.add_argument('-22', '--protocol22', action='store_true', default=False,
                                        help='''Signal Per Million Reads (Or Counts or Per Million whatever can be summed up in the 4th column). Bin_spmr = 1e6*Bin/Sum(Bins).
Use --SPXR to change scaling factor from 1e6 to X. If early (control) given, both are SPMR'd independently, then late/early (test/control).
When an early (control) sample is provided, you may also want to check the default pseudocount applied.''')

    parser_normalize_protocol.add_argument('-23', '--protocol23', action='store_true', default=False,
                                        help='''Rank standardize scores. First rank, then subtract and divide by middle: (r-M)/M, where r is a bins rank, and M is the theoretical middle rank: M=(min+max)/2. If early (control) sample given, first is control norm (late/early, test/control) followed by ranking. The fold-change calculation is not treated specially as it might be with other options. However, the rank should be the same or similar even if the t/c was not then median ratio normalized, for example, since that is just a scaling factor.''')

    parser_normalize_protocol.add_argument('-24', '--protocol24', action='store_true', default=False,
                                        help='''Pct difference from early (control): 100*(T-C)/abs(C). Usually for T and C values >= 0, but abs(C) allows both pos and neg values. This requires both test (late) and control (early) samples. It assumes samples are pre-prcoessed however you want. Other options below do pre-processing before this step.''')

    parser_normalize_protocol.add_argument('-25', '--protocol25', action='store_true', default=False,
                                        help='''Pct skew of late (test) vs early (control): 100*(T-C)/(abs(T)+abs(C)). Usually for T and C values >= 0, but (abs(T)+abs(C)) is an experimental way to allow both pos and neg values. This requires both test (late) and control (early) samples. It assumes samples are pre-prcoessed however you want. Other options below do pre-processing before this step.''')

    parser_normalize_protocol.add_argument('-26', '--protocol26', action='store_true', default=False,
                                        help='''Test Pct difference from early (control) after both test and control samples are transformed into Robust Z-scores: 100*(R_t-R_c)/abs(R_c). ''')


    parser_normalize_protocol.add_argument('-27', '--protocol27', action='store_true', default=False,
                                        help='''Pct skew given late (test) and early (control) samples, after both test and control samples are transformed into Robust Z-scores: 100*(R_t-R_c)/(abs(R_t)+abs(R_c)).''')


    parser_normalize_protocol.add_argument('-28', '--protocol28', action='store_true', default=False,
                                        help='''Test Pct difference from early (control) after both test and control samples are transformed into Ranks: 100*(R_t-R_c)/abs(R_c). ''')


    parser_normalize_protocol.add_argument('-29', '--protocol29', action='store_true', default=False,
                                        help='''Pct skew given late (test) and early (control) samples, after both test and control samples are transformed into Ranks: 100*(R_t-R_c)/(abs(R_t)+abs(R_c)).''')

    parser_normalize_protocol.add_argument('-30', '--protocol30', action='store_true', default=False,
                                        help='''Median ratio normalization. Late is normalized to early. Then those ratios are locally median normalized -- phycally local in the genome, with window size controlled by --halfwidth (default 10 bins to each side).
                                                This is similar to a what's used in DEseq2 (or TMM for EdgeR) -- on local genomic regions.
                                                A motivation to do this is to allow for non-linear corrections of ChIP vs control over a range of copy numbers, as is seen in intrachromosomal DNA amplification for example.
                                                The alternative is to use the global median ratio for everything, which assumes no distortion of ratios at different copy numbers.
                                                Note that if the experiment is to map relative copy numbers compared to a control DNA sample, the GLOBAL approach is what you would want to use.
                                                The local approach would make everything look like RCN=1 in that scenario.
                                                The local approach is aimed at eliminating effects on the ratios due to local biases such as copy number in order to leave only peaks due to ChIP, for example.''')


    parser_normalize_protocol.add_argument('-31', '--protocol31', action='store_true', default=False,
                                        help='''Chromosome-specific Median ratio normalization. Late is normalized to early. Then ratios on a given chromosome (sequence) are normalized to the median of that chromosome (sequence).
                                                This is similar to a what's used in DEseq2 (or TMM for EdgeR) -- on a chromosome by chromosome basis.
                                                A motivation to do this is that the global median is not always a safe bet for estimating the background of each chromosome or sequence in a file, especially if it is meta-genomic sample where each sequence may have a different copy number.
                                                An example from Sciara is that the mapping rate of the X chromosome is different from the autosomes for at least two scenarios: males that have one copy of the X, but two of each autosome; and X'X females that have one X and one variant of the X called X' (Xprime).
                                                In the latter scenario, either one is mapping X and X' data to the X if the X' sequence is not present and the X copy number may seem lower than autosomes when it is not; OR one is mapping to both X and X' sequences, putting each around haploid coverage compared to autosomes.
                                                In all these cases, the global median is not a good estimate of the background for the X chromosome.''')

    parser_normalize_protocol.add_argument('-32', '--protocol32', action='store_true', default=False,
                                        help='''Median smoothing. This is similar to protocol 30, but instead of median ratio normalization, it just smooths the values in the late stage sample (or late/early ratio early if present) with a local median defined by --halfwidth (default 10 bins to each side of a bin).
                                                Motivations to do this are any motivations for smoothing in general. This is useful for smoothing out local biases in the data.
                                                In particular, there are some spurious outlier bins that are not representative of the local region, and should be smoothed out.
                                                An option might be to set a cap on the maximum value in a bin, but that is not implemented yet.
                                                This is not a median ratio normalization.
                                                A desired pipeline might be: median ratio normalization with protocol 30 or 31, followed by local median smoothing with this protocol 32.''')


    parser_normalize_protocol.add_argument('-33', '--protocol33', action='store_true', default=False,
                                        help='''Trimmed Mean smoothing. Similar to protocol 32, but uses a trimmed mean instead of a median for smoothing.''')

    parser_normalize_protocol.add_argument('-34', '--protocol34', action='store_true', default=False,
                                        help='''Mean smoothing. Similar to protocol 32, but uses a mean instead of a median for smoothing.''')


    parser_normalize.add_argument('--stringcols', action='store_true', default=False,
                               help='''Just treat columns other than 4 as strings...''')

    parser_normalize.add_argument('--log2', action='store_true', default=False,
                               help='''Return log2 values. Default = False.''')
    parser_normalize.add_argument('--log10', action='store_true', default=False,
                               help='''Return log10 values. Default = False.''')
    parser_normalize.add_argument('--scalecov', type=float, default=1,
                               help='''Multiply coverage by this as part of protocol 17.''')    
    parser_normalize.add_argument('--SPXR', type=float, default=1e6,
                               help='''In essence, this is like --scalecov with a different default: 1e6.''')

    parser_normalize.add_argument('--halfwidth', type=int, default=10,
                               help='''In local operations (only protocol30, local med ratios atm), this is how many bins on each side of a central position is used to calculate a statistic.
                                        The total window size would therefore be:
                                            window = LHS + position + RHS = halfwidth + 1 + halfwidth = halfwidth*2 + 1.
                                        Default = 10.
                                        Note that the default has different spans depending on the input bin size.
                                        If using 500 bp bins, then 10 bins to each side equals 5 kb to each side (10.5 kb window), but just 1 kb (2.1 kb window) if using 100 bp bins.''') 


    # parser_normalize.add_argument('--safe', type=int, default=1,
    #                            help='''This only applies to protocols 30 and 32. When looking for local medians, if the median is 0, it tries to use trimmed mean; if that is 0, it tries to use mean; if that is 0, it returns the safe value.
    #                            This has been set to 1 by default, but can be set to 0 or any other value. 
    #                            For protocol 30, it makes sense to keep it as 1. This will result in those values being normalized to 1 (returned back as is, which was 0 in the first place).
    #                            For protocol 32, this safety should be turned off, since it is just smoothing, and 0s should be returned as 0s. Etc.''') 


    parser_normalize.add_argument('--pseudoZeroBins', action='store_true', default=False,
                               help='''Not to be confused with --pseudo. This option applies only to protocols 24-27 right now. It only need be used when there are zeros in the control (early) sample. In protocols 26 and 27, this is likely to happen from the robust z-score pre-processing. If an error is thrown, try --pseudoZeroBins or --addMinOtherPlusOneToBoth. --pseudoZeroBins adds min(abs(nonzero control values)) to bins in both samples that were 0 in contorl (early). --addMinOtherPlusOneToBoth shifts both distributions up by min(control values)+1, setting the min control value to 1. Both use minimum values for each chrom/contig independently rather than a global min. This is intended to reduce the effects of the modifications, but may introduce its own issues. These are not meant to be used together, but won't throw an error if they are.''')    
    parser_normalize.add_argument('--addMinOtherPlusOneToBoth', action='store_true', default=False,
                               help='''This option applies only to protocols 24-27 right now. It only need be used when there are zeros in the control (early) sample. In protocols 26 and 27, this is likely to happen from the robust z-score pre-processing. If an error is thrown, try --pseudoZeroBins or --addMinOtherPlusOneToBoth. --pseudoZeroBins adds min(abs(nonzero control values)) to bins in both samples that were 0 in contorl (early). --addMinOtherPlusOneToBoth shifts both distributions up by min(control values)+1, setting the min control value to 1. These are not meant to be used together, but won't throw an error if they are. NOTE: Both use minimum values for each chrom/contig independently rather than a global min. This is intended to reduce the effects of the modifications, but may introduce its own issues. ''')    
    parser_normalize.add_argument('--setToControlDist', action='store_true', default=False,
                               help='''This option applies only to protocols 24-27 right now. It resets the contorl RZ scores back to original, and scales the test R_z in same way. Z_c = (X_c-med_c)/mad_c ; X_c = mad_c * Z_c + med_c ;  X_t_c = mad_c * Z_t + med_c''')    


    parser_normalize.add_argument('-c', '--collapsed', action='store_true', default=False,
                               help='''Return collapsed variable-step bedGraph instead of expanded single-step bedGraph.
This is often a much smaller file.''')
    
    parser_normalize.add_argument('-ps', '--pseudo', type=float, default=0.1,
                               help=''' Before normalizing late to early, add this pseudocount to all counts in order to avoid division by zero.
Should be between 0 and 1.
Should be small enough to not change other values much,
but big enough such that numbers divided by 0+pseudo do not become massive.
Default: 0.1.''')
    
    parser_normalize.add_argument('-bw', '--bandwidth', type=int, default=2500,
                               help=''' If kernel smoothing, specify bandwidth (int).
Bandwidth should be bigger when no early stage normalization to try to smooth out sequencing biases, mappability biases, etc.
Default: 2500.''')
    parser_normalize.add_argument('--endsmoothing', action='store_true', default=False,
                               help=''' Add smoothing to the absolute end of any of the protocols for more flexibility here. This comes after log-transformation steps, for example, which optionally comes at the end of any protocol.''')
    parser_normalize.add_argument('--replaceNaN', action='store_true', default=False,
                               help=''' If kernel smoothing, NaNs can be generated. This option replaces those with local averages (see --localwinsize, default=5 bins). In cases where localaverages return NaN (very rare), it fills NaN with the global average for the given chrom/sequence (not whole geneome, so still local-ish).''')
    parser_normalize.add_argument('--localwinsize', type=int, default=5,
                               help=''' If kernel smoothing and/or using --replaceNan, this specifies the number of bins to use (centered on this so odd numbers preferable). Default = 5.''')
    parser_normalize.add_argument('--impute', type=int, default=False,
                               help=''' If imputing, specify bandwidth (int) for  kernel smoothing.
This bandwidth is generally longer than the one you would provide for regular smoothing.
Only bins with a count of 0 will take on smoothed (imputed) values.
Try: 10000.''')
        
##    parser_normalize.add_argument('--counts', type=str, default=False,
##                           help=''' Use this flag and specify an output prefix for the final normalized late stage bin counts bedGraph.''')


    parser_normalize.set_defaults(func=run_subtool)








###########GENERATE
        ## create sub-sommand for generate
    parser_generate = subparsers.add_parser('generate',
                                          help = '''Generate emitted_data and statepath bedGraphs.''')

    parser_generate.add_argument('-f','-b', '-i', '--bedgraph', type=str, required=True,
                                help='''Provide path to bedGraph that contains the intervals in first 3 columns to return with generated data.''')


    parser_generate.add_argument('-m', '--emodel', type=str, default='normal',
                               help='''Specify emissions model to assume for HMM. Options: normal, exponential. Default: normal.''')

    parser_generate.add_argument('--mu', type=str, default='1,2,4,8,16,32,64',
                               help=''' PuffCN has been optimized for mapping DNA puffs in the fungus fly.
The default state means were previously hard-coded.
This option allows some flexibility from the command-line to change the state means.
Default: 1,2,4,8,16,32,64
To change: Provide comma-seprated list of state means.
The number of states will be calculated from this list.
If changing state sigmas (used in normal model), it must have same number of states represented.

NOTE: If using exponential or geometric distribution, provide the expected mean RCN values of the states
    as you would for normal or poisson models. This script will automatically take their inverses to work
    in the exponential and geometric models.''')

    parser_generate.add_argument('--sigma', type=str, default=None,
                               help=''' PuffCN has been optimized for mapping DNA puffs in the fungus fly.
The default state sigmas (stdevs) were previously hard-coded.
This option allows some flexibility from the command-line to change the state sigmas.
Default: if not changed, defaults to square root of state means (Poisson-like).
To change: Provide comma-seprated list of state sigmas.
Alternatively: Use --mu_scale (default False) with a scaling factor multiplied against the MUs.
The number of states is calculated from this state mean list, which defaults to 7.
If changing state sigmas (used in normal model), it must have same number of states represented as state means.''')

    parser_generate.add_argument('--mu_scale', type=float, default=None,
                               help=''' See --sigma for more details on sigmas.
Use this to scale means (--mu) to use as stdevs (sigma) instead of taking square roots of means.
For example, --mu_scale 0.5 will use mu*0.5 as the stdev.''')


    parser_generate.add_argument('--special_idx', type=int, default=0,
                               help='''Only for use if you're very familiar with the program (and change defaults).
The default state means is 1,2,4,8,16,32,64.
The default index for the mean that represents copy number 1 is 0.
In this lingo - CN=1 is the special state, and the 0-based index of the special state in that list is 0.
If you were to change parameters that affect where the special state is in the list, make sure to change this index.
This index is only used to help construct initial probabilities and transition probabilies.
If understood, it can be used to designate any single special state (not necessarily the one that corresponds to CN=1).
The other parameters to use with this are:
--init_special (probability of starting in the special state (usually CN=1).
    The probabity of starting in a another state (usually copy number variant states) defaults to (1-init_special)/(nstates-1).
--prob_leave_special
--prob_stay_special
--prob_other_to_special
--prob_other_to_other
--prob_other_to_self

Alternative to, an initial probability vector can be given with --initialprobs 
''')

    parser_generate.add_argument('--init_special', type=float, default=0.997,
                               help='''Probability of starting in the 'special state' (usually copy number = 1). Default: 0.997.
The probabity of starting in a another state (usually copy number variant states) defaults to (1-init_special)/(nstates-1).
''')



    parser_generate.add_argument('--leave_special_state', type=float, default=0.001,
                               help='''Probability of leaving the 'special state' (usually copy number = 1).
Default: 0.001.
If number is betwen 0 and 1, it will be assumed a probability.
If number given is > 1, then it will be treated as the average length (number of bins) of the special state.
For example, if 1000 is given, it will be 1/1000 = 0.001.
In terms of bp lengths, one would need to multiply n_bins * bin_length OR divide bp length by bin_length
Thus, if you want to see a change every 500 kb w/ 500 bp bins, then 500kb/500 = 1 kb = 1000 -- which will be interpreted as 0.001.
Or as another example, if you expect to see a change every 2 Mb with 100 bp bins, then 2e6/1e2 = 1e4 = 10 kb = 10000, interpreted as 0.0001.

The probability of staying in this state is the complement: 1-p
''')

    parser_generate.add_argument('--leave_other', type=str, default=None,
                               help='''Probability of leaving one of the other states.
This defaults to --leave_special_state making all transition probabilities out of states the same (0.001 by default).

To change, provide a probability of leaving (p).

If the  first number is betwen 0 and 1, it will be assumed a probability.
If the first number given is > 1, then it will be treated as the average length (number of bins).
For example, if 1000 is given, it will be 1/1000 = 0.001.

If only 1 number is given, then that is assumed to be the probability of transitioning to all the other states.

You can also give a comma-separated pair of 2 probabilities:
    prob of leaving to special state
    prob of leaving to another 'non-special' state.
Make sure the probabilities sum to what you expect the overall probability of leaving the state is...
    which should be p_to_special + p_to_nonspecial * (num_non_special-1) = p_to_special + p_to_nonspecial (nstates-2)

For example, in a 7-state model:
    0.001,0.0002 --> 0.001 + 0.0002 * 5 = 0.001 + 0.001 = 0.002
    OR
    0.001,0.001 --> 0.001 + 0.001 * 5 = 0.006

If the second number is > 1, the same rules apply as to the first number.

For other analyses, I've used:
0.00001,0.000000000001
OR
0.001,0.0000000001

The probability of staying in these states is the complement: 1-p1-p2

NOTE: the program forces the transition probabilities of a given state to sum to 1.
''')

    parser_generate.add_argument('--initialprobs', type=str, default=None,
                               help='''PuffCN has been optimized for mapping DNA puffs in the fungus fly.
The default state means were previously hard-coded.
This option allows some flexibility from the command-line to change the state means.
Default: [0.997, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005]
The default will change with more or less states described w/ --mu and --sigma.
By default, the first state will start out as 0.997 as above, all other states will be (1-0.997)/n_other_states.
That behavior also changes with following parameters:
--special_idx -- determines which state (not necessarily first) will be given default 0.997 (OR other with --initcn1)
--init_special (probability of starting in the special state (usually CN=1).
    The probabity of starting in a another state (usually copy number variant states) defaults to (1-init_special)/(nstates-1).
--leave_special_state

--prob_other_to_special
--prob_other_to_other
--prob_other_to_self
To change the initial probs manually: Provide comma-separated list of initial probs -e.g.: '0.997,0.0005,0.0005,0.0005,0.0005,0.0005,0.0005'
This must have same number of states represented as state means (--mu; default 7).

''')

    
    parser_generate.set_defaults(func=run_subtool)





    ## create sub-command for filter/filterfish
    parser_filter = subparsers.add_parser('filter',
                                       help = '''Given a latest-stage sample (where all or most puffs have grown) and an optional earliest stage sample
    (for additional Fold-enrichment normalization), just return the late-stage sample with normalized values as specified by protocol options below.''')

    parser_filter.add_argument('--counts', type=str, default=False,
                            help=''' Use this flag and specify an output prefix for the final normalized late stage bin counts bedGraph.''')

    parser_filter_unit = parser_filter.add_mutually_exclusive_group()
    parser_filter_unit.add_argument('-sd1','--stdev_above', action='store_true', default=False, 
                             help='''Use value given as multiple of standard deviations above the mean.''')
    parser_filter_unit.add_argument('-sd2','--stdev_below', action='store_true', default=False, 
                             help='''Use value given as multiple of standard deviations BELOW the mean.''')
    parser_filter_unit.add_argument('-mu','--mean', action='store_true', default=False, 
                             help='''Use value given as multiple of the mean.''')
    parser_filter.add_argument('-V','--value', type=float, required=True,
                             help='''Value to filter on -- a float. Required.''')
    parser_filter.add_argument('-R','--relation', type=str, default=">",
                             help='''Relationship to value to filter on -- i.e. greater than, less than, etc. Accepted values are:
    gt, ge, lt, le, eq, ne -- respectively representing the relations >, >=, <, <=, ==, !=''')

    parser_filter.add_argument('-l','--latestage', type=str, required=True,
                             help='''Provide path to bedGraph (e.g. made from getcov) for a late stage sample.''')
    parser_filter.add_argument('-e','--earlystage', type=str, required=False, default=False,
                            help=''' Optional: Provide path to bedGraph (e.g. made from getcov) for an early stage sample. This is used after smoothing and median normalization to further normalize the late-stage sample (e.g. can correct for sequencing biases)''')

    parser_filter_protocol = parser_filter.add_mutually_exclusive_group(required=True)

    parser_filter_protocol.add_argument('-s','--skipnorm', action='store_true', default=False,
                             help='''Use provided bedGraph (late option) directly -- skip any normalization procedure.''')
    parser_filter_protocol.add_argument('-1', '--protocol1', action='store_true', default=False,
                                     help='''Late stage (and early stage if present) bin counts are median normalized.
    Then late stage is normalized to early stage if available..''')
    parser_filter_protocol.add_argument('-2', '--protocol2', action='store_true', default=False,
                                     help='''Late stage (and early stage if present) bin counts are first smoothed with bandwidth given by --bandwidth.
    Then they are median normalized.
    Then late stage is normalized to early stage if available.''')
    parser_filter_protocol.add_argument('-3', '--protocol3', action='store_true', default=False,
                                     help='''Late stage (and early stage if present) bin counts are first median normalized.
    Then they are smoothed with bandwidth given by --bandwidth.
    Then late stage is normalized to early stage if available.
    Note: if early is not present, this is same as protocol 4.''')
    parser_filter_protocol.add_argument('-4', '--protocol4', action='store_true', default=False,
                                     help='''Late stage (and early stage if present) bin counts are first median normalized.
    Then late stage is normalized to early stage if available.
    Then late/early is smoothed with bandwidth given by --bandwidth.
    Note: if early is not present, this is same as protocol 3.''')

    parser_filter_protocol.add_argument('-5', '--protocol5', action='store_true', default=False,
                                     help='''Late stage is normalized to early stage if available.
    Then late/early is smoothed with bandwidth given by --bandwidth. (i.e. L/E -> smooth).
    Note: if early is not present, this is same as protocol 6.''')

    parser_filter_protocol.add_argument('-6', '--protocol6', action='store_true', default=False,
                                     help='''Late stage (and early stage if present) bin counts are first smoothed with bandwidth given by --bandwidth.
    Then late stage is normalized to early stage if available. (i.e. smooth -> L/E).
    Note: if early is not present, this is same as protocol 5.''')

    parser_filter.add_argument('-c', '--collapsed', action='store_true', default=False,
                            help='''Return collapsed variable-step bedGraph instead of expanded single-step bedGraph.
    This is often a much smaller file.''')

    parser_filter.add_argument('-ps', '--pseudo', type=float, default=0.1,
                            help=''' Before normalizing late to early, add this pseudocount to all counts in order to avoid division by zero.
    Should be between 0 and 1.
    Should be small enough to not change other values much,
    but big enough such that numbers divided by 0+pseudo do not become massive.
    Default: 0.1.''')

    parser_filter.add_argument('-bw', '--bandwidth', type=int, default=2500,
                        help=''' If kernel smoothing, specify bandwidth (int).
    Bandwidth should be bigger when no early stage normalization to try to smooth out sequencing biases, mappability biases, etc.
    Default: 2500.''')
    parser_filter.add_argument('--impute', type=int, default=False,
                            help=''' If imputing, specify bandwidth (int) for  kernel smoothing.
    This bandwidth is generally longer than the one you would provide for regular smoothing.
    Only bins with a count of 0 will take on smoothed (imputed) values.
    Try: 10000.''')
     

    parser_filter.set_defaults(func=run_subtool)
 
 
 






    ## create sub-command for help
    parser_help = subparsers.add_parser('help', help=''' Gives more extensive guidance on using pufferfish.''')
    parser_help.set_defaults(func=run_subtool)


    
    ## parse the args and call the selected function
    args = parser.parse_args()


    


    ## check if args.quiet set (a default to all sub-modules from ArgumentParserWithDefaults class)
    if args.quiet:
        logger.setLevel(logging.ERROR)

    ## attempt to run args.func (which calls run_subtool() for all), catch errors
    try:
        args.func(parser, args)
    except IOError, e:
        ## often pipe will break for various reasons and will raise sigpipe error
        ## can import errno and do "if e.errno != errno.EPIPE:"
        ##   but errno.EPIPE = 32 -- so just use 32 here
        if e.errno != 32: ## ignore SIGPIPE
            raise


## Run main when this script used from command-line
if __name__ == "__main__":
    main()


## If needed parallelization
## Would have to parellize from here....
##    from joblib import Parallel, delayed
##    import time
##    from glob import glob
##    folder = "del"
##    files = glob('{}/*.txt'.format(folder))
##    def sleep(f):
##        print f
##        time.sleep(0.001)
##    ##    for f in files:
##    ##        sleep(f) #args.parallel
##    Parallel(n_jobs=2)(delayed(sleep)(f) for f in files)
