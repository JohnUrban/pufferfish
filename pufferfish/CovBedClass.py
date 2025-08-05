import os
from collections import defaultdict
import numpy as np
import pandas as pd
##np.seterr(divide='raise', invalid='raise')
np.seterr(divide='ignore', invalid='raise')
##np.seterr(divide='ignore', invalid='ignore')
import rpy2.robjects as robjects
ksmooth = robjects.r['ksmooth']
kmeans = robjects.r['kmeans']
intvec = robjects.IntVector
fltvec = robjects.FloatVector
matrixr = robjects.r.matrix
r = robjects.r
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage as stap
from puffR import *
puffR = stap(puffRstring, 'puffR')
import sys, datetime
from scipy.stats import spearmanr

import pandas as pd ## 2020 - potentially begin converting to pandas code where better
#from scipy.stats.mstats import winsorize



        


def s50(counts, stages, x=[50]):
        """
        counts,stages,x lists
        Returns sX for all x for a list of numbers "counts".
        Default: 50
        Assumes all values in list x are between 0 and 100.
        Interpretation: Returns stage at which >= 50% of reads is hit.
        """
        n = len(counts)
        x_to_stage = {e:0 for e in x}
        count_sum = sum(counts)
        total = 0
        i=0
        for e in sorted(x):
                target = count_sum*e/100.0
                while total < target and counts:
                        total += counts[i]
                        lastcount = counts[i]
                        laststage = stages[i]
                        i+=1
                try:
                    x_to_stage[e] = laststage
                except UnboundLocalError:
                    x_to_stage[e] = "."
        return x_to_stage


TMP_DIR = ".pufferfish_tmp_dir"
TMP_DIR = TMP_DIR[1:]

class CovBed(object):
    def __init__(self,
                 covbedfile,
                 count_only=False,
                 replace=False,
                 replace_with='0',
                 replace_this='.',
                 stringcols=False):
        ## "replace" means if you see the "replace_this" character in the count column, make it "replace_with"
        ## Made to deal with "." --> 0 by default when replace used.
        self.fopen = False
        self.connection = None
        self.file = covbedfile
        self.stringcols = stringcols
        self.start = {}
        self.end = {}
        self.count = {}
        self.chromosomes = set([])
        self.median = None
        self.mean = None
        self.sd = None
        self.mad = None
        self.rank = None
        self.rankstd = None
        self.sum = None
        self.nbins = None
        self.localmeds = None
        self.localmeans = None
        self.localtrimmedmeans = None
        self.count_only=count_only ## When False, start/end dicts are initialized, but remain empty: useful when comparing 2 bedgraphs of identical coords. See also MultiCovBed (though that currently requires a "stage file")
        self._extract_data(replace, replace_with, replace_this)
        
    def open(self):
        if self.fopen:
            self.close()        
        self.connection = open(self.file, 'r')
        self.fopen = True


    def close(self):
        self.connection.close()

    def _add_chromosome(self, chrom):
        self.chromosomes.add(chrom)
        self.start[chrom] = []
        self.end[chrom] = []
        self.count[chrom] = []

    def _update_data(self, chrom, start, end, count):
        if chrom not in self.chromosomes:
            self._add_chromosome(chrom)
        if not self.count_only: ## This allows the start/end dicts to be initialized, but remain empty
                self.start[chrom].append(start)
                self.end[chrom].append(end)
        self.count[chrom].append(count)

    def _finalize_data(self):
        ## convert set to list
        self.chromosomes = sorted(list(self.chromosomes))
        #convert lists to np arrays
        for chrom in self.chromosomes:
            self.start[chrom] = np.array(self.start[chrom])
            self.end[chrom] = np.array(self.end[chrom])
            self.count[chrom] = np.array(self.count[chrom])

                
    def _extract_data(self, replace=False, replace_with='0', replace_this='.'):
        self.open()
        for line in self.connection:
            chrom, start, end, count = line.strip().split()
            if replace and count == replace_this:
                count = replace_with
            if self.stringcols:
                self._update_data(chrom, start, end, float(count))
            else:
                self._update_data(chrom, int(start), int(end), float(count))
            ### JAN 9, 2018 -- I changed above int(float(count)) to float(count)
            ###         At this point, I don't know what it might break...
            ##          But now that I am using this more generally for signal data - not just counts - I need it as float
        self._finalize_data()
        self.close()


    def get_mean(self):
        if self.mean is None:
            counts = np.concatenate(self.count.values())
            self.mean = float(np.mean(counts))
            self.sd = float(np.std(counts))
        return self.mean

    def get_sd(self):
        if self.sd is None:
            counts = np.concatenate(self.count.values())
            self.mean = float(np.mean(counts))
            self.sd = float(np.std(counts,ddof=1))
        return self.sd

    def _get_median(self):
        counts = np.concatenate(self.count.values())
        self.median = float(np.median(counts))


    










    def _trimmed_mean_(self, x, extreme=0.25):
        ''' Returns the trimmed mean of x.
            Extreme is the proportion of values to trim from each end of sorted x.
            e.g. extreme=0.25 means 25% of values are trimmed from each end, so 50% of values are used to calculate the mean.
            If len(x) is 2 or 3, it returns the median of x.
            If len(x) < 2, it returns the only value in x.'''
        if len(x) < 2:
            return x[0]
        elif len(x) <= 3:
            return np.median(x)
        else:
            sort = sorted(x)
            l = len(sort)
            e = int(round(extreme*l))
            if e < (l-e):
                sort = sort[e:(l-e)] 
            ## This does have an else.....
            return np.mean(sort)


    def _safe_trimmed_mean_(self, x, unsafe=0, safe=1, meanFirst=True, extreme=0.25, play_it_safe=True):
        ''' if median is unsafe, make it safe.
            Optionally try chekcing and using safety of mean first'''
        ## Get trimmed mean of x.
        ans = self._trimmed_mean_(x, extreme=extreme)

        ## Some use cases may want to play with "unsafe" values, so this is an option (e.g. for median smoothing, there are no "unsafe" values).
        if not play_it_safe: ## i.e. if play_it_unsafe:
            ## RETURN!!! (Exit function here)
            return ans
        ## ELSE:
        if ans == unsafe and meanFirst:
            sys.stderr.write('Unsafe 2\n')
            ans = np.mean(x)
        ## Checks median or mean depending on safety of median and meanfirst option.
        if ans == unsafe:
            sys.stderr.write('Unsafe 3\n')
            ans = safe  ## This should be 1 if not log, 0 if log.
        return ans

    def _get_local_trimmed_means(self, halfwidth=10, play_it_safe=True, unsafe=0, safe=1, meanFirst=True, extreme=0.25):
        ''' halfwidth is the number of bins to each side.
            First position, A, is halfwidth
            Last position, Z, is L-halfwidth
            Positions 0 to A get local trimmed mean from 0-A+halfwidth+1.
            Positions Z+1 to L get local trimmed mean from Z-halfwidth to L.
            All others get position, P-haldwidth to P+halfwidth+1.
            Note: when play_it_safe=True, if the trimmed mean is the unsafe value, it is changed to a safe value.
                  when meanfirst=True, it will check if the mean gives a safe value, and use that.
                  Use case: local trimmed mean normalization, but 0s will give ZeroDivError.
                    Pseudocounts are sometimes used to prevent this, as well as eliminating 0 bins.
                    This is just a catch/check.
                  For trimming, the extreme proportion is substracted from both ends of sorted values before taking mean. (e.g. mean of bins 5:14 inclusive for 25% of 20 bins)'''
        self.localtrimmedmeans = {}
        for chrom in self.chromosomes:
            counts = np.array(self.count[chrom])
            L = len(self.count[chrom])
            A = halfwidth
            Z = L-halfwidth
            self.localtrimmedmeans[chrom] = np.zeros(L)
            if halfwidth > L:
                # Entire contig is just locally trimmed mean'd
                self.localtrimmedmeans[chrom][0:L] = self._safe_trimmed_mean_( x = counts,
                                                                 unsafe=unsafe, safe=safe,
                                                                 meanFirst=meanFirst, extreme=extreme, 
                                                                 play_it_safe=play_it_safe)
            else:
                # Init positions
                self.localtrimmedmeans[chrom][0:A] = self._safe_trimmed_mean_( x = counts[0:(A+halfwidth+1)],
                                                                 unsafe=unsafe, safe=safe,
                                                                 meanFirst=meanFirst, extreme=extreme, 
                                                                 play_it_safe=play_it_safe )
                # End positions
                self.localtrimmedmeans[chrom][(Z+1):L] = self._safe_trimmed_mean_( x = counts[(Z-halfwidth):L],
                                                                 unsafe=unsafe, safe=safe,
                                                                 meanFirst=meanFirst, extreme=extreme, 
                                                                 play_it_safe=play_it_safe )        
                if Z >= A:
                    ''' If Z < A, then entire contig has already been median norm'd with end positions overwriting overlapping start positions.
                    If Z > A, then there is at least one position in the middle.'''
                    # Middle positions:
                    for i in range(A,Z+1):
                        self.localtrimmedmeans[chrom][i] =  self._safe_trimmed_mean_( x = counts[(i-halfwidth):(i+halfwidth+1)],
                                                                         unsafe=unsafe, safe=safe,
                                                                         meanFirst=meanFirst, extreme=extreme, 
                                                                         play_it_safe=play_it_safe )



    def _get_local_means(self, halfwidth=10):
        ''' halfwidth is the number of bins to each side.
            First position, A, is halfwidth
            Last position, Z, is L-halfwidth
            Positions 0 to A get local trimmed mean from 0-A+halfwidth+1.
            Positions Z+1 to L get local trimmed mean from Z-halfwidth to L.
            All others get position, P-haldwidth to P+halfwidth+1.
            '''
        self.localmeans = {}
        for chrom in self.chromosomes:
            counts = np.array(self.count[chrom])
            L = len(self.count[chrom])
            A = halfwidth
            Z = L-halfwidth
            self.localmeans[chrom] = np.zeros(L)
            if halfwidth > L:
                # Entire contig is just mean'd
                self.localmeans[chrom][0:L] = np.mean(counts)
            else:
                # Init positions
                self.localmeans[chrom][0:A] = np.mean(counts[0:(A+halfwidth+1)])

                # End positions
                self.localmeans[chrom][(Z+1):L] = np.mean(counts[(Z-halfwidth):L])
                                                                                
                if Z >= A:
                    ''' If Z < A, then entire contig has already been median norm'd with end positions overwriting overlapping start positions.
                    If Z > A, then there is at least one position in the middle.'''
                    # Middle positions:
                    for i in range(A,Z+1):
                        self.localmeans[chrom][i] =  np.mean( counts[(i-halfwidth):(i+halfwidth+1)])





















    def _safe_median_(self, x, unsafe=0, safe=1, trimmedMeanFirst=True, meanFirst=True, extreme=0.25, play_it_safe=True):
        ''' if median is unsafe, make it safe.
            Optionally try chekcing and using safety of mean first'''
        ## Get median of x.
        ans = np.median(x)

        ## Some use cases may want to play with "unsafe" values, so this is an option (e.g. for median smoothing, there are no "unsafe" values).
        if not play_it_safe: ## i.e. if play_it_unsafe:
            ## RETURN!!! (Exit function here)
            return ans
        
        ## ELSE:
        ## Other use cases need to ensure that the "median" is safe (e.g. for local median ratio normalization, we try to avoid dividing by 0).
        if ans == unsafe and trimmedMeanFirst:
            sys.stderr.write('Unsafe 1\n')
            #sort = sorted(ans) ## This was here for some reason but ans is np.float64, not list... did I mean to sort x?
            sort = sorted(x)    ## 
            l = len(sort)
            if l > 3: ## other wise, depending on extreme it will exclude none, or sides (equiv to median), or all
                e = int(round(extreme*l))
                if e < (l-e):
                    sort = sort[e:(l-e)] 
            ans = np.mean(sort)
        if ans == unsafe and meanFirst:
            sys.stderr.write('Unsafe 2\n')
            ans = np.mean(x)
        ## Checks median or mean depending on safety of median and meanfirst option.
        if ans == unsafe:
            sys.stderr.write('Unsafe 3\n')
            ans = safe  ## This should be 1 if not log, 0 if log.
        return ans

    def _get_local_medians(self, halfwidth=10, play_it_safe=True, unsafe=0, safe=1, trimmedMeanFirst=True, meanFirst=True, extreme=0.25):
        ''' halfwidth is the number of bins to each side.
            First position, A, is halfwidth
            Last position, Z, is L-halfwidth
            Positions 0 to A get local median from 0-A+halfwidth+1.
            Positions Z+1 to L get local median from Z-halfwidth to L.
            All others get position, P-haldwidth to P+halfwidth+1.
            Note: when play_it_safe=True, if the median is the unsafe value, it is changed to a safe value.
                  when meanfirst=True and/or trimmedMeanFirst=True, it will check if the mean and/or trimmedMean gives a safe value, and use that.
                  trimmed mean takes precedent over mean; mean takes precedent over default safe value.
                  Use case: you likely plan to do local median normalization, but 0s will give ZeroDivError.
                    Pseudocounts are sometimes used to prevent this, as well as eliminating 0 bins.
                    This is just a catch/check.
                  For trimming, the extreme proportion is substracted from both ends of sorted values before taking mean. (e.g. mean of bins 5:14 inclusive for 25% of 20 bins)'''
        self.localmeds = {}
        for chrom in self.chromosomes:
            counts = np.array(self.count[chrom])
            L = len(self.count[chrom])
            A = halfwidth
            Z = L-halfwidth
            self.localmeds[chrom] = np.zeros(L)
            if halfwidth > L:
                # Entire contig is just locally median norm'd
                self.localmeds[chrom][0:L] = self._safe_median_( x = counts,
                                                                 unsafe=unsafe, safe=safe,
                                                                 trimmedMeanFirst=trimmedMeanFirst,
                                                                 meanFirst=meanFirst, extreme=extreme, 
                                                                 play_it_safe=play_it_safe)
            else:
                # Init positions
                self.localmeds[chrom][0:A] = self._safe_median_( x = counts[0:(A+halfwidth+1)],
                                                                 unsafe=unsafe, safe=safe,
                                                                 trimmedMeanFirst=trimmedMeanFirst,
                                                                 meanFirst=meanFirst, extreme=extreme, 
                                                                 play_it_safe=play_it_safe )
                # End positions
                self.localmeds[chrom][(Z+1):L] = self._safe_median_( counts[(Z-halfwidth):L],
                                                                 unsafe=unsafe, safe=safe,
                                                                 trimmedMeanFirst=trimmedMeanFirst,
                                                                 meanFirst=meanFirst, extreme=extreme, 
                                                                 play_it_safe=play_it_safe )
                if Z >= A:
                    ''' If Z < A, then entire contig has already been median norm'd with end positions overwriting overlapping start positions.
                    If Z > A, then there is at least one position in the middle.'''
                    # Middle positions:
                    for i in range(A,Z+1):
                        self.localmeds[chrom][i] =  self._safe_median_( x = counts[(i-halfwidth):(i+halfwidth+1)],
                                                                         unsafe=unsafe, safe=safe,
                                                                         trimmedMeanFirst=trimmedMeanFirst,
                                                                         meanFirst=meanFirst, extreme=extreme, 
                                                                         play_it_safe=play_it_safe )

    def _get_mad(self, relearn=False):
        if self.median is None or relearn:
            self._get_median()
        counts = np.concatenate(self.count.values())
        absdiffs = np.abs((counts - self.get_median()))
        self.mad = float(np.median(absdiffs))

    def _get_sum(self):
        counts = np.concatenate(self.count.values())
        self.sum = float(np.sum(counts))

    def _get_nbins(self):
        self.nbins = len(np.concatenate(self.count.values()))

    def _rank_data(self):
        counts = np.concatenate(self.count.values())
        ranks = np.array(pd.Series(counts).rank())
        assert len(counts) == len(ranks)
        self.rankdict = {counts[i]:ranks[i] for i in range(len(counts))}
        self.rank = {}
        for chrom in self.count.keys():
            self.rank[chrom] = [self.rankdict[e] for e in self.count[chrom]]

    def _rank_standardize_data(self):
        counts = np.concatenate(self.count.values())
        ranks = np.array(pd.Series(counts).rank())
        assert len(counts) == len(ranks)
        highest = self.get_nbins()
        lowest = 1.0
        M = (lowest+highest)/2.0
        self.rankdict = {counts[i]:((ranks[i]-M)/M) for i in range(len(counts))}
        self.rankstd = {}
        for chrom in self.count.keys():
            self.rankstd[chrom] = [self.rankdict[e] for e in self.count[chrom]]
        

        
##        counts = []
##        for chrom in self.chromosomes:
##            counts.append(self.count[chrom])
##        self.median = float(np.median(counts))
        
    def get_median(self, relearn=False):
        if self.median is None or relearn:
            self._get_median()
        return self.median

    def get_mad(self, relearn=False):
        if self.mad is None or relearn:
            self._get_mad()
        return self.mad

    def get_nbins(self):
        if self.nbins is None:
            self._get_nbins()
        return self.nbins

    def median_normalize_x(self, x, relearn=False):
        #x is np.array
        return x/self.get_median(relearn=relearn)

    def robust_z_normalize_x(self, x, relearn=False):
        #x is np.array
        return (x-self.get_median(relearn=relearn))/self.get_mad(relearn=relearn)

    def median_normalize_data(self, relearn=False):
        if relearn:
            self._get_median()
        for chrom in self.chromosomes:
            self.count[chrom] = self.median_normalize_x(self.count[chrom])

    def local_median_normalize_data(self, halfwidth=10, relearn=False):
        ''' halfwidth is the number of bins to each side.
            First position, A, is halfwidth
            Last position, Z, is L-halfwidth
            Positions 0 to A get local median from 0-A+halfwidth+1.
            Positions Z+1 to L get local median from Z-halfwidth to L.
            All others get position, P-haldwidth to P+halfwidth+1.'''
        
        if relearn or self.localmeds is None:
            self._get_local_medians(halfwidth=halfwidth)
        for chrom in self.chromosomes:
            self.count[chrom] = np.array(self.count[chrom]) / self.localmeds[chrom]

    def chromosome_median_normalize_data(self):
        ## DEVELOPING ON JULY 17, 2025; STILL NEEDS TESTING.
        '''Normalize the counts for a given chromosome by the median of that chromosome.'''
        for chrom in self.chromosomes:
            ## Get median count for this chromosome
            median_count = np.median(self.count[chrom])
            if median_count > 0:
                ## Normalize the counts for this chromosome by the median count
                self.count[chrom] = self.count[chrom] / median_count
            else:
                ## If median count is 0, Try the mean.
                mean_count = np.mean(self.count[chrom])
                if mean_count > 0:
                    self.count[chrom] = self.count[chrom] / mean_count
                else:
                    ## CAREFUL HERE: It is silently setting all counts to 1 to avoid division by zero in downstream analyses; but this may not be the desired behavior (e.g. if input is raw, raw ratio, or log2 ratio would have diff interpretations).
                    ## If both median and mean are 0, set all counts to 1; then they will be normalized to 1.
                    self.count[chrom] = np.ones_like(self.count[chrom])


    def local_median_smooth_data(self, halfwidth=10, relearn=False):
        ''' halfwidth is the number of bins to each side.
            First position, A, is halfwidth
            Last position, Z, is L-halfwidth
            Positions 0 to A get local median from 0-A+halfwidth+1.
            Positions Z+1 to L get local median from Z-halfwidth to L.
            All others get position, P-haldwidth to P+halfwidth+1.'''
        
        if relearn or self.localmeds is None:
            ## Get local medians around each bin in regions defined by halfwidth.
            ## Play_it_safe is False, so it will not attempt to change the median to a non-zero value (a useful behavior for local median ratio normalization, but not for local median smoothing).
            self._get_local_medians(halfwidth=halfwidth, play_it_safe=False) 
        for chrom in self.chromosomes:
            self.count[chrom] = self.localmeds[chrom]

    def local_mean_smooth_data(self, halfwidth=10, relearn=False):
        ''' halfwidth is the number of bins to each side.
            First position, A, is halfwidth
            Last position, Z, is L-halfwidth
            Positions 0 to A get local mean from 0-A+halfwidth+1.
            Positions Z+1 to L get local mean from Z-halfwidth to L.
            All others get position, P-haldwidth to P+halfwidth+1.'''
        if relearn or self.localmeans is None:
            ## Get local means around each bin in regions defined by halfwidth.
            self._get_local_means(halfwidth=halfwidth) 
        for chrom in self.chromosomes:
            self.count[chrom] = self.localmeans[chrom]

    def local_trimmed_mean_smooth_data(self, halfwidth=10, relearn=False):
        ''' halfwidth is the number of bins to each side.
            First position, A, is halfwidth
            Last position, Z, is L-halfwidth
            Positions 0 to A get local mean from 0-A+halfwidth+1.
            Positions Z+1 to L get local mean from Z-halfwidth to L.
            All others get position, P-haldwidth to P+halfwidth+1.'''
        if relearn or self.localtrimmedmeans is None:
            ## Get local trimmed means around each bin in regions defined by halfwidth.
            ## Play_it_safe is False, so it will not attempt to change the mean to a non-zero value (a useful behavior for local mean ratio normalization, but not for local mean smoothing).
            self._get_local_trimmed_means(halfwidth=halfwidth, play_it_safe=False) 
        for chrom in self.chromosomes:
            self.count[chrom] = self.localtrimmedmeans[chrom]


    def robust_z_normalize_data(self, relearn=False):
        if relearn:
            self._get_median()
            self._get_mad()
            
        for chrom in self.chromosomes:
            self.count[chrom] = self.robust_z_normalize_x(self.count[chrom])

    def rank_normalize_data(self, relearn=False):
        if self.rank is None or relearn:
            self._rank_data()
        for chrom in self.chromosomes:
            self.count[chrom] = self.rank[chrom]

    def rank_standardize_data(self, relearn=False):
        if self.rankstd is None or relearn:
            self._rank_standardize_data()
        for chrom in self.chromosomes:
            self.count[chrom] = self.rankstd[chrom]


    def spxr_normalize_data(self, x=1e6, relearn=False):
        ''' '''
        if self.sum is None or relearn:
            self._get_sum()
        scale_factor = float(x) / self.sum
        self.scale_data(scale=scale_factor)


    def scale_data(self, scale=1):
        for chrom in self.chromosomes:
            self.count[chrom] = scale*self.count[chrom]

    def log2_data(self, scale=1):
        for chrom in self.chromosomes:
            self.count[chrom] = np.log2(self.count[chrom])

    def log10_data(self, scale=1):
        for chrom in self.chromosomes:
            self.count[chrom] = np.log10(self.count[chrom])

            
    def expanded_bdg(self, bdg):
        ##bdg is just what should be in the 4th column
        string = ''
        for chrom in self.chromosomes:
            for i in range(len(self.start[chrom])):
                string += ('\t').join([chrom, str(self.start[chrom][i]),  str(self.end[chrom][i]),  str(bdg[chrom][i])]) + "\n"
        return string

    def expanded_bdg_two_cols(self, bdg1, bdg2):
        string = ''
        for chrom in self.chromosomes:
            for i in range(len( self.start[chrom] )):
                string += ('\t').join([chrom, str(self.start[chrom][i]), str(self.end[chrom][i]), str(bdg1[chrom][i]), str(bdg2[chrom][i])]) + "\n"
        return string
    
    def collapsed_bdg(self, bdg):
        ##bdg is just what should be in the 4th column
        string = ''
        for chrom in self.chromosomes:
            if len(self.start[chrom]) > 1:
                #init
                start = self.start[chrom][0]
                value = bdg[chrom][0]
                for i in range(1, len(self.start[chrom]) ):
                    if bdg[chrom][i] != value:
                        string += ('\t').join([chrom, str(start), str(self.end[chrom][i-1]), str(value)]) + "\n"
                        start = self.start[chrom][i]
                        value = bdg[chrom][i]
                ##finish chrom
                string += ('\t').join([chrom, str(start), str(self.end[chrom][i]), str(value)]) + "\n"
            else: #only 1 bin (very tiny contig)
                string += ('\t').join([chrom, str(self.end[chrom][0]), str(self.end[chrom][0]), str(bdg[chrom][0])]) + "\n"
        return string


    def get_bdg(self, bdg, collapsed=False):
        if not collapsed:
            return self.expanded_bdg(bdg)
        else:
            return self.collapsed_bdg(bdg)
    



    def filtered_bdg(self, relation = ">", value = 0, bdg=None):
        ##bdg is just what should be in the 4th column
        ## for this might typically be self.count
        if bdg is None:
            bdg = self.count
            string = ''
        if relation == "gt":
            keep = lambda x: x > value
        elif relation == "ge":
            keep = lambda x: x >= value
        elif relation == "lt":
            keep = lambda x: x < value
        elif relation == "le":
            keep = lambda x: x <= value
        elif relation == "eq":
            keep = lambda x: x == value
        elif relation == "ne":
            keep = lambda x: x != value
        for chrom in self.chromosomes:
            for i in range(len(self.start[chrom])):
                if keep(bdg[chrom][i]):
                    string += ('\t').join([chrom, str(self.start[chrom][i]),  str(self.end[chrom][i]),  str(bdg[chrom][i])]) + "\n"
        return string
    
    def __str__(self):
        return self.get_bdg(self.count)
    
    def get_chromosomes(self):
        return self.chromosomes

    def get_start_dict(self):
        return self.start
    def get_end_dict(self):
        return self.end
    def get_count_dict(self):
        return self.count

    def ksmooth_counts(self, bw=10000, rescueNaN=False, localWindow=5):
        for chrom in self.chromosomes:
            x = self.start[chrom]
            y = self.count[chrom]
            k = ksmooth(x = fltvec(x), y = fltvec(y), bandwidth = bw)
            self.count[chrom] = np.array(k[1])
            ## RESCUE NANs
            if rescueNaN:
                # 1. Compute global mean of self.count[chrom]
                mu = np.nanmean(self.count[chrom])
                # 2. Calculate local means (pd and convert to np)
                new = np.array(pd.DataFrame(self.count[chrom]).rolling(3, center=True, axis=0, min_periods=1).mean().transpose())[0]
                # 3. Replace any NaN in new with global mean
                new = np.where(np.isnan(new),mu,new)
                # 4. Replace any NaN in self.count[chrom] w/ corresponding value in new
                final = np.where(np.isnan(self.count[chrom]),new,self.count[chrom])
                # 5. Overwrite self.count[chrom] - I know this can be done in #4,  but felt like keeping it separate.
                self.count[chrom] = final


    def computeSkew(self):
        ''' Converts values, V, in counts to skew = (V[i]-V[i-1]) / (V[i]+V[i-1]).
        For i in 2:N. The first element is 0.'''
        for chrom in self.chromosomes:
            n = len(self.count[chrom])
            Next = self.count[chrom][1:n]
            Prev = self.count[chrom][:n-1]
            Diff = (Next - Prev)*100.0 ## 100.0 to ensure floats
            Sum = (Next + Prev)*1.0 ## 1.0 to ensure floats
            Skew = [0.0] + list(Diff / Sum)
            self.count[chrom] = np.array(Skew)

    def computePercentChange(self):
        ''' Converts values, V, in counts to skew = (V[i]-V[i-1]) / V[i-1].
        For i in 2:N. The first element is 0.'''
        for chrom in self.chromosomes:
            n = len(self.count[chrom])
            Next = self.count[chrom][1:n]
            Prev = self.count[chrom][:n-1]
            Diff = (Next - Prev)*100.0 ## 100.0 to ensure floats, and make as Pct
            PctChange = [0.0] + list(Diff / Prev)
            self.count[chrom] = np.array(PctChange)

    def computeSkewChange(self):
        self.computeSkew()
        for chrom in self.chromosomes:
            n = len(self.count[chrom])
            Next = self.count[chrom][1:n]
            Prev = self.count[chrom][:n-1]
            Diff = (Next - Prev)*1.0 ## 1.0 to ensure floats
            PctChange = [0.0] + list(Diff / 200.0) ## /200 b/c Skew was mult by 100 and I want to divide by 2 (similar to Hyrien segmentaiton of RFD)
            self.count[chrom] = np.array(PctChange)

    def computePercentChangeDerivative(self):
        ''' Converts values, V, in counts to skew = (V[i]-V[i-1]) / V[i-1].
        For i in 2:N. The first element is 0.'''
        self.computePercentChange()
        self.computePercentChange()
            
    def normalize_to_other(self, other, pseudocount=0.01):
        #other is another CovBed object with same bins from same genome
        for chrom in self.chromosomes:
            #print chrom
            self.count[chrom] = (np.array(self.count[chrom])+pseudocount)/(np.array(other.count[chrom])+pseudocount)


    def _opt_handle_zero_bins(self, other, chrom, setToControlDist=False, pseudoZeroBins=False, addMinOtherPlusOneToBoth=False):
        t = np.array(self.count[chrom])
        c = np.array(other.count[chrom])
        if setToControlDist:
            ## Assumes we are working with Z scores
            c = other.get_mad() * c + other.get_median()
            t = other.get_mad() * t + other.get_median()
        if pseudoZeroBins:
            ## slightly modify bins that have 0 in control (to avoid division)
            g = c == 0
            ng = c != 0
            m = np.abs(c[ng]).min() ## /10.0
            t[g] = t[g]+m
            c[g] = c[g]+m
        if addMinOtherPlusOneToBoth:
            ## Shift entire distro up -- not meant to be used with pseudoZeroBins, but won't throw error either.
            m = np.abs(c.min()) + 1
            t = t + m
            c = c + m
        return t, c
    
    def pct_diff_from_other(self, other, setToControlDist=False, pseudoZeroBins=False, addMinOtherPlusOneToBoth=False):
        #other is another CovBed object with same bins from same genome
        #I'm not supporting pseudo counts at this point -- in favor of removing 0 bins from both samples
        for chrom in self.chromosomes:
            t, c = self._opt_handle_zero_bins(other, chrom, setToControlDist, pseudoZeroBins, addMinOtherPlusOneToBoth)
            self.count[chrom] = 100.0*(t-c)/c

                

    def pct_skew_given_other(self, other, setToControlDist=False, pseudoZeroBins=False, addMinOtherPlusOneToBoth=False):
        #other is another CovBed object with same bins from same genome
        #I'm not supporting pseudo counts at this point -- in favor of removing 0 bins from both samples
        for chrom in self.chromosomes:
            t, c = self._opt_handle_zero_bins(other, chrom, setToControlDist, pseudoZeroBins, addMinOtherPlusOneToBoth)
            ###self.count[chrom] = 100.0*(np.array(self.count[chrom]) - np.array(other.count[chrom]))/(np.abs(np.array(self.count[chrom])) + np.abs(np.array(other.count[chrom])))
            self.count[chrom] = 100.0*(t - c)/(np.abs(t) + np.abs(c))
                
                
                

    def subtract_other(self, other, pseudocount=0.01):
        #other is another CovBed object with same bins from same genome
        for chrom in self.chromosomes:
            #print chrom
            self.count[chrom] = np.array(self.count[chrom]) - np.array(other.count[chrom])


    def normalize_with_glocalMedRatioNorm(self, other=None, pseudocount=0.01, globalweight=1, minlocalbins=3, minpropdata=0.1):
        #other is another CovBed object with same bins from same genome

        covd = self.create_local_medRatioNorm_dict(other=other,
                                              pseudocount=pseudocount)
        #print 'globalweight', globalweight
        covmeds = self.get_params_from_local_medRatioNorm_dict_with_glocalweighting(covd=covd,
                                                                               globalweight=globalweight,
                                                                               minlocalbins=minlocalbins,
                                                                               minpropdata=minpropdata)
        for chrom in self.chromosomes:
            norms = np.array(map(lambda x: covmeds[x], self.count[chrom]))
            if other is not None:
                self.count[chrom] = (1.0/norms)*(np.array(self.count[chrom])+pseudocount)/(np.array(other.count[chrom])+pseudocount)
            else:
                self.count[chrom] = norms*(np.array(self.count[chrom])+pseudocount)
        print("INDEV")
        #print(covd)
        #print()
        #print(covmeds)
        #print()
        #return (covd, covmeds)


    def impute_zeros(self, bw):
        '''When requiring mapq to be stringent, it leaves stretches of 0 that could benefit from being imputed.
        The 0 score often causes a state change in HMMs and can also lead to very inflated scores in FE after pseudocount added (if numerator is non-zero - e.g. 10/0.1 = 100)
        So this smooths the counts... and only uses the resulting smoothed values to substitute 0s.
        This means in very long 0 regions (e.g. contigs with no coverage), the score will remain 0 as desired.'''
        for chrom in self.chromosomes:
            x = self.start[chrom]
            y = self.count[chrom]
            k = puffR.impute_zeros(x = fltvec(x), y = fltvec(y), bw = bw)
            self.count[chrom] = np.array(k)
        
    def global_winsorize_counts(self, floor=0.00001, ceiling=0.99999):
        '''Set anything below the value of the floor quantile to that value.
           Set anything above the value of the ceiling quantile to that value.'''
        counts = np.concatenate(self.count.values())
        wange = np.quantile(counts, [floor, ceiling])
        for chrom in self.chromosomes:
            self.count[chrom][self.count[chrom] < wange[0]] = wange[0]
            self.count[chrom][self.count[chrom] > wange[1]] = wange[1]

    def local_winsorize_counts(self, floor=0.0001, ceiling=0.9999):
        '''Set anything below the value of the floor quantile to that value.
           Set anything above the value of the ceiling quantile to that value.'''
        counts = np.concatenate(self.count.values())
        gwange = np.quantile(counts, [floor, ceiling])
        for chrom in self.chromosomes:
            counts = np.array(self.counts[chrom])
            cwange = np.quantile(counts, [floor, ceiling])
            pass


    def create_local_medRatioNorm_dict(self, other=None, pseudocount=0.0):
        #other is another CovBed object with same bins from same genome
        covd = defaultdict(list)
        for chrom in self.chromosomes:
            if other is not None:
                ratios = (np.array(self.count[chrom])+pseudocount)/(np.array(other.count[chrom])+pseudocount)
            else:
                ratios = np.array(self.count[chrom])+pseudocount
            for i in range(len(self.count[chrom])):
                latecov = self.count[chrom][i]
                ratio = ratios[i]
                covd[latecov].append(ratio)
        
        return covd

    def get_params_from_local_medRatioNorm_dict_with_globalweighting(self, covd, globalweight=30):
        #covd = output from create_local_medRatioNorm_dict
        #globalweight = how many global_med values to add to a local bin
        global_med = [np.median(np.concatenate(covd.values()))]*globalweight
        covmeds = dict()
        for key in sorted(covd.keys()):
            ratios = np.concatenate([global_med,covd[key]])
            covmeds[key] = np.median(ratios)
        return covmeds

    def get_params_from_local_medRatioNorm_dict_with_glocalweighting(self, covd, globalweight=10, minlocalbins=3, minpropdata=0.1):
        #covd = output from create_local_medRatioNorm_dict
        #globalweight = how many global_med values to add to a local bin
        N = len(np.concatenate(covd.values()))
        pN = round(minpropdata*N)
        #print 'pN', pN
        global_med = [np.median(np.concatenate(covd.values()))]*globalweight
        #print 'G', global_med
        covmeds = dict()
        latecovs = sorted(covd.keys())
        halfbins = int((minlocalbins+1)//2) ## halves from odds rounded up
        n = len(latecovs)
        for i in range(n):
            ## FIRST SATISFY MIN BINS
            if minlocalbins >=3:
                l = max(0, i-halfbins)
                r = min(i+halfbins, n-1)
                nl = i-l
                nr = r - i
                if l == 0:
                    add = minlocalbins - nl
                    r += add
                    if r > n-1: ## CATCH -- this shouldn't happen
                        print "WARN: get_params_from_local_medRatioNorm_dict_with_glocalweighting #1"
                        r = n-1
                if r == n-1:
                    sub = minlocalbins - nr
                    l -= sub
                    if l < 0: ## CATCH -- this shouldn't happen
                        print "WARN: get_params_from_local_medRatioNorm_dict_with_glocalweighting #2"
                        l = 0
                bins =  [latecovs[j] for j in range(l,r)] 
            else:
                bins = [latecovs[i]] # [covd[j] for j in range(i,i+minlocalbins)]
            nbins = len(bins)
            #print "nbins post minbins", nbins

            ## SECOND, SATISFY MIN PROP OF DATA
            direction = 'R'
            localvals = list(np.concatenate([covd[latecov] for latecov in bins]))
            nvals = len(localvals)
            finiteloop=n
            while nvals < pN:
                if direction == 'R' and r <= n-1:
                    # use r before r change since it only went up to but did not include r before
                    localvals.append(covd[latecovs[r]])
                    r += 1
                    direction = 'L'
                    nvals = len(localvals)
                    finiteloop=n
                elif direction == 'L' and l > 0:
                    # use l AFTER l change since used l already
                    l-=1
                    localvals.append(covd[latecovs[l]])
                    direction = 'R'
                    nvals = len(localvals)
                    finiteloop=n
                else:
                    finiteloop-=1
                    if finiteloop <= 0:
                        print "WARN: get_params_from_local_medRatioNorm_dict_with_glocalweighting #3"
                        break


            ## GET GLOCAL MEDIAN RATIO
            ratios = np.concatenate([global_med, localvals])
            #print ratios
            covmeds[latecovs[i]] = np.median(ratios)
        return covmeds







class MultiCovBed(object):
    ## assumes all covbeds in list have same exact elements in 1st 3 colums and are all identically ordered
    ## this should be true if all are outputs from getcov using the same genome file (just different bams)
    def __init__(self, stagefile):
        # 'stagefile' is a FOFN-like file with 2 columns: 1=stage_integer,2=bincounts_filepaths
        ##  where stage integer is time points -- e.g. 1,2,3,4,5
        ##  stage_ints can be used on replicates (i.e. can use same int >1 time)
        self.stagefile = stagefile
        self.nfiles = 0
        self.nbins = 0
        self.covbeds = {}
        self.stages = {}
        self.stagelist = []
        self.files = {}
        self.count = {}
        self.median = {}
        self.start = {}
        self.end = {}
        self.chromosomes = set([])
        self._parse_stagefile()
        self._extract_data()
        self.corscores = {}
        self.rksmooth = None
        self.smooth_corscores = {}
        self.cor_states = {}
        self.s50 = {}
        self.a50 = {}
        self.filtered = {'start':{}, 'end':{}, 'count':{k:{} for k in range(self.nfiles)}}
        self.filteredchromosomes = []
        self.ntest = None


    ######################
    ##  INITIALIZATION
    ######################
    def _parse_stagefile(self):
        i = 0
        with open(self.stagefile) as f:
            for line in f:
                stage, fname = line.strip().split()
                self.stages[i] = int(stage)
                self.files[i] = fname
                self.stagelist.append(int(stage))
                i+=1
        self.nfiles = len(self.files.keys())
            
        
    def _extract_data(self):
        for i in sorted(self.files.keys()):
            self.add_covbed(findex = i)

    def add_covbed(self, findex):
        covbed = CovBed(self.files[findex])
        if not self.chromosomes:
            self._initialize(covbed)
        self.count[findex] = covbed.get_count_dict()
            

    def _initialize(self, covbed):
        self.start = covbed.get_start_dict()
        self.end = covbed.get_end_dict()
        self.chromosomes = sorted(list(covbed.get_chromosomes()))


    ######################
    ## Operations
    ######################
    def _get_median(self, findex):
        counts = np.concatenate(self.count[findex].values())
        self.median[findex] = float(np.median(counts))
        
    def get_median(self, findex, refresh=False):
        if refresh:
            self._get_median(findex)
        try:
            return self.median[findex]
        except KeyError as e:
            self._get_median(findex)
            return self.median[findex]

    def _refresh_medians(self):
        ## will re-calculate medians every time called
        for findex in range(self.nfiles):
            self._get_median(findex)

    def normalize(self, x, denom):
        #x is np.array
        # denom is float
        return x/denom

    def normalize_findex_by_x(self,findex, x):
        for chrom in self.chromosomes:
            self.count[findex][chrom] = self.normalize(self.count[findex][chrom], x)

    def nomalize_data_by_xdict(self,xdict):
        for findex in range(self.nfiles):
            self.normalize_findex_by_x(findex, xdict[findex])

    def median_normalize_findex(self,findex, refresh=False):
        self.normalize_findex_by_x(findex, x = self.get_median(findex, refresh))

    def median_normalize_data(self, refresh=False):
        for findex in range(self.nfiles):
            self.median_normalize_findex(findex, refresh)

    # def chromosome_median_normalize_data(self):
    #     ## DEVELOPING ON JULY 17, 2025; STILL NEEDS TESTING.
    #     '''Normalize the counts for a given chromosome by the median of that chromosome.'''
    #     for findex in range(self.nfiles):
    #         for chrom in self.chromosomes:
    #             ## Get median count for this chromosome
    #             median_count = np.median(self.count[findex][chrom])
    #             if median_count > 0:
    #                 ## Normalize the counts for this chromosome by the median count
    #                 self.count[findex][chrom] = self.normalize(self.count[findex][chrom], median_count)
    #             else:
    #                 ## If median count is 0, Try the mean.
    #                 mean_count = np.mean(self.count[findex][chrom])
    #                 if mean_count > 0:
    #                     self.count[findex][chrom] = self.normalize(self.count[findex][chrom], mean_count)
    #                 else:
    #                     ## CAREFUL HERE: It is silently setting all counts to 1 to avoid division by zero in downstream analyses; but this may not be the desired behavior (e.g. if input is raw, raw ratio, or log2 ratio would have diff interpretations).
    #                     ## If both median and mean are 0, set all counts to 1; then they will be normalized to 1.
    #                     self.count[findex][chrom] = np.ones_like(self.count[findex][chrom])

    def ksmooth_counts(self, bw=10000):
        for chrom in self.chromosomes:
            for findex in range(self.nfiles):
                x = self.start[chrom]
                y = self.count[findex][chrom]
                k = ksmooth(x = fltvec(x), y = fltvec(y), bandwidth = bw)
                self.count[findex][chrom] = np.array(k[1])

    def expanded_bdg(self, bdg):
        string = ''
        for chrom in self.chromosomes:
            for i in range(len( self.start[chrom] )):
                string += ('\t').join([chrom, str(self.start[chrom][i]), str(self.end[chrom][i]), str(bdg[chrom][i])]) + "\n"
        return string

    def expanded_bdg_two_cols(self, bdg1, bdg2):
        string = ''
        for chrom in self.chromosomes:
            for i in range(len( self.start[chrom] )):
                string += ('\t').join([chrom, str(self.start[chrom][i]), str(self.end[chrom][i]), str(bdg1[chrom][i]), str(bdg2[chrom][i])]) + "\n"
        return string
    
    def collapsed_bdg(self, bdg):
        string = ''
        for chrom in self.chromosomes:
            #init
            start = self.start[chrom][0]
            value = bdg[chrom][0]
            for i in range(1, len(self.start[chrom]) ):
                if bdg[chrom][i] != value:
                    string += ('\t').join([chrom, str(start), str(self.end[chrom][i-1]), str(value)]) + "\n"
                    start = self.start[chrom][i]
                    value = bdg[chrom][i]
            ##finish chrom
            string += ('\t').join([chrom, str(start), str(self.end[chrom][i]), str(value)]) + "\n"
        return string

    def get_bdg(self, bdg, collapsed=False):
        if not collapsed:
            return self.expanded_bdg(bdg)
        else:
            return self.collapsed_bdg(bdg)
        
    def find_slopes(self, stagelist=''):
        if not stagelist:
            stagelist = self.stagelist
        self.slopes = {}
        for chrom in self.chromosomes:
            self.slopes[chrom] = []
            for i in range(len(self.start[chrom])):
                counts = [self.count[j][chrom][i] for j in range(self.nfiles)]
                try:
                    slope = np.polyfit(x = stagelist, y = counts, deg = 1)[0]
                except FloatingPointError:
                    slope = 0
                self.slopes[chrom].append(slope)

    def get_slope_bdg(self, collapsed=False):
        ##assumes self.slopes already present
        return self.get_bdg(self.slopes, collapsed)


                
    def cor_score(self, stagelist=''):
        if not stagelist:
            stagelist = self.stagelist
        for chrom in self.chromosomes:
            self.corscores[chrom] = []
            for i in range(len(self.start[chrom])):
                counts = [self.count[j][chrom][i] for j in range(self.nfiles)]
                try:
##                    score = np.nan_to_num( np.corrcoef(x = [stagelist, counts])[0,1] )
                    score = np.corrcoef(x = [stagelist, counts])[0,1] 
                except FloatingPointError:
                    score = 0
                self.corscores[chrom].append(score)


    def get_corscore_bdg(self, collapsed=False):
        if not self.corscores:
            self.cor_score()
        return self.get_bdg(self.corscores, collapsed)

    def ksmooth_corscores(self, bw=10000):
        if not self.corscores:
            self.cor_score()
        for chrom in self.chromosomes:
            x = self.start[chrom]
            y = self.corscores[chrom]
            k = ksmooth(x = fltvec(x), y = fltvec(y), bandwidth = bw)
            self.smooth_corscores[chrom] = np.array(k[1])


    def get_smooth_corscore_bdg(self, collapsed=False):
        return self.get_bdg(self.smooth_corscores, collapsed)


    def get_cor_states(self, smoothed=False, emodel="normal"):
        if not self.corscores:
            self.cor_score()
        if smoothed and not self.smooth_corscores:
            sys.stderr.write("Smoothing cor scores with default bandwidth")
            self.ksmooth_corscores()
        if smoothed:
            scores = self.smooth_corscores
        else:
            scores = self.corscores
        for chrom in self.chromosomes:
##            sys.stderr.write( chrom + "\n" )
            v = puffR.viterbi_puff(emissions = puffR.emissions, transitions = puffR.transitions, initial = puffR.initial, states = intvec([1,2,3]), emitted_data = fltvec(scores[chrom]), emodel = emodel, logprobs=False)
##            f = puffR.forward_puff(emissions = puffR.emissions, transitions = puffR.transitions, initial = puffR.initial, states = intvec([-1,0,1]), emitted_data = fltvec(scores[chrom]), emodel = emodel, logprobs=False)
##            b = puffR.backward_puff(emissions = puffR.emissions, transitions = puffR.transitions, initial = puffR.initial, states = intvec([-1,0,1]), emitted_data = fltvec(scores[chrom]), emodel = emodel, logprobs=False)
            
            self.cor_states[chrom] = np.array(list(v[0]))


    def get_cor_state_bdg(self, collapsed=False):
        return self.get_bdg(self.cor_states, collapsed)

    

    
    def calc_s50(self, x=[50], stagelist=''):
        if not stagelist:
            stagelist = self.stagelist
        for chrom in self.chromosomes:
            self.s50[chrom] = []
            for i in range(len(self.start[chrom])):
                counts = [self.count[j][chrom][i] for j in range(self.nfiles)]
                ans = s50(counts,stagelist,x=x)
                ans = ans[x[0]]
                self.s50[chrom].append(ans)

                
    def get_s50_bdg(self):
        if not self.s50:
            self.calc_s50()
        string = ''
        for chrom in self.chromosomes:
            for i in range(len( self.start[chrom] )):
                if self.s50[chrom][i] != ".":
                    string += ('\t').join([chrom, str(self.start[chrom][i]), str(self.end[chrom][i]), str(self.s50[chrom][i])]) + "\n"
        return string

    def analyze_state0_bins(self):
        ## assumes median normalized (or otherwise)
        if not self.corscores:
            self.cor_score()
        if not self.cor_states:
            self.get_cor_states()
        counts = {k:[] for k in range(self.nfiles)}
        for chrom in self.chromosomes:
            for i in range(len( self.start[chrom] )):
                if self.cor_states[chrom][i] == 2:#states are 1,2,3 for -1,0,1
                    for j in range(self.nfiles):
                        counts[j].append( self.count[j][chrom][i] )
        ## mean would be skewed above 1 since even if the count of numbers < 1 equals the count of numbers > 1, the mean will be >1.
        total_median = [np.median(counts.values())]
        medians = {findex:np.median(counts[findex]) for findex in range(self.nfiles)}
        self.state0_medians = {0:medians,1:total_median}

    def _need_more_cycles(self):
        for m in self.state0_medians[0].values():
            if m != 1:
                return True
        return False

    def _normalize_data_by_state0_median(self):
        for findex in range(self.nfiles):
            if self.state0_medians[0][findex] != 1:
                self.normalize_findex_by_x(findex, self.state0_medians[0][findex])

    def find_cn1(self, n_iter=10, max_empty_bin_pct=0.4, max_offending_samples_pct=0.4, verbose=True):
        if verbose:
            sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 - getting medians before and after filtering..\n")
        self._refresh_medians() ## get medians for first time
        if verbose:
            sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 - medians before filtering were %s...\n" % (str(self.median)))
            sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 - filtering..\n")
        self.filter_null_contigs(max_empty_bin_pct, max_offending_samples_pct) ## filter bad contigs
        self._refresh_medians() ## get medians after filtration
        if verbose:
            sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 - medians after filtering were %s...\n" % (str(self.median)))
            sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 - median normalizing..\n")
        self.median_normalize_data() # median norm
        self._refresh_medians() ## get medians after normalization
        if verbose:
            sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 - medians (of all bins) after normalizing should be 1 and were %s...\n" % (str(self.median)))
            sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 - getting bin correlations \n")
        self.cor_score() # give correlation score to each bin
        if verbose:
            sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 - getting state path \n")
        self.get_cor_states() # get state path through bins
        self.analyze_state0_bins() # get median of 0-correlation (most likely cn=1) bins
        for i in range(n_iter):
            if verbose:
                sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 iter %d starting 0-corr medians: %s...\n" % (i+1, str(self.state0_medians)))
            if self._need_more_cycles(): # if 0-corr median from any sample != 1, re-normalize by 0-corr median, est new correlations, get updated statepath, get new 0corr medians
                self._normalize_data_by_state0_median()
                self.cor_score()
                self.get_cor_states()
                self.analyze_state0_bins()
            else:
                if verbose:
                    sys.stderr.write(str(datetime.datetime.now()) +": ..CN=1 iter %d canceled - 0-corr medians all = 1...\n" % (i+1))
                break
        
    def filter_null_contigs(self, max_empty_bin_pct=0.4, max_offending_samples_pct=0.4):
        ##default: if >40% of samples have >40% empty bins on a chrom, remove it
        sample_thres = max_offending_samples_pct*self.nfiles
        for chrom in self.chromosomes:
            nbins = len(self.start[chrom])
            bin_thres = max_empty_bin_pct*nbins
            n_bad_samp = 0
            for findex in range(self.nfiles):
                n_bad_bins = len(self.count[findex][chrom][self.count[findex][chrom] == 0])
                if n_bad_bins > bin_thres:
                    n_bad_samp += 1
            if n_bad_samp > sample_thres:
                self.filtered['start'][chrom] = self.start.pop(chrom)
                self.filtered['end'][chrom] = self.end.pop(chrom)
                for findex in range(self.nfiles):
                    counts = (self.count[findex]).pop(chrom)
                    self.filtered['count'][findex][chrom] = counts
                self.filteredchromosomes.append(chrom)
                self.chromosomes.remove(chrom)


            
    def pct_state(self, state=3):
        total = 0
        nstate = 0
        for chrom in self.chromosomes:
            total += len(self.cor_states[chrom])
            nstate += sum(self.cor_states[chrom][self.cor_states[chrom] == state])
        return 100.0*nstate/total

    def n_state(self, state=3):
        nstate = 0
        for chrom in self.chromosomes:
            nstate += sum(self.cor_states[chrom][self.cor_states[chrom] == state])
        return nstate
    
    def eFDR1(self, stage=3, n_iter=10):
        ## shuffles stages -- not best way since many permutations uphold correlations...
        ## assumes find_cn1 already run
        if self.ntest is None:
            self.ntest = self.n_state(3)
        controls = []
        for i in range(n_iter):
            stagelist = self.stagelist[:] #clone it
            np.random.shuffle(stagelist)
            sys.stderr.write(str(stagelist)+"\n") ## PRINT TEMP
            self.cor_score(stagelist)
            self.get_cor_states()
            controls.append(self.n_state(3))
        return self.ntest, controls, 100.0*np.array(controls)/self.ntest

    def eFDR2(self):
        ## 
        pass

    def pval(self):
        # for each bin, shuffle scores, take cor, store -- get 1000 of these, p ~ n_gt_cor/N
        # BH correct p-values - only keep bins with q < 0.1
        ## OR -- take state+ bins, combine counts in all for each stage, do this
        ## that would be less tests and it would not treat all bins independently - would treat regions independently
        ## hmmm.... maybe not worth doing this at all...
        pass

    def discretize_cor_values(self):
        self.dcorscores = {}
        for chrom in self.chromosomes:
            self.dcorscores[chrom] = map(round, np.array(self.corscores[chrom])*4+5) #9-sided dice emission symbols 1-9 (for Rindexing) where 1-4 neg, 5=0, 6-9pos
    ######################
    ## Printing etc
    ######################
    def get_counts_bdg(self):
        string = ''
        for chrom in self.chromosomes:
            for i in range(len( self.start[chrom] )):
                string += ('\t').join([chrom, str(self.start[chrom][i]), str(self.end[chrom][i])] + [str(self.count[j][chrom][i]) for j in range(self.nfiles)]) + "\n"
        return string

    def get_filtered_contigs_bdg(self):
        string = ''
        for chrom in self.filteredchromosomes:
            for i in range(len( self.filtered['start'][chrom] )):
                string += ('\t').join([chrom, str(self.filtered['start'][chrom][i]), str(self.filtered['end'][chrom][i])] + [str(self.filtered['count'][j][chrom][i]) for j in range(self.nfiles)]) + "\n"
        return string
    
    def __str__(self):
        return self.get_counts_bdg()

    
            









## DOES NOT BEHAVE IN USEFUL WAY - i.e. results not more useful than s50 (probably less so).
##    def calc_a50(self,x=[50]):
##        ## ASSUMES COUNTS ARE MEDIAN NORMALIZED
##        ## for all cn=1 areas, s50=3 -- i.e. 50% of normalized read counts seen halfway through
##        ## a50 is trying to mask cn=1 to highlight when 50% of amplification is done
##        ## it does this by subtracting 1 from the median normalized counts (and taking max of that or 0 to avoid negatives)
##        ##  --> and later ignoring anything that is "."
##        ## this is experimental - not totally sure it will give what I want
##        for chrom in self.chromosomes:
##            self.a50[chrom] = []
##            for i in range(len(self.start[chrom])):
##                counts = [max([self.count[j][chrom][i]-1,0]) for j in range(self.nfiles)]
##                ans = s50(counts,self.stagelist,x=x)
##                ans = ans[x[0]]
##                self.a50[chrom].append(ans)
##                
##    def get_a50_bdg(self):
##        if not self.a50:
##            self.calc_a50()
##        string = ''
##        for chrom in self.chromosomes:
##            for i in range(len( self.start[chrom] )):
##                if self.a50[chrom][i] != ".":
##                    string += ('\t').join([chrom, str(self.start[chrom][i]), str(self.end[chrom][i]), str(self.a50[chrom][i])]) + "\n"
####                if i > 20: break
####            if i > 20: break
##        return string
