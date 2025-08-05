import sys, datetime
from CovBedClass import *
from pk2txt import bdgmsg, newmsg


class NormalizeProtocol(object):
    def __init__(self, args):
        #latestage, protocol=1, earlystage=False, pseudo=0.1, bandwidth=2500, quiet=False, impute=False, replace=False, replace_with='0', replace_this='.'
        self.args = args
        self._normalize()
        self._finalize()
        
    def _set_protocol(self):
        if self.args.protocol1:
            self.protocol = 1
            self._run_protocol = self._protocol1
        elif self.args.protocol2:
            self.protocol = 2
            self._run_protocol = self._protocol2
        elif self.args.protocol3:
            self.protocol = 3
            self._run_protocol = self._protocol3
        elif self.args.protocol4:
            self.protocol = 4
            self._run_protocol = self._protocol4
        elif self.args.protocol5:
            self.protocol = 5
            self._run_protocol = self._protocol5
        elif self.args.protocol6:
            self.protocol = 6
            self._run_protocol = self._protocol6
        elif self.args.protocol7:
            self.protocol = 7
            self._run_protocol = self._protocol7
        elif self.args.protocol8:
            self.protocol = 8
            self._run_protocol = self._protocol8
        elif self.args.protocol9:
            self.protocol = 9
            self._run_protocol = self._protocol9
        elif self.args.protocol10:
            self.protocol = 10
            self._run_protocol = self._protocol10
        elif self.args.protocol11:
            self.protocol = 11
            self._run_protocol = self._protocol11
        elif self.args.protocol12:
            self.protocol = 12
            self._run_protocol = self._protocol12
        elif self.args.protocol13:
            self.protocol = 13
            self._run_protocol = self._protocol13
        elif self.args.protocol14:
            self.protocol = 14
            self._run_protocol = self._protocol14
        elif self.args.protocol15:
            self.protocol = 15
            self._run_protocol = self._protocol15
        elif self.args.protocol16:
            self.protocol = 16
            self._run_protocol = self._protocol16
        elif self.args.protocol17:
            self.protocol = 17
            self._run_protocol = self._protocol17
        elif self.args.protocol18:
            self.protocol = 18
            self._run_protocol = self._protocol18
        elif self.args.protocol19:
            self.protocol = 19
            self._run_protocol = self._protocol19
        elif self.args.protocol20:
            self.protocol = 20
            self._run_protocol = self._protocol20
        elif self.args.protocol21:
            self.protocol = 21
            self._run_protocol = self._protocol21
        elif self.args.protocol22:
            self.protocol = 22
            self._run_protocol = self._protocol22
        elif self.args.protocol23:
            self.protocol = 23
            self._run_protocol = self._protocol23
        elif self.args.protocol24:
            self.protocol = 24
            self._run_protocol = self._protocol24
        elif self.args.protocol25:
            self.protocol = 25
            self._run_protocol = self._protocol25
        elif self.args.protocol26:
            self.protocol = 26
            self._run_protocol = self._protocol26
        elif self.args.protocol27:
            self.protocol = 27
            self._run_protocol = self._protocol27
        elif self.args.protocol28:
            self.protocol = 28
            self._run_protocol = self._protocol28
        elif self.args.protocol29:
            self.protocol = 29
            self._run_protocol = self._protocol29
        elif self.args.protocol30:
            self.protocol = 30
            self._run_protocol = self._protocol30
        elif self.args.protocol31:
            self.protocol = 31
            self._run_protocol = self._protocol31
        elif self.args.protocol32:
            self.protocol = 32
            self._run_protocol = self._protocol32
        elif self.args.protocol33:
            self.protocol = 33
            self._run_protocol = self._protocol33
        elif self.args.protocol34:
            self.protocol = 34
            self._run_protocol = self._protocol34

            
    def _protocol1(self):
        self.late.median_normalize_data()
        if self.early:
            self.early.median_normalize_data()
            self.late.normalize_to_other(self.early,
                                              self.args.pseudo)

    def _protocol2(self):
        #print self.late.get_median()
        self.late.ksmooth_counts(bw=self.args.bandwidth,
                                      rescueNaN=self.args.replaceNaN,
                                      localWindow=self.args.localwinsize)
        self.late.median_normalize_data()
        if self.early:
            self.early.ksmooth_counts(bw=self.args.bandwidth,
                                      rescueNaN=self.args.replaceNaN,
                                      localWindow=self.args.localwinsize)
            self.early.median_normalize_data()
            self.late.normalize_to_other(self.early,
                                              self.args.pseudo)

    def _protocol3(self):
        self.late.median_normalize_data()
        self.late.ksmooth_counts(bw=self.args.bandwidth,
                                      rescueNaN=self.args.replaceNaN,
                                      localWindow=self.args.localwinsize)
        if self.early:
            self.early.median_normalize_data()
            self.early.ksmooth_counts(bw=self.args.bandwidth,
                                      rescueNaN=self.args.replaceNaN,
                                      localWindow=self.args.localwinsize)
            self.late.normalize_to_other(self.early,
                                         self.args.pseudo)

    def _protocol4(self):
        self.late.median_normalize_data()
        if self.early:
            self.early.median_normalize_data()
            self.late.normalize_to_other(self.early,
                                         self.args.pseudo)
        self.late.ksmooth_counts(bw=self.args.bandwidth,
                                      rescueNaN=self.args.replaceNaN,
                                      localWindow=self.args.localwinsize)

    def _protocol5(self):
        #smoothing only --  late/early before smooth
        if self.early:
            self.late.normalize_to_other(self.early,
                                         self.args.pseudo)
        self.late.ksmooth_counts(bw=self.args.bandwidth,
                                      rescueNaN=self.args.replaceNaN,
                                      localWindow=self.args.localwinsize)

    def _protocol6(self):
        #smoothing only --  late/early AFTER smooth
        self.late.ksmooth_counts(bw=self.args.bandwidth,
                                      rescueNaN=self.args.replaceNaN,
                                      localWindow=self.args.localwinsize)
        if self.early:
            self.early.ksmooth_counts(bw=self.args.bandwidth,
                                      rescueNaN=self.args.replaceNaN,
                                      localWindow=self.args.localwinsize)
            self.late.normalize_to_other(self.early,
                                         self.args.pseudo)

    def _protocol7(self):
        #no smoothing, no median norm
        # late:early only if early present
        if self.early:
            self.late.normalize_to_other(self.early,
                                         self.args.pseudo)

    def _protocol8(self):
        self.late.computeSkew()

    def _protocol9(self):
        self.late.computePercentChange()

    def _protocol10(self):
        self.late.computeSkewChange()

    def _protocol11(self):
        self.late.computePercentChange()

    def _protocol12(self):
        # median ratio norm - e.g. similar to DEseq (or TMM in EdgeR) -- global version
        # if no "early" sample present, then this just returns median norm late -- like protocol 1
        # step 1: get ratios.
        # step 2: divide all ratios by median ratio.
        # This is analogous to subtracting the median log2 ratio from each log2 ratio.
        if self.early:
            self.late.normalize_to_other(self.early,
                                         self.args.pseudo)
        self.late.median_normalize_data()

    def _protocol13(self):
        # median ratio norm w/ pre-ratio smoothing - e.g. similar to DEseq (or TMM in EdgeR) -- global version
        self.late.ksmooth_counts(bw=self.args.bandwidth,
                                      rescueNaN=self.args.replaceNaN,
                                      localWindow=self.args.localwinsize)
        if self.early:
            self.early.ksmooth_counts(bw=self.args.bandwidth,
                                      rescueNaN=self.args.replaceNaN,
                                      localWindow=self.args.localwinsize)
            self.late.normalize_to_other(other=self.early,
                                         pseudocount=self.args.pseudo)
        self.late.median_normalize_data()

    def _protocol14(self):
        # median ratio norm w/ post-ratio smoothing - e.g. similar to DEseq (or TMM in EdgeR) -- global version
        if self.early:
            self.late.normalize_to_other(self.early,
                                         self.args.pseudo)
        self.late.ksmooth_counts(bw=self.args.bandwidth,
                                      rescueNaN=self.args.replaceNaN,
                                      localWindow=self.args.localwinsize)
        self.late.median_normalize_data()
        

    def _protocol15(self):
        # median ratio norm w/ end-smoothing (after both ratio and median norm) - e.g. similar to DEseq (or TMM in EdgeR) -- global version
        if self.early:
            self.late.normalize_to_other(self.early,
                                         self.args.pseudo)
        self.late.median_normalize_data()
        self.late.ksmooth_counts(bw=self.args.bandwidth,
                                      rescueNaN=self.args.replaceNaN,
                                      localWindow=self.args.localwinsize)

    def _protocol16(self):
        # experimental glocal median ratio norm w/ end-smoothing - 
        self.late.normalize_with_glocalMedRatioNorm(other=self.early,
                                                    pseudocount=self.args.pseudo,
                                                    globalweight=100,
                                                    minlocalbins=5,
                                                    minpropdata=0.333) ## LAST 3 HARDCODED FOR NOW -- less local
    def _protocol17(self):
        ## First median normalizes both samples
        ## Then scales to a target median bin cov -- default 1000.
        ## Then adds pseudocount - default 1. (so default pseudo is 0.1% median, 1000-fold less than median)
        ## Then normalizes late to early
        ## Then median ratio normalizes a la DEseq/EdgeR
        self.late.median_normalize_data()
        self.late.scale_data(scale=self.args.scalecov)
        if self.early:
            self.early.median_normalize_data()
            self.early.scale_data(scale=self.args.scalecov)
            
            self.late.normalize_to_other(self.early,
                                              self.args.pseudo)
        self.late.median_normalize_data(relearn=True)


    def _protocol18(self):
        ## Robust Z scores
        if self.early:
            self.late.normalize_to_other(self.early,
                                              self.args.pseudo)
        self.late.robust_z_normalize_data()

    def _protocol19(self):
        ## Rank scores
        if self.early:
            self.late.normalize_to_other(self.early,
                                              self.args.pseudo)
        self.late.rank_normalize_data()

    def _protocol20(self):
        ## Robust Z score difference
        if not self.early:
            sys.stderr.write("This requires late (test) and early (control) files, not just late (test).....\n")
            quit()
        self.late.robust_z_normalize_data()
        self.early.robust_z_normalize_data()
        self.late.subtract_other(self.early)

    def _protocol21(self):
        ## Rank difference
        if not self.early:
            sys.stderr.write("This requires late (test) and early (control) files, not just late (test).....\n")
            quit()
        self.late.rank_normalize_data()
        self.early.rank_normalize_data()
        self.late.subtract_other(self.early)

    def _protocol22(self):
        ## SPMR - all counts are summed up, all bins divided by sum and multipled by 1e6 (default)
        self.late.spxr_normalize_data(x=self.args.SPXR)
        if self.early:
            self.early.spxr_normalize_data(x=self.args.SPXR)
            self.late.normalize_to_other(self.early,
                                              self.args.pseudo)
    def _protocol23(self):
        ## Std Rank: Rank scores, subtract theoretical middle rank (min+max)/2, fold-normalize by middle rank
        if self.early:
            self.late.normalize_to_other(self.early,
                                              self.args.pseudo)
        self.late.rank_standardize_data()

    def _protocol24(self):
        ## pct difference: (Z_t-Z_c)/abs(Z_c)
        ## t < c; both neg
        ## (-1 - -0.5)/(-0.5) = -0.5/-0.5 = 1
        ## (-1 - -0.5)/abs(-0.5) = -0.5/0.5 = -1
        ## t > c; both neg
        ## (-0.5 - -1) /(-1) = 0.5/-1 = -1.5
        ## (-0.5 - -1) /abs(-1) = 0.5/1 = 0.5
        if not self.early:
            sys.stderr.write("This requires late (test) and early (control) files, not just late (test).....\n")
            quit()
        self.late.pct_diff_from_other(self.early, setToControlDist=self.args.setToControlDist, pseudoZeroBins=self.args.pseudoZeroBins, addMinOtherPlusOneToBoth=self.args.addMinOtherPlusOneToBoth)

    def _protocol25(self):
        ## pct skew: (Z_t - Z_c)/(abs(Z_t) + abs(Z_c))
        ## t < c; both neg
        ## (-1 - -0.5)/(-1 + -0.5) = -0.5/-1.5 = 0.333
        ## (-1 - -0.5)/abs(-1 + -0.5) = -0.5/1.5 = -0.333
        ## (-1 - -0.5)/(abs(-1) + abs(-0.5)) = -0.5/1.5 = -0.333
        ## t > c; both neg
        ## (-0.5 - -1) /(-0.5 + -1) = 0.5/-1.5 = -0.333
        ## (-0.5 - -1) /abs(-0.5 + -1) = 0.5/1.5 = 0.333
        ## (-0.5 - -1) /(abs(-0.5) + abs(-1)) = 0.5/1.5 = 0.333
        ##
        ## t < c; pos and neg
        ## (-1 - 0.5)/(-1 + 0.5) = -1.5/-0.5 = 3
        ## (-1 - 0.5)/abs(-1 + 0.5) = -1.5/0.5 = -3
        ## (-1 - 0.5)/(abs(-1) + abs(0.5)) = -1.5/1.5 = -1
        ## t > c; pos and neg
        ## (0.5 - -1) /(0.5 + -1) = 1.5/-0.5 = -0.333
        ## (0.5 - -1) /abs(0.5 + -1) = 1.5/0.5 = 0.333
        ## (0.5 - -1) /(abs(0.5) + abs(-1)) = 1.5/1.5 = 1
        ##
        ## t < c; both pos
        ## (0.5 - 1)/(0.5 + 1) = -0.5/1.5 = -0.333
        ## (0.5 - 1)/abs(0.5 + 1) = -0.5/1.5 = -0.333
        ## (0.5 - 1)/(abs(0.5) + abs(1)) = -0.5/1.5 = -0.333
        ## t > c; both pos
        ## (1 - 0.5) /(1 + 0.5) = 0.5/1.5 = 0.333
        ## (1 - 0.5) /abs(1 + 0.5) = 0.5/1.5 = 0.333
        ## (1 - 0.5) /(abs(1) + abs(0.5)) = 0.5/1.5 = 0.333
        ##
        ## inputs should really just be positive.
        ## However, the "most correct" way when anticipating pos and neg is to sum the abs vals of each.
        
        if not self.early:
            sys.stderr.write("This requires late (test) and early (control) files, not just late (test).....\n")
            quit()
        self.late.pct_skew_given_other(self.early, setToControlDist=self.args.setToControlDist, pseudoZeroBins=self.args.pseudoZeroBins, addMinOtherPlusOneToBoth=self.args.addMinOtherPlusOneToBoth)





    def _protocol26(self):
        ## Robust Z score pct difference: (Z_t-Z_c)/abs(Z_c)
        if not self.early:
            sys.stderr.write("This requires late (test) and early (control) files, not just late (test).....\n")
            quit()
        self.late.robust_z_normalize_data()
        self.early.robust_z_normalize_data()
        self.late.pct_diff_from_other(self.early, setToControlDist=self.args.setToControlDist, pseudoZeroBins=self.args.pseudoZeroBins, addMinOtherPlusOneToBoth=self.args.addMinOtherPlusOneToBoth)

    def _protocol27(self):
        ## Robust Z score pct skew: (Z_t - Z_c)/(abs(Z_t) + abs(Z_c))
        if not self.early:
            sys.stderr.write("This requires late (test) and early (control) files, not just late (test).....\n")
            quit()
        self.late.robust_z_normalize_data()
        self.early.robust_z_normalize_data()
        self.late.pct_skew_given_other(self.early, setToControlDist=self.args.setToControlDist, pseudoZeroBins=self.args.pseudoZeroBins, addMinOtherPlusOneToBoth=self.args.addMinOtherPlusOneToBoth)



    def _protocol28(self):
        ## Rank pct difference: (Z_t-Z_c)/abs(Z_c)
        if not self.early:
            sys.stderr.write("This requires late (test) and early (control) files, not just late (test).....\n")
            quit()
        self.late.rank_normalize_data()
        self.early.rank_normalize_data()
        self.late.pct_diff_from_other(self.early)

    def _protocol29(self):
        ## Rank pct skew: (Z_t - Z_c)/(abs(Z_t) + abs(Z_c))
        if not self.early:
            sys.stderr.write("This requires late (test) and early (control) files, not just late (test).....\n")
            quit()
        self.late.rank_normalize_data()
        self.early.rank_normalize_data()
        self.late.pct_skew_given_other(self.early)

    def _protocol30(self):
        # median ratio norm for local windows - e.g. similar to DEseq (or TMM in EdgeR) -- PHYSICALLY-LOCAL version
        # if no "early" sample present, then this just returns phys-local median norm late 
        if self.early:
            self.late.normalize_to_other(self.early,
                                         self.args.pseudo)
        self.late.local_median_normalize_data(halfwidth=self.args.halfwidth)

    def _protocol31(self):
        # chromosome-specific median ratio norm
        # if no "early" sample present, then this just returns chromosome-specific median norm late.
        # step 1: get ratios, if both late and early present.
        # step 2: divide all ratios on a given chromosome by the median ratio of that chromosome.
        # This is analogous to subtracting the chromosome-specific median log2 ratio from each log2 ratio on that chromosome.
        if self.early:
            self.late.normalize_to_other(self.early,
                                         self.args.pseudo)
        self.late.chromosome_median_normalize_data()

    def _protocol32(self):
        # Median smoothing in windows defined by halfwidth.
        if self.early:
            self.late.normalize_to_other(self.early,
                                         self.args.pseudo)
        self.late.local_median_smooth_data(halfwidth=self.args.halfwidth)

    def _protocol33(self):
        # Trimmed Mean smoothing in windows defined by halfwidth.
        if self.early:
            self.late.normalize_to_other(self.early,
                                         self.args.pseudo)
        self.late.local_trimmed_mean_smooth_data(halfwidth=self.args.halfwidth)

    def _protocol34(self):
        # Mean smoothing in windows defined by halfwidth.
        if self.early:
            self.late.normalize_to_other(self.early,
                                         self.args.pseudo)
        self.late.local_mean_smooth_data(halfwidth=self.args.halfwidth)

    def _normalize(self):
        if not self.args.quiet:
            newmsg("loading late stage file")
        self.late = CovBed(self.args.latestage,
                      replace=self.args.replace,
                      replace_with=self.args.replace_with,
                      replace_this=self.args.replace_this,
                      stringcols=self.args.stringcols)
        if self.args.impute:
            if not self.args.quiet:
                newmsg("imputing late stage bins with missing data")
            self.late.impute_zeros(bw=self.args.impute)
                
        if self.args.earlystage:
            if not self.args.quiet:
                newmsg("loading early stage file")
            self.early = CovBed(self.args.earlystage,
                                replace=self.args.replace,
                                replace_with=self.args.replace_with,
                                replace_this=self.args.replace_this,
                                stringcols=self.args.stringcols)
            if self.args.impute:
                if not self.args.quiet:
                    newmsg("imputing early stage bins with missing data")
                self.early.impute_zeros(bw=self.args.impute)
        else:
            self.early = False
        ## todo -- add filtering out 0 contigs option...
        

        if not self.args.quiet:
            if self.args.earlystage:
                emsg = ' with additional early stage normalization'
            else:
                emsg = ' without additional early stage normalization'
            
        self._set_protocol()
        newmsg("following normalization protocol "+str(self.protocol)+emsg)
        self._run_protocol()

    def _finalize(self):
        ## Optional Log Transforms
        if self.args.log10:
            self.late.log10_data()
        elif self.args.log2:
            self.late.log2_data()

        ## Optional end smoothing
        if self.args.endsmoothing:
            self.late.ksmooth_counts(bw=self.args.bandwidth,
                                      rescueNaN=self.args.replaceNaN,
                                      localWindow=self.args.localwinsize)


##TODO - allow imputing values locally when a bin is 0 -- i.e. if surrounded by 2 bins with values >0, take average.
## various 0 spots are causing short state changes in the CN hmm.
## perhaps do ksmoothing with something like 10-50kb bandwidth -- then go back to raw signal, and wherever 0 is, substitute with Ksmoothed value.
## 0 spots lead to 2 problems:
##   1. these become fe=1 when both early and late are 0 -- leads to state drops and even bad smoothing
##   2. when not 0 in late, this becomes (count+0.1)/0.1 in fe step (after pseudocount)
##    --> which could mean a really high bin
## best stage to impute would be prior to ANY normalization
## i.e. impute late and early separately --> then median norm --> then FE --> then smooth if want
## If one does protocol 2 or 3... technically imputing is done just by smoothing...
##   only problem is it still shows a significant drop at 0 spots... whereas would be less the case if pre-imputed even before this smoothing

## NOTES: imputing in this way creates its own problems.
## e.g. if late needs to be imputed, but early does not
##  then late will get value close to other values, but early will remain near 0 -- so you get massive spike
##  would want to change same bins in both samples...
##  didnt seem to change state path at all

        

def run(parser, args):
    protocol = NormalizeProtocol(args)
    sys.stdout.write(
        protocol.late.get_bdg(
            protocol.late.count,
            args.collapsed
            )
        )
