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
        if self.early:
            self.late.normalize_to_other(self.early,
                                         self.args.pseudo)
        self.late.median_normalize_data()

    def _protocol13(self):
        # median ratio norm w/ pre-smoothing - e.g. similar to DEseq (or TMM in EdgeR) -- global version
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
        # median ratio norm w/ post-smoothing - e.g. similar to DEseq (or TMM in EdgeR) -- global version
        if self.early:
            self.late.normalize_to_other(self.early,
                                         self.args.pseudo)
        self.late.ksmooth_counts(bw=self.args.bandwidth,
                                      rescueNaN=self.args.replaceNaN,
                                      localWindow=self.args.localwinsize)
        self.late.median_normalize_data()
        

    def _protocol15(self):
        # median ratio norm w/ end-smoothing - e.g. similar to DEseq (or TMM in EdgeR) -- global version
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
        ## Similar to 17, but only adds pseudo to bins with 0....
        ##
        ## First median normalizes both samples
        ## Then scales to a target median bin cov -- default 1000.
        ## Then adds pseudocount - default 1. (so default pseudo is 0.1% median, 1000-fold less than median)
        ## Then normalizes late to early
        ## Then median ratio normalizes a la DEseq/EdgeR
        pass

        
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
