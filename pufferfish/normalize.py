import sys, datetime
from CovBedClass import *
from pk2txt import bdgmsg, newmsg


class NormalizeProtocol(object):
    def __init__(self, args):
        #latestage, protocol=1, earlystage=False, pseudo=0.1, bandwidth=2500, quiet=False, impute=False, replace=False, replace_with='0', replace_this='.'
        self.args = args
        self._normalize()
        
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

    def _protocol1(self):
        self.late.median_normalize_data()
        if self.early:
            self.early.median_normalize_data()
            self.late.normalize_to_other(self.early,
                                              self.args.pseudocount)

    def _protocol2(self):
        print self.late.get_median()
        self.late.ksmooth_counts(bw=self.args.bandwidth)
        self.late.median_normalize_data()
        if self.early:
            self.early.ksmooth_counts(bw=self.args.bandwidth)
            self.early.median_normalize_data()
            self.late.normalize_to_other(self.early,
                                              self.args.pseudocount)

    def _protocol3(self):
        self.late.median_normalize_data()
        self.late.ksmooth_counts(bw=self.args.bandwidth)
        if self.early:
            self.early.median_normalize_data()
            self.early.ksmooth_counts(bw=self.args.bw)
            self.late.normalize_to_other(self.early,
                                         self.args.pseudocount)

    def _protocol4(self):
        self.late.median_normalize_data()
        if self.early:
            self.early.median_normalize_data()
            self.late.normalize_to_other(self.early,
                                         self.args.pseudocount)
        self.late.ksmooth_counts(bw=self.args.bandwidth)

    def _protocol5(self):
        #smoothing only --  late/early before smooth
        if self.early:
            self.late.normalize_to_other(self.early,
                                         self.args.pseudocount)
        self.late.ksmooth_counts(bw=self.args.bandwidth)

    def _protocol6(self):
        #smoothing only --  late/early AFTER smooth
        self.late.ksmooth_counts(bw=self.args.bandwidth)
        if self.early:
            self.early.ksmooth_counts(bw=self.args.bandwidth)
            self.late.normalize_to_other(self.early,
                                         self.args.pseudocount)

    def _protocol7(self):
        #no smoothing, no median norm
        # late:early only if early present
        if self.early:
            self.late.normalize_to_other(self.early, pseudocount)

    def _protocol8(self):
        self.late.computeSkew()

    def _protocol9(self):
        self.late.computePercentChange()

    def _protocol10(self):
        self.late.computeSkewChange()

    def _protocol11(self):
        self.late.computePercentChange()

        
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
