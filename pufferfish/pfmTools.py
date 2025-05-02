import random, sys
import pandas as pd
import numpy as np
from scipy.stats import geom
from Bio import SeqIO

def killmsg():
    print "Kill Error: Reached a sequence of different length than previous sequences."
    print "Sequence length established as:", seqLen
    print "Sequence name:", currentHeader
    print "Sequence number:", seqCount
    print "Sequence length:", len(sequence)
    quit()
    
def addToPfm(seqLen, sequence, posFreqMatDict, seqCount):
    if len(sequence) == seqLen:
        sequence = sequence.upper() ## make all bases same case
        ## for each position along sequence
        for i in range(seqLen):
            ## Take the base at position i
            b_i = sequence[i]
            try:
                ## Add 1 count at position i for b_i = b
                posFreqMatDict[i][b_i] += 1
            except KeyError:
                try: ## perhaps position i is there, but not yet this base
                    posFreqMatDict[i][b_i] = 1
                except KeyError: ## this means position i was not yet there
                    posFreqMatDict[i] = {}
                    posFreqMatDict[i][b_i] = 1
            
    else:
        killmsg()

    return posFreqMatDict






def seqlencheck(seqLen, sequence, posFreqMatDict, pseudo=0):
    if seqLen == None:
        ## This should only occur once
        ## establish seqLen 
        seqLen = len(sequence)
        ## establish what the minimum dictionary should be
        for i in range(seqLen):
            posFreqMatDict[i] = {}
            for b in 'ACGTN':
                posFreqMatDict[i][b] = pseudo
    return posFreqMatDict, seqLen




def pfmDictFromFasta(fastaFile):
    """Fasta File should contain N sequences all of the same length.
    If a sequence of a length != to the length of the first sequence analyzed is encountered,
    it will raise an error.
    Note, there should be no empty lines in the fasta file except for the end of file.
    Use: "grep -v ^$ fasta file > new.fa" to fix if you need."""

    ## open connection to the fasta file
    fasta = open(fastaFile, 'r')

    ## initialize a variable to establish sequence length
    seqLen = None

    ## initialize a variable to store current sequences
    sequence = ''

    ## initalize the PFM dictionary
    posFreqMatDict = {}

    ## Do line#1 outside of for loop
    currentHeader = fasta.readline().strip()[1:]
    seqCount = 1

    ## For each line in fasta, determine if it is an empty line (end of file), a header line, or a sequence line and do appropriate actions.
    for line in fasta:
        if line.startswith('>') == False: ## sequence line
            sequence += line.strip() ## Add all but \n at end
        else: ##line starts with '>'
            ## Then there must be a sequence that has just been put together
            ## analyze previous sequence before moving on

            ## Make sure stored sequence is same size as all previous sequences (or establish the seqLen now)
            posFreqMatDict, seqLen = seqlencheck(seqLen, sequence, posFreqMatDict)

            # Update
            posFreqMatDict = addToPfm(seqLen, sequence, posFreqMatDict, seqCount)
            currentHeader = line[1:-1]
            seqCount += 1
            sequence = ''
    # Final update
    posFreqMatDict = addToPfm(seqLen, sequence, posFreqMatDict, seqCount)
    fasta.close()
    return posFreqMatDict




def pfmDictFromSeqList(seqlist, pseudo=0):
    """seqlist is a list object that should contain N sequences all of the same length.
    If a sequence of a length != to the length of the first sequence analyzed is encountered,
    it will raise an error."""

    ## initialize a variable to establish sequence length
    seqLen = None

    ## initialize a variable to store current sequences
    sequence = ''

    ## initalize the PFM dictionary
    posFreqMatDict = {}

    ## Init
    seqCount = 1

    ## Iter.
    for sequence in seqlist:

        ## Make sure stored sequence is same size as all previous sequences (or establish the seqLen now)
        posFreqMatDict, seqLen = seqlencheck(seqLen, sequence, posFreqMatDict, pseudo)

        # Update
        posFreqMatDict = addToPfm(seqLen, sequence, posFreqMatDict, seqCount)
        seqCount += 1
    return posFreqMatDict



def constructPwmDFfromSeqList(seqlist,  distributeN=False, marginalizeOutN=False, pseudo=0):
    posFreqMatDict = pfmDictFromSeqList(seqlist, pseudo)
    return pfmDictToPwmDF(posFreqMatDict, distributeN, marginalizeOutN)
    

def pfmDictToPfmDF(posFreqMatDict):
    return pd.DataFrame(posFreqMatDict).transpose()


def pfmDictToPwmDF(posFreqMatDict, distributeN=False, marginalizeOutN=False):
    ''' distributeN and marginalizeOutN are mutually exclusive'''
    # Check
    assert not (distributeN and marginalizeOutN)

    # Get
    freqs = pfmDictToPfmDF(posFreqMatDict)
    
    # Pre-process
    if distributeN:
        nfrac = freqs['N'].map(lambda x: x/4.0)
        Nfrac = pd.DataFrame( {e:nfrac for e in 'ACGT'} ) ## will add this much to each column
        
    if  marginalizeOutN or distributeN:
        freqs.drop('N', axis=1, inplace=True)

    if distributeN:
        freqs += Nfrac

    # Return - (z.transpose() / z.sum(axis=1)).transpose()
    #return freqs.divide( freqs.sum(axis=1) )
    return (freqs.transpose() / freqs.sum(axis=1)).transpose()



def pfmFromPfmDict(posFreqMatDict):
    """Takes in PFM dictionary object -- e.g. output of pfmFromFasta
    Writes it out to file ...."""
    
    ##1. Find the full set of symbols
    bases = set()
    for i in posFreqMatDict:
        for base in posFreqMatDict[i].keys():
            bases.add(base) ## "bases" is a set variable -- thus, if the element 'base' is in the set, it will not be repeated

    ##2. Ensure that each position in dictionary has all bases:
    for i in posFreqMatDict:
        for b in bases:
            try:
                posFreqMatDict[i][b] += 0
            except KeyError:
                posFreqMatDict[i][b] = 0
            
    ## 3. Write out header
    header = ''
    orderedBases = ''
    for e in bases:
        orderedBases += e
        header += e+"\t"
    header = header[:-1]
    print header

    ## 4. write out such that each line is equal to current position
    for i in range(len(posFreqMatDict.keys())):
        line = ''
        for b in orderedBases:
            line += str(posFreqMatDict[i][b])+"\t"
        line = line[:-1]
        print line

def pfmConsensus(posFreqMatDict, verbose=False):
    ''' Returns max of A,C,G,T -- if tie, it randomly selects one'''
    bases = ['A','C','G','T']
    baseDict = {0:'A',1:'C',2:'G',3:'T'}
    cSeq = ''
    numRands = 0
    randIndexes = []
    for i in range(len(posFreqMatDict.keys())):
        options = [posFreqMatDict[i]['A'], posFreqMatDict[i]['C'], posFreqMatDict[i]['G'], posFreqMatDict[i]['T']]
        maxOpt = max(options)
        indexes = [k for k, j in enumerate(options) if j == maxOpt]
        if len(indexes) > 1:
            indexes = [indexes[random.randint(0,len(indexes)-1)]]
            numRands += 1
            randIndexes.append(i)
        cSeq += baseDict[indexes[0]]
        
    sys.stdout.write(cSeq+"\n")
    if verbose:
        sys.stderr.write(str(numRands) + " positions had ties for which a base was randomly selected.\n")
        sys.stderr.write("Positions: " + (", ").join([str(e) for e in randIndexes]) + '\n')
    
    
#def stackMultAln(fastaFile, blocksize=50):






def transseqlencheck(seqLen, sequence, posTransFreqMatDict):
    if seqLen == None:
        ## This should only occur once
        ## establish seqLen 
        seqLen = len(sequence)
        ## establish what the minimum dictionary should be
        transdf = pd.DataFrame({b:np.zeros(5) for b in 'ACGTN'}, index='A C G T N'.split())
        for i in range(seqLen):
            posTransFreqMatDict[i] = transdf.copy()
    return posTransFreqMatDict, seqLen

def ptfmDictFromSeqList(seqlist):
    """
    ptfm = position transition frequency matrix
    
    seqlist is a list object that should contain N sequences all of the same length.
    If a sequence of a length != to the length of the first sequence analyzed is encountered,
    it will raise an error."""
    ## initialize a variable to establish sequence length
    seqLen = None

    ## initialize a variable to store current sequences
    sequence = ''

    ## initalize the PFM dictionary
    posTransFreqMatDict = {}

    ## Init
    seqCount = 1

    ## Iter.
    for sequence in seqlist:

        ## Make sure stored sequence is same size as all previous sequences (or establish the seqLen now)
        posTransFreqMatDict, seqLen = seqlencheck(seqLen, sequence, posTransFreqMatDict)

        # Update
        posFreqMatDict = addToPtfm(seqLen, sequence, posTransFreqMatDict, seqCount)
        seqCount += 1
    return posTransFreqMatDict






















def addToPtfm(seqLen, sequence, posTransFreqMatDict, letters):
    if len(sequence) == seqLen:
        N_in_ltrs = 'N' in letters
        sequence = sequence.upper() ## make all bases same case
        ## for each position along sequence
        for i in range(1, seqLen):
            ## Take the base at position i
            b_prev = sequence[i-1]
            b_i = sequence[i]

            try:
                ## Add 1 count
                if b_i == 'N' and not N_in_ltrs:
                    for b in letters:
                        posTransFreqMatDict[i][b].loc[b_prev] += 0.25
                elif b_prev == 'N' and not N_in_ltrs:
                    for b in letters:
                        posTransFreqMatDict[i][b_i].loc[b] += 0.25
                else:
                    posTransFreqMatDict[i][b_i].loc[b_prev] += 1
            except:
                print("ERROR IN addToPtfm")
                quit()
            
    else:
        killmsg()
    return posTransFreqMatDict


def transseqlencheck(seqLen, sequence, posTransFreqMatDict, letters, pseudo):
    if seqLen == None:
        ## This should only occur once
        ## establish seqLen 
        seqLen = len(sequence)
        ## establish what the minimum dictionary should be
        transdf = pd.DataFrame({b:np.zeros(len(letters))+pseudo for b in letters}, index=letters)
        for i in range(1, seqLen):
            posTransFreqMatDict[i] = transdf.copy()
    ## Else it just passes variables back through.
    return posTransFreqMatDict, seqLen




def ptfmDictFromSeqList(seqlist, includeN=True, includeDash=False, pseudo=0):
    """
    ptfm = position transition frequency matrix

    Will collect trans probs for range(1,N) -- position 0 has impicit trans from nothing as p=1
    
    seqlist is a list object that should contain N sequences all of the same length.
    If a sequence of a length != to the length of the first sequence analyzed is encountered,
    it will raise an error."""

    # Define letters
    letters = 'A C G T'.split()
    if includeN:
        letters.append( 'N' )
    if includeDash:
        letters.append( '-' )

    ## initialize a variable to establish sequence length
    seqLen = None

    ## initialize a variable to store current sequences
    sequence = ''

    ## initalize the PFM dictionary
    posTransFreqMatDict = {}


    ## Iter.
    for sequence in seqlist:

        ## Make sure stored sequence is same size as all previous sequences (or establish the seqLen now)
        posTransFreqMatDict, seqLen = transseqlencheck(seqLen, sequence, posTransFreqMatDict, letters, pseudo)

        # Update
        posFreqMatDict = addToPtfm(seqLen, sequence, posTransFreqMatDict, letters)

    # End of def
    return posTransFreqMatDict


def ptwmDictFromSeqList(seqlist, includeN=True, includeDash=False, pseudo=0):
    posTransFreqMatDict = ptfmDictFromSeqList(seqlist, includeN, includeDash, pseudo)
    for i in range(1, max(posTransFreqMatDict.keys())+1):
        posTransFreqMatDict[i] = posTransFreqMatDict[i] / posTransFreqMatDict[i].sum()
    return posTransFreqMatDict






def define_regex_from_pwm(pwm):
    regex = ''
    letters = list(pwm.columns)
    for i in range(pwm.shape[0]):
        pos='['
        for ltr in letters:
            if pwm[ltr][i] > 0:
                if ltr == 'N':
                    pos = '[AaCcGgTtN'
                    break
                else:
                    pos += ltr.upper() + ltr.lower()
        pos += ']'
        regex += pos
    return regex





class tripartiteProfileModel(object):
    '''
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
            
            
        

        '''
    def __init__(self, trainfile, distributeN=False, marginalizeOutN=False, includeN=True, includeDash=False, pseudo=1e-2):
        ## Store filename
        self.fh = trainfile

        ## Create dataframe
        self.train = pd.read_csv(trainfile, sep="\t")

        ## LHS and RHS lengths
        self.lhs_len = len(self.train['left'][0])
        self.rhs_len = len(self.train['right'][0])
        
        ## Create PWMs for LHS and RHS Emission Prob
        self.emit_lhs = constructPwmDFfromSeqList(list(self.train['left']),  distributeN, marginalizeOutN, pseudo)
        self.emit_rhs = constructPwmDFfromSeqList(list(self.train['right']),  distributeN, marginalizeOutN, pseudo)

        ## Define LHS and RHS log10 emission prob sums
        self.train['lhs_emit_prob'] = self.train['left'].map(lambda x: np.array([
            np.log10(self.emit_lhs[x[i]][i]) for i in range(self.emit_lhs.shape[0])]).sum())
        self.train['rhs_emit_prob'] = self.train['right'].map(lambda x: np.array([
            np.log10(self.emit_rhs[x[i]][i]) for i in range(self.emit_rhs.shape[0])]).sum())
        

        ## Create transition prob matrices for LHS and RHS
        self.trans_lhs = ptwmDictFromSeqList(list(self.train['left']), includeN, includeDash, pseudo)
        self.trans_rhs = ptwmDictFromSeqList(list(self.train['right']), includeN, includeDash, pseudo)

        ## Define LHS trans prob product - a[1]['G'].loc['T']
        self.train['lhs_trans_prob'] = self.train['left'].map(lambda x: np.array([
            np.log10(self.trans_lhs[i][x[i]].loc[x[i-1]]) for i in range(1, max(self.trans_lhs.keys())+1)]).sum())
        self.train['rhs_trans_prob'] = self.train['right'].map(lambda x: np.array([
            np.log10(self.trans_rhs[i][x[i]].loc[x[i-1]]) for i in range(1, max(self.trans_rhs.keys())+1)]).sum())

        
        ## Define Length probs of middle seqs - planned to model with geom, but just doing discrete freqs found in train now
        self.train['mid_seq_len'] = self.train['middle'].map(lambda x: 0 if x == "-" else len(x) )
        self.mid_seq_mean_len = self.train['mid_seq_len'].mean()
        self.mid_seq_min_len = self.train['mid_seq_len'].min()
        self.mid_seq_max_len = self.train['mid_seq_len'].max()
        self.min_motif_len = self.lhs_len + self.mid_seq_min_len + self.rhs_len
        self.max_motif_len = self.lhs_len + self.mid_seq_max_len + self.rhs_len
        
        self.mid_seq_length_probs = np.zeros(self.train['mid_seq_len'].max()+1) + pseudo
        for l in self.train['mid_seq_len']:
            self.mid_seq_length_probs[l] += 1
        self.mid_seq_length_probs = self.mid_seq_length_probs / self.mid_seq_length_probs.sum()
        self.train['mid_seq_len_prob'] = self.train['mid_seq_len'].map(lambda x: np.log10(self.mid_seq_length_probs[x] ))
        #self.mid_seq_len_prob = 1/float( self.mid_seq_mean_len )
        #self.train['mid_seq_len_prob'] = self.train['mid_seq_len'].map(lambda x: np.log10(geom.pmf(x, self.mid_seq_len_prob)) )


        ## Define composition probs of middle seqs:
        ##      - get count of all symbols (ACGTN-)
        self.midsymbolcounts = pd.DataFrame({b:self.train['middle'].map(lambda x: x.count(b)) for b in 'A C G T N -'.split(' ')})
        ##      - define probs by freq/sum(freqs)
        self.midsymbolprobs = self.midsymbolcounts.sum() / self.midsymbolcounts.sum().sum()
        ##      - get mean prob of sequence
        self.train['mid_seq_comp_prob'] = self.train['middle'].map(lambda x: np.log10(np.array([self.midsymbolprobs[e] for e in x]).mean()))
        
        
        ## Define transition prob from LHS to mid-seq and from mid-seq to RHS:
        lhs_mid = []
        mid_rhs = []
        for i in self.train.index:
            lhs_mid.append( self.train['left'][i][-1] +  self.train['middle'][i][0] )
            mid_rhs.append( self.train['middle'][i][-1] +  self.train['right'][i][0] )
        self.train['lhs_mid'] = lhs_mid
        self.train['mid_rhs'] = mid_rhs
        self.trans_lhs_mid = ptwmDictFromSeqList(list(self.train['lhs_mid']), includeDash=True, includeN=includeN, pseudo=pseudo)
        self.trans_mid_rhs = ptwmDictFromSeqList(list(self.train['mid_rhs']), includeDash=True, includeN=includeN, pseudo=pseudo)
        self.train['lhs_mid_trans_prob'] = self.train['lhs_mid'].map(lambda x: np.array([
            np.log10(self.trans_lhs_mid[i][x[i]].loc[x[i-1]]) for i in range(1, max(self.trans_lhs_mid.keys())+1)]).sum())
        self.train['mid_rhs_trans_prob'] = self.train['mid_rhs'].map(lambda x: np.array([
            np.log10(self.trans_mid_rhs[i][x[i]].loc[x[i-1]]) for i in range(1, max(self.trans_mid_rhs.keys())+1)]).sum())


        ## Define total seq prob - sum of log10 probs
        bipartite_prob = []
        tripartite_prob = []
        indiv_bi_probs = ['lhs_emit_prob', 'rhs_emit_prob', 'lhs_trans_prob', 'rhs_trans_prob']
        indiv_tri_probs = indiv_bi_probs + ['mid_seq_len_prob', 'mid_seq_comp_prob', 'lhs_mid_trans_prob', 'mid_rhs_trans_prob']
        for i in self.train.index:
            bipartite_prob.append( np.array([self.train[p][i] for p in indiv_bi_probs]).sum()  )
            tripartite_prob.append( np.array([self.train[p][i] for p in indiv_tri_probs]).sum()  )
        self.train['bipartite_prob'] = bipartite_prob
        self.train['tripartite_prob'] = tripartite_prob

        ## Define the regular expression





class triPartiteSequenceSearch(object):
    ''' To search the Reverse Complement Strand, simply provide it as the sequence, and use the reverse option for writing Wigs and Beds.'''
    def __init__(self, sequence, name, model, bipcutoff='auto', tripcutoff='auto', minMidSeqLen=0, verbose=True):
        ''' cutoffs can be any number or the following strings: auto, min, max, mean, median. '''
        ## VARS
        self.sequence = sequence
        self.name = name
        self.model = model
        self.seqlen = len(sequence)
        self.store_lhs = {}
        self.store_rhs = {}
        self.store_mid_comp = {}
        self.store_lhs_mid = {}
        self.store_mid_rhs = {}
        self.store_mid = {}
        self.lhs_likelihood = np.zeros(self.seqlen-self.model.lhs_len+1)
        self.rhs_likelihood = np.zeros(self.seqlen-self.model.rhs_len+1)
        self.bip_likelihood = np.zeros(self.seqlen - self.model.max_motif_len + 1)
        self.trip_likelihood = np.zeros(self.seqlen - self.model.max_motif_len + 1)
        self.bipcalls = []
        self.tripcalls = []
        self.progress = 0
        self.updateprog = 0.1
        self.verbose = verbose
        self.minMidSeqLen = self.model.mid_seq_min_len if minMidSeqLen == "auto" else int(minMidSeqLen)
        
        
        ## EXEC
        ## Todo -- can compute LHS and RHS in one pass...
        self.__stderr("Computing LHS likelihoods.")
        self.__lhs_likelihood()
        self.__stderr("Computing RHS likelihoods.")
        self.__rhs_likelihood()
        self.__stderr("Computing Bipartite/Tripartite likelihoods.")
        # Get cutoffs
        self.bipcutoff = bipcutoff
        self.tripcutoff = tripcutoff
        if bipcutoff in ('auto', 'min', 'max', 'mean', 'median'):
            if bipcutoff in ('auto', 'min'):
                self.bipcutoff = self.model.train['bipartite_prob'].min()
            elif bipcutoff == 'max':
                self.bipcutoff = self.model.train['bipartite_prob'].max()
            elif bipcutoff == 'mean':
                self.bipcutoff = self.model.train['bipartite_prob'].mean()
            elif bipcutoff == 'median':
                self.bipcutoff = np.median(self.model.train['bipartite_prob'])
        if tripcutoff in ('auto', 'min', 'max', 'mean', 'median'):
            if tripcutoff in ('auto', 'min'):
                self.tripcutoff = self.model.train['tripartite_prob'].min()
            elif tripcutoff == 'max':
                self.tripcutoff = self.model.train['tripartite_prob'].max()
            elif tripcutoff == 'mean':
                self.tripcutoff = self.model.train['tripartite_prob'].mean()
            elif tripcutoff == 'median':
                self.tripcutoff = np.median(self.model.train['tripartite_prob'])

        self.__bi_and_tripartite_likelihoods(self.bipcutoff, self.tripcutoff)
    def __stderr(self, msg):
        if self.verbose:
            sys.stderr.write('tripartiteSequenceSearch INFO ::: ' + msg + '\n')
    def to_stderr(self, msg):
        self.__stderr(msg)
    def __progmsg(self, j):
        pct = j/float(self.seqlen) 
        if j/float(self.seqlen)  > self.progress:
            self.__stderr(str(round(100*pct)) + "% of " + str(self.seqlen) + "bp sequence processed...")
            self.progress += self.updateprog
        if self.progress >= 1:
            self.progress=0
    def __side_likelihood(self, likelihood, sidelen, emitmat, transmat, store):
        for j in range(self.seqlen-sidelen+1):
            self.__progmsg(j)
                
            x = self.sequence[j:j+sidelen]

            ## try using stored value; if not stored, then calculate it and store
            try:
                likelihood[j] = store[x]
            except:
                ## Emissions
                store[x] = (
                    np.array([
                        np.log10(emitmat[x[i]][i])
                        for i in range(sidelen)
                        ]).sum()
                    )
                ## ADD Transitions
                store[x] += (
                    np.array([
                        np.log10(transmat[i][x[i]].loc[x[i-1]])
                        for i in range(1, max(transmat.keys())+1)
                        ]).sum()
                    )
                likelihood[j] = store[x]
        ## Reset progress
        self.__stderr("Finished step...")
        self.progress = 0
        ###return likelihood
            
    def __lhs_likelihood(self):
        self.__side_likelihood(likelihood=self.lhs_likelihood,
                                sidelen=self.model.lhs_len,
                                emitmat=self.model.emit_lhs,
                                transmat=self.model.trans_lhs,
                                store=self.store_lhs)
    
    def __rhs_likelihood(self):
        self.__side_likelihood(likelihood=self.rhs_likelihood,
                                sidelen=self.model.rhs_len,
                                emitmat=self.model.emit_rhs,
                                transmat=self.model.trans_rhs,
                                store=self.store_rhs)
    
    def __bi_and_tripartite_likelihoods(self, bipcutoff, tripcutoff, ndigits=2):
        ''' For each position, gather the:
            max bipartite likelihood array
            max bipartite BED coords if passes cutoff.
            max tripartite likelihood array
            max tripartite BED coords if passes cutoff'''

        ## Bip and Trip can be calculated from sum of those 2 or 3 stored values
        for i in range(self.seqlen - self.model.max_motif_len + 1):
            self.__progmsg(i)
            #mid_rhs = float('-inf')
            lhs = self.lhs_likelihood[i]      # lhs + mid + rhs
            rhs = float('-inf')
            trip = float('-inf')
            # identify max relative rhs pos
            for midlen in range(self.minMidSeqLen, self.model.mid_seq_max_len + 1): ## min is 0 by default, minlen seen if auto specified, or w/e int given
                cand_rhs_pos = i + self.model.lhs_len + midlen
                cand_rhs = self.rhs_likelihood[cand_rhs_pos]
                if cand_rhs > rhs:
                    rhs = cand_rhs
                    bip = lhs + rhs
                    biplen = cand_rhs_pos + self.model.rhs_len - i

                ## also calculate tripartite likelihood
                ## GET MID PROB
                ## Define Middle sequence
                mid_seq = "-" if midlen == 0 else self.sequence[(i + self.model.lhs_len):cand_rhs_pos]
                ## Define Middle sequence flanked by single bases from LHS and RHS
                last_mid_base = "-" if midlen == 0 else self.sequence[cand_rhs_pos - 1]
                first_rhs_base = "-" if midlen == 0 else self.sequence[cand_rhs_pos]
                last_lhs_base = self.sequence[i + self.model.lhs_len - 1]
                first_mid_base = "-" if midlen == 0 else self.sequence[i + self.model.lhs_len]
                flanked_mid_seq = last_lhs_base + mid_seq + first_rhs_base

                ## Check if this mid value is already stored. If not compute:
                try:
                    mid = self.store_mid[flanked_mid_seq]
                except:
                    ## Get prob of trans from LHS to Mid
                    x = last_lhs_base + first_mid_base
                    try:
                        lhs_to_mid = self.store_lhs_mid[x]
                    except:    
                        self.store_lhs_mid[x] = np.array([
                            np.log10(self.model.trans_lhs_mid[k][x[k]].loc[x[k-1]])
                            for k in range(1, max(self.model.trans_lhs_mid.keys())+1)
                            ]).sum()
                        lhs_to_mid = self.store_lhs_mid[x]
                    ## Get prob of mid seq len (already stored)
                    mid_len = np.log10(self.model.mid_seq_length_probs[midlen])

                    ## Get prob of mid seq nt comp
                    try:
                        mid_comp = self.store_mid_comp[mid_seq]
                    except:
                        self.store_mid_comp[mid_seq] = np.log10(np.array([
                            self.model.midsymbolprobs[e] for e in x
                            ]).mean() + 1e-100) ### PSEUDO ADDED HERE TO AVOID DIV ZERO ERROR.............................................
                        mid_comp = self.store_mid_comp[mid_seq]

                    ## Get prob of trans from Mid to RHS
                    x = last_mid_base + first_rhs_base
                    try:
                        mid_to_rhs = self.store_mid_rhs[x]
                    except:
                        self.store_mid_rhs[x] = np.array([
                            np.log10(self.model.trans_mid_rhs[k][x[k]].loc[x[k-1]])
                            for k in range(1, max(self.model.trans_mid_rhs.keys())+1)
                            ]).sum()
                        mid_to_rhs = self.store_mid_rhs[x]
                    ## Put it all together for the midseq prob
                    self.store_mid[flanked_mid_seq] = lhs_to_mid + mid_len + mid_comp + mid_to_rhs
                    mid = self.store_mid[flanked_mid_seq]

                ## CHECK
                cand_trip = lhs + mid + rhs
                if cand_trip > trip:
                    trip = cand_trip
                    triplen = cand_rhs_pos + self.model.rhs_len - i


            # store max bip
            self.bip_likelihood[i] = round(bip, ndigits=ndigits)
            # check if bip exceeds cutoff, and store if so
            if bip >= bipcutoff:
                # collect bed coords
                self.bipcalls.append( (i,i+biplen) )
            # store max trip
            self.trip_likelihood[i] = round(trip, ndigits=ndigits)
            # check if trip exceeds cutoff, and store if so
            if trip >= tripcutoff:
                # collect bed coords
                self.tripcalls.append( (i,i+triplen) )
        ## END FXN
        ## Reset progress
        self.__stderr("Finished step...")
        self.progress = 0

    def wig_header(self, start=1):             
        return 'fixedStep chrom=' + self.name + ' start=' + str(start) + ' step=1\n'

    def lhs_wig(self, reverse=False, start=1):
        if reverse:
            return self.wig_header(start=start + self.model.lhs_len - 1) + '\n'.join([str(e) for e in np.flip(self.lhs_likelihood)])              
        else:
            return self.wig_header(start=start) + '\n'.join([str(e) for e in self.lhs_likelihood])
    
    def rhs_wig(self, reverse=False, start=1):
        if reverse:
            return self.wig_header(start=start + self.model.rhs_len - 1) + '\n'.join([str(e) for e in np.flip(self.rhs_likelihood)])
        else:
            return self.wig_header(start=start) + '\n'.join([str(e) for e in self.rhs_likelihood])

    def bipartite_wig(self, adjust=True, reverse=False, start=1):
        add = abs(self.bipcutoff) if adjust else 0
        if reverse:
            return self.wig_header(start=start + self.model.max_motif_len - 1) + '\n'.join([str(e) for e in np.flip(self.bip_likelihood + add)])
        else:
            return self.wig_header(start=start) + '\n'.join([str(e) for e in self.bip_likelihood + add])

    def tripartite_wig(self, adjust=True, reverse=False, start=1):
        add = abs(self.tripcutoff) if adjust else 0
        if reverse:
            return self.wig_header(start=start + self.model.max_motif_len - 1) + '\n'.join([str(e) for e in np.flip(self.trip_likelihood + add)])
        else:
            return self.wig_header(start=start) + '\n'.join([str(e) for e in self.trip_likelihood + add])
        #return self.wig_header() + '\n'.join([str(e) for e in 10**self.trip_likelihood])

    def bipartite_bedgraph(self, adjust=True, reverse=False, start=0):
##        str(self.seqlen - self.bipcalls[i][1] + 1),
##         str(self.seqlen - self.bipcalls[i][0] + 1)
        add = abs(self.bipcutoff) if adjust else 0
        
        if reverse:
            return ('\n').join([
                '\t'.join([
                    self.name,
                    str(start + self.seqlen - self.bipcalls[i][0] + 1),
                    str(start + self.seqlen - self.bipcalls[i][0] + 1 + self.bipcalls[i][1]-self.bipcalls[i][0]),
                    str(np.flip(self.bip_likelihood)[ self.seqlen -self.model.max_motif_len - self.bipcalls[i][0]  ]  + add )
                    ]) 
                for i in range(len(self.bipcalls))
                ]) 
        else:
            return ('\n').join([
                '\t'.join([
                    self.name,
                    str(start + self.bipcalls[i][0]),
                    str(start + self.bipcalls[i][1]),
                    str(self.bip_likelihood[ self.bipcalls[i][0] ] + add)
                    ])
                for i in range(len(self.bipcalls))
                ]) 

    def tripartite_bedgraph(self, adjust=True, reverse=False, start = 0):
        add = abs(self.tripcutoff) if adjust else 0
        if reverse:
            return ('\n').join([
                '\t'.join([
                    self.name,
                    str(start + self.seqlen - self.tripcalls[i][0] + 1),
                    str(start + self.seqlen - self.tripcalls[i][0] + 1 + self.tripcalls[i][1]-self.tripcalls[i][0]),
                    str(np.flip(self.trip_likelihood)[ self.seqlen -self.model.max_motif_len - self.tripcalls[i][0]  ]  + add)
                    ]) 
                    
                for i in range(len(self.tripcalls))
                ]) 
        else:
            return ('\n').join([
                '\t'.join([
                    self.name,
                    str(start + self.tripcalls[i][0]),
                    str(start + self.tripcalls[i][1]),
                    str(self.trip_likelihood[ self.tripcalls[i][0] ]  + add)
                    ])
                for i in range(len(self.tripcalls))
                ]) 

    

## Uncomment for development purposes.
#self = tripartiteProfileModel('/Users/johnurban/Documents/data/sciara/EcRE/Table1-ver5-April16-MikeF-or-Yutaka-EcRE-analysis.transcribed.txt')
#seq = triPartiteSequenceSearch('GAGGTCAATGACCTCAAAAAAAAAAGGGTGCGATGAATCAGGGGGGGGGGGGGG', 'exampleseq', self)

#self = tripartiteProfileModel('/Users/johnurban/GoogleDrives/Carnegie/spradling_lab/mypapers/miiko/yeastnsseq/data/oridb/ACS/tripartite/acs-1992/train-table-1992.txt')
#seq = triPartiteSequenceSearch('GATTACATTTATATTTAGATTACA', 'exampleseq', self)
