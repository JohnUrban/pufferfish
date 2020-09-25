import sys, datetime
from CovBedClass import *
from pk2txt import bdgmsg, newmsg
from normalize import NormalizeProtocol #protocol1, protocol2, protocol3, protocol4, protocol5, protocol6, normalize
from hmm_fxns_for_R import *



## TODO:
## ALLOW OPTION TO USE KMEANS CLUSTERING TO FIND MEAN OF K STATES (w/ or w/o normalization)
## --kmeans 5 -- this will over-ride the normal emissions params given or used and will learn them...
## might also want to give option to randomly sample from data if kmeans takes to long on all data..



###########################################################################
'''SUPPLEMENTAL FUNCTIONS'''
###########################################################################

def report(args, late, statepath, i):
    ## FIGURE OUT OUTPUT NAME OR STDOUT
    if args.outpfx is None:
        out = sys.stdout
    else:
        out = open(args.outpfx + '_iter' + str(i) + '.bedGraph', 'w')
    ## OPTION: If opted for, get the state means/levels
    if args.levels:
        newmsg("Getting levels for alternative reporting...")
        levels = get_levels(statepath)

    ## TALK
    if not args.quiet:
        bdgmsg("REPORTING STATE PATH ITER " + str(i), args.collapsed)

    ## REPORT
    if args.levels:
        out.write(late.get_bdg(levels, args.collapsed))
    else:
        out.write(late.get_bdg(statepath, args.collapsed))

    ## CLOSE FILE
    if args.outpfx is not None:
        out.close()


def report_counts(args, late):
    fn = args.counts + ".bedGraph"
    if not args.quiet:
        bdgmsg("REPORTING:: Final normalized late stage counts", args.collapsed)
        newmsg("Written to: " + fn)
    o = open(fn, 'w')
    o.write(late.get_bdg(late.count, args.collapsed))
    o.close()

def get_initial_probs(args, late):
    if not args.quiet:
        newmsg("Constructing probability matrices...")

    ## CONSTRUCT EMISSIONS/TRANSITION/INITIAL PROBABILITY MATRIXES FOR R
    if args.kmeans is not None:
        nstates = args.kmeans
        data = []
        for chrom in late.count:
            data += list(late.count[chrom])
            
        r_eprobs, r_tprobs, r_iprobs, np_eprobs, np_tprobs, np_iprobs = help_get_state_emissions_from_kmeans(data, args.kmeans)
    else:
        if args.emodel == "discrete":
            nstates = len( args.mu.strip().strip('\\').split(';') )
        else:
            nstates = len( args.mu.strip().strip('\\').split(',') )
        r_eprobs, r_tprobs, r_iprobs, np_eprobs, np_tprobs, np_iprobs = help_get_prob_matrices_from_params(args.mu,
                                                                                                           args.sigma,
                                                                                                           args.mu_scale,
                                                                                                           args.leave_special_state,
                                                                                                           args.leave_other,
                                                                                                           args.special_idx,
                                                                                                           args.init_special,
                                                                                                           args.initialprobs,
                                                                                                           args.transprobs,
                                                                                                           args.discrete)

    return nstates, r_eprobs, r_tprobs, r_iprobs, np_eprobs, np_tprobs, np_iprobs




def do_hmm_iter_steps(args, late, nstates, r_eprobs, r_tprobs, r_iprobs, np_eprobs, np_tprobs, np_iprobs, converged=1e-9):
    #
    
    ## HIDDEN MARKOV MODEL: Find most probable path through states
    ## STATE PATH LEARNING AND RETURN:
    ## This is "Viterbi Training" when using the Viterbi Decoded path.
    ## It is a modified Baum-Welch when using posterior decoding -- it uses the posterior state path as the 100% solution rather than propagating probabilities.
    ## In both cases, the current algo uses a predetermined number of iters (default 1).
    ##      ...rather than looking for a small change in log likelihood of the model
    nstates = np_tprobs.shape[0]
    log10probs = np.zeros(args.iters)
    for i in range(0,args.iters):
        if not args.quiet:
            newmsg("PARAMETERS: iter " + str(i))
            newmsg("\nTransition Probs:\n"+str(r_tprobs))
            newmsg("\nEmission Probs:\n"+str(r_eprobs))
            newmsg("\nInitial Probs:\n"+str(r_iprobs))
        
        if not args.quiet:
            newmsg("FINDING STATE PATH: iter " + str(i))

        ## STEP 1: FIND STATE PATH WITH CURRENT PARAMETERS.
        statepath = hmmR(late, args.path, args.emodel,
                         eprobs=r_eprobs, tprobs=r_tprobs, iprobs=r_iprobs)
        
        report(args, late,
               statepath, i)
        
        ## STEP 1.5: GET LOG PROB OF STATEPATH
        if not args.quiet:
            newmsg("Computing LOG10 PROB STATEPATH: iter " + str(i) + ".........")
        log10probs[i] = log10_prob_state_path(late, statepath, np_eprobs, np_tprobs, np_iprobs, args.emodel)
        if not args.quiet:
            newmsg("LOG10 PROB STATEPATH: iter " + str(i) + " = " + str(log10probs[i]))
        ## CONVERGED?                                                                                   nstates) #r_eprobs, r_tprobs, r_iprobs, np_eprobs, np_tprobs, np_iprobs)
        if args.iters > 1 and i > 0:
            abslogdiff = abs(log10probs[i]-log10probs[i-1])
            if not args.quiet:
                newmsg("DIFFERENCE FROM LAST ITER: " + " = " + str(abslogdiff))
            if abslogdiff <= converged:
                ## FILL REST OF LOG10 WITH DUMMIES OF BREAKING ITER
                for i_leftover in range(i+1, args.iters):
                    log10probs[i_leftover] = log10probs[i]
                if not args.quiet:
                    newmsg("Model has converged at iter: " + str(i) + "\nLog10 State path given model = " + str(log10probs[i]) + "\n... STOPPING. ")
                break
                
        ## STEP 2: UPDATE PARAMETERS WITH CURRENT STATE PATH.
        if args.iters > 1 and i < args.iters-1: ## Don't need to do if only 1 iter, or if last iter.
            if not args.quiet and args.iters > 1:
                newmsg("UPDATING PARAMETERS: iter " + str(i))
            ## Update parameters with current state path.
            r_eprobs, r_tprobs, r_iprobs, np_eprobs, np_tprobs, np_iprobs = updateparameters(late,
                                                                                             statepath,
                                                                                             nstates,
                                                                                             old_np_eprobs=np_eprobs,
                                                                                             learnpseudo=args.learnpseudo,
                                                                                             emodel=args.emodel,
                                                                                             emitpseudo=args.emitpseudo,
                                                                                             constrainEmit=args.constrainEmit)
        
    return log10probs

            







###########################################################################
'''RUN/MAIN FUNCTION'''
###########################################################################

def run(parser, args):
    ## Don't have kmeans for the discrete option
    assert not (args.emodel == "discrete" and args.kmeans is not None)

    ## Set args.discrete variable based on emodel
    args.discrete = True if args.emodel == "discrete" else False
    ##sys.stderr.write("DEBUG: " + str(args.discrete)+"\n")

    
    ## Warning about iterations and output prefixes
    if args.iters > 1 and args.outpfx is None:
        newmsg("\n\tWARNING WARNING WARNING!!!!!!!!!!!!\n\tIters > 1, but no output prefix specified!\n\tIgnore if that was intentional.\n\tElse restart with --outpfx!\n\tWARNING WARNING WARNING!!!!!!!!!!!!")

    ## Transform data with specified normalization protocol (incl no normalization)
    protocol = NormalizeProtocol(args)
    late = protocol.late #normalize(latestage=args.latestage, protocol=protocol, earlystage=args.earlystage, pseudo=args.pseudo, bandwidth=args.bandwidth, quiet=args.quiet, impute=args.impute, replace=args.replace, replace_with=args.replace_with, replace_this=args.replace_this)

    ## OPTIONAL REPORTING OF TRANSFORMED DATA
    if args.counts:
        report_counts(args, late)

    ## INITIALIZE PARAMETERS
    nstates, r_eprobs, r_tprobs, r_iprobs, np_eprobs, np_tprobs, np_iprobs = get_initial_probs(args, late)

    ## ITERATE
    log10probs = do_hmm_iter_steps(args, late, nstates,
                      r_eprobs, r_tprobs, r_iprobs,
                      np_eprobs, np_tprobs, np_iprobs,
                                   args.converge)

    ## REPORTING LOGPROBS FOR ALL ITERS
    for i in range(args.iters):
        newmsg(str(i) + "\t" + str(log10probs[i]))

    
        
        
        
