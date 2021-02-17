# -*- coding: utf-8 -*-
"""

Created on Wed May 13 11:42:50 2020

@author: barthele

Functions :
    dosimuls
    simulate
    sample
    params
    getAbund

"""
    # TODO : doc !!!
import msprime
import random
import numpy as np
import sys
from ete3 import Tree
import pandas as pd
#from loguniform import LogUniform
from scipy.stats import loguniform

from . import pastdemo
from . import islmodel
from . import toPhylo

def dosimuls(nsim, sample_size, comprior, muprior, lim_mrca = None, sstype="SFS",
             prior_distrib = "uniform", npop=1, withmigr=False, init_ratesprior=None,
             init_sizeprior=None, pastprior=None, changetime = None,
             nsplit=None, massprior=None, migrfrom=None, migrto=None,
             verbose=False, savetrees= False, saveto = "", seed = None):

    # TODO : doc !!!
    # TODO : idiotproof
 # CHECKS HERE FOR IDIOT-PROOFING
    
    if prior_distrib not in ("uniform", "log_unif"):
        sys.exit("Ecophylo only supports uniform or log-uniform prior distributions for the moment")
    
    comments = 'Simulating eco-evolutionary dynamics over the following parameters ranges:'
    df = pd.DataFrame()

    # CREATING PARAMETER DATAFRAME
    
    # sample size
    if isinstance(sample_size, list):
        df['sampsize'] = ""
        comments += (f'\n -sample size in {sample_size}')
    else:
        samp = sample_size
    
    # mutation rate
    if muprior is not None:
        df['mu'] = ""
        comments += (f'\n -mutation rate in {muprior}')
   
    # community size (Ne)
    if comprior is not None:
        df['comsize'] = ""
        comments += (f'\n -community size in {comprior}')
    
    # past size(s)
    if pastprior is not None:
        nepoch = len(changetime)
        colnames = ['pastsize{}'.format(i) for i in range(1, nepoch+1)]
        for i in range(nepoch):
            df[colnames[i]] = ""
        comments += (f'\n -{nepoch} past epoch(s) in which community size varies in {pastprior}')
    else:
        past_sizes = None

    # populations sizes and growth rates
    if init_sizeprior is not None and init_ratesprior is not None and npop is not None:
        comments += (f'\n -{npop} population(s) in which initial community sizes vary in {init_sizeprior} and growth rates vary in {init_ratesprior}')
        colnames = ['initrate{}'.format(i) for i in range(1, npop+1)]
        for i in range(npop):
            df[colnames[i]] = ""
        colnames = ['initsize{}'.format(i) for i in range(1, npop+1)]
        for i in range(npop):
            df[colnames[i]] = ""
    else:
        init_sizes = None
        init_rates = None

    # migration between populations
    if withmigr and npop is not None:
        df['m'] = ""
        comments += ('\n -there is migration between populations')
    else:
        m = 0

    # past split dates
    if massprior is not None:
        colnames = ['splitdate{}'.format(i) for i in range(1, nsplit+1)]
        for i in range(nsplit):
            df[colnames] = ""
        comments += (f'\n -populations {migrfrom} and {migrto} split at dates ranging in {massprior}')
    else:
        split_dates = None

    comments += f"\nEstimating the associated model's {len(df.columns)} parameter(s) based on the {sstype} summary statistic(s)"
    
    if savetrees and saveto is not None:
        comments += f"\n  //!\\\ Phylogenies will be saved to the specified file"
        trees = ""

    if verbose:
        print(comments)
    
    # FILL PARAMETER DATAFRAME AND SIMULATE PHYLOGENIES
    safecount = 0
    i = 0
    ss = list()

    # While loop for calling simulate
    while i < nsim:
        if safecount > 100:
            sys.exit("Too many simulations have failed")

        try:
            df = df.append(pd.Series(), ignore_index=True)
            # sample parameters from prior for simulation

            if isinstance(sample_size, list):
                samp = sample(sample_size[0], sample_size[1], seed = seed)
                df['sampsize'] = samp
        
            if muprior is not None:
                if len(muprior) == 1:
                    mu = muprior[0]
                else:
                    mu = sample(muprior[0], muprior[1], distr = prior_distrib, seed = seed)
                df['mu'].iloc[i,] = mu

            if comprior is not None:
                if len(comprior) == 1:
                    com_size = comprior[0]
                else:
                    com_size = sample(comprior[0], comprior[1], distr = prior_distrib, typ = "int", seed = seed)
                df['comsize'].iloc[i,] = com_size
            
            if pastprior is not None:
                past_sizes = [sample(pastprior[0], pastprior[1], distr = prior_distrib, typ = "int", seed = seed) for _ in range(nepoch)]
                colnames = ['pastsize{}'.format(i) for i in range(1, nepoch+1)]
                for j in range(nepoch):
                    df[colnames[j]].iloc[i,] = past_sizes[j]

            if init_sizeprior is not None and init_ratesprior is not None and npop is not None:
                init_rates = [sample(init_ratesprior[0], init_ratesprior[1], distr = prior_distrib, typ = "float", seed = seed) for _ in range(npop)][0]
                init_sizes = [sample(init_sizeprior[0], init_sizeprior[1], distr = prior_distrib, typ = "float", seed = seed) for _ in range(npop)][0]
                colnames1 = ['initrate{}'.format(i) for i in range(1, npop+1)]
                colnames2 = ['initsize{}'.format(i) for i in range(1, npop+1)]
                for j in range(npop):
                    df[colnames1[j]].iloc[i,] = init_rates[j][0]
                    df[colnames2[j]].iloc[i,] = init_sizes[j][0]

            if withmigr and npop is not None:
                m = sample(0, 1, distr = prior_distrib, typ = "float", seed = seed)
                df['m'].iloc[i,] = m

            if massprior is not None:
                # probably better to fix split dates same as changetime
                split_dates = [params(massprior, 1, distrib = prior_distrib) for _ in range(nsplit)][0]
                colnames = ['splitdate{}'.format(i) for i in range(1, nsplit+1)]
                for j in range(nsplit):
                    df[colnames[j]].iloc[i,] = split_dates[j][0]
            # simulate phylogeny
            phylo = simulate(sample_size = samp,
                             com_size = com_size,
                             mu = mu,
                             mrca = lim_mrca,
                             npop = npop,
                             changetime = changetime,
                             m = m,
                             init_rates = init_rates,
                             init_sizes = init_sizes,
                             past_sizes = past_sizes,
                             split_dates = split_dates,
                             migrfrom = migrfrom,
                             migrto = migrto,
                             verbose = verbose, seed = seed)

            if sstype == 'SFS':
                ss.append(getAbund(phylo, samp))

            if savetrees:
                trees += phylo.write() + "\n"
        except: 
            df = df[:-1]
        else:
            i += 1
        finally: 
            safecount+=1

    if sstype == 'SFS':
        SFS = np.zeros((len(ss), max(len(x) for x in ss)))
        for k, j in enumerate(ss):
            SFS[k][0:len(j)] = j
    ssdf = pd.DataFrame(SFS)
    
    if savetrees:
        saved = df.to_string()
        if sstype == 'SFS':
            saved += "\n###\n" + ssdf.to_string()
        saved += "\n###\n" + trees

        f = open(saveto, "w")
        f.write(saved)
        f.close()
    return df, ssdf

def simulate(sample_size, com_size, mu, mrca = None, npop = 1, 
             m = 0, init_rates = None, init_sizes = None, 
             past_sizes = None, changetime = None, split_dates = None, 
             migrfrom = None, migrto = None, verbose = False, seed = None):
    """
    Simulate a phylogeny with msprime
    
    Parameters
    ----------
    sample_size : TYPE
        DESCRIPTION
    com_size : TYPE
        DESCRIPTION
    mu : TYPE
        DESCRIPTION
    mrca = None : TYPE
        DESCRIPTION
    npop = 1 : TYPE
        DESCRIPTION
    m = 0 : TYPE
        DESCRIPTION
    init_rates = None : TYPE
        DESCRIPTION
    init_sizes = None : TYPE
        DESCRIPTION
    past_sizes = None : TYPE
        DESCRIPTION
    changetime = None : TYPE
        DESCRIPTION
    split_dates = None : TYPE
        DESCRIPTION
    migrfrom = None : TYPE
        DESCRIPTION
    migrto = None : TYPE
        DESCRIPTION
    verbose = False : TYPE
        DESCRIPTION
    seed = None
    
    tree : TYPE
        DESCRIPTION.
    mu : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    Examples
    --------
    >>> t = simulate(10, 1e5, 0.03, seed = 42)
    >>> print(t)
    <BLANKLINE>
          /-1
       /-|
      |  |   /-8
      |   \-|
    --|      \-0
      |
      |   /-7
      |  |
       \-|      /-5
         |   /-|
          \-|   \-3
            |
             \-6
    >>> simulate(1e5, 1e5, 0.03, seed = 42)
    """

    # TODO : doc !!!
    # TODO : idiotproof
    # TODO : more examples
    # do dummy checks here --> try to make code stupid-proof
    if sample_size >= com_size:
        sys.exit("Sample size should not exceed community size")
    
    popchange = []
    massmigration = []
    popconfig = None
    migration = None

    # make past demographic changes between different time frames
    if past_sizes is not None and changetime is not None:
        if len(past_sizes) != len(changetime):
            sys.exit("There should be as many sizes as there are past epochs")
        popchange = demographic_events(changetime, past_sizes)

    # make island model
    if npop > 1:
        if init_sizes is None or init_rates is None:
            sys.exit("Initial population sizes and growth rates should be provided when there are more than one population (npop>1)")
        migration = islmodel.migration_matrix(npop, m)
        samples = np.ones(npop, dtype=int)*sample_size
        popconfig = islmodel.population_configurations(samples, init_sizes, init_rates)

        # possible mass migration between populations
        if split_dates is not None:
            # implement option later for limited mass dispersal
            M = 1
            massmigration = islmodel.mass_migrations(split_dates, migrfrom, migrto, M)

    demography = popchange + massmigration
    if len(demography) == 0:
        demography = None

    # if verbose should print the demography debugger - only for debugging purposes!!! 
    if verbose:
        dd = msprime.DemographyDebugger(Ne = com_size, 
                                        demographic_events= demography, 
                                        migration_matrix= migration, 
                                        population_configurations= popconfig)
        dd.print_history()
    
    if npop > 1:
         treeseq = msprime.simulate(Ne = com_size,
                                    random_seed= seed,
                                    population_configurations = popconfig,
                                    migration_matrix = migration,
                                    demographic_events = demography)
    else: 
        treeseq = msprime.simulate(sample_size= sample_size,
                                   Ne = com_size,
                                   random_seed= seed,
                                   demographic_events = demography)

    tree = treeseq.first()
    if verbose: print(tree.draw(format = 'unicode'))
    if mrca is not None:
        if tree.time(tree.root) > mrca : 
            raise Exception(f"Simulated MRCA ({tree.time(tree.root)}) predates"+
                             " fixed limit ({mrca})")
    #print(tree.draw(format="unicode"))
    node_labels = {u: str(u) for u in tree.nodes() if tree.is_sample(u)}
    tree = Tree(tree.newick(node_labels = node_labels))
    phylo = toPhylo.toPhylo(tree, mu, seed = seed)

    return phylo   


def sample(lower, upper, distr = "uniform", typ = "float", seed = None):
    random.seed(seed)
    if upper == lower :
        p = upper
    
    if distr == "uniform":
        if typ == "int":
            p = random.randrange(lower,upper)
        else:
            p = random.uniform(lower,upper)

    if distr == "log_unif":
        p =   loguniform(lower, upper).rvs()

    return p 


def params(lim, nsim, distrib = "uniform", typ = "float", seed = None):
    """
    Make a list of parameters to test from prior distribution limits.

    Parameters
    ----------
    lim: list
        the inferior and superior limits of the prior distribution.
    nsim: int
        number of replicates.
    distrib: str, optional
        the type of prior distribution. The default is "uniform".

    Returns
    -------
    params: list
        a list of nsim parameter values.

    Examples
    --------
    """
    # TODO : idiotproof
    # TODO : more examples
    if len(lim) == 1:
        p = np.repeat([lim[0]], nsim)
    else:
        np.random.seed(seed)
        if distrib == "uniform":
            if typ == "int":
                p = np.random.randint(lim[0], lim[1], size = nsim, dtype=np.int64)
            if typ == "float":
                p = np.random.uniform(lim[0], lim[1], size = nsim)
        if distrib == "logunif":
            p = np.exp(np.random.uniform(lim[0], lim[1], size = nsim))
    return p


def getAbund(tree, samp_size):
    """
    

    Parameters
    ----------
    tree : TYPE
        DESCRIPTION.

    Returns
    -------
    sfs : TYPE
        DESCRIPTION.

    Examples
    --------
    """
    # TODO : idiotproof
    # TODO : more examples
    sfs = list()
    abund = list()
    for leaf in tree.iter_leaves():
        try:
            inds = leaf.mergedInd.lstrip(" ")
            inds = list(inds.split(" "))
            abund.append(len(inds))
        except AttributeError:
            abund.append(1)
    sfs.extend(abund)
    # thin about catching error when phylogeny has only 1 sp
    if sum(sfs) != samp_size:
        raise Exception(f"Simulated phylogeny has only one species!")
    return sfs

if __name__ == "__main__":
        import doctest
        doctest.testmod()
