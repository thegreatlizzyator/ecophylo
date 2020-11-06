# -*- coding: utf-8 -*-
"""
 ALL ECOPHYLO


Created on Wed May 13 11:42:50 2020

@author: barthele
"""
import msprime
import numpy as np
import sys
from ete3 import Tree
import pandas as pd

from . import toPhylo

def dosimuls(nsim, sample_size, comprior, muprior, lim_mrca = None, sstype="SFS",
             npop=1, nepoch=1, withmigr=False, init_ratesprior=None,
             init_sizeprior=None, pastprior=None, maxtime=None,
             nsplit=None, massprior=None, migrfrom=None, migrto=None,
             verbose=False):

    # CHECKS HERE FOR IDIOT-PROOFING
    if maxtime is None and pastprior is not None:
        sys.exit("Time of past demographic change should be provided")
    
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
        colnames = ['pastsize{}'.format(i) for i in range(1, nepoch)]
        for i in range(nepoch-1):
            df[colnames[i]] = ""
        comments += (f'\n -{nepoch-1} past epoch(s) in which community size varies in {pastprior}')
    else:
        past_sizes = None

    # date of past demographic change
    if maxtime is not None and isinstance(maxtime, list):
        df['changetime'] = ""
    else:
        changetime = None

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
    if verbose:
        print(comments)
    
    # FILL PARAMETER DATAFRAME AND SIMULATE PHYLOGENIES
    safecount = 0
    i = 0
    ss = list()

    while i < nsim:
        if safecount > 10000:
            if i > 0:
                print(f"Many simulations have failed. Returning fewer simulations than {nsim}.")
                if sstype == 'SFS':
                    SFS = np.zeros((nsim, max(len(x) for x in ss)))
                    for i, j in enumerate(ss):
                        SFS[i][0:len(j)] = j
                ssdf = pd.DataFrame(SFS)
                return df, ssdf
            else:
                sys.exit("Too many simulations have failed")

        try:
            df = df.append(pd.Series(), ignore_index=True)
            # sample parameters from prior for simulation
            if isinstance(sample_size, list):
                samp = params(sample_size, nsim, typ = "int")
                df['sampsize'] = samp
        
            if muprior is not None:
                mu = params(muprior, 1)[0]
                df['mu'].iloc[i,] = mu

            if comprior is not None:
                com_size = params(comprior, 1, typ="int")[0]
                df['comsize'].iloc[i,] = com_size
            
            if pastprior is not None:
                past_sizes = [params(pastprior, 1, typ ="int") for _ in range(nepoch-1)]
                colnames = ['pastsize{}'.format(i) for i in range(1, nepoch)]
                for j in range(nepoch-1):
                    df[colnames[j]].iloc[i,] = past_sizes[j][0]

            if maxtime is not None and isinstance(maxtime, list):
                changetime = params(maxtime, 1, typ = "int")[0]
                df['changetime'].iloc[i,] = changetime

            if init_sizeprior is not None and init_ratesprior is not None and npop is not None:
                init_rates = [params(init_ratesprior, 1) for _ in range(npop)][0]
                init_sizes = [params(init_sizeprior, 1) for _ in range(npop)][0]
                colnames1 = ['initrate{}'.format(i) for i in range(1, npop+1)]
                colnames2 = ['initsize{}'.format(i) for i in range(1, npop+1)]
                for j in range(npop):
                    df[colnames1[j]].iloc[i,] = init_rates[j][0]
                    df[colnames2[j]].iloc[i,] = init_sizes[j][0]

            if withmigr and npop is not None:
                m = params([0, 1], 1)[0]
                df['m'].iloc[i,] = m

            if massprior is not None:
                split_dates = [params(massprior, 1) for _ in range(nsplit)][0]
                colnames = ['splitdate{}'.format(i) for i in range(1, nsplit+1)]
                for j in range(nsplit):
                    df[colnames[j]].iloc[i,] = split_dates[j][0]

            # simulate phylogeny
            phylo = simulate(sample_size = samp,
                             com_size = com_size,
                             mu = mu,
                             mrca = lim_mrca,
                             npop = npop,
                             nepoch = nepoch,
                             m = m,
                             init_rates = init_rates,
                             init_sizes = init_sizes,
                             past_sizes = past_sizes,
                             maxtime = changetime,
                             split_dates = split_dates,
                             migrfrom = migrfrom,
                             migrto = migrto,
                             verbose = False)

            if sstype == 'SFS':
                ss.append(getSFS(phylo))
        except: 
            df = df[:-1]
        else:
            i = i+1
        finally: 
            safecount+=1
    
    if sstype == 'SFS':
        SFS = np.zeros((nsim, max(len(x) for x in ss)))
        for i, j in enumerate(ss):
            SFS[i][0:len(j)] = j
    ssdf = pd.DataFrame(SFS)

    return df, ssdf

def simulate(sample_size, com_size, mu, mrca = None, npop = 1, nepoch = 1, m = 0, init_rates = None, init_sizes = None, past_sizes = None, maxtime = None, split_dates = None, migrfrom = None, migrto = None, seed = 1, verbose = False):

    # do dummy checks here --> try to make code stupid-proof
    if sample_size >= com_size:
        sys.exit("Sample size should not exceed community size")
    
    popchange = []
    massmigration = []
    popconfig = None
    migration = None

    # make past demographic changes between different time frames
    if nepoch > 1:
        if len(past_sizes) != nepoch - 1 :
            sys.exit("There should be as many sizes as there are epochs")
        if nepoch > 2: 
            times = timeframes(nepoch-1, maxtime, 0.05)
        else:
            times = [maxtime]
        popchange = demographic_events(times, past_sizes)

    # make island model
    if npop > 1:
        if init_sizes is None or init_rates is None:
            sys.exit("Initial population sizes and growth rates should be provided when there are more than one population (npop>1)")
        migration = migration_matrix(npop, m)
        samples = np.ones(npop, dtype=int)*sample_size
        popconfig = population_configurations(samples, init_sizes, init_rates)

        # possible mass migration between populations
        if split_dates is not None:
            # implement option later for limited mass dispersal
            M = 1
            massmigration = mass_migrations(split_dates, migrfrom, migrto, M)

    demography = popchange + massmigration

    if len(demography) == 0:
        demography = None

    # if verbose should print the demography debugger - only for debugging purposes!!! 
    if verbose:
        dd = msprime.DemographyDebugger(Ne = com_size, demographic_events= demography, migration_matrix= migration, population_configurations= popconfig)
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
    if mrca is not None:
        if tree.time(tree.root) > mrca : 
            raise Exception(f"Simulated MRCA ({tree.time(tree.root)}) predates fixed limit ({mrca})")
    #print(tree.draw(format="unicode"))
    node_labels = {u: str(u) for u in tree.nodes() if tree.is_sample(u)}
    tree = Tree(tree.newick(node_labels = node_labels))
    phylo = toPhylo.toPhylo(tree, mu)

    return phylo   

def params(lim, nsim, distrib = "uniform", typ = "float"):
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

    """
    if len(lim) == 1:
        p = np.repeat([lim[0]], nsim)
    else:
        #only supports uniform prior distributions at the moment
        if distrib == "uniform":
            if typ == "int":
                p = np.random.randint(lim[0], lim[1], size = nsim, dtype=np.int64)
            else:
                p = np.random.uniform(lim[0], lim[1], size = nsim)
    return p


def getSFS(tree):
    """
    

    Parameters
    ----------
    tree : TYPE
        DESCRIPTION.

    Returns
    -------
    sfs : TYPE
        DESCRIPTION.

    """
    sfs = list()
    abund = list()
    for leaf in tree.iter_leaves():
        try:
            inds = leaf.mergedInd.lstrip(" ")
            inds = list(inds.split(" "))
            abund.append(len(inds))
        except AttributeError:
            abund.append(1)
    abund.sort()
    sfs.extend(abund)
    return sfs

