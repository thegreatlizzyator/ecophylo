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
    getDeme

"""
    # TODO : doc !!!
import msprime
import random
import numpy as np
import sys
import warnings
from ete3 import Tree
import pandas as pd
#from loguniform import LogUniform
from scipy.stats import loguniform

from ecophylo import pastdemo
from ecophylo import islmodel
from ecophylo import phylogen

def dosimuls(nsim, sample_size, comprior, muprior, lim_mrca = None, sstype="SFS",
             prior_distrib = "uniform", npop=1, withmigr=False, init_ratesprior=None,
             init_sizeprior=None, pastprior=None, changetime = None,
             nsplit=None, massprior=None, migrfrom=None, migrto=None,
             verbose=False, savetrees= False, saveto = "", seed = None):
    """
  
    nsim : int
    comprior : list of int 
        lenght = 2
    muprior : list of float
        length = 2
    sstype : str
    prior_distrib
    withmigr : bool
    inits_ratesprior : TYPE
    nsplit : TYPE
    massprior : TYPE
    pastprior : TYPE
    savetrees : TYPE
    saveto : TYPE
  
    # TODO : find simulate parameters
    # TODO : rename lim_mrca with mrca
    
    Examples
    --------
    """

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
    sample_size : int
        number of individual in the community
        Sample size should not exceed community size
        # TODO : renommer init_size
    com_size : int
        taille de la meta-communauté sensu hubbel 2001
    mu : float
        mutation rate, must be comprised between between 0 and 1.
    mrca = None : int
        number of generation
        # TODO : test if it is not a float
        DESCRIPTION
    npop = 1 : int
        > 0
        DESCRIPTION
    m = 0 : float
        DESCRIPTION
        # TODO : rename migr: float taux de migration
        # TODO : test between 0 and 1
        overall symetric migration rate. Default is 0. 
    init_rates = None : list of float$
        list of length = npop
        DESCRIPTION
        # TODO : rename with ilsmodel rates : list of float
        # TODO : test rates limites (real, -1:1 or -inf:inf)
        The initial population growth rates
    init_sizes = None : list of int
        DESCRIPTION
        # TODO : same as in islmodels init_sizes : list of int
        positive values
        The initial population sizes
    past_sizes = None : list of int
        Population past_sizes at the different time periods. Must be > 0.
        Must be of same length as npop and same length as changetime.
    changetime = None : list of int
        When the demographic changes have occured. Must be >= 0.
        Must be of same length as past_sizes
    split_dates = None : list of int
        # TODO : check if it work with list length == 1
        DESCRIPTION
    migrfrom = None : TYPE
        DESCRIPTION
        # TODO : rename accordingly to sources in islmodel
    migrto = None : TYPE
        DESCRIPTION
        # TODO : rename accordingly to destinations in islmodel
    verbose = False : bool
        DESCRIPTION
    seed = None : int
        An integer used to set the seed in all the random events in msprime
        and in the mutation on the phylogeny.

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
    Traceback (most recent call last):
      File "/usr/lib/python3.6/doctest.py", line 1330, in __run
        compileflags, 1), test.globs)
      File "<doctest ecophylo.dosimulate.simulate[2]>", line 1, in <module>
        simulate(1e5, 1e5, 0.03, seed = 42)
      File "/home/maxime/Bureau/BEE/ecophylo/ecophylo/dosimulate.py", line 289, in simulate
        sys.exit("Sample size should not exceed community size")
    SystemExit: Sample size should not exceed community size
    """
    
    if not isinstance(seed, (int,float)):
        sys.exit('seed must be an integer')
    if isinstance(seed, float):
      seed = int(seed)

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
        popchange = pastdemo.demographic_events(changetime, past_sizes)

    # make island model
    if npop > 1:
        if init_sizes is None or init_rates is None:
            sys.exit("Initial population sizes and growth rates should be provided when there are more than one population (npop>1)")
        
        # subpops =  npop
        
        migration = islmodel.migration_matrix(npop, m)
        samples = np.ones(npop, dtype=int)*sample_size
        # TODO : allow differential sampling in pop (provide sample list same length as npop)
        popconfig = islmodel.population_configurations(init_sizes, samples, init_rates)

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
        dd.print_history(output=sys.stderr) # this will get deprecated with msprime 1.0
    
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
    # if verbose: print(tree.draw(format = 'unicode'))
    if mrca is not None:
        if tree.time(tree.root) > mrca : 
            raise Exception(f"Simulated MRCA ({tree.time(tree.root)}) predates"+
                             " fixed limit ({mrca})")
    #print(tree.draw(format="unicode"))
    node_labels = {u: str(u) for u in tree.nodes() if tree.is_sample(u)}
    tree = Tree(tree.newick(node_labels = node_labels))
    phylo = phylogen.toPhylo(tree, mu, seed = seed)

    return phylo


def simulate_dolly(sample_size, com_size, mu, init_rates = None, 
                   changetime = None, mrca = None, 
                   migr = 1, migr_time = None, verbose = False, seed = None):
    """
    This function implements the simulation algorithm described in Barthelemy
    et al. 2021 in which (i) the shared co-ancestry of present individuals is
    simulated backward in time using coalescent theory (ii) speciation events
    are sprinkled over the simulated genealogy conditionally to its topology
    and branch lengths and (iii) the phylogenetic relationships amongst
    individuals and their abundances are finally obtained by merging
    paraphyletic clades into single species. Coalescent reconstruction of the
    genealogy of individuals can be simulated to represent past demographic
    fluctuations due to varying habitat availability, or include multiple demes
    linked by migration events and/or vicariance.

    Phylogenies are returned in Newick format given the desired parameter
    combinations accounting for the demographic history of Jm

    Parameters
    ----------
    sample_size : int or list of int
        number of sampled individuals in an assemblage for which the shared
        co-ancestry should be reconstructed. If multiple demes are to be
        simulated, this should be a list of sample sizes for each deme. These
        should not exceed the present or past assemblage sizes.
        # TODO : rename to samples
    com_size : int or nested list of ints
        the size of Jm for each deme at each given period. Should be a nested
        list containing for each deme, a list of past Jm sizes in which the
        first element is the current size of the assemblage and the nth element
        is the size of Jm at epoch n 
        # TODO: rename to Jm ?
    mu : float
        the point mutation rate must be comprised between between 0 and 1.
    init_rates: float or nested list of floats
        the growth rates for each deme at each given period. Should be a nested
        list containing for each deme, a list of past growth rates in which the
        first element is the current growth rate and the nth element is the
        growth rate at epoch n. If no growth rates are given, then changes in
        deme sizes occur instantenously following sizes provided in com_size at
        the different times given in changetime
        # TODO: rename to gr_rates ?
    changetime: list of int or nested list of int
        the times (in generation before present) at which either growth rates or
        the size of the assemblages Jm have changed. If multiple demes are to be
        simulated, should be a nested list containing for each deme, a list of
        times at which changes occured in which the first element is 0.
        # TODO: rename epochs? times?
    mrca = None : int
        # TODO : document this when it is implemented
    migr = 1 : int, float or list of int,float or nested lists of int,float
        the migration rates between pairs of demes at each given period. Can be
        an int or float comprised between 0 and 1, in which case constant
        symmetric migration is assumed between all demes for all epochs.

        If migration rate are to change then migr should be a list of ints or
        floats comprised between 0 and 1 containing the different symmetric
        migration rates at each given time period in which the first element is
        the current symmetric migration rate and the nth element is the migration
        rate at epoch n.

        For non-symmetric migration rates, migr should be a list of migration
        matrices M of size dxd where d is the number of demes. Migr should then
        contain as many matrices M as there are time periods in migr_time where
        M[j,k] is the rate at which individuals move from deme j to deme k in
        the coalescent process, backwards in time. Individuals that move from
        deme j to k backwards in time actually correspond to individuals
        migrating from deme k to j forwards in time.
    migr_time = None: list of ints
        the times (in generation before present) at which migration rates have
        changed in which the first element is 0
    verbose = False : bool
        whether or not to print a summary of the demographic history and the
        resulting genealogy to be passed to a phylogeny
    seed = None : int
        set seed for entire simulation
    
    Examples
    --------
    >>> t = simulate_dolly(sample_size = [10], com_size = [[1e5]], mu = 0.03, seed = 42)
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
    >>> t = simulate_dolly(sample_size = [5, 5], com_size = [[1e5], [1e5]], mu = 0.03, migr = 1, seed = 42)
    >>> print(t)
    <BLANKLINE>
          /-7
         |
       /-|      /-4
      |  |   /-|
      |   \-|   \-8
    --|     |
      |      \-0
      |
      |   /-3
       \-|
          \-1
    >>> t = simulate_dolly(sample_size = [5], com_size = [[1e3]], mu = 0.03, migr = 1, seed = 42)
    >>> print(t)
    <BLANKLINE>
          /-4
       /-|
    --|   \-0
      |
       \-2
    >>> t = simulate_dolly(sample_size = [5], com_size = [[1e3, 2e3]], changetime = [[0, 50]], mu = 0.03, migr = 1, seed = 42)
    >>> print(t)
    <BLANKLINE>
          /-3
       /-|
    --|   \-0
      |
       \-2
    >>> t = simulate_dolly(sample_size = [5, 5], com_size = [[1e5], [1e5]], mu = 0.03, migr = 1, seed = 42)
    >>> print(t)
    <BLANKLINE>
          /-7
         |
       /-|      /-4
      |  |   /-|
      |   \-|   \-8
    --|     |
      |      \-0
      |
      |   /-3
       \-|
          \-1
    >>> t = simulate_dolly(sample_size = [5, 5], com_size = [[1e3, 2e3], [1e3, 5e2]], changetime = [[0, 50],[0, 30]], mu = 0.03, migr = 1, seed = 42)
    >>> print(t)
    <BLANKLINE>
          /-2
       /-|
      |   \-0
      |
    --|         /-6
      |      /-|
      |   /-|   \-3
      |  |  |
       \-|   \-1
         |
          \-7
    >>> t = simulate_dolly(sample_size = [5, 5], com_size = [[1e3, 2e3], [1e3, 5e2]], changetime = [[0, 500],[0, 300]], mu = 0.03, migr = [0, 1], migr_time = [0, 200], seed = 42)
    >>> print(t)
    <BLANKLINE>
          /-3
       /-|
      |  |   /-1
      |   \-|
    --|      \-0
      |
      |      /-6
      |   /-|
       \-|   \-2
         |
          \-5
    """         
    # # parameters that will be used later when mass migration will be coded
    # split_dates = None # won't be used
    # migrfrom = None # won't be used
    # migrto = None # won't be used
    
    # Idiotproof
    # check sample_size
    if not isinstance(sample_size, list):
        sample_size = [sample_size]
    isint_samp = [isinstance(s, int) for s in sample_size]
    if not all(isint_samp):
        sys.exit("sample_size should all be ints")
    ispos_samp = [s>0 for s in sample_size]
    if not all(ispos_samp):
        sys.exit("sample_size should all be positive")
    # compute number of populations
    npop = len(sample_size)

        # check changetime
    if changetime is not None:
        if not isinstance(changetime, list):
            if isinstance(changetime, (int,float)) : 
                if changetime >= 0 :
                    changetime = [[changetime]]
                else :
                    sys.exit('changetime must be positive values')
            else :
                sys.exit('changetime must be int, list of int or'+
                ' nested list of int')
        else :
            for x in changetime:
                if isinstance(x, list):
                    if not all(isinstance(y, (float, int)) for y in x) : 
                        sys.exit('changetime must be int, list of int or'+
                                 ' nested list of int')
                    if any(y < 0 for y in x[1:]) : 
                        sys.exit('changetime must be positive values')
                    if x[0] != 0:
                        sys.exit('first element of changetime for a Deme'+
                        ' must be equal to 0')
                    if len(set(x)) != len(x) :
                        sys.exit('Duplicated times in changetime for a Deme' +
                                 ' are not possible')
                else :
                    if len(set(changetime)) != len(changetime) :
                        sys.exit('Duplicated times in changetime are not possible')
                    if not isinstance(x, (float, int)):
                        sys.exit('changetime must be int, list of int or'+
                                 ' nested list of int')
                    if x < 0 :
                        sys.exit('changetime must be positive values')
                    if changetime[0] != 0:
                        sys.exit('first element of changetime'+
                        ' must be equal to 0')
            if not isinstance(x, list):
                changetime = [changetime]
        if len(changetime) != npop :
            sys.exit("there should be as many past sizes as there " + 
            "are epochs in changetime")
    else :
        changetime = [[0]] * npop

    # check com_size 
    if changetime is not None  :
        isint_com = True
        sampl_com = True
        if not isinstance(com_size, list) :
            if isinstance(com_size, (int,float)) and com_size > 0 : 
                if com_size < sample_size[0] : sampl_com = False
                com_size = [[int(com_size)]] * npop
            else :
                isint_com = False
        else :
            for i in range(len(com_size)):
                if isinstance(com_size[i], list):
                    if not all(isinstance(y, (float, int)) for y in com_size[i]) :
                        isint_com = False
                    else :
                        com_size[i] = [int(x) if isinstance(x, float) else x for x in com_size[i]]
                    if len(com_size) != npop :
                        sys.exit("there should be as many elements in"+
                                 " com_size as there are demes")
                    if len(com_size[i]) !=  len(changetime[i]) :
                        sys.exit("there should be as many past "+
                        "com_size as there are epochs in changetime")
                    if isint_com and any([s <= 0 for s in com_size[i]]):
                        sys.exit("all past sizes should be strictly positive")
                    if isint_com and any([x < sample_size[i] for x in com_size[i]]):
                        sampl_com = False 
                else :
                    if len(com_size) != npop :
                        sys.exit("there should be as many elements in"+
                                 " com_size as there are demes")
                    if not isinstance(com_size[i], (float, int)):
                        isint_com = False 
                    if isint_com and com_size[i] <= 0 :
                        sys.exit("all past sizes should be strictly positive")
                    if isint_com and com_size[i] < sample_size[i] :
                        sampl_com = False 
                    com_size[i] = [com_size[i]]
        if not isint_com:
            sys.exit("community sizes should be strictly positive int")
        if not sampl_com:
            sys.exit('com_size must be superior to samples')
        
    # if isinstance(com_size, float) :
    #     com_size = [[int(com_size)]] * npop
    # check mu
    if not isinstance(mu, (int,float)) or mu < 0 or mu > 1 :
        sys.exit('mu must be a float between 0 and 1')
    # check init_rates
    if init_rates is not None and changetime is not None:
        isint_rates = True
        if not isinstance(init_rates, list):
            if isinstance(init_rates, (int,float)) : 
                init_rates = [[init_rates]] * npop
            else :
                isint_rates = False
        else :
            for i in range(len(init_rates)):
                if isinstance(init_rates[i], list):
                    if not all(isinstance(y, (float, int)) for y in init_rates[i]) :
                        isint_rates = False
                    if len(init_rates[i]) !=  len(changetime[i]) :
                        sys.exit("there should be as many past growth "+
                        "init_rates as there are epochs in changetime")
                else :
                    if len(init_rates) != npop :
                        sys.exit("there should be as many elements in"+
                        " init_rates as there are demes")
                    if not isinstance(init_rates[i], (float, int)):
                        isint_rates = False 
                    init_rates[i] <- [init_rates[i]]
        if not isint_rates:
            sys.exit('init_rates must be float, list of float or'+
                                 ' nested list of float')
        if len(init_rates) != npop :
            sys.exit("there should be as many past sizes as there " + 
            "are epochs in init_rates")
    else :
        init_rates = [[0]] * npop

    # check mrca
    if mrca is not None :
        print('prout') # TODO : do this
    # check migr
    if migr is not None :
        if npop == 1 :
            # warnings.warn("no migration matrix is needed for a single deme")
            migr = None
    if migr is not None :
        m = np.array(migr)
        dim = m.shape

        if np.sum(m) == 0 :
            sys.exit("migration matrices cannot all be empty")

        if len(dim) == 0:
            migr = [migr]
            m = np.array(migr)
            dim = m.shape
        if len(dim) > 1 :
            if dim[1] != dim[2] or dim[1] != npop or dim[2] != npop:
                sys.exit("custom migration matrices should be of size ndeme" +
                " x ndeme")
            for mat in migr:
                isint_migr = [isinstance(i, (float,int)) for i in mat]
                if not all(isint_migr) :
                    sys.exit("found custom migration matrix that is not made" +
                    " of ints or floats")
                ispos_migr = [i>=0 and i<= 1 for i in mat]
                if not all(ispos_migr):
                    sys.exit("found custom migration matrix with negative" +
                    " migration rates or greater than 1")
        else :
            isint_migr = [isinstance(i, (float,int)) for i in migr]
            if not all(isint_migr):
                sys.exit("migration rates should be either ints or floats")
            ispos_migr = [i>=0 and i<= 1 for i in migr]
            if not all(ispos_migr): 
                sys.exit("migration rate should be positive (or zero) and" + 
                " not exceed 1")
    
    # check migr_time
    if migr is not None and len(migr) > 1 :
        if len(migr) != len(migr_time):
            sys.exit("there should be as many migration rates or matrices as" +
            " there are times in migr_time")
    # check verbose
    if not isinstance(verbose, bool):
        sys.exit('verbose must be a boolean') 
    # check seed
    if seed is not None and not isinstance(seed, (int,float)):
        sys.exit('seed must be an integer')
    if seed is not None and isinstance(seed, float):
        seed = int(seed)
  
    # compute number of populations
    npop = len(sample_size)

    # initialise demography object for msprime.sim_ancestry()
    demography = msprime.Demography()
    
    # Build samples
    samples = {"pop_"+str(i):sample_size[i] for i in range(len(sample_size))}
    pop_ids = ["pop_" + str(p) for p in range(npop)]

    # initialize population configurations
    for pop in range(npop):
        demography.add_population(
            initial_size= com_size[pop][0],       
            growth_rate=init_rates[pop][0])
        
        # if population sizes have fluctuated in the past:
        if len(changetime[pop]) > 1 and len(com_size[pop]) > 1:
            for i in range(len(changetime[pop][1:])):
                demography.add_population_parameters_change(
                    time = changetime[pop][i+1] , 
                    initial_size=com_size[pop][i+1], 
                    population= pop_ids[pop])
        
        # if population growth rates have fluctuated in the past:
        if len(changetime[pop]) > 1 and len(init_rates[pop]) > 1:
            for i in range(len(changetime[pop][1:])):
                demography.add_population_parameters_change(
                    time = changetime[pop][i+1] , 
                    growth_rate=init_rates[pop][i+1], 
                    population= pop_ids[pop])

    ## MIGRATION
    if migr is not None :
        if len(dim) == 1 :
            # symmetric migration matrix
            demography.set_symmetric_migration_rate(
                populations = range(npop), rate = migr[0])

            # if symmatric migration rate has changed in the past:
            if len(migr) > 1:
                for m in range(len(migr[1:])):
                    demography.add_migration_rate_change(
                        time = migr_time[m+1], rate = migr[m+1])

        if len(dim) > 2 :
            # custom migration matrix
            for row in range(npop):
                for col in range(npop):
                    if migr[0][row][col] == 0 :
                        continue
                    demography.set_migration_rate(
                        source = pop_ids[row], dest = pop_ids[col], 
                        rate = migr[0][row][col])

            # if migration matrix has changed in the past
            if len(migr)>1:
                for t in range(len(migr_time)):
                    for row in range(npop):
                        for col in range(npop):
                            if migr[t][row][col] == 0 :
                                continue
                            demography.add_migration_rate_change(
                                time = migr_time[t], rate = migr[t][row][col], 
                                source=pop_ids[row], dest=pop_ids[col])

    # sort events chronologically
    demography.sort_events()

    # if verbose should print the demography debugger - only for debugging purposes!!! 
    if verbose:
        print(demography.debug())
        
    treeseq = msprime.sim_ancestry(samples = samples, 
          demography=demography, random_seed=seed, ploidy = 1)
    
    # Work on the result tree
    tree = treeseq.first()
    if verbose: print(tree.draw(format = 'ascii'))
    if mrca is not None:
        if tree.time(tree.root) > mrca : 
            raise Exception(f"Simulated MRCA ({tree.time(tree.root)}) predates"+
                             " fixed limit ({mrca})")
    #print(tree.draw(format="unicode"))
    node_labels = {u: str(u) for u in tree.nodes() if tree.is_sample(u)}
    tree = Tree(tree.newick(node_labels = node_labels))
    phylo = phylogen.toPhylo(tree, mu, seed = seed)

    return phylo  


def sample(lower, upper, distr = "uniform", typ = "float", seed = None):
    """
    

    Parameters
    ----------
    lower : int or float
        the inferior and superior limits of the prior distribution.
        length == 2
    upper : int or float
        number of replicates.
    distr : str
        distribution law choice are uniform and log_unif
    typ : str
        float or int, typ of lim elements.
    seed : int
        None by default, set the seed for mutation random events.

    Returns
    -------
    params: 
      # TODO : test what type is returned

    Examples
    --------
    """
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
    lim: list of int or float
        the inferior and superior limits of the prior distribution.
        length == 2
    nsim: int
        number of replicates.
    distrib: str, optional
        the type of prior distribution. The default is "uniform".
        choice are uniform, logunif
    typ : str
        float or int, typ of lim elements.
    seed : int
        None by default, set the seed for mutation random events.

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
            # TODO : add case for int input
            p = np.exp(np.random.uniform(lim[0], lim[1], size = nsim))
    return p


def getAbund(tree, sample_size = None):
    """
    
    Parameters
    ----------
    tree : (ete3 class)
        Phylogeny with attributes on leafs. This attributes is a character 
        string containing all names of the species individual (mean to use 
        topPhylo result). Names is formated like this :
          " name1 name2 name3"
    sample_size : int
        number of individual in the community
        # TODO : rename sample_size
        # TODO : set this check as optionnal

    Returns
    -------
    sfs : TYPE
        DESCRIPTION.

    Examples
    --------
    >>> from ete3 import Tree
    >>> tree = Tree('(((A:5,(B:3, C:3))1:2,(D:2, E:2)1:5)1:2, (F:3, G:3)1:6);')
    >>> import ecophylo as eco
    >>> phylo = eco.toPhylo(tree, 0.5, seed = 42)
    >>> print(phylo)
    <BLANKLINE>
          /-A
       /-|
    --|   \-D
      |
       \-F
    >>> getAbund(phylo, 7)
    [3, 2, 2]
    """
    # Idiot proof
    if tree.__class__.__name__ != 'TreeNode' :
        sys.exit('tree must have a class TreeNode')
    if sample_size != None:
      if not isinstance(sample_size, int):
          sys.exit('sample_size must be an integer')

    sfs = list()
    abund = list()
    for leaf in tree.iter_leaves():
        try:
            inds = leaf.mergedInd.lstrip(" ") # remove 1st space
            inds = list(inds.split(" ")) # strip on spaces
            abund.append(len(inds)) # get length 
        except AttributeError:
            abund.append(1)
    sfs.extend(abund)
    # think about catching error when phylogeny has only 1 sp
    
    if sample_size != None and sum(sfs) != sample_size:
        raise Exception(f"Simulated phylogeny has only one species!")
        # TODO : modify error with a better check here
    return sfs

def getDeme(tree, div = False):
    """
    
    Parameters
    ----------
    tree : (ete3 class)
        Phylogeny with attributes on leafs. This attributes is a character 
        string containing all names of the species individual (mean to use 
        topPhylo result). Names is formated like this :
          " name1 name2 name3"
    div : bool
        Option to simplify the matrix to a simple list of Deme species diversity.

    Returns
    -------
    indiv : nest list of int
        site/species matrix of the tree

    Examples
    --------
    >>> from ete3 import Tree
    >>> tree = Tree('(((A_0:5,(B_0:3, C_1:3))1:2,(D_1:2, E_1:2)1:5)1:2, (F_2:3, G_2:3)1:6);')
    >>> import ecophylo as eco
    >>> phylo = eco.toPhylo(tree, 0.5, seed = 42)
    >>> print(phylo)
    <BLANKLINE>
          /-A
       /-|
    --|   \-D
      |
       \-F
    >>> getDeme(phylo)
    [[2, 1, 0], [0, 2, 0], [0, 0, 2]]
    >>> getDeme(phylo, div = True)
    [1, 2, 1]
    """
    # Idiot proof
    if tree.__class__.__name__ != 'TreeNode' :
        sys.exit('tree must have a class TreeNode')
    
    indiv = list()
    for leaf in tree.iter_leaves():
        try:
            indiv.append(leaf.popInd)
        except AttribueError :
            indiv.append(1)
    if div:
        indiv = np.array(indiv)
        indiv = [sum(indiv[:,i] > 0) for i in range(indiv.shape[1])]
    return indiv

if __name__ == "__main__":
        import doctest
        doctest.testmod()
        #simulate_dolly(sample_size = [5, 5], com_size = [[500], [500]], mu = 0.05, migr = 1, verbose = True, seed = 42)
        # simulate_dolly(sample_size = [10], com_size = [[500, 1000]], mu = 0.05, changetime= [[0,100]], seed = 42, verbose = True )
        # t = simulate_dolly(sample_size = [5, 5], com_size = [[1e3, 2e3], [1e3, 5e2]], changetime = [[0, 50],[0, 30]], mu = 0.03, migr = 2, verbose = True)
        # print(t)

        ## SINGLE POP
        # stat discret
        # t = simulate_dolly(sample_size = [5], com_size = [[1e3]], mu = 0.03, migr = 2, seed = 42, verbose = True)
        # stat continue 
        # t = simulate_dolly(sample_size = [5], com_size = [[1e3]], stable_pop = False, mu = 0.03, migr = 2, seed = 42, verbose = True)
        
        # fluct discret
        # t = simulate_dolly(sample_size = [5], com_size = [[1e3, 2e3]], changetime = [[0, 50]], mu = 0.03, migr = 2, seed = 42, verbose = True)
        # fluct continue # TODO : marche avec certaines valeurs cheloues
        # t = simulate_dolly(sample_size = [5], com_size = [[10, 2e5]], changetime = [[0, 800]], stable_pop = False, mu = 0.03, migr = 2, seed = 42, verbose = True)
        
        ## MULT POP

        # stat discret
        # t = simulate_dolly(sample_size = [5, 5], com_size = [[1e3], [1e3]], mu = 0.03, migr = 1, seed = 42, verbose = True)
        # stat continue # 
        # t = simulate_dolly(sample_size = [5, 5], com_size = [[1e3], [1e3]], stable_pop = False, mu = 0.03, migr = 2, seed = 42, verbose = True)

        # fluct discret
        # t = simulate_dolly(sample_size = [5, 5], com_size = [[1e3, 2e3], [1e3, 5e2]], changetime = [[0, 50],[0, 30]], mu = 0.03, migr = 2, seed = 42, verbose = True)
        # fluct continue # TODO : marche bizarrement
        # t = simulate_dolly(sample_size = [5, 5], com_size = [[1e3, 2e3], [1e3, 5e2]], changetime = [[0, 50],[0, 30]], stable_pop = False, mu = 0.03, migr = 2, seed = 42, verbose = True)

        # t = simulate_dolly(sample_size = [5, 5], com_size = [[1e3, 2e3], [1e3, 5e2]], changetime = [[0, 500],[0, 300]], mu = 0.03, migr = [0, 1], migr_time = [0, 200], seed = 42, verbose = True)
        # print(t)


