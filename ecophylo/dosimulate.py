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
import warnings
from ete3 import Tree
import pandas as pd
from scipy.stats import loguniform
import collections

from ecophylo import pastdemo
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
        raise ValueError("Ecophylo only supports uniform or log-uniform prior distributions for the moment")
    
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
            raise ValueError("Too many simulations have failed")

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


def simulate(samples, com_size, mu, init_rates = None, 
                   changetime = None, mrca = None, 
                   migr = 1, migr_time = None, vic_events = None,
                   verbose = False, seed = None):
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
    samples : int or list of int
        number of sampled individuals in an assemblage for which the shared
        co-ancestry should be reconstructed. If multiple demes are to be
        simulated, this should be a list of sample sizes for each deme. These
        should not exceed the present or past assemblage sizes.
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
    vic_events = nested list of ints
        a nested list detailing the different split events that should be 
        included in the simulation. Each element of vic_events should be a list
        specifying, in order: the date (in generations before present) at which
        the split occured, the demes resulting from the split (as a list of ints)
        and finally the ancestral deme number. For instance, if deme 1 splits 
        into deme 0 and deme 1 then vic_events =  [[time01, [0,1], 1]]
        Note that time01 should appear in changetime. Also, user should specify
        in com_size (at the correct position i.e to the size of the ancestral 
        deme at time01) the size of the ancestral deme when the split occurs. 
    verbose = False : bool
        whether or not to print a summary of the demographic history and the
        resulting genealogy to be passed to a phylogeny
    seed = None : int
        set seed for entire simulation
    
    Examples
    --------
    >>> t = simulate(samples = [10], com_size = [[1e5]], mu = 0.03, seed = 42)
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
    >>> t = simulate(samples = [5, 5], com_size = [[1e5], [1e5]], 
    ... mu = 0.03, migr = 1, seed = 42)
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
    """         
    # # parameters that will be used later when mass migration will be coded
    # split_dates = None # won't be used
    # migrfrom = None # won't be used
    # migrto = None # won't be used
    
    # Idiotproof
    # check samples
    if not isinstance(samples, list):
        samples = [samples]
    isint_samp = [isinstance(s, int) for s in samples]
    if not all(isint_samp):
        raise ValueError("samples should all be ints")
    ispos_samp = [s>0 for s in samples]
    if not all(ispos_samp):
        raise ValueError("samples should all be positive")
    # compute number of populations
    npop = len(samples)

    # check changetime
    if changetime is not None:
        if not isinstance(changetime, list):
            if isinstance(changetime, (int,float)) : 
                if changetime >= 0 :
                    changetime = [[changetime]]
                else :
                    raise ValueError('changetime must be positive values')
            else :
                raise ValueError('changetime must be int, list of int or'+
                ' nested list of int')
        else :
            for x in changetime:
                if isinstance(x, list):
                    if not all(isinstance(y, (float, int)) for y in x) : 
                        raise ValueError('changetime must be int, list of int or'+
                                 ' nested list of int')
                    if any(y < 0 for y in x[1:]) : 
                        raise ValueError('changetime must be positive values')
                    if x[0] != 0:
                        raise ValueError('first element of changetime for a Deme'+
                        ' must be equal to 0')
                    if len(set(x)) != len(x) :
                        raise ValueError('Duplicated times in changetime for a Deme' +
                                 ' are not possible')
                else :
                    if len(set(changetime)) != len(changetime) :
                        raise ValueError('Duplicated times in changetime are not possible')
                    if not isinstance(x, (float, int)):
                        raise ValueError('changetime must be int, list of int or'+
                                 ' nested list of int')
                    if x < 0 :
                        raise ValueError('changetime must be positive values')
                    if changetime[0] != 0:
                        raise ValueError('first element of changetime'+
                        ' must be equal to 0')
            if not isinstance(x, list):
                changetime = [changetime]
        if len(changetime) != npop :
            raise ValueError("there should be as many past sizes as there " + 
            "are epochs in changetime")
    else :
        changetime = [[0]] * npop

    # check com_size 
    if changetime is not None  :
        isint_com = True
        sampl_com = True
        if not isinstance(com_size, list) :
            if isinstance(com_size, (int,float)) and com_size > 0 : 
                if com_size < samples[0] : sampl_com = False
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
                        raise ValueError("there should be as many elements in"+
                                 " com_size as there are demes")
                    if len(com_size[i]) !=  len(changetime[i]) :
                        raise ValueError("there should be as many past "+
                        "com_size as there are epochs in changetime")
                    if isint_com and any([s <= 0 for s in com_size[i]]):
                        raise ValueError("all past sizes should be strictly positive")
                    if isint_com and any([x < samples[i] for x in com_size[i]]):
                        sampl_com = False 
                else :
                    if len(com_size) != npop :
                        raise ValueError("there should be as many elements in"+
                                 " com_size as there are demes")
                    if not isinstance(com_size[i], (float, int)):
                        isint_com = False 
                    if isint_com and com_size[i] <= 0 :
                        raise ValueError("all past sizes should be strictly positive")
                    if isint_com and com_size[i] < samples[i] :
                        sampl_com = False 
                    com_size[i] = [com_size[i]]
        if not isint_com:
            raise ValueError("community sizes should be strictly positive int")
        if not sampl_com:
            raise ValueError('com_size must be superior to samples')
        
    # if isinstance(com_size, float) :
    #     com_size = [[int(com_size)]] * npop
    # check mu
    if not isinstance(mu, (int,float)) or mu < 0 or mu > 1 :
        raise ValueError('mu must be a float between 0 and 1')
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
                        raise ValueError("there should be as many past growth "+
                        "init_rates as there are epochs in changetime")
                else :
                    if len(init_rates) != npop :
                        raise ValueError("there should be as many elements in"+
                        " init_rates as there are demes")
                    if not isinstance(init_rates[i], (float, int)):
                        isint_rates = False 
                    init_rates[i] <- [init_rates[i]]
        if not isint_rates:
            raise ValueError('init_rates must be float, list of float or'+
                                 ' nested list of float')
        if len(init_rates) != npop :
            raise ValueError("there should be as many past sizes as there " + 
            "are epochs in init_rates")
    else :
        init_rates = [[0]] * npop

    # check mrca
    # if mrca is not None :
    #     # TODO : do this
    # check migr_time
    if migr_time is not None and migr is not None:
        if not isinstance(migr_time, list):
            raise ValueError("migr_time should be a list of int")
        if not all([isinstance(x, (int,float)) for x in migr_time]):
            raise ValueError("migr_time should be a list of int")
    else :
        migr_time = [0]
    # check migr & migr_time
    if migr is not None :
        if npop == 1 :
            # warnings.warn("no migration matrix is needed for a single deme")
            migr = None
    if migr is not None :
        if not isinstance(migr, list) : # case 'a
            if not isinstance(migr, (int,float)) :
                raise ValueError("migration rate must be a float or an int.")
            if migr < 0 or migr > 1 :
                raise ValueError("migration rate should be positive (or zero) and" + 
                " not exceed 1")
            # migr = np.ones((npop,npop))*migr
            # np.fill_diagonal(migr, 0)
            migr = [migr]
        else :
            for i in range(len(migr)):
                if not isinstance(migr[i], list): # case ['a, ... ,'b]
                    if len(migr) != len(migr_time):
                        raise ValueError("there should be as many migration rates" + 
                            " or matrices as there are times in migr_time")
                    if not isinstance(migr[i], (int,float)) :
                        raise ValueError("migration rate must be a float or an int.")
                    if migr[i] < 0 or migr[i] > 1 :
                        raise ValueError("migration rate should be positive (or zero)" + 
                                 " and not exceed 1")
                    # check len of migr is done with migr_time
                    migr[i] = np.ones((npop,npop))*migr[i]
                    np.fill_diagonal(migr[i], 0)
                else :
                    if not isinstance(migr[i][0], list) : # case [[0,'a], ['a, 0]]
                        if len(migr[i]) != len(migr) or len(migr) != npop:
                            raise ValueError("custom migration matrices should be of" + 
                                     " size ndeme x ndeme")
                        if not all([ isinstance(r, (float,int)) for r in migr[i]]) :
                            raise ValueError("found custom migration matrix that is" + 
                                     " not made of ints or floats")
                        if any([r<0 or r> 1 for r in migr[i]]):
                            raise ValueError("found custom migration matrix with" + 
                                " negative migration rates or greater than 1")
                    else: # case [[[0, 'a], ['a, 0]], [[0, 'b], ['b, 0]]]
                        if len(migr) != len(migr_time):
                            raise ValueError("there should be as many migration rates" + 
                                " or matrices as there are times in migr_time")
                        for j in range(len(migr[i])) :
                            if len(migr[i][j]) != len(migr[i]) or len(migr[i]) != npop:
                                raise ValueError("custom migration matrices should be" + 
                                    " of size ndeme x ndeme")
                            if any ([not isinstance(r, (float,int)) for r in migr[i][j]]) :
                                raise ValueError("found custom migration matrix that" + 
                                    " is not made of ints or floats")
                            if any ([r<0 or r> 1 for r in migr[i][j]]):
                                raise ValueError("found custom migration matrix with" + 
                                 " negative migration rates or greater than 1")

    m = np.array(migr)
    dim = m.shape
    if np.sum(m) == 0 :
        raise ValueError("migration matrices cannot all be empty")
    # check vic_events
    if vic_events is not None:
        if not isinstance(vic_events, list) or any([not isinstance(x, list) for x in vic_events]):
            raise ValueError("vic_events should be a nested list of list with "+
            "a length 3")
        if any([len(v)!=3 for v in vic_events]) :
            raise ValueError("all elements in vic_events should be lists of"+
            " lenght 3")
        
        if any([len(v[1])!=2 for v in vic_events]):
            raise ValueError("second element of vic_events should be a list"+
            " of 2 deme ids")
        # TODO : cath float and format them
        if not all([isinstance(v, int) for v in flatten(vic_events)]):
             raise ValueError("all elements of vic_events should be ints")
        
        if any([t<0 for t in [v[0] for v in vic_events]]):
            raise ValueError("all times in _vic_events should be strictly"+
            " positive")
        
        if any([test not in flatten(changetime) for test in [v[0] for v in vic_events]]):
            raise ValueError("split times in vic_events should also appear"+
            " in changetime")
        
        if not all([x in range(npop) for x in sum([flatten(v[1:]) for v in vic_events],[]) ]):
            raise ValueError("Split events do not match provided deme"+
            " information")
        
        vic_dates = [v[0] for v in vic_events]
        if vic_dates != sorted(vic_dates):
            raise ValueError("Split dates should be provided in chronological"+
            " order")
        
        if any([v[2] not in v[1] for v in vic_events]):
            raise ValueError("Splits events of two demes should be defined"+
            " relative to one of the demes' id")
        
        for i in range(len(vic_events)) :
            if i == 0 :
                continue
            if vic_events[i][1][0] in vic_events[i-1][1] and vic_events[i][1][0] != vic_events[i-1][2]:
                raise ValueError("Trying to merge with inactive deme")
            if vic_events[i][1][1] in vic_events[i-1][1] and vic_events[i][1][1] != vic_events[i-1][2]:
                raise ValueError("Trying to merge with inactive deme")
    # check verbose
    if not isinstance(verbose, bool):
        raise ValueError('verbose must be a boolean') 
    # check seed
    if seed is not None and not isinstance(seed, (int,float)):
        raise ValueError('seed must be an integer')
    if seed is not None and isinstance(seed, float):
        seed = int(seed)
  
    ############################################################################
    ####                    Now we can compute something                    ####
    ############################################################################
    # compute number of populations
    npop = len(samples)

    # initialise demography object for msprime.sim_ancestry()
    demography = msprime.Demography()
    
    # Build samples
    samples = {"pop_" + str(i):samples[i] for i in range(len(samples))}
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

    ## VICARIANCE EVENTS
    if vic_events is not None:
        # extract some informations about 
        nvic = len(vic_events)
        ancestrals = [str(vic_events[i][2]) for i in range(nvic)] 
        count = {}
        for i, anc in enumerate(ancestrals):
            cnt = count.get(anc, 0)
            count[anc] = cnt + 1
            ancestrals[i] += chr(ord('a') + cnt)

        ancestrals = ["pop_" + a for a in ancestrals]
        derived = []
        # add the ancestral pop and their size change
        for v in range(nvic):
            demography.add_population(
                name = ancestrals[v],
                initial_size= com_size[vic_events[v][2]][changetime[vic_events[v][2]].index(vic_dates[v])])
            tmp = changetime[vic_events[v][2]].index(vic_dates[v]) + 1
            an_changetime = changetime[vic_events[v][2]][tmp:]
            an_com_size = com_size[vic_events[v][2]][tmp:]
            for i in range(len(an_changetime)):
                demography.add_population_parameters_change(
                    population=ancestrals[v],
                    time = an_changetime[i],
                    initial_size= an_com_size[i])
            derived.extend(vic_events[v][1])

        #set up split events
        derived = [str(o) for o in derived]
        count = {}
        for i, o in enumerate(derived):
            cnt = count.get(o, 0)
            count[o] = cnt + 1
            if cnt > 0:
                derived[i] += chr(ord('a') + cnt - 1) 
        derived = ["pop_" + d for d in derived]
        derived = [derived[i*len(derived) // nvic: (i+1)*len(derived) // nvic] for i in range(nvic)] 
        
        for v in range(nvic):
            demography.add_population_split(time = vic_events[v][0], 
                                            derived = derived[v], 
                                            ancestral = ancestrals[v])

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
    node_labels = {u: str(u)+'_'+str(tree.population(u)) for u in tree.nodes() if tree.is_sample(u)}
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


def getAbund(tree, samples = None):
    """
    
    Parameters
    ----------
    tree : (ete3 class)
        Phylogeny with attributes on leafs. This attributes is a character 
        string containing all names of the species individual (mean to use 
        topPhylo result). Names is formated like this :
          " name1 name2 name3"
    samples : int
        number of individual in the community
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
        raise ValueError('tree must have a class TreeNode')
    if samples != None:
      if not isinstance(samples, int):
          raise ValueError('samples must be an integer')

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
    
    if samples != None and sum(sfs) != samples:
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
        raise ValueError('tree must have a class TreeNode')
    
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


def flatten(x):
    """
    lizzy need to document this because thx SO
    """
    if isinstance(x, collections.Iterable) and not isinstance(x, str):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]


if __name__ == "__main__":
        import doctest
        doctest.testmod()
        #simulate(samples = [5, 5], com_size = [[500], [500]], mu = 0.05, migr = 1, verbose = True, seed = 42)
        # simulate(samples = [10], com_size = [[500, 1000]], mu = 0.05, changetime= [[0,100]], seed = 42, verbose = True )
        # t = simulate(samples = [5, 5], com_size = [[1e3, 2e3], [1e3, 5e2]], changetime = [[0, 50],[0, 30]], mu = 0.03, migr = 2, verbose = True)
        # print(t)

        ## SINGLE POP
        # stat discret
        # t = simulate(samples = [5], com_size = [[1e3]], mu = 0.03, migr = 2, seed = 42, verbose = True)
        # stat continue 
        # t = simulate(samples = [5], com_size = [[1e3]], stable_pop = False, mu = 0.03, migr = 2, seed = 42, verbose = True)
        
        # fluct discret
        # t = simulate(samples = [5], com_size = [[1e3, 2e3]], changetime = [[0, 50]], mu = 0.03, migr = 2, seed = 42, verbose = True)

        ## MULT POP

        # stat discret
        # t = simulate(samples = [5, 5], com_size = [[1e3], [1e3]], mu = 0.03, migr = 1, seed = 42, verbose = True)
        # stat continue # 
        # t = simulate(samples = [5, 5], com_size = [[1e3], [1e3]], stable_pop = False, mu = 0.03, migr = 2, seed = 42, verbose = True)

        # fluct discret
        # t = simulate(samples = [5, 5], com_size = [[1e3, 2e3], [1e3, 5e2]], changetime = [[0, 50],[0, 30]], mu = 0.03, migr = 2, seed = 42, verbose = True)

        # t = simulate(samples = [5, 5], com_size = [[1e3, 2e3], [1e3, 5e2]], changetime = [[0, 500],[0, 300]], mu = 0.03, migr = [0, 1], migr_time = [0, 200], seed = 42, verbose = True)
        # print(t)

        # t = simulate(samples = [5, 5], com_size = [[1e3, 2e3], [1e3, 5e2]], changetime = [[0, 50],[0, 30]], mu = 0.03, migr = 1, seed = 42, verbose = True)
        # print(t)

