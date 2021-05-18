# -*- coding: utf-8 -*-
"""

Created on Wed May 13 11:42:50 2020

@author : Maxime Jaunatre <maxime.jaunatre@yahoo.fr>
@author : Elizabeth Bathelemy <barthelemy.elizabeth@gmail.com>

Functions :
    dosimuls
    simulate
    sample
    getAbund
    getDeme
    check_params

"""
    # TODO : general doc !!!
import msprime
import random
import numpy as np
from ete3 import Tree
import pandas as pd
from scipy.stats import loguniform
import collections

from ecophylo import pastdemo
from ecophylo import phylogen
from ecophylo import sumstat

def dosimuls(nsim, samples, deme_sizes, mu, tau = 0, gr_rates = None, 
             changetimes = None, mrca = None, migr = 1, migr_times = None,
             splits = None, 
             verbose = False,
             output = ['Params'], # Params, Sumstat, Tree
             file_name = None, seed = None):
    """
    This function allows simulating large datasets over wide ranges of eco-
    evolutive parameters by repeatedly calling the simulate function and 
    retreiving summary statistics generated for different parameter values 
    drawn from specified distributions.
    
    This function returns a table of sampled parameter values used for the 
    simulations, a table of summary statistics and (if specified) exports the 
    simulated phylogeneties.

    Parameters
    ----------
    nsim : a positive int
        the number of simulations to run
    output = ['Params'] : list of str
        A list specifying which elements to return among the following options:
            - Params :  returns the table of parameter values used for the 
                        simulations 
            - Sumstats: returns the table of relative abundances as well as 
                        alpha diversity metrics per deme if multiple demes 
                        are simulated
            - Tree:     returns the simulated phylogenies to a specified file
    file_name = None :str
        a string specifying the file path to save the output
    seed = None : int
        set seed for entire simulation. Note that if seed is specified, all 
        nsim simulation will be identical.
        
    all other parameters are documented in the simulate function 

    Notes
    -----
        Whether a given parameter for the simulation should be drawn from a 
        prior distribution should be specified by remplacing the parameter 
        value by a list containing in order: the parameter bounds as well as a
        string specifying the shape of the distribution, as follows:
        [min_prior, max_prior, "sample_law"] 
        
        Implemented distributions include:
            - "uniform" for a uniform distibution
            - "log_unif" for a log-uniform distribution
      
    Examples
    --------
    #ecophylo.dosimuls(nsim = 5, samples = [10, 9], com_size = [[500, 1000], [2000, [1500, 5000, "uniform"], 6000]], mu = 0.001, migr = 0.5, changetime = [[0, 200], [0, 100, 500]])
    
    
    """
    # Idiotproof dosimul parameters
    # nsim
    if not isinstance(nsim, (float,int)) or nsim < 1 :
        raise ValueError("nsim must be a int value above 0")
    if isinstance(nsim, float) : nsim = int(nsim)
    # output
    if not isinstance(output, list) or any([not isinstance(x, str) for x in output]) :
        raise ValueError("output must be a list of string")
    if len(output) > 3 or any([ x not in ["Params", "Sumstat", "Trees"] for x in output]):
        raise ValueError("output must contain only 'Params', 'Sumstat' or/and 'Trees'")
    # file_name
    if file_name is not None and not isinstance(file_name, str) :
        raise ValueError("file_name must be string value.")
            

    # register None values before they are initiated
    Nonedef = [gr_rates, changetimes, mrca, migr_times, splits]
    for i in range(len(Nonedef)):
        if Nonedef[i] is not None : Nonedef[i] = True
        else : Nonedef[i] = False
    
    # IDIOTPROOF and Register locations
    samples, deme_sizes, mu, gr_rates, changetimes, mrca, migr, migr_times, \
                    splits, verbose, seed, prior_locate = check_params(
        samples = samples, deme_sizes = deme_sizes, mu = mu, gr_rates = gr_rates, 
        changetimes = changetimes, mrca = mrca, migr = migr, 
        migr_times = migr_times, splits = splits, 
        verbose = verbose, seed = seed, prior_locate = 'naive'
    )
    if seed == 324:
        raise Warning("JE NE MANGE PAS DE GRAINES !")

    npop = len(samples)
    priors = [prior_locate[p][0] for p in range(len(prior_locate))]

    ##########################################################################
    ####                    SETUP DATAFRAME                               ####
    ##########################################################################
    # Creating parameter dataframe
    # if prior locate is empty deal with this later
    comments = ("Simulating eco-evolutionary dynamics over the following" +
               " parameters ranges:")
    params = pd.DataFrame()

    prior_names = []
    for prior in prior_locate :
        if prior[0] == "samples": prior_names.append(f'samples_pop{prior[1]}')
        if prior[0] == "mu": prior_names.append('mu')
        if prior[0] == "deme_sizes":
            prior_names.append(f'deme_sizes_pop{prior[1]}_t{prior[2]}')
        if prior[0] == "changetimes":
            prior_names.append(f'time_pop{prior[1]}_t{prior[2]}')
        if prior[0] == "gr_rates":
            prior_names.append(f'rate_pop{prior[1]}_t{prior[2]}')
        if prior[0] == "migr": prior_names.append(f'migr_t{migr_times[prior[1]]}')
    ## samples
    col_samples = ['samples_pop{}'.format(i) for i in range(len(samples))]
    for j in range(len(samples)):
        if col_samples[j] in prior_names: params[col_samples[j]] =  [None]*nsim
        else: params[col_samples[j]] = [samples[j]]*nsim
    ## mu
    if "mu" in priors: params["mu"] = [None]*nsim 
    else: params["mu"] = [mu]*nsim
    ## deme_sizes
    for j in range(len(deme_sizes)):
        for jj in range(len(deme_sizes[j])):
            tmp_name = f'deme_sizes_pop{j}_t{jj}'
            if tmp_name in prior_names: params[tmp_name] = [None] * nsim
            else : params[tmp_name] = [deme_sizes[j][jj]] * nsim
    ## changetimes
    if Nonedef[1]:
        for i in range(npop):
            for ii in range(len(changetimes[i])):
                if ii != 0:
                    tmp_name = f'time_pop{i}_t{ii}'
                    if tmp_name in prior_names: params[tmp_name] = [None] * nsim
                    else : params[tmp_name] = [changetimes[i][ii]] * nsim
    ## gr_rates
    if Nonedef[0]:
        for i in range(npop):
            for ii in range(len(gr_rates[i])):
                tmp_name = f'rates_pop{i}_t{ii}'
                if tmp_name in prior_names: params[tmp_name] = [None] * nsim
                else : params[tmp_name] = [gr_rates[i][ii]] * nsim
    ## migr and migr_times
    if migr is not None and isinstance(migr[0], (int, float)):
        for i in range(len(migr)):
            tmp_name = f'migr_t{migr_times[i]}'
            if tmp_name in prior_names: params[tmp_name] =  [None]*nsim
            else: params[tmp_name] = [migr[i]]*nsim
        # # This should work for matrices but we don't want to fuck our brains
        # col_migr = []
        # migr_val = []
        # for mt in range(len(migr)):
        #     if isinstance(migr[mt], (int, float)):
        #         migr[mt] = np.ones((npop,npop))*migr[mt]
        #         np.fill_diagonal(migr[mt], 0)
        #     x = np.array(migr[mt])
        #     XX,YY = np.meshgrid(np.arange(x.shape[1]),np.arange(x.shape[0]))
        #     table = np.vstack((x.ravel(),XX.ravel(),YY.ravel())).T
        #     col_migr.extend([f'migr_pop{int(t[1])}_pop{int(t[2])}_t{mt}' for t in table if np.unique(t[1:]).size != 1])
        #     migr_val.extend([t[0] for t in table if np.unique(t[1:]).size != 1])
        # for j in range(len(col_migr)):
        #     if 'migr' in priors:
        #         params[col_migr[j]] =  [None]*nsim
        #     else: 
        #         params[col_migr[j]] = [migr_val[j]]*nsim

    # print(params)

    result = list()
    abund = list()
    diver = list()
    trees = ""
    print("\n****SIMULATING***************\n\n")
    failed = 0
    i = 0
    while i < nsim and failed < 100:
        # SIMULATE
        phylo = simulate(
            samples = samples, deme_sizes = deme_sizes, mu = mu, tau = tau,
            gr_rates = gr_rates, changetimes = changetimes, 
            mrca = mrca, migr = migr, migr_times = migr_times, 
            splits = splits, verbose = False, seed = seed, force = True
        )
        #################################################################
        ####                    CHECK TREE                           ####
        #################################################################
        # check tree
        sp  = phylo.get_leaf_names()
        if len(sp) > 2 and len(set(sp)) > 1 :
            for ii in range(len(prior_locate)):
                tmp_p = prior_locate[ii]
                if tmp_p[0] == "samples":
                    print('not done') # TODO : samples 
                if tmp_p[0] == "mu":
                    params.loc[i,'mu'] = mu
                if tmp_p[0] == "deme_sizes":
                    params.loc[i,(f'deme_sizes_pop{tmp_p[1]}_t{tmp_p[2]}')] = \
                        deme_sizes[tmp_p[1]][tmp_p[2]]
                if tmp_p[0] == "changetimes":
                    params.loc[i,(f'time_pop{tmp_p[1]}_t{tmp_p[2]}')] = \
                        changetimes[tmp_p[1]][tmp_p[2]]
                if tmp_p[0] == "gr_rates":
                    params.loc[i,(f'rates_pop{tmp_p[1]}_t{tmp_p[2]}')] = \
                        gr_rates[tmp_p[1]][tmp_p[2]]
                if tmp_p[0] == 'migr':
                    params.loc[i,(f'migr_t{migr_times[tmp_p[1]]}')] = \
                        migr[tmp_p[1]]
            # Save sumstat
            abund.append(sumstat.getAbund(phylo))
            diver.append(sumstat.getDeme(phylo, div = True))
            # Save tree
            trees += phylo.write() + "\n"
            failed = 0
            i+=1
        else :
            failed += 1      
        #################################################################
        ####                    RESAMPLING                           ####
        #################################################################
        samples, deme_sizes, mu, gr_rates, changetimes, mrca, migr, migr_times, \
                    splits, verbose, seed, prior_locate = check_params(
            samples = samples, deme_sizes = deme_sizes, mu = mu, gr_rates = gr_rates, 
            changetimes = changetimes, mrca = mrca, migr = migr, 
            migr_times = migr_times, splits = splits, 
            verbose = verbose, seed = seed, prior_locate = prior_locate
        )

    if failed >= 100 : raise ValueError("Too many simulations have failed")

    if verbose:
        print(comments)

    # abund
    abund_pd = np.zeros((nsim, max(len(x) for x in abund)))
    for k, j in enumerate(abund):
        abund_pd[k][0:len(j)] = j
    abund_pd = pd.DataFrame(abund_pd)
    # diver
    diver_pd = pd.DataFrame(diver)
    print(['alpha'+ str(x) for x in range(diver_pd.shape[1] )])
    diver_pd.columns = ['alpha'+ str(x) for x in range(diver_pd.shape[1] )]
    print(diver_pd)
    sumstats = pd.concat([abund_pd.reset_index(drop=True), diver_pd], axis=1)

    print(sumstats)

    if file_name is not None :
        saved = ""
        if "Params" in output:
            saved += params.to_string()
        if "Sumstat" in output :
            saved += "\n###\n" + sumstats.to_string()
        if "Trees" in output : 
            saved += "\n###\n" + trees

        print('\nSimulation saved in :', file_name, '\n')
        f = open(file_name, "w")
        f.write(saved)
        f.close()

    result = [params, sumstats]
    return result

def check_params(samples, deme_sizes, mu, gr_rates = None, 
                   changetimes = None, mrca = None, 
                   migr = 1, migr_times = None, splits = None,
                   verbose = False, seed = None, prior_locate = None):
    """
    Internal function used to check parameters of simulate and dosimulate 
    functions. Parameters thus need to have the same format as described in 
    simulate function documentation.
    Can be used to sample priors, see dosimulate documentation

    Parameters
    ----------
    prior_locate = None : str in "naive" or list with specific format
        If none, will check parameters format
        If naive will sample priors and locate them in itself
        If lsit will use info stored inside for resampling

    """
    # Idiotproof prior_locate 
    if prior_locate is not None and prior_locate != "naive" and \
                                     not isinstance(prior_locate, list):
        prior_locate = None
        # TODO : warning here maybe and more idiotproof

    # prior_locate will pass this if educated
    if isinstance(prior_locate, list) :
        # draw and return prior values
        for prior in prior_locate :
            if prior[0] == "samples": # draw prior for qamples
                samples[prior[1]] = round(
                    sample(prior[3][0], prior[3][1], prior[3][2], seed = seed)
                )
            if prior[0] == "mu": # draw prior for mu
                mu = sample(prior[3][0], prior[3][1], prior[3][2], seed = seed)
            if prior[0] == "deme_sizes": # draw prior for deme_sizes
                deme_sizes[prior[1]][prior[2]] = round(
                    sample(prior[3][0], prior[3][1], prior[3][2], seed = seed)
                )
            if prior[0] == "changetimes": # draw prior for changetimes
                changetimes[prior[1]][prior[2]] = round(
                    sample(prior[3][0], prior[3][1], prior[3][2], seed = seed)
                )
            if prior[0] == "gr_rates": # draw prior for gr_rates
                gr_rates[prior[1]][prior[2]] = sample(
                    prior[3][0], prior[3][1], prior[3][2], seed = seed
                )
            if prior[0] == "migr":
                migr[prior[1]] = sample(
                    prior[3][0], prior[3][1], prior[3][2], seed = seed
                )

    #################################################################
    ####                    CHECK ALL AND SAMPLE                 ####
    #################################################################
    else :
        if prior_locate == "naive":
            prior_locate = []
            # else prior_locate is None
        # Idiotproof
        # check samples
        if not isinstance(samples, list):
            samples = [samples]
        for i in range(len(samples)):
            if prior_locate is not None and isinstance(samples[i], list):
                prior_locate.append(["samples", i, 0, samples[i]])
                samples[i] = round(sample(
                    samples[i][0], samples[i][1], samples[i][2], seed = seed
                ))
            if not isinstance(samples[i], int):
                raise ValueError("samples should all be ints")
            if samples[i] <= 0:
                raise ValueError("samples should all be positive")
        # compute number of populations
        npop = len(samples)

        # check changetimes
        if changetimes is not None:
            if not isinstance(changetimes, list):
                if isinstance(changetimes, (int,float)) : 
                    if changetimes >= 0 :
                        changetimes = [[changetimes]]
                    else :
                        raise ValueError("changetimes must be positive values")
                else :
                    raise ValueError("changetimes must be int, list of int or"+
                    " nested list of int")
            else :
                for i in range(len(changetimes)):
                    if isinstance(changetimes[i], list):
                        for ii in range(len(changetimes[i])):
                            # draw priors
                            if prior_locate is not None and isinstance(changetimes[i][ii], list):
                                prior_locate.append(["changetimes", i, ii, changetimes[i][ii]])
                                changetimes[i][ii] = round(sample(
                                    changetimes[i][ii][0], changetimes[i][ii][1],
                                    changetimes[i][ii][2], seed = seed
                                ))
                        if not all(isinstance(y, (float, int)) for y in changetimes[i]) :
                            raise ValueError("changetimes must be int, list of"+
                            " int or nested list of int")
                        if any(y < 0 for y in changetimes[i][1:]) : 
                            raise ValueError("changetimes must be positive values")
                        if changetimes[i][0] != 0:
                            raise ValueError("first element of changetimes for a Deme"+
                            " must be equal to 0")
                        if len(set(changetimes[i])) != len(changetimes[i]) :
                            raise ValueError("Duplicated times in changetimes for a Deme" +
                                     " are not possible")
                    else :
                        if len(set(changetimes)) != len(changetimes) :
                            raise ValueError("Duplicated times in changetimes are not possible")
                        if not isinstance(changetimes[i], (float, int)):
                            raise ValueError("changetimes must be int, list of int or"+
                                     " nested list of int")
                        if changetimes[i] < 0 :
                            raise ValueError("changetimes must be positive values")
                        if changetimes[0] != 0:
                            raise ValueError("first element of changetimes"+
                            " must be equal to 0")
                if not isinstance(changetimes[i], list):
                    changetimes = [changetimes]
            if len(changetimes) != npop :
                raise ValueError("there should be as many past sizes as there " + 
                "are epochs in changetimes")
        else :
            changetimes = [[0]] * npop
    
        # check deme_sizes 
        if changetimes is not None  :
            isint_com = True
            sampl_com = True
            if not isinstance(deme_sizes, list) :
                if isinstance(deme_sizes, (int,float)) and deme_sizes > 0 : 
                    if deme_sizes < samples[0] : sampl_com = False
                    deme_sizes = [[int(deme_sizes)]] * npop
                else :
                    isint_com = False
            else :
                for i in range(len(deme_sizes)):
                    if isinstance(deme_sizes[i], list):
                        for ii in range(len(deme_sizes[i])): 
                            # draw priors
                            if prior_locate is not None and isinstance(deme_sizes[i][ii], list):
                                prior_locate.append(["deme_sizes", i, ii, deme_sizes[i][ii]])
                                deme_sizes[i][ii] = round(sample(
                                    deme_sizes[i][ii][0], deme_sizes[i][ii][1],
                                    deme_sizes[i][ii][2], seed = seed
                                ))

                        if not all(isinstance(y, (float, int)) for y in deme_sizes[i]) :
                            isint_com = False
                        else :
                            deme_sizes[i] = [int(x) if isinstance(x, float) else x for x in deme_sizes[i]]
                        if len(deme_sizes) != npop :
                            raise ValueError("there should be as many elements in"+
                                     " deme_sizes as there are demes")
                        if len(deme_sizes[i]) !=  len(changetimes[i]) :
                            raise ValueError("there should be as many past "+
                            "deme_sizes as there are epochs in changetimes")
                        if isint_com and any([s <= 0 for s in deme_sizes[i]]):
                            raise ValueError("all past sizes should be strictly positive")
                        if isint_com and any([x < samples[i] for x in deme_sizes[i]]):
                            sampl_com = False 
                    else :
                        if len(deme_sizes) != npop :
                            raise ValueError("there should be as many elements in"+
                                     " deme_sizes as there are demes")
                        if not isinstance(deme_sizes[i], (float, int)):
                            isint_com = False 
                        if isint_com and deme_sizes[i] <= 0 :
                            raise ValueError("all past sizes should be strictly positive")
                        if isint_com and deme_sizes[i] < samples[i] :
                            sampl_com = False 
                        deme_sizes[i] = [deme_sizes[i]]
            if not isint_com:
                raise ValueError("community sizes should be strictly positive int")
            if not sampl_com:
                raise ValueError("deme_sizes must be superior to samples")
            
        # if isinstance(deme_sizes, float) :
        #     deme_sizes = [[int(deme_sizes)]] * npop
        # check mu
        # draw prior
        if isinstance(mu, list) and len(mu) == 3 and prior_locate is not None and \
             all([x > 0 and x < 1 for x in mu[:2]]):
            prior_locate.append(["mu", 0, 0, mu])
            mu = sample(mu[0], mu[1], mu[2])
        if not isinstance(mu, (int,float)) or mu < 0 or mu > 1 :
            raise ValueError("mu must be a float between 0 and 1")
        # check gr_rates
        if gr_rates is not None and changetimes is not None:
            isint_rates = True
            if not isinstance(gr_rates, list):
                if isinstance(gr_rates, (int,float)) : 
                    gr_rates = [[gr_rates]] * npop
                else :
                    isint_rates = False
            else :
                for i in range(len(gr_rates)):
                    if isinstance(gr_rates[i], list):
                        for ii in range(len(gr_rates[i])): 
                            # draw priors
                            if prior_locate is not None and isinstance(gr_rates[i][ii], list):
                                prior_locate.append(["gr_rates", i, ii, deme_sizes[i][ii]])
                                gr_rates[i][ii] = round(sample(
                                    gr_rates[i][ii][0], gr_rates[i][ii][1],
                                    gr_rates[i][ii][2], seed = seed
                                ))
                        if not all(isinstance(y, (float, int)) for y in gr_rates[i]) :
                            isint_rates = False
                        if len(gr_rates[i]) !=  len(changetimes[i]) :
                            raise ValueError("there should be as many past growth "+
                            "gr_rates as there are epochs in changetimes")
                    else :
                        if len(gr_rates) != npop :
                            raise ValueError("there should be as many elements in"+
                            " gr_rates as there are demes")
                        if not isinstance(gr_rates[i], (float, int)):
                            isint_rates = False 
                        gr_rates[i] <- [gr_rates[i]]
            if not isint_rates:
                raise ValueError("gr_rates must be float, list of float or"+
                                     " nested list of float")
            if len(gr_rates) != npop :
                raise ValueError("there should be as many past sizes as there " + 
                "are epochs in gr_rates")
        else :
            gr_rates = [[0]] * npop
    
        # check mrca
        # if mrca is not None :
        #     # TODO : do this
        # check migr_times
        if migr_times is not None and migr is not None:
            if not isinstance(migr_times, list):
                raise ValueError("migr_times should be a list of int")
            if not all([isinstance(x, (int,float)) for x in migr_times]):
                raise ValueError("migr_times should be a list of int")
        else :
            migr_times = [0]
        # check migr & migr_times
        if migr is not None :
            if npop == 1 :
                # raise Warning("no migration matrix is needed for a single deme")
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
                        if len(migr) != len(migr_times):
                            raise ValueError("there should be as many migration rates" + 
                                " or matrices as there are times in migr_times")
                        if not isinstance(migr[i], (int,float)) :
                            raise ValueError("migration rate must be a float or an int.")
                        if migr[i] < 0 or migr[i] > 1 :
                            raise ValueError("migration rate should be positive (or zero)" + 
                                     " and not exceed 1")
                        # check len of migr is done with migr_times
                        # migr[i] = np.ones((npop,npop))*migr[i]
                        # np.fill_diagonal(migr[i], 0)
                    else :
                        # draw prior
                        if prior_locate is not None and len(migr[i]) == 3 and isinstance(migr[i][2], str) :
                            prior_locate.append(["migr", i, 0, migr[i]])
                            migr[i] = sample(
                                migr[i][0], migr[i][1], migr[i][2], seed = seed
                            )
                        elif not isinstance(migr[i][0], list) : # case [[0,'a], ['a, 0]]
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
                            if len(migr) != len(migr_times):
                                raise ValueError("there should be as many migration rates" + 
                                    " or matrices as there are times in migr_times")
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
        # check splits
        if splits is not None:
            if not isinstance(splits, list):
                raise ValueError("splits should be a nested list of list with "+
                "a length 3")
            for i in range(len(splits)) :
                if not isinstance(splits[i], list):
                    raise ValueError("splits should be a nested list of list with "+
                                 "a length 3")
                if len(splits[i]) != 3 :
                    raise ValueError("all elements in splits should be lists of"+
                                         " lenght 3")
                if len(splits[i][0])!=2:
                    raise ValueError("first element of splits should be a list"+
                    " of 2 deme ids")
                # # TODO : cath float and format them
                if not all([isinstance(x, int) for x in flatten(splits[i])]):
                     raise ValueError("all elements of splits should be ints")
                if splits[i][2] < 0 or any([ x < 0 for x in splits[i][0]]):
                    raise ValueError("all times in _splits should be strictly"+
                    " positive")
                
                if any(splits[i][0]) > npop or splits[i][1] > npop  :
                    raise ValueError("Split events do not match provided deme"+
                    " information")

                if splits[i][2] >= len(changetimes[splits[i][1]]):
                    raise ValueError("split times in splits should also appear"+
                                       " in changetimes")
                if splits[i][1] not in splits[i][0]:
                    raise ValueError("Splits events of two demes should be defined"+
                    " relative to one of the demes' id")
            
            vic_dates = [v[2] for v in splits]
            if vic_dates != sorted(vic_dates):
                raise ValueError("Split dates should be provided in chronological"+
                " order")
            
            for i in range(len(splits)) :
                if i == 0 :
                    continue
                if splits[i][0][0] in splits[i-1][0] and splits[i][0][0] != splits[i-1][1]:
                    raise ValueError("Trying to merge with inactive deme")
                if splits[i][0][1] in splits[i-1][0] and splits[i][0][1] != splits[i-1][1]:
                    raise ValueError("Trying to merge with inactive deme")

        # check coalesc
        if migr is None and splits is None and npop > 1:
            raise ValueError("Multiple demes must be linked by either migration"+
            " or vicariance events in order to coalesce")
        # check verbose
        if not isinstance(verbose, bool):
            raise ValueError("verbose must be a boolean") 
        # check seed
        if seed is not None and not isinstance(seed, (int,float)):
            raise ValueError("seed must be an integer")
        if seed is not None and isinstance(seed, float):
            seed = int(seed)

    return samples, deme_sizes, mu, gr_rates, changetimes, mrca, migr, \
                    migr_times, splits, verbose, seed, prior_locate



def simulate(samples, deme_sizes, mu, tau = 0, spmodel = "SGD",
             gr_rates = None, changetimes = None, mrca = None, 
             migr = 1, migr_times = None, splits = None,
             verbose = False, seed = None, force = False):
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
    deme_sizes : int or nested list of ints
        the size of Jm for each deme at each given period. Should be a nested
        list containing for each deme, a list of past Jm sizes in which the
        first element is the current size of the assemblage and the nth element
        is the size of Jm at epoch n 
    mu : float
        the point mutation rate must be comprised between between 0 and 1.
    gr_rates: float or nested list of floats
        the growth rates for each deme at each given period. Should be a nested
        list containing for each deme, a list of past growth rates in which the
        first element is the current growth rate and the nth element is the
        growth rate at epoch n. If no growth rates are given, then changes in
        deme sizes occur instantenously following sizes provided in deme_sizes at
        the different times given in changetimes
    changetimes: list of int or nested list of int
        the times (in generation before present) at which either growth rates or
        the size of the assemblages Jm have changed. If multiple demes are to be
        simulated, should be a nested list containing for each deme, a list of
        times at which changes occured in which the first element is 0.
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
        contain as many matrices M as there are time periods in migr_times where
        M[j,k] is the rate at which individuals move from deme j to deme k in
        the coalescent process, backwards in time. Individuals that move from
        deme j to k backwards in time actually correspond to individuals
        migrating from deme k to j forwards in time.
    migr_times = None: list of ints
        the times (in generation before present) at which migration rates have
        changed in which the first element is 0
    splits = nested list of ints # TODO : change doc for splits
        a nested list detailing the different split events that should be 
        included in the simulation. Each element of splits should be a list
        specifying, in order: the date (in generations before present) at which
        the split occured, the demes resulting from the split (as a list of ints)
        and finally the ancestral deme number. For instance, if deme 1 splits 
        into deme 0 and deme 1 then splits =  [[time01, [0,1], 1]]
        Note that time01 should appear in changetimes. Also, user should specify
        in deme_sizes (at the correct position i.e to the size of the ancestral 
        deme at time01) the size of the ancestral deme when the split occurs. 
    verbose = False : bool
        whether or not to print a summary of the demographic history and the
        resulting genealogy to be passed to a phylogeny
    seed = None : int
        set seed for entire simulation
    
    Examples
    --------
    >>> t = simulate(samples = [10], deme_sizes = [[1e5]], mu = 0.03, seed = 42)
    >>> print(t)
    <BLANKLINE>
          /-sp1
       /-|
      |  |   /-sp2
      |   \-|
    --|      \-sp3
      |
      |   /-sp4
      |  |
       \-|      /-sp5
         |   /-|
          \-|   \-sp6
            |
             \-sp7
    >>> t = simulate(samples = [5, 5], deme_sizes = [[1e5], [1e5]], 
    ... mu = 0.03, migr = 1, seed = 42)
    >>> print(t)
    <BLANKLINE>
          /-sp1
         |
       /-|      /-sp2
      |  |   /-|
      |   \-|   \-sp3
    --|     |
      |      \-sp4
      |
      |   /-sp5
       \-|
          \-sp6
    """ 
    # Idiotproof
    if not force :
        samples, deme_sizes, mu, gr_rates, changetimes, mrca, migr, migr_times,\
                        splits, verbose, seed, prior_locate = check_params(
            samples = samples, deme_sizes = deme_sizes, mu = mu,  
            gr_rates = gr_rates, changetimes = changetimes,
            mrca = mrca, migr = migr, migr_times = migr_times,
            splits = splits, verbose = verbose, seed = seed, 
            prior_locate = None
        )
  
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
            initial_size= deme_sizes[pop][0],       
            growth_rate=gr_rates[pop][0])
        
        # if population sizes have fluctuated in the past:
        if len(changetimes[pop]) > 1 and len(deme_sizes[pop]) > 1:
            for i in range(len(changetimes[pop][1:])):
                demography.add_population_parameters_change(
                    time = changetimes[pop][i+1] , 
                    initial_size=deme_sizes[pop][i+1], 
                    population= pop_ids[pop])
        
        # if population growth rates have fluctuated in the past:
        if len(changetimes[pop]) > 1 and len(gr_rates[pop]) > 1:
            for i in range(len(changetimes[pop][1:])):
                demography.add_population_parameters_change(
                    time = changetimes[pop][i+1] , 
                    growth_rate=gr_rates[pop][i+1], 
                    population= pop_ids[pop])

    ## VICARIANCE EVENTS
    if splits is not None:
        vic_dates = [changetimes[v[1]][v[2]] for v in splits]
        # extract some informations about 
        nvic = len(splits)
        ancestrals = [str(splits[i][1]) for i in range(nvic)] 
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
                initial_size= deme_sizes[splits[v][1]]\
                           [changetimes[splits[v][1]].index(vic_dates[v])])
            tmp = changetimes[splits[v][1]].index(vic_dates[v]) + 1
            an_changetimes = changetimes[splits[v][1]][tmp:]
            an_deme_sizes = deme_sizes[splits[v][1]][tmp:]
            for i in range(len(an_changetimes)):
                demography.add_population_parameters_change(
                    population=ancestrals[v],
                    time = an_changetimes[i],
                    initial_size= an_deme_sizes[i])
            derived.extend(splits[v][0])

        #set up split events
        derived = [str(o) for o in derived]
        count = {}
        for i, o in enumerate(derived):
            cnt = count.get(o, 0)
            count[o] = cnt + 1
            if cnt > 0:
                derived[i] += chr(ord('a') + cnt - 1) 
        derived = ["pop_" + d for d in derived]
        derived = [derived[i*len(derived) // nvic: (i+1)*len(derived) // nvic]\
                                                          for i in range(nvic)] 
        
        for v in range(nvic):
            demography.add_population_split(time = vic_dates[v], 
                                            derived = derived[v], 
                                            ancestral = ancestrals[v])

    ## MIGRATION
    if migr is not None :
        m = np.array(migr)
        dim = m.shape
        if len(dim) == 1 :
            # symmetric migration matrix
            demography.set_symmetric_migration_rate(
                populations = range(npop), rate = migr[0])

            # if symmatric migration rate has changed in the past:
            if len(migr) > 1:
                for m in range(len(migr[1:])):
                    demography.add_migration_rate_change(
                        time = migr_times[m+1], rate = migr[m+1])

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
                for t in range(len(migr_times)):
                    for row in range(npop):
                        for col in range(npop):
                            if migr[t][row][col] == 0 :
                                continue
                            demography.add_migration_rate_change(
                                time = migr_times[t], rate = migr[t][row][col], 
                                source=pop_ids[row], dest=pop_ids[col])

    # sort events chronologically
    demography.sort_events()

    # if verbose should print the demography debugger - only for debugging !!! 
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
    node_labels = {u: str(u)+'_'+str(tree.population(u)) for u in tree.nodes()\
                                                           if tree.is_sample(u)}
    tree = Tree(tree.newick(node_labels = node_labels))
    phylo = phylogen.toPhylo(
        tree= tree, mu= mu, tau= tau, spmodel = spmodel, seed= seed
    )

    return phylo


def sample(lower, upper, distr = "uniform", seed = None):
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
      float single values

    Examples
    --------
    >>> sample(lower=1, upper=6, distr="uniform", seed = 42)
    4.197133992289419
    """
    if not isinstance(lower, (int,float)):
            raise ValueError("lower must be an integer or a float")
    if not isinstance(upper, (int,float)):
            raise ValueError("upper must be an integer or a float")
    if not isinstance(distr, str):
            raise ValueError("distr must be a string in this list :"+
            " uniform"+
            ", log_unif"+
            ".")
    if seed is not None and not isinstance(seed, (int,float)):
            raise ValueError("seed must be an integer")

    random.seed(seed)
    if upper == lower :
        return upper

    if lower > upper :
        tmp = lower
        lower = upper
        upper = tmp

    if distr == "uniform":
        p = random.uniform(lower,upper)
    elif distr == "log_unif":
        p = loguniform(lower, upper).rvs()
    else :
        raise ValueError("This distribution is not implemented")
    return p 


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

