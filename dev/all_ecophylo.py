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
import random
from loguniform import LogUniform

import os
import datetime
import h5py

def dosimuls(nsim, sample_size, comprior, muprior, lim_mrca = None, sstype="SFS",
             prior_distrib = "uniform", npop=1, withmigr=False, init_ratesprior=None,
             init_sizeprior=None, pastprior=None, changetime = None,
             nsplit=None, massprior=None, migrfrom=None, migrto=None,
             verbose=False, savetrees= False, saveto = ""):

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

    while i < nsim:
        if safecount > 10000:
            sys.exit("Too many simulations have failed")

        try:
            df = df.append(pd.Series(), ignore_index=True)
            # sample parameters from prior for simulation
            if isinstance(sample_size, list):
                samp = sample(sample_size[0], sample_size[1])
                df['sampsize'] = samp
        
            if muprior is not None:
                if len(muprior) == 1:
                    mu = muprior[0]
                else:
                    mu = sample(muprior[0], muprior[1], distr = prior_distrib)
                df['mu'].iloc[i,] = mu

            if comprior is not None:
                if len(comprior) == 1:
                    com_size = comprior[0]
                else:
                    com_size = sample(comprior[0], comprior[1], distr = prior_distrib, typ = "int")
                df['comsize'].iloc[i,] = com_size
            
            if pastprior is not None:
                past_sizes = [sample(pastprior[0], pastprior[1], distr = prior_distrib, typ = "int") for _ in range(nepoch)]
                colnames = ['pastsize{}'.format(i) for i in range(1, nepoch+1)]
                for j in range(nepoch):
                    df[colnames[j]].iloc[i,] = past_sizes[j]

            if init_sizeprior is not None and init_ratesprior is not None and npop is not None:
                init_rates = [sample(init_ratesprior[0], init_ratesprior[1], distr = prior_distrib, typ = "float") for _ in range(npop)][0]
                init_sizes = [sample(init_sizeprior[0], init_sizeprior[1], distr = prior_distrib, typ = "float") for _ in range(npop)][0]
                colnames1 = ['initrate{}'.format(i) for i in range(1, npop+1)]
                colnames2 = ['initsize{}'.format(i) for i in range(1, npop+1)]
                for j in range(npop):
                    df[colnames1[j]].iloc[i,] = init_rates[j][0]
                    df[colnames2[j]].iloc[i,] = init_sizes[j][0]

            if withmigr and npop is not None:
                m = sample(0, 1, distr = prior_distrib, typ = "float")
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
                             verbose = False)
            
            if sstype == 'SFS':
                ss.append(getSFS(phylo, samp))

            if savetrees:
                trees += phylo.write() + "\n"
        except: 
            df = df[:-1]
        else:
            i = i+1
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

def simulate(sample_size, com_size, mu, mrca = None, npop = 1, m = 0, init_rates = None, init_sizes = None, past_sizes = None, changetime = None, split_dates = None, migrfrom = None, migrto = None, seed = 1, verbose = False):

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
    phylo = toPhylo(tree, mu)

    return phylo 

def sample(lower, upper, distr = "uniform", typ = "float"):
    if upper == lower :
        p = upper
    
    if distr == "uniform":
        if typ == "int":
            p = random.randrange(lower,upper)
        else:
            p = random.uniform(lower,upper)

    if distr == "log_unif":
        p =   LogUniform(lower, upper).rvs()

    return p 


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
        if distrib == "uniform":
            if typ == "int":
                p = np.random.randint(lim[0], lim[1], size = nsim, dtype=np.int64)
            if typ == "float":
                p = np.random.uniform(lim[0], lim[1], size = nsim)
        if distrib == "logunif":
            p = np.exp(np.random.uniform(lim[0], lim[1], size = nsim))
    return p


def timeframes(I, T, a):
    """
    Compute window frames.

    cf: Boitard et al 2016 Plos Genetics

    Parameters
    ----------
    I: int
        Number of time windows
    T: int
        Maximum time in generation time
    a: float
        Resolution focus on a particular period

    Returns
    -------
    a list of dates (in generation time) corresponding to the
    different time windows

    """
    I = I + 1
    times = [(np.exp((np.log(1+a*T)*i)/(I-1))-1)/a for i in range(1, I)]
    return(times)


def toPhylo(tree, mu, spmodel = "SGD", force_ultrametric = True):
    """
    

    Parameters
    ----------
    tree : TYPE
        DESCRIPTION.
    mu : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    innerNodeIndex = 0
    nIndsORI = 0
    spID = 0
    
    for node in tree.traverse("preorder"):
        try:
            node.sp
        except AttributeError:
            node.add_features(sp=1)

        if not node.is_leaf():
            node.name = "n%d" % innerNodeIndex
            innerNodeIndex += 1
        else:
            nIndsORI += 1

        if not node.is_leaf():
            umut = ubranch_mutation(node, mu)
            if umut:
                # print(f"Speciation event @ node {node.name}")
                spID += 1
                node.sp = spID
                for leaf in node:
                    try:
                        leaf.sp = spID
                    except AttributeError:
                        leaf.add_features(sp=1)
            # print(f"node {innerNodeIndex} --> sp: {node.sp}")

    if spmodel == "NTB" :
        traversedNodes = set()
        for node in tree.traverse("postorder"):
            if not node.is_leaf():
                children = node.get_children()
                csp = [i.sp for i in children]
                if csp.count(csp[0]) == len(csp):
                    mergedLeaves = ""
                    for childnode in node.traverse():
                        traversedNodes.add(childnode)
                        if childnode.is_leaf():
                            mergedLeaves = mergedLeaves+ " "+ childnode.name

                    children[1].mergedInd = mergedLeaves
                    children[0].delete()
        
        
        
    if spmodel == "SGD" : 
        traversedNodes = set()
        for node in tree.traverse("preorder"):
            if node not in traversedNodes:
                if not node.is_leaf():
                    children = node.get_children()
                    if len(children) != 2:
                        sys.exit("The algorithm does not know how to deal with non dichotomic trees!")
                    csp1 = set()
                    for j in children[0].iter_leaves():
                        csp1.add(j.sp)
                    csp2 = set()
                    for j in children[1].iter_leaves():
                        csp2.add(j.sp)
                    common = csp1.intersection(csp2)  # compares species label between children
                    if len(common) > 0:  # if paraphyletic
                        if not node.is_root():
                            upNode = node.up  # parent
                            newLeaf = node.get_farthest_leaf()  # finds new leaf
                            newDist = newLeaf[1] + node.dist

                            mergedLeaves = ""

                            for childnode in node.traverse():
                                traversedNodes.add(childnode)
                                if childnode.is_leaf():
                                    mergedLeaves = mergedLeaves+" "+childnode.name
                            node.detach()
                            upNode.add_child(newLeaf[0], newLeaf[0].name, newDist)
                            newLeaf[0].mergedInd = mergedLeaves
                    
                        else:
                            # populate "mergedInd" feature for future SFS
                            mergedLeaves = ""
                            for l in node.iter_leaves():
                                mergedLeaves = mergedLeaves+" "+l.name
                            # collapse the subtree
                            newLeaf = tree.get_farthest_leaf()
                            for child in tree.get_children():
                                child.detach()
                                node.add_child(newLeaf[0], newLeaf[0].name, newLeaf[1])
    
                            # actualize "mergedInd" feature of new leaf
                            newLeaf[0].mergedInd = mergedLeaves
    
    if force_ultrametric:
        tree_dist = tree.get_farthest_leaf()[1]
        for leaf in tree.iter_leaves():
            dst = tree.get_distance(leaf)
            if dst != tree_dist:
                leaf.dist += tree_dist - dst
    
    return tree


def ubranch_mutation(node, mu):
    """
    Draw mutations following a poisson process.

    Parameters
    ----------
    node : ete3.coretype.tree.TreeNode
        node from which to compute branch length
    mu : float
        mutation rate

    Returns
    -------
    bool
        whether or not a mutation should appear on the tree at this node

    """
    lambd = node.dist * mu
    rb = np.random.poisson(lambd)
    if rb >= 1:
        return True
    else:
        return False


def getSFS(tree, samp_size):
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
    sfs.extend(abund)
    # think about catching error when phylogeny has only 1 species! 
    if sum(sfs) != samp_size:
        raise Exception(f"Simulated phylogeny has only one species!")
    return sfs

def timeframes(I, T, a):
    """
    Compute window frames.

    cf: Boitard et al 2016 Plos Genetics

    Parameters
    ----------
    I: int
        Number of time windows
    T: int
        Maximum time in generation time
    a: float
        Resolution focus on a particular period

    Returns
    -------
    a list of dates (in generation time) corresponding to the
    different time windows

    """
    I = I + 1
    times = [(np.exp((np.log(1+a*T)*i)/(I-1))-1)/a for i in range(1, I)]
    return(times)


def demographic_events(epochs, sizes):
    """
    Change the demography of populations over given time periods.

    Parameters
    ----------
    epochs: list
        When the demographic changes have occured.
    sizes: list
        Population sizes at the different time periods.

    Returns
    -------
    a list object to be passed into msprime.simulate to change the demography 
    over given times periods. Changes affect Population sizes at the different
    times, growth rates are automatically set to 0 (constant size).Should later
    implement an option to change a specific population for now the changes 
    affect all populations simultaneously.


    """
    dc = [msprime.PopulationParametersChange(time=t, initial_size =s) for t, s in zip(epochs, sizes)]
    return dc


def population_configurations(samples, init_sizes, rates):
    """
    Set up the initial population configurations.

    Parameters
    ----------
    samples : list
        The different population sample sizes
    init_sizes : list
        The initial population sizes
    rates : list
        The initial population growth rates

    Returns
    -------
    a list object that can be passed into msprime.simulate to indicate initial 
    population configurations

    """
    pc = [msprime.PopulationConfiguration(sample_size = s, initial_size = i, growth_rate = g) for s, i, g in zip(samples, init_sizes, rates)]
    return pc

def migration_matrix(subpops, migr = 0):
    """
    Set up the migration matrix between sub-populations.
    
    Only supports symmetric island model for the moment, 
    migration rate is identical between all pairs of subpopulations. 
    
    Parameters
    ----------
    subpops: int
        Number of sub-populations. Should be at least 2.
    migr: float
        overall symetric migration rate. Default is 0. 

    Returns
    -------
    Given N populations, an NxN numpy array of between-subpopulation 
    migration rates. 

    """
    if subpops < 1:
        sys.exit("there should be at least 2 populations")
    
    m = migr / (2 * (subpops - 1))

    # symmetric island model (later - implement other types of models)
    migration_matrix = np.ones((subpops,subpops))*m
    np.fill_diagonal(migration_matrix, 0)

    return migration_matrix

def mass_migrations(times, sources, destinations, migr = 1):
    """
    Include mass migration events.

    Parameters
    ----------
    times: int or list
        One or more dates (in generation time) when mass migrations have occured
    sources : int or list
        The source population(s) from the different mass migration events
    destinations : int or list
        The destination population(s) from the different mass migration events
    migr : int or list
        The probability(ies) of immigrating from the source population(s) to 
        the destination populations(s). Default is 1 (merge/split populations)

    Returns
    -------
    a list object that can be passed into msprime.simulate to introduce mass
    migration events

    """
    # if only one mass migration event should still be of type list
    if np.isscalar(times):
        times = [times]
    if np.isscalar(sources):
        sources = [sources]
    if np.isscalar(destinations):
        destinations = [destinations]
    if np.isscalar(migr):
        migr = [migr]
    if len(migr) != len(times):
        migr = migr*len(times)
    if not len(times) == len(sources) == len(destinations) == len(migr):
        sys.exit("Argument lengths should match")
    
    M = [msprime.MassMigration(time = t, source = s, dest = d, proportion = m) for t, s, d, m in zip(times, sources, destinations, migr)]
    return M
