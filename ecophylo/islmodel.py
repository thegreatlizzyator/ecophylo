# -*- coding: utf-8 -*-
"""
islmodel.py file 

Module with functions creating msprime (dependence) objects to be used 
by the ecophylo.simulate function.

Created on Wed May 13 11:42:50 2020

@author: barthele

Functions :
    population_configurations
    migration_matrix
    mass_migrations

"""
# Dependencies

import msprime
import numpy as np
import sys
import math


def population_configurations(samples, init_sizes, rates) :
    """
    Set up the initial population configurations.

    Parameters
    ----------
    samples : list of int
        positive values
        The different population sample sizes
    init_sizes : list of int
        positive vlaues
        The initial population sizes
    rates : list of float
        # TODO : test rates limites (real, -1:1 or -inf:inf)
        The initial population growth rates

    Returns
    -------
    a list object that can be passed into 
    msprime.simulate to indicate initial 
    population configurations
    
    Examples
    -------
    >>> print("test")
    "no"
    """

    # TODO : examples
    
    if not isinstance(samples, list) :
        sys.exit('samples must be a list')
    if not isinstance(init_sizes, list) :
        sys.exit('init_sizes must be a list')
    if not isinstance(rates, list) :
        sys.exit('rates must be a list')
    if not len(samples) == len(init_sizes) == len(rates) :
        sys.exit('all parameters must be list of same lengths')
        
    # if not all(isinstance(x, int) for x in past_sizes) :
    #     sys.exit('past_sizes must be a list of int')

    pc = [msprime.PopulationConfiguration(sample_size = s, initial_size = i, growth_rate = g) for s, i, g in zip(samples, init_sizes, rates)]
    return pc

def sizes2rates(init_size, past_sizes, times):
    """
    Compute the growth rates corresponding to a set of past sizes at different 
    times for a single population with a given initial size

    Parameters
    ----------
    init_size : int
        positive value
        the initial size of the population
    past_sizes: list of int
        positive values
        the sizes of the population in the past
        # TODO : should have the same length as time 
    times : list of int
        times at which the population size changes
        # TODO : should have the same length as t past_sizes

    Returns
    -------
    a list object containing the different growth rates of a given population
    at each given time period
    
    Examples
    -------
    >>> sizes2rates(2, [2, 4, 2, 5], [10, 20, 30, 40])
    [0.0, 0.06931471805599453, -0.06931471805599453, 0.0916290731874155]
    >>> print("test")
    "no"
    """
    stime = 0
    sprev = init_size
    rates = []
    for i in range(len(past_sizes)):
        rates.append(math.log(past_sizes[i]/sprev)/(times[i] - stime))
        sprev = past_sizes[i]
        stime = times[i]
    return rates


def migration_matrix(npop, m = 0):
    """
    Set up the migration matrix between sub-populations.
    
    Only supports symmetric island model for the moment, 
    migration rate is identical between all pairs of subpopulations. 
    
    Parameters
    ----------
    npop: int
        Number of sub-populations. Should be at least 2.
    m: float
        overall symetric migration rate. Default is 0, maximum is 1.
        Is the percentage of the population to be immigrants.

    Returns
    -------
    Given N populations, an NxN numpy array of between-subpopulation 
    migration rates. 
    
    Examples
    --------
    >>> migration_matrix(npop=2, m=0.5)
    array([[0.  , 0.25],
           [0.25, 0.  ]])
    >>> migration_matrix(npop=2)
    array([[0., 0.],
           [0., 0.]])
    >>> migration_matrix(npop=1)
    Traceback (most recent call last):
      ...
    SystemExit: npop must be an integer, with a minimum value of 2
    >>> migration_matrix(npop=3, m=1.2)
    Traceback (most recent call last):
      ...
    SystemExit: migr must be a float between 0 and 1 (both included
    """
    
    if not isinstance(npop, int) or npop < 2 :
        sys.exit('npop must be an integer, with a minimum value of 2')
    
    if not isinstance(m, (int,float)) or m < 0 or m > 1 :
        sys.exit('migr must be a float between 0 and 1 (both included')
        
    m = m / (2 * (npop - 1))

    # symmetric island model (later - implement other types of models)
    migration_matrix = np.ones((npop,npop))*m
    np.fill_diagonal(migration_matrix, 0)

    return migration_matrix


def mass_migrations(times, sources, destinations, migr = 1):
    """
    Include mass migration events.

    Parameters
    ----------
    times: int or list
        One or more dates (in generation time) when 
        mass migrations have occured
    sources : int or list
        The source population(s) from the different 
        mass migration events
    destinations : int or list
        The destination population(s) from the different 
        mass migration events
    migr : int or list
        The probability(ies) of immigrating from 
        the source population(s) to 
        the destination populations(s). 
        Default is 1 (merge/split populations)

    Returns
    -------
    a list object that can be passed into 
    msprime.simulate to introduce mass
    migration events

    Examples
    -------
    >>> print("test")
    "no"
    """
    # TODO : idiotproof
    # TODO : examples

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

if __name__ == "__main__":
        import doctest
        doctest.testmod()
