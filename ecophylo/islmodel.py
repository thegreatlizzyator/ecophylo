# -*- coding: utf-8 -*-
"""
 ALL ECOPHYLO


Created on Wed May 13 11:42:50 2020

@author: barthele
"""
import msprime
import numpy as np
import sys

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
    a list object that can be passed into 
    msprime.simulate to indicate initial 
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
    >>> migration_matrix(subpops=2, migr=0.5)
    array([[0.  , 0.25],
           [0.25, 0.  ]])
    >>> migration_matrix(subpops=2)
    array([[0., 0.],
           [0., 0.]])
    """
    if subpops < 1:
        sys.exit("there should be at least 2 populations")
    # TODO: incoherence between > 1 and at least 2 statement !
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

if __name__ == "__main__":
        import doctest
        doctest.testmod()
