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

def sizes2rates(init_size, past_sizes, changetime):
    """
    Compute the growth rates corresponding to a set of past sizes at different 
    times for a single population with a given initial size

    Parameters
    ----------
    init_size : int
        positive value
        the initial size of the population
    past_sizes : list of int
        positive values
        the sizes of the population in the past
        should have the same length as time 
    changetime : list of int
        times at which the population size changes
        should have the same length as t past_sizes

    Returns
    -------
    a list object containing the different growth rates of a given population
    at each given time period
    
    Examples
    -------
    >>> sizes2rates(2, [2, 4, 2, 5], [10, 20, 30, 40])
    [0.0, 0.069, -0.069, 0.092]
    """
    
    # Idiot proof
    if len(changetime) != len(past_sizes) :
        sys.exit('changetime and past_sizes list must be of same length')
    if len(set(changetime)) != len(changetime) :
        sys.exit('Duplicated times in changetime are not possible')
    
    if not all(isinstance(x, int) for x in past_sizes) :
        sys.exit('past_sizes must be a list of int')
    if not all(isinstance(x, int) for x in changetime) :
        sys.exit('changetime must be a list of int')
    if not all((x > 0) for x in changetime) :
        sys.exit('changetime must be strict positive values')
    if not all((x > 0) for x in past_sizes) :
        sys.exit('past_sizes must be strict positive values')
        
    if not isinstance(init_size, int) or init_size <= 0 :
        sys.exit('init_size must be a single strict positive value')
    
    stime = 0
    sprev = init_size
    rates = []
    for i in range(len(past_sizes)):
        val = math.log(past_sizes[i]/sprev)/(changetime[i] - stime)
        rates.append(round(val, 3))
        sprev = past_sizes[i]
        stime = changetime[i]
    return rates


def mergesizes2rates(past_values, changetime, init_size = None , sizevalues = True):
    """
    Compute the growth rates corresponding to a set of past sizes at different 
    times for a single population with a given initial size

    Parameters
    ----------
    past_values : list of int or list of list of int
        positive values
        the sizes of the population in the past
        # TODO : should have the same length as time
    init_size : list of int
        positive value
        the initial size of the population
        # TODO : same len as past_value matrix
    changetime : list of int
        None by default
        times at which the population size changes
        # TODO : should have the same length as t past_sizes
    sizesvalues : bool
        indicate if the past_values list is a matrix of sizes or growth rates.

    Returns
    -------
    a list object containing the different growth rates of a given population
    at each given time period
    
    Examples
    -------
    >>> mergesizes2rates([[3, 5, 3, 4]], [[10, 20, 30, 40]], [2], True)
    [[10, 20, 30, 40], [0.041, 0.026, -0.017, 0.007]]
    
    >>> mergesizes2rates([[0.2, 0.1, -0.09, 0.009]], [[10, 20, 30, 40]], [2], False)
    [[10, 20, 30, 40], [0.2, 0.1, -0.09, 0.009]]
    
    >>> init_size = [2, 3]
    >>> past_values = [[4, 3, 10, 7], [1, 2, 5, 2,3]]
    >>> changetime = [[10, 20, 30, 40], [10, 15, 21, 60, 70]]
    >>> mergesizes2rates(past_values, changetime, init_size, True)
    [[10, 15, 20, 21, 30, 40, 60, 70], [0.069, -0.014, -0.014, 0.04, 0.04, -0.009, -0.009, -0.009], [-0.11, 0.046, 0.044, 0.044, -0.015, -0.015, -0.015, 0.006]]
    
    >>> init_size = [2, 3]
    >>> past_values = [[0.2, 0.1, -0.09, 0.009], [0.42, -0.4, 0.07, 0.42,0.01]]
    >>> changetime = [[10, 20, 30, 40], [10, 15, 21, 60, 70]]
    >>> mergesizes2rates(past_values, changetime, init_size, False)
    [[10, 15, 20, 21, 30, 40, 60, 70], [0.2, 0.1, 0.1, -0.09, -0.09, 0.009, 0.009, 0.009], [0.42, -0.4, 0.07, 0.07, 0.42, 0.42, 0.42, 0.01]]
    
    """
    
    # Idiot proof
    #  past_values : list of int or list of list of int
    #     positive values
    #     the sizes of the population in the past
    #     # TODO : should have the same length as time
    # init_size : list of int
    #     positive value
    #     the initial size of the population
    #     # TODO : same len as past_value matrix
    # changetime : list of int
    #     None by default
    #     times at which the population size changes
    #     # TODO : should have the same length as t past_sizes
    # sizesvalues : bool
    #     indicate if the past_values list is a matrix of sizes or growth rates.
    
    if len(changetime) != len(past_values) :
        sys.exit('changetime and past_values list must be of same length')
    # if len(set(changetime)) != len(changetime) :
    #     sys.exit('Duplicated times in changetime are not possible')
    
    # if not all(isinstance(x, int) for x in past_values) :
    #     sys.exit('past_values must be a list of int')
    # if not all(isinstance(x, int) for x in changetime) :
    #     sys.exit('changetime must be a list of int')
    # if not all((x > 0) for x in changetime) :
    #     sys.exit('changetime must be strict positive values')
    # if not all((x > 0) for x in past_values) :
    #     sys.exit('past_values must be strict positive values')
        
    # if not isinstance(init_size, int) or init_size <= 0 :
    #     sys.exit('init_size must be a single strict positive value')
        
    # TODO : check if sizevalue,  if past values are int !
    # TODO : check if changetime[1] same size as pastvalue[1] etc
    
    npop = len(past_values)
    # aggregate the times.
    nrow = changetime[0]
    for i in range(1,len(changetime)):
      nrow = sorted(nrow + list(set(changetime[i]) - set(nrow) ))
    
    matrix = [nrow]
    for i in range(npop) :
      tmpsize = list(past_values[i])
      tmptime = list(changetime[i])
      res = [] # init empty list
      for ii in  nrow:
        
        if len(tmptime) == 0 : # if no time left for this pop
          res.append('NA') 
          continue
        
        if ii == tmptime[0]: # when time is defined for a pop
          val = tmpsize.pop(0)
          tmptime.remove(tmptime[0])
          if sizevalues : # compute directly the growrate
            tmp = val
            val = sizes2rates(init_size[i], 
                                  past_sizes = [val], 
                                  changetime = [ii])[0]
            init_size[i] = tmp
          
          res.append(val) # add the value
          # remove previous NA (only if its grates)
          place = len(res)-1
          while res[place-1] == 'NA' and place >= 0 :
            place-=1
          
          # change NA with values
          for iii in range(place, len(res)) :
            res[iii] = val
        else: # if nothing defined for this pop at this time
          res.append('NA')
      
      if(res[-1] == 'NA'): # if ending with NA, expand last known rate
        place = len(res)-1
        while res[place-1] == 'NA' and place >= 0 :
          place-=1
          # change NA with values
        for iii in range(place, len(res)) :
          res[iii] = res[place-1]
      matrix.append(res) # append this pop list in matrix
    return matrix


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
