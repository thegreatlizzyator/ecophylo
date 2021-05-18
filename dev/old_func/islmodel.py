# -*- coding: utf-8 -*-
"""
islmodel.py file 

Module with functions creating msprime (dependence) objects to be used 
by the ecophylo.simulate function.

Created on Wed May 13 11:42:50 2020

@author: barthele
@Author : Maxime Jaunatre <maxime.jaunatre@yahoo.fr>

Functions :
    population_configurations
    sizes2rates
    mergesizes2rates
    migration_matrix
    mass_migrations

"""
# Dependencies

import msprime
import numpy as np
import sys
import math


def population_configurations(init_sizes, samples, rates = None) :
    """
    Set up the initial population configurations.
    Parameters
    ----------
    init_sizes : list of int
        positive values
        the initial size of the population
    samples : int or list of int
        Positive values.
        The number of individuals simulated in each population. 
        Must be a list with as many elements as there are populations. If only one
        value is provided, all populations are assumed to have the same sample size
        and the value is duplicated npop times. 
    rates = None : list of float
        The initial community growth rates.
        Must be a list with as many elements as populations containing the initial
        growth rates.
  
    Returns
    -------
    a list object that can be passed into 
    msprime.simulate to indicate initial 
    population configurations
    """
    # TODO : examples
    
    # Idiot proof
    # check init_sizes
    if not isinstance(init_sizes, list):
        if isinstance(init_sizes, (int,float)) : 
            init_sizes = [int(init_sizes)]
        else :
            sys.exit('init_sizes must be an int or a list of int. '+
                     'Float will be rounded to int')
    if not all(isinstance(x, (int, float)) for x in init_sizes) :
        sys.exit('init_sizes must be a list of int')
    if not all(x > 0 for x in init_sizes) :
        sys.exit('init_sizes must be positive values')
    else : 
        for x in range(len(init_sizes)) : 
            init_sizes[x] = int(init_sizes[x])
    # check samples
    if not isinstance(samples, list):
        if isinstance(samples, (int,float)) : 
            samples = [int(samples)]
        else :
            sys.exit('samples must be an int or a list of int. '+
                     'Float will be rounded to int')
    if not all(isinstance(x, (int, float)) for x in samples) :
        sys.exit('samples must be a list of int')
    if not all(x > 0 for x in samples) :
        sys.exit('samples must be positive values')
    else : 
        for x in range(len(samples)) : 
            samples[x] = int(samples[x])
    if len(samples) == 1:
        samples = samples*len(init_sizes)
    # check rates
    if rates != None :
        if not isinstance(rates, list) or not all(isinstance(y, (float, int)) for y in rates):
            sys.exit("rates must be a list of floats")
    else : 
        rates = [0]*len(init_sizes)
    # check lenghts
    if len(init_sizes) != len(samples) or ( rates is not None and len(init_sizes) != len(rates) ):
        sys.exit('init_sizes, samples and rates must have the same'+
                 ' number of elements')

    pc = [msprime.PopulationConfiguration(sample_size = s, initial_size = i, growth_rate = g) for s, i, g in zip(samples, init_sizes, rates)]
    return pc

# TODO : stripe is the second mogwai born of Gizmo, be gentle with him
def population_configurations_stripe(init_sizes, past_sizes, changetime, samples, stable_pop = True, rates = None, demo = None):
    """
    Set up the initial population configurations and past demographic events.

    Parameters
    ----------
    init_sizes : list of int
        positive values
        the initial size of the population
    past_sizes : list of list of int
        positive values
        The past community sizes. 
        Must be a list with as many elements as populations each containing a list
        of past sizes.
    changetime : list of list of int
        The times at which populations have changed sizes.
        Must be a list with as many elements as populations each containing a list
        of times.
    samples : int or list of int
        Positive values.
        The number of individuals simulated in each population. 
        Must be a list with as many elements as there are populations. If only one
        value is provided, all populations are assumed to have the same sample size
        and the value is duplicated npop times. 
    stable_pop = True : bool
        Whether or not changes in population sizes are instantenous or progressive.
        This parameter is deprecated if growth rates are provided.
    rates = None : list of float
        The initial community growth rates.
        Must be a list with as many elements as populations containing the initial
        growth rates.
     
    Notes
    -----
        All parameters must have the same length, with an exception for stable_pop.
        If a population change its size, you need to indicate at least one population 
        change for other ones. For example a constant population should be coded with 
        at least a past_size = init_size at changetime = 0.

    Returns
    -------
        two list objects that can be passed into msprime.simulate to indicate
        initial community states and past demographic changes.
    
    Examples
    -------
    >>> import msprime
    >>> demo = population_configurations_stripe(init_sizes= 500, past_sizes= None, changetime= None, samples = 10 )
    >>> demo = population_configurations_stripe(init_sizes= 500, past_sizes= [[500]], changetime= [[0]], samples = 10 )
    >>> demo
    Demography(populations=[Population(initial_size=500, growth_rate=0, name='pop_0', description='', extra_metadata={}, default_sampling_time=None, initially_active=None, id=0)], events=[], migration_matrix=array([[0.]]))
    >>> demo = population_configurations_stripe(init_sizes= [500], past_sizes= [[1000]], changetime=[[100]], samples = [10] )
    >>> demo
    Demography(populations=[Population(initial_size=500, growth_rate=0, name='pop_0', description='', extra_metadata={}, default_sampling_time=None, initially_active=None, id=0)], events=[PopulationParametersChange(time=100, initial_size=1000, growth_rate=0, population=0)], migration_matrix=array([[0.]]))
    >>> demo = population_configurations_stripe(init_sizes= [500], past_sizes= [[1000]], changetime=[[100]], stable_pop = False, samples = [10] )
    >>> demo
    Demography(populations=[Population(initial_size=500, growth_rate=0.007, name='pop_0', description='', extra_metadata={}, default_sampling_time=None, initially_active=None, id=0)], events=[PopulationParametersChange(time=100, initial_size=1000, growth_rate=0.007, population=0)], migration_matrix=array([[0.]]))
    >>> demo = population_configurations_stripe(init_sizes= [500, 300], past_sizes= [[1000], [500]], changetime=[[100], [40]], samples = [10] )
    >>> demo
    Demography(populations=[Population(initial_size=500, growth_rate=0, name='pop_0', description='', extra_metadata={}, default_sampling_time=None, initially_active=None, id=0), Population(initial_size=300, growth_rate=0, name='pop_1', description='', extra_metadata={}, default_sampling_time=None, initially_active=None, id=1)], events=[PopulationParametersChange(time=40, initial_size=500, growth_rate=0, population=0), PopulationParametersChange(time=40, initial_size=500, growth_rate=0, population=1), PopulationParametersChange(time=100, initial_size=1000, growth_rate=0, population=0), PopulationParametersChange(time=100, initial_size=500, growth_rate=0, population=1)], migration_matrix=array([[0., 0.],
           [0., 0.]]))
    >>> #print(demo.debug(), file=sys.stderr)
    """
    # use a simplier function !
    if past_sizes is None or changetime is None:
        pc = population_configurations(init_sizes, samples, rates)
        return pc, None

    # Idiot proof
    # check init_sizes
    if not isinstance(init_sizes, list):
        if isinstance(init_sizes, (int,float)) : 
            init_sizes = [int(init_sizes)]
        else :
            sys.exit('init_sizes must be an int or a list of int. '+
                     'Float will be rounded to int')
    if not all(isinstance(x, (int, float)) for x in init_sizes) :
         sys.exit('init_sizes must be a list of int')
    if not all(x > 0 for x in init_sizes) :
        sys.exit('init_sizes must be positive values')
    else : 
        for x in range(len(init_sizes)) : 
            init_sizes[x] = int(init_sizes[x])
    # check past_sizes
    if not isinstance(past_sizes, list):
        if not isinstance(past_sizes, (int,float)) : 
            sys.exit('past_sizes must be int, list of int or'+
            ' nested list of int')
    else :
        for x in past_sizes:
            if isinstance(x, list):
                if not all(isinstance(y, (float, int)) for y in x) : 
                    sys.exit('past_sizes must be int, list of int or'+
                             ' nested list of int')
            else :
                if not isinstance(x, (float, int)):
                    sys.exit('past_sizes must be int, list of int or'+
                    ' nested list of int')
                past_sizes = [past_sizes]
    # check changetime
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
                if any(y < 0 for y in x) : 
                    sys.exit('changetime must be positive values')
            else :
                if not isinstance(x, (float, int)):
                    sys.exit('changetime must be int, list of int or'+
                             ' nested list of int')
                if x < 0:
                    sys.exit('changetime must be positive values')
                changetime = [changetime]
    # stable_pop
    if not isinstance(stable_pop, bool):
        sys.exit('stable_pop must be a boolean')
    # check rates
    if rates != None :
        if not isinstance(rates, list) or not all(isinstance(y, (float, int)) for y in rates):
            sys.exit("rates must be a list of floats")
    # check samples
    if not isinstance(samples, list):
        if isinstance(samples, (int,float)) : 
            samples = [int(samples)]
        else :
            sys.exit('samples must be an int or a list of int. '+
                     'Float will be rounded to int')
    if not all(isinstance(x, (int, float)) for x in samples) :
        sys.exit('samples must be a list of int')
    if not all(x > 0 for x in samples) :
        sys.exit('samples must be positive values')
    else : 
        for x in range(len(samples)) : 
            samples[x] = int(samples[x])
    # check lenghts
    if len(past_sizes) != len(changetime) or len(changetime) != len(init_sizes):
        sys.exit('past_sizes, changetime, rates and init_sizes must have the same'+
                 ' number of elements')
    if rates != None and len(past_sizes) != len(rates):
        sys.exit('past_sizes, changetime, rates and init_sizes must have the same'+
                 ' number of elements')
    if len(samples) == 1:
        samples = samples*len(init_sizes)
    
    # need to check :
    # sample is uniq, then duplicate
    # same with stable_pop
    
    # sizes = list()
    # changetimen = list()
    # tmp_past_sizes = list()
    # tmp_init_sizes = list()
    # tmp_changetime = list()
    # for i in range(len(past_sizes)):
    #     sizes.append( [init_sizes[i]] + past_sizes[i] )
    #     changetimen.append([0] + changetime[i])
    
    # check insiders lengths
    for i in range(len(past_sizes)):
        if len(past_sizes[i]) != len(changetime[i]):
            sys.exit('each element of past_sizes and changetime lists must be '
              +'of same lenghts')
        # # extract init_values
        # tmp_init_sizes.append(sizes[i][0])
        # if len(sizes[i]) == 1:
        #     tmp_past_sizes.append(sizes[i][0])
        #     tmp_changetime.append([0])
        # else:
        #     tmp_past_sizes.append(sizes[i][1:])
        #     tmp_changetime.append(changetimen[i][1:])
    
    
    if rates != None :
        stable_pop = False
    if sum([sum(i) for i in changetime]) == 0 : 
        stable_pop = True
    npop = len(init_sizes)
    
    if demo is None:
        demo = msprime.Demography()
    sizes = mergesizes2rates(past_sizes, changetime, init_sizes, False)
    times = sizes[0]
    
    # working on rates
    if stable_pop: # no rates
        rates = [ [0]*len(times) for _ in range(npop)]
        rates = [times] + rates
    
    elif rates != None : # if rates are provided, fixed rates for all the epochs
        if len(rates) == npop and all(isinstance(x, (int,float)) for x in rates):
            tmp_rates = [times]
            for i in range(npop):
                tmp_rates = tmp_rates + [[rates[i]]*len(times)]
            rates = tmp_rates
    #     else : # not homogenous rates in time...too difficult to code right now
    #         rates = eco.mergesizes2rates(rates, changetime, init_sizes, True)
    else : # no rates provided so computed them from past_sizes and changetime
        rates = mergesizes2rates(past_sizes, changetime, init_sizes, True)
    
    for i in range(npop): # initiate every pop
        demo.add_population(name="pop_"+str(i), 
                            initial_size=init_sizes[i], 
                            growth_rate=rates[i+1][0])
        
        rates[i+1].append(rates[i+1][-1]) # duplicate late rate for infinity
        rates[i + 1].remove(rates[i+1][0]) # remove first rate
    
    for i in range(len(times)): # later populations parameter changes
        for ii in range(npop):
            if times[i] == 0:
                continue
            demo.add_population_parameters_change(
                time = times[i], initial_size = sizes[ii+1][i], 
                growth_rate = rates[ii+1][i], population = ii
            )
    
    return demo


def sizes2rates(sizes, changetime):
    """
    Compute the growth rates corresponding to a set of past sizes at different 
    times for a single population with a given initial size

    Parameters
    ----------
    sizes : list of int
        positive values
        the sizes of the population in the past
        should have the same length as time and
        first one is the initial size of the population
    changetime : list of int
        times at which the population size changes
        should have the same length as t past_sizes
        First one should be 0.

    Returns
    -------
    a list object containing the different growth rates of a given population
    at each given time period
    
    Examples
    -------
    >>> sizes2rates([2, 2, 4, 2, 5], [0,10, 20, 30, 40])
    [0.0, 0.069, -0.069, 0.092]
    """
    
    # Idiot proof
    if len(changetime) != len(sizes) :
        sys.exit('changetime and sizes list must be of same length')
    if len(set(changetime)) != len(changetime) :
        sys.exit('Duplicated times in changetime are not possible')
    
    if not all(isinstance(x, (int, float)) for x in sizes) :
        sys.exit('sizes must be a list of int')
    if not all(isinstance(x, (int, float)) for x in changetime) :
        sys.exit('changetime must be a list of int')
    if not all((x >= 0) for x in changetime) :
        sys.exit('changetime must be strict positive values')
    if not all((x > 0) for x in sizes) :
        sys.exit('sizes must be strict positive values')

    sprev = sizes[0]
    stime = changetime[0]
    rates = []
    for i in range(1,len(sizes)):
        val = math.log(sizes[i]/sprev)/(changetime[i] - stime)
        rates.append(round(val, 3))
        sprev = sizes[i]
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
    changetime : list of int or list of list of int
        None by default
        times at which the population size changes
    init_size : list of int
        positive values
        the initial size of the population
    sizesvalues : bool
        indicate if the past_values list is a matrix of sizes or growth rates.
        will be overridden if any past_values is a float below one.
        
    Notes
    -----
        sizes in float will be changed to int.
        past_values, changetime and init_sizes must have the same lengths
        every element of past_values and changetime need to have the same length
        ex : [3, 2] & [3, 2] where the first element is past_values and the second
        is changetime. The first element of past_values is of size 3 and 
        changetime too.

    Returns
    -------
    a list object containing the different growth rates of a given population
    at each given time period
    
    Examples
    -------
    >>> mergesizes2rates([[1000], [500]], [[100], [40]], [500, 300], False)
    [[40, 100], [500, 1000], [500, 500]]
    >>> mergesizes2rates([[3, 5, 3, 4]], [[10, 20, 30, 40]], [2], True)
    [[10, 20, 30, 40], [0.041, 0.051, -0.051, 0.029]]
    >>> mergesizes2rates([[3, 5, 3, 4]], [[10, 20, 30, 40]], 2, True)
    [[10, 20, 30, 40], [0.041, 0.051, -0.051, 0.029]]
    
    >>> mergesizes2rates([[0.2, 0.1, -0.09, 0.009]], [[10, 20, 30, 40]], [2], False)
    [[10, 20, 30, 40], [0.2, 0.1, -0.09, 0.009]]
    
    >>> init_size = [2, 3]
    >>> past_values = [[4, 3, 10, 7], [1, 2, 5, 2,3]]
    >>> changetime = [[10, 20, 30, 40], [10, 15, 21, 60, 70]]
    >>> mergesizes2rates(past_values, changetime, init_size, True)
    [[10, 15, 20, 21, 30, 40, 60, 70], [0.069, -0.029, -0.029, 0.12, 0.12, -0.036, -0.036, -0.036], [-0.11, 0.139, 0.153, 0.153, -0.023, -0.023, -0.023, 0.041]]
    >>> mergesizes2rates(past_values, changetime, init_size, False)
    [[10, 15, 20, 21, 30, 40, 60, 70], [4, 4, 3, 3, 10, 7, 7, 7], [1, 2, 2, 5, 5, 5, 2, 3]]
    
    >>> init_size = [2, 3]
    >>> past_values = [[0.2, 0.1, -0.09, 0.009], [0.42, -0.4, 0.07, 0.42,0.01]]
    >>> changetime = [[10, 20, 30, 40], [10, 15, 21, 60, 70]]
    >>> mergesizes2rates(past_values, changetime, init_size, True)
    [[10, 15, 20, 21, 30, 40, 60, 70], [0.2, 0.1, 0.1, -0.09, -0.09, 0.009, 0.009, 0.009], [0.42, -0.4, 0.07, 0.07, 0.42, 0.42, 0.42, 0.01]]
    
    """
    # Idiot proof
    ## Init_size 
    if not isinstance(init_size, list):
        if isinstance(init_size, (int,float)) : 
            init_size = [int(init_size)]
        else :
            sys.exit('init_size must be an int or a list of int. '+
                     'Float will be rounded to int')
    if not all(isinstance(x, (int, float)) for x in init_size) :
        sys.exit('init_size must be a list of int')
    if not all(x > 0 for x in init_size) :
        sys.exit('init_size must be positive values')
    else : 
        for x in range(len(init_size)) : 
            init_size[x] = int(init_size[x])
    # changetime : list of int or list of list of int
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
                if any(y < 0 for y in x) : 
                    sys.exit('changetime must be positive values')
            else :
                if not isinstance(x, (float, int)):
                    sys.exit('changetime must be int, list of int or'+
                             ' nested list of int')
                if x < 0:
                    sys.exit('changetime must be positive values')
                changetime = [changetime]
    # past_values
    if not isinstance(past_values, list):
        if not isinstance(past_values, (int,float)) : 
            sys.exit('past_values must be int, list of int or'+
                     ' nested list of int')
    else :
        for x in past_values:
            if isinstance(x, list):
                if not all(isinstance(y, (float, int)) for y in x) : 
                    sys.exit('past_values must be int, list of int or'+
                             ' nested list of int')
            else :
                if not isinstance(x, (float, int)):
                    sys.exit('past_values must be int, list of int or'+
                             ' nested list of int')
                past_values = [past_values]
    # size_values
    if not isinstance(sizevalues, bool):
        sys.exit('sizevalues must be a boolean')
    # check lengths
    if len(past_values) != len(changetime) or len(changetime) != len(init_size):
        sys.exit('past_values, changetime and init_sizes must have the same'+
        ' number of elements')
    
    for i in range(len(past_values)):
        if len(past_values[i]) != len(changetime[i]):
            sys.exit('each element of past_values and changetime lists must be '
            +'of same lenghts')
    
    is_rate = False
    for i in range(len(past_values)):
        for ii in past_values[i]:
            if ii  < 1 :
                is_rate = True
                break
    first_size = list(init_size)
    
    npop = len(past_values)
    # aggregate the times.
    nrow = changetime[0]
    for i in range(1,len(changetime)):
        nrow = sorted(nrow + list(set(changetime[i]) - set(nrow) ))
    
    matrix = [nrow]
    for i in range(npop) :
        tmpsize = list(past_values[i])
        tmptime = list(changetime[i])
        init_time = 0
        res = [] # init empty list
        for ii in nrow:
            if len(tmptime) == 0 : # if no time left for this pop
                res.append('NA') 
                continue
            
            if ii == tmptime[0]: # when time is defined for a pop
                val = tmpsize.pop(0)
                tmptime.remove(tmptime[0])
                if sizevalues and not is_rate : # compute directly the growrate
                    tmp = val
                    # print("sizes :",[first_size[i], val])
                    # print("changetime :", [init_time,ii])
                    val = -sizes2rates(sizes = [first_size[i], val], 
                                      changetime = [init_time,ii])[0]
                    first_size[i] = tmp
                    init_time = ii

                res.append(val) # add the value
                # remove previous NA (only if its grates)
                place = len(res)-1
                while res[place-1] == 'NA' and place >= 0 :
                    place-=1

                if not sizevalues : # if it is sizes, get last sizes and apply it
                    if place == 0:
                        for iii in range(len(res)-1):
                            res[iii] = first_size[i]
                    else :  
                        for iii in range(place, len(res)-1):
                            res[iii] = res[place - 1]
                else: # apply NA with lateest rates
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


def migration_configuration(npop, migr = 1, migr_time = None):
    if isinstance(migr, (int, float)) :
        migration = islmodel.migration_matrix(npop = npop, migr = migr)
    else : 
        migration = [[0., 0.5], [0.5, 0.]]


def migration_matrix(npop, migr = 2):
    """
    Set up the migration matrix between sub-populations.
    
    Only supports symmetric island model for the moment, 
    migration rate is identical between all pairs of subpopulations. 
    
    Parameters
    ----------
    npop: int
        Number of sub-populations. Should be at least 2.
    migr: float
        overall symetric migration rate. Default is 0, maximum is npop.
        Is the percentage of the population to be immigrants.
    
    Returns
    -------
    Given N populations, an NxN numpy array of between-subpopulation 
    migration rates. 
    
    Examples
    --------
    >>> migration_matrix(npop=2, migr=0.5)
    array([[0.  , 0.25],
           [0.25, 0.  ]])
    >>> migration_matrix(npop=2)
    array([[0., 1.],
           [1., 0.]])
    >>> migration_matrix(npop=1)
    Traceback (most recent call last):
      ...
    SystemExit: npop must be an integer, with a minimum value of 2
    >>> migration_matrix(npop=3, migr=3.2)
    Traceback (most recent call last):
      ...
    SystemExit: migr must be a float between 0 and npop (both included
    """
    
    if not isinstance(npop, int) or npop < 2 :
        sys.exit('npop must be an integer, with a minimum value of 2')
    
    if not isinstance(migr, (int,float)) or migr < 0 or migr > npop :
        sys.exit('migr must be a float between 0 and npop (both included')
        
    migr = migr / (2 * (npop - 1))
    
    # symmetric island model (later - implement other types of models)
    migration_matrix = np.ones((npop,npop))*migr
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
    # pc, config = population_configurations_stripe(init_sizes= [500, 300], past_sizes= [[1000], [500]], changetime=[[100], [0]], samples = [10] )
    # population_configurations_stripe(init_sizes= [500], past_sizes= [[1000]], changetime=[[100]], samples = [10] )
