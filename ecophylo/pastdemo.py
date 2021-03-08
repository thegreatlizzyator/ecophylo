# -*- coding: utf-8 -*-
"""


Created on Wed May 13 11:42:50 2020

@author: barthele

Functions :
    timeframes
    demographic_events

"""
# TODO : more info in pastdemo

import msprime
import numpy as np
import sys


# TODO : make an example in rmd about timeframes
def timeframes(I, T, a):
    """
    Compute window frames.

    cf: Boitard et al 2016 Plos Genetics

    Parameters
    ----------
    I : int
        Number of time windows. Minimum is 1.
    T : int
        Maximum time in generation time. Must be in R+.
    a : float
        Resolution focus on a particular period. Must be in R+.

    Returns
    -------
    a list of dates (in generation time) corresponding to the
    different time windows

    Examples
    --------
    >>> timeframes(I=2, T=2, a=0.3)
    [0.8830368802245059, 2.0000000000000004]
    
    >>> timeframes(I=3, T=2, a=0.3)
    [0.5653569842838218, 1.226602524471192, 2.0000000000000004]

    >>> timeframes(I=3, T=0.5, a=0.3)
    [0.15896517723882417, 0.32551133275002, 0.4999999999999997]

    >>> timeframes(I=-2, T=2, a=0.3)
    Traceback (most recent call last):
      ...
    SystemExit: I number of time windows must be an integer superior to 0.
    
    >>> timeframes(I=2, T=-2, a=0.3)
    Traceback (most recent call last):
      ...
    SystemExit: T maximum time in generation time must be a strict positive float.
    
    >>> timeframes(I=3, T=0, a=0.3)
    Traceback (most recent call last):
      ...
    SystemExit: T maximum time in generation time must be a strict positive float.
    
    >>> timeframes(I=2, T=2, a=-0.3)
    Traceback (most recent call last):
      ...
    SystemExit: The resolution a must be a float superior to 0.
    """
    # Idiotproof
    if not isinstance(I, int) or I <= 0 :
        sys.exit('I number of time windows must be an integer superior to 0.')
    if not isinstance(T, (int,float)) or T <= 0 :
        sys.exit('T maximum time in generation time must be a strict positive float.')
    if not isinstance(a, (int, float)) or  a <= 0 :
        sys.exit("The resolution a must be a float superior to 0.")
    
    I += 1
    times = [(np.exp((np.log(1+a*T)*i)/(I-1))-1)/a for i in range(1, I)]
    
    return(times)


def demographic_events(changetime, past_sizes):
    """
    Change the demography of populations over given time periods.

    Parameters
    ----------
    changetime : list of int
        When the demographic changes have occured. Must be >= 0.
    past_sizes : list of int
        Population past_sizes at the different time periods. Must be > 0.

    Returns
    -------
    a list object to be passed into msprime.simulate to change 
    the demography over given times periods. 
    Changes affect Population past_sizes at the different times, 
    growth rates are automatically set to 0 (constant size).
    Should later implement an option to change a specific population 
    for now the changes affect all populations simultaneously.

    Examples
    --------
    >>> demographic_events([1, 2], [42, 9000])
    [{'type': 'population_parameters_change', 'time': 1, 'growth_rate': None, 'initial_size': 42, 'population': -1}, {'type': 'population_parameters_change', 'time': 2, 'growth_rate': None, 'initial_size': 9000, 'population': -1}]
    
    >>> demographic_events([1], [42])
    [{'type': 'population_parameters_change', 'time': 1, 'growth_rate': None, 'initial_size': 42, 'population': -1}]
    
    >>> demographic_events([1, 2], [42])
    Traceback (most recent call last):
      ...
    SystemExit: changetime and past_sizes list must be of same length
    
    >>> demographic_events([1, '2'], [42, 9000])
    Traceback (most recent call last):
      ...
    SystemExit: changetime must be a list of int
    >>> demographic_events([1, 2], [42, '9000'])
    Traceback (most recent call last):
      ...
    SystemExit: past_sizes must be a list of int
    
    >>> demographic_events([1, -2], [42, 9000])
    Traceback (most recent call last):
      ...
    SystemExit: changetime must be strict positive values
    >>> demographic_events([1, 2], [42, -9000])
    Traceback (most recent call last):
      ...
    SystemExit: past_sizes must be strict positive values
    """
    # Idiot proof
    if len(changetime) != len(past_sizes) :
        sys.exit('changetime and past_sizes list must be of same length')
    if not all(isinstance(x, int) for x in changetime) :
        sys.exit('changetime must be a list of int')
    if not all(isinstance(x, int) for x in past_sizes) :
        sys.exit('past_sizes must be a list of int')
    if not all((x > 0) for x in changetime) :
        sys.exit('changetime must be strict positive values')
    if not all((x > 0) for x in past_sizes) :
        sys.exit('past_sizes must be strict positive values')

    dc = [msprime.PopulationParametersChange(time=t, initial_size =s) for t, s in zip(changetime, past_sizes)]
    return dc

if __name__ == "__main__":
        import doctest
        doctest.testmod()
