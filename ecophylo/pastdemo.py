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


# TODO : figure out what to do with timeframes function
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

    Examples
    --------
    >>> timeframes(I=2, T=2, a=0.3)
    [0.8830368802245059, 2.0000000000000004]
    
    >>> timeframes(I=3, T=2, a=0.3)
    [0.5653569842838218, 1.226602524471192, 2.0000000000000004]

    >>> timeframes(I=3, T=0, a=0.3)
    [0.0, 0.0, 0.0]

    >>> timeframes(I=3, T=0.5, a=0.3)
    [0.15896517723882417, 0.32551133275002, 0.4999999999999997]

    >>> timeframes(I=3, T=-0.7, a=0.3)
    [-0.2518881782078394, -0.48474206011941934, -0.7]
    """
    # TODO : idiotproof the T value
    if not isinstance(I, int) or I <= 0 :
        sys.exit('The number of time windows must be an integer superior to 0')
    if not isinstance(T, (int,float)) :
        sys.exit('The maximum time in generation time must be a float')
    if not isinstance(a, (int, float)) or  a <= 0 :
        sys.exit("The resolution a must be a float superior to 0.")
    
    I += 1
    times = [(np.exp((np.log(1+a*T)*i)/(I-1))-1)/a for i in range(1, I)]
    
    return(times)


def demographic_events(epochs, sizes):
    """
    Change the demography of populations over given time periods.

    Parameters
    ----------
    epochs: list of int
        When the demographic changes have occured.
    sizes: list of int
        Population sizes at the different time periods.

    Returns
    -------
    a list object to be passed into msprime.simulate to change 
    the demography over given times periods. 
    Changes affect Population sizes at the different times, 
    growth rates are automatically set to 0 (constant size).
    Should later implement an option to change a specific population 
    for now the changes affect all populations simultaneously.

    Examples
    --------
    >>> demographic_events([1,2], [42, 9000])
    [{'type': 'population_parameters_change', 'time': 1, 'growth_rate': None, 'initial_size': 42, 'population': -1}, {'type': 'population_parameters_change', 'time': 2, 'growth_rate': None, 'initial_size': 9000, 'population': -1}]
    """
    # TODO : test if input is float
    if not isinstance(epochs, list) :
        sys.exit('epochs must be a list')
    if not isinstance(sizes, list) :
        sys.exit('sizes must be a list')
    if len(epochs) != len(sizes) :
        sys.exit('epochs and sizes list must be of same length')

    dc = [msprime.PopulationParametersChange(time=t, initial_size =s) for t, s in zip(epochs, sizes)]
    return dc

if __name__ == "__main__":
        import doctest
        doctest.testmod()
