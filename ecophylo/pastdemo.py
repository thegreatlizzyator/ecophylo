# -*- coding: utf-8 -*-
"""
 ALL ECOPHYLO


Created on Wed May 13 11:42:50 2020

@author: barthele
"""
import msprime
import numpy as np

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

