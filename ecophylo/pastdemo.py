# -*- coding: utf-8 -*-
"""


Created on Wed May 13 11:42:50 2020

@author : Maxime Jaunatre <maxime.jaunatre@yahoo.fr>
@author : Elizabeth Bathelemy <barthelemy.elizabeth@gmail.com>

Functions :
    timeframes

"""
# TODO : more info in pastdemo

import msprime
import numpy as np


# TODO : make an example in rmd about timeframes
def timeframes(I, time, a):
    """
    Compute window frames.

    cf: Boitard et al 2016 Plos Genetics

    Parameters
    ----------
    I : int
        Number of time windows. Minimum is 1.
    time : int
        Maximum time in generation time. Must be in R+.
    a : float
        Resolution focus on a particular period. Must be in R+.

    Returns
    -------
    a list of dates (in generation time) corresponding to the
    different time windows

    Examples
    --------
    >>> timeframes(I=2, time=2, a=0.3)
    [0.8830368802245059, 2.0000000000000004]
    
    >>> timeframes(I=3, time=2, a=0.3)
    [0.5653569842838218, 1.226602524471192, 2.0000000000000004]

    >>> timeframes(I=3, time=0.5, a=0.3)
    [0.15896517723882417, 0.32551133275002, 0.4999999999999997]
    """
    # TODO : change for new comsize
    # Idiotproof
    if not isinstance(I, int) or I <= 0 : # TODO : accept and change to float
        raise ValueError('I number of time windows must be an integer superior to 0.')
    if not isinstance(time, (int,float)) or time <= 0 :
        raise ValueError('time maximum time in generation time must be a strict positive float.')
    if not isinstance(a, (int, float)) or  a <= 0 :
        raise ValueError("The resolution a must be a float superior to 0.")
    
    I += 1
    times = [(np.exp((np.log(1+a*time)*i)/(I-1))-1)/a for i in range(1, I)]
    
    return(times)


if __name__ == "__main__":
        import doctest
        doctest.testmod()
