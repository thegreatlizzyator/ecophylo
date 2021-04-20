# -*- coding: utf-8 -*-
"""
islmodel_mugger.py file 

Temporary script for creating migration matrix configuration functions

Created on Wed May 13 11:42:50 2020

@author: barthele
@Author : Maxime Jaunatre <maxime.jaunatre@yahoo.fr>

Functions :
    migration_matrix

"""
# Dependencies
import numpy as np

def migration_matrix(npop, m = 0, migrtime = None:
  """
  """
    if isinstance(m, (int,float)):
      # if m is a number, build single migration matrix of size npop x npop
    
    if isinstance(m, list):
        m = np.array(m)
        dim = m.shape
        if np.sum(m) == 0 :
            #ERROR: migration matrices cannot all be empty
        if len(dim) == 2:
        # when a matrix is provided by user but doesn't change in time
            if dim[0]=!dim[1] or dim[0] =! npop or dim[1] =! npop:
                # ERROR: if a single migration matrix is provided, it should be 
                # symetric and be of size npop
            else:
                diag = np.diagonal(m)
                if np.all(diag != diag[0]):
                    #ERROR: diagonal elements of migration matrix must be 0 
        if len(dim) == 3:
        # when a list of matixes are provided
              if dim[0] =! migrtime : 
              # ERROR : if a list of migration matrixes is provided, there should be
              # as many matrixes as there are migrtime
    else :
        # ERROR
        
  
  
  
  

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
