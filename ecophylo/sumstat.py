# -*- coding: utf-8 -*-
"""

Created on Wed May 13 11:42:50 2020

@author : Maxime Jaunatre <maxime.jaunatre@yahoo.fr>
@author : Elizabeth Bathelemy <barthelemy.elizabeth@gmail.com>

Functions :
    getAbund
    getDeme

"""

import numpy as np

def getAbund(tree, samples = None):
    """
    
    Parameters
    ----------
    tree : (ete3 class)
        Phylogeny with attributes on leafs. This attributes is a character 
        string containing all names of the species individual (mean to use 
        topPhylo result). Names is formated like this :
          " name1 name2 name3"
    samples : int
        number of individual in the community
        # TODO : set this check as optionnal

    Returns
    -------
    sfs : TYPE
        DESCRIPTION.

    Examples
    --------
    >>> from ete3 import Tree
    >>> tree = Tree('(((A:5,(B:3, C:3))1:2,(D:2, E:2)1:5)1:2, (F:3, G:3)1:6);')
    >>> import ecophylo as eco
    >>> phylo = eco.toPhylo(tree, 0.5, seed = 42)
    >>> print(phylo)
    <BLANKLINE>
          /-sp1
       /-|
    --|   \-sp2
      |
       \-sp3
    >>> getAbund(phylo, 7)
    [3, 2, 2]
    """
    # Idiot proof
    if tree.__class__.__name__ != 'TreeNode' :
        raise ValueError('tree must have a class TreeNode')
    if samples != None:
      if not isinstance(samples, int):
          raise ValueError('samples must be an integer')

    sfs = list()
    abund = list()
    for leaf in tree.iter_leaves():
        try:
            inds = leaf.mergedInd.lstrip(" ") # remove 1st space
            inds = list(inds.split(" ")) # strip on spaces
            abund.append(len(inds)) # get length 
        except AttributeError:
            abund.append(1)
    sfs.extend(abund)
    # think about catching error when phylogeny has only 1 sp
    
    if samples != None and sum(sfs) != samples:
        raise Exception(f"Simulated phylogeny has only one species!")
        # TODO : modify error with a better check here
    return sfs

def getDeme(tree, div = False):
    """
    
    Parameters
    ----------
    tree : (ete3 class)
        Phylogeny with attributes on leafs. This attributes is a character 
        string containing all names of the species individual (mean to use 
        topPhylo result). Names is formated like this :
          " name1 name2 name3"
    div : bool
        Option to simplify the matrix to a simple list of Deme species diversity.

    Returns
    -------
    indiv : nest list of int
        site/species matrix of the tree

    Examples
    --------
    >>> from ete3 import Tree
    >>> tree = Tree('(((A_0:5,(B_0:3, C_1:3))1:2,(D_1:2, E_1:2)1:5)1:2, (F_2:3, G_2:3)1:6);')
    >>> import ecophylo as eco
    >>> phylo = eco.toPhylo(tree, 0.5, seed = 42)
    >>> print(phylo)
    <BLANKLINE>
          /-sp1
       /-|
    --|   \-sp2
      |
       \-sp3
    >>> getDeme(phylo)
    [[2, 1, 0], [0, 2, 0], [0, 0, 2]]
    >>> getDeme(phylo, div = True)
    [1, 2, 1]
    """
    # Idiot proof
    if tree.__class__.__name__ != 'TreeNode' :
        raise ValueError('tree must have a class TreeNode')
    
    indiv = list()
    for leaf in tree.iter_leaves():
        try:
            indiv.append(leaf.popInd)
        except AttribueError :
            indiv.append(1)
    if div:
        indiv = np.array(indiv)
        indiv = [sum(indiv[:,i] > 0) for i in range(indiv.shape[1])]
        # TODO : t(indiv) and compute div 
    return indiv