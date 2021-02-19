# -*- coding: utf-8 -*-
"""
tophylo
Created on Fri Nov 6 13:20:00 2020

@author: barthele

Functions : 
    toPhylo
    ubranch_mutation

"""
# TODO : more info in help toPhylo

import numpy as np
import sys

def toPhylo(tree, mu, spmodel = "SGD", force_ultrametric = True, seed = None):
    """
    Merge branches of genealogy following speciation model of the user choice 
    after applying mutation to the tree.

    Parameters
    ----------
    tree : TreeNode (ete3 class)
        A tree with all the individuals.
    mu : float
        mutation rate, must be comprised between between 0 and 1.
    spmodel = "SGD" : string
        # TODO :DESCRIPTION
        choices SGD, NTD ; type of model of speciation wanted
        NTD -> hubbel speciation (point mutation)
        SGD -> manseau & al 2015 (with a tau != 1 is different model)
    force_ultrametric = True : bool
        msprime tree are not ultrametric by default so here is the correction.
    seed = None : int
        None by default, set the seed for mutation random events.

    Returns
    -------
    Tree Node (ete3 class)
      Individuals are merged based on the mutation present to represent species.
      A species population size can be assessed by the getAbund function.

    Examples
    --------
    >>> from ete3 import Tree
    >>> tree = Tree('(((A:5,(B:3, C:3))1:2,(D:2, E:2)1:5)1:2, (F:3, G:3)1:6);')
    >>> print(tree)
    <BLANKLINE>
             /-A
          /-|
         |  |   /-B
         |   \-|
       /-|      \-C
      |  |
      |  |   /-D
    --|   \-|
      |      \-E
      |
      |   /-F
       \-|
          \-G
    >>> phylo = toPhylo(tree, 0.5, seed = 42)
    >>> print(phylo)
    <BLANKLINE>
          /-A
       /-|
    --|   \-D
      |
       \-F
    >>> import ecophylo as eco
    >>> eco.getAbund(phylo, 7)
    [3, 2, 2]
    
    >>> toPhylo(tree = 'bamboo', mu = 0.5, seed = 42)
    Traceback (most recent call last):
      ...
    SystemExit: tree must have a class TreeNode
    
    >>> toPhylo(tree = tree, mu = -1, seed = 42)
    Traceback (most recent call last):
      ...
    SystemExit: mu must be a float between 0 and 1
    
    >>> toPhylo(tree = tree, mu = 0.5, spmodel = 'creationism', seed = 42)
    Traceback (most recent call last):
      ...
    SystemExit: creationism is not a correct model. spmodel must be either "SGD" or "NTD" string
    
    >>> toPhylo(tree = tree, mu = 0.5, force_ultrametric = 'ça fait plaisir', seed = 42)
    Traceback (most recent call last):
      ...
    SystemExit: force_ultrametric must be a boolean
    
    >>> tree_er = Tree('((A:4,(B:4,C:3)1:4)1:5, Err:1, (D:1,E:1)1:1);')
    >>> print(tree_er)
    <BLANKLINE>
          /-A
       /-|
      |  |   /-B
      |   \-|
      |      \-C
    --|
      |--Err
      |
      |   /-D
       \-|
          \-E
    >>> eco.toPhylo(tree_er, 0.5, seed = 42)
    Traceback (most recent call last):
      ...
    SystemExit: The algorithm does not know how to deal with non dichotomic trees!
    """
    # Idiot proof
    if tree.__class__.__name__ != 'TreeNode' :
        sys.exit('tree must have a class TreeNode')
    if mu < 0 or mu > 1 or not isinstance(mu, (int,float)):
        sys.exit('mu must be a float between 0 and 1')
    if not spmodel in ['SGD', 'NTD']:
        sys.exit(spmodel+' is not a correct model. '+
                'spmodel must be either "SGD" or "NTD" string')
    if not isinstance(force_ultrametric, bool):
        sys.exit('force_ultrametric must be a boolean')
    if not isinstance(seed, int):
        sys.exit('seed must be an integer')

    # init some parameters
    innerNodeIndex = 0
    nIndsORI = 0
    spID = 0
    
    # mutation model on branches
    for node in tree.traverse("preorder"): # traverse les noeuds
        try:
            node.sp
        except AttributeError:
            node.add_features(sp=1)

        if not node.is_leaf():
            node.name = "n%d" % innerNodeIndex
            innerNodeIndex += 1
        else:
            nIndsORI += 1

        if not node.is_leaf():
            umut = ubranch_mutation(node, mu, seed = seed)
            if umut:
                # print(f"Speciation event @ node {node.name}")
                spID += 1
                node.sp = spID
                for leaf in node:
                    try:
                        leaf.sp = spID
                    except AttributeError:
                        leaf.add_features(sp=1)
            # print(f"node {innerNodeIndex} --> sp: {node.sp}")

    # merging the branches with different models
    if spmodel == "NTB" :
        traversedNodes = set()
        for node in tree.traverse("postorder"):
            if not node.is_leaf():
                children = node.get_children()
                csp = [i.sp for i in children]
                if csp.count(csp[0]) == len(csp):
                    mergedLeaves = ""
                    for childnode in node.traverse():
                        traversedNodes.add(childnode)
                        if childnode.is_leaf():
                            mergedLeaves = mergedLeaves+ " "+ childnode.name

                    children[1].mergedInd = mergedLeaves
                    children[0].delete()
        
        
        
    if spmodel == "SGD" : 
        traversedNodes = set()
        for node in tree.traverse("preorder"):
            if node not in traversedNodes:
                if not node.is_leaf():
                    children = node.get_children()
                    if len(children) != 2:
                        sys.exit("The algorithm does not know how to deal with non dichotomic trees!")
                    csp1 = set()
                    for j in children[0].iter_leaves():
                        csp1.add(j.sp)
                    csp2 = set()
                    for j in children[1].iter_leaves():
                        csp2.add(j.sp)
                    common = csp1.intersection(csp2)  # compares species label between children
                    if len(common) > 0:  # if paraphyletic
                        if not node.is_root():
                            upNode = node.up  # parent
                            newLeaf = node.get_farthest_leaf()  # finds new leaf
                            newDist = newLeaf[1] + node.dist

                            mergedLeaves = ""

                            for childnode in node.traverse():
                                traversedNodes.add(childnode)
                                if childnode.is_leaf():
                                    mergedLeaves = mergedLeaves+" "+childnode.name
                            node.detach()
                            upNode.add_child(newLeaf[0], newLeaf[0].name, newDist)
                            newLeaf[0].mergedInd = mergedLeaves
                    
                        else:
                            # populate "mergedInd" feature for future SFS
                            mergedLeaves = ""
                            for l in node.iter_leaves():
                                mergedLeaves = mergedLeaves+" "+l.name
                            # collapse the subtree
                            newLeaf = tree.get_farthest_leaf()
                            for child in tree.get_children():
                                child.detach()
                                node.add_child(newLeaf[0], newLeaf[0].name, newLeaf[1])
    
                            # actualize "mergedInd" feature of new leaf
                            newLeaf[0].mergedInd = mergedLeaves
    
    if force_ultrametric: # TODO : add is.ultramtric from ete3
        tree_dist = tree.get_farthest_leaf()[1]
        for leaf in tree.iter_leaves():
            dst = tree.get_distance(leaf)
            if dst != tree_dist:
                leaf.dist += tree_dist - dst
    
    return tree


def ubranch_mutation(node, mu, seed = None):
    """
    Draw mutations following a poisson process.
    # TODO : add tau parameter for speciation
    # TODO : 

    Parameters
    ----------
    node : ete3.coretype.tree.TreeNode
        node from which to compute branch length
    mu : float
        mutation rate, must be comprised between between 0 and 1.
    seed : int
        None by default, set the seed for mutation random events.
        
    Returns
    -------
    bool
        whether or not a mutation should appear on the tree at this node

    Examples
    -------
    >>> from ete3 import Tree
    >>> tree = Tree('((A:1,(B:1,C:1)1:1)1:5,(D:1,E:1)1:1);')
    >>> node = tree.children[0] # first non-root node
    
    >>> ubranch_mutation(node = node, mu = 0, seed = 42)
    False
    
    >>> ubranch_mutation(node = node, mu = 1, seed = 42)
    True
    
    >>> ubranch_mutation(node = node, mu = 0.5, seed = 42)
    True
    
    >>> ubranch_mutation(node = 'bamboo', mu = -1, seed = 42)
    Traceback (most recent call last):
      ...
    SystemExit: node must have a class TreeNode
    
    >>> ubranch_mutation(node = node, mu = -1, seed = 42)
    Traceback (most recent call last):
      ...
    SystemExit: mu must be a float between 0 and 1
    """
    # Idiot proof
    if node.__class__.__name__ != 'TreeNode' :
        sys.exit('node must have a class TreeNode')
    if mu < 0 or mu > 1 or not isinstance(mu, (int,float)):
        sys.exit('mu must be a float between 0 and 1')
    if not isinstance(seed, int):
        sys.exit('seed must be an integer')
    
    lambd = node.dist * mu # modify node.dist -> max((node.dist - tau), 0)
    # set the seed
    np.random.seed(seed)
    rb = np.random.poisson(lambd) 
    return rb >= 1 # parametrize the 1 by a value n


if __name__ == "__main__":
        import doctest
        doctest.testmod()
