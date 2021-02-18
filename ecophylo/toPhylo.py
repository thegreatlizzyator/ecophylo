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
    Merge branches of genealogy following speciation model of the user choice.

    Parameters
    ----------
    tree : TreeNode (ete3 class)
        # TODO : DESCRIPTION.
        Is a genealogy
    mu : float
        # TODO :DESCRIPTION.
        0 : 1
    spmodel : string
        # TODO :DESCRIPTION
        choices SGD, NTD ; type of model of speciation wanted
        NTD -> hubbel speciation (point mutation)
        SGD -> manseau & al 2015 (with a tau != 1 is different model)
    force_ultrametric : bool
        True by default
        # TODO : check if option is usefull
    seed : int
        None by default, set the seed for mutation random events.

    Returns
    -------
    Tree Node (ete3 class)
    Is a species tree

    Examples
    --------
    >>> print("test")
    "no"
    """
    # TODO : example toPhylo
    # TODO : idiot proof toPhylo
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
    
    if force_ultrametric:
        tree_dist = tree.get_farthest_leaf()[1]
        for leaf in tree.iter_leaves():
            dst = tree.get_distance(leaf)
            if dst != tree_dist:
                leaf.dist += tree_dist - dst
    
    return tree


def ubranch_mutation(node, mu, seed = None):
    """
    Draw mutations following a poisson process.

    Parameters
    ----------
    node : ete3.coretype.tree.TreeNode
        node from which to compute branch length
    mu : float
        # TODO :DESCRIPTION.
        0 : 1
    seed : int
        None by default, set the seed for mutation random events.
    # TODO : add tau parameter for speciation

    Returns
    -------
    bool
        whether or not a mutation should appear on the tree at this node

    Examples
    -------
    >>> from ete3 import Tree
    >>> tree = Tree()
    >>> tree.populate(5)
    >>> node = tree.get_tree_root()
    
    >>> ubranch_mutation(node = node, mu = 0, seed = 42)
    False
    
    >>> ubranch_mutation(node = node, mu = 1, seed = 42)
    False
    
    >>> ubranch_mutation(node = node, mu = 0.2, seed = 42)
    False
    
    >>> ubranch_mutation(node = node, mu = -1, seed = 42)
    sys.exit('mu must be a float between 0 and 1')
    """
    
    if mu < 0 or mu > 1 or not isinstance(mu, float):
        sys.exit('mu must be a float between 0 and 1')
    
    # TODO : example ubranch_mutation
    # TODO : idiot proof ubranch_mutation
    lambd = node.dist * mu # modify node.dist -> max((node.dist - tau), 0)
    # set the seed
    np.random.seed(seed)
    rb = np.random.poisson(lambd) 
    return rb >= 1 # parametrize the 1 by a value n


if __name__ == "__main__":
        import doctest
        doctest.testmod()
