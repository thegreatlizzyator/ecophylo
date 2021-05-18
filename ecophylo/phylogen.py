# -*- coding: utf-8 -*-
"""
tophylo
Created on Fri Nov 6 13:20:00 2020

@author : Maxime Jaunatre <maxime.jaunatre@yahoo.fr>
@author : Elizabeth Bathelemy <barthelemy.elizabeth@gmail.com>

Functions : 
    toPhylo
    ubranch_mutation

"""
# TODO : more info in help toPhylo

import numpy as np

def toPhylo(tree, mu, tau = 0, spmodel = "SGD", 
            force_ultrametric = True, seed = None):
    """
    Merge branches of genealogy following speciation model of the user choice 
    after sprinkling mutation events over the branches of simulated genealogies
    depending on branch lengths. 
    
    Mutation events are sprinkled over the branches of simulated genealogies 
    depending on branch lengths, so that the number of mutations over a branch 
    follows a Poisson distribution with parameter 𝜇·𝐵 where 𝜇 is the point 
    mutation rate and 𝐵 is the length of the branch. 
    
    The descendants stemming from a branch with at least one mutation define
    a genetically distinct clade. Since an extant species should be a 
    monophyletic genetic clade distinct from other species, all paraphyletic 
    clades of haplotypes at present are merged to form a single species. 
    
    Parameters
    ----------
    tree : TreeNode (ete3 class)
        A tree representing the genealogy of simulated individuals.
    mu : float
        point mutation rate, must be comprised between between 0 and 1. 
    tau = 0 : float
        The minimum number of generations monophyletic lineages have to be 
        seperated for to be considered distinct species
    spmodel = "SGD" : string
        the type of speciation model to implement. Default if "SGD" and 
        corresponds to a generalisation of the Speciation by Genetic 
        Differentiation. Note that setting tau to 1 will equate to the SGD 
        model as described in Manceau et al. 2015.
        "NTB" corresponds to the speciation model as described in Hubbell 2001, 
        in which point mutations instantenously give rise to new species.
    force_ultrametric = True : bool
        Whether or note to force phylogenetic tree ultrametry 
    seed = None : int
        None by default, set the seed for mutation random events.

    Returns
    -------
    Tree Node (ete3 class)
        A phylogeny representing the phylogenetic relationships among species 
        as well as the number of individuals descending from a speciation event
        in the genealogy, which defined the species abundance in the sample at 
        present (abundances can be retrived using the getAbund function).

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
          /-sp1
       /-|
    --|   \-sp2
      |
       \-sp3
    >>> import ecophylo as eco
    >>> eco.getAbund(phylo, 7)
    [3, 2, 2]
    """
    # Idiot proof
    if tree.__class__.__name__ != 'TreeNode' :
        raise ValueError('tree must have a class TreeNode')
    if mu < 0 or mu > 1 or not isinstance(mu, (int,float)):
        raise ValueError('mu must be a float between 0 and 1')
    if not spmodel in ['SGD', 'NTD']:
        raise ValueError(spmodel+' is not a correct model. '+
                'spmodel must be either "SGD" or "NTD" string')
    if not isinstance(force_ultrametric, bool):
        raise ValueError('force_ultrametric must be a boolean')
    if seed is not None and not isinstance(seed, int):
        raise ValueError('seed must be an integer')

    # init some parameters
    innerNodeIndex = 0
    nIndsORI = 0
    spID = 0
    demeID = 0
    ndeme = 0
    
    # mutation model on branches
    for node in tree.traverse("preorder"): # traverse les noeuds
        try:
            node.sp
        except AttributeError:
            node.add_features(sp=1)
        try:
            node.deme
        except AttributeError:
            node.add_features(deme=1)

        if not node.is_leaf():
            node.name = "n%d" % innerNodeIndex
            innerNodeIndex += 1
        else:
            nIndsORI += 1
            name_deme = node.name.split("_")
            if len(name_deme) == 1: # if no population added, only one deme
                name_deme.append('0')
            if(int(name_deme[1]) > ndeme):
                ndeme += 1
            node.deme = int(name_deme[1])
            node.name = name_deme[0]

        if not node.is_leaf():
            umut = ubranch_mutation(node= node, mu= mu, tau= tau, seed= seed)
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
                    popInd = [0] * (ndeme + 1)
                    for childnode in node.traverse():
                        traversedNodes.add(childnode)
                        if childnode.is_leaf():
                            mergedLeaves = mergedLeaves+ " "+ childnode.name
                            popInd = [a+b for a, b in zip(childnode.popInd, popInd)]

                    children[1].mergedInd = mergedLeaves
                    children[0].delete()
            else :
                popInd = [0] * (ndeme + 1)
                popInd[node.deme] += 1
                node.popInd = popInd

        
        
        
    if spmodel == "SGD" : 

        for leaf in tree.iter_leaves():
            popInd = [0] * (ndeme + 1)
            popInd[leaf.deme] += 1
            leaf.popInd = popInd

        traversedNodes = set()
        for node in tree.traverse("preorder"):
            if node not in traversedNodes:
                if not node.is_leaf():
                    children = node.get_children()
                    if len(children) != 2:
                        raise ValueError("The algorithm does not know how to"+
                                         " deal with non dichotomic trees!")
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
                            popInd = [0] * (ndeme + 1)

                            for childnode in node.traverse():
                                traversedNodes.add(childnode)
                                if childnode.is_leaf():
                                    mergedLeaves = mergedLeaves+" "+childnode.name
                                    popInd[childnode.deme] += 1
                            node.detach()
                            upNode.add_child(newLeaf[0], newLeaf[0].name, newDist)
                            newLeaf[0].mergedInd = mergedLeaves
                            newLeaf[0].popInd = popInd
                    
                        else:
                            # populate "mergedInd" feature for future SFS
                            mergedLeaves = ""
                            popInd = [0] * (ndeme + 1)
                            for l in node.iter_leaves():
                                mergedLeaves = mergedLeaves+" "+l.name
                                popInd = [a+b for a, b in zip(l.popInd, popInd)]
                            # collapse the subtree
                            newLeaf = tree.get_farthest_leaf()
                            for child in tree.get_children():
                                child.detach()
                                node.add_child(newLeaf[0], newLeaf[0].name, newLeaf[1])
    
                            # actualize "mergedInd" feature of new leaf
                            newLeaf[0].mergedInd = mergedLeaves
                            newLeaf[0].popInd = popInd
    
    if force_ultrametric: # TODO : add is.ultramtric from ete3
        tree_dist = tree.get_farthest_leaf()[1]
        for leaf in tree.iter_leaves():
            dst = tree.get_distance(leaf)
            if dst != tree_dist:
                leaf.dist += tree_dist - dst

    nsp = 1
    for leaf in tree.iter_leaves():
        leaf.name = "sp"+str(nsp)
        nsp += 1
    
    return tree


def ubranch_mutation(node, mu, tau = 0, seed = None):
    """
    Draw mutations following a poisson process with parameter 
    max((B - tau), 0)*mu where mu is the point mutation rate, B is the 
    length of the branch at a given node and tau 
>>>
    Parameters
    ----------
    node : ete3.coretype.tree.TreeNode
        node from which to compute branch length
    mu : float
        point mutation rate, must be comprised between between 0 and 1. 
    tau = 1 : float
        The minimum number of generations monophyletic lineages have to be 
        seperated for to be considered distinct species
    seed : int
        None by default, set the seed for mutation random events.
        
    Returns
    -------
    bool
        whether or not a at least one mutation should appear on the tree at 
        this node

    Examples
    -------
    TODO: Examples with tau ?
    >>> from ete3 import Tree
    >>> tree = Tree('((A:1,(B:1,C:1)1:1)1:5,(D:1,E:1)1:1);')
    >>> node = tree.children[0] # first non-root node
    
    >>> ubranch_mutation(node = node, mu = 0, seed = 42)
    False
    
    >>> ubranch_mutation(node = node, mu = 1, seed = 42)
    True
    
    >>> ubranch_mutation(node = node, mu = 0.5, seed = 42)
    True
    """
    # Idiot proof
    if node.__class__.__name__ != 'TreeNode' :
        raise ValueError('node must have a class TreeNode')
    if mu < 0 or mu > 1 or not isinstance(mu, (int,float)):
        raise ValueError('mu must be a float between 0 and 1')
    if tau < 0 or not isinstance(tau, (int,float)):
        raise ValueError('tau must be a float superior or equal to 0')
    if seed is not None and not isinstance(seed, int):
        raise ValueError('seed must be an integer')
    
    lambd = max((node.dist - tau), 0) * mu 
    # set the seed
    np.random.seed(seed)
    rb = np.random.poisson(lambd) 
    return rb >= 1 # parametrize the 1 by a value n


if __name__ == "__main__":
        import doctest
        doctest.testmod()
