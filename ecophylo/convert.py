# -*- coding: utf-8 -*-
#
# This file is a part of the ecophylo package.
# 
# ecophylo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the CeCILL FREE SOFTWARE LICENSE AGREEMENT
# along with ecophylo.

"""
Module responsible for converting genealogies into phylogenies.

@author: Elizabeth Barthelemy

"""

import numpy as np
import sys

def toPhylo(tree, mu, spmodel = "SGD"):
    """
    

    Parameters
    ----------
    tree : TYPE
        DESCRIPTION.
    mu : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    innerNodeIndex = 0
    nIndsORI = 0
    spID = 0
    
    for node in tree.traverse("preorder"):
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
            umut = ubranch_mutation(node, mu)
            if umut:
                # print(f"Speciation event @ node {innerNodeIndex}")
                spID += 1
                node.sp = spID
                for leaf in node:
                    try:
                        leaf.sp = spID
                    except AttributeError:
                        leaf.add_features(sp=1)
            # print(f"node {innerNodeIndex} --> sp: {node.sp}")

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
    return tree

def ubranch_mutation(node, mu):
    """
    Draw mutations following a poisson process.

    Parameters
    ----------
    node : ete3.coretype.tree.TreeNode
        node from which to compute branch length
    mu : float
        mutation rate

    Returns
    -------
    bool
        whether or not a mutation should appear on the tree at this node

    """
    lambd = node.dist * mu
    rb = np.random.poisson(lambd)
    if rb >= 1:
        return True
    else:
        return False