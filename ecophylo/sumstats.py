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
Module responsible for computing summary statistics on the simulated 
phylogenies

@author: Elizabeth Barthelemy

"""
def getSFS(tree):
    """
    

    Parameters
    ----------
    tree : TYPE
        DESCRIPTION.

    Returns
    -------
    sfs : TYPE
        DESCRIPTION.

    """
    sfs = list()
    abund = list()
    for leaf in tree.iter_leaves():
        try:
            inds = leaf.mergedInd.lstrip(" ")
            inds = list(inds.split(" "))
            abund.append(len(inds))
        except AttributeError:
            abund.append(1)
    abund.sort()
    sfs.extend(abund)
    return sfs
