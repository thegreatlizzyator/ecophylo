# -*- coding: utf-8 -*-

# import msprime
# from ete3 import Tree
# import ecophylo

# sample_size = [10]
# com_size = [[1e5]]
# seed = 42

# #demography = None

# ecophylo.mergesizes2rates([[1000], [500]], [[100], [40]]; [500, 300], False)

# demography = [msprime.PopulationParametersChange(time=0, initial_size = com_size[0][0], growth_rate = 0, population_id= 0)]

# treeseq = msprime.simulate(sample_size= sample_size[0],
#                  Ne = com_size[0][0],
#                  random_seed= seed,
#                  demographic_events = demography)

# tree = treeseq.first()
# print(tree.draw(format = 'ascii'))

# import ecophylo
# from ete3 import Tree
# tree = Tree('(((A:5,(B:3, C:3))1:2,(D:2, E:2)1:5)1:2, (F:3, G:3)1:6);')
# print(tree)
    # <BLANKLINE>
    #          /-A
    #       /-|
    #      |  |   /-B
    #      |   \-|
    #    /-|      \-C
    #   |  |
    #   |  |   /-D
    # --|   \-|
    #   |      \-E
    #   |
    #   |   /-F
    #    \-|
    #       \-G
# phylo = ecophylo.toPhylo(tree, 0.5, seed = 42)
# print(phylo)
    # <BLANKLINE>
    #       /-A
    #    /-|
    # --|   \-D
    #   |
    #    \-F
# ecophylo.getAbund(phylo, 7)

# import msprime as ms
# from ete3 import Tree
# import ecophylo as eco
# seed = 42
# # pc = [ms.PopulationConfiguration(sample_size = s, initial_size = i, growth_rate = g) for s, i, g in zip([5], [500], [0])]
# # cas avec très petit taux de migr
# # treeseq=ms.simulate(population_configurations = pc, random_seed= seed)

# pc = [ms.PopulationConfiguration(sample_size = s, initial_size = i, growth_rate = g) for s, i, g in zip([5, 3, 4], [500, 500, 500], [0,0, 0])]
# # cas avec taux de migr
# treeseq=ms.simulate(population_configurations = pc, random_seed= seed, migration_matrix = [[0, 0.0001, 0.0001], [0.0001, 0, 0.0001], [0.0001, 0.0001, 0]])

# tree = treeseq.first()
# print(tree.draw(format = 'ascii'))

# # tree.population(18)
# # print( [tree.population(i) for i in range(10)] )
# node_labels = {u: str(u)+'_'+str(tree.population(u)) for u in tree.nodes() if tree.is_sample(u)}
# tree = Tree(tree.newick(node_labels = node_labels))
# print(tree)

# phylo = eco.toPhylo(tree, 0.5, seed = 42)

# print(phylo)

# print(eco.getAbund(phylo))

# print(eco.getDeme(phylo))

import numpy as np
migr = [[[0,0.1],
        [0.2,0]]
        ,
        [[0,0.3],
        [0.4,0]]
        ,
        [[0,0.3],
        [0.4,0]]
        ,
        [[0,0.3],
        [0.4,0]]
        ,
        [[0,0.5],
        [0.6,0]]]

m = np.array(migr)
print(m.shape)

npop = 2
migr = 0.4
print(np.ones((npop,npop))*migr)

migr = [[0, 1],[1, 0, 0]]
migr_time = [0]
npop = 2

migr = 1

import sys

if True :
    # check migr
    if migr is not None :
        if npop == 1 :
            # warnings.warn("no migration matrix is needed for a single deme")
            migr = None
    if migr is not None :
        if not isinstance(migr, list) : # case 'a
            if not isinstance(migr, (int,float)) :
                sys.exit("migration rate must be a float or an int.")
            if migr < 0 or migr > 1 :
                sys.exit("migration rate should be positive (or zero) and" + 
                " not exceed 1")
            # migr = np.ones((npop,npop))*migr
            # np.fill_diagonal(migr, 0)
        else :
            for i in range(len(migr)):
                if not isinstance(migr[i], list): # case ['a, ... ,'b]
                    if len(migr) != len(migr_time):
                        sys.exit("there should be as many migration rates" + 
                            " or matrices as there are times in migr_time")
                    if not isinstance(migr[i], (int,float)) :
                        sys.exit("migration rate must be a float or an int.")
                    if migr[i] < 0 or migr[i] > 1 :
                        sys.exit("migration rate should be positive (or zero)" + 
                                 " and not exceed 1")
                    # check len of migr is done with migr_time
                    migr[i] = np.ones((npop,npop))*migr[i]
                    np.fill_diagonal(migr[i], 0)
                else :
                    if not isinstance(migr[i][0], list) : # case [[0,'a], ['a, 0]]
                        if len(migr[i]) != len(migr) or len(migr) != npop:
                            sys.exit("custom migration matrices should be of" + 
                                     " size ndeme x ndeme")
                        if not all([ isinstance(r, (float,int)) for r in migr[i]]) :
                            sys.exit("found custom migration matrix that is" + 
                                     " not made of ints or floats")
                        if any([r<0 or r> 1 for r in migr[i]]):
                            sys.exit("found custom migration matrix with" + 
                                " negative migration rates or greater than 1")
                    else: # case [[[0, 'a], ['a, 0]], [[0, 'b], ['b, 0]]]
                        if len(migr) != len(migr_time):
                            sys.exit("there should be as many migration rates" + 
                                " or matrices as there are times in migr_time")
                        for j in range(len(migr[i])) :
                            if len(migr[i][j]) != len(migr[i]) or len(migr[i]) != npop:
                                sys.exit("custom migration matrices should be" + 
                                    " of size ndeme x ndeme")
                            if any ([not isinstance(r, (float,int)) for r in migr[i][j]]) :
                                sys.exit("found custom migration matrix that" + 
                                    " is not made of ints or floats")
                            if any ([r<0 or r> 1 for r in migr[i][j]]):
                                sys.exit("found custom migration matrix with" + 
                                 " negative migration rates or greater than 1")

    m = np.array(migr)
    dim = m.shape
    if np.sum(m) == 0 :
        sys.exit("migration matrices cannot all be empty")

print(dim)
print(len(dim))
