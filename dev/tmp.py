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

import numpy
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

m = numpy.array(migr)
print(m.shape)