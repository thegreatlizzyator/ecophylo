# -*- coding: utf-8 -*-

import msprime
import ecophylo
from ete3 import Tree

# seed = 42
# mu = 0.03

# demography = msprime.Demography()
# demography.add_population(name="A", initial_size=100_000)
# print(demography.debug())
# print(dir(demography))
# treeseq = msprime.sim_ancestry(samples = {"A": 10}, 
#     demography=demography, random_seed=seed, ploidy = 1)
# # treeseq = msprime.sim_ancestry(2, ploidy=1)

# tree = treeseq.first()
# # print(tree.draw(format = 'ascii'))

# node_labels = {u: str(u)+'_'+str(tree.population(u)) for u in tree.nodes() if tree.is_sample(u)}
# tree = Tree(tree.newick(node_labels = node_labels))
# # print(tree)
# phylo = ecophylo.toPhylo(tree, mu, seed = seed)
# # print(phylo)
        
# t = ecophylo.simulate_dolly(
#     sample_size = [10], 
#     com_size = [[1e5]], 
#     mu = mu, seed = seed, verbose = True)
# # print(t)

demography = msprime.Demography()
demography.add_population(name="A", initial_size=10_000)
demography.add_population(name="B", initial_size=5_000)
demography.add_population(name="C", initial_size=1_000)
demography.add_population_split(time=1000, derived=["A", "B"], ancestral="C")
ts = msprime.sim_ancestry(samples={"A": 1, "B": 1}, demography=demography, random_seed=12)
ts
print(ts.first().draw(format = 'ascii'))