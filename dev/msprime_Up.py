# -*- coding: utf-8 -*-

import msprime
import ecophylo

seed = 42
mu = 0.03

treeseq = msprime.sim_ancestry(
        samples=10,
        population_size=100_000,
        random_seed=seed)

tree = treeseq.first()
print(tree.draw(format = 'ascii'))

node_labels = {u: str(u) for u in tree.nodes() if tree.is_sample(u)}
tree = Tree(tree.newick(node_labels = node_labels))
phylo = phylogen.toPhylo(tree, mu, seed = seed)
print(phylo)
        
t = ecophylo.simulate_dolly(
    sample_size = [10], 
    com_size = [[1e5]], 
    mu = mu, seed = seed)
print(t)