
import ecophylo as eco

import msprime
import random
import numpy as np
import sys
from ete3 import Tree
import pandas as pd
from loguniform import LogUniform

sample_size = 10
com_size = 1000
mu = 0.03

mrca = None ; npop = 1 ; m = 0 ; init_rates = None ; init_sizes = None ; 
past_sizes = None ; changetime = None ; split_dates = None ; migrfrom = None ;
migrto = None ; verbose = False ; seed = 42

# TODO : doc !!!
# TODO : idiotproof
# TODO : more examples
# do dummy checks here --> try to make code stupid-proof
if sample_size >= com_size:
    sys.exit("Sample size should not exceed community size")

popchange = []
massmigration = []
popconfig = None
migration = None

# make past demographic changes between different time frames
if past_sizes is not None and changetime is not None:
    if len(past_sizes) != len(changetime):
        sys.exit("There should be as many sizes as there are past epochs")
    popchange = demographic_events(changetime, past_sizes)

# make island model
if npop > 1:
    if init_sizes is None or init_rates is None:
        sys.exit("Initial population sizes and growth rates should be provided when there are more than one population (npop>1)")
    migration = islmodel.migration_matrix(npop, m)
    samples = np.ones(npop, dtype=int)*sample_size
    popconfig = islmodel.population_configurations(samples, init_sizes, init_rates)

    # possible mass migration between populations
    if split_dates is not None:
        # implement option later for limited mass dispersal
        M = 1
        massmigration = islmodel.mass_migrations(split_dates, migrfrom, migrto, M)

demography = popchange + massmigration
if len(demography) == 0:
    demography = None

# if verbose should print the demography debugger - only for debugging purposes!!! 
if verbose:
    dd = msprime.DemographyDebugger(Ne = com_size, 
                                    demographic_events= demography, 
                                    migration_matrix= migration, 
                                    population_configurations= popconfig)
    dd.print_history()

if npop > 1:
     treeseq = msprime.simulate(Ne = com_size,
                                random_seed= seed,
                                population_configurations = popconfig,
                                migration_matrix = migration,
                                demographic_events = demography)
else: 
    treeseq = msprime.simulate(sample_size= sample_size,
                               Ne = com_size,
                               random_seed= seed,
                               demographic_events = demography)
                               
tree = treeseq.first()
if verbose: print(tree.draw(format = 'unicode'))
if mrca is not None:
    if tree.time(tree.root) > mrca : 
        raise Exception(f"Simulated MRCA ({tree.time(tree.root)}) predates"+
                         " fixed limit ({mrca})")
#print(tree.draw(format="unicode"))
node_labels = {u: str(u) for u in tree.nodes() if tree.is_sample(u)}
tree = Tree(tree.newick(node_labels = node_labels))
phylo = eco.toPhylo(tree, mu, seed = seed)

print(phylo)



