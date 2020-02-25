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
Module responsible for simulating a phylogeny from a given set of parameters.

@author: Elizabeth Barthelemy

notes: in the simplest of cases (with default values): one population, no past changes - constant size - SGD? 

"""

import msprime
import numpy as np
import sys
from ete3 import Tree

import convert
import generate


def simulate(sample_size, com_size, mu, npop = 1, nepoch = 1, m = 0, init_rates = None, init_sizes = None, past_sizes = None, maxtime = None, split_dates = None, migrfrom = None, migrto = None, seed = 1, verbose = False):

    # do dummy checks here --> try to make code stupid-proof
    if sample_size >= com_size:
        sys.exit("Sample size should not exceed community size")
    
    popchange = []
    massmigration = []
    popconfig = None
    migration = None

    # make past demographic changes between different time frames
    if nepoch > 1:
        if len(past_sizes) != nepoch :
            sys.exit("There should be as many past growth rates as there are epochs")
        times = generate.timeframes(nepoch, maxtime, 0.05)
        popchange = generate.demographic_events(times, past_sizes)

    # make island model
    if npop > 1:
        if init_sizes is None or init_rates is None:
            sys.exit("Initial population sizes and growth rates should be provided when there are more than one population (npop>1)")
        migration = generate.migration_matrix(npop, m)
        samples = np.ones(npop, dtype=int)*sample_size
        popconfig = generate.population_configurations(samples, init_sizes, init_rates)

        # possible mass migration between populations
        if split_dates is not None:
            # implement option later for limited mass dispersal
            M = 1
            massmigration = generate.mass_migrations(split_dates, migrfrom, migrto, M)

    demography = popchange + massmigration

    if len(demography) == 0:
        demography = None

    # if verbose should print the demography debugger - only for debugging purposes!!! 
    if verbose:
        dd = msprime.DemographyDebugger(demographic_events= demography, migration_matrix= migration, population_configurations= popconfig)
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
    #print(tree.draw(format="unicode"))
    tree = Tree(tree.newick())
    phylo = convert.toPhylo(tree, mu)

    return phylo   
