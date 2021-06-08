Simulation API
==============
We introduce the Python package ecophylo dedicated to coalescent-based simulation of eco-evolutionary dynamics. The model rests on coalescent theory to simulate a shared ancestry of co-occurring individuals, under the influence of past demographic fluctuations due to for example, habitat fluctuations, fragmentation and/or migration events among separate areas. Mutations occur over time in the genealogy, and divergent genotypes represent distinct extant species.

First, we start by simulating the genealogies of individuals in an assemblage experiencing past demographic fluctuations and/or linked by vicariance or migration events using the ms coalescent simulator, without allotting a species label to each individual. Second, we sprinkle mutation events over the branches of simulated genealogies depending on branch lengths. Since an extant species should be a monophyletic genetic clade distinct from other species, all paraphyletic clades of haplotypes at present were merged to form a single species.

The package includes tools to simulate large numbers of datasets and associated summary statistics, so that Approximate Bayesian Computation methods can be used to estimate parameter values for these processes. Diverse patterns of taxonomic and phylogenetic compositions can be generated.


Core simulation algorithm
-------------------------
The ecophylo.simulate function implements the core simulation algorithm mentionned above. It yields a phylogeny of species belonging to an extant assemblage in Newick format, along with corresponding species abundances, for a combination of parameter values representing the past biogeographic history of Jm(t). Each biogeographic event accounting for the history of Jm(t) occurs at a specific time that users must supply in changetime as a list in the order in which they occur.

.. autofunction:: ecophylo.simulate


Running multiple simulations
----------------------------
Wether to conduct virtual experiments examining how different eco-evulionary scenarios have shaped patterns of diversity (Barthelemy et al. 2021) or to carry out inference the ecophylo.dosimul function can be used to call ecophylo.simulate and retreive summary statistics generated for different parameter values drawn from specified distributions.

.. autofunction:: ecophylo.dosimuls


